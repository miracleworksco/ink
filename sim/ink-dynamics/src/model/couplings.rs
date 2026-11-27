use super::{Input, M3, Model, Moment, Num, ROUNDING, Rides, V3, Wrench, mul, mv, rot, split, tr, vee};
use ink_core::eval::Ctx;
use ink_core::v3::{cross, dot, norm};

/// Where a coupling's loss is charged at a moment (Model::couple): a joint
/// a coupling or a mode passes power through, the moment or force `held`
/// its lines exert along it as the ideal couplings pass it and its rate,
/// and the efficiency `eff` of what passes from it where the power enters
/// the coupling there (held times rate below 0): a coupling's A end E, its
/// B end EB, each part of a mode E. Where the power leaves it charges
/// nothing, so a coupling is charged at the end power enters it, whichever
/// end is written first. Its loss is friction in that joint (losses):
/// where the joint balances couplings' tensions (an intermediate shaft, a
/// passive output driving back), that friction is in its balance, so the
/// stage driving it carries it, `gain` (k, g): the moment this site's lines
/// gain per unit of site k's friction, for each k that moves it (itself
/// too, where its own balance gives its lines). The effort along
/// coordinate k is less each site's friction times weights[k]. `total`:
/// the moment its lines exert with the losses past it, whose sign with the
/// rate's says whether power enters there. `worst`: the most its loss can
/// be of the power it passes (1/E - 1 at A, 1 - EB at B).
#[derive(Clone, Debug)]
pub struct Site {
    pub held: f64,
    /// Roundoff in the products summed at this site, before opposing line tensions cancel.
    pub held_roundoff: f64,
    pub rate: f64,
    pub eff: f64,
    pub worst: f64,
    pub weights: Vec<f64>,
    pub gain: Vec<(usize, f64)>,
    pub total: f64,
    /// Its friction where it is charged in a joint balancing the tensions
    /// of several couplings besides its own (a head two belts drive): it
    /// is charged there to the actuators, not carried by those stages, so
    /// the loss it would add upstream is bounded instead (loss_bound).
    pub shared: bool,
    pub unfed: f64,
    /// What passes power together here (a mode, else a coupling): a stage.
    pub stage: usize,
}

/// What a site's friction acts on (Model::couple): the joint it is along
/// (turn, unit axis and origin in the world), the parts it acts on (their
/// ends: coupling, 0 for A and 1 for B), every end of its stage, and the
/// end it pairs with, a coupling's (the other end, which it drives).
struct Acts {
    joint: (bool, V3, V3),
    parts: Vec<PartActs>,
    all: Vec<(usize, usize)>,
    partner: Option<(usize, usize)>,
}

/// A part a site's friction acts on: its ends, its lines' coefficients in
/// the site's held moment (its share of the friction is its share of
/// that), and its axle where the site's joint is one carrying its own (the
/// friction then a force there along the way that joint moves it; None: a
/// moment about the joint's axis, or a force along it through its origin).
type PartActs = (Vec<(usize, usize)>, Vec<(usize, f64)>, Option<V3>);

/// A coupling's friction at a moment as a wrench on parts (Model::couple):
/// the ends it acts on, every end of its stage, its force (N) and its
/// moment about the world's origin (N mm), and whether its opposite acts
/// on the stage's other part (paired), or on no part a reaction separates
/// (a routed belt's: how its parts' supports share it is not found).
#[derive(Clone, Debug)]
pub struct Friction {
    on: Vec<(usize, usize)>,
    all: Vec<(usize, usize)>,
    pub force: V3,
    pub moment: V3,
    pub paired: bool,
}

/// The wrench (force, moment about the world's origin) a friction x along
/// P's joint makes on the part Q it drives, its opposite on P, Newton's
/// third law: each joint (turn, unit axis, origin), Q's where it rides one.
/// Of those that give P x along its joint and Q nothing along its own, the
/// least force: a moment alone where one does (a screw's thread friction
/// on its nut, a rack's on the rack), else a force through Q's axis across
/// the line between the two axes (two wheels on parallel shafts: the
/// moment over their distance apart). None where none does (two coaxial
/// wheels, two slides along one line).
fn pair(p: (bool, V3, V3), q: Option<(bool, V3, V3)>, x: f64, level: Option<V3>) -> Option<(V3, V3)> {
    let (pt, a, op) = p;
    // A point on an axis (unit u through o) level with the lines.
    let at = |u: V3, o: V3| -> V3 {
        match level {
            Some(l) => {
                let s = (0..3).map(|k| (l[k] - o[k]) * u[k]).sum::<f64>();
                std::array::from_fn(|k| o[k] + s * u[k])
            }
            None => o,
        }
    };
    let scale = |v: V3, k: f64| v.map(|c| k * c);
    let add = |u: V3, v: V3| -> V3 { std::array::from_fn(|i| u[i] + v[i]) };
    // Two unit directions: the least c a + d b with u.a = -x, u.b = 0.
    let two = |a: V3, b: V3| -> Option<V3> {
        let c = dot(a, b);
        (1.0 - c * c > 1e-12).then(|| add(scale(a, -x / (1.0 - c * c)), scale(b, x * c / (1.0 - c * c))))
    };
    match (pt, q) {
        (true, None) | (true, Some((false, _, _))) => Some(([0.0; 3], scale(a, -x))),
        (true, Some((true, b, oq))) => match two(a, b) {
            Some(m) => Some(([0.0; 3], m)),
            None => {
                let d: V3 = std::array::from_fn(|i| oq[i] - op[i]);
                let across: V3 = add(d, scale(a, -dot(d, a)));
                let n2 = dot(across, across);
                if n2 <= 1e-18 * (1.0 + dot(d, d)) {
                    return None;
                }
                let f = scale(cross(a, across), -x / n2);
                Some((f, cross(at(b, oq), f)))
            }
        },
        (false, None) => {
            let f = scale(a, -x);
            Some((f, cross(op, f)))
        }
        (false, Some((true, b, oq))) => {
            let f = scale(a, -x);
            Some((f, cross(at(b, oq), f)))
        }
        (false, Some((false, b, _))) => two(a, b).map(|f| (f, cross(op, f))),
    }
}

/// What an actuator's couplings' loss may miss, per unit of its worst
/// compounded efficiency's share (Moment::worst): its share `p` of what the
/// unseparated bodies need along each coordinate, and of each friction
/// charged unfed (Site::unfed), which the stages driving it would carry.
/// `friction`: each site's at another state than its own (another answer's),
/// its unfed friction that where it is `shared`.
pub fn unfed_need(p: &[f64], apart: &[f64], sites: &[Site], friction: Option<&[f64]>) -> f64 {
    let unfed = |j: usize, s: &Site| match friction {
        Some(x) if s.shared => x[j],
        Some(_) => 0.0,
        None => s.unfed,
    };
    let fed: f64 = sites.iter().enumerate().filter(|(_, s)| s.shared).map(|(j, s)| (unfed(j, s) * p.iter().zip(&s.weights).map(|(p, w)| p * w).sum::<f64>()).abs()).sum();
    p.iter().zip(apart).map(|(p, a)| (p * a).abs()).sum::<f64>() + fed
}

/// Each site's friction (Site) where their lines' ideal moments are
/// `held`, from T, the moment its lines exert with the friction of the
/// other sites in the balances that give its tensions: none where power
/// does not enter there (T times its rate not below 0); else, with
/// efficiency e and g its own gain, F = (1 - e) T / (e - (1 - e) g), so that
/// what its lines pass, T + g F, is e of what enters: (1/e - 1) T where its
/// lines are set elsewhere (g = 0: a motor driving a load), (1 - e) T where
/// its own balance sets them (g = -1: a load driving the motor back; a
/// self-locking screw's e <= 0 then has the motor drive it). Each stage's
/// loss is so on the power entering it, the losses of the stages it drives
/// included, and stages in series multiply. With each site's T. Err where
/// the frictions do not settle (couplings driving each other round a loop).
pub fn losses(sites: &[Site], held: &[f64]) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = sites.len();
    let own = |j: usize| sites[j].gain.iter().find(|g| g.0 == j).map_or(0.0, |g| g.1);
    let friction = |j: usize, t: f64| {
        if t * sites[j].rate >= 0.0 {
            return 0.0;
        }
        let (e, g) = (sites[j].eff, own(j));
        let d = e - (1.0 - e) * g;
        // Its loss as a share of what enters where nothing settles it
        // (a self-locking end whose lines are set elsewhere).
        if d <= 0.0 { (1.0 - e) * t } else { (1.0 - e) * t / d }
    };
    // No site feeding another (a CoreXY's belts, a single mesh): each site's
    // friction on its own lines alone.
    if sites.iter().enumerate().all(|(j, s)| s.gain.iter().all(|g| g.0 == j)) {
        return Ok(((0..n).map(|j| friction(j, held[j])).collect(), held.to_vec()));
    }
    let mut x = vec![0.0; n];
    // Each pass carries the friction one stage further up a train, so a
    // train of n stages settles in n passes.
    let others = |j: usize, x: &[f64]| held[j] + sites[j].gain.iter().filter(|g| g.0 != j).map(|&(k, g)| g * x[k]).sum::<f64>();
    for _ in 0..2 * n + 4 {
        let next: Vec<f64> = (0..n).map(|j| friction(j, others(j, &x))).collect();
        let settled = (0..n).all(|j| (next[j] - x[j]).abs() <= 1e-12 * (held[j].abs() + next[j].abs() + x[j].abs()));
        x = next;
        if settled {
            let total = (0..n).map(|j| others(j, &x)).collect();
            return Ok((x, total));
        }
    }
    Err("the couplings' losses do not settle: couplings pass power round a loop, each one's loss charged on what the next passes back".into())
}

/// Whether two points are one, to the rounding of where they were worked
/// out (a span and a belt's path drawn by the same route).
fn same_point(a: V3, b: V3) -> bool {
    norm(std::array::from_fn(|k| a[k] - b[k])) <= 1e-6 * (1.0 + norm(a))
}

/// What each piece of a carried belt's path lies on, from its ends and the
/// spans' (coupling, line, A's end then B's): a line on the span whose two
/// ends are its own, an arc on the wheel of the part a span meets at one of
/// its ends. None where a piece lies on none: the belt is not separated.
pub(super) fn belt_pieces(pieces: &[[V3; 2]], arcs: &[bool], spans: &[(usize, usize, [V3; 2])]) -> Option<Vec<Rides>> {
    pieces
        .iter()
        .zip(arcs)
        .map(|(&[p0, p1], &arc)| {
            if arc {
                spans.iter().find_map(|&(c, l, e)| (0..2).find(|&k| same_point(p0, e[k]) || same_point(p1, e[k])).map(|k| Rides::Wrap(c, l, k)))
            } else {
                spans.iter().find(|(_, _, [a, b])| (same_point(p0, *a) && same_point(p1, *b)) || (same_point(p0, *b) && same_point(p1, *a))).map(|&(c, l, _)| Rides::Span(c, l))
            }
        })
        .collect()
}

impl Model {
    /// What the couplings pull with on the part at `on` (an end of one, or a part within an end) at moment m:
    /// the largest of their lines' tensions (N; a screw's thrust, a belt's tight span, a tooth's force) and
    /// the coupling it is of. None: no coupling has an end there. Err: the tensions are not found there, or a
    /// moving body with no mass is on a coordinate they rest on.
    /// How far `pull`'s tension can move for what the bodies' measures are known to (Moment::mass_t).
    pub fn pull_mass(&self, on: &str, m: &Moment) -> f64 {
        let within = |end: &str| on == end || on.strip_prefix(end).is_some_and(|rest| rest.starts_with(['/', ' ']));
        let lines = |c: usize| self.couplings[c].first..self.couplings[c].first + self.couplings[c].lines.len();
        let tightest = |c: usize| lines(c).max_by(|&a, &b| m.tensions[a].total_cmp(&m.tensions[b]));
        (0..self.couplings.len())
            .filter(|&c| self.couplings[c].ends.iter().any(|e| within(&e.path)))
            .filter_map(|c| tightest(c).map(|l| (m.tensions[l], l)))
            .max_by(|a, b| a.0.total_cmp(&b.0))
            .and_then(|(_, l)| m.mass_t.get(l).copied())
            .unwrap_or(0.0)
    }

    pub fn pull(&self, on: &str, m: &Moment) -> Option<Result<(f64, &str), String>> {
        let within = |end: &str| on == end || on.strip_prefix(end).is_some_and(|rest| rest.starts_with(['/', ' ']));
        let ours: Vec<usize> = (0..self.couplings.len()).filter(|&c| self.couplings[c].ends.iter().any(|e| within(&e.path))).collect();
        if ours.is_empty() {
            return None;
        }
        if let Some(why) = &m.coupling {
            return Some(Err(why.clone()));
        }
        let without: Vec<&str> = ours.iter().flat_map(|&c| self.blocks_coupling(c)).collect();
        if !without.is_empty() {
            return Some(Err(format!("these moving bodies its tension rests on have no mass: {}", self.massless_said(&without))));
        }
        let tightest = |c: usize| m.tensions[self.couplings[c].first..self.couplings[c].first + self.couplings[c].lines.len()].iter().fold(0.0f64, |p, x| p.max(*x));
        let c = ours.into_iter().max_by(|&a, &b| tightest(a).total_cmp(&tightest(b))).expect("a coupling");
        Some(Ok((tightest(c), self.couplings[c].name.as_str())))
    }

    /// Each coupling: its name (A~B (kind)) and its lines' place among a
    /// moment's tensions.
    pub fn couplings(&self) -> Vec<(&str, std::ops::Range<usize>)> {
        self.couplings.iter().map(|c| (c.name.as_str(), c.first..c.first + c.lines.len())).collect()
    }

    /// What the couplings with one end on what reaction r's path carries and the other off it exert on it at
    /// moment m (N, the world's axes): each line's tension along it, and their frictions, as the reaction
    /// counts them.
    pub fn coupling_pull(&self, r: usize, m: &Moment) -> Option<V3> {
        let mut pull = [0.0; 3];
        if self.couplings.iter().all(|c| c.side[r] == 0) && !m.frictions.iter().any(|fr| self.bears(fr, r)) {
            return Some(pull);
        }
        if m.coupling.is_some() || m.tensions.is_empty() {
            return None;
        }
        let ctx = self.ctx(m.t).ok()?;
        for c in self.couplings.iter().filter(|c| c.side[r] != 0) {
            for (i, (_, u, _, _)) in c.lines.iter().enumerate() {
                let u = u.each_ref().map(|e| e.at(&ctx));
                let t = m.tensions[c.first + i] / norm(u);
                (0..3).for_each(|k| pull[k] += f64::from(c.side[r]) * t * u[k]);
            }
        }
        for fr in m.frictions.iter().filter(|fr| self.bears(fr, r)) {
            (0..3).for_each(|k| pull[k] += fr.force[k]);
        }
        Some(pull)
    }

    /// Whether reaction r's path carries what friction fr acts on, and
    /// not its whole stage (whose friction is then its own).
    pub(super) fn bears(&self, fr: &Friction, r: usize) -> bool {
        let has = |&(c, e): &(usize, usize)| self.couplings[c].carries[r][e];
        fr.on.iter().all(has) && !fr.all.iter().all(has)
    }

    /// Whether reaction r may carry a coupling's friction whose opposite
    /// is on no part (Friction::paired): an end of a routed belt (one part
    /// of several couplings) giving an efficiency, carried without the
    /// whole belt.
    pub fn bears_unpaired(&self, r: usize) -> bool {
        self.couplings.iter().any(|c| {
            c.efficiency.is_some()
                && (0..2).any(|e| {
                    c.carries[r][e]
                        && self.couplings.iter().any(|d| !std::ptr::eq(c, d) && d.ends.iter().any(|x| x.path == c.ends[e].path))
                        && !self.couplings.iter().filter(|d| d.ends.iter().any(|x| x.path == c.ends[e].path)).all(|d| d.carries[r] == [true, true])
                })
        })
    }

    /// Whether any coupling states an efficiency: the actuators then supply more than the motion takes.
    pub fn lossy(&self) -> bool {
        self.couplings.iter().any(|c| c.efficiency.is_some())
    }

    /// The couplings passing power for actuator a at moment m that give no
    /// efficiency: its effort leaves their loss out (Model::refusal).
    pub fn unknown_losses(&self, a: usize, m: &Moment) -> Vec<&str> {
        let act = &self.actuators[a];
        m.lossless.iter().map(|&c| &self.couplings[c]).filter(|c| c.moves(&act.coupled)).map(|c| c.name.as_str()).collect()
    }

    /// Whether actuator a's couplings' loss rests on tensions that leave
    /// unseparated bodies out: a coupling with an efficiency passes power
    /// along its coordinates, and some body moves with the balances
    /// unseparated.
    pub fn loss_leaves_out(&self, a: usize) -> bool {
        let act = &self.actuators[a];
        !self.unseparated.is_empty() && self.couplings.iter().any(|c| c.efficiency.is_some() && c.moves(&act.coupled))
    }

    /// Whether a coupling reaches what reaction r's path carries (one end
    /// carried, the other not).
    pub fn coupled(&self, r: usize) -> bool {
        self.couplings.iter().any(|c| c.side[r] != 0)
    }

    /// The lines of the coupling relating instance paths a and b (either
    /// way), as their place among a moment's tensions.
    pub fn coupling_of(&self, a: &str, b: &str) -> Option<std::ops::Range<usize>> {
        let c = &self.couplings[self.coupling_index(a, b)?];
        Some(c.first..c.first + c.lines.len())
    }

    /// The index of the coupling relating instance paths a and b (either
    /// way).
    pub fn coupling_index(&self, a: &str, b: &str) -> Option<usize> {
        self.couplings.iter().position(|c| (c.ends[0].path == a && c.ends[1].path == b) || (c.ends[0].path == b && c.ends[1].path == a))
    }

    /// The force and moment (about the world's origin) the couplings exert
    /// on instance path `path` at moment m, where some has an end there:
    /// what a support that is one of their parts carries besides its share.
    pub(crate) fn coupled_on(&self, path: &str, m: &Moment) -> Option<(V3, V3)> {
        if m.tensions.is_empty() {
            return None;
        }
        let ctx = self.ctx(m.t).ok()?;
        let (mut force, mut moment, mut any) = ([0.0; 3], [0.0; 3], false);
        for c in &self.couplings {
            let s = f64::from(i8::from(c.ends[1].path == path) - i8::from(c.ends[0].path == path));
            if s == 0.0 {
                continue;
            }
            any = true;
            for (i, (p, u, h, _)) in c.lines.iter().enumerate() {
                let (p, u) = (p.each_ref().map(|e| e.at(&ctx)), u.each_ref().map(|e| e.at(&ctx)));
                let t = m.tensions[c.first + i] / norm(u);
                let f = u.map(|x| s * t * x);
                let mo = cross(p, f);
                for k in 0..3 {
                    force[k] += f[k];
                    moment[k] += mo[k] + h.at(&ctx) * f[k];
                }
            }
        }
        // Their frictions on it alone (a screw's thread friction on its nut).
        for fr in m.frictions.iter().filter(|fr| fr.on.iter().all(|&(c, e)| self.couplings[c].ends[e].path == path)) {
            any = true;
            for k in 0..3 {
                force[k] += fr.force[k];
                moment[k] += fr.moment[k];
            }
        }
        any.then_some((force, moment))
    }

    /// Why a carried belt separated where the model was built is not at
    /// moment ctx, if it is not: a piece of its path is off the span or
    /// wheel it lay on (Model::belts).
    pub(super) fn belts_off(&self, ctx: &Ctx) -> Option<String> {
        for s in self.strips.iter().filter(|s| !s.on.is_empty()) {
            for (i, rides) in s.on.iter().enumerate() {
                let [p0, p1] = s.ends[i].each_ref().map(|p| p.each_ref().map(|n| n.at(ctx)));
                let (c, l) = match *rides {
                    Rides::Span(c, l) | Rides::Wrap(c, l, _) => (c, l),
                };
                let cp = &self.couplings[c];
                let (q, u, _, len) = &cp.lines[l];
                let (q, u, len) = (q.each_ref().map(|n| n.at(ctx)), u.each_ref().map(|n| n.at(ctx)), len.at(ctx));
                let n = norm(u);
                let e: [V3; 2] = [std::array::from_fn(|k| q[k] + u[k] / n * len), q];
                let on = match *rides {
                    Rides::Span(..) => (same_point(p0, e[0]) && same_point(p1, e[1])) || (same_point(p0, e[1]) && same_point(p1, e[0])),
                    Rides::Wrap(_, _, k) => same_point(p0, e[k]) || same_point(p1, e[k]),
                };
                if !on {
                    return Some(format!(
                        "piece {} of the carried belt {} lies off {} here, which it lies on where the model is built: its material's push on the parts \
                         is not separated",
                        i + 1,
                        s.label,
                        cp.name
                    ));
                }
            }
        }
        None
    }

    /// Each coupling's ends, A and B, by instance path, in `couplings`'
    /// order.
    pub fn coupling_ends(&self) -> Vec<[&str; 2]> {
        self.couplings.iter().map(|c| [c.ends[0].path.as_str(), c.ends[1].path.as_str()]).collect()
    }

    /// Each coupling line at moment t, in the order of a moment's tensions:
    /// its points nearest A's and B's origins (where each is placed: the
    /// same line of action) and its unit direction in the world, along
    /// which its tension T pulls B (and A the other way), and its pitch h
    /// (mm: a screw's lead / 2 pi), its moment about the direction h T. Err
    /// where the solves it reads do not solve at t.
    pub fn lines_at(&self, t: f64) -> Result<Vec<([V3; 2], V3, f64)>, String> {
        let ctx = self.ctx(t)?;
        let mut out = vec![];
        for c in &self.couplings {
            let o: [V3; 2] = std::array::from_fn(|e| std::array::from_fn(|k| c.ends[e].on[9 + k].at(&ctx)));
            for (p, u, h, _) in &c.lines {
                let (p, u) = (p.each_ref().map(|e| e.at(&ctx)), u.each_ref().map(|e| e.at(&ctx)));
                let u = u.map(|x| x / norm(u));
                let nearest = |o: V3| -> V3 {
                    let s: f64 = (0..3).map(|k| (o[k] - p[k]) * u[k]).sum();
                    std::array::from_fn(|k| p[k] + s * u[k])
                };
                out.push(([nearest(o[0]), nearest(o[1])], u, h.at(&ctx)));
            }
        }
        Ok(out)
    }

    /// Every coupling's efficiency set to 1: its efforts the ideal
    /// mechanism's, as references of inertia alone state them.
    #[cfg(test)]
    pub(crate) fn without_losses(&mut self) {
        for c in &mut self.couplings {
            c.efficiency = Some((Num::C(1.0), Num::C(1.0)));
        }
        self.moments.lock().unwrap().clear();
        self.chosen.lock().unwrap().clear();
    }

    /// The balances giving the couplings' tensions: each one's reaction (by
    /// its place among `reactions`) and whether it is about a turn's axis.
    #[cfg(test)]
    pub(crate) fn balancing(&self) -> Vec<(usize, bool)> {
        self.balances.iter().map(|b| (b.reaction, b.turn)).collect()
    }

    /// Why coupling c's loss is not found though it gives an efficiency: an end of it rides a carrier
    /// that moves (`couple`).
    pub(super) fn loss_on_a_carrier(&self, c: usize) -> String {
        format!(
            "{} passes power between parts one of which rides a joint that moves (a planet on its carrier, a pulley on a carriage): its loss is on its parts' turning in \
             that carrier's frame, which is not found here",
            self.couplings[c].name
        )
    }

    /// Whether the coupling named `name` gives an efficiency (one in `unknown_losses` that does has its
    /// loss unfound for riding a moving carrier).
    pub fn gives_efficiency(&self, name: &str) -> bool {
        self.couplings.iter().any(|c| c.name == name && c.efficiency.is_some())
    }

    /// Whether coupling c's lines or frictions reach what reaction r
    /// carries: it carries an end of its train and not every one.
    pub(super) fn reaches(&self, c: usize, r: usize) -> bool {
        let ends = || (0..self.couplings.len()).filter(|&d| self.trains[d] == self.trains[c]).flat_map(|d| self.couplings[d].carries[r]);
        self.couplings[c].carries[r].iter().any(|&x| x) && !ends().all(|x| x)
    }

    /// The couplings at a moment (couple.rs). Each line must do no work along
    /// any coordinate, its two ends' points on it moving alike along it as
    /// the placements move them (else the coupling does not match the
    /// joints' motion: refused). Each line's tension comes from the balances
    /// of the passive joints they reach, and is then taken from the
    /// reactions carrying one end of it (what they get through it is the
    /// coupling's, not the rest of the model's). The loss of what passes
    /// power together, a mode (a belt's spans, a mesh's flanks) or else a
    /// coupling's lines, is friction (losses) on the part power enters it
    /// by, along the joints moving that part: a coupling's along each end's
    /// own joint (E at A, EB at B); a mode's along each part's, its own and
    /// those carrying it, but one moving every part. Each friction acts on
    /// its part as a wrench (with a coupling's, its pair on the other end:
    /// Friction), so the efforts, the reactions carrying it and the
    /// balances read it alike, the tensions solved again with it in them.
    /// Its tensions; why none are found; the couplings with no efficiency
    /// given passing power here.
    pub(super) fn couple(
        &self,
        ctx: &Ctx,
        reactions: &mut [Wrench],
        f: &mut [f64],
        scale: &mut [f64],
        frictions: &mut Vec<Friction>,
        tension_noise: &mut Vec<f64>,
        changes: &mut Vec<Vec<(usize, V3, V3)>>,
        mass_t: &mut Vec<f64>,
    ) -> (Vec<f64>, Option<String>, Vec<usize>, Vec<Site>) {
        if self.couplings.is_empty() {
            return (Vec::new(), None, Vec::new(), Vec::new());
        }
        let v3 = |n: &[Num; 3]| -> V3 { n.each_ref().map(|e| e.at(ctx)) };
        let vals = |n: &[Num]| -> Vec<f64> { n.iter().map(|e| e.at(ctx)).collect() };
        let lines: Vec<(V3, V3, f64, f64)> = self
            .couplings
            .iter()
            .flat_map(|c| &c.lines)
            .map(|(p, u, h, l)| {
                let u = v3(u);
                let n = norm(u);
                (v3(p), u.map(|x| x / n), h.at(ctx), l.at(ctx))
            })
            .collect();
        // A line's wrench on B, per N of tension, along a joint's motion
        // (axis a through o): its moment about it for a turn, its force
        // along it for a slide.
        let along = |(p, u, h, _): &(V3, V3, f64, f64), o: V3, a: V3, turn: bool| {
            if turn { dot(cross(std::array::from_fn(|k| p[k] - o[k]), *u), a) + h * dot(*u, a) } else { dot(*u, a) }
        };
        let column = |frame: &[f64], k: usize| -> V3 { [frame[k], frame[3 + k], frame[6 + k]] };
        let name = |l: usize| {
            let c = self.couplings.iter().rfind(|c| c.first <= l).unwrap();
            format!("line {} of {}", l - c.first + 1, c.name)
        };
        let mut why = None;
        for c in &self.couplings {
            if let Some(w) = c.unsound(ctx) {
                why = Some(w);
            }
        }
        // Each line does no work along any coordinate: along it, the point
        // of A on it and the point of B move alike, and they turn alike
        // about it times its pitch, to 1e-3 of the terms each end's motion
        // is made of (a point at rest for both ends, a planet's contact on
        // its sun, has terms that cancel; a radius rounded to a catalog's
        // digits is off by 1e-4; a crossed belt or a sign is off by its
        // whole speed).
        for c in &self.couplings {
            let ends: Vec<(M3, V3, Vec<(usize, M3, V3)>)> = c
                .ends
                .iter()
                .map(|e| {
                    let (r, o) = split(&vals(&e.on));
                    (r, o, e.dq.iter().map(|(k, d)| (*k, split(&vals(d)))).map(|(k, (dr, dp))| (k, dr, dp)).collect())
                })
                .collect();
            let mut coords: Vec<usize> = ends.iter().flat_map(|e| e.2.iter().map(|x| x.0)).collect();
            coords.sort_unstable();
            coords.dedup();
            for (i, l) in lines[c.first..c.first + c.lines.len()].iter().enumerate() {
                let (p, u, h, _) = *l;
                for &k in &coords {
                    // An end's velocity at P and turning rate along k, and
                    // the size of the terms they are made of.
                    let moving = |e: usize| -> (V3, V3, f64) {
                        let (r, o, dq) = &ends[e];
                        match dq.iter().find(|x| x.0 == k) {
                            None => ([0.0; 3], [0.0; 3], 0.0),
                            Some((_, dr, dp)) => {
                                let w = mul(dr, &tr(r));
                                let arm: V3 = std::array::from_fn(|j| p[j] - o[j]);
                                let turn = mv(&w, &arm);
                                (std::array::from_fn(|j| turn[j] + dp[j]), vee(&w), norm(turn) + norm(*dp))
                            }
                        }
                    };
                    let ((va, wa, sa), (vb, wb, sb)) = (moving(0), moving(1));
                    let work = dot(u, std::array::from_fn(|j| vb[j] - va[j])) + h * dot(u, std::array::from_fn(|j| wb[j] - wa[j]));
                    let size = sa + sb + h.abs() * (norm(wa) + norm(wb));
                    if work.abs() > 1e-3 * size {
                        // One part alone moves with this coordinate: no sign or ratio is wrong, the two are
                        // written as moving apart (a nut on a slide of its own beside the screw's angle).
                        let moves = |e: usize| ends[e].2.iter().any(|x| x.0 == k);
                        let cause = if !(moves(0) && moves(1)) {
                            "only one of its two parts moves with that coordinate. A coupling's parts move together by one coordinate: write the driven \
                             part's motion as a function of the driver's (a nut on a screw turned by `angle`: (slide z (* (/ angle 360) lead))), not \
                             as a coordinate of its own passed in step"
                        } else {
                            "the coupling does not match the joints' motion (a crossed belt written as an open one, a ratio's, a turn's or a lead's sign)"
                        };
                        why = Some(format!(
                            "line {} of {} separates its parts at {work:.6} mm per unit of {}: {cause}, so it would do work the motion does not",
                            i + 1,
                            c.name,
                            self.coords[k].name
                        ));
                    }
                }
            }
        }
        // Each balancing joint's balance along its motion, from what the
        // reactions hold before the lines are taken out.
        let balances = |reactions: &[Wrench]| -> Vec<crate::couple::Balance> {
            self.balances
                .iter()
                .map(|b| {
                    let w = &reactions[b.reaction];
                    let (a, o) = (column(&w.frame, b.axis), [w.frame[9], w.frame[10], w.frame[11]]);
                    let need = if b.turn { dot(std::array::from_fn(|k| w.moment[k] - cross(o, w.force)[k]), a) } else { dot(w.force, a) };
                    let noise = ROUNDING * if b.turn { w.noise.1 + norm(o) * w.noise.0 } else { w.noise.0 };
                    let mut row = vec![0.0; lines.len()];
                    for c in &self.couplings {
                        for i in 0..c.lines.len() {
                            row[c.first + i] = f64::from(c.side[b.reaction]) * along(&lines[c.first + i], o, a, b.turn);
                        }
                    }
                    crate::couple::Balance { joint: self.reactions[b.reaction].path.clone(), row, need, noise }
                })
                .collect()
        };
        let eqs = balances(reactions);
        // Each line's rounding: each balance's in its own units over the line's coefficient in it.
        *tension_noise = (0..lines.len()).map(|l| eqs.iter().filter(|q| q.row[l] != 0.0).map(|q| q.noise / q.row[l].abs()).sum()).collect();
        // What the balances read is unfound where a reaction they rest on
        // is (a body with no mass, one whose shape moves).
        if let Some(b) = self.balances.iter().find(|b| self.static_refusal(b.reaction).is_some()) {
            why = Some(format!(
                "the balance of {} gives its couplings' tensions, and its reaction is not found: {}",
                self.reactions[b.reaction].path,
                self.static_refusal(b.reaction).unwrap()
            ));
        }
        let read: Vec<crate::couple::Line> = self
            .couplings
            .iter()
            .flat_map(|c| {
                let t0 = c.pretension.at(ctx);
                lines[c.first..c.first + c.lines.len()].iter().map(move |l| crate::couple::Line { t0, length: l.3, held: c.held })
            })
            .collect();
        let (tensions, modes) = match why {
            Some(why) => return (Vec::new(), Some(why), Vec::new(), Vec::new()),
            None => match crate::couple::tensions(&eqs, &read, &name) {
                Ok(t) => t,
                Err(why) => return (Vec::new(), Some(why), Vec::new(), Vec::new()),
            },
        };
        // What each measured moment's change moves the tensions by (couple::tension_moves), and so the
        // reactions at the lines' ends: a line's linear move with the change itself, each set's lift as a
        // change of its own (its sign not known).
        {
            let dneeds: Vec<Vec<f64>> = changes
                .iter()
                .map(|ch| {
                    self.balances
                        .iter()
                        .map(|b| {
                            let w = &reactions[b.reaction];
                            let (a, o) = (column(&w.frame, b.axis), [w.frame[9], w.frame[10], w.frame[11]]);
                            ch.iter()
                                .filter(|c| c.0 == b.reaction)
                                .map(|(_, df, dm)| if b.turn { dot(std::array::from_fn(|k| dm[k] - cross(o, *df)[k]), a) } else { dot(*df, a) })
                                .sum()
                        })
                        .collect()
                })
                .collect();
            let moves = crate::couple::tension_moves(&eqs, &read, &dneeds);
            *mass_t = vec![0.0; lines.len()];
            // A line's pull on the reactions at its ends, per N of tension.
            let pull = |l: usize| -> (V3, V3) {
                let (p, u, h, _) = lines[l];
                (u, std::array::from_fn(|k| cross(p, u)[k] + h * u[k]))
            };
            let mut lifted: Vec<Vec<(usize, V3, V3)>> = Vec::new();
            for (ch, (dx, lifts)) in changes.iter_mut().zip(&moves) {
                for c in &self.couplings {
                    for l in c.first..c.first + c.lines.len() {
                        mass_t[l] += dx[l].abs();
                        let (u, m) = pull(l);
                        for r in (0..reactions.len()).filter(|&r| c.side[r] != 0) {
                            let s = f64::from(c.side[r]) * dx[l];
                            ch.push((r, u.map(|x| -s * x), m.map(|x| -s * x)));
                        }
                    }
                }
                for (lift, ls) in lifts.iter().filter(|x| x.0 > 0.0) {
                    let mut moved = Vec::new();
                    for &(l, c) in ls {
                        mass_t[l] += lift * c.abs();
                        let owner = self.couplings.iter().find(|cp| (cp.first..cp.first + cp.lines.len()).contains(&l)).expect("a line's coupling");
                        let (u, m) = pull(l);
                        for r in (0..reactions.len()).filter(|&r| owner.side[r] != 0) {
                            let s = f64::from(owner.side[r]) * lift * c;
                            moved.push((r, u.map(|x| -s * x), m.map(|x| -s * x)));
                        }
                    }
                    lifted.push(moved);
                }
            }
            changes.extend(lifted);
        }
        let owner: Vec<usize> = self.couplings.iter().enumerate().flat_map(|(i, c)| std::iter::repeat_n(i, c.lines.len())).collect();
        // A joint's frame, axis, moment unit and rate (rad or mm a second).
        let joint = |inp: &Input| {
            let frame = vals(&inp.frame);
            let unit = if inp.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
            (column(&frame, inp.axis), [frame[9], frame[10], frame[11]], unit, inp.rate.at(ctx) * unit)
        };
        // Each place a loss may be charged, with the coefficient of each
        // line's tension in its held moment, and the balance its joint
        // gives, if it gives one.
        let mut sites: Vec<Site> = Vec::new();
        let mut reads: Vec<(Vec<(usize, f64)>, Option<usize>, usize, bool)> = Vec::new();
        // What each site's friction acts on (Acts).
        let mut acts: Vec<Acts> = Vec::new();
        let mut site = |inp: &Input, coef: Vec<(usize, f64)>, eff: f64, worst: f64, stage: usize, solo: bool, act: Acts| {
            acts.push(act);
            let unit = if inp.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
            let mut weights = vec![0.0; self.coords.len()];
            for (k, dk) in &inp.along {
                weights[*k] = unit * dk.at(ctx);
            }
            let held = coef.iter().map(|&(l, c)| c * tensions[l]).sum();
            let held_roundoff = f64::EPSILON * coef.len() as f64 * coef.iter().map(|&(l, c)| (c * tensions[l]).abs()).sum::<f64>();
            let balance = self.balances.iter().position(|b| self.reactions[b.reaction].path == inp.joint);
            reads.push((coef, balance, stage, solo));
            sites.push(Site { held, held_roundoff, rate: joint(inp).3, eff, worst, weights, gain: Vec::new(), total: held, shared: false, unfed: 0.0, stage });
        };
        // Each joint moving the ends of couplings `cs` (an end's own, then
        // those carrying it), with the ends it moves.
        let moving = |cs: &[usize]| -> Vec<(&Input, Vec<(usize, usize)>)> {
            let mut out: Vec<(&Input, Vec<(usize, usize)>)> = Vec::new();
            for &c in cs {
                for e in 0..2 {
                    for inp in self.couplings[c].ends[e].chain() {
                        match out.iter_mut().find(|(i, _)| i.joint == inp.joint) {
                            Some((_, es)) => es.push((c, e)),
                            None => out.push((inp, vec![(c, e)])),
                        }
                    }
                }
            }
            out
        };
        // What passes power together: a mode, else a coupling's lines.
        let mut groups: Vec<((bool, usize), Vec<usize>)> = Vec::new();
        for l in 0..lines.len() {
            let key = modes[l].map_or((false, owner[l]), |m| (true, m));
            match groups.iter_mut().find(|g| g.0 == key) {
                Some(g) => g.1.push(l),
                None => groups.push((key, vec![l])),
            }
        }
        let mut lossless = Vec::new();
        for (stage, (_, group)) in groups.iter().enumerate() {
            let mut cs: Vec<usize> = group.iter().map(|&l| owner[l]).collect();
            cs.dedup();
            if let [ci] = cs.as_slice() {
                // Each end that rides a joint: A, whose lines pull it the
                // other way, then B.
                // Each end that rides a joint: A, whose lines pull it the
                // other way, then B. Its loss is on the power through that
                // joint alone, the frame its axle rides (Willis's: a
                // planet's mesh loses on its turning on the carrier, not on
                // the carrier's), and its friction pairs with the other end,
                // which it drives.
                let c = &self.couplings[*ci];
                let eff = c.efficiency.as_ref().map(|(e, b)| (e.at(ctx), b.at(ctx)));
                // An end whose axle is carried by a joint that moves (a planet on its carrier, a pulley
                // on a carriage): the mesh's power is its ends' turning in that carrier's frame, which no
                // end's own joint gives (a sun held still has none; one turning slower than the carrier
                // gives it the other way). Its loss is not found there, and is said so, as one not given.
                let own = |e: usize| c.ends[e].input.as_ref().map(|i| i.joint.as_str());
                let rides = (0..2).any(|e| c.ends[e].chain().any(|inp| Some(inp.joint.as_str()) != own(e) && joint(inp).3 != 0.0));
                // Willis's frame: where a turning carrier carries an end's axle (a planet's), the coupling's
                // power is its ends' turning in the carrier's frame. An end the carrier carries turns there
                // by its own joint; another, coaxial with the carrier, by its own turn less the carrier's
                // (a sun held still turns back at the carrier's rate). Charged there, its loss is found.
                let carrier: Option<&Input> = (0..2).flat_map(|e| c.ends[e].chain().filter(move |inp| Some(inp.joint.as_str()) != own(e))).find(|inp| joint(inp).3 != 0.0);
                let carried = |e: usize| carrier.is_some_and(|cr| c.ends[e].chain().any(|inp| inp.joint == cr.joint));
                let relative: Option<Vec<Option<Input>>> = carrier.filter(|cr| cr.turn).and_then(|cr| {
                    let (ca, co, _, _) = joint(cr);
                    (0..2)
                        .map(|e| {
                            if carried(e) {
                                return Some(None);
                            }
                            let mine = c.ends[e].input.as_ref();
                            // Coaxial with the carrier, or held (no joint: the mesh's centre on its axis).
                            if let Some(m) = mine {
                                let (ma, mo, _, _) = joint(m);
                                let off: V3 = std::array::from_fn(|k| mo[k] - co[k]);
                                let across = norm(std::array::from_fn(|k| off[k] - dot(off, ca) * ca[k]));
                                if !m.turn || dot(ma, ca).abs() < 1.0 - 1e-9 || across > 1e-6 * (1.0 + norm(off)) {
                                    return None;
                                }
                            }
                            let sign = mine.map_or(1.0, |m| dot(joint(m).0, ca).signum());
                            let mut along: Vec<(usize, f64)> = cr.along.iter().map(|(k, d)| (*k, -d.at(ctx))).collect();
                            let mut rate = -cr.rate.at(ctx);
                            if let Some(m) = mine {
                                rate += sign * m.rate.at(ctx);
                                for (k, d) in &m.along {
                                    match along.iter_mut().find(|x| x.0 == *k) {
                                        Some(x) => x.1 += sign * d.at(ctx),
                                        None => along.push((*k, sign * d.at(ctx))),
                                    }
                                }
                            }
                            let frame: [Num; 12] = vals(&cr.frame).try_into().map(|f: [f64; 12]| f.map(Num::C)).ok()?;
                            Some(Some(Input {
                                joint: mine.map_or_else(|| format!("{} (in its carrier)", c.ends[e].path), |m| m.joint.clone()),
                                frame,
                                axis: cr.axis,
                                turn: true,
                                rate: Num::C(rate),
                                along: along.into_iter().map(|(k, d)| (k, Num::C(d))).collect(),
                            }))
                        })
                        .collect()
                });
                let rides = rides && relative.is_none();
                // (One that loses nothing has no loss to find.)
                let unfound = rides && eff.is_some_and(|(e, b)| e != 1.0 || b != 1.0);
                // Ends on one line (two slides along it, two wheels on one axle) move as one along it: the
                // coupling passes its power between them with nothing turning or bending in it, and loses
                // none, whatever it states (a cord between two carts on a rail). Nor could its friction be
                // paired: a force along that line on the other does work there.
                let joints: Vec<Option<(bool, V3, V3)>> = (0..2)
                    .map(|e| {
                        c.ends[e].input.as_ref().map(|i| {
                            let (a, o, _, _) = joint(i);
                            (i.turn, a, o)
                        })
                    })
                    .collect();
                let one_line = matches!((joints[0], joints[1]), (Some(p), q @ Some(_)) if pair(p, q, 1.0, None).is_none());
                let eff = eff.filter(|_| !unfound && !one_line);
                let mut passes = unfound && group.iter().any(|&l| tensions[l] != 0.0);
                for (end, sign) in [(0, -1.0), (1, 1.0)] {
                    let Some(inp) = relative.as_ref().and_then(|r| r[end].as_ref()).or(c.ends[end].input.as_ref()) else { continue };
                    let (a, o, _, rate) = joint(inp);
                    let coef: Vec<(usize, f64)> = group.iter().map(|&l| (l, sign * along(&lines[l], o, a, inp.turn))).collect();
                    let held: f64 = coef.iter().map(|&(l, c)| c * tensions[l]).sum();
                    passes |= held * rate != 0.0;
                    if let Some((e, b)) = eff {
                        let (eff, worst) = if end == 0 { (e, 1.0 / e - 1.0) } else { (b, 1.0 - b) };
                        let act =
                            Acts { joint: (inp.turn, a, o), parts: vec![(vec![(*ci, end)], coef.clone(), None)], all: vec![(*ci, 0), (*ci, 1)], partner: Some((*ci, 1 - end)) };
                        site(inp, coef, eff, worst, stage, true, act);
                    }
                }
                if eff.is_none() && passes && !one_line {
                    lossless.push(*ci);
                }
                continue;
            }
            // Its parts, each with the moment the mode's lines exert on it
            // along its joint, and each line's coefficient in it.
            // Each joint moving its parts (their own, and those carrying
            // them: a gantry's slide pulling its idlers' axles), with the
            // moment the mode's lines exert along it on the parts it moves,
            // and each line's coefficient in it. One moving every part
            // passes nothing.
            cs.sort_unstable();
            cs.dedup();
            let all: Vec<(usize, usize)> = cs.iter().flat_map(|&c| [(c, 0), (c, 1)]).collect();
            // Each part it moves: its ends, its lines' coefficients, and
            // where a friction along a joint carrying its own acts on it:
            // at its axle, along the way that joint moves it.
            let mut parts: Vec<(&Input, f64, Vec<(usize, f64)>, Vec<PartActs>)> = Vec::new();
            for (inp, ends) in moving(&cs) {
                let (a, o, _, _) = joint(inp);
                let paths: Vec<&str> = ends.iter().map(|&(c, e)| self.couplings[c].ends[e].path.as_str()).collect();
                if all.iter().all(|&(c, e)| paths.contains(&self.couplings[c].ends[e].path.as_str())) {
                    continue;
                }
                let mut each: Vec<PartActs> = Vec::new();
                for &(c, e) in &ends {
                    let end = &self.couplings[c].ends[e];
                    let coef: Vec<(usize, f64)> =
                        group.iter().filter(|&&l| owner[l] == c).map(|&l| (l, if e == 0 { -1.0 } else { 1.0 } * along(&lines[l], o, a, inp.turn))).collect();
                    let axle = (end.input.as_ref().is_none_or(|i| i.joint != inp.joint)).then(|| match &end.input {
                        Some(i) => joint(i).1,
                        None => std::array::from_fn(|k| end.on[9 + k].at(ctx)),
                    });
                    match each.iter_mut().find(|p| self.couplings[p.0[0].0].ends[p.0[0].1].path == end.path) {
                        Some(p) => {
                            p.0.push((c, e));
                            p.1.extend(coef);
                        }
                        None => each.push((vec![(c, e)], coef, axle)),
                    }
                }
                let coef: Vec<(usize, f64)> = each.iter().flat_map(|p| p.1.iter().copied()).collect();
                parts.push((inp, coef.iter().map(|&(l, c)| c * tensions[l]).sum(), coef, each));
            }
            let given = parts.iter().any(|p| p.1 * joint(p.0).3 < 0.0);
            let effs: Vec<Option<(f64, f64)>> = cs.iter().map(|&c| self.couplings[c].efficiency.as_ref().map(|(e, b)| (e.at(ctx), b.at(ctx)))).collect();
            if let Some(Some((e, _))) = effs.first().filter(|_| effs.iter().all(Option::is_some)) {
                for p in &parts {
                    site(
                        p.0,
                        p.2.clone(),
                        *e,
                        1.0 / e - 1.0,
                        stage,
                        false,
                        Acts { joint: (p.0.turn, joint(p.0).0, joint(p.0).1), parts: p.3.clone(), all: all.clone(), partner: None },
                    );
                }
            }
            if !given {
                continue;
            }
            if effs.iter().any(Option::is_none) {
                lossless.extend(&cs);
                continue;
            }
            let e = effs[0].unwrap().0;
            if effs.iter().flatten().any(|&(x, b)| x != e || b != e) {
                let named: Vec<&str> = cs.iter().map(|&c| self.couplings[c].name.as_str()).collect();
                return (tensions, Some(format!("{} tighten together as one belt and give different efficiencies: give each its drive's one", named.join(", "))), lossless, sites);
            }
        }
        // Friction x in balance b's joint leaves its lines need - x, fed to
        // the stage driving that joint where one is (an intermediate shaft:
        // its balance reads one other stage's lines), and to the coupling's
        // own lines at that end (a passive end driving back). A stage's loss
        // is charged once, on what enters it: its ends and parts gain
        // nothing from each other's (a self-locking screw lowered is driven
        // by its motor through its EB alone). A coupling feeds; a mode of
        // several (a routed belt) does not, where it meets another stage.
        // Where several other stages meet in the joint (a head two belts
        // drive) the friction is not fed: `unfed`.
        let mut stage_of = vec![0; lines.len()];
        for (g, (_, ls)) in groups.iter().enumerate() {
            for &l in ls {
                stage_of[l] = g;
            }
        }
        // Whether each site's friction is fed (None: nothing to feed).
        let feeds: Vec<Option<bool>> = reads
            .iter()
            .map(|(_, b, g, solo)| {
                let b = (*b)?;
                let mut others: Vec<usize> = (0..lines.len()).filter(|&l| eqs[b].row[l] != 0.0 && stage_of[l] != *g).map(|l| stage_of[l]).collect();
                others.sort_unstable();
                others.dedup();
                // A mode's part (a routed belt's) meeting another stage (a
                // head two belts drive) is not fed: it is not a stage in
                // series but beside one. Alone in its joint (an idler, a
                // carriage the belt alone pulls), its friction is in its
                // own balance.
                Some(if *solo { others.len() <= 1 } else { others.is_empty() })
            })
            .collect();
        // Each site's gain: the moment its lines exert per unit of friction
        // in a balancing joint, from the tensions with that balance's need
        // moved. They are linear in it while no line slackens or tightens:
        // a step that changes which lines are slack is shrunk until none
        // does. Only where a stage's friction is fed: the others' tensions
        // are not solved again.
        let mut slopes: Vec<Option<Vec<f64>>> = vec![None; eqs.len()];
        for b in reads.iter().zip(&feeds).filter(|(_, f)| **f == Some(true)).filter_map(|(r, _)| r.1) {
            if slopes[b].is_some() {
                continue;
            }
            let mut step = 1e-3 * eqs.iter().map(|e| e.need.abs()).fold(1.0, f64::max);
            for _ in 0..6 {
                let moved: Vec<crate::couple::Balance> = eqs
                    .iter()
                    .enumerate()
                    .map(|(i, e)| crate::couple::Balance { joint: e.joint.clone(), row: e.row.clone(), need: e.need + if i == b { step } else { 0.0 }, noise: e.noise })
                    .collect();
                let t = match crate::couple::tensions(&moved, &read, &name) {
                    Ok((t, _)) => t,
                    Err(why) => return (tensions, Some(why), lossless, sites),
                };
                slopes[b] = Some(t.iter().zip(&tensions).map(|(a, b)| (a - b) / step).collect());
                // Slack: within rounding of 0.
                let tiny = 1e-9 * tensions.iter().chain(&t).fold(1.0, |m: f64, x| m.max(x.abs()));
                if t.iter().zip(&tensions).all(|(a, b)| (a.abs() <= tiny) == (b.abs() <= tiny)) {
                    break;
                }
                step *= 1e-2;
            }
        }
        let gains: Vec<Vec<(usize, f64)>> = reads
            .iter()
            .enumerate()
            .map(|(j, (coef, _, g, _))| {
                reads
                    .iter()
                    .zip(&feeds)
                    .enumerate()
                    .filter_map(|(k, ((_, b, h, _), fed))| match b.and_then(|b| slopes[b].as_ref()) {
                        _ if *fed != Some(true) || (h == g && k != j) => None,
                        Some(dt) => Some((k, -coef.iter().map(|&(l, c)| c * dt[l]).sum::<f64>())).filter(|x| x.1 != 0.0),
                        None => None,
                    })
                    .collect()
            })
            .collect();
        for (s, g) in sites.iter_mut().zip(gains) {
            s.gain = g;
        }
        let held: Vec<f64> = sites.iter().map(|s| s.held).collect();
        let x = match losses(&sites, &held) {
            Ok((x, total)) => {
                for (((s, x), total), fed) in sites.iter_mut().zip(&x).zip(total).zip(&feeds) {
                    s.total = total;
                    s.shared = *fed == Some(false);
                    s.unfed = if s.shared { *x } else { 0.0 };
                    for (k, w) in s.weights.iter().enumerate() {
                        f[k] -= x * w;
                        scale[k] += (x * w).abs();
                    }
                }
                x
            }
            Err(why) => return (tensions, Some(why), lossless, sites),
        };
        // Each friction acts on the parts it is charged at, along the joint
        // power enters by (a moment about a turn's axis, a force along a
        // slide's through its origin), as the coupling's own: in the
        // reactions carrying them, and in a balancing joint's balance, so
        // the stages driving it pass it. At an end's own joint it pairs with
        // the other end (Newton's third law: pair).
        for (j, a) in acts.iter().enumerate().filter(|(j, _)| x[*j] != 0.0) {
            let (turn, axis, o) = a.joint;
            let other = a.partner.and_then(|(c, e)| self.couplings[c].ends[e].input.as_ref()).map(|i| {
                let (b, p, _, _) = joint(i);
                (i.turn, b, p)
            });
            // A coupling whose lines have a pitch (a screw): its friction is
            // its thread's, a moment about its axis on both parts whichever
            // end power enters by (Rexroth's drive torque at the nut's load
            // either way), M U on B and -M U on A, doing the work x does
            // along the site's joint: M (wB - wA).U = x w along each
            // coordinate.
            if let Some((c, _)) = a.partner {
                let cp = &self.couplings[c];
                if let Some(l) = (0..cp.lines.len()).find(|&i| lines[cp.first + i].2 != 0.0) {
                    let u = lines[cp.first + l].1;
                    let turning = |e: usize| -> Vec<(usize, V3)> {
                        let r = rot(&vals(&cp.ends[e].on)[..9]);
                        cp.ends[e].dq.iter().map(|(k, d)| (*k, vee(&mul(&rot(&vals(d)[..9]), &tr(&r))))).collect()
                    };
                    let (wa, wb) = (turning(0), turning(1));
                    let rel = |k: usize| -> f64 {
                        let w = |ws: &[(usize, V3)]| ws.iter().find(|x| x.0 == k).map_or([0.0; 3], |x| x.1);
                        dot(std::array::from_fn(|i| w(&wb)[i] - w(&wa)[i]), u)
                    };
                    let k = (0..sites[j].weights.len()).max_by(|&p, &q| rel(p).abs().total_cmp(&rel(q).abs()));
                    if let Some(k) = k.filter(|&k| rel(k).abs() > 1e-12) {
                        let m = x[j] * sites[j].weights[k] / rel(k);
                        let moment = u.map(|v| m * v);
                        frictions.push(Friction { on: vec![(c, 1)], all: a.all.clone(), force: [0.0; 3], moment, paired: true });
                        frictions.push(Friction { on: vec![(c, 0)], all: a.all.clone(), force: [0.0; 3], moment: moment.map(|v| -v), paired: true });
                        continue;
                    }
                }
            }
            // Where the pair's force acts along Q's axis: level with the
            // coupling's lines (a belt's plane, a mesh's face).
            let level = a.partner.map(|(c, _)| lines[self.couplings[c].first].0);
            if let Some((q, (force, moment))) = a.partner.zip(pair((turn, axis, o), other, x[j], level)) {
                frictions.push(Friction { on: a.parts[0].0.clone(), all: a.all.clone(), force: force.map(|v| -v), moment: moment.map(|v| -v), paired: true });
                frictions.push(Friction { on: vec![q], all: a.all.clone(), force, moment, paired: true });
                continue;
            }
            // Unpaired (a routed belt's parts on a joint carrying them): each part its share, by the size
            // of what its lines pull it with along the joint. The shares sum to 1 and none passes the
            // friction: by their signed moments over the joint's net, where the parts pull against each
            // other (a gantry's idlers, each by its pretension) and the net is small, each took the loss
            // of its own pull, the gantry's idlers 3.5 N either way at 0.95 and 15 lbf, appearing and
            // going as the gantry's power changed sign.
            let pulls: Vec<f64> = a.parts.iter().map(|(_, coef, _)| coef.iter().map(|&(l, c)| c * tensions[l]).sum::<f64>().abs()).collect();
            let all: f64 = pulls.iter().sum();
            for ((on, _, axle), pull) in a.parts.iter().zip(&pulls) {
                let share = if all > 0.0 { pull / all } else { 1.0 / a.parts.len() as f64 };
                let xi = x[j] * share;
                if xi == 0.0 {
                    continue;
                }
                let (force, moment) = match axle {
                    // Along a joint carrying its own: a force at its axle,
                    // F . v = xi, v how that joint moves the axle.
                    Some(p) => {
                        let v = if turn { cross(axis, std::array::from_fn(|k| p[k] - o[k])) } else { axis };
                        let vv = dot(v, v);
                        if vv <= 0.0 {
                            continue;
                        }
                        let f = v.map(|c| xi * c / vv);
                        (f, cross(*p, f))
                    }
                    None if turn => ([0.0; 3], axis.map(|v| xi * v)),
                    None => {
                        let f = axis.map(|v| xi * v);
                        (f, cross(o, f))
                    }
                };
                frictions.push(Friction { on: on.clone(), all: a.all.clone(), force, moment, paired: false });
            }
        }
        for fr in frictions.iter() {
            for (r, w) in reactions.iter_mut().enumerate() {
                if self.bears(fr, r) {
                    for k in 0..3 {
                        w.force[k] -= fr.force[k];
                        w.moment[k] -= fr.moment[k];
                    }
                    w.noise.0 += norm(fr.force);
                    w.noise.1 += norm(fr.moment);
                }
            }
        }
        // The tensions with the frictions in the balances.
        let tensions = if frictions.is_empty() {
            tensions
        } else {
            match crate::couple::tensions(&balances(reactions), &read, &name) {
                Ok((t, _)) => t,
                Err(why) => return (tensions, Some(format!("with the couplings' losses in the balances: {why}")), lossless, sites),
            }
        };
        // Each line taken from the reactions carrying one end of it: what
        // they get through it is the coupling's, not the rest of the
        // model's.
        for c in &self.couplings {
            for (l, t) in lines[c.first..c.first + c.lines.len()].iter().zip(&tensions[c.first..c.first + c.lines.len()]) {
                let (p, u, h, _) = *l;
                let force = u.map(|x| t * x);
                let moment: V3 = std::array::from_fn(|k| cross(p, force)[k] + h * force[k]);
                for (r, w) in reactions.iter_mut().enumerate() {
                    let s = f64::from(c.side[r]);
                    if s != 0.0 {
                        for k in 0..3 {
                            w.force[k] -= s * force[k];
                            w.moment[k] -= s * moment[k];
                        }
                        w.noise.0 += norm(force);
                        w.noise.1 += norm(moment);
                    }
                }
            }
        }
        (tensions, None, lossless, sites)
    }
}
