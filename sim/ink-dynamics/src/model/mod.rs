//! Inverse dynamics by projection onto the mechanism's coordinates (Kane's
//! method, the virtual-power form of d'Alembert's principle). Every body's
//! placement is an expression of t and of the coordinates (the arguments
//! the target's instances pass, built as value + a param at 0:
//! Doc::with_coordinates), and of what solves set from them (solved again
//! at each moment), differentiated exactly: along t for each body's
//! acceleration and angular velocity and acceleration, along each
//! coordinate for its partial velocities. The effort along coordinate k is
//! f_k = sum over bodies of Jv_k . m (a - g) + Jw_k . (I alpha + omega x I
//! omega); the actuators supply it, B^T tau = f, B the rate of each
//! actuator's drive along each coordinate (radians for a turn), split by
//! least norm weighted by each one's span of limits. Loads (`(load ...)`)
//! enter through the partial velocity of the point they act at. A path's
//! reaction, the wrench the rest of the model exerts on what it carries, is
//! the same bodies' m (a - g) and I alpha + omega x I omega summed (Newton
//! and Euler on them), less the loads and couplings' lines on them; the
//! lines' tensions come from the balance of the passive joints they reach
//! (couple.rs), and a coupling's loss is friction it exerts on its parts,
//! in the efforts, the reactions and the balances alike (Model::couple).

mod bodies;
mod bounded;
mod build;
mod couplings;
mod drives;
mod inverse;
mod kinematics;
mod refusal;
mod sampling;
mod spin;

/// Natural frequencies at a pose: a child module, reading the model's parts.
#[path = "../modes.rs"]
pub mod modes;

/// Springs in a simulation: the modes' model stepped forward.
#[path = "../springs.rs"]
pub mod springs;

#[cfg(test)]
mod tests;

pub(crate) use self::bounded::Plan;
pub use self::build::riding;
pub(crate) use self::build::{carried, joint_axis};
pub use self::couplings::{Friction, Site, losses, unfed_need};
pub use self::drives::shares;
pub(crate) use self::drives::{eigen, pinv};
#[cfg(test)]
pub(crate) use self::sampling::CHOSEN;
pub use self::sampling::{OFF_LINE, REACTION_STEPS, STEPS, Serve};

use self::build::{path_frame, rides};
use self::drives::dot_n;
use self::kinematics::nums;
use self::sampling::{Chosen, Tracked};
use ink_core::eval::{Ctx, eval_point};
use ink_core::graph::{Graph, NodeId, Var};
use ink_core::interval::I;
use ink_core::tape::Tape;
use ink_core::v3::{cross, dot, norm};
use ink_lang::view::ModelView;
use ink_lang::{Doc, InkError};
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

/// A number of the model: a constant, or a field of params and t; one
/// reading what a solve on t sets (a placement through a loop's closure)
/// with its rate, exact through the solve (Doc::derivative), since a
/// solved value has no bound over a step.
#[derive(Clone)]
pub(crate) enum Num {
    C(f64),
    T(Arc<Tape>),
    S(Arc<Tape>, Arc<Tape>),
}

impl Num {
    pub(crate) fn of(g: &Graph, n: NodeId) -> Num {
        match g.constant_value(n) {
            Some(c) => Num::C(c),
            None => Num::T(Arc::new(Tape::new(g, n))),
        }
    }

    /// `of`, or S where n reads what a solve on t sets.
    fn timed(d: &mut Doc, n: NodeId) -> Num {
        if !d.solves_read(&[n]).iter().any(|&i| d.solves[i].on_t) {
            return Num::of(&d.g, n);
        }
        let r = d.derivative(n, Var::T);
        Num::S(Arc::new(Tape::new(&d.g, n)), Arc::new(Tape::new(&d.g, r)))
    }

    pub(crate) fn at(&self, ctx: &Ctx) -> f64 {
        match self {
            Num::C(c) => *c,
            Num::T(t) | Num::S(t, _) => eval_point(t, ctx, [0.0; 3]),
        }
    }

    /// Its values at points `pts` (along a carried body's path).
    fn at_points(&self, ctx: &Ctx, pts: &[[f64; 3]], out: &mut Vec<f64>) {
        match self {
            Num::C(c) => {
                out.clear();
                out.resize(pts.len(), *c);
            }
            Num::T(t) | Num::S(t, _) => ink_core::eval::eval_points(t, ctx, pts, out),
        }
    }

    /// Its bounds over ctx's scope (a step: its t_span); through a solve,
    /// none.
    fn over(&self, ctx: &Ctx) -> I {
        match self {
            Num::C(c) => I::point(*c),
            Num::T(t) => ink_core::eval::eval_interval(t, ctx, [I::point(0.0); 3], None),
            Num::S(..) => I::ENTIRE,
        }
    }
}

type M3 = [[f64; 3]; 3];
type V3 = [f64; 3];

fn mul(a: &M3, b: &M3) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}

fn tr(a: &M3) -> M3 {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}

fn mv(a: &M3, v: &V3) -> V3 {
    std::array::from_fn(|i| (0..3).map(|k| a[i][k] * v[k]).sum())
}

/// The axial vector of a matrix's skew part.
fn vee(w: &M3) -> V3 {
    [0.5 * (w[2][1] - w[1][2]), 0.5 * (w[0][2] - w[2][0]), 0.5 * (w[1][0] - w[0][1])]
}

/// A rotation's 9 numbers (row by row) as R.
fn rot(e: &[f64]) -> M3 {
    [[e[0], e[1], e[2]], [e[3], e[4], e[5]], [e[6], e[7], e[8]]]
}

/// A placement's 12 numbers (R row by row, then its translation) as R and t.
fn split(e: &[f64]) -> (M3, V3) {
    (rot(e), [e[9], e[10], e[11]])
}

/// What an effort reads of a placement: its rotation, the rotation's rate
/// and its translation, its numbers' second derivatives along t, and their
/// derivatives along each coordinate that moves it.
struct Placement {
    r: [Num; 9],
    dr: [Num; 9],
    p: [Num; 3],
    dp: [Num; 3],
    dt2: [Num; 12],
    dq: Vec<(usize, [Num; 12])>,
}

/// A rigid body that moves along some coordinate, or that a reaction asked
/// for carries: its mass properties in its own frame, its placement, and
/// the reactions it is of (by index).
struct Body {
    /// Its instance path.
    label: String,
    mass: f64,
    /// Centre of mass, in the body's frame (mm).
    center: V3,
    /// About the centre of mass, in the body's frame (g mm^2).
    inertia: M3,
    at: Placement,
    /// Its group (Group), where it is a member, and the axis in its frame
    /// (a point on it, its direction): read round about it where the group
    /// is read with its spin at rest.
    spin: Option<(usize, (V3, V3))>,
    /// What of its mass off its spin's axis, and of its second moment not round about it, was read as
    /// none within what its measure is known to (known_round): g mm, g mm^2. Not known, so not none.
    unknown: [f64; 2],
    /// What each of its ten measured moments' error (as Model::compile
    /// reads it) moves its mass, centre and inertia by.
    /// An effort is linear in the moments, so each moves it by what its own
    /// change does, whatever the others' (Moment::mass).
    measure: Vec<(f64, V3, M3)>,
    of: Vec<usize>,
    /// What moves with a balancing joint that no accounted joint carries
    /// (Model::unseparated): its effort is also summed apart.
    apart: bool,
    /// The joints released for springs carrying it (springs.rs), by index.
    released: Vec<usize>,
}

/// `(load NAME PATH (force fx fy fz) (torque tx ty tz) (at FRAME))`, as the
/// model reads it: the force and torque in the world's axes, the point it
/// acts at (FRAME's origin, in the world), the placement of what PATH
/// places (which the point moves with) with its partials along the
/// coordinates, and the reactions it is of.
struct Applied {
    force: [Num; 3],
    torque: [Num; 3],
    at: [Num; 3],
    on: [Num; 12],
    dq: Vec<(usize, [Num; 12])>,
    of: Vec<usize>,
    /// It acts on an unseparated body (Body::apart).
    apart: bool,
    /// The path it acts on, and the joints released for springs carrying
    /// it (springs.rs), by index.
    path: String,
    released: Vec<usize>,
}

/// A relation's coupling (`(couple ...)`), as the model reads it: its ends'
/// paths, its lines (a point, a direction and a pitch, in the world), its
/// pretension and efficiencies (forward and back; None: not given), on
/// which side each reaction has it (+1: it carries B alone, -1: A alone, 0:
/// both or neither), and the joint A rides, where the coupling's loss is
/// charged.
struct Coupled {
    name: String,
    /// Each line's point, direction, pitch and length (mm; 1 where not
    /// given).
    lines: Vec<([Num; 3], [Num; 3], Num, Num)>,
    /// Its first line's index among every coupling's.
    first: usize,
    pretension: Num,
    held: bool,
    efficiency: Option<(Num, Num)>,
    side: Vec<i8>,
    /// Whether each reaction carries A, and B.
    carries: Vec<[bool; 2]>,
    /// A and B.
    ends: [End; 2],
}

impl Coupled {
    /// Whether a joint one of its ends rides moves along any of `coords`.
    fn moves(&self, coords: &[usize]) -> bool {
        self.ends.iter().filter_map(|e| e.input.as_ref()).any(|i| i.along.iter().any(|(k, _)| coords.contains(k)))
    }

    /// Why its numbers at `ctx` are refused, if they are: its pretension is
    /// a tension, at least 0; its efficiency forward in (0, 1], back at most
    /// 1 (an effort's and a mode's alike).
    pub(crate) fn unsound(&self, ctx: &Ctx) -> Option<String> {
        let (t0, eff) = (self.pretension.at(ctx), self.efficiency.as_ref().map(|(e, b)| (e.at(ctx), b.at(ctx))));
        if !(t0 >= 0.0) {
            return Some(format!("the pretension of {} is {t0}: it is a tension, at least 0", self.name));
        }
        let (e, b) = eff.filter(|&(e, b)| !(e > 0.0 && e <= 1.0 && b <= 1.0))?;
        Some(format!("the efficiency of {} is {e} forward and {b} back: forward it is in (0, 1], back at most 1", self.name))
    }
}

/// A coupling's end: its instance path, its placement in the world with
/// its partials along the coordinates (a point of it moves as the
/// placement does: a line's work), and the joint it rides.
struct End {
    path: String,
    on: [Num; 12],
    dq: Vec<(usize, [Num; 12])>,
    input: Option<Input>,
    /// The joints carrying the one it rides, out to the target (a
    /// gantry's slide carrying an idler's axle): power the coupling passes
    /// through them is power it passes too.
    carriers: Vec<Input>,
}

impl End {
    /// The joints moving it: its own, then those carrying that.
    fn chain(&self) -> impl Iterator<Item = &Input> {
        self.input.iter().chain(&self.carriers)
    }
}

/// The joint a coupling's end rides: its frame in the world, its axis and
/// kind, its drive's rate (deg/s or mm/s) and partial along each coordinate
/// moving it.
struct Input {
    joint: String,
    frame: [Num; 12],
    axis: usize,
    turn: bool,
    rate: Num,
    along: Vec<(usize, Num)>,
}

/// A passive joint a coupling reaches: its reaction (index) and its motion's
/// axis and kind. It passes nothing along its motion: its couplings do.
struct Balancing {
    reaction: usize,
    axis: usize,
    turn: bool,
}

/// A path whose reaction is asked for: the bodies under it and those placed
/// on what it carries (Instance::on), every one counted (static ones by
/// their weight). Its frame in the world (R row major, then T): the joint's
/// (its WHERE moved by its motion) where it has one, else its placement;
/// the target's, the world's.
pub struct Reaction {
    pub path: String,
    frame: [Num; 12],
    /// The bodies (item labels) it carries, measured or not: those under
    /// it and those placed on what it carries (`carried`).
    pub carries: Vec<String>,
    /// Bodies it carries that have no mass, or that are not proven
    /// bounded (measured over the window alone): it is not found.
    pub without_mass: Vec<String>,
    pub unbounded: Vec<String>,
    /// Bodies it carries whose shape moves (deforming), with their mass:
    /// their inertia is left out.
    pub leaves_out: Vec<String>,
    /// The coordinates moving what it carries, by index.
    coords: Vec<usize>,
    /// Bodies moving along the coordinates the drive of the joint it rides
    /// reads (its own; else the nearest it is under or placed on) that it
    /// does not carry (not placed on its frames) and that do not support it,
    /// and bodies placed partly on what it carries: what they exert on it
    /// is in it, unseparated, so where any is listed the reaction is not
    /// the joint's load alone.
    pub moves_with: Vec<String>,
    /// Bodies it carries placed by a formula of t other than through a
    /// coordinate or a joint (Model::held_by_nothing).
    pub moved_by_t: Vec<String>,
    /// Carried belts whose material lies on parts it carries, and which it
    /// does not carry: their inertia is in it (Model::belts).
    pub belts: Vec<String>,
}

/// What the rest of the model exerts on what a path carries at a moment,
/// in the world: its force (N), its moment about the world's origin (N
/// mm), the path's frame there (Reaction), and the rounding's scales of
/// the force and the moment.
#[derive(Clone)]
pub struct Wrench {
    pub force: V3,
    pub moment: V3,
    pub frame: [f64; 12],
    noise: (f64, f64),
}

impl Wrench {
    /// One with no rounding, in `frame`: for a split by hand, or a change of one.
    pub(crate) fn of(force: V3, moment: V3, frame: [f64; 12]) -> Wrench {
        Wrench { force, moment, frame, noise: (0.0, 0.0) }
    }

    /// Its force and moment in `frame` (a point p of it is at R p + T in
    /// the world: R row major, then T), about its origin, in its axes.
    pub fn in_frame(&self, frame: &[f64]) -> (V3, V3) {
        let (r, o) = split(frame);
        let oxf = cross(o, self.force);
        let m: V3 = std::array::from_fn(|k| self.moment[k] - oxf[k]);
        (mv(&tr(&r), &self.force), mv(&tr(&r), &m))
    }

    /// In its path's frame.
    pub fn local(&self) -> (V3, V3) {
        self.in_frame(&self.frame)
    }

    /// In a frame placed in its path's (`local`: R row major, then T, as `in_frame`'s), about that frame's
    /// origin, in its axes: a frame a part on the joint declares.
    pub fn in_local(&self, local: &[f64; 12]) -> (V3, V3) {
        let ((rj, tj), (rl, tl)) = (split(&self.frame), split(local));
        let r = mul(&rj, &rl);
        let t = mv(&rj, &tl);
        let mut world = [0.0; 12];
        for i in 0..3 {
            world[3 * i..3 * i + 3].copy_from_slice(&r[i]);
            world[9 + i] = t[i] + tj[i];
        }
        self.in_frame(&world)
    }

    /// What a moment's reaction is known to, in its force and in its moment
    /// about the world's origin: ROUNDING of what they sum (`Model::at`
    /// leaves it so scaled).
    pub(crate) fn rounding(&self) -> (f64, f64) {
        self.noise
    }

    /// What its arithmetic can leave in its force and in its moment about
    /// the world's origin, in any direction: FLOAT of what they sum.
    pub(crate) fn float(&self) -> (f64, f64) {
        (self.noise.0 / ROUNDING * FLOAT, self.noise.1 / ROUNDING * FLOAT)
    }
}

/// A spinning body read round within what its measure is known to (`Model::made_round`).
pub struct MadeRound {
    pub body: String,
    /// What its first and second moments across its axis are known to (g mm, g mm^2), each where
    /// that part of it was read round (0 where it was read as measured).
    pub first: f64,
    pub second: f64,
}

/// A spinning body whose centre is not exact: how far its measure leaves it (mm), and the reactions carrying it.
#[derive(Clone, Debug)]
pub struct SpinCentre {
    pub body: String,
    pub within: f64,
    pub of: Vec<usize>,
}

/// A body whose material moves as its `carried` says (../kernel/docs/IR.md, Paths):
/// material station m rides the body's placement and its path at station
/// m + SHIFT, its mass on the path line. Its effort is its material's,
/// integrated along it: Gauss-Legendre on each piece of the path the
/// material covers at the moment.
struct Strip {
    label: String,
    /// Mass per mm along the material (g): its profile measured cut to its
    /// length, spread evenly along it.
    mass: f64,
    /// How far its mass per mm can be from `mass` (g): its measure's proven
    /// bounds, at what its density is known to. Infinite: not known.
    err: f64,
    shift: Num,
    length: Num,
    closed: bool,
    total: Num,
    /// Per piece of the path: its start station and length, its sweep
    /// (an arc's), and its point x mm along it.
    starts: Vec<Num>,
    lengths: Vec<Num>,
    sweeps: Vec<Option<Num>>,
    pieces: Vec<Point>,
    /// Each piece's two ends in the world (its station 0 and its length).
    ends: Vec<[[Num; 3]; 2]>,
    /// Each piece's coupling parts (Model::belts): the span it lies on, or
    /// the wheel it wraps; empty where its material is not separated.
    on: Vec<Rides>,
    of: Vec<usize>,
    apart: bool,
}

/// What a piece of a carried belt's path lies on (Model::belts): a span, a
/// coupling's line with a length (its material shared between the line's
/// two ends by where it is along it, as a string's load is between its
/// ends), or a wheel it wraps (all of it on the end a span meets there).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Rides {
    /// Coupling, line.
    Span(usize, usize),
    /// Coupling, line, end (0: A, 1: B) where the wrap meets the line.
    Wrap(usize, usize, usize),
}

/// What an effort reads of a point mass: its position and acceleration,
/// and its partial velocity along each coordinate that moves it.
struct Point {
    pos: [Num; 3],
    acc: [Num; 3],
    dq: Vec<(usize, [Num; 3])>,
}

/// A coordinate: its name (`INSTANCE.ARG`) and its value as passed, with
/// the value's rate, acceleration and jerk.
pub struct Coordinate {
    pub name: String,
    q: Num,
    dq: Num,
    ddq: Num,
    /// Its jerk, whose rate bounds steps (Model::steps), and the jerk's
    /// rate (a stop's direction where the lower ones vanish: path.rs).
    dddq: Num,
    d4q: Num,
    /// Its value as passed reads t: it moves (a constant is a pose held).
    moves: bool,
}

/// A joint with effort limits under the target, by its instance path: its
/// kind, its drive's rate along each coordinate, its rate in t, and its
/// limits.
pub struct Actuator {
    pub path: String,
    pub turn: bool,
    along: Vec<Num>,
    rate: Num,
    lo: Num,
    hi: Num,
    /// The coordinates its effort rests on: those it moves, and those
    /// coupled with them (moving one body together, or read by one
    /// actuator or one deforming body), by index.
    coupled: Vec<usize>,
    /// What its effort leaves out: the deforming bodies moving with those
    /// coordinates, each with its mass.
    pub leaves_out: Vec<String>,
    /// Where its limits are written (file:line).
    pub source: String,
    /// Its `(command PATH EFFORT)`, the effort it supplies in a simulation
    /// (clamped to its limits there: simulate.rs); None: none.
    command: Option<Num>,
}

/// A model's coordinates' params at a state (Doc::with_states): each one's
/// position, and where it moves its rate's and acceleration's; and the
/// moment's, t0.
struct States {
    of: Vec<(usize, Option<(usize, usize)>)>,
    t0: usize,
}

/// Over a step of the written motion, a coordinate's: bounds on its rate,
/// its rate's rate and that one's (None: it may jump in the step: a select
/// or floor on t switches there), and its rate's and acceleration's
/// bounds over the step (Model::path_step).
pub struct PathStep {
    pub speed: Option<f64>,
    pub accel: Option<f64>,
    pub jerk: Option<f64>,
    pub rate: I,
    pub acc: I,
}

/// A joint's `(speed-limits LO HI)` under the target, by its instance path:
/// its drive's rate and its limits (deg/s or mm/s).
pub struct Speed {
    pub path: String,
    pub(crate) rate: Num,
    pub(crate) lo: Num,
    pub(crate) hi: Num,
}

/// Bodies spinning with one link, not all shown round about its axis (a
/// rotor's D-cut shaft, the pulley and grub screw on its shaft frame): one
/// rigid mass turning in the frame the link turns in. Their effort reads
/// the spin's angle only through their combined departure from round
/// about the axis (Model::departure), so where that can move no effort by
/// more than the steps' share of its peak, the steps read them with the
/// spin at rest (Model::choose).
struct Group {
    /// The link's instance path.
    link: String,
    /// Each body not shown round: its path, its index, and its placement in
    /// the frame at rest (a body point p is R p + T there).
    members: Vec<(String, usize, M3, V3)>,
    /// The axis the link turns about, of the frame it turns in.
    axis: usize,
    /// That frame at rest in the world (R row major, then T), with its
    /// derivatives along t (twice) and along each coordinate that moves it.
    frame: [Num; 12],
    dt: [Num; 12],
    dt2: [Num; 12],
    dq: Vec<(usize, [Num; 12])>,
    /// The spin's rate and acceleration (radians), and its partial along
    /// each coordinate that turns it.
    rate: Num,
    acc: Num,
    along: Vec<(usize, Num)>,
    /// Its bodies move along the coordinates (an effort reads them); else
    /// only reactions carry them (t spins them).
    effort: bool,
}

/// The mechanism under a target, ready to answer at any moment.
pub struct Model {
    /// m/s^2.
    pub gravity: [f64; 3],
    pub coords: Vec<Coordinate>,
    pub actuators: Vec<Actuator>,
    bodies: Vec<Body>,
    strips: Vec<Strip>,
    loads: Vec<Applied>,
    /// The paths whose reactions were asked for, in the order asked, then
    /// the passive joints couplings reach (Balancing), not asked.
    pub reactions: Vec<Reaction>,
    /// The relations' couplings under the target, and the passive joints
    /// whose balance gives their tensions.
    couplings: Vec<Coupled>,
    balances: Vec<Balancing>,
    /// What moves with the balances' coordinates that no accounted joint
    /// carries (Reaction::moves_with of each): the tensions leave out what
    /// it exerts.
    pub unseparated: Vec<String>,
    /// Each coupling's train: those whose tensions rest on one another's
    /// (a balance reads both, or they share a part), by its least member.
    trains: Vec<usize>,
    /// Groups never read with their spin at rest: carried by a joint whose
    /// balance gives a coupling's tension, which their departure from round
    /// (where it is not none) would move by no bound found.
    behind: Vec<bool>,
    /// What an effort rests on that moves, each once: every placement's
    /// rotation, its translation's acceleration, its partial velocities, a
    /// spinning body's with its spin at rest and its spin's partials (Spin;
    /// Model::steps keeps each within 1% inside a step).
    tracked: Vec<Tracked>,
    groups: Vec<Group>,
    /// The params its efforts read, by name (its coordinates' own aside),
    /// and those only its reactions read besides.
    pub reads: Vec<String>,
    pub reads_carried: Vec<String>,
    /// The file's span (its `(t-span ...)`), over which its effort claims
    /// read each moment (Model::read).
    pub span: Option<(f64, f64)>,
    /// The parameters it is read at (its coordinates' at 0).
    params: Vec<f64>,
    /// Bodies placed by t other than through a coordinate: no actuator is
    /// charged for that motion.
    pub driven_by_nothing: Vec<String>,
    /// Those placed by a formula of t riding no joint: nothing the model
    /// declares holds them (Model::moved_by_t_note).
    pub held_by_nothing: Vec<String>,
    /// Passive joints moved by an actuator's coordinate through no relation (a belt, a cord, a mesh) and
    /// not riding its moving frame, each with those actuators: nothing between them is in the model.
    pub unlinked: Vec<(String, Vec<String>)>,
    /// Spinning bodies not shown round whose departure from round about
    /// their axis lies, in some part, within what their measure is known to
    /// (known_round): that part is read as none, and what it may be is not
    /// known, which the reactions carrying the body say (`imbalance`).
    pub made_round: Vec<MadeRound>,
    /// Spinning bodies whose centre is not exact (spin_centre_note).
    pub spin_centres: Vec<SpinCentre>,
    /// Bodies whose shape moves (reads t, or a coordinate that moves):
    /// their inertia is not counted. Each with its mass (g) at the moment
    /// it was read, where it has one.
    pub deforming: Vec<(String, Option<f64>)>,
    /// Bodies whose material moves as their `carried` says, counted so,
    /// each with its mass (g), measured from its profile.
    pub carried: Vec<(String, f64)>,
    /// Moving bodies whose mass is measured, and those of them weighed from
    /// a published mass, or whose spin is published (`spin-inertia`).
    pub measured: usize,
    pub published: usize,
    pub spin_inertia: Vec<String>,
    /// The finest and coarsest cells the masses were measured at (mm), and
    /// those that did not settle, why: of the bodies an effort reads, and of
    /// those only reactions carry.
    pub res: (f64, f64),
    pub unsettled: Vec<String>,
    pub unsettled_carried: Vec<String>,
    /// Built leaving out the moving bodies that have no mass (`compile_kinematic`): each one's label, its
    /// definition and the coordinates it moves along. Its motion is whole (every coordinate's value and
    /// rate, each moment read); an effort or a reaction such a body would change is refused, naming it
    /// (`blocked`, `blocks_reaction`), and no other is changed by leaving it out.
    pub massless: Vec<(String, String, Vec<usize>)>,
    /// The set each coordinate is coupled in (those moving one body together, read by one actuator or one
    /// deforming body, or joined by a coupling): two of one number rest on each other's efforts.
    sets: Vec<usize>,
    /// The solves it reads that vary (index in Doc::solves), in their
    /// order, set again at each moment from `params` (those on t followed
    /// along t: Follow).
    solves: Vec<usize>,
    follow: Option<Arc<ink_lang::Follow>>,
    /// In a model of a load with states (Doc::with_states): each
    /// coordinate's params at a state, and every joint's speed limits.
    states: Option<States>,
    /// Each loop's root at the last state a simulation read (state_near),
    /// by its solve's index: where the next is found from.
    roots: Mutex<HashMap<usize, Vec<f64>>>,
    /// The coordinates passed as `(simulate X0 [V0])`, in the order those
    /// are written: each one's index, and the params its value reads (its
    /// position and first four derivatives at a moment, and the moment).
    pub simulated: Vec<(usize, [usize; 6])>,
    /// The bodies whose shape moves (deforming, with a mass), each with
    /// the coordinates it moves with: what an effort along those leaves out.
    leaves: Vec<(String, Vec<usize>)>,
    pub speeds: Vec<Speed>,
    /// Its springs, in a simulation that releases them (springs.rs).
    springs: Option<springs::Springs>,
    /// Moments, by their t and the groups read round in them.
    moments: Mutex<HashMap<(u64, Vec<bool>), Arc<Moment>>>,
    /// Steps chosen, by span and whom for: the efforts, or a claim on a
    /// reaction (Serve::claim).
    chosen: Mutex<HashMap<(u64, u64, Option<usize>), Arc<Chosen>>>,
    /// What the bodies and loads sum into over a step, by the step's ends
    /// (bounded.rs); None where something read has no bound over it.
    swept: Mutex<HashMap<(u64, u64), Option<Arc<bounded::Swept>>>>,
    /// What its efforts read includes a param a solve on t sets at each
    /// moment, which no bound over a step here follows (Model::bounded).
    reads_moment: bool,
    /// Where it is a load's model at its written moments: the file and the
    /// target to build its model of states from, and that model once built
    /// (where an actuator's limits cross, where along its speed they begin
    /// to: Model::crossed).
    pub twin: Option<(ink_lang::Reload, String, Arc<std::sync::OnceLock<Option<Model>>>)>,
}

/// The model at one moment.
#[derive(Clone)]
pub struct Moment {
    pub t: f64,
    /// Per coordinate: its value, rate and acceleration as passed, the
    /// effort along it (N mm per unit of it), and, along one driven or moving,
    /// what of that the actuators do not supply (`unpowered`, where above
    /// rounding).
    pub q: Vec<[f64; 3]>,
    pub f: Vec<f64>,
    pub unpowered: Vec<Option<f64>>,
    /// Per actuator: its effort (N mm a turn, N a slide), its drive's rate
    /// (deg/s or mm/s), its limits there, and the inertia reflected to it
    /// (g mm^2 a turn, g a slide).
    pub effort: Vec<f64>,
    pub rate: Vec<f64>,
    pub lo: Vec<f64>,
    pub hi: Vec<f64>,
    pub reflected: Vec<f64>,
    /// The rounding scale of each coordinate's effort (1e-9 of the terms
    /// summed into it) and of each actuator's (what those carry into it):
    /// a difference below it is no difference.
    noise_f: Vec<f64>,
    noise: Vec<f64>,
    /// Each actuator's effort per unit of each coordinate's (tau = P f).
    pm: Vec<Vec<f64>>,
    /// Read over a span (Model::over), the groups read with their spin at
    /// rest read round: how far each coordinate's effort and each
    /// actuator's may be from the real one, whatever the spins' angles
    /// (Group::moves over the step holding t); 0 at one moment (Model::at).
    pub bound_f: Vec<f64>,
    pub bound: Vec<f64>,
    /// Per actuator, the most its effort can be off by for what its bodies'
    /// measured mass, centre and inertia are known to (Body::measure),
    /// through its couplings' worst loss.
    pub mass: Vec<f64>,
    /// What each body's measured moments, each anywhere in its proven bounds, move each reaction by: per
    /// reaction, one (force, moment about the world's origin) change per measured moment of a body it
    /// carries, through the couplings' tensions too; summed component by component in a claim's frame,
    /// the reaction's band for its masses (Moment::mass_reach). Each change alone, so that each component
    /// sums only its own.
    pub mass_w: Vec<Vec<(V3, V3)>>,
    /// Per coupling line, how far its tension can move for what the bodies' measures are known to.
    pub mass_t: Vec<f64>,
    /// Per actuator, the most its couplings' loss can be off by (Moment::
    /// compute): the tensions it rests on leave the unseparated bodies out.
    pub loss_bound: Vec<f64>,
    /// |I - B^T P|: how much of each coordinate's effort no actuator
    /// carries into each one's unpowered effort; and, read over a span, how
    /// far each unpowered effort may be from the one read (through it).
    uncovered: Vec<Vec<f64>>,
    pub bound_u: Vec<f64>,
    /// Each reaction asked for (then those of the joints couplings reach):
    /// what the rest of the model exerts on what it carries, its couplings'
    /// lines and frictions aside.
    pub reactions: Vec<Wrench>,
    /// Each coupling line's tension (N), in the couplings' order; each
    /// coupling's friction on its parts (Friction); why none is found,
    /// where it is not; and the couplings whose efficiency is not given
    /// that pass power here (their loss unknown).
    pub tensions: Vec<f64>,
    /// Each line's tension's rounding (N): each balance giving it, in its
    /// own units (N mm about a turn's axis, N along a slide), over the
    /// line's coefficient there (its arm about that axis), summed.
    pub tension_noise: Vec<f64>,
    pub frictions: Vec<Friction>,
    pub coupling: Option<String>,
    pub lossless: Vec<usize>,
    /// How the efforts are made, for an answer that varies the state: the
    /// effort along each coordinate before the
    /// couplings' losses, what the unseparated bodies need along each,
    /// each actuator's worst efficiency's share (its loss bound is that
    /// times |P apart|), each actuator's drive along each coordinate (B,
    /// radians for a turn), and where the losses were charged.
    pub f_ideal: Vec<f64>,
    pub apart: Vec<f64>,
    pub worst: Vec<f64>,
    pub bm: Vec<Vec<f64>>,
    pub sites: Vec<Site>,
    /// The generalized inertia along the coordinates (g mm^2 per unit^2 of
    /// them: an effort is 1e-6 of it times their accelerations, plus what
    /// their speeds and gravity need), and each actuator's command there
    /// (None: it has none): what a simulation steps by (simulate.rs).
    pub gm: Vec<Vec<f64>>,
    /// The heft of each of gm's rows, what its diagonal sums: each body's
    /// mass times its speed squared per unit of the row, and its turning
    /// squared times its inertia's diagonal there. A row whose inertia is
    /// within ROUNDING of it carries none: an idler given none, its centre
    /// measured a hair off its axis, sums m r^2 there and m r^2 less in its
    /// own inertia to rounding, of either sign (simulate::forces).
    pub heft: Vec<f64>,
    pub command: Vec<Option<f64>>,
    /// Along each joint released for springs (springs.rs), in order: the
    /// effort what it carries needs along its deviation (N a slide, N mm per
    /// radian a turn: its bodies' motion and weight, the loads on them, as
    /// along a coordinate), and its rounding; `gm` runs over the
    /// coordinates, then those joints. The springs at this moment.
    pub fr: Vec<f64>,
    pub noise_r: Vec<f64>,
    pub springs: Option<Arc<springs::At>>,
    /// Why nothing is found at this moment, if nothing is: a solve the
    /// placements rest on does not solve here, or its rates are not
    /// defined (a dead point). Its numbers are then NaN.
    pub unsolved: Option<String>,
}

impl Moment {
    /// The rounding of coordinate k's effort: 1e-9 of the terms summed
    /// into it.
    pub fn rounding(&self, k: usize) -> f64 {
        self.noise_f[k]
    }

    /// A fixed-effort rating is an equality. A nonzero discrepancy within
    /// the rounding of this actuator's own load cannot witness its failure.
    pub(super) fn fixed_margin_unresolved(&self, a: usize) -> bool {
        self.lo[a] == self.hi[a] && self.margin(a) < 0.0 && -self.margin(a) <= self.noise[a]
    }

    /// min(effort - lo, hi - effort) of actuator a.
    /// How far reaction r's six components in `frame` (force, then moment about the frame's origin, in its
    /// axes, as Wrench::in_frame) can be off for what the bodies it carries are measured to: each body's
    /// measured moment's change, component by component (`mass_w`).
    pub fn mass_reach(&self, r: usize, frame: &[f64; 12]) -> [f64; 6] {
        let (rot, o) = split(frame);
        let rt = tr(&rot);
        self.mass_w[r].iter().fold([0.0; 6], |mut out, (df, dm)| {
            let odf = cross(o, *df);
            let (f, m) = (mv(&rt, df), mv(&rt, &std::array::from_fn(|k| dm[k] - odf[k])));
            for k in 0..3 {
                out[k] += f[k].abs();
                out[3 + k] += m[k].abs();
            }
            out
        })
    }

    /// `mass_reach` through a linear map from reaction r (force, then moment about its frame's origin,
    /// world axes: split::split's) to each row's value.
    pub fn mass_reach_mapped(&self, r: usize, rows: &[[f64; 6]]) -> Vec<f64> {
        let o = [self.reactions[r].frame[9], self.reactions[r].frame[10], self.reactions[r].frame[11]];
        rows.iter()
            .map(|row| {
                self.mass_w[r]
                    .iter()
                    .map(|(df, dm)| {
                        let odf = cross(o, *df);
                        (0..3).map(|k| row[k] * df[k] + row[3 + k] * (dm[k] - odf[k])).sum::<f64>().abs()
                    })
                    .sum()
            })
            .collect()
    }

    pub fn margin(&self, a: usize) -> f64 {
        (self.effort[a] - self.lo[a]).min(self.hi[a] - self.effort[a])
    }
}

/// What rounding is, of a sum: this much of the magnitudes summed.
pub(crate) const ROUNDING: f64 = 1e-9;

/// What the arithmetic itself can leave of nothing, as a share of the
/// magnitudes it sums: f64's 2.2e-16 over the few thousand operations a
/// reaction and a solve on it pass through. Where a part of an answer is
/// zero by statics (what no support carries), a large number elsewhere in
/// the problem hides no more than this of itself in it.
pub(crate) const FLOAT: f64 = 4096.0 * f64::EPSILON;

/// Why bodies have no mass, and what gives them one.
pub const NO_MASS: &str = "give each a material or density in the file defining it, or the published mass of an assembly placing it, (attr mass M): \
    an effort rests on every moving body's inertia, so none is guessed";

impl Coordinate {
    /// Its value as passed reads t: it moves (else a pose held).
    pub fn moves(&self) -> bool {
        self.moves
    }
}

/// Disjoint sets of coordinates, by index (union-find).
struct Sets(Vec<usize>);

impl Sets {
    fn find(&mut self, mut k: usize) -> usize {
        while self.0[k] != k {
            self.0[k] = self.0[self.0[k]];
            k = self.0[k];
        }
        k
    }

    fn join(&mut self, ks: &[usize]) {
        for w in ks.windows(2) {
            let (a, b) = (self.find(w[0]), self.find(w[1]));
            self.0[a] = b;
        }
    }
}
