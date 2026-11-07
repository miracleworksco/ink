//! A field's turn symmetry about the z axis of its own coordinates, shown
//! from how the field is built: never sampled.
//!
//! `turn_symmetry(g, f)` gives the largest order n it can show: f, as a
//! function of real numbers (each node its operation on reals, under the
//! core's value rules: x / 0 = 0, sqrt of a negative 0, ...; each constant
//! the float it is), is unchanged by a turn of 2 pi / n about z through the
//! origin of f's own coordinates. `Order::Every` is every turn (a function
//! of x^2 + y^2, z, t and params); `Of(1)` is nothing shown. For every
//! point p, integer k and every t and param value:
//!
//!   f(R(2 pi k / n) p) = f at p with each angular fold's argument moved by
//!   at most `slack` radians (`Symmetry::slack`; 0 with no fold).
//!
//! A fold is `mod(a, P)` (std's `(- a (* P (floor (/ a P))))`) of an angle
//! a about z (atan2 of a turning pair, plus anything invariant), with P a
//! constant within 1e-9 rad of 2 pi / m: P is a float, not 2 pi / m, so m
//! folds of it close the circle short or long by |m P - 2 pi|, and a turn
//! by 2 pi k / n crossing atan2's seam lands that much off a whole number
//! of folds. `slack` is the |m P - 2 pi| (rounded up) of the fold the
//! field reads: about 1e-15 rad. It reads one fold (one angle, one pitch,
//! however many times), or nothing is shown: two folds are each off by
//! their own error after a turn, and what combines them (a select of
//! their difference) may change by any amount, so no one rotation of the
//! point bounds the change. A caller that needs a distance turns
//! it into one. With one fold, moving its argument by e is turning the
//! point by e, but across atan2's seam, which costs another slack: so
//! f(R(2 pi k / n) p) = f(R(e) p) with |e| <= 2 slack, a point moved at
//! most 2 r slack along its circle (3e-13 mm at r = 100); a field with a
//! slope bound L changes by L times that.
//!
//! Each rule is an identity in real numbers:
//! - constants, params, z, t: every turn;
//! - x, y: linear forms A x + B y (coefficients products of invariant
//!   values, compared as written, never computed: a product of constants
//!   and params alone is its float, one value, not its factors); a turn
//!   about z, a mirror or a scale of them stays one;
//! - two forms u = A x + B y, v = -B x + A y (or v = B x - A y): u^2 + v^2
//!   = (A^2 + B^2)(x^2 + y^2), every turn their coefficients keep (so
//!   `len2 x y`, a cylinder, a rotated one);
//! - atan2 v u of such a pair: an angle, which a turn by a shifts by +-a
//!   (mod 2 pi; at u = v = 0 it is 0 and stays 0);
//! - mod(angle, P) with P within 1e-9 of 2 pi / m: order m (above);
//! - any operation whose operands are all invariant: the gcd of their
//!   orders (`Every` neutral): a function of invariant values is
//!   invariant (of one fold, above: a 12-fold and an 8-fold part together
//!   show nothing); a fold's turned coordinates r cos a, r sin a (std's repeat-polar) are
//!   m-fold, and so is anything warped by them.
//! - tag, lipschitz, a derivative's mark and a sampled field's read are
//!   their values.
//!
//! Anything else (a form plus a constant: a move off the axis; |x|; x y;
//! a pitch in degrees or not 2 pi / m; a fold of a param) is order 1.
//! Copies at float angles (std's pattern-polar, a bolt circle) are not
//! exactly symmetric in reals: order 1.
//!
//! Whether a placed body turns about this axis (its joint's axis is its own
//! z through its origin) is the caller's question: this is about the field
//! in the coordinates it is given in.
//!
//! `carried(g, f, centre, drive)` shows what of f every turn about the
//! vertical axis through `centre` leaves unchanged with an input, the
//! drive, moved with it: a belt's teeth round its pulley, each tooth's
//! station its angle about the centre times the pitch radius, less the
//! belt's shift, which its drive moves. Its rules are the ones above,
//! and:
//! - x - cx and y - cy are the linear forms (bare x where cx is 0); an
//!   angle a phase reads with a drive value is of forms with nonzero
//!   constant coefficients (one that is a value may be 0 somewhere, where
//!   atan2 is 0 and does not turn);
//! - the drive, and any function of it and values every turn keeps, is a
//!   drive value: the turn does not keep it;
//! - a phase is k angle + sum_i c_i S_i + a value every turn keeps (an
//!   angle is k = +-1), S_i drive values, k and c_i bounds of the real
//!   coefficients the field's constants make (products and quotients by
//!   constants, rounded outward); two angles are not added (each crosses
//!   atan2's seam where the other does not);
//! - cos or sin of a phase that reads a drive value, with k within 1e-9 of
//!   a whole n != 0, is a function of the phase alone (`Coupling`). Across
//!   atan2's seam the angle jumps by 2 pi and the phase by 2 pi k, 2 pi n
//!   give or take 2 pi |k - n|: `Coupling::seam`, radians of the phase per
//!   crossing. One phase is read so (one node): two would each be off by
//!   their own seam.
//!
//! So each node of `Carried::invariant` reads the angle and the drive only
//! through the phase: for every point p off the axis, b and values S, S' of the drive
//! values, it is the same at (R_c(b) p, S') as at (p, S) with the phase
//! moved by k b + sum_i c_i (S'_i - S_i), and by at most `seam` per
//! crossing of atan2's seam the turn makes; and the phase moved by e is
//! the point turned by e / k. With no drive (None) nothing reads one, and
//! the nodes are those `turn_invariant` gives about the centre.

use crate::graph::{BinOp, Graph, Node, NodeId, UnOp, Var};
use crate::interval::{self, I};
use std::collections::HashMap;

/// Turns about z that leave a field unchanged: every multiple of 2 pi / n
/// (`Of(n)`, n >= 1), or every turn.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    Of(u64),
    Every,
}

impl Order {
    /// The turns both keep.
    pub fn gcd(self, o: Order) -> Order {
        match (self, o) {
            (Order::Every, o) | (o, Order::Every) => o,
            (Order::Of(mut a), Order::Of(mut b)) => {
                while b != 0 {
                    (a, b) = (b, a % b);
                }
                Order::Of(a)
            }
        }
    }
}

const NONE: Order = Order::Of(1);

/// What `turn_symmetry` shows (see the module's contract).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Symmetry {
    pub order: Order,
    /// Radians each fold's argument may be off after a turn of the order
    /// (0 without folds, or of order 1).
    pub slack: f64,
    /// The angle the field's one fold reads (its argument a in mod(a, P)),
    /// where it has one: every turn of the order is the field with that
    /// angle moved by at most `slack` (and its multiple of the pitch), so a
    /// caller covers the turns by ranging a shift added to it.
    pub fold: Option<NodeId>,
}

/// `root`'s turn symmetry about z in its own coordinates (the module's
/// contract).
pub fn turn_symmetry(g: &Graph, root: NodeId) -> Symmetry {
    let mut a = Analysis::new(g, [0.0; 2], None);
    for id in value_nodes(g, root) {
        let c = a.classify(id);
        a.class.insert(id, c);
    }
    let order = match a.class[&root] {
        // Two folds (of other angles or pitches) are each off by their own
        // error after a turn, and what combines them (a select of their
        // difference) may change by any amount: no one rotation of the
        // point bounds the change, so nothing is shown.
        Class::Inv(_) if a.folds.len() > 1 => NONE,
        Class::Inv(o) => o,
        _ => NONE,
    };
    let fold = match a.folds[..] {
        [(angle, _)] if order != NONE => Some(angle),
        _ => None,
    };
    Symmetry { order, slack: if order == NONE { 0.0 } else { a.slack }, fold }
}

/// The nodes of `root`'s value that read x or y and that every turn about
/// z leaves unchanged (`Order::Every`, by the module's rules: a function of
/// x^2 + y^2, z, t and params): each is the same at R p as at p for every
/// turn R about z, exactly in real numbers.
pub fn turn_invariant(g: &Graph, root: NodeId) -> std::collections::HashSet<NodeId> {
    carried(g, root, [0.0; 2], None).invariant
}

/// What `carried` shows (the module's contract).
#[derive(Clone, Debug)]
pub struct Carried {
    /// The nodes reading x or y that every turn about the axis leaves
    /// unchanged but for the one phase they read (the module's contract).
    pub invariant: std::collections::HashSet<NodeId>,
    /// That phase, where they read the drive.
    pub coupling: Option<Coupling>,
}

/// The one phase k angle + sum_i c_i S_i + c that `Carried::invariant`
/// reads periodically: k and each drive value S_i (a node reading the
/// drive and values every turn keeps) with c_i, bounds of the real
/// coefficients, and `seam`, radians the phase is off a whole turn of it
/// per crossing of atan2's seam.
#[derive(Clone, Debug)]
pub struct Coupling {
    pub k: I,
    pub drives: Vec<(NodeId, I)>,
    pub seam: f64,
}

/// What of `root` turns about the vertical axis through `centre` with
/// `drive` (the module's contract; None: no input is one).
pub fn carried(g: &Graph, root: NodeId, centre: [f64; 2], drive: Option<Var>) -> Carried {
    let mut a = Analysis::new(g, centre, drive);
    let ids = value_nodes(g, root);
    for &id in &ids {
        let c = a.classify(id);
        a.class.insert(id, c);
    }
    let invariant = ids.into_iter().filter(|id| matches!(a.class[id], Class::Inv(Order::Every)) && g.depends_on(*id, |v| matches!(v, Var::X | Var::Y))).collect();
    Carried { invariant, coupling: a.coupling.map(|c| c.1) }
}

/// The centres `root`'s parts are written about, on each axis (x - c and
/// y - c, a bare x or y about 0), ascending: where `carried` may find
/// parts every turn about a vertical axis leaves unchanged.
pub fn centres(g: &Graph, root: NodeId) -> [Vec<f64>; 2] {
    let mut out = [vec![0.0], vec![0.0]];
    for id in value_nodes(g, root) {
        if let Some((i, c)) = about(g, id) {
            out[i].push(c);
        }
    }
    for c in &mut out {
        c.sort_by(f64::total_cmp);
        c.dedup();
    }
    out
}

/// x - c or y - c (x + -c ...): the axis (0 or 1) and c.
fn about(g: &Graph, id: NodeId) -> Option<(usize, f64)> {
    let axis = |v: NodeId| match g.nodes[v as usize] {
        Node::Var(Var::X) => Some(0),
        Node::Var(Var::Y) => Some(1),
        _ => None,
    };
    match g.nodes[id as usize] {
        Node::Bin(BinOp::Sub, v, c) => Some((axis(v)?, g.as_const(c)?)),
        Node::Bin(BinOp::Add, u, v) => match (axis(u), axis(v)) {
            (Some(i), _) => Some((i, -g.as_const(v)?)),
            (_, Some(i)) => Some((i, -g.as_const(u)?)),
            _ => None,
        },
        _ => None,
    }
}

/// The operands a node's value is made of (a promise's frame and a
/// derivative's field only bound it).
fn value_children(n: Node) -> ([NodeId; 3], usize) {
    match n {
        Node::Tag(_, a) | Node::Lipschitz(a, ..) | Node::Deriv(a, ..) | Node::Read(_, _, a, _) | Node::Real(a) => ([a, 0, 0], 1),
        Node::Frame(_) => ([0; 3], 0),
        n => n.children(),
    }
}

/// The nodes `root`'s value reads, ascending (operands first).
fn value_nodes(g: &Graph, root: NodeId) -> Vec<NodeId> {
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![root];
    while let Some(n) = stack.pop() {
        if seen.insert(n) {
            let (ch, k) = value_children(g.nodes[n as usize]);
            stack.extend(&ch[..k]);
        }
    }
    let mut out: Vec<NodeId> = seen.into_iter().collect();
    out.sort_unstable();
    out
}

/// Whether `id` is a constants-and-params value, taken in floats
/// (Tape::fixed): every node it reads is a constant, a param, or an
/// operation or tag of those.
fn fixed(g: &Graph, id: NodeId) -> bool {
    g.reachable(&[id]).iter().all(|&n| matches!(g.nodes[n as usize], Node::Const(_) | Node::Var(Var::Param(_)) | Node::Un(..) | Node::Bin(..) | Node::Select(..) | Node::Tag(..)))
}

/// A factor of a coefficient: a constant's magnitude (its bits) or an
/// invariant node, or (`true`) its reciprocal under the value rules (0 at 0);
/// or pi / 180 exactly (UnOp::Rad; `true`: 180 / pi, UnOp::Deg).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Factor {
    Num(u64, bool),
    Node(NodeId, bool),
    Pi180(bool),
}

/// A real coefficient as written: its sign times the product of its
/// factors (sorted, so equal products compare equal).
#[derive(Clone, Debug, PartialEq, Eq)]
struct Coef {
    neg: bool,
    factors: Vec<Factor>,
}

impl Coef {
    fn one() -> Coef {
        Coef { neg: false, factors: Vec::new() }
    }
    fn negated(&self) -> Coef {
        Coef { neg: !self.neg, factors: self.factors.clone() }
    }
}

/// A x + B y (`None`: a coefficient 0), whose coefficients every turn of
/// `order` keeps.
#[derive(Clone, Debug, PartialEq)]
struct Lin {
    a: Option<Coef>,
    b: Option<Coef>,
    order: Order,
}

fn neg_coef(c: &Option<Coef>) -> Option<Coef> {
    c.as_ref().map(Coef::negated)
}

impl Lin {
    fn negated(&self) -> Lin {
        Lin { a: neg_coef(&self.a), b: neg_coef(&self.b), order: self.order }
    }
    /// Whether (u, v) = (self, v) turns with the point (v = -B x + A y) or
    /// against it (v = B x - A y): either way u^2 + v^2 = (A^2 + B^2)(x^2 +
    /// y^2), and atan2 v u moves by the turn's angle, one way or the other.
    fn pairs(&self, v: &Lin) -> bool {
        (self.a.is_some() || self.b.is_some()) && ((v.a == neg_coef(&self.b) && v.b == self.a) || (v.a == self.b && v.b == neg_coef(&self.a)))
    }
}

#[derive(Clone, Debug)]
enum Class {
    /// Unchanged by the turns of the order.
    Inv(Order),
    Lin(Lin),
    /// A value unchanged by the turns of the order, plus sign * lin^2 for
    /// each term (never empty: what pairs up is invariant).
    Squares(Order, Vec<(bool, Lin)>),
    /// k angle + ks S + a value unchanged by the turns of the order
    /// (`Phase`); an angle about the axis is k = +-1 with no S.
    Phase(Phase),
    /// x or y where the centre's is not 0: only x - cx is a form.
    Raw,
    /// A drive value: the drive, or a function of it and values every turn
    /// keeps.
    Drive,
}

/// k angle + sum_i c_i S_i + an invariant value: a turn by b moves it by
/// k b, and across atan2's seam by 2 pi k more; each drive value S_i moved
/// moves it by c_i times that. k and c_i bound the real coefficients.
/// `fixed`: its angle's forms have nonzero constant coefficients (one
/// that is a value could be 0 somewhere, where atan2 is 0 and turns not).
#[derive(Clone, Debug)]
struct Phase {
    k: I,
    drives: Vec<(NodeId, I)>,
    order: Order,
    fixed: bool,
}

impl Phase {
    fn scaled(self, c: I) -> Phase {
        let m = |x: I| interval::mul_o(x, c);
        Phase { k: m(self.k), drives: self.drives.into_iter().map(|(s, ks)| (s, m(ks))).collect(), ..self }
    }

    /// Plus drive value `s` times c.
    fn plus(mut self, s: NodeId, c: I) -> Phase {
        match self.drives.iter_mut().find(|d| d.0 == s) {
            Some(d) => d.1 = interval::add_o(d.1, c),
            None => self.drives.push((s, c)),
        }
        self
    }

    /// An angle (k = +-1, no drive value): what a fold reads.
    fn angle(&self) -> bool {
        self.drives.is_empty() && self.k.lo == self.k.hi && self.k.lo.abs() == 1.0
    }
}

struct Analysis<'g> {
    g: &'g Graph,
    class: HashMap<NodeId, Class>,
    slack: f64,
    /// The folds read: each one's angle and pitch (its bits).
    folds: Vec<(NodeId, u64)>,
    /// The axis's point (x, y), and the input that is a drive (`carried`).
    centre: [f64; 2],
    drive: Option<Var>,
    /// The phase read periodically (its node) and its coupling.
    coupling: Option<(NodeId, Coupling)>,
}

impl Analysis<'_> {
    fn new(g: &Graph, centre: [f64; 2], drive: Option<Var>) -> Analysis<'_> {
        Analysis { g, class: HashMap::new(), slack: 0.0, folds: Vec::new(), centre, drive, coupling: None }
    }

    fn classify(&mut self, id: NodeId) -> Class {
        let node = self.g.nodes[id as usize];
        let form = |i: usize| Class::Lin(Lin { a: (i == 0).then(Coef::one), b: (i == 1).then(Coef::one), order: Order::Every });
        match node {
            Node::Var(Var::X) => return if self.centre[0] == 0.0 { form(0) } else { Class::Raw },
            Node::Var(Var::Y) => return if self.centre[1] == 0.0 { form(1) } else { Class::Raw },
            Node::Var(v) if Some(v) == self.drive => return Class::Drive,
            Node::Const(_) | Node::Var(_) | Node::Frame(_) => return Class::Inv(Order::Every),
            Node::Tag(..) | Node::Lipschitz(..) | Node::Deriv(..) | Node::Read(..) | Node::Real(..) => {
                return self.class[&value_children(node).0[0]].clone();
            }
            _ => {}
        }
        // x - cx: the form about the centre.
        if let Some((i, c)) = about(self.g, id)
            && c == self.centre[i]
        {
            return form(i);
        }
        // A function of invariant values; of them and drive values, a drive
        // value (of values every turn keeps: a phase reads it whole).
        let (ch, k) = value_children(node);
        let (mut all, mut drive) = (Some(Order::Every), false);
        for c in &ch[..k] {
            all = match (all, &self.class[c]) {
                (Some(o), Class::Inv(p)) => Some(o.gcd(*p)),
                (Some(o), Class::Drive) => {
                    drive = true;
                    Some(o)
                }
                _ => None,
            };
        }
        match all {
            Some(Order::Every) if drive => return Class::Drive,
            Some(_) if drive => return Class::Inv(NONE),
            Some(o) => return Class::Inv(o),
            None => {}
        }
        let class = |n: NodeId| self.class[&n].clone();
        match node {
            Node::Un(op, a) => match (op, class(a)) {
                (UnOp::Neg, Class::Lin(l)) => Class::Lin(l.negated()),
                (UnOp::Neg, Class::Squares(o, t)) => Class::Squares(o, t.into_iter().map(|(s, l)| (!s, l)).collect()),
                (UnOp::Neg, Class::Phase(p)) => Class::Phase(p.scaled(I::point(-1.0))),
                (UnOp::Square, Class::Lin(l)) => Class::Squares(Order::Every, vec![(false, l)]),
                (UnOp::Cos | UnOp::Sin, Class::Phase(p)) => self.periodic(a, p),
                // An exact scaling, one factor for every form it scales.
                (UnOp::Rad | UnOp::Deg, Class::Lin(l)) => {
                    let f = Factor::Pi180(op == UnOp::Deg);
                    let times = |c: Option<Coef>| {
                        c.map(|mut c| {
                            c.factors.push(f);
                            c.factors.sort_unstable();
                            c
                        })
                    };
                    Class::Lin(Lin { a: times(l.a), b: times(l.b), order: l.order })
                }
                (UnOp::Rad, Class::Phase(p)) => Class::Phase(p.scaled(interval::RAD)),
                (UnOp::Deg, Class::Phase(p)) => Class::Phase(p.scaled(interval::DEG)),
                _ => Class::Inv(NONE),
            },
            Node::Bin(op, a, b) => self.bin(op, a, b),
            _ => Class::Inv(NONE),
        }
    }

    fn bin(&mut self, op: BinOp, a: NodeId, b: NodeId) -> Class {
        let (ca, cb) = (self.class[&a].clone(), self.class[&b].clone());
        let constant = |n: NodeId| self.g.as_const(n).map(I::point);
        match op {
            BinOp::Add | BinOp::Sub => {
                if op == BinOp::Sub
                    && let Some(c) = self.fold(a, b)
                {
                    return c;
                }
                // A phase plus a drive value.
                let sign = I::point(if op == BinOp::Sub { -1.0 } else { 1.0 });
                match (ca, cb) {
                    (Class::Phase(p), Class::Drive) => Class::Phase(p.plus(b, sign)),
                    (Class::Drive, Class::Phase(q)) => Class::Phase(q.scaled(sign).plus(a, I::point(1.0))),
                    (ca, cb) => add(ca, if op == BinOp::Sub { negate(cb) } else { cb }),
                }
            }
            BinOp::Mul => match (ca, cb) {
                (Class::Lin(l), Class::Lin(_)) if a == b => Class::Squares(Order::Every, vec![(false, l)]),
                (Class::Lin(l), Class::Inv(o)) => self.scaled(l, b, o, false),
                (Class::Inv(o), Class::Lin(l)) => self.scaled(l, a, o, false),
                (Class::Phase(p), Class::Inv(_)) if constant(b).is_some() => Class::Phase(p.scaled(constant(b).unwrap())),
                (Class::Inv(_), Class::Phase(p)) if constant(a).is_some() => Class::Phase(p.scaled(constant(a).unwrap())),
                _ => Class::Inv(NONE),
            },
            BinOp::Div => match (ca, cb) {
                (Class::Lin(l), Class::Inv(o)) => self.scaled(l, b, o, true),
                (Class::Phase(p), Class::Inv(_)) => match constant(b) {
                    Some(c) if c.lo != 0.0 => Class::Phase(p.scaled(interval::bin(BinOp::Div, I::point(1.0), c).0)),
                    _ => Class::Inv(NONE),
                },
                _ => Class::Inv(NONE),
            },
            // atan2 v u of forms that turn together: an angle, moved by the
            // turn (k = 1) or against it (-1).
            BinOp::Atan2 => match (ca, cb) {
                (Class::Lin(v), Class::Lin(u)) if u.pairs(&v) => {
                    let with = v.a == neg_coef(&u.b) && v.b == u.a;
                    let fixed = [&u.a, &u.b]
                        .into_iter()
                        .flatten()
                        .all(|c| c.factors.iter().all(|f| matches!(f, Factor::Num(x, _) if f64::from_bits(*x) != 0.0) || matches!(f, Factor::Pi180(_))));
                    Class::Phase(Phase { k: I::point(if with { 1.0 } else { -1.0 }), drives: Vec::new(), order: u.order.gcd(v.order), fixed })
                }
                _ => Class::Inv(NONE),
            },
            _ => Class::Inv(NONE),
        }
    }

    /// cos or sin of phase `arg` (`p`): a function of the phase alone, where
    /// it reads a drive value, its angle's coefficient is a whole number
    /// but for 1e-9 and its invariant part keeps every turn; the one phase
    /// so read (`Coupling`).
    fn periodic(&mut self, arg: NodeId, p: Phase) -> Class {
        let n = p.k.mid().round();
        let off = (p.k.lo - n).abs().max((p.k.hi - n).abs());
        let finite = |x: I| x.lo.is_finite() && x.hi.is_finite();
        if p.drives.is_empty() || !p.fixed || p.order != Order::Every || n == 0.0 || !(off <= 1e-9) || !p.drives.iter().all(|d| finite(d.1)) {
            return Class::Inv(NONE);
        }
        match &self.coupling {
            None => {
                let seam = interval::mul_up(std::f64::consts::TAU.next_up(), off);
                self.coupling = Some((arg, Coupling { k: p.k, drives: p.drives, seam }));
                Class::Inv(Order::Every)
            }
            Some((a, _)) if *a == arg => Class::Inv(Order::Every),
            _ => Class::Inv(NONE),
        }
    }

    /// l times (or divided by) the invariant k.
    fn scaled(&self, l: Lin, k: NodeId, order: Order, recip: bool) -> Class {
        let mut neg = false;
        let mut fs = Vec::new();
        self.factors(k, recip, &mut neg, &mut fs);
        let scale = |c: Option<Coef>| {
            c.map(|mut c| {
                c.neg ^= neg;
                c.factors.extend(&fs);
                c.factors.sort_unstable();
                c
            })
        };
        Class::Lin(Lin { a: scale(l.a), b: scale(l.b), order: l.order.gcd(order) })
    }

    /// k as a sign and factors: constants by magnitude, negations and
    /// products of real numbers taken apart (not under a reciprocal: 1 /
    /// (a b) is written as one factor). A product of constants and params
    /// alone is its f64 value (the core's value rules), which reordering
    /// its factors can change (one order underflows to 0, another not):
    /// it is one factor, equal only to itself.
    fn factors(&self, k: NodeId, recip: bool, neg: &mut bool, out: &mut Vec<Factor>) {
        match self.g.nodes[k as usize] {
            Node::Const(v) => {
                *neg ^= v < 0.0;
                if v.abs() != 1.0 {
                    out.push(Factor::Num(v.abs().to_bits(), recip));
                }
            }
            Node::Un(UnOp::Neg, a) => {
                *neg ^= true;
                self.factors(a, recip, neg, out);
            }
            Node::Bin(BinOp::Mul, a, b) if !recip && !fixed(self.g, k) => {
                self.factors(a, false, neg, out);
                self.factors(b, false, neg, out);
            }
            _ => out.push(Factor::Node(k, recip)),
        }
    }

    /// `a - m`, where m is P floor(a / P): mod(a, P), of an angle a and a
    /// constant P within 1e-9 of 2 pi / n, is a function of a periodic
    /// in P, so of order n (with the slack of the module's contract).
    fn fold(&mut self, a: NodeId, m: NodeId) -> Option<Class> {
        let o = match &self.class[&a] {
            Class::Phase(p) if p.angle() => p.order,
            _ => return None,
        };
        let Node::Bin(BinOp::Mul, m1, m2) = self.g.nodes[m as usize] else { return None };
        let n = &self.g.nodes;
        let pitch = [(m1, m2), (m2, m1)].into_iter().find_map(|(p, f)| {
            let Node::Un(UnOp::Floor, d) = n[f as usize] else { return None };
            matches!(n[d as usize], Node::Bin(BinOp::Div, a2, p2) if a2 == a && p2 == p).then_some(())?;
            self.g.as_const(p)
        })?;
        let p = pitch.abs();
        let turns = (std::f64::consts::TAU / p).round();
        if !(turns >= 1.0 && turns < 2f64.powi(52)) {
            return None;
        }
        // |n P - 2 pi|, rounded up: n P rounds by at most half an ulp of
        // 2 pi, TAU is within half of one, and their difference is exact.
        let off = (turns * p - std::f64::consts::TAU).abs() + 2.0 * f64::EPSILON * std::f64::consts::TAU;
        if off > 1e-9 {
            return None;
        }
        self.slack = self.slack.max(off);
        if !self.folds.contains(&(a, p.to_bits())) {
            self.folds.push((a, p.to_bits()));
        }
        Some(Class::Inv(Order::Of(turns as u64).gcd(o)))
    }
}

fn negate(c: Class) -> Class {
    match c {
        Class::Lin(l) => Class::Lin(l.negated()),
        Class::Squares(o, t) => Class::Squares(o, t.into_iter().map(|(s, l)| (!s, l)).collect()),
        Class::Phase(p) => Class::Phase(p.scaled(I::point(-1.0))),
        c => c, // an invariant's negation is one
    }
}

/// The sum of two classes (a difference: the second negated).
fn add(a: Class, b: Class) -> Class {
    match (a, b) {
        (Class::Lin(p), Class::Lin(q)) => {
            // Coefficients are not added: one side's, where the other has none.
            let slot = |x: Option<Coef>, y: Option<Coef>| match (x, y) {
                (Some(_), Some(_)) => Err(()),
                (x, y) => Ok(x.or(y)),
            };
            match (slot(p.a, q.a), slot(p.b, q.b)) {
                (Ok(a), Ok(b)) => Class::Lin(Lin { a, b, order: p.order.gcd(q.order) }),
                _ => Class::Inv(NONE),
            }
        }
        (Class::Phase(q), Class::Inv(p)) | (Class::Inv(p), Class::Phase(q)) => Class::Phase(Phase { order: q.order.gcd(p), ..q }),
        (Class::Squares(o, t), Class::Inv(p)) | (Class::Inv(p), Class::Squares(o, t)) => squares(o.gcd(p), t),
        (Class::Squares(o, mut t), Class::Squares(p, u)) => {
            t.extend(u);
            squares(o.gcd(p), t)
        }
        _ => Class::Inv(NONE),
    }
}

/// Pair up terms of one sign whose forms turn together: each pair is
/// invariant under the turns its coefficients keep.
fn squares(mut order: Order, terms: Vec<(bool, Lin)>) -> Class {
    let mut left: Vec<(bool, Lin)> = Vec::new();
    for (s, l) in terms {
        match left.iter().position(|(s2, l2)| *s2 == s && l2.pairs(&l)) {
            Some(i) => {
                let (_, l2) = left.swap_remove(i);
                order = order.gcd(l.order).gcd(l2.order);
            }
            None => left.push((s, l)),
        }
    }
    if left.is_empty() { Class::Inv(order) } else { Class::Squares(order, left) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, TAU};

    fn xy(g: &mut Graph) -> (NodeId, NodeId) {
        (g.var(Var::X), g.var(Var::Y))
    }

    fn len2(g: &mut Graph, u: NodeId, v: NodeId) -> NodeId {
        let (u2, v2) = (g.un(UnOp::Square, u), g.un(UnOp::Square, v));
        let s = g.add(u2, v2);
        g.un(UnOp::Sqrt, s)
    }

    /// std's (mod a m).
    fn modulo(g: &mut Graph, a: NodeId, m: NodeId) -> NodeId {
        let d = g.bin(BinOp::Div, a, m);
        let f = g.un(UnOp::Floor, d);
        let k = g.mul(m, f);
        g.sub(a, k)
    }

    /// std's repeat-polar of f(x, y, z): f at the point turned into the
    /// first of n wedges.
    fn repeat_polar(g: &mut Graph, f: NodeId, n: f64) -> NodeId {
        let (x, y) = xy(g);
        let p = g.constant(2.0 * PI / n);
        let half = g.constant(0.5 * (2.0 * PI / n));
        let at = g.bin(BinOp::Atan2, y, x);
        let shifted = g.add(at, half);
        let m = modulo(g, shifted, p);
        let a = g.sub(m, half);
        let r = len2(g, x, y);
        let (c, s) = (g.un(UnOp::Cos, a), g.un(UnOp::Sin, a));
        let (rx, ry) = (g.mul(r, c), g.mul(r, s));
        g.subst(f, &[(Var::X, rx), (Var::Y, ry)])
    }

    /// A box 2 x 1 centred at (5, 0): no symmetry about z of its own.
    fn tooth(g: &mut Graph) -> NodeId {
        let (x, y) = xy(g);
        let five = g.constant(5.0);
        let dx = g.sub(x, five);
        let (ax, ay) = (g.un(UnOp::Abs, dx), g.un(UnOp::Abs, y));
        let (one, half) = (g.constant(1.0), g.constant(0.5));
        let (qx, qy) = (g.sub(ax, one), g.sub(ay, half));
        g.max(qx, qy)
    }

    fn order(g: &Graph, f: NodeId) -> Order {
        turn_symmetry(g, f).order
    }

    #[test]
    fn round_fields_keep_every_turn() {
        let mut g = Graph::new();
        let (x, y) = xy(&mut g);
        let z = g.var(Var::Z);
        let r = len2(&mut g, x, y);
        let three = g.constant(3.0);
        let dr = g.sub(r, three);
        let dz = g.un(UnOp::Abs, z);
        let cyl = g.bin(BinOp::Emax, dr, dz);
        assert_eq!(order(&g, cyl), Order::Every);
        // Turned about z by 30 degrees (std's rotate-z) and mirrored, still round.
        let (c, s) = (g.constant((PI / 6.0).cos()), g.constant((PI / 6.0).sin()));
        let (cx, sy, cy, sx) = (g.mul(c, x), g.mul(s, y), g.mul(c, y), g.mul(s, x));
        let (u, v) = (g.add(cx, sy), g.sub(cy, sx));
        let nu = g.neg(u);
        let turned = g.subst(cyl, &[(Var::X, nu), (Var::Y, v)]);
        assert_eq!(order(&g, turned), Order::Every);
        // Moved off the axis, or squashed: not.
        let moved = g.sub(x, three);
        let off = g.subst(cyl, &[(Var::X, moved)]);
        assert_eq!(order(&g, off), Order::Of(1));
        let squashed = g.mul(three, x);
        let oval = g.subst(cyl, &[(Var::X, squashed)]);
        assert_eq!(order(&g, oval), Order::Of(1));
        // x or y alone, and a D-flat (the cylinder cut by x < 2).
        assert_eq!(order(&g, x), Order::Of(1));
        let two = g.constant(2.0);
        let flat = g.sub(x, two);
        let d = g.max(cyl, flat);
        assert_eq!(order(&g, d), Order::Of(1));
    }

    #[test]
    fn a_fixed_product_is_its_float_not_its_factors() {
        // u = x ((p a) b) c and v = y (p b) a c: the same factors, but
        // (p a) b and p b are floats of params and constants, and at p =
        // 1e-200 the first underflows to 0 while the second is 1e100. The
        // strip sqrt(u^2 + v^2) is then |v|: not round.
        let mut g = Graph::new();
        let (x, y) = xy(&mut g);
        let p = g.var(Var::Param(0));
        let (a, b, c) = (g.constant(1e-200), g.constant(1e300), g.constant(1e100));
        let pa = g.mul(p, a);
        let pab = g.mul(pa, b);
        let xk = g.mul(x, pab);
        let u = g.mul(xk, c);
        let pb = g.mul(p, b);
        let yk = g.mul(y, pb);
        let yka = g.mul(yk, a);
        let v = g.mul(yka, c);
        let r = len2(&mut g, u, v);
        assert_eq!(order(&g, r), Order::Of(1));
        // Products of real numbers (a form times each factor in turn) may
        // be reordered, and one fixed product is equal to itself.
        let (xp, yp) = (g.mul(x, p), g.mul(y, p));
        let (xpa, ypb) = (g.mul(xp, a), g.mul(yp, b));
        let (u, v) = (g.mul(xpa, b), g.mul(ypb, a));
        let r = len2(&mut g, u, v);
        assert_eq!(order(&g, r), Order::Every);
        let (u, v) = (g.mul(x, pab), g.mul(pab, y));
        let r = len2(&mut g, u, v);
        assert_eq!(order(&g, r), Order::Every);
    }

    #[test]
    fn folds_give_their_count_and_unions_the_gcd() {
        let mut g = Graph::new();
        let t = tooth(&mut g);
        assert_eq!(order(&g, t), Order::Of(1));
        for n in [1.0, 2.0, 7.0, 12.0, 45.0, 360.0] {
            let f = repeat_polar(&mut g, t, n);
            let s = turn_symmetry(&g, f);
            assert_eq!(s.order, Order::Of(n as u64), "{n}");
            if n > 1.0 {
                // At least the real |n P - 2 pi| of the float pitch, and tiny.
                let p = 2.0 * PI / n;
                assert!(s.slack >= (n * p - TAU).abs() && s.slack < 1e-13, "{n}: {}", s.slack);
            }
        }
        let (twelve, eight) = (repeat_polar(&mut g, t, 12.0), repeat_polar(&mut g, t, 8.0));
        // Two folds: each is off by its own error after a turn, which no
        // one rotation of the point bounds.
        let both = g.min(twelve, eight);
        assert_eq!(order(&g, both), Order::Of(1));
        // A fold whose pitch closes no circle (2 pi / 7.5) shows nothing.
        let odd = repeat_polar(&mut g, t, 7.5);
        assert_eq!(order(&g, odd), Order::Of(1));
        // The angle itself turns with the point, and so does the fold's
        // cell index (which tooth): only the fold is invariant.
        let (x, y) = xy(&mut g);
        let at = g.bin(BinOp::Atan2, y, x);
        let p = g.constant(2.0 * PI / 12.0);
        let cell = g.bin(BinOp::Div, at, p);
        let cell = g.un(UnOp::Floor, cell);
        let folded = modulo(&mut g, at, p);
        assert_eq!([at, cell, folded].map(|f| order(&g, f)), [Order::Of(1), Order::Of(1), Order::Of(12)]);
    }

    /// A ring of n teeth of pitch 2 about (cx, 0) whose stations run with
    /// drive value `s` (`speed` mm of station per mm of s), as std's
    /// path-repeat-at places a belt's teeth round a wrap: each a band of
    /// the ring where cos of the station's phase is near 1, cut by the
    /// wedge x < cx (the wrap's half). The cut ring and its teeth.
    fn ring(g: &mut Graph, n: f64, cx: f64, speed: f64, s: NodeId) -> (NodeId, NodeId) {
        let (x, y) = xy(g);
        let c = g.constant(cx);
        let qx = g.sub(x, c);
        let a = g.bin(BinOp::Atan2, y, qx);
        let r = g.constant(n * 2.0 / TAU);
        let station = g.mul(r, a);
        let k = g.constant(speed);
        let moved = g.mul(k, s);
        let d = g.sub(station, moved);
        let w = g.constant(PI);
        let ph = g.mul(w, d);
        let cos = g.un(UnOp::Cos, ph);
        let rho = len2(g, qx, y);
        let dr = g.sub(rho, r);
        let ar = g.un(UnOp::Abs, dr);
        let half = g.constant(0.5);
        let radial = g.sub(ar, half);
        let nine = g.constant(0.9);
        let ac = g.un(UnOp::Neg, cos);
        let along = g.add(ac, nine);
        let tooth = g.max(radial, along);
        (g.max(tooth, qx), tooth)
    }

    /// The teeth of a ring that runs with its drive are carried by a turn
    /// about its centre with the drive moved at the rate their phase gives
    /// (the pitch radius, mm per radian); about another
    /// centre, at a pitch that closes no circle, still, or with no drive,
    /// nothing is.
    #[test]
    fn a_ring_of_teeth_running_with_its_drive_is_carried() {
        let mut g = Graph::new();
        let drive = Some(Var::Param(0));
        let u = g.var(Var::Param(0));
        let three = g.constant(3.0);
        let s = g.mul(three, u);
        let r = 40.0 / TAU;
        for cx in [0.0, 430.0] {
            let (cut, tooth) = ring(&mut g, 20.0, cx, 1.0, s);
            let c = carried(&g, cut, [cx, 0.0], drive);
            let k = c.coupling.expect("carried");
            // Phase 20 a - pi s: s moves r per radian of turn.
            let ([(d, cs)], n) = (&k.drives[..], 20.0) else { panic!("{:?}", k.drives) };
            let rate = -k.k.mid() / cs.mid();
            assert!(g.depends_on(*d, |v| v == Var::Param(0)) && k.k.lo <= n && n <= k.k.hi && (rate - r).abs() < 1e-12, "{k:?}");
            assert!(k.seam > 0.0 && k.seam < 1e-12, "{}", k.seam);
            // Turned by b about the centre with u moved by rate b / 3: the
            // same, in floats to their rounding.
            let tape = crate::tape::Tape::new(&g, tooth);
            let at = |u: f64| crate::eval::Ctx { t: 0.0, params: vec![u], t_span: None, p_spans: Default::default() };
            for i in 0..200 {
                let (b, u0) = (0.37 * i as f64 - 30.0, 0.11 * i as f64);
                let p = [cx + r * (0.05 * i as f64).cos(), r * (0.05 * i as f64).sin(), 0.0];
                let q = [cx + (p[0] - cx) * b.cos() - p[1] * b.sin(), (p[0] - cx) * b.sin() + p[1] * b.cos(), 0.0];
                let (f, h) = (crate::eval::eval_point(&tape, &at(u0), p), crate::eval::eval_point(&tape, &at(u0 + rate * b / 3.0), q));
                assert!((f - h).abs() < 1e-9, "{b} {u0}: {f} {h}");
            }
            assert!(carried(&g, cut, [cx + 1.0, 0.0], drive).coupling.is_none());
            assert!(carried(&g, cut, [cx, 0.0], None).invariant.iter().all(|&n| n != tooth));
        }
        let (odd, _) = ring(&mut g, 20.5, 0.0, 1.0, s);
        assert!(carried(&g, odd, [0.0, 0.0], drive).coupling.is_none());
        let zero = g.constant(0.0);
        let (still, _) = ring(&mut g, 20.0, 0.0, 1.0, zero);
        assert!(carried(&g, still, [0.0, 0.0], drive).coupling.is_none());
        // cos(20 atan2(A y, A x) - s): carried for a constant A; not for A a
        // value (a param), 0 somewhere, where atan2 is 0 and does not turn.
        let (x, y) = xy(&mut g);
        let teeth = |g: &mut Graph, a: NodeId| {
            let (ax, ay) = (g.mul(a, x), g.mul(a, y));
            let angle = g.bin(BinOp::Atan2, ay, ax);
            let twenty = g.constant(20.0);
            let k = g.mul(twenty, angle);
            let ph = g.sub(k, s);
            g.un(UnOp::Cos, ph)
        };
        let (two, param) = (g.constant(2.0), g.var(Var::Param(1)));
        let (by_two, by_param) = (teeth(&mut g, two), teeth(&mut g, param));
        assert!(carried(&g, by_two, [0.0, 0.0], drive).coupling.is_some() && carried(&g, by_param, [0.0, 0.0], drive).coupling.is_none());
    }
}
