//! A tape's screw symmetry about z: the values a turn about z through the
//! origin by any angle a, with a move along z by L a / (2 pi) (a screw
//! motion of lead L), leaves unchanged, shown from how the tape is built:
//! never sampled. Two threads of one lead on one axis, in any phase, are
//! such values, so where a field is one it is the same at every point as
//! at a point of the half-plane y = 0, x > 0, at a height within one lead:
//! a question of the cross-section (`Screw`).
//!
//! `screw(tape)` classifies each instruction, as a function of real
//! numbers (each its operation on reals, under the core's value rules;
//! each constant the float it is). An instruction is `invariant` where,
//! for every point p off the axis (x^2 + y^2 > 0), every angle a, and the
//! screw motion M_a of the tape's one lead L, its value at M_a p is its
//! value at p (at every t and param value). The rules, each an identity in
//! real numbers:
//! - constants, params, t, a promise's frame (value 0): every motion;
//! - x, y and c x, c y (c a constant): linear forms A x + B y + C, each
//!   coefficient one constant or none, C the sum of the constants added,
//!   bounded outward (a frame moved off and back, x - -10 - 10, is exactly
//!   x); two forms u = A x + B y and v = -B x + A y (or B x - A y), each C
//!   exactly 0, turn together: u^2 + v^2 is (A^2 + B^2)(x^2 + y^2), kept by
//!   every motion;
//! - atan2 v u of such a pair: an angle, which M_a moves by a or -a, and by
//!   2 pi more or less where it crosses atan2's cut (u < 0, v = 0);
//! - z, k z + c, scaled and negated: an axial value, which M_a moves by k
//!   L a / (2 pi);
//! - a phase: an angle times k_a plus an axial value k_z z, plus values the
//!   motions keep: M_a moves it by (k_a + k_z L / (2 pi)) a, 0 where L is
//!   -2 pi k_a / k_z, and across the cut by J = 2 pi k_a; the coefficients
//!   are kept as written, products of constants and of pi or its reciprocal
//!   (UnOp::Deg's 180 / pi), so two phases have one lead only where their
//!   leads, so written, are the same product;
//! - q - floor(q) of a phase q whose J is a whole number (shown exactly):
//!   kept by every motion (q moves by whole numbers only);
//! - any operation of values every motion keeps.
//!
//! Anything else (a form plus a constant, a move off the axis; a phase read
//! other than through its fractional part; a second lead; a sampled
//! field, a derivative) is kept by none. The field then has period L along
//! z (M_{2 pi} is that move), so its value at (r cos b, r sin b, z) is its
//! value at (r, 0, z - L b / (2 pi) + n L) for every whole n.

use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// A factor of a coefficient as written: a constant's magnitude (its
/// bits) or its reciprocal (`true`), or pi or its reciprocal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Factor {
    Num(u64, bool),
    Pi(bool),
}

impl Factor {
    fn recip(self) -> Factor {
        match self {
            Factor::Num(x, r) => Factor::Num(x, !r),
            Factor::Pi(r) => Factor::Pi(!r),
        }
    }
}

/// A real coefficient: its sign times the product of its factors, sorted,
/// a factor and its reciprocal cancelled (so equal reals written alike
/// compare equal).
#[derive(Clone, Debug, PartialEq, Eq)]
struct Coef {
    neg: bool,
    factors: Vec<Factor>,
}

impl Coef {
    fn one() -> Coef {
        Coef { neg: false, factors: Vec::new() }
    }
    fn of(c: f64) -> Option<Coef> {
        (c != 0.0 && c.is_finite()).then(|| Coef { neg: c < 0.0, factors: vec![Factor::Num(c.abs().to_bits(), false)] })
    }
    fn negated(mut self) -> Coef {
        self.neg = !self.neg;
        self
    }
    fn times(mut self, f: Factor) -> Coef {
        match self.factors.iter().position(|&g| g == f.recip()) {
            Some(i) => {
                self.factors.remove(i);
            }
            None => {
                self.factors.push(f);
                self.factors.sort_unstable();
            }
        }
        self
    }
    fn mul(self, o: &Coef) -> Coef {
        let neg = self.neg != o.neg;
        let mut c = o.factors.iter().fold(self, |c, &f| c.times(f));
        c.neg = neg;
        c
    }
    fn recip(self) -> Coef {
        Coef { neg: self.neg, factors: self.factors.into_iter().map(Factor::recip).collect::<Vec<_>>() }.sorted()
    }
    fn sorted(mut self) -> Coef {
        self.factors.sort_unstable();
        self
    }
    /// Its value, rounded outward (exact where each step is); None where
    /// a factor of pi is left (no float is it).
    fn value(&self) -> Option<I> {
        let mut v = I::point(if self.neg { -1.0 } else { 1.0 });
        for f in &self.factors {
            v = match *f {
                Factor::Num(x, false) => interval::mul_o(v, I::point(f64::from_bits(x))),
                Factor::Num(x, true) => interval::div_o(v, I::point(f64::from_bits(x))),
                Factor::Pi(_) => return None,
            };
        }
        Some(v)
    }
    /// Its value rounded outward, a factor of pi taken as pi's bounds.
    fn bound(&self) -> I {
        let pi = I::new(std::f64::consts::PI, std::f64::consts::PI.next_up());
        let mut v = I::point(if self.neg { -1.0 } else { 1.0 });
        for f in &self.factors {
            v = match *f {
                Factor::Num(x, false) => interval::mul_o(v, I::point(f64::from_bits(x))),
                Factor::Num(x, true) => interval::div_o(v, I::point(f64::from_bits(x))),
                Factor::Pi(false) => interval::mul_o(v, pi),
                Factor::Pi(true) => interval::div_o(v, pi),
            };
        }
        v
    }
}

/// A x + B y + C, each coefficient one constant (its bits) or none (0), C
/// bounded.
#[derive(Clone, Copy, Debug)]
struct Lin {
    a: Option<u64>,
    b: Option<u64>,
    c: I,
}

fn neg_bits(c: Option<u64>) -> Option<u64> {
    c.map(|x| (-f64::from_bits(x)).to_bits())
}

impl Lin {
    fn of(a: Option<f64>, b: Option<f64>) -> Lin {
        Lin { a: a.map(f64::to_bits), b: b.map(f64::to_bits), c: I::point(0.0) }
    }
    fn negated(self) -> Lin {
        Lin { a: neg_bits(self.a), b: neg_bits(self.b), c: I::new(-self.c.hi, -self.c.lo) }
    }
    /// Through the axis: C exactly 0.
    fn about_axis(&self) -> bool {
        self.c.lo == 0.0 && self.c.hi == 0.0
    }
    /// Whether (self, v), both through the axis, turn together (v = -B x +
    /// A y) or against (v = B x - A y): Some(1) or Some(-1).
    fn pairs(&self, v: &Lin) -> Option<f64> {
        if self.a.is_none() && self.b.is_none() || !self.about_axis() || !v.about_axis() {
            return None;
        }
        if v.a == neg_bits(self.b) && v.b == self.a {
            Some(1.0)
        } else if v.a == self.b && v.b == neg_bits(self.a) {
            Some(-1.0)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug)]
enum Class {
    /// Kept by every motion of the lead.
    Inv,
    Lin(Lin),
    /// A kept value plus sign * lin^2 for each term.
    Squares(Vec<(bool, Lin)>),
    /// k z plus a kept value.
    Axial(Coef),
    /// k_a angle (+ k_z z) plus a kept value.
    Phase(Coef, Option<Coef>),
    /// floor of the phase at this instruction.
    Floor(usize),
    /// Kept by no motion shown.
    No,
}

/// What `screw` shows of a tape (the module's contract).
#[derive(Clone, Debug)]
pub struct Screw {
    /// The lead, bounded outward.
    pub lead: I,
    invariant: Vec<bool>,
}

impl Screw {
    /// Whether instruction `i` (an ssa index) is kept by every screw motion
    /// of the lead (the module's contract).
    pub fn invariant(&self, i: usize) -> bool {
        self.invariant[i]
    }
}

/// The tape's screw symmetry about z (the module's contract): None where
/// it reads no phase of one lead (nothing helical is shown).
pub fn screw(tape: &Tape) -> Option<Screw> {
    let ssa = &tape.ssa;
    let mut class: Vec<Class> = Vec::with_capacity(ssa.len());
    let mut lead: Option<Coef> = None;
    let constant = |i: u32| match ssa[i as usize].op {
        Op::Const(c) => Some(c),
        _ => None,
    };
    for ins in ssa {
        let (args, k) = ins.args();
        let of = |n: usize| &class[args[n] as usize];
        let c = match ins.op {
            Op::Const(_) | Op::Var(Var::T | Var::Param(_)) | Op::Frame => Class::Inv,
            Op::Var(Var::X) => Class::Lin(Lin::of(Some(1.0), None)),
            Op::Var(Var::Y) => Class::Lin(Lin::of(None, Some(1.0))),
            Op::Var(Var::Z) => Class::Axial(Coef::one()),
            Op::Tag(_) | Op::Lipschitz => of(0).clone(),
            Op::Grid(..) | Op::Read | Op::Deriv(_) | Op::Real | Op::Within(_) => Class::No,
            _ if args[..k].iter().all(|&a| matches!(class[a as usize], Class::Inv)) => Class::Inv,
            Op::Un(op) => match (op, of(0).clone()) {
                (UnOp::Neg, Class::Lin(l)) => Class::Lin(l.negated()),
                (UnOp::Neg, Class::Squares(t)) => Class::Squares(t.into_iter().map(|(s, l)| (!s, l)).collect()),
                (UnOp::Neg, Class::Axial(c)) => Class::Axial(c.negated()),
                (UnOp::Neg, Class::Phase(a, z)) => Class::Phase(a.negated(), z.map(Coef::negated)),
                (UnOp::Square, Class::Lin(l)) => Class::Squares(vec![(false, l)]),
                (UnOp::Floor, Class::Phase(..)) => Class::Floor(ins.a as usize),
                (UnOp::Deg | UnOp::Rad, cl @ (Class::Axial(_) | Class::Phase(..))) => {
                    // 180 / pi, or pi / 180, exactly.
                    let s = |c: Coef| {
                        if op == UnOp::Deg {
                            c.times(Factor::Num(180f64.to_bits(), false)).times(Factor::Pi(true))
                        } else {
                            c.times(Factor::Pi(false)).times(Factor::Num(180f64.to_bits(), true))
                        }
                    };
                    match cl {
                        Class::Axial(c) => Class::Axial(s(c)),
                        Class::Phase(a, z) => Class::Phase(s(a), z.map(s)),
                        _ => unreachable!(),
                    }
                }
                _ => Class::No,
            },
            Op::Bin(op) => {
                let (ca, cb) = (of(0).clone(), of(1).clone());
                let sub = op == BinOp::Sub;
                match op {
                    // The fractional part of a phase whose jump at the cut
                    // is whole: kept, the lead now required.
                    BinOp::Sub if matches!(cb, Class::Floor(q) if q == ins.a as usize) => match ca {
                        Class::Phase(ka, Some(kz)) if whole(&ka) => {
                            // L = -2 pi k_a / k_z.
                            let l = ka.mul(&kz.recip()).times(Factor::Num(2f64.to_bits(), false)).times(Factor::Pi(false)).negated();
                            match &lead {
                                Some(m) if *m != l => Class::No,
                                _ => {
                                    lead = Some(l);
                                    Class::Inv
                                }
                            }
                        }
                        _ => Class::No,
                    },
                    BinOp::Add | BinOp::Sub => {
                        let flip = |c: Class| if sub { negate(c) } else { c };
                        // A form plus a constant: its C moved.
                        let moved = |l: Lin, k: u32, s: f64| match constant(k) {
                            Some(c) => Class::Lin(Lin { c: interval::add_o(l.c, I::point(s * c)), ..l }),
                            None => Class::No,
                        };
                        match (ca.clone(), cb.clone()) {
                            (Class::Lin(l), Class::Inv) => return_class(moved(l, ins.b, if sub { -1.0 } else { 1.0 })),
                            (Class::Inv, Class::Lin(l)) => return_class(moved(if sub { l.negated() } else { l }, ins.a, 1.0)),
                            _ => None,
                        }
                        .unwrap_or_else(|| match (ca, flip(cb)) {
                            (Class::Inv, c) | (c @ (Class::Axial(_) | Class::Phase(..) | Class::Squares(_)), Class::Inv) => match c {
                                c @ (Class::Axial(_) | Class::Phase(..) | Class::Squares(_)) => c,
                                _ => Class::No,
                            },
                            (Class::Phase(a, None), Class::Axial(z)) | (Class::Axial(z), Class::Phase(a, None)) => Class::Phase(a, Some(z)),
                            (Class::Lin(u), Class::Lin(v)) => match (u.a.xor(v.a), u.b.xor(v.b)) {
                                (a, b) if (u.a.is_none() || v.a.is_none()) && (u.b.is_none() || v.b.is_none()) => Class::Lin(Lin { a, b, c: interval::add_o(u.c, v.c) }),
                                _ => Class::No,
                            },
                            (Class::Squares(mut s), Class::Squares(t)) => {
                                s.extend(t);
                                squares(s)
                            }
                            _ => Class::No,
                        })
                    }
                    BinOp::Mul => {
                        let by = |cl: Class, c: Option<f64>| -> Class {
                            match (cl, c.and_then(Coef::of)) {
                                (Class::Axial(k), Some(c)) => Class::Axial(k.mul(&c)),
                                (Class::Phase(a, z), Some(c)) => Class::Phase(a.mul(&c), z.map(|z| z.mul(&c))),
                                (Class::Lin(l), Some(_)) if l.about_axis() => {
                                    let c = c.unwrap();
                                    let one = Some(1f64.to_bits());
                                    match (l.a, l.b) {
                                        (a, None) if a == one => Class::Lin(Lin::of(Some(c), None)),
                                        (None, b) if b == one => Class::Lin(Lin::of(None, Some(c))),
                                        _ => Class::No,
                                    }
                                }
                                _ => Class::No,
                            }
                        };
                        match (ca, cb) {
                            (Class::Lin(l), Class::Lin(_)) if ins.a == ins.b => Class::Squares(vec![(false, l)]),
                            (Class::Inv, c) => by(c, constant(ins.a)),
                            (c, Class::Inv) => by(c, constant(ins.b)),
                            _ => Class::No,
                        }
                    }
                    BinOp::Div => match (ca, constant(ins.b).and_then(Coef::of)) {
                        (Class::Axial(k), Some(c)) => Class::Axial(k.mul(&c.recip())),
                        (Class::Phase(a, z), Some(c)) => {
                            let r = c.recip();
                            Class::Phase(a.mul(&r), z.map(|z| z.mul(&r)))
                        }
                        _ => Class::No,
                    },
                    BinOp::Atan2 => match (ca, cb) {
                        (Class::Lin(v), Class::Lin(u)) => match u.pairs(&v) {
                            Some(s) => Class::Phase(if s > 0.0 { Coef::one() } else { Coef::one().negated() }, None),
                            None => Class::No,
                        },
                        _ => Class::No,
                    },
                    _ => Class::No,
                }
            }
            _ => Class::No,
        };
        class.push(c);
    }
    let lead = lead?.bound();
    (lead.lo > 0.0 || lead.hi < 0.0)
        .then(|| Screw { lead: I::new(lead.lo.abs().min(lead.hi.abs()), lead.lo.abs().max(lead.hi.abs())), invariant: class.iter().map(|c| matches!(c, Class::Inv)).collect() })
}

/// A class decided, or None to try the other rules.
fn return_class(c: Class) -> Option<Class> {
    match c {
        Class::No => None,
        c => Some(c),
    }
}

fn negate(c: Class) -> Class {
    match c {
        Class::Inv => Class::Inv,
        Class::Lin(l) => Class::Lin(l.negated()),
        Class::Squares(t) => Class::Squares(t.into_iter().map(|(s, l)| (!s, l)).collect()),
        Class::Axial(k) => Class::Axial(k.negated()),
        Class::Phase(a, z) => Class::Phase(a.negated(), z.map(Coef::negated)),
        _ => Class::No,
    }
}

/// A sum of signed squares: kept where it is two of one sign whose forms
/// turn together (u^2 + v^2), else still a sum.
fn squares(t: Vec<(bool, Lin)>) -> Class {
    match t[..] {
        [(s, u), (r, v)] if s == r && u.pairs(&v).is_some() => Class::Inv,
        _ => Class::Squares(t),
    }
}

/// Whether a phase's jump at atan2's cut, 2 pi k_a, is a whole number
/// other than 0, shown exactly.
fn whole(ka: &Coef) -> bool {
    let j = ka.clone().times(Factor::Num(2f64.to_bits(), false)).times(Factor::Pi(false));
    j.value().is_some_and(|v| v.lo == v.hi && v.lo != 0.0 && v.lo.fract() == 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Graph, NodeId};

    /// thread-iso's phase and fractional part, written as std writes it:
    /// deg(atan2(y', x')) / 360 - z' / p, then its distance from a crest's
    /// centre, about the frame (x', y', z').
    fn crest(g: &mut Graph, at: [NodeId; 3], p: f64) -> NodeId {
        let a = g.bin(BinOp::Atan2, at[1], at[0]);
        let d = g.un(UnOp::Deg, a);
        let c360 = g.constant(360.0);
        let turns = g.bin(BinOp::Div, d, c360);
        let pc = g.constant(p);
        let zz = g.bin(BinOp::Div, at[2], pc);
        let q = g.sub(turns, zz);
        let f = g.un(UnOp::Floor, q);
        let w = g.sub(q, f);
        let r2 = {
            let (x2, y2) = (g.un(UnOp::Square, at[0]), g.un(UnOp::Square, at[1]));
            g.add(x2, y2)
        };
        let r = g.un(UnOp::Sqrt, r2);
        g.add(w, r)
    }

    fn vars(g: &mut Graph) -> [NodeId; 3] {
        [g.var(Var::X), g.var(Var::Y), g.var(Var::Z)]
    }

    #[test]
    fn two_threads_of_one_lead_are_kept() {
        let mut g = Graph::new();
        let v = vars(&mut g);
        let a = crest(&mut g, v, 0.5);
        // The mate moved along z and turned half over about x: y and z
        // negated.
        let h = g.constant(5.5);
        let z = g.sub(v[2], h);
        let (ny, nz) = (g.neg(v[1]), g.neg(z));
        let b = crest(&mut g, [v[0], ny, nz], 0.5);
        let root = g.max(a, b);
        let t = Tape::new(&g, root);
        let s = screw(&t).expect("a lead");
        assert!(s.invariant(t.ssa.len() - 1));
        assert!(s.lead.lo <= 0.5 && 0.5 <= s.lead.hi && s.lead.hi - s.lead.lo < 1e-15);
    }

    #[test]
    fn two_leads_or_a_second_axis_are_not() {
        for case in 0..3 {
            let mut g = Graph::new();
            let v = vars(&mut g);
            let a = crest(&mut g, v, 0.5);
            let b = match case {
                // Another pitch.
                0 => crest(&mut g, v, 0.7),
                // A left-hand thread: the angle against the turn.
                1 => {
                    let ny = g.neg(v[1]);
                    crest(&mut g, [v[0], ny, v[2]], 0.5)
                }
                // Off the axis.
                _ => {
                    let c = g.constant(0.01);
                    let x = g.sub(v[0], c);
                    crest(&mut g, [x, v[1], v[2]], 0.5)
                }
            };
            let root = g.max(a, b);
            let t = Tape::new(&g, root);
            let kept = screw(&t).is_some_and(|s| s.invariant(t.ssa.len() - 1));
            assert!(!kept, "case {case}");
        }
    }

    /// A mate whose frame is moved off the axis and back (x - -10 - 10,
    /// as a tapped part placed at its frame writes it) is on the axis, its
    /// height moved too: kept.
    #[test]
    fn a_frame_moved_off_and_back_is_kept() {
        let mut g = Graph::new();
        let v = vars(&mut g);
        let a = crest(&mut g, v, 0.5);
        let (m10, p10, m18, p13) = (g.constant(-10.0), g.constant(10.0), g.constant(-18.0), g.constant(13.0));
        let at = [v[0], v[1]].map(|u| {
            let off = g.sub(u, m10);
            g.sub(off, p10)
        });
        let z = g.sub(v[2], m18);
        let z = g.sub(z, p13);
        let b = crest(&mut g, [at[0], at[1], z], 0.5);
        let root = g.max(a, b);
        let t = Tape::new(&g, root);
        assert!(screw(&t).is_some_and(|s| s.invariant(t.ssa.len() - 1)));
    }

    #[test]
    fn a_phase_read_whole_is_not() {
        // The phase itself jumps by 1 at atan2's cut: only its fractional
        // part is kept.
        let mut g = Graph::new();
        let v = vars(&mut g);
        let a = g.bin(BinOp::Atan2, v[1], v[0]);
        let d = g.un(UnOp::Deg, a);
        let c360 = g.constant(360.0);
        let turns = g.bin(BinOp::Div, d, c360);
        let pc = g.constant(0.5);
        let zz = g.bin(BinOp::Div, v[2], pc);
        let q = g.sub(turns, zz);
        let t = Tape::new(&g, q);
        assert!(screw(&t).is_none_or(|s| !s.invariant(t.ssa.len() - 1)));
    }
}
