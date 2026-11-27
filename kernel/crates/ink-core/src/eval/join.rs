//! Exact equality at clock switches, including those between adjacent
//! floats. Fixed subexpressions have the language's float value. The real
//! expressions use rational arithmetic within a bounded work budget.
use super::{Ctx, fixed};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::I;
use crate::tape::{Op, Tape};

use num_rational::BigRational as Q;

const BUDGET: usize = 256;

fn zero() -> Q {
    Q::from_integer(0.into())
}

fn small(q: &Q) -> bool {
    q.numer().bits().max(q.denom().bits()) <= 4096
}

fn fixed_value(tape: &Tape, ctx: &Ctx, b: [I; 3], i: usize) -> Option<Q> {
    let iv = super::eval_interval(&tape.cone(i), ctx, b, None);
    (iv.lo == iv.hi).then(|| Q::from_float(iv.lo)).flatten()
}

/// m*t + b. At the join t is supplied as a constant, so both branch
/// expressions may be nonlinear there.
#[derive(Clone, PartialEq, Eq)]
struct Line(Q, Q);

struct At<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    b: [I; 3],
    time: Line,
    budget: usize,
}

impl At<'_> {
    fn read(&mut self, i: u32) -> Option<Line> {
        let line = self.expression(i)?;
        // Bound the cost of repeated squaring. Both operands of the next
        // operation have already passed this check.
        (small(&line.0) && small(&line.1)).then_some(line)
    }

    fn expression(&mut self, i: u32) -> Option<Line> {
        self.budget = self.budget.checked_sub(1)?;
        let ins = self.tape.ssa[i as usize];
        let constant = |q| Some(Line(zero(), q));
        match ins.op {
            Op::Const(x) => constant(Q::from_float(x)?),
            Op::Var(Var::T) => Some(self.time.clone()),
            Op::Var(Var::Param(p)) if self.ctx.p_spans.get(p).is_none() => constant(Q::from_float(self.ctx.params[p as usize])?),
            Op::Tag(_) | Op::Real => self.read(ins.a),
            _ if fixed(self.tape, self.ctx, i as usize) => constant(fixed_value(self.tape, self.ctx, self.b, i as usize)?),
            Op::Un(UnOp::Neg) => {
                let Line(m, b) = self.read(ins.a)?;
                Some(Line(-m, -b))
            }
            Op::Un(UnOp::Square) => {
                let Line(m, b) = self.read(ins.a)?;
                if m != zero() {
                    return None;
                }
                constant(&b * &b)
            }
            Op::Bin(op) => {
                let (Line(m, a), Line(n, b)) = (self.read(ins.a)?, self.read(ins.b)?);
                match op {
                    BinOp::Add => Some(Line(m + n, a + b)),
                    BinOp::Sub => Some(Line(m - n, a - b)),
                    BinOp::Mul if m == zero() || n == zero() => Some(Line(m * &b + n * &a, a * b)),
                    BinOp::Div if n == zero() && b != zero() => Some(Line(m / &b, a / b)),
                    _ => None,
                }
            }
            Op::Select => {
                let Line(m, b) = self.read(ins.a)?;
                if m != zero() {
                    return None;
                }
                if b == zero() {
                    // At a coincident switch a point value need not be
                    // either limit. Require both continuations to agree.
                    let left = self.read(ins.b)?;
                    return (left == self.read(ins.c)?).then_some(left);
                }
                self.read(if b < zero() { ins.b } else { ins.c })
            }
            _ => None,
        }
    }
}

pub(super) fn equal(tape: &Tape, ctx: &Ctx, b: [I; 3], i: usize) -> Option<()> {
    let select = tape.ssa[i];
    if ctx.t_span.is_none() || tape.fixed[select.a as usize] {
        return None;
    }
    let mut at = At { tape, ctx, b, time: Line(Q::from_integer(1.into()), zero()), budget: BUDGET };
    let Line(m, c) = at.read(select.a)?;
    if m == zero() {
        return None;
    }
    at.time = Line(zero(), -c / m);
    let (Line(_, left), Line(_, right)) = (at.read(select.b)?, at.read(select.c)?);
    (left == right).then_some(())
}

/// A rational times an integer power of pi. Degree conversion keeps pi
/// exact, and sine and cosine at quarter turns give rationals back.
#[derive(Clone, PartialEq, Eq)]
struct Exact {
    q: Q,
    pi: i16,
}

impl Exact {
    fn new(q: Q, pi: i16) -> Self {
        let pi = if q == zero() { 0 } else { pi };
        Self { q, pi }
    }

    fn integer(n: i32) -> Self {
        Self::new(Q::from_integer(n.into()), 0)
    }

    fn bin(&self, op: BinOp, b: &Self) -> Option<Self> {
        let (p, q) = (&self.q, &b.q);
        Some(match op {
            BinOp::Add | BinOp::Sub if q == &zero() => self.clone(),
            BinOp::Add if p == &zero() => b.clone(),
            BinOp::Sub if p == &zero() => Self::new(-q, b.pi),
            BinOp::Add if self.pi == b.pi => Self::new(p + q, self.pi),
            BinOp::Sub if self.pi == b.pi => Self::new(p - q, self.pi),
            BinOp::Mul => Self::new(p * q, self.pi.checked_add(b.pi)?),
            BinOp::Div if q != &zero() => Self::new(p / q, self.pi.checked_sub(b.pi)?),
            _ => return None,
        })
    }

    fn un(&self, op: UnOp) -> Option<Self> {
        let a = &self.q;
        Some(match op {
            UnOp::Neg => Self::new(-a, self.pi),
            UnOp::Abs => Self::new(if a < &zero() { -a } else { a.clone() }, self.pi),
            UnOp::Square => self.bin(BinOp::Mul, self)?,
            UnOp::Rad => Self::new(a / Self::integer(180).q, self.pi.checked_add(1)?),
            UnOp::Deg => Self::new(a * Self::integer(180).q, self.pi.checked_sub(1)?),
            UnOp::Sin | UnOp::Cos => {
                let turns = if a == &zero() {
                    0
                } else {
                    let quarters = a * Self::integer(2).q;
                    if self.pi != 1 || !quarters.is_integer() {
                        return None;
                    }
                    let four = Self::integer(4).q;
                    let rem = ((quarters % &four) + &four) % &four;
                    usize::try_from(rem.to_integer()).ok()?
                };
                let phase = (turns + usize::from(op == UnOp::Cos)) % 4;
                Self::integer([0, 1, 0, -1][phase])
            }
            _ => return None,
        })
    }
}

/// An expression continuous on one open interval. Its endpoint limits
/// can be unknown (e.g. a sine between quarter turns). Where affine and
/// both limits are exact, they also determine every value inside.
#[derive(Clone)]
struct Limits {
    ends: [Option<Exact>; 2],
    affine: bool,
}

impl Limits {
    fn constant(x: Exact) -> Self {
        Self { ends: [Some(x.clone()), Some(x)], affine: true }
    }

    fn constant_value(&self) -> Option<&Exact> {
        let a = self.ends[0].as_ref()?;
        (self.affine && self.ends[1].as_ref() == Some(a)).then_some(a)
    }

    fn affine_ends(&self) -> Result<[&Exact; 2], Unproved> {
        match &self.ends {
            [Some(a), Some(b)] if self.affine => Ok([a, b]),
            _ => Err(Unsupported),
        }
    }

    fn bin_ends(&self, op: BinOp, b: &Self) -> [Option<Exact>; 2] {
        std::array::from_fn(|k| {
            let (a, b) = (self.ends[k].as_ref(), b.ends[k].as_ref());
            // Even an unknown limit here is finite: the open piece was
            // built only out of finite, continuous operations.
            if op == BinOp::Mul && a.into_iter().chain(b).any(|x| x.q == zero()) {
                return Some(Exact::integer(0));
            }
            a?.bin(op, b?)
        })
    }
}

enum Unproved {
    Split(Q),
    Unsupported,
}

use Unproved::{Split, Unsupported};

/// One run over an open interval, or at a single moment if its ends are
/// equal. The open run uses the branch taken inside the interval even if
/// a different branch is taken at an endpoint.
struct Piece<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    b: [I; 3],
    time: &'a [Q; 2],
    seen: Vec<Option<Limits>>,
    budget: &'a mut usize,
    /// Only at a switch can the point value differ from the open limit.
    checks: [bool; 2],
}

impl<'a> Piece<'a> {
    fn new(tape: &'a Tape, ctx: &'a Ctx, b: [I; 3], time: &'a [Q; 2], budget: &'a mut usize) -> Self {
        Self { tape, ctx, b, time, seen: vec![None; tape.ssa.len()], budget, checks: [false; 2] }
    }

    fn read(&mut self, i: u32) -> Result<Limits, Unproved> {
        let i = i as usize;
        if let Some(v) = &self.seen[i] {
            return Ok(v.clone());
        }
        *self.budget = self.budget.checked_sub(1).ok_or(Unsupported)?;
        let mut v = self.expression(i)?;
        // An end too large to compare is dropped, not the piece: a switch on
        // this value, or a check at that end, then fails as unproved.
        for end in &mut v.ends {
            if end.as_ref().is_some_and(|x| !small(&x.q)) {
                *end = None;
                v.affine = false;
            }
        }
        self.seen[i] = Some(v.clone());
        Ok(v)
    }

    /// Does this affine value lie strictly below zero inside the piece?
    fn negative(&mut self, v: &Limits) -> Result<bool, Unproved> {
        let [left, right] = v.affine_ends()?;
        let (a, b) = (&left.q, &right.q);
        let z = zero();
        if (a < &z && b > &z) || (a > &z && b < &z) {
            if left.pi != right.pi {
                return Err(Unsupported);
            }
            return Err(Split(&self.time[0] + (&self.time[1] - &self.time[0]) * a / (a - b)));
        }
        self.checks[0] |= a == &z && b != &z;
        self.checks[1] |= b == &z && a != &z;
        Ok(a < &z || b < &z)
    }

    fn expression(&mut self, i: usize) -> Result<Limits, Unproved> {
        let ins = self.tape.ssa[i];
        match ins.op {
            Op::Var(Var::T) => Ok(Limits { ends: self.time.clone().map(|q| Some(Exact::new(q, 0))), affine: true }),
            Op::Tag(_) | Op::Real | Op::Deriv(3) => self.read(ins.a),
            _ if fixed(self.tape, self.ctx, i) => Ok(Limits::constant(Exact::new(fixed_value(self.tape, self.ctx, self.b, i).ok_or(Unsupported)?, 0))),
            Op::Un(op) => {
                let v = self.read(ins.a)?;
                if op == UnOp::Floor {
                    let [a, b] = v.affine_ends()?;
                    if a.pi != 0 || b.pi != 0 {
                        return Err(Unsupported);
                    }
                    let floor = (&a.q).min(&b.q).floor();
                    let next = Exact::new(&floor + Exact::integer(1).q, 0);
                    let ends = v.ends.map(|x| x?.bin(BinOp::Sub, &next));
                    self.negative(&Limits { ends, affine: true })?;
                    return Ok(Limits::constant(Exact::new(floor, 0)));
                }
                let affine = match op {
                    UnOp::Neg | UnOp::Rad | UnOp::Deg => v.affine,
                    UnOp::Abs | UnOp::Square | UnOp::Sin | UnOp::Cos => v.constant_value().is_some(),
                    _ => return Err(Unsupported),
                };
                Ok(Limits { ends: v.ends.map(|x| x?.un(op)), affine })
            }
            Op::Bin(op) => {
                let (a, b) = (self.read(ins.a)?, self.read(ins.b)?);
                let affine = a.affine && b.affine;
                if matches!(op, BinOp::Min | BinOp::Max) {
                    let ends = a.bin_ends(BinOp::Sub, &b);
                    let below = self.negative(&Limits { ends, affine })?;
                    return Ok(if below == (op == BinOp::Min) { a } else { b });
                }
                let affine = match op {
                    BinOp::Add | BinOp::Sub => affine,
                    BinOp::Mul => a.constant_value().is_some() && b.affine || b.constant_value().is_some() && a.affine,
                    BinOp::Div if b.constant_value().is_some_and(|x| x.q != zero()) => a.affine,
                    _ => return Err(Unsupported),
                };
                let ends = a.bin_ends(op, &b);
                Ok(Limits { ends, affine })
            }
            Op::Select => {
                let c = self.read(ins.a)?;
                let branch = if self.negative(&c)? { ins.b } else { ins.c };
                self.read(branch)
            }
            _ => Err(Unsupported),
        }
    }
}

/// Prove continuity of time-only expressions that switch on affine
/// pieces of the clock. Compare the exact left limit, point value and
/// right limit of the *whole expression* wherever its formula switches.
/// The same check covers drives a time derivative reads.
pub(super) fn continuous(tape: &Tape, ctx: &Ctx, b: [I; 3], span: (f64, f64)) -> Option<()> {
    let mut time = vec![[Q::from_float(span.0)?, Q::from_float(span.1)?]];
    // A switch takes separate runs on each side and at their endpoints.
    let mut budget = 4 * BUDGET;
    let roots = || std::iter::once(tape.ssa.len() as u32 - 1).chain(tape.ssa.iter().filter(|ins| matches!(ins.op, Op::Deriv(3))).map(|ins| ins.b));
    while let Some(ends) = time.pop() {
        let mut open = Piece::new(tape, ctx, b, &ends, &mut budget);
        match roots().map(|i| open.read(i)).collect::<Result<Vec<_>, _>>() {
            Ok(limits) => {
                let checks = open.checks;
                for (k, t) in ends.iter().enumerate() {
                    if !checks[k] {
                        continue;
                    }
                    let point = [t.clone(), t.clone()];
                    let mut at = Piece::new(tape, ctx, b, &point, &mut budget);
                    for (i, limit) in roots().zip(&limits) {
                        if at.read(i).ok()?.ends[0].as_ref()? != limit.ends[k].as_ref()? {
                            return None;
                        }
                    }
                }
            }
            Err(Split(t)) if small(&t) => time.extend([[ends[0].clone(), t.clone()], [t, ends[1].clone()]]),
            Err(_) => return None,
        }
    }
    Some(())
}
