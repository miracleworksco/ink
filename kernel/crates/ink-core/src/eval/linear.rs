//! Linear interval bounds over boxes and far sectors.

use super::context::Ctx;
use super::fork::{Fork, forked_live};
use super::gradient::Slopes;
use super::promise::Frames;
use super::range::interval_ins;
use crate::V3;
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I, add_o, mul_o, sub_o};
use crate::tape::{Op, Tape};

/// A value linear in x, y, z: coefficients and constant term, as intervals
/// that hold the real-number ones (a coefficient may be the range of
/// another value, as in a rotation by an angle that varies).
#[derive(Clone, Copy)]
struct Linear {
    a: [I; 3],
    d: I,
}

/// The region of `eval_interval_linear`: a box, or a far sector: the points
/// c + r (s q) (componentwise) for s in `s` and q in `q`, which lie on the
/// faces of the cube [-1, 1]^3. A sector reaches infinity when s does.
#[derive(Clone, Copy, Debug)]
pub enum Region {
    Box([I; 3]),
    Sector { c: V3, r: V3, q: [I; 3], s: I },
}

impl Region {
    /// The coordinates' ranges.
    pub fn bounds(&self) -> [I; 3] {
        match *self {
            Region::Box(b) => b,
            Region::Sector { c, r, q, s } => std::array::from_fn(|k| add_o(I::point(c[k]), mul_o(I::point(r[k]), mul_o(s, q[k])))),
        }
    }
    /// Over a sector with s > 0, 1 / s: the factor that takes a value
    /// bounded as s^k w to one of lower degree (`Far`).
    fn inverse(&self) -> Option<I> {
        match *self {
            Region::Sector { s, .. } if s.lo > 0.0 => Some(I::new((1.0 / s.hi).next_down().max(0.0), (1.0 / s.lo).next_up())),
            _ => None,
        }
    }
}

/// A value over a far sector bounded by its growth with the sector's
/// distance s: at every point at distance s it lies in s^k w (k >= 1).
/// Where the linear forms stop (a radius, sqrt(x^2 + y^2); an angle, atan2
/// of two coordinates) these go on: the radius of a sector is s times the
/// radius of its directions, its angle that of its directions, both narrow
/// however far it reaches, so a polar fold turns a narrow sector into a
/// narrow one. A value of lower degree is lifted by 1 / s (its range times
/// (1 / s)^j); one of degree 0 is its range.
#[derive(Clone, Copy)]
struct Far {
    k: u8,
    w: I,
}

/// The highest degree kept: a product of more coordinates is bounded by its
/// range alone.
const FAR_DEGREE: u8 = 4;

impl Far {
    /// A bound of degree `k` from one of degree `from` (`w`): w (1 / s)^(k - from).
    fn lift(w: I, from: u8, k: u8, inv: I) -> I {
        (from..k).fold(w, |w, _| mul_o(w, inv))
    }
    /// Its range over distances `s`: s^k w.
    fn range(&self, s: I) -> I {
        (0..self.k).fold(self.w, |w, _| mul_o(w, s))
    }
}

impl Linear {
    fn constant(d: I) -> Linear {
        Linear { a: [I::point(0.0); 3], d }
    }
    fn is_constant(&self) -> bool {
        self.a.iter().all(|c| c.lo == 0.0 && c.hi == 0.0)
    }
    fn map(&self, f: impl Fn(I) -> I) -> Linear {
        Linear { a: self.a.map(&f), d: f(self.d) }
    }
    fn zip(&self, o: &Linear, f: impl Fn(I, I) -> I) -> Linear {
        Linear { a: std::array::from_fn(|k| f(self.a[k], o.a[k])), d: f(self.d, o.d) }
    }
    /// Its range over the region. Over a sector, d + a.c + s (a.(r q)): the
    /// distance s multiplies the direction's sum once, so a rotated
    /// coordinate is one-signed on a narrow sector however far it reaches.
    fn range(&self, region: &Region) -> I {
        match *region {
            Region::Box(b) => (0..3).fold(self.d, |acc, k| add_o(acc, mul_o(self.a[k], b[k]))),
            Region::Sector { c, r, q, s } => {
                let near = (0..3).fold(self.d, |acc, k| add_o(acc, mul_o(self.a[k], I::point(c[k]))));
                let dir = (0..3).fold(I::point(0.0), |acc, k| add_o(acc, mul_o(mul_o(self.a[k], I::point(r[k])), q[k])));
                add_o(near, mul_o(s, dir))
            }
        }
    }
}

/// `eval_interval` over a region, where values linear in x, y, z (sums and
/// multiples of the coordinates: rotations of rotations, projections) are
/// also bounded from their collected coefficients, and each value takes the
/// tighter bound. Plain interval arithmetic loses the dependency between
/// such terms, which over an unbounded region loses everything (the far side
/// of a rotated part); collected, and over a far sector, a rotation is exact
/// up to the sector's width. Over a far sector each value is also bounded
/// by its growth with the distance (`Far`), so a radius and an angle stay
/// narrow too. Sound: every bound holds for the real-number value.
pub fn eval_interval_linear(tape: &Tape, ctx: &Ctx, region: &Region, choices: Option<&mut [u8]>) -> I {
    let far = region.inverse().map(|inv| (inv, vec![None; tape.nslots]));
    let mut run = Lin {
        tape,
        ctx,
        region,
        b: region.bounds(),
        s: vec![I::point(0.0); tape.nslots],
        lin: vec![None; tape.nslots],
        far,
        slopes: Slopes::of(tape, tape.nslots),
        frames: Frames::default(),
    };
    forked_live(tape, &mut run, choices);
    run.s[tape.out_slot()]
}

/// eval_interval_linear's state.
struct Lin<'a> {
    tape: &'a Tape,
    ctx: &'a Ctx,
    region: &'a Region,
    b: [I; 3],
    s: Vec<I>,
    lin: Vec<Option<Linear>>,
    /// Over a far sector, 1 / s and each slot's `Far` bound (None: degree
    /// 0, its range).
    far: Option<(I, Vec<Option<Far>>)>,
    slopes: Option<Slopes>,
    frames: Frames,
}

impl Lin<'_> {
    /// Slot `k`'s bound as degree and w (its range at degree 0).
    fn far_of(&self, k: u32) -> (u8, I) {
        match &self.far {
            Some((_, f)) => f[k as usize].map_or((0, self.s[k as usize]), |f| (f.k, f.w)),
            None => (0, self.s[k as usize]),
        }
    }
    /// Slots `a` and `b`'s bounds lifted to one degree, the higher.
    fn far_pair(&self, inv: I, a: u32, b: u32) -> (u8, I, I) {
        let ((ka, wa), (kb, wb)) = (self.far_of(a), self.far_of(b));
        let k = ka.max(kb);
        (k, Far::lift(wa, ka, k, inv), Far::lift(wb, kb, k, inv))
    }
    /// Instruction `ins`'s `Far` bound from its arguments' (`l` its linear
    /// form), or the bound of degree 0 it gives where its degree falls to
    /// 0 (an angle, a quotient). Read before its slot is written.
    fn far_ins(&self, ins: &crate::tape::Ins, l: Option<&Linear>) -> (Option<Far>, Option<I>) {
        let (Some((inv, _)), Region::Sector { c, r, q, .. }) = (&self.far, self.region) else { return (None, None) };
        let inv = *inv;
        // A linear form: d + a.c + s a.(r q), so w = (d + a.c) / s + a.(r q),
        // met with the rule's (d may be the range of a value of degree 1).
        let linear = l.filter(|l| !l.is_constant()).map(|l| {
            let near = (0..3).fold(l.d, |acc, k| add_o(acc, mul_o(l.a[k], I::point(c[k]))));
            let dir = (0..3).fold(I::point(0.0), |acc, k| add_o(acc, mul_o(mul_o(l.a[k], I::point(r[k])), q[k])));
            add_o(mul_o(near, inv), dir)
        });
        let (far, angle) = self.far_rule(ins, inv);
        match (far, linear) {
            (None, Some(w)) => (Some(Far { k: 1, w }), angle),
            (Some(f), Some(w)) if f.k == 1 => (Some(Far { k: 1, w: interval::meet(f.w, w) }), angle),
            _ => (far, angle),
        }
    }
    /// `far_ins` by the instruction's own rule, its linear form aside.
    fn far_rule(&self, ins: &crate::tape::Ins, inv: I) -> (Option<Far>, Option<I>) {
        let (k, w) = match ins.op {
            Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => self.far_of(ins.a),
            Op::Un(op @ (UnOp::Neg | UnOp::Abs)) => {
                let (k, w) = self.far_of(ins.a);
                (k, interval::un(op, w))
            }
            Op::Un(UnOp::Square) => {
                let (k, w) = self.far_of(ins.a);
                (2 * k, interval::un(UnOp::Square, w))
            }
            // Of an odd degree, lifted to the even one above.
            Op::Un(UnOp::Sqrt) => {
                let (k, w) = self.far_of(ins.a);
                let e = k + k % 2;
                (e / 2, interval::un(UnOp::Sqrt, Far::lift(w, k, e, inv)))
            }
            Op::Bin(BinOp::Mul) if ins.a == ins.b => {
                let (k, w) = self.far_of(ins.a);
                (2 * k, interval::un(UnOp::Square, w))
            }
            Op::Bin(BinOp::Mul) => {
                let ((ka, wa), (kb, wb)) = (self.far_of(ins.a), self.far_of(ins.b));
                (ka + kb, mul_o(wa, wb))
            }
            // A divisor one-signed: a quotient of lower degree is lifted.
            Op::Bin(BinOp::Div) => {
                let ((ka, wa), (kb, wb)) = (self.far_of(ins.a), self.far_of(ins.b));
                if !(wb.lo > 0.0 || wb.hi < 0.0) {
                    return (None, None);
                }
                let w = interval::div_o(wa, wb);
                if ka >= kb { (ka - kb, w) } else { (0, Far::lift(w, 0, kb - ka, inv)) }
            }
            // s^k > 0: a sum, a min or max, and an Euclidean max of two
            // bounds of one degree are those of their w; atan2's angle is
            // that of its w, of degree 0.
            Op::Bin(op @ (BinOp::Add | BinOp::Sub | BinOp::Min | BinOp::Max | BinOp::Emax | BinOp::Atan2)) => {
                let (k, wa, wb) = self.far_pair(inv, ins.a, ins.b);
                if k == 0 {
                    return (None, None);
                }
                (if op == BinOp::Atan2 { 0 } else { k }, interval::bin(op, wa, wb).0)
            }
            _ => return (None, None),
        };
        match k {
            0 => (None, Some(w)),
            k if k > FAR_DEGREE => (None, None),
            k => (Some(Far { k, w }), None),
        }
    }
}

impl Fork for Lin<'_> {
    type Out = (I, Option<Linear>, Option<interval::G>, Option<Far>);
    type Args = (I, I);
    /// An angle's coordinates over a far sector as their `Far` bounds of
    /// one degree: atan2 does not see a positive scale, so its jump's sides
    /// are those of the sector's directions.
    fn args(&self, i: usize) -> (I, I) {
        let ins = &self.tape.code[i];
        if let (Op::Bin(BinOp::Atan2), Some((inv, _))) = (ins.op, &self.far) {
            let (k, y, x) = self.far_pair(*inv, ins.a, ins.b);
            if k > 0 {
                return (y, x);
            }
        }
        (self.s[ins.a as usize], self.s[ins.b as usize])
    }
    fn bounds(args: &(I, I)) -> (I, I) {
        *args
    }
    fn ctx(&self) -> &Ctx {
        self.ctx
    }
    fn code(&self) -> &[crate::tape::Ins] {
        &self.tape.code
    }
    fn value(&self, k: u32) -> I {
        self.s[k as usize]
    }
    #[inline(always)]
    fn step(&mut self, i: usize, choices: Option<&mut [u8]>) -> I {
        let (tape, ctx, region, b) = (self.tape, self.ctx, self.region, self.b);
        let (s, lin) = (&self.s, &self.lin);
        let ins = &tape.code[i];
        let mut v = interval_ins(tape, ctx, i, &tape.code, s, b, choices, Some(&mut self.frames), true);
        if let Some(sl) = &mut self.slopes {
            v = sl.tighten(tape, ctx, i, &tape.code, b, v, true);
        }
        let arg = |k: u32| lin[k as usize];
        // A value's linear form, or its range as a constant term (to scale
        // another's coefficients by).
        let scale = |k: u32| lin[k as usize].filter(|l| l.is_constant()).map_or(s[k as usize], |l| l.d);
        // Time is a constant term too: it does not vary over the region.
        let l = if tape.fixed[i] || matches!(ins.op, Op::Const(_) | Op::Within(_) | Op::Var(Var::Param(_) | Var::T)) {
            Some(Linear::constant(v))
        } else {
            match ins.op {
                Op::Var(Var::X) | Op::Var(Var::Y) | Op::Var(Var::Z) => {
                    let k = match ins.op {
                        Op::Var(Var::X) => 0,
                        Op::Var(Var::Y) => 1,
                        _ => 2,
                    };
                    let mut l = Linear::constant(I::point(0.0));
                    l.a[k] = I::point(1.0);
                    Some(l)
                }
                Op::Un(UnOp::Neg) => arg(ins.a).map(|l| l.map(|c| I::new(-c.hi, -c.lo))),
                Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => arg(ins.a),
                Op::Bin(BinOp::Add) => arg(ins.a).zip(arg(ins.b)).map(|(p, q)| p.zip(&q, add_o)),
                Op::Bin(BinOp::Sub) => arg(ins.a).zip(arg(ins.b)).map(|(p, q)| p.zip(&q, sub_o)),
                // A linear value times another value: the other's range
                // scales its coefficients (exact for a constant).
                Op::Bin(BinOp::Mul) => match (arg(ins.a), arg(ins.b)) {
                    (Some(p), q) if q.is_none_or(|q| q.is_constant()) => Some(p.map(|c| mul_o(c, scale(ins.b)))),
                    (p, Some(q)) if p.is_none_or(|p| p.is_constant()) => Some(q.map(|c| mul_o(c, scale(ins.a)))),
                    _ => None,
                },
                Op::Bin(BinOp::Div) => match arg(ins.a) {
                    Some(p) if arg(ins.b).is_none_or(|q| q.is_constant()) && (scale(ins.b).lo > 0.0 || scale(ins.b).hi < 0.0) => {
                        Some(p.map(|c| interval::bin(BinOp::Div, c, scale(ins.b)).0))
                    }
                    _ => None,
                },
                // One-signed over the region: the form itself, or negated.
                Op::Un(UnOp::Abs) => arg(ins.a).and_then(|l| match s[ins.a as usize] {
                    x if x.lo >= 0.0 => Some(l),
                    x if x.hi <= 0.0 => Some(l.map(|c| I::new(-c.hi, -c.lo))),
                    _ => None,
                }),
                // A member that wins everywhere over the region: its form.
                Op::Bin(BinOp::Min | BinOp::Max) => {
                    let (p, q) = (s[ins.a as usize], s[ins.b as usize]);
                    let first = match ins.op {
                        Op::Bin(BinOp::Min) => (p.hi <= q.lo, q.hi <= p.lo),
                        _ => (p.lo >= q.hi, q.lo >= p.hi),
                    };
                    match first {
                        (true, _) => arg(ins.a),
                        (_, true) => arg(ins.b),
                        _ => None,
                    }
                }
                _ => None,
            }
        };
        // A value that does not depend on x, y, z (of time, say) is a
        // constant term.
        let l = l.or_else(|| {
            let (args, k) = ins.args();
            (k > 0 && args[..k].iter().all(|&a| lin[a as usize].is_some_and(|l| l.is_constant()))).then(|| Linear::constant(v))
        });
        if let Some(l) = l {
            v = interval::meet(v, l.range(region));
        }
        let (far, angle) = self.far_ins(ins, l.as_ref());
        if let (Some(f), Region::Sector { s, .. }) = (far, region) {
            v = interval::meet(v, f.range(*s));
        }
        if let Some(a) = angle {
            v = interval::meet(v, a);
        }
        self.s[ins.out as usize] = v;
        self.lin[ins.out as usize] = l;
        if let Some((_, f)) = &mut self.far {
            f[ins.out as usize] = far;
        }
        v
    }
    fn side(&mut self, i: usize, _: &(I, I), _: usize, side: I) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = side;
        // A floor's side is one number over the region.
        self.lin[o] = (side.lo == side.hi).then(|| Linear::constant(side));
        if let Some((_, f)) = &mut self.far {
            f[o] = None;
        }
        if let Some(sl) = &mut self.slopes {
            sl.0[o].v = interval::meet(sl.0[o].v, side);
        }
    }
    fn out(&self, i: usize) -> Self::Out {
        let o = self.tape.code[i].out as usize;
        (self.s[o], self.lin[o], self.slopes.as_ref().map(|sl| sl.0[o]), self.far.as_ref().and_then(|(_, f)| f[o]))
    }
    fn put(&mut self, i: usize, (v, l, g, far): Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = v;
        self.lin[o] = l;
        if let Some((_, f)) = &mut self.far {
            f[o] = far;
        }
        if let (Some(sl), Some(g)) = (&mut self.slopes, g) {
            sl.0[o] = g;
        }
    }
    fn join(&mut self, i: usize, a: Self::Out, b: Self::Out) {
        let o = self.tape.code[i].out as usize;
        self.s[o] = interval::hull(a.0, b.0);
        // Coefficients that hold either form's (a value is one or the other).
        self.lin[o] = a.1.zip(b.1).map(|(p, q)| p.zip(&q, interval::hull));
        if let Some((inv, f)) = &mut self.far {
            // Each side's bound lifted to the higher degree.
            let (ka, wa) = a.3.map_or((0, a.0), |f| (f.k, f.w));
            let (kb, wb) = b.3.map_or((0, b.0), |f| (f.k, f.w));
            let k = ka.max(kb);
            f[o] = (k > 0).then(|| Far { k, w: interval::hull(Far::lift(wa, ka, k, *inv), Far::lift(wb, kb, k, *inv)) });
        }
        if let (Some(sl), Some(ga), Some(gb)) = (&mut self.slopes, a.2, b.2) {
            sl.0[o] = ga.hull(&gb);
        }
    }
}
