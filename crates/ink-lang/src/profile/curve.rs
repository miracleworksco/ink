//! Points of a drafted profile that move on curves as its walls lean: a
//! line wall meeting an arc's (a plane and a cone) or two arcs' meeting, and
//! where a hull line touches an arc (from a point, or common to two arcs).
//! Each is a closed form in the span's parameter t (0 at z0, 1 at z1), one
//! expression evaluated in floats, in outward-rounded bounds over parts of
//! the span (`Sl`: the load-time proof) and lowered into the field, so that
//! the proof and the field read the same curve.

use super::*;
use ink_core::interval::{self, I};
use std::rc::Rc;

/// An expression in t: a constant, a + b t, or an operation on two, or a
/// square root. A curve's own expression folds constants and linear terms
/// in floats as it is built (so that intervals of a difference of two
/// nearly equal linear terms stay tight): that is its definition, which the
/// field and the proof both read. Its derivative (`Ex::d`) must be the
/// derivative of exactly that, so it folds only where the result is exact.
#[derive(Debug)]
pub(crate) enum Ex {
    K(f64),
    Lin(f64, f64),
    Bin(BinOp, X, X),
    Sqrt(X),
}

pub(crate) type X = Rc<Ex>;
/// A moving point.
pub(crate) type XP = [X; 2];

pub(crate) fn kx(v: f64) -> X {
    Rc::new(Ex::K(v))
}

/// a + s t.
fn lin(a: f64, s: f64) -> X {
    if s == 0.0 { kx(a) } else { Rc::new(Ex::Lin(a, s)) }
}

fn bx(op: BinOp, a: &X, b: &X) -> X {
    make(op, a, b, true)
}

/// a op b; `round`: folding may round (a curve's definition), else only
/// exact folds (a derivative).
fn make(op: BinOp, a: &X, b: &X, round: bool) -> X {
    // x op y folded: in floats, or where it is exact (its outward interval
    // a point).
    let f = |op: BinOp, x: f64, y: f64| {
        if round {
            return Some(f64::bin(op, x, y));
        }
        let r = interval::bin(op, I::point(x), I::point(y)).0;
        (r.lo == r.hi).then_some(r.lo)
    };
    let k = |x: &Ex| if let Ex::K(v) = *x { Some(v) } else { None };
    // Exact identities: x + 0, x - 0, x 1, x 0.
    match (op, k(a), k(b)) {
        (BinOp::Add, Some(0.0), _) => return b.clone(),
        (BinOp::Add | BinOp::Sub, _, Some(0.0)) => return a.clone(),
        (BinOp::Mul, Some(0.0), _) | (BinOp::Mul, _, Some(0.0)) => return kx(0.0),
        (BinOp::Mul, Some(1.0), _) => return b.clone(),
        (BinOp::Mul | BinOp::Div, _, Some(1.0)) => return a.clone(),
        _ => {}
    }
    if let (Some(x), Some(y)) = (k(a), k(b))
        && let Some(v) = f(op, x, y)
    {
        return kx(v);
    }
    // A linear term's coefficients: a + s t.
    let coef = |x: &Ex| match *x {
        Ex::K(v) => Some((v, 0.0)),
        Ex::Lin(a, s) => Some((a, s)),
        _ => None,
    };
    if let (Some((a0, a1)), Some((b0, b1))) = (coef(a), coef(b)) {
        let folded = match op {
            BinOp::Add => f(BinOp::Add, a0, b0).zip(f(BinOp::Add, a1, b1)),
            BinOp::Sub => f(BinOp::Sub, a0, b0).zip(f(BinOp::Sub, a1, b1)),
            BinOp::Mul if a1 == 0.0 => f(BinOp::Mul, a0, b0).zip(f(BinOp::Mul, a0, b1)),
            BinOp::Mul if b1 == 0.0 => f(BinOp::Mul, a0, b0).zip(f(BinOp::Mul, a1, b0)),
            _ => None,
        };
        if let Some((c0, c1)) = folded {
            return lin(c0, c1);
        }
    }
    Rc::new(Ex::Bin(op, a.clone(), b.clone()))
}
pub(crate) fn xadd(a: &X, b: &X) -> X {
    bx(BinOp::Add, a, b)
}
pub(crate) fn xsub(a: &X, b: &X) -> X {
    bx(BinOp::Sub, a, b)
}
pub(crate) fn xmul(a: &X, b: &X) -> X {
    bx(BinOp::Mul, a, b)
}
pub(crate) fn xdiv(a: &X, b: &X) -> X {
    bx(BinOp::Div, a, b)
}
/// The operations a condition on the paths is built from (walls.rs): only
/// exact folds (`make`), so that its bounds hold of the real value of the
/// paths it reads, which the field lowers as they are.
pub(crate) mod exact {
    use super::*;
    pub(crate) use super::{kx, pk};
    pub(crate) fn xadd(a: &X, b: &X) -> X {
        make(BinOp::Add, a, b, false)
    }
    pub(crate) fn xsub(a: &X, b: &X) -> X {
        make(BinOp::Sub, a, b, false)
    }
    pub(crate) fn xmul(a: &X, b: &X) -> X {
        make(BinOp::Mul, a, b, false)
    }
    pub(crate) fn xdiv(a: &X, b: &X) -> X {
        make(BinOp::Div, a, b, false)
    }
    pub(crate) fn psub(a: &XP, b: &XP) -> XP {
        [xsub(&a[0], &b[0]), xsub(&a[1], &b[1])]
    }
    pub(crate) fn padd(a: &XP, b: &XP) -> XP {
        [xadd(&a[0], &b[0]), xadd(&a[1], &b[1])]
    }
    pub(crate) fn pscale(a: &XP, s: &X) -> XP {
        [xmul(&a[0], s), xmul(&a[1], s)]
    }
    pub(crate) fn pdot(a: &XP, b: &XP) -> X {
        xadd(&xmul(&a[0], &b[0]), &xmul(&a[1], &b[1]))
    }
    pub(crate) fn pcross(a: &XP, b: &XP) -> X {
        xsub(&xmul(&a[0], &b[1]), &xmul(&a[1], &b[0]))
    }
    pub(crate) fn pnorm(a: &XP) -> X {
        Rc::new(Ex::Sqrt(pdot(a, a)))
    }
}

pub(crate) fn xsqrt(a: &X) -> X {
    if let Ex::K(v) = **a {
        return kx(f64::sqrt(v));
    }
    Rc::new(Ex::Sqrt(a.clone()))
}
/// From a at t = 0 to b at t = 1: fixed only where they are equal. (Motion
/// small beside a coordinate or a radius is real: a wall 1e12 from the
/// origin moving 0.5, a radius of 1e8 growing 5e-5. Read as rounding, the
/// paths froze while the walls moved.)
pub(crate) fn xlin(a: f64, b: f64) -> X {
    lin(a, b - a)
}
pub(crate) fn plin(a: P, b: P) -> XP {
    [xlin(a[0], b[0]), xlin(a[1], b[1])]
}
pub(crate) fn pk(a: P) -> XP {
    [kx(a[0]), kx(a[1])]
}
pub(crate) fn psub(a: &XP, b: &XP) -> XP {
    [xsub(&a[0], &b[0]), xsub(&a[1], &b[1])]
}
pub(crate) fn padd(a: &XP, b: &XP) -> XP {
    [xadd(&a[0], &b[0]), xadd(&a[1], &b[1])]
}
pub(crate) fn pscale(a: &XP, s: &X) -> XP {
    [xmul(&a[0], s), xmul(&a[1], s)]
}
pub(crate) fn pdot(a: &XP, b: &XP) -> X {
    xadd(&xmul(&a[0], &b[0]), &xmul(&a[1], &b[1]))
}
/// a x b.
pub(crate) fn pcross(a: &XP, b: &XP) -> X {
    xsub(&xmul(&a[0], &b[1]), &xmul(&a[1], &b[0]))
}
pub(crate) fn pnorm(a: &XP) -> X {
    xsqrt(&pdot(a, a))
}

/// A value an expression evaluates to: a float, or an interval holding the
/// real result (rounded outward).
pub(crate) trait Val: Copy {
    fn k(v: f64) -> Self;
    fn bin(op: BinOp, a: Self, b: Self) -> Self;
    fn sqrt(a: Self) -> Self;
}

impl Val for f64 {
    fn k(v: f64) -> f64 {
        v
    }
    fn bin(op: BinOp, a: f64, b: f64) -> f64 {
        match op {
            BinOp::Add => a + b,
            BinOp::Sub => a - b,
            BinOp::Mul => a * b,
            BinOp::Div => a / b,
            _ => unreachable!("a curve's expression adds, subtracts, multiplies and divides"),
        }
    }
    fn sqrt(a: f64) -> f64 {
        a.max(0.0).sqrt()
    }
}

impl Val for I {
    fn k(v: f64) -> I {
        I::point(v)
    }
    fn bin(op: BinOp, a: I, b: I) -> I {
        interval::bin(op, a, b).0
    }
    fn sqrt(a: I) -> I {
        interval::un(UnOp::Sqrt, a)
    }
}

/// A node's key in a memo: an expression shares its parts (a contact's
/// formula reads its circle's centre several times, a partner's curve
/// inside it), so each is evaluated once.
fn key(x: &X) -> usize {
    Rc::as_ptr(x) as usize
}

impl Ex {
    pub(crate) fn at<V: Val>(self: &Rc<Self>, t: V) -> V {
        self.at_in(t, &mut HashMap::new())
    }

    fn at_in<V: Val>(self: &Rc<Self>, t: V, memo: &mut HashMap<usize, V>) -> V {
        if let Some(&v) = memo.get(&key(self)) {
            return v;
        }
        let v = match &**self {
            Ex::K(v) => V::k(*v),
            Ex::Lin(a, b) => V::bin(BinOp::Add, V::k(*a), V::bin(BinOp::Mul, V::k(*b), t)),
            Ex::Bin(op, a, b) => {
                let (a, b) = (a.at_in(t, memo), b.at_in(t, memo));
                V::bin(*op, a, b)
            }
            Ex::Sqrt(a) => V::sqrt(a.at_in(t, memo)),
        };
        memo.insert(key(self), v);
        v
    }

    /// Its derivative in t.
    pub(crate) fn d(x: &X) -> X {
        Ex::d_in(x, &mut HashMap::new())
    }

    fn d_in(x: &X, memo: &mut HashMap<usize, X>) -> X {
        if let Some(d) = memo.get(&key(x)) {
            return d.clone();
        }
        let m = |op: BinOp, a: &X, b: &X| make(op, a, b, false);
        let d = match &**x {
            Ex::K(_) => kx(0.0),
            Ex::Lin(_, s) => kx(*s),
            Ex::Bin(op, a, b) => {
                let (da, db) = (Ex::d_in(a, memo), Ex::d_in(b, memo));
                match op {
                    BinOp::Add => m(BinOp::Add, &da, &db),
                    BinOp::Sub => m(BinOp::Sub, &da, &db),
                    BinOp::Mul => m(BinOp::Add, &m(BinOp::Mul, &da, b), &m(BinOp::Mul, a, &db)),
                    BinOp::Div => m(BinOp::Div, &m(BinOp::Sub, &m(BinOp::Mul, &da, b), &m(BinOp::Mul, a, &db)), &m(BinOp::Mul, b, b)),
                    _ => unreachable!(),
                }
            }
            Ex::Sqrt(a) => m(BinOp::Div, &Ex::d_in(a, memo), &m(BinOp::Mul, &kx(2.0), x)),
        };
        memo.insert(key(x), d.clone());
        d
    }

    /// The quantities under its square roots (each must stay positive).
    pub(crate) fn radicands(x: &X, out: &mut Vec<X>) {
        match &**x {
            Ex::K(_) | Ex::Lin(..) => {}
            Ex::Bin(_, a, b) => {
                Ex::radicands(a, out);
                Ex::radicands(b, out);
            }
            Ex::Sqrt(a) => {
                out.push(a.clone());
                Ex::radicands(a, out);
            }
        }
    }
}

/// Whether an expression is linear in t.
pub(crate) fn linear(x: &X) -> bool {
    matches!(**x, Ex::K(_) | Ex::Lin(..))
}

pub(crate) fn at<V: Val>(p: &XP, t: V) -> [V; 2] {
    let mut memo = HashMap::new();
    [p[0].at_in(t, &mut memo), p[1].at_in(t, &mut memo)]
}

/// Where a line wall (unit normal n, n . p = k(t)) meets a circle wall
/// (centre c, radius r(t)), on side s along n's quarter turn from their
/// foot.
pub(crate) fn line_circle(n: P, k: [f64; 2], c: P, r: [f64; 2], s: f64) -> XP {
    let h = xlin(k[0] - dot(n, c), k[1] - dot(n, c));
    let r = xlin(r[0], r[1]);
    let q = xsqrt(&xmul(&xsub(&r, &h), &xadd(&r, &h)));
    let foot = padd(&pk(c), &pscale(&pk(n), &h));
    padd(&foot, &pscale(&pk([-n[1] * s, n[0] * s]), &q))
}

/// Where two circle walls (fixed centres, radii r1(t), r2(t)) meet, on side
/// s along the quarter turn from c1 to c2. The square of the half-chord,
/// r1^2 - along^2, is written as its four linear factors over (2 d)^2, so
/// that circles nearly tangent read their small factor (d - |r1 - r2|, or
/// r1 + r2 - d) directly, not as a difference of squares of their size.
pub(crate) fn circle_circle(c1: P, r1: [f64; 2], c2: P, r2: [f64; 2], s: f64) -> XP {
    let e = sub(c2, c1);
    let d = norm(e);
    let u = mul(e, 1.0 / d);
    let (r1, r2) = (xlin(r1[0], r1[1]), xlin(r2[0], r2[1]));
    let (sum, dif, kd) = (xadd(&r1, &r2), xsub(&r1, &r2), kx(d));
    let along = xdiv(&xadd(&kx(d * d), &xmul(&dif, &sum)), &kx(2.0 * d));
    let q4 = xmul(&xmul(&xsub(&sum, &kd), &xadd(&sum, &kd)), &xmul(&xsub(&kd, &dif), &xadd(&kd, &dif)));
    let q = xsqrt(&xdiv(&q4, &kx(4.0 * d * d)));
    padd(&padd(&pk(c1), &pscale(&pk(u), &along)), &pscale(&pk([-u[1] * s, u[0] * s]), &q))
}

/// Where the common tangent of circles (c, r) and (c2, r2) touches the
/// first (a point is a circle of radius 0), on side s: its normal n has
/// n . (c2 - c) = r - r2, so n = ((r - r2) d + s w d') / |d|^2, d = c2 - c, d'
/// its quarter turn, w = sqrt(|d|^2 - (r - r2)^2), read as (|d| - (r - r2))
/// (|d| + (r - r2)) so that circles nearly tangent read their small factor
/// directly; the point is c + r n.
///
/// Built of exact operations only (`exact`), so that what follows holds of
/// the expressions the field lowers, not just of the formula: |d'|^2 = |d|^2
/// and d . d' = 0, so n . n = 1 and n . d = r - r2. The point is on its
/// circle (|p - c| = |r|), and the line from it to c2 + r2 n (a point: c2;
/// or the other circle's contact, built so on side -s, whose n is this n)
/// is at right angles to n: tangent to both, at distance |c2 + r2 n - p| =
/// w apart. Where w > 0 (its radicands' condition) n and the line move
/// continuously, so which side of the line the circle lies on is fixed.
///
/// The graph folds an operation on two constants in floats as it is
/// lowered (`Graph::bin`, `Graph::un`): so is the point here, as the field
/// will read it (`lowered`), and whether any such fold rounded (then the
/// identities above hold of the formula but not of what the field reads).
pub(crate) fn touch(c: &XP, r: &X, c2: &XP, r2: &X, s: f64) -> (XP, bool) {
    use exact::{padd, pdot, pscale, psub, xadd, xdiv, xmul, xsub};
    let xsqrt = |a: &X| Rc::new(Ex::Sqrt(a.clone()));
    let pperp = |a: &XP| [xsub(&kx(0.0), &a[1]), a[0].clone()];
    let d = psub(c2, c);
    let dd = pdot(&d, &d);
    let del = xsub(r, r2);
    let len = xsqrt(&dd);
    let w = xsqrt(&xmul(&xsub(&len, &del), &xadd(&len, &del)));
    let n = pscale(&padd(&pscale(&d, &del), &pscale(&pperp(&d), &xmul(&kx(s), &w))), &xdiv(&kx(1.0), &dd));
    let p = padd(c, &pscale(&n, r));
    let (mut memo, mut rounded) = (HashMap::new(), false);
    let p = [lowered(&p[0], &mut memo, &mut rounded), lowered(&p[1], &mut memo, &mut rounded)];
    (p, rounded)
}

/// An expression with its operations on constants folded in floats, as the
/// graph folds them when it is lowered (`Graph::bin`, `Graph::un`: but an
/// undefined case, which it keeps); `rounded` set where one is folded (an
/// exact expression keeps such an operation only where its fold rounds).
fn lowered(x: &X, memo: &mut HashMap<usize, X>, rounded: &mut bool) -> X {
    if let Some(y) = memo.get(&key(x)) {
        return y.clone();
    }
    let k = |x: &X| if let Ex::K(v) = **x { Some(v) } else { None };
    let y = match &**x {
        Ex::K(_) | Ex::Lin(..) => x.clone(),
        Ex::Bin(op, a, b) => {
            let (a, b) = (lowered(a, memo, rounded), lowered(b, memo, rounded));
            match (k(&a), k(&b)) {
                (Some(u), Some(v)) if ink_core::graph::undefined_bin(*op, u, v).is_none() => {
                    *rounded = true;
                    kx(ink_core::graph::bin_f64(*op, u, v))
                }
                _ => Rc::new(Ex::Bin(*op, a, b)),
            }
        }
        Ex::Sqrt(a) => {
            let a = lowered(a, memo, rounded);
            match k(&a) {
                Some(v) => {
                    *rounded = true;
                    kx(ink_core::graph::un_f64(UnOp::Sqrt, v))
                }
                _ => Rc::new(Ex::Sqrt(a)),
            }
        }
    };
    memo.insert(key(x), y.clone());
    y
}

/// The side of two that puts a curve at t nearest p, and which.
pub(crate) fn nearer(f: impl Fn(f64) -> XP, p: P, t: f64) -> (XP, f64) {
    let (a, b) = (f(1.0), f(-1.0));
    let off = |x: &XP| norm(sub(at(x, t), p));
    if off(&a) <= off(&b) { (a, 1.0) } else { (b, -1.0) }
}

/// An expression's bounds over a part of the span in second-order slope
/// form: its value v and derivative d at the part's middle m, and a second
/// slope s with f(t) = f(m) + f'(m) (t - m) + s (t - m)^2 for every t in
/// the part (t - m is dt), and its bounds there, met from its operands'
/// bounds and from that form. Their width falls with the part's square
/// times the expression's own second derivative, where plain bounds' falls
/// with the part times the spread its terms' dependence adds: a quantity
/// read as a difference of nearly equal terms (a point barely off a circle,
/// circles nearly tangent), or through a square root whose radicand is
/// small, settles in a few parts where its value moves smoothly.
#[derive(Clone, Copy)]
pub(crate) struct Sl {
    v: I,
    d: I,
    s: I,
    pub(crate) e: I,
    dt: I,
}

fn ib(op: BinOp, a: I, b: I) -> I {
    interval::bin(op, a, b).0
}

impl Sl {
    /// t itself over the part.
    pub(crate) fn of(t: I) -> Sl {
        let m = I::point(t.mid());
        Sl { v: m, d: I::point(1.0), s: I::point(0.0), e: t, dt: ib(BinOp::Sub, t, m) }
    }

    /// Its first slope over the part: d + s dt.
    fn slope(&self) -> I {
        ib(BinOp::Add, self.d, ib(BinOp::Mul, self.s, self.dt))
    }

    fn with(v: I, d: I, s: I, plain: I, dt: I) -> Sl {
        let sq = I::new(0.0, interval::mag(ib(BinOp::Mul, dt, dt)));
        let form = ib(BinOp::Add, ib(BinOp::Add, v, ib(BinOp::Mul, d, dt)), ib(BinOp::Mul, s, sq));
        Sl { v, d, s, e: interval::meet(plain, form), dt }
    }
}

impl Val for Sl {
    fn k(v: f64) -> Sl {
        let (p, z) = (I::point(v), I::point(0.0));
        Sl { v: p, d: z, s: z, e: p, dt: z }
    }
    fn bin(op: BinOp, a: Sl, b: Sl) -> Sl {
        let dt = if a.dt.lo != 0.0 || a.dt.hi != 0.0 { a.dt } else { b.dt };
        let (mul, add, sub, div) = (BinOp::Mul, BinOp::Add, BinOp::Sub, BinOp::Div);
        let v = ib(op, a.v, b.v);
        let (d, s) = match op {
            BinOp::Add | BinOp::Sub => (ib(op, a.d, b.d), ib(op, a.s, b.s)),
            // gh - g(m)h(m) - (g'h + gh')(m) dt = (g(m) sh + h(m) sg + G H) dt^2,
            // G and H the first slopes.
            BinOp::Mul => (ib(add, ib(mul, a.d, b.v), ib(mul, a.v, b.d)), ib(add, ib(add, ib(mul, a.v, b.s), ib(mul, b.v, a.s)), ib(mul, a.slope(), b.slope()))),
            // g/h: f' = (g' - f h') / h; the rest over dt^2 is
            // ((sg - f(m) sh) h(m) - (g' - f h')(m) H) / (h(t) h(m)).
            BinOp::Div => {
                let am = ib(sub, a.d, ib(mul, v, b.d));
                let num = ib(sub, ib(mul, ib(sub, a.s, ib(mul, v, b.s)), b.v), ib(mul, am, b.slope()));
                (ib(div, am, b.v), ib(div, num, ib(mul, b.e, b.v)))
            }
            _ => unreachable!("a curve's expression adds, subtracts, multiplies and divides"),
        };
        Sl::with(v, d, s, ib(op, a.e, b.e), dt)
    }
    fn sqrt(a: Sl) -> Sl {
        let (mul, sub, div) = (BinOp::Mul, BinOp::Sub, BinOp::Div);
        let v = <I as Val>::sqrt(a.v);
        let plain = <I as Val>::sqrt(a.e);
        // Where the radicand's bounds stay above 0: f' = g' / 2f, and the
        // rest over dt^2 is (2 f(m) sg - g'(m) F) / (2 f(m) (f(t) + f(m))),
        // F = G / (f(t) + f(m)) the first slope.
        if !(v.lo > 0.0 && a.e.lo > 0.0) {
            return Sl::with(v, I::ENTIRE, I::ENTIRE, plain, a.dt);
        }
        let two_f = ib(mul, I::point(2.0), v);
        let sum = ib(BinOp::Add, plain, v);
        let first = ib(div, a.slope(), sum);
        let s = ib(div, ib(sub, ib(mul, two_f, a.s), ib(mul, a.d, first)), ib(mul, two_f, sum));
        Sl::with(v, ib(div, a.d, two_f), s, plain, a.dt)
    }
}

/// Expressions flattened once (each shared part once, in order), to be
/// evaluated often: a condition over the parts of a span.
pub(crate) struct Tape {
    ops: Vec<Op>,
    roots: Vec<usize>,
}

enum Op {
    K(f64),
    Lin(f64, f64),
    Bin(BinOp, usize, usize),
    Sqrt(usize),
}

impl Tape {
    pub(crate) fn new(xs: &[X]) -> Tape {
        fn put(x: &X, ops: &mut Vec<Op>, at: &mut HashMap<usize, usize>) -> usize {
            if let Some(&i) = at.get(&key(x)) {
                return i;
            }
            let op = match &**x {
                Ex::K(v) => Op::K(*v),
                Ex::Lin(a, b) => Op::Lin(*a, *b),
                Ex::Bin(op, a, b) => {
                    let (a, b) = (put(a, ops, at), put(b, ops, at));
                    Op::Bin(*op, a, b)
                }
                Ex::Sqrt(a) => Op::Sqrt(put(a, ops, at)),
            };
            ops.push(op);
            at.insert(key(x), ops.len() - 1);
            ops.len() - 1
        }
        let (mut ops, mut at) = (Vec::new(), HashMap::new());
        let roots = xs.iter().map(|x| put(x, &mut ops, &mut at)).collect();
        Tape { ops, roots }
    }

    /// Its expressions' values at t.
    pub(crate) fn at<V: Val>(&self, t: V) -> Vec<V> {
        let mut v: Vec<V> = Vec::with_capacity(self.ops.len());
        for op in &self.ops {
            let x = match *op {
                Op::K(c) => V::k(c),
                Op::Lin(a, b) => V::bin(BinOp::Add, V::k(a), V::bin(BinOp::Mul, V::k(b), t)),
                Op::Bin(op, a, b) => V::bin(op, v[a], v[b]),
                Op::Sqrt(a) => V::sqrt(v[a]),
            };
            v.push(x);
        }
        self.roots.iter().map(|&r| v[r]).collect()
    }
}

/// How far from t0 toward t1 (either way) `holds` is proven, given whether
/// it holds over all of a part of the span (an interval of t): parts it
/// does not settle are halved, nearest t0 first, and where it fails at a
/// part's middle, the part toward t0 is halved on to find where. The
/// furthest t it is proven to, and where it was found to fail (both t1
/// where it holds throughout); a part too narrow to halve that it does not
/// settle fails there.
pub(crate) fn reach(holds: impl Fn(I) -> bool, t0: f64, t1: f64) -> (f64, f64) {
    let part = |a: f64, b: f64| I::new(a.min(b), a.max(b));
    let (mut upto, mut fail) = (t0, t1);
    let mut todo = vec![(t0, t1)];
    while let Some((a, b)) = todo.pop() {
        if holds(part(a, b)) {
            upto = b;
            continue;
        }
        let m = 0.5 * (a + b);
        if m == a || m == b {
            return (upto, b);
        }
        if holds(I::point(m)) {
            todo.push((m, b));
            todo.push((a, m));
        } else {
            todo.clear();
            todo.push((a, m));
            fail = m;
        }
    }
    (upto, fail)
}

/// The bounds of f over t in [a, b], given its bounds over a part of it
/// (in slope form, tight to the part's square): over the whole, halved where
/// a bound is unbounded (a derivative divides by a square root whose
/// radicand's bound reaches 0 there, though the radicand does not). None
/// where it is unbounded at a point.
pub(crate) fn hull(f: impl Fn(I) -> I, a: f64, b: f64) -> Option<I> {
    let mut out: Option<I> = None;
    let mut todo = vec![(a, b)];
    while let Some((p, q)) = todo.pop() {
        let v = f(I::new(p, q));
        if v.lo.is_finite() && v.hi.is_finite() {
            out = Some(out.map_or(v, |o| I::new(o.lo.min(v.lo), o.hi.max(v.hi))));
            continue;
        }
        let mid = 0.5 * (p + q);
        let m = f(I::point(mid));
        if mid <= p || mid >= q || !(m.lo.is_finite() && m.hi.is_finite()) {
            return None;
        }
        todo.push((mid, q));
        todo.push((p, mid));
    }
    out
}

/// The most of |f| over t in [a, b] (`hull`).
pub(crate) fn most(f: impl Fn(I) -> I, a: f64, b: f64) -> Option<f64> {
    hull(f, a, b).map(interval::mag)
}

/// A moving point's bounds over t in [a, b] (`hull`; unbounded where it
/// has none).
pub(crate) fn bounds(p: &XP, a: f64, b: f64) -> [I; 2] {
    [0, 1].map(|k| hull(|t| p[k].at(Sl::of(t)).e, a, b).unwrap_or(I::ENTIRE))
}

impl Doc {
    /// An expression in t, as a field of the node t (`memo`: its parts
    /// lowered already).
    pub(super) fn lower_ex(&mut self, x: &X, t: NodeId, memo: &mut HashMap<usize, NodeId>) -> NodeId {
        if let Some(&n) = memo.get(&key(x)) {
            return n;
        }
        let n = match &**x {
            Ex::K(v) => self.g.constant(*v),
            Ex::Lin(a, b) => self.ramp(*a, *b, t),
            Ex::Bin(op, a, b) => {
                let (a, b) = (self.lower_ex(a, t, memo), self.lower_ex(b, t, memo));
                self.g.bin(*op, a, b)
            }
            Ex::Sqrt(a) => {
                let a = self.lower_ex(a, t, memo);
                self.g.un(UnOp::Sqrt, a)
            }
        };
        memo.insert(key(x), n);
        n
    }
}
