//! Interval arithmetic: sound bounds of every op over a box.
//!
//! Directed rounding: every inexact result is widened outward by whole ulps
//! (one for correctly rounded +, -, *, /, sqrt, on the side its exact
//! rounding error lies, so exact results stay exact; more for library
//! transcendentals, assumed accurate to a few ulps). Min, max, abs, negation
//! and floor are exact. So a bound holds for the real-number result of the
//! evaluated expression, not just its rounded value.

use crate::graph::{BinOp, UnOp, emax_f64};
use crate::tape::{BOTH, FIRST, SECOND};
use std::f64::consts::{FRAC_PI_2, PI, TAU};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct I {
    pub lo: f64,
    pub hi: f64,
}

impl I {
    pub const ENTIRE: I = I { lo: f64::NEG_INFINITY, hi: f64::INFINITY };

    pub fn new(lo: f64, hi: f64) -> I {
        if lo.is_nan() || hi.is_nan() { I::ENTIRE } else { I { lo, hi } }
    }

    pub fn point(v: f64) -> I {
        I::new(v, v)
    }

    pub fn contains(&self, v: f64) -> bool {
        self.lo <= v && v <= self.hi
    }

    /// The middle, in the interval and finite for finite ends: halved
    /// before the sum (0.5 (lo + hi) overflows for ends near 1e308) and
    /// clamped into it (the halvings round below 2^-1021). NaN for an
    /// interval with both ends infinite.
    pub fn mid(&self) -> f64 {
        (0.5 * self.lo + 0.5 * self.hi).clamp(self.lo, self.hi)
    }

    fn contains_zero(&self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }
}

// Directed rounding of one correctly rounded operation: r is the float
// result and e = exact - r, found exactly by an error-free transformation
// (two-sum, or fma's residual), or NaN where it is not (overflow, and
// results near underflow, whose error need not be representable). A bound
// moves one ulp out only on the side the exact result lies: an exact
// result stays exact (3.175 - 3.175 is 0, not a range around it), and an
// unknown error widens both sides. The error's sign is as good as random,
// so the side is chosen without a branch.

fn down(r: f64, e: f64) -> f64 {
    std::hint::select_unpredictable(e >= 0.0, r, r.next_down())
}

fn up(r: f64, e: f64) -> f64 {
    std::hint::select_unpredictable(e <= 0.0, r, r.next_up())
}

/// Results smaller than this may have lost bits to underflow: their
/// rounding error is not taken as exact.
const TINY: f64 = f64::MIN_POSITIVE * (1u64 << 54) as f64;

/// x + y and its rounding error (fast two-sum, the larger magnitude
/// first: exact in round-to-nearest; past overflow an infinite error of
/// the right sign, or NaN).
fn sum_err(x: f64, y: f64) -> (f64, f64) {
    let s = x + y;
    let (big, small) = if x.abs() >= y.abs() { (x, y) } else { (y, x) };
    (s, small - (s - big))
}

/// x y and its rounding error. 0 times anything is 0: as a limit (the real
/// numbers in a range are finite, and 0 x = 0) and as a value.
fn prod_err(x: f64, y: f64) -> (f64, f64) {
    if x == 0.0 || y == 0.0 {
        return (0.0, 0.0);
    }
    let p = x * y;
    (p, if p.abs() >= TINY { x.mul_add(y, -p) } else { f64::NAN })
}

/// x / y and the sign of its rounding error (as a value of that sign):
/// x - q y is exact by fma when the quotient is normal and x is far from
/// underflow, and x / y - q = (x - q y) / y.
fn quot_err(x: f64, y: f64) -> (f64, f64) {
    if x == 0.0 {
        return (x / y, 0.0);
    }
    let q = x / y;
    let ok = x.abs() >= TINY && q.abs() >= f64::MIN_POSITIVE && q.is_finite() && y.is_finite();
    (q, if ok { (-q).mul_add(y, x) * y.signum() } else { f64::NAN })
}

/// The product: the least and greatest of the four endpoint pairs'
/// nearest products, each rounded outward by its own error (every pair
/// that ties for an end: equal rounded products can err either way).
fn mul(a: I, b: I) -> I {
    // A point times a range (a constant, a zero or unit partial): two
    // products, ordered by the point's sign.
    if a.lo == a.hi || b.lo == b.hi {
        let (c, r) = if a.lo == a.hi { (a.lo, b) } else { (b.lo, a) };
        if c == 0.0 {
            return I::point(0.0);
        }
        let (x, y) = if c > 0.0 { (r.lo, r.hi) } else { (r.hi, r.lo) };
        let ((p, ep), (q, eq)) = (prod_err(c, x), prod_err(c, y));
        return I::new(down(p, ep), up(q, eq));
    }
    let f = [(a.lo, b.lo), (a.lo, b.hi), (a.hi, b.lo), (a.hi, b.hi)].map(|(x, y)| prod_err(x, y));
    if f.iter().any(|v| v.0.is_nan()) {
        return I::ENTIRE;
    }
    I::new(f.iter().map(|&(p, e)| down(p, e)).fold(f64::INFINITY, f64::min), f.iter().map(|&(p, e)| up(p, e)).fold(f64::NEG_INFINITY, f64::max))
}

/// The quotient, each end one division rounded outward. Not via the
/// reciprocal: 1 / b can overflow where a / b is fine (u / u for the
/// smallest subnormal u), and two roundings exceed one ulp. With b > 0 the
/// extremes are a.lo over b.hi (a.lo >= 0) or b.lo, and a.hi over b.lo
/// (a.hi >= 0) or b.hi; b < 0 is (-a) / (-b).
fn div(a: I, b: I) -> I {
    // x / 0 = 0.
    if b.lo == 0.0 && b.hi == 0.0 {
        return I::point(0.0);
    }
    if b.contains_zero() {
        return I::ENTIRE;
    }
    if b.hi < 0.0 {
        return div(neg_i(a), neg_i(b));
    }
    let (lo, el) = quot_err(a.lo, if a.lo >= 0.0 { b.hi } else { b.lo });
    let (hi, eh) = quot_err(a.hi, if a.hi >= 0.0 { b.lo } else { b.hi });
    I::new(down(lo, el), up(hi, eh))
}

fn abs(a: I) -> I {
    if a.lo >= 0.0 {
        a
    } else if a.hi <= 0.0 {
        I::new(-a.hi, -a.lo)
    } else {
        I::new(0.0, (-a.lo).max(a.hi))
    }
}

fn sin(a: I) -> I {
    if !(a.lo.is_finite() && a.hi.is_finite()) || a.hi - a.lo >= TAU {
        return I::new(-1.0, 1.0);
    }
    let (s0, s1) = (libm::sin(a.lo), libm::sin(a.hi));
    let mut lo = s0.min(s1);
    let mut hi = s0.max(s1);
    // Maxima at pi/2 + 2k*pi, minima at -pi/2 + 2k*pi.
    // Extremum detection with slack, so rounding can only add an extremum.
    let slack = 8.0 * (mag(a) + TAU).next_up().max(1.0) * f64::EPSILON;
    let k = ((a.lo - FRAC_PI_2 - slack) / TAU).ceil();
    if FRAC_PI_2 + TAU * k <= a.hi + slack {
        hi = 1.0;
    }
    let k = ((a.lo + FRAC_PI_2 - slack) / TAU).ceil();
    if -FRAC_PI_2 + TAU * k <= a.hi + slack {
        lo = -1.0;
    }
    I::new(lo, hi)
}

fn tan(a: I) -> I {
    if !(a.lo.is_finite() && a.hi.is_finite()) || a.hi - a.lo >= PI {
        return I::ENTIRE;
    }
    // Poles at pi/2 + k pi, detected with slack (as sin's extrema), so
    // rounding of the constants and of k pi can only add a pole.
    let slack = 8.0 * (mag(a) + PI).next_up().max(1.0) * f64::EPSILON;
    let k = ((a.lo - FRAC_PI_2 - slack) / PI).ceil();
    if FRAC_PI_2 + PI * k <= a.hi + slack {
        return I::ENTIRE;
    }
    I::new(libm::tan(a.lo), libm::tan(a.hi))
}

/// atan, as `un` bounds a library function: exact at 0, else four ulps
/// out (ulps_un). atan is no op of its own (std's is atan2 of y over 1);
/// atan2's bounds are built on it.
fn atan(a: I) -> I {
    if a.lo == 0.0 && a.hi == 0.0 {
        return I::point(0.0);
    }
    widen(I::new(libm::atan(a.lo), libm::atan(a.hi)), 4)
}

fn ln(a: I) -> I {
    // ln 0 = -inf; ln of a negative is 0.
    if a.hi < 0.0 {
        return I::point(0.0);
    }
    let f = |v: f64| if v > 0.0 { libm::log(v) } else { f64::NEG_INFINITY };
    let r = I::new(f(a.lo), f(a.hi));
    if a.lo < 0.0 { hull(r, I::point(0.0)) } else { r }
}

fn exp(a: I) -> I {
    I::new(libm::exp(a.lo), libm::exp(a.hi))
}

/// atan2 by half-plane, composed of the outward-rounded rules (division,
/// atan, subtraction), with pi/2 and pi enclosed by the floats either side
/// (the f64 constants lie below them): every step bounds the real result
/// of the one before, so the whole holds for real numbers. Not a few ulps
/// of the result added at the end: near the +x axis pi/2 - atan(x / y)
/// cancels to a small angle, and its rounding errors are ulps of pi/2
/// (2e-16), billions of ulps of an angle of 1e-10.
fn atan2(y: I, x: I) -> I {
    let half_pi = I::new(FRAC_PI_2, FRAC_PI_2.next_up());
    let atan_of = |a: I, b: I| atan(bin(BinOp::Div, a, b).0);
    if x.lo > 0.0 {
        atan_of(y, x)
    } else if y.lo > 0.0 {
        sub_o(half_pi, atan_of(x, y))
    } else if y.hi < 0.0 {
        sub_o(neg_i(half_pi), atan_of(x, y))
    } else if y.lo >= 0.0 && x.hi < 0.0 {
        // One zero: the negative x axis is +pi (Value rules), so y >= 0
        // left of the y axis stays on the upper side of the cut.
        sub_o(I::new(PI, PI.next_up()), atan_of(y, neg_i(x)))
    } else {
        I::new(-PI.next_up(), PI.next_up())
    }
}

/// atan2's branch on one side of its cut, continued across it: over y and
/// x < 0 throughout, the branch it takes where y >= 0 (`up`: pi + atan(y /
/// x)) or where y < 0 (atan(y / x) - pi), each continuous over the whole
/// region, its partials atan2's formula (d_atan2). A bound first order in
/// the scope takes one side of a cut as this (eval::forked).
pub fn atan2_branch(y: I, x: I, up: bool) -> I {
    let (pi, a) = (I::new(PI, PI.next_up()), atan(bin(BinOp::Div, y, x).0));
    if up { add_o(pi, a) } else { sub_o(a, pi) }
}

/// m^n for m >= 0 by squaring, every product rounded down (`up` false) or
/// up: a bound whatever n. A fixed allowance of ulps on the nearest result
/// is not one: the rounding errors of repeated squaring grow with n (31 ulps
/// at n = 64).
fn powi_dir(m: f64, n: u64, up: bool) -> f64 {
    let step = |x: f64, y: f64| {
        if x == 0.0 || y == 0.0 {
            return 0.0;
        }
        // The residual x y - p, exact by fma away from underflow and
        // overflow (elsewhere NaN: move a step): a product that is exact
        // stays so, and only an end the rounding moved inward moves out.
        let p = x * y;
        let e = if p.is_finite() && p >= f64::MIN_POSITIVE * 2f64.powi(54) { x.mul_add(y, -p) } else { f64::NAN };
        if up {
            if e <= 0.0 { p } else { p.next_up() }
        } else if e >= 0.0 {
            p
        } else {
            p.next_down().max(0.0)
        }
    };
    let (mut r, mut b, mut n) = (None, m, n);
    loop {
        if n & 1 == 1 {
            r = Some(r.map_or(b, |r| step(r, b)));
        }
        n >>= 1;
        if n == 0 {
            return r.unwrap_or(1.0);
        }
        b = step(b, b);
    }
}

/// Integer powers, rounded outward.
fn powi(a: I, n: i64) -> I {
    if n == 0 {
        return I::point(1.0);
    }
    if n < 0 {
        // pow(0, n < 0) is +inf (not 1 / 0, which the value rules make 0).
        let p = powi(a, -n);
        if p.lo == 0.0 && p.hi == 0.0 {
            return I::point(f64::INFINITY);
        }
        if p.contains_zero() {
            return I::ENTIRE;
        }
        return div_o(I::point(1.0), p);
    }
    let n = n as u64;
    // Odd powers keep the sign: -(|v|^n), its bounds swapped.
    let signed = |v: f64, up: bool| if v >= 0.0 { powi_dir(v, n, up) } else { -powi_dir(-v, n, !up) };
    if n.is_multiple_of(2) {
        let m = abs(a);
        I::new(powi_dir(m.lo, n, false), powi_dir(m.hi, n, true))
    } else {
        I::new(signed(a.lo, false), signed(a.hi, true))
    }
}

fn pow(a: I, b: I) -> I {
    if b.lo == b.hi && b.lo.fract() == 0.0 && b.lo.abs() < 1e9 {
        return powi(a, b.lo as i64);
    }
    // |a|^b for bases m >= 0 (each step rounded outward: ln, exp and the
    // product are not exact).
    let of = |m: I| {
        let r = widen(exp(widen(mul(b, widen(ln(m), 4)), 2)), 4);
        if m.lo == 0.0 { I::new(0.0f64.min(r.lo), r.hi) } else { r }
    };
    let pos = (a.hi >= 0.0).then(|| of(I::new(a.lo.max(0.0), a.hi)));
    if a.lo >= 0.0 {
        return pos.unwrap();
    }
    // Negative bases: 0 where the exponent is not a whole number, and
    // +-|a|^n at whole numbers n in b.
    let mut r = pos.map_or(I::point(0.0), |p| hull(p, I::point(0.0)));
    if b.lo.ceil() <= b.hi {
        let m = of(I::new(if a.hi < 0.0 { -a.hi } else { 0.0 }, -a.lo)).hi;
        r = hull(r, I::new(-m, m));
    }
    r
}

/// Widen outward by k ulps (infinities stay).
pub fn widen(i: I, k: u32) -> I {
    let (mut lo, mut hi) = (i.lo, i.hi);
    for _ in 0..k {
        lo = lo.next_down();
        hi = hi.next_up();
    }
    I::new(lo, hi)
}

/// Rounding error allowance (ulps) per op, added to a result whose every
/// end is one rounding of exact inputs: none for exact ops, one for a
/// correctly rounded one, four for one call of a library function (the
/// `libm` crate's, musl's: FreeBSD's routines, documented as under one ulp
/// for exp and ln and "nearly rounded" for sin and tan, with no bound
/// stated for atan, so four leaves room; with exp, ln and atan monotone, an
/// end is the function at an end, and sin and tan find their extrema and
/// poles with slack). An op of several roundings composes these rules
/// instead (atan2, pow; cos's shift has its own allowance, in `un`): an
/// allowance at the end covers the last rounding, not the earlier ones'
/// errors carried through it.
fn ulps_un(op: UnOp) -> u32 {
    match op {
        // Sqrt and square are rounded by their errors (un_rounded).
        UnOp::Neg | UnOp::Abs | UnOp::Floor | UnOp::Sqrt | UnOp::Square | UnOp::Rad | UnOp::Deg => 0,
        UnOp::Exp | UnOp::Ln | UnOp::Sin | UnOp::Cos | UnOp::Tan => 4,
    }
}

/// The value of a library function at the one float where it is a float:
/// sin, tan and atan of 0 are 0, cos and exp of 0 are 1, ln of 1 is 0. At
/// every other float (an algebraic number) each value is transcendental
/// (Lindemann-Weierstrass), so a point's bound there is never exact.
fn exact_at(op: UnOp, x: f64) -> Option<f64> {
    match op {
        UnOp::Sin | UnOp::Tan if x == 0.0 => Some(0.0),
        UnOp::Cos | UnOp::Exp if x == 0.0 => Some(1.0),
        UnOp::Ln if x == 1.0 => Some(0.0),
        _ => None,
    }
}

pub fn un(op: UnOp, a: I) -> I {
    // Exact results stay exact (as `bin`'s): a point on the +x axis has
    // angle 0, not a few ulps either side of a wrap at 0 (std's mod of
    // atan2, which then spans its whole period).
    if let Some(v) = exact_at(op, a.lo).filter(|_| a.lo == a.hi) {
        return I::point(v);
    }
    let r = widen(un_rounded(op, a), ulps_un(op));
    match op {
        // cos(a) is computed as sin(a + pi/2): allow for the rounded shift.
        // At each end the float sum a + FRAC_PI_2 is within half an ulp of
        // m of the real one (m is at least the sum's magnitude, so its ulp
        // is at least the sum's), and FRAC_PI_2 within 6.2e-17 of pi/2
        // (0.28 ulp of FRAC_PI_2 <= m): sin, 1-Lipschitz, moves by under
        // 0.78 ulp of m. Subtracting e rounds by half an ulp of a result
        // of magnitude at most 1 + e, within half an ulp of m: 1.28 ulps of
        // m in all, under the 4 allowed.
        UnOp::Cos => {
            let m = mag(a) + FRAC_PI_2;
            let e = 4.0 * (m.next_up() - m);
            I::new((r.lo - e).max(-1.0), (r.hi + e).min(1.0))
        }
        // Never negative: rounding must not take the bound below 0.
        UnOp::Sqrt | UnOp::Square => I::new(r.lo.max(0.0), r.hi),
        _ => r,
    }
}

fn un_rounded(op: UnOp, a: I) -> I {
    match op {
        UnOp::Neg => I::new(-a.hi, -a.lo),
        UnOp::Abs => abs(a),
        UnOp::Sqrt => {
            // sqrt(x) - r has the sign of x - r^2, exact by fma for x far
            // from underflow.
            let root = |x: f64| {
                let r = x.max(0.0).sqrt();
                (
                    r,
                    if x == 0.0 {
                        0.0
                    } else if x >= TINY {
                        (-r).mul_add(r, x)
                    } else {
                        f64::NAN
                    },
                )
            };
            let ((lo, el), (hi, eh)) = (root(a.lo), root(a.hi));
            I::new(down(lo, el), up(hi, eh))
        }
        UnOp::Square => {
            let m = abs(a);
            let ((lo, el), (hi, eh)) = (prod_err(m.lo, m.lo), prod_err(m.hi, m.hi));
            I::new(down(lo, el), up(hi, eh))
        }
        UnOp::Exp => exp(a),
        UnOp::Ln => ln(a),
        UnOp::Sin => sin(a),
        UnOp::Cos => sin(I::new(a.lo + FRAC_PI_2, a.hi + FRAC_PI_2)),
        UnOp::Tan => tan(a),
        UnOp::Floor => I::new(a.lo.floor(), a.hi.floor()),
        UnOp::Rad => mul_o(a, RAD),
        UnOp::Deg => mul_o(a, DEG),
    }
}

/// pi / 180 between adjacent floats: the nearest, below it, and the next
/// (UnOp::Rad's factor, holding both the real one and the float's).
pub const RAD: I = I { lo: 0.017453292519943295, hi: 0.0174532925199433 };
/// 180 / pi likewise: the float before the nearest, and the nearest,
/// above it (UnOp::Deg's factor).
pub const DEG: I = I { lo: 57.295779513082316, hi: 57.29577951308232 };

/// Returns the bound and, for min/max, which argument always wins.
/// Every result is rounded outward only where it is inexact (see `down`):
/// exact ones stay exact, so a field that is exactly 0 on a face bounds
/// to 0 there, a rotation's partials (a cosine times 1) stay points, and
/// opposite faces' slopes cancel exactly.
pub fn bin(op: BinOp, a: I, b: I) -> (I, u8) {
    let zero = |i: I| i.lo == 0.0 && i.hi == 0.0;
    let sum = |x: I, y: I| {
        if zero(y) {
            return x;
        }
        if zero(x) {
            return y;
        }
        let ((lo, el), (hi, eh)) = (sum_err(x.lo, y.lo), sum_err(x.hi, y.hi));
        I::new(down(lo, el), up(hi, eh))
    };
    match op {
        BinOp::Add => (sum(a, b), BOTH),
        BinOp::Sub => (sum(a, neg_i(b)), BOTH),
        BinOp::Mul => (mul(a, b), BOTH),
        BinOp::Div => (div(a, b), BOTH),
        // Decided only when strictly apart: where the ranges touch the two
        // may tie, and there the derivative is the average of both sides.
        BinOp::Min => {
            if a.hi < b.lo {
                (a, FIRST)
            } else if b.hi < a.lo {
                (b, SECOND)
            } else {
                (I::new(a.lo.min(b.lo), a.hi.min(b.hi)), BOTH)
            }
        }
        BinOp::Max => {
            if a.lo > b.hi {
                (a, FIRST)
            } else if b.lo > a.hi {
                (b, SECOND)
            } else {
                (I::new(a.lo.max(b.lo), a.hi.max(b.hi)), BOTH)
            }
        }
        BinOp::Atan2 => (atan2(a, b), BOTH),
        BinOp::Pow => (pow(a, b), BOTH),
        // Nondecreasing in both arguments; with one side <= 0 it is max,
        // exact, and otherwise one library call (hypot, which neither
        // overflows nor underflows on the way), allowed 4 ulps (`ulps_un`).
        BinOp::Emax => {
            let c = if b.hi <= 0.0 && b.hi < a.lo {
                FIRST
            } else if a.hi <= 0.0 && a.hi < b.lo {
                SECOND
            } else {
                BOTH
            };
            let hypot = |x: f64, y: f64| x > 0.0 && y > 0.0;
            let lo = emax_f64(a.lo, b.lo);
            let hi = emax_f64(a.hi, b.hi);
            let lo = if hypot(a.lo, b.lo) { widen(I::point(lo), 4).lo } else { lo };
            let hi = if hypot(a.hi, b.hi) { widen(I::point(hi), 4).hi } else { hi };
            (I::new(lo, hi), c)
        }
    }
}

pub fn select(c: I, a: I, b: I) -> (I, u8) {
    if c.hi < 0.0 {
        (a, FIRST)
    } else if c.lo >= 0.0 {
        (b, SECOND)
    } else {
        (I::new(a.lo.min(b.lo), a.hi.max(b.hi)), BOTH)
    }
}

// ---------------------------------------------------------------- interval gradients
//
// Sound bounds on a field's gradient over a box (interval forward-mode AD).
// A gradient bound is a Lipschitz bound, which turns field values into
// guaranteed distances. Discontinuities (floor across an integer, or an
// undecided select, whose branches may disagree) give unbounded derivatives.

/// A value interval with interval partials d/dx, d/dy, d/dz, and `l`, an
/// upper bound on the gradient's length carried as a scalar. The two bound
/// |grad| in different ways: the partials are tight under rotations, the
/// scalar is tight where min/max switch between branches (a rim between a
/// face and a bore: the partials' hull has length sqrt 2, the scalar 1).
/// `lipschitz` takes the smaller.
#[derive(Clone, Copy, Debug)]
pub struct G {
    pub v: I,
    pub d: [I; 3],
    pub l: f64,
    /// When the value is a sum of squares u1^2 + u2^2 + ... (and constants
    /// at least 0): the sum of the squared bounds of the u_i (infinite where
    /// one is unbounded). Its square root bounds the gradient of
    /// sqrt(value) by Cauchy-Schwarz, even where the value touches 0 and
    /// sqrt's own derivative is unbounded: the corners of exact distance
    /// fields (len2, len3) stay 1-Lipschitz-ish. None for any other value.
    pub s2: Option<f64>,
    /// With `s2`: the Gram matrix sum of grad u_i grad u_i^T, as intervals.
    /// The gradient of sqrt(value) is at most the square root of its
    /// largest eigenvalue (|sum u_i grad u_i| <= sigma_max |u|): exactly 1
    /// for len2 and len3 of coordinates, rotated or not. Its diagonal bounds
    /// each partial alone, whatever the others.
    pub m: [[I; 3]; 3],
}

/// No Gram matrix (see `G::m`).
pub const NO_GRAM: [[I; 3]; 3] = [[I { lo: 0.0, hi: 0.0 }; 3]; 3];

/// Bound on the largest eigenvalue of a symmetric matrix with entries in
/// `m`, and positive semidefinite: the smaller of Gershgorin and the trace.
fn gram_bound(m: &[[I; 3]; 3]) -> f64 {
    let gersh = (0..3).map(|i| sum_up((0..3).map(|j| if j == i { m[i][i].hi.max(0.0) } else { mag(m[i][j]) }))).fold(0.0, f64::max);
    let trace = sum_up((0..3).map(|i| m[i][i].hi.max(0.0)));
    gersh.min(trace)
}

/// The largest |x| over `a`.
pub fn mag(a: I) -> f64 {
    a.lo.abs().max(a.hi.abs())
}

fn neg_i(a: I) -> I {
    I::new(-a.hi, -a.lo)
}
/// The least interval holding both.
#[inline]
pub fn hull(a: I, b: I) -> I {
    I::new(a.lo.min(b.lo), a.hi.max(b.hi))
}

/// Two bounds of one value, met: their intersection. Both hold the value,
/// so they overlap; where they do not, one rule is broken, and which one is
/// unknown here, so the hull, which holds the value if either does.
pub fn meet(a: I, b: I) -> I {
    let r = I::new(a.lo.max(b.lo), a.hi.min(b.hi));
    if r.lo <= r.hi { r } else { hull(a, b) }
}

/// The hull of a min's or max's slopes, with room for its tie's mean as a
/// derivative formula writes it (graph::Graph::d_branch: da/2 + db/2),
/// since a derivative's range meets these bounds (eval::Slopes). Where the
/// sides are constants and params the halvings are taken in floats, and a
/// halving below 2^-1021 rounds by up to half the least subnormal: the mean
/// then lies within the least subnormal of the hull, and beyond it only at
/// an end of magnitude below 2^-1021, never across 0. So an end that is not
/// 0 moves out by the least subnormal, rounded to nearest: exact there, and
/// the end itself at any larger magnitude.
fn tie_slack(h: I) -> I {
    let tiny = f64::from_bits(1);
    I::new(if h.lo == 0.0 { 0.0 } else { h.lo - tiny }, if h.hi == 0.0 { 0.0 } else { h.hi + tiny })
}

impl G {
    /// A bound of a value that `self` or `o` bounds: their hull, the larger
    /// slope bound, no sum of squares (the two need not be one).
    pub fn hull(&self, o: &G) -> G {
        G { v: hull(self.v, o.v), d: std::array::from_fn(|k| hull(self.d[k], o.d[k])), l: self.l.max(o.l), s2: None, m: NO_GRAM }
    }

    pub fn constant(v: f64) -> G {
        // A constant c >= 0 is a (trivial) sum of squares: adding it to one
        // only lowers the slope of its square root.
        G { v: I::point(v), d: [I::point(0.0); 3], l: 0.0, s2: (v >= 0.0).then_some(0.0), m: NO_GRAM }
    }
    /// A coordinate: value interval, unit partial along `axis`.
    pub fn coordinate(v: I, axis: usize) -> G {
        G { v, d: std::array::from_fn(|k| I::point(if k == axis { 1.0 } else { 0.0 })), l: 1.0, s2: None, m: NO_GRAM }
    }
    /// A value with unknown (possibly unbounded) derivatives.
    pub fn unbounded(v: I) -> G {
        G { v, d: [I::ENTIRE; 3], l: f64::INFINITY, s2: None, m: NO_GRAM }
    }
    /// Unbounded partials (a jump, or an infinite slope), except along an
    /// axis none of `args` varies along (each one's partial exactly 0 over
    /// the box): a function of values constant along it is constant along
    /// it, so its partial is exactly 0 too (atan2 across its cut, of
    /// operands that do not move along z).
    fn entire_d(v: I, args: &[&G]) -> G {
        let flat = |k: usize| args.iter().all(|a| a.d[k].lo == 0.0 && a.d[k].hi == 0.0);
        G { d: std::array::from_fn(|k| if flat(k) { I::point(0.0) } else { I::ENTIRE }), ..G::unbounded(v) }
    }
}

pub fn g_un(op: UnOp, a: G) -> G {
    let g = g_un_nearest(op, a);
    // The operand's tightest bound: its partials or its scalar.
    let la = lipschitz(&a);
    // |f(a)|' <= sup |f'| over a's interval times a's bound.
    if op == UnOp::Sqrt
        && let Some(s2) = a.s2
    {
        // sqrt of a sum of squares: |grad| <= sqrt(sum of the squared bounds),
        // and <= sqrt(the Gram matrix's largest eigenvalue).
        let l = sqrt_up(s2.min(gram_bound(&a.m)));
        let direct = if a.v.lo > 0.0 { mul_up(div_up(0.5, a.v.lo.sqrt().next_down()), la) } else { f64::INFINITY };
        // Each partial too: |sum u_i du_i/dx_k| / |u| <= sqrt(sum (du_i/dx_k)^2),
        // the Gram matrix's diagonal (bounded even where the sum is 0).
        let d = std::array::from_fn(|k| {
            let r = sqrt_up(a.m[k][k].hi.max(0.0));
            let x = g.d[k];
            if s2 == 0.0 { I::point(0.0) } else { I::new(x.lo.max(-r), x.hi.min(r)) }
        });
        return G { v: g.v, d, l: l.min(direct), s2: None, m: NO_GRAM };
    }
    let fac = match op {
        UnOp::Neg | UnOp::Abs => 1.0,
        UnOp::Sqrt => {
            if a.v.lo > 0.0 {
                div_up(0.5, a.v.lo.sqrt().next_down())
            } else {
                f64::INFINITY
            }
        }
        UnOp::Square => 2.0 * mag(a.v),
        UnOp::Exp => mag(g.v),
        UnOp::Ln => {
            if a.v.lo > 0.0 {
                div_up(1.0, a.v.lo)
            } else {
                f64::INFINITY
            }
        }
        UnOp::Sin | UnOp::Cos => 1.0,
        UnOp::Rad => RAD.hi,
        UnOp::Deg => DEG.hi,
        UnOp::Tan => add_up(1.0, mul_up(mag(g.v), mag(g.v))),
        UnOp::Floor => {
            if a.v.lo.floor() == a.v.hi.floor() {
                0.0
            } else {
                f64::INFINITY
            }
        }
    };
    let l = if fac == 0.0 { 0.0 } else { mul_up(fac, la) };
    // A square is a one-term sum of squares: la^2.
    let s2 = (op == UnOp::Square).then(|| mul_up(la, la));
    // A square's Gram matrix: grad u grad u^T.
    let m = if op == UnOp::Square { sym(|i, j| if i == j { un(UnOp::Square, a.d[i]) } else { mul_o(a.d[i], a.d[j]) }) } else { NO_GRAM };
    G { v: g.v, d: g.d, l: if l.is_nan() { f64::INFINITY } else { l }, s2, m }
}

fn g_un_nearest(op: UnOp, a: G) -> G {
    let v = un(op, a.v);
    match d_un(op, a.v, v, a.d) {
        Some(d) => G { v, d, l: f64::INFINITY, s2: None, m: NO_GRAM },
        None => G::entire_d(v, &[&a]),
    }
}

// The derivative rules, shared by the gradient bounds (G: d/dx, d/dy, d/dz)
// and the rates over the scope (eval::eval_interval_rates: d/dt, d/dp). Each
// gives bounds on the partials of op's result from its arguments' values and
// partials over the box, or None where they are unbounded: an infinite slope,
// or a jump (floor across an integer, a division whose denominator may be 0,
// atan2 across its cut, an undecided select whose branches may disagree).
// Callers rely on None for continuity: a mean value bound needs it. Every
// step rounds outward.

/// Partials of `op(a)`, where `a` has value bounds `a` and partials `da`, and
/// the result has value bounds `v`.
pub fn d_un<const N: usize>(op: UnOp, a: I, v: I, da: [I; N]) -> Option<[I; N]> {
    let scale = |f: I| da.map(|x| mul_o(f, x));
    Some(match op {
        UnOp::Neg => da.map(neg_i),
        // At 0 itself the derivative is 0 (the average of both sides).
        UnOp::Abs => {
            if a.lo > 0.0 {
                da
            } else if a.hi < 0.0 {
                da.map(neg_i)
            } else if a.lo == 0.0 {
                da.map(|x| hull(x, I::point(0.0)))
            } else if a.hi == 0.0 {
                da.map(|x| hull(neg_i(x), I::point(0.0)))
            } else {
                da.map(|x| hull(x, neg_i(x)))
            }
        }
        UnOp::Sqrt => {
            if a.lo <= 0.0 {
                return None;
            }
            scale(div_o(I::point(0.5), v))
        }
        UnOp::Square => scale(mul_o(I::point(2.0), a)),
        UnOp::Exp => scale(v),
        UnOp::Ln => {
            if a.lo <= 0.0 {
                return None;
            }
            scale(div_o(I::point(1.0), a))
        }
        UnOp::Sin => scale(un(UnOp::Cos, a)),
        UnOp::Cos => scale(neg_i(un(UnOp::Sin, a))),
        UnOp::Tan => scale(add_o(I::point(1.0), un(UnOp::Square, v))),
        UnOp::Floor => {
            if a.lo.floor() == a.hi.floor() {
                [I::point(0.0); N]
            } else {
                return None;
            }
        }
        UnOp::Rad => scale(RAD),
        UnOp::Deg => scale(DEG),
    })
}

/// Partials of `op(a, b)` (values `a`, `b`, partials `da`, `db`), where the
/// result has value bounds `v` and `choice` is the branch bin() decided.
pub fn d_bin<const N: usize>(op: BinOp, a: I, b: I, v: I, choice: u8, da: [I; N], db: [I; N]) -> Option<[I; N]> {
    Some(match op {
        BinOp::Add => std::array::from_fn(|k| add_o(da[k], db[k])),
        BinOp::Sub => std::array::from_fn(|k| add_o(da[k], neg_i(db[k]))),
        BinOp::Mul => std::array::from_fn(|k| add_o(mul_o(da[k], b), mul_o(a, db[k]))),
        BinOp::Div => {
            if b.lo <= 0.0 && b.hi >= 0.0 {
                return None;
            }
            // (1 / b) da - ((a / b) / b) db, as graph::Graph::diff writes
            // it: each float operation of the formula is one of these on
            // values they hold, rounded outward (eval::Slopes). It loses
            // that b is one value in both terms; quotient_slope keeps it.
            let (fa, fb) = (div_o(I::point(1.0), b), div_o(div_o(a, b), b));
            std::array::from_fn(|k| {
                let r = add_o(mul_o(fa, da[k]), neg_i(mul_o(fb, db[k])));
                // Both hold the partial, so they meet.
                quotient_slope(da[k], mul_o(a, db[k]), fa).map_or(r, |q| meet(q, r))
            })
        }
        // min/max are continuous: an undecided choice has a derivative in the
        // hull (with room for the tie's mean, tie_slack).
        BinOp::Min | BinOp::Max => match choice {
            FIRST => da,
            SECOND => db,
            _ => std::array::from_fn(|k| tie_slack(hull(da[k], db[k]))),
        },
        BinOp::Atan2 => {
            if atan2_cut(a, b) {
                return None;
            }
            return d_atan2(a, b, da, db);
        }
        // A combination w_a da + w_b db with weights in [0, 1]: where both
        // sides are above 0 it is their hypot, the weights a / v and b / v
        // (their bounds' quotients hold them, as v's bounds hold the hypot
        // of any a and b in theirs). Elsewhere it may be either side's, as
        // max is.
        BinOp::Emax => match choice {
            FIRST => da,
            SECOND => db,
            _ if a.lo > 0.0 && b.lo > 0.0 && v.lo > 0.0 => {
                let unit = I::new(0.0, 1.0);
                let (wa, wb) = (meet(div_o(a, v), unit), meet(div_o(b, v), unit));
                std::array::from_fn(|k| add_o(mul_o(wa, da[k]), mul_o(wb, db[k])))
            }
            _ => {
                let w = I::new(0.0, 1.0);
                std::array::from_fn(|k| add_o(mul_o(w, da[k]), mul_o(w, db[k])))
            }
        },
        BinOp::Pow => {
            // Constant exponent: d a^n = n a^(n-1) da. Otherwise only for a > 0.
            if b.lo == b.hi && db.iter().all(|x| x.lo == 0.0 && x.hi == 0.0) {
                let n = b.lo;
                if a.lo <= 0.0 && n < 1.0 {
                    return None;
                }
                // n - 1 as a point where it is exact (Knuth's two-sum error
                // is 0: whole exponents), else rounded outward.
                let m = n - 1.0;
                let (bn, bm) = (m - n, m - (m - n));
                let m = if (n - bm) + (-1.0 - bn) == 0.0 { I::point(m) } else { sub_o(I::point(n), I::point(1.0)) };
                let (dv, _) = bin(BinOp::Pow, a, m);
                let f = mul_o(I::point(n), dv);
                da.map(|x| mul_o(f, x))
            } else if a.lo > 0.0 {
                let la = un(UnOp::Ln, a);
                std::array::from_fn(|k| mul_o(v, add_o(mul_o(db[k], la), mul_o(b, div_o(da[k], a)))))
            } else {
                return None;
            }
        }
    })
}

/// Partials of atan2(y, x) where it is continuous (off its cut, or one
/// side's branch continued across it: atan2_branch): d atan2 = (x / r2) dy
/// + (-y / r2) dx, r2 = x^2 + y^2. None where r2 may be 0.
pub fn d_atan2<const N: usize>(y: I, x: I, dy: [I; N], dx: [I; N]) -> Option<[I; N]> {
    let r2 = add_o(un(UnOp::Square, y), un(UnOp::Square, x));
    if r2.lo <= 0.0 {
        return None;
    }
    let (fy, fx) = (div_o(x, r2), div_o(neg_i(y), r2));
    Some(std::array::from_fn(|k| add_o(mul_o(fy, dy[k]), mul_o(fx, dx[k]))))
}

/// The range of d s - c s^2 over d, c and s in their bounds, rounded
/// outward: the quotient rule da / b - a db / b^2 in s = 1 / b, c = a db,
/// with b one value in both terms. It holds the partial both as a real
/// number (s = 1 / b in `s`, the bound of 1 / b, and c in a db's) and as
/// the formula's value with its parts of constants and params taken in
/// floats (eval::Slopes): those parts hold b, so b is constant along the
/// axis, db is the constant 0, and the formula is (1 / b) da with 1 / b
/// rounded or not, s = 1 / b or its float, both in `s`, and c = 0 in
/// a db's bound; the product rounded or not is inside too (its float
/// lies between the float ends of a bound of its real value).
///
/// For fixed s the value is affine in (d, c), so its extremes over the
/// box are at the four corners of (d, c); at each, a quadratic in s
/// whose extremes are at the ends of `s` or where its slope d - 2 c s is
/// 0, there d^2 / (4 c), counted unless that slope has one strict sign
/// at both ends. None where a bound is unbounded.
fn quotient_slope(d: I, c: I, s: I) -> Option<I> {
    if ![d.lo, d.hi, c.lo, c.hi, s.lo, s.hi].iter().all(|v| v.is_finite()) {
        return None;
    }
    let p = I::point;
    let mut v = Vec::with_capacity(12);
    for d in [d.lo, d.hi] {
        for c in [c.lo, c.hi] {
            // Written so a part overflows only where the value does.
            let cs = |s: f64| mul_o(p(c), p(s));
            v.extend([s.lo, s.hi].map(|s| mul_o(p(s), add_o(p(d), neg_i(cs(s))))));
            if c != 0.0 {
                let slope = |s: f64| add_o(p(d), neg_i(mul_o(p(2.0), cs(s))));
                let (g0, g1) = (slope(s.lo), slope(s.hi));
                if !(g0.lo > 0.0 && g1.lo > 0.0 || g0.hi < 0.0 && g1.hi < 0.0) {
                    v.push(mul_o(div_o(p(d), mul_o(p(4.0), p(c))), p(d)));
                }
            }
        }
    }
    // An overflow's inf - inf is ENTIRE (I::new): no bound.
    v.into_iter().reduce(hull)
}

/// Partials of select(c, a, b): the decided branch's, or with the branches
/// the same field (`same`), either's; else unbounded (they may disagree
/// where the condition changes sign).
pub fn d_select<const N: usize>(choice: u8, same: bool, da: [I; N], db: [I; N]) -> Option<[I; N]> {
    match choice {
        FIRST => Some(da),
        SECOND => Some(db),
        _ if same => Some(da),
        _ => None,
    }
}

/// Partials of component `m` of a sampled field at coordinates with bounds
/// `c` and partials `dc`, by the chain rule through the components with one
/// more axis (which jump across cell faces along an axis the component
/// already has: unbounded there, 0 within a cell). Beyond the grid's ends a
/// coordinate clamps and the slope along it is 0: where the box may reach
/// an end, 0 joins the slope's bounds (without it, slopes along two axes
/// can cancel where one of them is really 0).
pub fn d_grid<const N: usize>(g: &crate::grid::GridData, m: u8, c: [I; 3], dc: [[I; N]; 3]) -> [I; N] {
    let slope: [I; 3] = std::array::from_fn(|a| {
        if m >> a & 1 == 0 {
            let s = g.component_bound(c, m | 1 << a);
            if g.may_clamp(c[a], a) { hull(s, I::point(0.0)) } else { s }
        } else if g.crosses_face(c[a], a) {
            I::ENTIRE
        } else {
            I::point(0.0)
        }
    });
    std::array::from_fn(|k| (0..3).fold(I::point(0.0), |acc, a| add_o(acc, mul_o(slope[a], dc[a][k]))))
}

pub fn g_bin(op: BinOp, a: G, b: G) -> G {
    let (g, choice) = g_bin_nearest(op, a, b);
    let (la, lb) = (lipschitz(&a), lipschitz(&b));
    let d = g.d;
    // |b| |grad a| + |a| |grad b|, the numerator of the product and quotient rules.
    let cross = || add_up(mul_up(mag(b.v), la), mul_up(mag(a.v), lb));
    let l = match op {
        BinOp::Add | BinOp::Sub => add_up(la, lb),
        BinOp::Mul => cross(),
        BinOp::Div => {
            let bmin = if b.v.lo > 0.0 {
                b.v.lo
            } else if b.v.hi < 0.0 {
                -b.v.hi
            } else {
                0.0
            };
            if bmin > 0.0 { div_up(cross(), (bmin * bmin).next_down().max(0.0)) } else { f64::INFINITY }
        }
        // The point of the scalar: a switch between branches is no steeper
        // than the steeper branch.
        BinOp::Min | BinOp::Max => match choice {
            FIRST => la,
            SECOND => lb,
            _ => la.max(lb),
        },
        BinOp::Atan2 => {
            let r2 = (un(UnOp::Square, a.v).lo + un(UnOp::Square, b.v).lo).next_down();
            if r2 > 0.0 && !atan2_cut(a.v, b.v) { div_up(cross(), r2) } else { f64::INFINITY }
        }
        BinOp::Pow => {
            // From the partials (their length is a valid bound too).
            norm(d.map(mag))
        }
        // Where max decides, the steeper branch; where both are positive,
        // |w_a grad a + w_b grad b| over unit weights: the Gram matrix's
        // largest eigenvalue (Gershgorin), or Cauchy-Schwarz.
        BinOp::Emax => match choice {
            FIRST => la,
            SECOND => lb,
            _ if a.v.hi <= 0.0 || b.v.hi <= 0.0 => la.max(lb),
            _ => {
                let (gaa, gbb, gab) = (dot_o(&a.d, &a.d).hi.min(mul_up(la, la)), dot_o(&b.d, &b.d).hi.min(mul_up(lb, lb)), mag(dot_o(&a.d, &b.d)));
                let round = sqrt_up(add_up(gaa.max(gbb), gab)).min(sqrt_up(add_up(mul_up(la, la), mul_up(lb, lb))));
                la.max(lb).max(round)
            }
        },
    };
    // Sums of sums of squares stay sums of squares.
    let s2 = if op == BinOp::Add { a.s2.zip(b.s2).map(|(p, q)| add_up(p, q)) } else { None };
    let m = if s2.is_some() { sym(|i, j| add_o(a.m[i][j], b.m[i][j])) } else { NO_GRAM };
    G { v: g.v, d, l: if l.is_nan() { f64::INFINITY } else { l }, s2, m }
}

fn g_bin_nearest(op: BinOp, a: G, b: G) -> (G, u8) {
    let (v, choice) = bin(op, a.v, b.v);
    let g = match d_bin(op, a.v, b.v, v, choice, a.d, b.d) {
        Some(d) => G { v, d, l: f64::INFINITY, s2: None, m: NO_GRAM },
        None => G::entire_d(v, &[&a, &b]),
    };
    (g, choice)
}

/// Component `m` of a sampled field (graph::Node::Grid) at coordinates `c`:
/// its range, and its partials by the chain rule through the components
/// with one more axis (which jump across cell faces along an axis the
/// component already has: unbounded there, 0 within a cell).
pub fn g_grid(g: &crate::grid::GridData, m: u8, c: [G; 3]) -> G {
    let b = [c[0].v, c[1].v, c[2].v];
    let v = g.component_bound(b, m);
    let d = d_grid(g, m, b, [c[0].d, c[1].d, c[2].d]);
    let l = norm(d.map(mag));
    G { v, d, l: if l.is_nan() { f64::INFINITY } else { l }, s2: None, m: NO_GRAM }
}

/// Whether atan2(y, x) jumps by 2 pi inside the box: across the negative x
/// axis (y below 0 and reaching it where x < 0; on the axis it is +pi, the
/// upper side's).
fn atan2_cut(y: I, x: I) -> bool {
    y.lo < 0.0 && y.hi >= 0.0 && x.lo < 0.0
}

/// `same`: the branches are the same field (perhaps under different names),
/// so the result is continuous where the condition changes sign.
pub fn g_select(c: G, a: G, b: G, same: bool) -> G {
    let (v, choice) = select(c.v, a.v, b.v);
    match choice {
        FIRST => G { v, d: a.d, l: a.l, s2: a.s2, m: a.m },
        SECOND => G { v, d: b.d, l: b.l, s2: b.s2, m: b.m },
        // Otherwise the branches may disagree where the condition changes
        // sign (equal bounds do not mean equal values: min and max of the
        // same arguments have the same bounds).
        _ if same => G { v, ..a },
        _ => G::entire_d(v, &[&c, &a, &b]),
    }
}

/// The coordinates of a Lipschitz promise (graph::Node::Frame). The value is
/// 0; `l` carries the frame's stretch: a bound on the largest singular value
/// of its Jacobian J, the smaller of Frobenius (from the rows' bounds) and
/// the square root of J J^T's largest eigenvalue (gram_bound, from the
/// partials: 1 for any rigid motion).
pub fn g_frame(c: [G; 3]) -> G {
    let frobenius = norm(std::array::from_fn(|k| lipschitz(&c[k])));
    let jjt = sym(|i, j| dot_o(&c[i].d, &c[j].d));
    let l = frobenius.min(sqrt_up(gram_bound(&jjt)));
    G { v: I::point(0.0), d: [I::point(0.0); 3], l: if l.is_nan() { f64::INFINITY } else { l }, s2: None, m: NO_GRAM }
}

/// `a` promised k-Lipschitz in `frame`'s coordinates: its gradient is at most
/// k times the frame's stretch. A k that is not positive promises nothing (a
/// library's promise for only some of its arguments).
pub fn g_lipschitz(a: G, k: G, frame: G) -> G {
    let promised = if k.v.lo > 0.0 { mul_up(k.v.hi, frame.l) } else { f64::INFINITY };
    G { l: lipschitz(&a).min(if promised.is_nan() { f64::INFINITY } else { promised }), ..a }
}

// Upper bounds of nonnegative quantities (slopes, radii, norms), rounded up
// at every step. A nearest result is within half a step of the exact one,
// subnormal or not, so the next float up is at least the exact result: a
// positive product never rounds to 0, which a relative margin (x (1 + e))
// cannot repair. A zero operand is exact (0 times infinity stays NaN: no
// bound).

/// a b, rounded up.
pub fn mul_up(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 { a * b } else { (a * b).next_up() }
}

/// a / b, rounded up (b >= 0, itself rounded down).
fn div_up(a: f64, b: f64) -> f64 {
    if a == 0.0 { 0.0 } else { (a / b).next_up() }
}

/// a + b, rounded up.
pub fn add_up(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 { a + b } else { (a + b).next_up() }
}

/// A sum rounded up: at least the exact sum of the terms.
pub fn sum_up(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(0.0, add_up)
}

/// The square root, rounded up.
pub fn sqrt_up(a: f64) -> f64 {
    if a == 0.0 { 0.0 } else { a.sqrt().next_up() }
}

/// A float s with s |v| <= 1, proven, |v| the Euclidean norm of a vector
/// whose components are at most `v`'s in size: a reciprocal length rounded
/// down, the scale that keeps a field along that direction (a line through
/// its normal, a leaf whose slopes are bounded by v) at most 1-Lipschitz in
/// real numbers, where a quotient or a norm rounded to nearest can leave it
/// an ulp over. The norm is rounded up (scaled by the largest component, so
/// squares neither underflow nor overflow), and s = 1 / norm is checked by
/// an exact fma; so s is exactly 1 / |v| where that is a float (one nonzero
/// component, or an exact norm). Infinite for the zero vector.
pub fn inv_norm_down(v: &[f64]) -> f64 {
    let nz: Vec<f64> = v.iter().map(|x| x.abs()).filter(|&x| x != 0.0).collect();
    let l = match nz.as_slice() {
        [] => return f64::INFINITY,
        [one] => *one,
        _ => {
            let m = nz.iter().fold(0.0f64, |a, &x| a.max(x));
            mul_up(
                m,
                sqrt_up(sum_up(nz.iter().map(|&x| {
                    // x / m rounded up: its square at least the exact one.
                    let q = if x == m { 1.0 } else { (x / m).next_up() };
                    mul_up(q, q)
                }))),
            )
        }
    };
    let s = 1.0 / l;
    // s l - 1, exactly rounded: above 0 only where s l > 1.
    if s.mul_add(l, -1.0) > 0.0 { s.next_down() } else { s }
}

/// Interval sum, product, quotient and dot product, rounded outward where
/// inexact (as in bin).
pub fn add_o(a: I, b: I) -> I {
    bin(BinOp::Add, a, b).0
}
pub fn sub_o(a: I, b: I) -> I {
    bin(BinOp::Sub, a, b).0
}
pub fn mul_o(a: I, b: I) -> I {
    bin(BinOp::Mul, a, b).0
}
pub fn div_o(a: I, b: I) -> I {
    bin(BinOp::Div, a, b).0
}
fn dot_o(x: &[I; 3], y: &[I; 3]) -> I {
    (0..3).fold(I::point(0.0), |acc, k| add_o(acc, mul_o(x[k], y[k])))
}

/// Whether every vector in box `x` points away from every vector in box
/// `y` (no two of them, one from each, are parallel and alike): their dot
/// product is below the product of their least lengths throughout. False
/// where either box holds 0.
pub fn apart(x: &[I; 3], y: &[I; 3]) -> bool {
    // Each box scaled by a power of 2 (exact) to its largest bound near 1:
    // directions are kept, and products neither underflow nor overflow.
    let scaled = |v: &[I; 3]| -> [I; 3] {
        let m = v.iter().map(|i| i.lo.abs().max(i.hi.abs())).fold(0.0, f64::max);
        if !(m > 0.0 && m.is_finite()) {
            return *v;
        }
        // In two steps: 2^e itself may not be a float (e up to 1074).
        let e = -m.log2().round() as i32;
        let (a, b) = (2f64.powi(e / 2), 2f64.powi(e - e / 2));
        v.map(|i| I::new(i.lo * a * b, i.hi * a * b))
    };
    let (x, y) = (&scaled(x), &scaled(y));
    let least = |v: &[I; 3]| {
        v.iter()
            .map(|i| {
                if i.lo > 0.0 {
                    i.lo
                } else if i.hi < 0.0 {
                    -i.hi
                } else {
                    0.0
                }
            })
            .map(|m| m * m)
            .sum::<f64>()
    };
    let (lx, ly, d) = (least(x), least(y), dot_o(x, y).hi);
    lx > 0.0 && ly > 0.0 && (d < 0.0 || d * d < lx * ly * (1.0 - 1e-12))
}

/// A symmetric matrix from its entries on and above the diagonal.
fn sym(f: impl Fn(usize, usize) -> I) -> [[I; 3]; 3] {
    let (a, b, c) = (f(0, 1), f(0, 2), f(1, 2));
    [[f(0, 0), a, b], [a, f(1, 1), c], [b, c, f(2, 2)]]
}

/// Largest gradient magnitude allowed by the bounds (a Lipschitz constant).
pub fn lipschitz(g: &G) -> f64 {
    norm(g.d.map(mag)).min(g.l)
}

/// Euclidean norm, rounded up, and scaled by the largest term m so squares
/// neither underflow (a slope of 1e-200 is not 0) nor overflow. The scaled
/// sum is at least 1, so its rounding errors are relative (a quotient that
/// underflows loses less than one ulp of it), and 8 ulps cover them; the
/// final rescale by m, which can fall below the normal range, rounds up.
/// Three terms: the 8 ulps are not enough for a long sum.
pub fn norm(v: [f64; 3]) -> f64 {
    let m = v.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    if m == 0.0 || !m.is_finite() {
        return m;
    }
    mul_up(m, v.iter().map(|x| (x / m).powi(2)).sum::<f64>().sqrt() * (1.0 + 8.0 * f64::EPSILON))
}

#[cfg(test)]
mod tests {

    #[test]
    fn inv_norm_down_is_at_most_the_reciprocal_length() {
        // s^2 (x^2 + y^2) <= 1 in integers, s = m 2^e: m^2 n <= 2^(-2e).
        for (x, y) in (1..60u64).flat_map(|x| (0..60u64).map(move |y| (x, y))) {
            let s = inv_norm_down(&[x as f64, y as f64]);
            let (m, e) = {
                let b = s.to_bits();
                ((b & ((1 << 52) - 1)) | (1 << 52), ((b >> 52) & 0x7ff) as i32 - 1075)
            };
            let n = (x * x + y * y) as u128;
            assert!((m as u128) * (m as u128) * n <= 1u128 << (-2 * e), "({x}, {y}): {s}");
        }
        assert_eq!(inv_norm_down(&[0.0, 4.0]), 0.25);
        assert_eq!(inv_norm_down(&[0.0, 0.0]), f64::INFINITY);
    }

    use super::*;

    #[test]
    fn pi_constants_lie_below_pi() {
        // atan2 encloses pi/2 and pi by FRAC_PI_2 and PI and the floats
        // above them: sin(PI) = sin(pi - PI) > 0 says PI < pi (by 1.2e-16,
        // far beyond sin's error there), and FRAC_PI_2 is PI halved, exactly.
        assert!(PI.sin() > 0.0 && FRAC_PI_2 * 2.0 == PI);
    }

    #[test]
    fn upward_helpers_bound_subnormal_results() {
        // Exact results known: each rounded-up result is at least it, where
        // nearest rounding (and a relative margin) would fall below.
        let u = f64::from_bits(1);
        assert!(mul_up(u, 0.5) >= u, "u/2 rounds to 0");
        assert!(mul_up(1e-200, 1e-130) > 0.0, "1e-330 underflows");
        assert!(mul_up(3.0 * u, 0.5) >= 2.0 * u, "1.5 u rounds to 2 u (even), the next up 3 u");
        // The clearance claim's offset L d (a relative margin would give 0 and u).
        assert!(mul_up(u, 1.25) >= 2.0 * u && mul_up(0.0, 1.0) == 0.0);
        assert!(add_up(0.1, 0.2) > 0.1 + 0.2);
        // norm [u, u, 0] = sqrt(2) u: the next float up from u is 2 u.
        assert!(norm([u, u, 0.0]) >= 2.0 * u);
        assert!(norm([3.0 * u, 4.0 * u, 0.0]) >= 5.0 * u);
        assert!(norm([1e-200, 1e-200, 1e-200]) >= 1e-200 * 3f64.sqrt());
        // Exact zeros stay 0 (a constant's slope), and 0 x infinity is no bound.
        assert_eq!((mul_up(0.0, 5.0), add_up(0.0, 0.0), sum_up([0.0; 3].into_iter()), norm([0.0; 3])), (0.0, 0.0, 0.0, 0.0));
        assert!(mul_up(0.0, f64::INFINITY).is_nan());
    }

    #[test]
    fn product_partials_survive_cancellation() {
        // a = (1 + e) x - 1, b = x + 1 - e at x = 0: (ab)' = (1 + e)(1 - e) - 1
        // = -e^2, which a sum of nearest products cancels to 0.
        let e = 2f64.powi(-27);
        let z = I::point(0.0);
        let a = G { d: [I::point(1.0 + e), z, z], ..G::unbounded(I::point(-1.0)) };
        let b = G { d: [I::point(1.0), z, z], ..G::unbounded(I::point(1.0 - e)) };
        let r = g_bin(BinOp::Mul, a, b);
        let exact = -2f64.powi(-54);
        assert!(r.d[0].contains(exact), "{r:?}");
        assert!(lipschitz(&r) >= exact.abs(), "{r:?}");
        // The quotient a / b at a = x, x = 0, b = 2^600: b^2 overflows, the
        // partial is 2^-600.
        let a = G::coordinate(I::point(0.0), 0);
        let b = G::constant(2f64.powi(600));
        let r = g_bin(BinOp::Div, a, b);
        assert!(r.d[0].contains(2f64.powi(-600)), "{r:?}");
    }

    /// x = m 2^e with m in [0.5, 1) (0 for 0), for finite x.
    fn frexp(x: f64) -> (f64, i32) {
        if x == 0.0 {
            return (0.0, 0);
        }
        let (x, k) = if x.abs() < f64::MIN_POSITIVE { (x * 2f64.powi(64), -64) } else { (x, 0) };
        let bits = x.to_bits();
        (f64::from_bits(bits & !(0x7ff << 52) | 1022 << 52), ((bits >> 52) & 0x7ff) as i32 - 1022 + k)
    }

    /// p q against n in real numbers, for finite p, q, n. p and q are
    /// scaled to [0.5, 1) first, so the error of their product (which a
    /// fused multiply-add gives exactly) cannot underflow, and the product
    /// cannot overflow or underflow, however large or small p q is.
    fn cmp_product(p: f64, q: f64, n: f64) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        let ((mp, ep), (mq, eq), (mn, en)) = (frexp(p), frexp(q), frexp(n));
        let m = mp * mq;
        if m == 0.0 || n == 0.0 {
            return m.partial_cmp(&0.0).unwrap().then(0.0f64.partial_cmp(&n).unwrap());
        }
        if (m < 0.0) != (n < 0.0) {
            return m.partial_cmp(&0.0).unwrap();
        }
        // Same sign: |m| in (0.25, 1), |mn| in [0.5, 1), so a gap of more
        // than 2 in the exponents decides.
        let (e, k, s) = (mp.mul_add(mq, -m), ep + eq - en, if m < 0.0 { Less } else { Greater });
        if k.abs() > 2 {
            return if k > 0 { s } else { s.reverse() };
        }
        // The error is under half an ulp of m, so m against the float t
        // decides unless they are equal.
        let t = mn * 2f64.powi(-k);
        m.partial_cmp(&t).unwrap().then(e.partial_cmp(&0.0).unwrap())
    }

    /// x <= p q <= y in real numbers (an infinite end holds).
    fn holds_product(x: f64, y: f64, p: f64, q: f64) -> bool {
        use std::cmp::Ordering::*;
        (x == f64::NEG_INFINITY || cmp_product(p, q, x) != Less) && (y == f64::INFINITY || cmp_product(p, q, y) != Greater)
    }

    #[test]
    fn products_compare_exactly_at_any_scale() {
        use std::cmp::Ordering::*;
        // The product underflows to 0: 2^-1200 > 0.
        assert_eq!(cmp_product(2f64.powi(-300), 2f64.powi(-900), 0.0), Greater);
        assert!(!holds_product(0.0, 0.0, 2f64.powi(-300), 2f64.powi(-900)));
        // Overflows: 2^1200 > MAX; and an exact product of subnormals' scale.
        assert_eq!(cmp_product(2f64.powi(600), 2f64.powi(600), f64::MAX), Greater);
        assert_eq!(cmp_product(3.0 * 2f64.powi(-600), 2f64.powi(-474), 3.0 * 2f64.powi(-1074)), Equal);
        // (1 + 2^-52)^2 = 1 + 2^-51 + 2^-104 rounds to 1 + 2^-51, below it.
        let a = 1.0 + f64::EPSILON;
        assert_eq!(cmp_product(a, a, a * a), Greater);
        assert_eq!(cmp_product(-a, a, -(a * a)), Less);
        assert_eq!(cmp_product(-2.0, 3.0, -6.0), Equal);
        assert_eq!(cmp_product(0.0, 3.0, -1.0), Greater);
    }

    #[test]
    fn quotient_slopes_hold_the_partial() {
        // Random bounds of integers scaled by powers of two, so that at
        // points inside, n = da b - a db and d = b^2 are exact floats (a
        // db and da b share a scale): the partial n / d lies in the rule's
        // bound, checked exactly at any scale (r.lo d <= n <= r.hi d, cmp_product). With db = 0, so
        // does (1 / b) da with 1 / b rounded, its product real or rounded.
        // Tighter than the formula's own bounds and (da b - a db) / b^2, but
        // for their rounding, and often by much.
        let seed = std::cell::Cell::new(0x9e3779b97f4a7c15u64);
        let next = || {
            seed.set(seed.get().wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407));
            seed.get() >> 11
        };
        let (mut tighter, cases) = (0, 20000);
        for _ in 0..cases {
            let int = |m: u64| (next() % (2 * m + 1)) as f64 - m as f64;
            let range = |m: u64| {
                let (x, y) = (int(m), int(m));
                (x.min(y), x.max(y))
            };
            // Narrow ranges (near points) and wide ones.
            let m = if next() % 3 == 0 { 3 } else { 1 << 20 };
            let (a, da, db) = (range(m), range(m), range(m));
            let (x, y) = (int(m).abs() + 1.0, int(m).abs() + 1.0);
            let b = if next() % 2 == 0 { (x.min(y), x.max(y)) } else { (-x.max(y), -x.min(y)) };
            let k = || (next() % 601) as i32 - 300;
            let (ka, kb, kdb) = (k(), k(), k());
            let kda = ka + kdb - kb;
            let sc = |v: f64, k: i32| v * 2f64.powi(k);
            let iv = |(lo, hi): (f64, f64), k| I::new(sc(lo, k), sc(hi, k));
            let (ai, bi, dai, dbi) = (iv(a, ka), iv(b, kb), iv(da, kda), iv(db, kdb));
            let r = d_bin(BinOp::Div, ai, bi, I::ENTIRE, BOTH, [dai], [dbi]).unwrap()[0];
            // The rule, and quotient_slope's bound alone (not met with the
            // formula's).
            let q = quotient_slope(dai, mul_o(ai, dbi), div_o(I::point(1.0), bi)).unwrap_or(I::ENTIRE);
            let q0 = quotient_slope(dai, I::point(0.0), div_o(I::point(1.0), bi)).unwrap_or(I::ENTIRE);
            for _ in 0..8 {
                let pick = |(lo, hi): (f64, f64)| lo + (next() % ((hi - lo) as u64 + 1)) as f64;
                let (x, y, dx, dy) = (pick(a), pick(b), pick(da), pick(db));
                let (n, d) = (sc(dx * y - x * dy, ka + kdb), sc(y * y, 2 * kb));
                for r in [r, q] {
                    let lo = r.lo == f64::NEG_INFINITY || cmp_product(r.lo, d, n) != std::cmp::Ordering::Greater;
                    let hi = r.hi == f64::INFINITY || cmp_product(r.hi, d, n) != std::cmp::Ordering::Less;
                    assert!(lo && hi, "a {ai:?} b {bi:?} da {dai:?} db {dbi:?}: {r:?} misses {n} / {d}");
                }
            }
            let r0 = d_bin(BinOp::Div, ai, bi, I::ENTIRE, BOTH, [dai], [I::point(0.0)]).unwrap()[0];
            for (y, dx) in [(b.0, da.0), (b.1, da.1), (b.0, da.1)] {
                let (f, dx) = (1.0 / sc(y, kb), sc(dx, kda));
                for r0 in [r0, q0] {
                    assert!(holds_product(r0.lo, r0.hi, f, dx) && r0.lo <= f * dx && f * dx <= r0.hi, "{r0:?} misses {f} {dx}");
                }
            }
            let formula = add_o(mul_o(div_o(I::point(1.0), bi), dai), neg_i(mul_o(div_o(div_o(ai, bi), bi), dbi)));
            let old = div_o(add_o(mul_o(dai, bi), neg_i(mul_o(ai, dbi))), mul_o(bi, bi));
            for o in [formula, old] {
                let slack = 1e-12 * o.lo.abs().max(o.hi.abs());
                assert!(r.lo >= o.lo - slack && r.hi <= o.hi + slack, "a {ai:?} b {bi:?} da {dai:?} db {dbi:?}: {r:?} looser than {o:?}");
            }
            if r.hi - r.lo < 0.9 * (formula.hi - formula.lo).min(old.hi - old.lo) {
                tighter += 1;
            }
        }
        assert!(tighter > cases / 10, "{tighter} of {cases} tighter by 10%");
        // a = b = [1, 2], da = db = 1: over the box, da / b - a db / b^2
        // spans exactly [-1, 1/4].
        let (ab, one) = (I::new(1.0, 2.0), [I::point(1.0)]);
        let r = d_bin(BinOp::Div, ab, ab, I::ENTIRE, BOTH, one, one).unwrap()[0];
        assert!(r.lo <= -1.0 && r.lo > -1.0 - 1e-12 && r.hi >= 0.25 && r.hi < 0.25 + 1e-12, "{r:?}");
        // With b constant along the axis (db = 0) the formula is (1 / b) da
        // with 1 / b a float where b is constants and params: at b = 3,
        // da = 3 that is 1 - 2^-54, below the real 1.
        let r = d_bin(BinOp::Div, I::new(0.0, 1.0), I::point(3.0), I::ENTIRE, BOTH, [I::point(3.0)], [I::point(0.0)]).unwrap()[0];
        assert!(holds_product(r.lo, r.hi, 1.0 / 3.0, 3.0) && r.hi >= 1.0 && r.hi - r.lo < 1e-15, "{r:?}");
    }

    #[test]
    fn division_handles_subnormal_denominators() {
        let u = f64::from_bits(1);
        for a in [u, -u] {
            let r = bin(BinOp::Div, I::point(a), I::point(u)).0;
            assert!(r.contains(if a > 0.0 { 1.0 } else { -1.0 }), "{r:?}");
            let r = bin(BinOp::Div, I::point(1.0), I::point(a)).0;
            assert!(r.lo <= 1.0 / a && 1.0 / a <= r.hi, "{r:?}");
        }
        // The derivative rules' quotient too (d ln a = da / a).
        let r = d_un(UnOp::Ln, I::point(u), un(UnOp::Ln, I::point(u)), [I::point(u)]).unwrap();
        assert!(r[0].contains(1.0), "{r:?}");
        // x / 0 = 0; a denominator touching 0 gives no bound.
        assert!(bin(BinOp::Div, I::point(1.0), I::point(0.0)).0.contains(0.0));
        for b in [I::new(0.0, u), I::new(-u, 0.0), I::new(-u, u)] {
            assert_eq!(bin(BinOp::Div, I::point(1.0), b).0, I::ENTIRE);
        }
        assert!(bin(BinOp::Div, I::point(u), I::point(2.0)).0.hi >= u);
        // Unbounded ranges: a / b over a, b >= 1 is (0, inf).
        let r = bin(BinOp::Div, I::new(1.0, f64::INFINITY), I::new(1.0, f64::INFINITY)).0;
        assert!(r.lo <= 0.0 && r.hi == f64::INFINITY && r.lo > -1e-300, "{r:?}");
        let r = bin(BinOp::Div, I::new(-3.0, 2.0), I::new(-4.0, -0.5)).0;
        assert!(r.lo <= -4.0 && r.lo > -4.0001 && r.hi >= 6.0 && r.hi < 6.0001, "{r:?}");
    }

    #[test]
    fn power_slopes_round_the_exponent_outward() {
        // d a^n / da = n a^(n - 1) with n = 2^-54: n - 1 rounds to -1, and
        // a^(n - 1) at a = 2^500 is 2^-500 (1 + 1.9e-14), 87 ulps above a^-1.
        let n = 2f64.powi(-54);
        let a = I::point(2f64.powi(500));
        let v = bin(BinOp::Pow, a, I::point(n)).0;
        let r = d_bin(BinOp::Pow, a, I::point(n), v, BOTH, [I::point(1.0)], [I::point(0.0)]).unwrap();
        assert!(r[0].hi >= 2f64.powi(-554) * (1.0 + 1.9e-14), "{r:?}");
        // A whole exponent keeps the exact power rule (x^3: 3 x^2).
        let r = d_bin(BinOp::Pow, I::new(-2.0, 1.0), I::point(3.0), I::new(-8.0, 1.0), BOTH, [I::point(1.0)], [I::point(0.0)]).unwrap();
        assert!(r[0].lo > -1e-12 && r[0].hi < 12.0 + 1e-12, "{r:?}");
    }

    #[test]
    fn integer_powers_round_outward_for_any_exponent() {
        // x^n by repeated squaring, nearest, is off by 31 ulps at n = 64 and
        // 731 at n = 1000 for these x (exact values from rational arithmetic:
        // between the two doubles given).
        for (x, n, lo, hi) in [
            (1.0020677346905054, 64.0, 1.1413347280818533, 1.1413347280818535),
            (1.0036629244665904, 100.0, 1.4414120588965522, 1.4414120588965524),
            (1.000654916230334, 1000.0, 1.9245686530402346, 1.9245686530402348),
        ] {
            let r = bin(BinOp::Pow, I::point(x), I::point(n)).0;
            assert!(r.lo <= lo && hi <= r.hi && r.hi - r.lo < 1e-12, "{x}^{n}: {r:?}");
            // Odd and negative exponents: the same bounds, mirrored and inverted.
            let r = bin(BinOp::Pow, I::point(-x), I::point(n + 1.0)).0;
            assert!(r.lo <= -x * hi && -x * lo <= r.hi && r.hi < 0.0, "{r:?}");
            let r = bin(BinOp::Pow, I::point(x), I::point(-n)).0;
            assert!(r.lo <= 1.0 / hi && 1.0 / lo <= r.hi && r.lo > 0.0, "{r:?}");
        }
        // Small powers stay exact where they are.
        assert_eq!(bin(BinOp::Pow, I::new(-3.0, 2.0), I::point(2.0)).0, I::new(0.0, 9.0));
        assert_eq!(bin(BinOp::Pow, I::new(-3.0, 2.0), I::point(3.0)).0, I::new(-27.0, 8.0));
        assert_eq!(bin(BinOp::Pow, I::point(0.0), I::point(-2.0)).0, I::point(f64::INFINITY));
    }

    #[test]
    fn tangent_poles_are_found_through_rounding() {
        // Adjacent doubles around the poles 22.5 pi and 332.5 pi (from
        // rational arithmetic with 100 digits of pi): the rounded pole
        // FRAC_PI_2 + k PI falls outside them, so tan's ends may come out
        // with opposite signs, an inverted range.
        for (a, b) in [(70.68583470577035, 70.68583470577036), (1044.5795573186062, 1044.5795573186065)] {
            assert_eq!(b, f64::next_up(a));
            assert_eq!(un(UnOp::Tan, I::new(a, b)), I::ENTIRE, "{a} {b}");
        }
        for p in [FRAC_PI_2, 1000.0 * PI + FRAC_PI_2, -FRAC_PI_2] {
            for (a, b) in [(p, p.next_up()), (p.next_down(), p), (p.next_down().next_down(), p.next_up().next_up())] {
                assert_eq!(un(UnOp::Tan, I::new(a, b)), I::ENTIRE, "{a} {b}");
            }
        }
        let r = un(UnOp::Tan, I::new(-1.0, 1.0));
        assert!(r.lo < -1.557 && r.hi > 1.557 && r.hi < 1.5575, "{r:?}");
    }

    #[test]
    fn grid_slopes_vanish_where_coordinates_clamp() {
        // value(i, j, k) = i + j on a 2x2x2 unit grid, at x = p - 1, y = 0.5 - p,
        // z = 0.5 near p = 0: x clamps to 0, so the value is 0.5 - p, slope -1.
        let values = (0..8).map(|n| ((n & 1) + (n >> 1 & 1)) as f32).collect();
        let g = crate::grid::GridData::new([0.0; 3], 1.0, [2, 2, 2], values);
        let (z, one) = (I::point(0.0), I::point(1.0));
        let cx = G { v: I::new(-1.1, -0.9), d: [one, z, z], l: 1.0, s2: None, m: NO_GRAM };
        let cy = G { v: I::new(0.4, 0.6), d: [I::point(-1.0), z, z], l: 1.0, s2: None, m: NO_GRAM };
        let r = g_grid(&g, 0, [cx, cy, G::constant(0.5)]);
        assert!(r.v.contains(0.5), "{r:?}");
        assert!(r.d[0].contains(-1.0), "{r:?}");
        assert!(lipschitz(&r) >= 1.0, "{r:?}");
    }

    #[test]
    fn exact_point_results_stay_points() {
        // A coefficient times 1 or -1, and 0 minus one, are exact: so a
        // rotation's partials stay exact and opposite faces' cancel.
        let c = (30f64).to_radians().cos();
        assert_eq!(bin(BinOp::Mul, I::point(c), I::point(1.0)).0, I::point(c));
        assert_eq!(bin(BinOp::Mul, I::point(-1.0), I::new(c, 2.0)).0, I::new(-2.0, -c));
        assert_eq!(bin(BinOp::Sub, I::point(0.0), I::point(c)).0, I::point(-c));
        // Others are widened past the real result: 0.1 + 0.2, 0.1 x 3, and
        // products that underflow.
        // The real sum of the floats 0.1 and 0.2, 0.30000000000000001665...,
        // lies strictly between the floats 0.3 and 0.1 + 0.2 (rounded to
        // nearest, 0.30000000000000004): a bound must hold both.
        let r = bin(BinOp::Add, I::point(0.1), I::point(0.2)).0;
        assert!(r.lo <= 0.3 && 0.1 + 0.2 <= r.hi, "{r:?}");
        let r = bin(BinOp::Mul, I::point(0.1), I::point(3.0)).0;
        assert!(r.lo < r.hi, "{r:?}");
        let u = f64::from_bits(1);
        for (x, y) in [(u, 0.5), (3.0 * u, 0.5), (1e-200, 1e-150), (1e-160, 1e-160)] {
            // Not a point, and around the nearest product (the real one is
            // within half an ulp of it).
            let r = bin(BinOp::Mul, I::point(x), I::point(y)).0;
            assert!(r.lo < x * y && x * y < r.hi, "{x:e} {y:e}: {r:?}");
        }
        // The pair bound's combination is exact only where every step is.
        let comb = |x: f64, y: f64| add_o(mul_o(I::point(0.5), I::point(x)), mul_o(I::point(0.5), I::point(y)));
        assert_eq!(comb(c, -c), I::point(0.0));
        assert!(comb(0.1, 0.2).lo < comb(0.1, 0.2).hi);
        assert!(comb(u, 0.0).lo < comb(u, 0.0).hi);
        // Any exact result stays a point, wherever it comes from: a face's
        // field on its own plane (3.175 - 3.175) is 0, not a range around
        // it, so a ray along the face is not inside it by rounding.
        assert_eq!(bin(BinOp::Sub, I::point(3.175), I::point(3.175)).0, I::point(0.0));
        assert_eq!(bin(BinOp::Div, I::point(1.0), I::point(4.0)).0, I::point(0.25));
        assert_eq!(un(UnOp::Sqrt, I::point(2.25)), I::point(1.5));
        assert_eq!(un(UnOp::Square, I::point(-1.5)), I::point(2.25));
        assert_eq!(bin(BinOp::Emax, I::point(-1.0), I::point(0.0)).0, I::point(0.0));
        // Inexact ends move one ulp, on the side of the real result.
        let third = bin(BinOp::Div, I::point(1.0), I::point(3.0)).0;
        assert!(third.lo < 1.0 / 3.0 || third.hi > 1.0 / 3.0);
        assert!(third.hi == third.lo.next_up(), "{third:?}");
        let r = un(UnOp::Sqrt, I::point(2.0));
        assert!(r.lo * r.lo < 2.0 && r.hi * r.hi > 2.0 && r.hi == r.lo.next_up(), "{r:?}");
        // Exact ops stay exact; sin 1 is widened, tightly.
        let (m, _) = bin(BinOp::Min, I::point(1.5), I::point(2.5));
        assert_eq!((m.lo, m.hi), (1.5, 1.5));
        let s = un(UnOp::Sin, I::point(1.0));
        assert!(s.lo < 1f64.sin() && 1f64.sin() < s.hi && s.hi - s.lo < 1e-14);
        // A library function at its one float with a float value.
        let at = |op, x| un(op, I::point(x));
        let exact = [at(UnOp::Sin, 0.0), at(UnOp::Tan, -0.0), atan(I::point(0.0)), at(UnOp::Cos, 0.0), at(UnOp::Exp, 0.0), at(UnOp::Ln, 1.0)];
        assert_eq!(exact.map(|i| (i.lo, i.hi)), [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (0.0, 0.0)]);
    }
}
