//! Scalar semantics shared by constant folding and every evaluator.

use super::{BinOp, UnOp};

// The core's value rules (docs/CORE.md): the IEEE-754 result, with one zero
// (-0 reads as 0), x / 0 = 0, and anything else undefined (NaN) = 0. So
// sqrt and ln of a negative are 0, pow of a negative base with an exponent
// that is not a whole number is 0; ln 0 is -inf. conformance/ tests every
// evaluator against an independent reference.

/// Where each operation is undefined (docs/IR.md, Fields, Value rules: its
/// value is 0, and it is reported): the one table of operator domains, read
/// over operands' ranges [lo, hi] (a point for a value). The case an op
/// may meet there, if any. Every check reads it: the evaluators' report and
/// constant folding (at points: `undefined_un`, `undefined_bin`), and
/// whether a promise's constants are real numbers (eval::promise_valid).
/// sqrt of a negative is no such case: it is sqrt(max(a, 0)), continuous.
pub fn undefined_un_over(op: UnOp, [lo, hi]: [f64; 2]) -> Option<&'static str> {
    match op {
        UnOp::Ln if lo < 0.0 => Some("ln of a negative"),
        UnOp::Sin | UnOp::Cos | UnOp::Tan if lo.is_infinite() || hi.is_infinite() => Some("sin, cos or tan of infinity"),
        UnOp::Neg | UnOp::Abs | UnOp::Square | UnOp::Exp | UnOp::Ln | UnOp::Sqrt | UnOp::Sin | UnOp::Cos | UnOp::Tan | UnOp::Floor | UnOp::Rad | UnOp::Deg => None,
    }
}

/// See `undefined_un_over`. A NaN bound may be anything.
pub fn undefined_bin_over(op: BinOp, a: [f64; 2], b: [f64; 2]) -> Option<&'static str> {
    let nan = a.iter().chain(&b).any(|x| x.is_nan());
    let (zero, infinite) = (|r: [f64; 2]| r[0] <= 0.0 && r[1] >= 0.0, |r: [f64; 2]| r[0] == f64::NEG_INFINITY || r[1] == f64::INFINITY);
    match op {
        BinOp::Div if zero(b) => Some("division by zero"),
        BinOp::Div if nan || (infinite(a) && infinite(b)) => Some("infinity divided by infinity"),
        BinOp::Add if nan || (a[1] == f64::INFINITY && b[0] == f64::NEG_INFINITY) || (a[0] == f64::NEG_INFINITY && b[1] == f64::INFINITY) => Some("infinity minus infinity"),
        BinOp::Sub if nan || (a[1] == f64::INFINITY && b[1] == f64::INFINITY) || (a[0] == f64::NEG_INFINITY && b[0] == f64::NEG_INFINITY) => Some("infinity minus infinity"),
        BinOp::Mul if nan || (zero(a) && infinite(b)) || (zero(b) && infinite(a)) => Some("zero times infinity"),
        BinOp::Pow if a[0] < 0.0 && !(b[0] == b[1] && b[0].fract() == 0.0) => Some("a negative base to a power that is not a whole number"),
        BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Div | BinOp::Pow | BinOp::Min | BinOp::Max | BinOp::Atan2 | BinOp::Emax => None,
    }
}

/// The undefined case an op meets at these arguments (`undefined_un_over`
/// at a point).
pub fn undefined_un(op: UnOp, a: f64) -> Option<&'static str> {
    undefined_un_over(op, [a, a])
}

/// See `undefined_un`.
pub fn undefined_bin(op: BinOp, a: f64, b: f64) -> Option<&'static str> {
    undefined_bin_over(op, [a, a], [b, b])
}

/// The core's value for an IEEE result: NaN (undefined) is 0. Library
/// functions (exp, ln, sin, cos, tan, atan, atan2, pow, hypot) are the
/// `libm` crate's, plain f64 code with no CPU dispatch, not the system's:
/// a model means the same, bit for bit, on every machine, so its answers
/// do too (and the result cache needs no maths library in its keys).
#[inline]
pub fn defined(v: f64) -> f64 {
    if v.is_nan() { 0.0 } else { v }
}

/// `a / b` with `x / 0 = 0`.
#[inline]
pub fn div_f64(a: f64, b: f64) -> f64 {
    if b == 0.0 { 0.0 } else { defined(a / b) }
}

#[inline]
pub fn pow_f64(a: f64, b: f64) -> f64 {
    // + 0.0: one zero, so pow(-0, -1) is +inf like pow(0, -1).
    if a < 0.0 && b.fract() != 0.0 { 0.0 } else { libm::pow(a + 0.0, b) }
}

#[inline]
pub fn ln_f64(a: f64) -> f64 {
    if a > 0.0 {
        libm::log(a)
    } else if a == 0.0 {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

#[inline]
pub fn atan2_f64(y: f64, x: f64) -> f64 {
    // + 0.0: one zero, so the negative x axis is +pi from either side of 0.
    libm::atan2(y + 0.0, x + 0.0)
}

#[inline]
pub fn un_f64(op: UnOp, a: f64) -> f64 {
    match op {
        UnOp::Neg => -a,
        UnOp::Abs => a.abs(),
        UnOp::Sqrt => a.max(0.0).sqrt(),
        UnOp::Square => a * a,
        UnOp::Exp => libm::exp(a),
        UnOp::Ln => ln_f64(a),
        UnOp::Sin => defined(libm::sin(a)),
        UnOp::Cos => defined(libm::cos(a)),
        UnOp::Tan => defined(libm::tan(a)),
        UnOp::Floor => a.floor(),
        UnOp::Rad => a * (std::f64::consts::PI / 180.0),
        UnOp::Deg => a * (180.0 / std::f64::consts::PI),
    }
}

#[inline]
pub fn bin_f64(op: BinOp, a: f64, b: f64) -> f64 {
    match op {
        BinOp::Add => defined(a + b),
        BinOp::Sub => defined(a - b),
        BinOp::Mul => defined(a * b),
        BinOp::Div => div_f64(a, b),
        BinOp::Min => a.min(b),
        BinOp::Max => a.max(b),
        BinOp::Atan2 => atan2_f64(a, b),
        BinOp::Pow => pow_f64(a, b),
        BinOp::Emax => emax_f64(a, b),
    }
}

#[inline]
pub fn emax_f64(a: f64, b: f64) -> f64 {
    if a > 0.0 && b > 0.0 { libm::hypot(a, b) } else { a.max(b) }
}

/// Below this, a product's or quotient's rounding error may itself round
/// away, so its float is not shown exact.
const EXACT_TINY: f64 = 1e-270;

/// `op` of `a` where its float is its real value, shown so (an exact
/// square root, a square without rounding, sin 0), else None: what a
/// rewrite may fold an operation that read a spatial value into
/// (`Graph::un_real`).
pub fn exact_un(op: UnOp, a: f64) -> Option<f64> {
    if !a.is_finite() || undefined_un(op, a).is_some() {
        return None;
    }
    let v = un_f64(op, a);
    let exact = match op {
        UnOp::Neg | UnOp::Abs | UnOp::Floor => true,
        UnOp::Square => a == 0.0 || (v.is_finite() && v >= EXACT_TINY && a.mul_add(a, -v) == 0.0),
        UnOp::Sqrt => a <= 0.0 || (a >= EXACT_TINY && v.mul_add(v, -a) == 0.0),
        UnOp::Exp | UnOp::Cos => a == 0.0,
        UnOp::Sin | UnOp::Tan => a == 0.0,
        UnOp::Ln => a == 1.0,
        UnOp::Rad | UnOp::Deg => a == 0.0,
    };
    exact.then_some(v)
}

/// `exact_un` for two operands: a sum or difference with no rounding
/// error, a product or quotient that is exact, min and max.
pub fn exact_bin(op: BinOp, a: f64, b: f64) -> Option<f64> {
    if !a.is_finite() || !b.is_finite() || undefined_bin(op, a, b).is_some() {
        return None;
    }
    let v = bin_f64(op, a, b);
    if !v.is_finite() {
        return None;
    }
    // The sum's rounding error (Knuth's two-sum), which is a float.
    let two_sum = |a: f64, b: f64| {
        let s = a + b;
        let bb = s - a;
        (a - (s - bb)) + (b - bb)
    };
    let exact = match op {
        BinOp::Add => two_sum(a, b) == 0.0,
        BinOp::Sub => two_sum(a, -b) == 0.0,
        BinOp::Mul => a == 0.0 || b == 0.0 || (v.abs() >= EXACT_TINY && a.mul_add(b, -v) == 0.0),
        BinOp::Div => a == 0.0 || (v.abs() >= EXACT_TINY && a.abs() >= EXACT_TINY && v.mul_add(b, -a) == 0.0),
        BinOp::Min | BinOp::Max => true,
        BinOp::Emax => !(a > 0.0 && b > 0.0),
        BinOp::Atan2 | BinOp::Pow => false,
    };
    exact.then_some(v)
}

// Derivatives: each formula evaluated under the value rules, on the branch
// the value took; at an exact tie between branches (min, max, abs at 0) the
// average of both sides, so the result does not depend on argument order.

/// Value and derivative d/da.
#[inline]
pub fn dun(op: UnOp, a: f64) -> (f64, f64) {
    match op {
        UnOp::Neg => (-a, -1.0),
        UnOp::Abs => (
            a.abs(),
            if a < 0.0 {
                -1.0
            } else if a > 0.0 {
                1.0
            } else {
                0.0
            },
        ),
        UnOp::Sqrt => {
            let v = a.max(0.0).sqrt();
            (v, if v > 0.0 { 0.5 / v } else { 0.0 })
        }
        UnOp::Square => (a * a, 2.0 * a),
        UnOp::Exp => {
            let e = libm::exp(a);
            (e, e)
        }
        UnOp::Ln => (ln_f64(a), if a > 0.0 { 1.0 / a } else { 0.0 }),
        UnOp::Sin => (defined(libm::sin(a)), defined(libm::cos(a))),
        UnOp::Cos => (defined(libm::cos(a)), defined(-libm::sin(a))),
        UnOp::Tan => {
            let t = defined(libm::tan(a));
            (t, 1.0 + t * t)
        }
        UnOp::Floor => (a.floor(), 0.0),
        UnOp::Rad => {
            let k = std::f64::consts::PI / 180.0;
            (a * k, k)
        }
        UnOp::Deg => {
            let k = 180.0 / std::f64::consts::PI;
            (a * k, k)
        }
    }
}

/// Value and partials d/da, d/db.
#[inline]
pub fn dbin(op: BinOp, a: f64, b: f64) -> (f64, f64, f64) {
    match op {
        BinOp::Add => (defined(a + b), 1.0, 1.0),
        BinOp::Sub => (defined(a - b), 1.0, -1.0),
        BinOp::Mul => (defined(a * b), b, a),
        BinOp::Div => {
            if b == 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                // -(a / b) / b, not -a / b^2: b^2 underflows to 0 for tiny b.
                (defined(a / b), 1.0 / b, defined(-(a / b) / b))
            }
        }
        // At exact ties the average of the two branch derivatives: surface
        // points often land exactly on ties (e.g. min(d, 0) on a box face),
        // and a sharp edge gets the direction between its two faces.
        BinOp::Min => {
            if a == b {
                (a, 0.5, 0.5)
            } else if a < b {
                (a, 1.0, 0.0)
            } else {
                (b, 0.0, 1.0)
            }
        }
        BinOp::Max => {
            if a == b {
                (a, 0.5, 0.5)
            } else if a > b {
                (a, 1.0, 0.0)
            } else {
                (b, 0.0, 1.0)
            }
        }
        BinOp::Atan2 => {
            let r2 = a * a + b * b;
            if r2 > 0.0 { (atan2_f64(a, b), defined(b / r2), defined(-a / r2)) } else { (atan2_f64(a, b), 0.0, 0.0) }
        }
        BinOp::Pow => {
            let v = pow_f64(a, b);
            // 0 where the value is the constant 0 (a negative base, an
            // exponent that is not a whole number); b - 1 can round to a
            // whole number, so test b itself.
            let fa = if b == 0.0 || (a < 0.0 && b.fract() != 0.0) { 0.0 } else { defined(b * pow_f64(a, b - 1.0)) };
            let fb = if a > 0.0 { defined(v * libm::log(a)) } else { 0.0 };
            (v, fa, fb)
        }
        BinOp::Emax => {
            if a > 0.0 && b > 0.0 {
                // Each weight at most 1 (hypot may round below a or b).
                let e = libm::hypot(a, b);
                (e, defined(a / e).min(1.0), defined(b / e).min(1.0))
            } else {
                dbin(BinOp::Max, a, b)
            }
        }
    }
}

/// Whether `v` is outside ink's range: finite, not 0, and smaller than
/// 1e-100 or larger than 1e100 in size.
pub fn wild(v: f64) -> bool {
    v != 0.0 && v.is_finite() && !(1e-100..=1e100).contains(&v.abs())
}

/// Whether floats lost `r`, computed from `args` (by an operation that
/// `scales`: a product, a quotient, a power, exp): a subnormal, 0 from
/// finite numbers that are not, or an infinity from them (an overflow; ln 0
/// is exactly -inf).
pub fn lost(r: f64, args: &[f64], scales: bool) -> bool {
    let real = args.iter().all(|&a| a != 0.0 && a.is_finite());
    (r != 0.0 && r.abs() < f64::MIN_POSITIVE) || (real && (r.is_infinite() || (scales && r == 0.0)))
}
