//! What a constant is in real numbers, beside its float: an exact fraction
//! where its arithmetic has one (a decimal literal, a unit written as one,
//! their sums, products and quotients), else the interval its rounding
//! cannot leave. A value claim on constants reads it (claims.rs), so that a
//! comparison within rounding of its limit is not decided by the float, and
//! one the fractions decide is exact (12 in to the foot), not by luck.

use crate::interval::I;

/// A fraction n / d, d > 0, in lowest terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ratio {
    pub n: i128,
    pub d: i128,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    (a, b) = (a.abs(), b.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

impl Ratio {
    pub fn new(n: i128, d: i128) -> Option<Ratio> {
        if d == 0 {
            return None;
        }
        let g = gcd(n, d).max(1);
        let s = if d < 0 { -1 } else { 1 };
        Some(Ratio { n: s * (n / g), d: s * (d / g) })
    }

    /// A decimal as written ("25.4", "-1.5e3", "7"): its exact value, where it fits.
    pub fn of_decimal(s: &str) -> Option<Ratio> {
        let (m, e) = match s.find(['e', 'E']) {
            Some(i) => (&s[..i], s[i + 1..].parse::<i32>().ok()?),
            None => (s, 0),
        };
        let (neg, m) = match m.strip_prefix('-') {
            Some(r) => (true, r),
            None => (false, m.strip_prefix('+').unwrap_or(m)),
        };
        let (int, frac) = m.split_once('.').unwrap_or((m, ""));
        if int.is_empty() && frac.is_empty() || !int.bytes().chain(frac.bytes()).all(|c| c.is_ascii_digit()) {
            return None;
        }
        let digits: String = format!("{int}{frac}");
        let n: i128 = digits.trim_start_matches('0').parse().unwrap_or(0);
        let exp = e - frac.len() as i32;
        let pow = |k: i32| 10i128.checked_pow(k.unsigned_abs());
        let r = if exp >= 0 { Ratio::new(n.checked_mul(pow(exp)?)?, 1)? } else { Ratio::new(n, pow(exp)?)? };
        Some(if neg { r.negate() } else { r })
    }

    /// A float's exact value, where it fits.
    pub fn of_f64(v: f64) -> Option<Ratio> {
        if !v.is_finite() {
            return None;
        }
        if v == 0.0 {
            return Some(Ratio { n: 0, d: 1 });
        }
        let bits = v.to_bits();
        let (sign, exp, frac) = ((bits >> 63) as i128, ((bits >> 52) & 0x7ff) as i32, (bits & ((1 << 52) - 1)) as i128);
        let (m, e) = if exp == 0 { (frac, -1074) } else { (frac | 1 << 52, exp - 1075) };
        let r = if e >= 0 {
            Ratio::new(m.checked_mul(1i128.checked_shl(e as u32).filter(|&p| p > 0)?)?, 1)?
        } else {
            Ratio::new(m, 1i128.checked_shl((-e) as u32).filter(|&p| p > 0 && -e < 127)?)?
        };
        Some(if sign == 1 { r.negate() } else { r })
    }

    pub fn negate(self) -> Ratio {
        Ratio { n: -self.n, d: self.d }
    }

    pub fn plus(self, o: Ratio) -> Option<Ratio> {
        Ratio::new(self.n.checked_mul(o.d)?.checked_add(o.n.checked_mul(self.d)?)?, self.d.checked_mul(o.d)?)
    }

    pub fn minus(self, o: Ratio) -> Option<Ratio> {
        self.plus(o.negate())
    }

    pub fn times(self, o: Ratio) -> Option<Ratio> {
        Ratio::new(self.n.checked_mul(o.n)?, self.d.checked_mul(o.d)?)
    }

    pub fn over(self, o: Ratio) -> Option<Ratio> {
        Ratio::new(self.n.checked_mul(o.d)?, self.d.checked_mul(o.n)?)
    }

    /// Its sign against `o`: -1, 0 or 1, where the products fit.
    pub fn cmp_to(self, o: Ratio) -> Option<i32> {
        let (a, b) = (self.n.checked_mul(o.d)?, o.n.checked_mul(self.d)?);
        Some(a.cmp(&b) as i32)
    }

    /// The interval of floats that holds it: its float where that is it exactly, else a few steps either
    /// side of the nearest (n and d each rounded, then their quotient).
    pub fn band(self) -> I {
        let v = self.n as f64 / self.d as f64;
        if Ratio::of_f64(v) == Some(self) {
            return I::point(v);
        }
        I::new(v.next_down().next_down().next_down(), v.next_up().next_up().next_up())
    }
}

/// A constant in real numbers: its exact fraction where known, and the interval that holds it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Exact {
    pub ratio: Option<Ratio>,
    pub band: I,
}

impl Exact {
    /// Of a fraction.
    pub fn of(r: Ratio) -> Exact {
        Exact { ratio: Some(r), band: r.band() }
    }

    /// Two accounts of one constant (two literals with one float): what both allow.
    pub fn merge(self, o: Exact) -> Exact {
        Exact { ratio: self.ratio.filter(|r| Some(*r) == o.ratio), band: I::new(self.band.lo.min(o.band.lo), self.band.hi.max(o.band.hi)) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decimals_and_floats_are_their_exact_values() {
        assert_eq!(Ratio::of_decimal("25.4"), Ratio::new(127, 5));
        assert_eq!(Ratio::of_decimal("-1.5e3"), Ratio::new(-1500, 1));
        assert_eq!(Ratio::of_decimal("0.10"), Ratio::new(1, 10));
        assert_eq!(Ratio::of_f64(0.5), Ratio::new(1, 2));
        assert!(Ratio::of_f64(0.1).is_some_and(|r| r != Ratio::new(1, 10).unwrap()));
        let foot = Ratio::of_decimal("304.8").unwrap();
        let inch = Ratio::of_decimal("25.4").unwrap();
        assert_eq!(foot.over(inch.times(Ratio::new(12, 1).unwrap()).unwrap()), Ratio::new(1, 1));
        assert!(Ratio::new(1, 10).unwrap().band().contains(0.1) && Ratio::new(1, 3).unwrap().band().lo < 1.0 / 3.0);
    }
}
