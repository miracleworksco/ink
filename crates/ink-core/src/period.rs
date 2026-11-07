//! A tape's period in one of its variables, shown from how it is built:
//! never sampled.
//!
//! `period(tapes, v, ctx)` gives a base a > 0, an exact rational, such
//! that each tape, as the kernel evaluates it in ctx, is unchanged by v ->
//! v + 2 pi m / a, for every integer m and every value of the ranged
//! params and the coordinates. As the kernel evaluates it: an instruction
//! the evaluator takes in floats in ctx (eval::fixed: constants and params
//! no ranged one among them, a select's branches included) is the float it
//! computes (eval::fixed_value); any other is its operation on the real
//! numbers it reads (`(real a)` is a's value, a real number). v is t or a
//! param ctx ranges (else what reads it would be floats).
//! The rules, each an identity in real numbers:
//! - v is linear in v: c v + w, c an exact rational times a whole power
//!   of pi (deg->rad and rad->deg scale by pi / 180 and 180 / pi, pi
//!   exact: a turn of k degrees a second repeats every 360 / k), w anything not
//!   reading v. Sums and differences of such forms, a negation, and a
//!   product with or quotient by a value whose real number is known
//!   exactly (a float, or exact arithmetic of floats) stay linear, the
//!   coefficients computed exactly; where an exact value would not fit
//!   (a 128-bit fraction) nothing is shown;
//! - sin and cos of c v + w are unchanged by a shift s with c s a multiple
//!   of 2 pi; tan, with c s a multiple of pi;
//! - any operation of values unchanged and values not reading v is
//!   unchanged.
//!
//! Then a is the greatest rational each coefficient inside sin and cos (2
//! c for tan) is a whole multiple of, and the shift 2 pi / a is enclosed,
//! every conversion to floats rounded outward. Any other use of v (a
//! product of two forms, a form in a select's condition or at the output;
//! a tape reading a form, which may read v inside) shows no period: None.
//! A tape that does not read v is unchanged by every shift
//! (`Shown::Every`).
//!
//! `quotient(tapes, ctx)` asks one period of each range longer than one:
//! a search over the rest asks the same values again.

use crate::eval::{Ctx, fixed, fixed_value, span_mid};
use crate::graph::{BinOp, UnOp, Var};
use crate::interval::{self, I};
use crate::tape::{Op, Tape};

/// An exact rational n / d: d > 0, in lowest terms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Q {
    n: i128,
    d: i128,
}

impl Q {
    /// n / d in lowest terms; None where d is 0 or a sign overflows.
    fn new(n: i128, d: i128) -> Option<Q> {
        if d == 0 {
            return None;
        }
        let g = gcd(n.unsigned_abs(), d.unsigned_abs());
        let (n, d) = (n / g as i128, d / g as i128);
        Some(if d < 0 { Q { n: n.checked_neg()?, d: d.checked_neg()? } } else { Q { n, d } })
    }

    /// The real number the float f is; None where it is not finite or
    /// does not fit.
    fn of(f: f64) -> Option<Q> {
        if !f.is_finite() {
            return None;
        }
        if f == 0.0 {
            return Some(Q { n: 0, d: 1 });
        }
        let bits = f.to_bits();
        let e = ((bits >> 52) & 0x7ff) as i32;
        let frac = (bits & ((1 << 52) - 1)) as i128;
        let (mut m, mut e) = if e == 0 { (frac, -1074) } else { (frac | 1 << 52, e - 1075) };
        while m & 1 == 0 {
            m >>= 1;
            e += 1;
        }
        let m = if f < 0.0 { -m } else { m };
        match e {
            0.. if e < 64 => Some(Q { n: m.checked_mul(1 << e)?, d: 1 }),
            ..0 if e > -126 => Some(Q { n: m, d: 1 << -e }),
            _ => None,
        }
    }

    fn add(self, o: Q) -> Option<Q> {
        Q::new(self.n.checked_mul(o.d)?.checked_add(o.n.checked_mul(self.d)?)?, self.d.checked_mul(o.d)?)
    }

    fn neg(self) -> Option<Q> {
        Some(Q { n: self.n.checked_neg()?, d: self.d })
    }

    fn mul(self, o: Q) -> Option<Q> {
        let (g1, g2) = (gcd(self.n.unsigned_abs(), o.d.unsigned_abs()) as i128, gcd(o.n.unsigned_abs(), self.d.unsigned_abs()) as i128);
        Q::new((self.n / g1).checked_mul(o.n / g2)?, (self.d / g2).checked_mul(o.d / g1)?)
    }

    fn div(self, o: Q) -> Option<Q> {
        self.mul(Q::new(o.d, o.n)?)
    }

    /// An enclosure of it.
    pub fn value(&self) -> I {
        interval::bin(BinOp::Div, whole(self.n), whole(self.d)).0
    }
}

/// An exact coefficient: a rational times a whole power of pi (deg->rad
/// and rad->deg scale by pi / 180 and 180 / pi with pi exact: UnOp::Rad,
/// UnOp::Deg).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coef {
    q: Q,
    pi: i8,
}

impl Coef {
    /// An enclosure of 2 pi / it, the shift it is the base of: 2 pi^(1 -
    /// k) / q, each factor enclosed.
    pub fn shift(&self) -> I {
        let pi = I::new(std::f64::consts::PI, std::f64::consts::PI.next_up());
        let e = 1 - self.pi as i32;
        let power = (0..e.unsigned_abs()).fold(I::point(1.0), |a, _| mul_o(a, pi));
        let power = if e < 0 { interval::bin(BinOp::Div, I::point(1.0), power).0 } else { power };
        interval::bin(BinOp::Div, mul_o(I::point(2.0), power), self.q.value()).0
    }

    /// n times it.
    pub fn times(self, n: i128) -> Option<Coef> {
        self.by(Q::new(n, 1)?, false)
    }

    /// The sum, of one power of pi: None where the powers differ or it
    /// does not fit; Some(None) where it is 0.
    pub fn plus(self, o: Coef) -> Option<Option<Coef>> {
        if self.pi != o.pi {
            return None;
        }
        let q = self.q.add(o.q)?;
        Some((q.n != 0).then_some(Coef { q, pi: self.pi }))
    }

    /// Times or over the rational x.
    fn by(self, x: Q, over: bool) -> Option<Coef> {
        Some(Coef { q: if over { self.q.div(x)? } else { self.q.mul(x)? }, ..self })
    }
}

fn mul_o(a: I, b: I) -> I {
    interval::bin(BinOp::Mul, a, b).0
}

/// An enclosure of the whole number n: its float, or the floats either
/// side where that is rounded.
fn whole(n: i128) -> I {
    let f = n as f64;
    // f's own value, exactly, where it fits (|f| < 2^127): n's float is n.
    let exact = f.abs() < 1.7e38 && f as i128 == n;
    if exact { I::point(f) } else { I::new(f.next_down(), f.next_up()) }
}

/// What `period` shows.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Shown {
    /// The tapes do not read v: every shift keeps them.
    Every,
    /// Every shift by a multiple of 2 pi / base keeps them.
    Base(Coef),
}

/// What an instruction is, as a function of v.
#[derive(Clone, Copy)]
enum Is {
    /// It does not read v: its real value, where known exactly.
    Free(Option<Q>),
    /// c v plus something not reading v.
    Lin(Coef),
    /// Unchanged by the shifts the coefficients `seen` allow.
    Kept,
}

/// The period of `tapes` in `v` as evaluated in `ctx` (the module's
/// contract). None: no period shown.
pub fn period(tapes: &[&Tape], v: Var, ctx: &Ctx) -> Option<Shown> {
    period_kept(tapes, v, ctx, &[])
}

/// `period`, where each var of `kept` stands for a value that every shift
/// by a whole multiple of 2 pi / its coefficient keeps (None: every shift),
/// and that is otherwise unknown: a shape turned with v whose own turn
/// symmetry absorbs the turn. Each coefficient joins the ones sin and cos
/// read, so the base is a whole divisor of it too.
pub fn period_kept(tapes: &[&Tape], v: Var, ctx: &Ctx, kept: &[(Var, Option<Coef>)]) -> Option<Shown> {
    if let Var::Param(k) = v
        && ctx.p_spans.get(k).is_none()
    {
        return None;
    }
    let mut seen: Vec<Coef> = kept.iter().filter_map(|k| k.1).collect();
    for tape in tapes {
        if matches!(read(tape, v, ctx, kept, &mut seen)?, Is::Lin(_)) {
            return None;
        }
    }
    if seen.is_empty() {
        return Some(Shown::Every);
    }
    // One power of pi, and the greatest rational each is a whole multiple
    // of: the gcd of their numerators over a common denominator, over it.
    let pi = seen[0].pi;
    if seen.iter().any(|c| c.pi != pi) {
        return None;
    }
    let d = seen.iter().try_fold(1i128, |l, c| l.checked_mul(c.q.d / gcd(l.unsigned_abs(), c.q.d.unsigned_abs()) as i128))?;
    let n = seen.iter().try_fold(0u128, |g, c| Some(gcd(g, c.q.n.checked_mul(d / c.q.d)?.unsigned_abs())))?;
    Some(Shown::Base(Coef { q: Q::new(i128::try_from(n).ok()?, d)?, pi }))
}

/// A tape's value as a function of v: c v plus what does not read v
/// (Some(Some(c)), c never 0), or not reading v (Some(None)). None: any
/// other.
pub fn form(tape: &Tape, v: Var, ctx: &Ctx) -> Option<Option<Coef>> {
    match read(tape, v, ctx, &[], &mut Vec::new())? {
        Is::Lin(c) => Some(Some(c)),
        Is::Free(_) => Some(None),
        Is::Kept => None,
    }
}

/// What the tape's output is as a function of v (the module's rules),
/// each coefficient sin and cos read pushed to `seen`. None: a use of v
/// the rules do not read.
fn read(tape: &Tape, v: Var, ctx: &Ctx, kept: &[(Var, Option<Coef>)], seen: &mut Vec<Coef>) -> Option<Is> {
    if !tape.forms.is_empty() {
        return None;
    }
    let mut is: Vec<Is> = Vec::with_capacity(tape.ssa.len());
    // The floats the evaluator takes its fixed instructions as.
    let mut floats: Vec<f64> = vec![f64::NAN; tape.ssa.len()];
    for (i, ins) in tape.ssa.iter().enumerate() {
        let (args, n) = Tape::reads(&tape.ssa, i);
        let arg = |j: usize| is[args[j] as usize];
        let exact = |j: usize| match arg(j) {
            Is::Free(q) => q,
            _ => None,
        };
        let reads = (0..n).any(|j| !matches!(arg(j), Is::Free(..)));
        let out = match ins.op {
            Op::Var(w) if same(w, v) => Is::Lin(Coef { q: Q { n: 1, d: 1 }, pi: 0 }),
            // Unchanged by the shifts its coefficient allows: no known value.
            Op::Var(w) if kept.iter().any(|k| same(k.0, w)) => Is::Free(None),
            // The float the evaluator computes, and so its value.
            _ if fixed(tape, ctx, i) => {
                floats[i] = fixed_value(ins, |a| floats[a as usize], ctx);
                Is::Free(Q::of(floats[i]))
            }
            _ if !reads => Is::Free(match ins.op {
                Op::Un(UnOp::Neg) => exact(0).and_then(Q::neg),
                Op::Un(UnOp::Square) => exact(0).and_then(|a| a.mul(a)),
                Op::Bin(BinOp::Add) => exact(0).zip(exact(1)).and_then(|(a, b)| a.add(b)),
                Op::Bin(BinOp::Sub) => exact(0).zip(exact(1).and_then(Q::neg)).and_then(|(a, b)| a.add(b)),
                Op::Bin(BinOp::Mul) => exact(0).zip(exact(1)).and_then(|(a, b)| a.mul(b)),
                // x / 0 is 0, no quotient.
                Op::Bin(BinOp::Div) => exact(0).zip(exact(1)).and_then(|(a, b)| a.div(b)),
                Op::Tag(_) | Op::Real => exact(0),
                _ => None,
            }),
            _ if !(0..n).any(|j| matches!(arg(j), Is::Lin(_))) => Is::Kept,
            Op::Un(UnOp::Neg) => Is::Lin(lin(arg(0))?.by(Q { n: -1, d: 1 }, false)?),
            // Degrees in radians, and back: pi exact.
            Op::Un(UnOp::Rad) => {
                let c = lin(arg(0))?;
                Is::Lin(Coef { q: c.q.div(Q { n: 180, d: 1 })?, pi: c.pi.checked_add(1)? })
            }
            Op::Un(UnOp::Deg) => {
                let c = lin(arg(0))?;
                Is::Lin(Coef { q: c.q.mul(Q { n: 180, d: 1 })?, pi: c.pi.checked_sub(1)? })
            }
            Op::Tag(_) | Op::Real => arg(0),
            Op::Bin(BinOp::Add | BinOp::Sub) => {
                let b = match arg(1) {
                    Is::Lin(c) if matches!(ins.op, Op::Bin(BinOp::Sub)) => Is::Lin(c.by(Q { n: -1, d: 1 }, false)?),
                    b => b,
                };
                match (arg(0), b) {
                    (Is::Lin(c), Is::Free(..)) | (Is::Free(..), Is::Lin(c)) => Is::Lin(c),
                    (Is::Lin(c), Is::Lin(e)) => match c.plus(e)? {
                        None => Is::Free(None),
                        Some(c) => Is::Lin(c),
                    },
                    _ => return None,
                }
            }
            Op::Bin(BinOp::Mul) => match (arg(0), arg(1)) {
                (Is::Lin(c), Is::Free(Some(x))) | (Is::Free(Some(x)), Is::Lin(c)) if x.n != 0 => Is::Lin(c.by(x, false)?),
                _ => return None,
            },
            Op::Bin(BinOp::Div) => match (arg(0), arg(1)) {
                (Is::Lin(c), Is::Free(Some(x))) if x.n != 0 => Is::Lin(c.by(x, true)?),
                _ => return None,
            },
            Op::Un(UnOp::Sin | UnOp::Cos) => {
                seen.push(lin(arg(0))?);
                Is::Kept
            }
            Op::Un(UnOp::Tan) => {
                seen.push(lin(arg(0))?.by(Q { n: 2, d: 1 }, false)?);
                Is::Kept
            }
            _ => return None,
        };
        is.push(out);
    }
    Some(is.last().copied().unwrap_or(Is::Free(None)))
}

/// A form's coefficient (None: not a form).
fn lin(a: Is) -> Option<Coef> {
    match a {
        Is::Lin(c) => Some(c),
        _ => None,
    }
}

fn same(a: Var, b: Var) -> bool {
    matches!((a, b), (Var::X, Var::X) | (Var::Y, Var::Y) | (Var::Z, Var::Z) | (Var::T, Var::T)) || matches!((a, b), (Var::Param(p), Var::Param(q)) if p == q)
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.max(1)
}

/// ctx with each independent ranged param, and t's span, longer than a
/// period the tapes show in it (`period`; none where a derived param they
/// read rests on it) narrowed to one period from its low end (rounded up),
/// its band with it: every value the tapes take over ctx they take over
/// the result, at values inside ctx's spans (a moment past the period is
/// one in it moved by whole periods). None: nothing narrowed.
pub fn quotient(tapes: &[&Tape], ctx: &Ctx) -> Option<Ctx> {
    let read: Vec<u32> = tapes.iter().flat_map(|t| t.ssa.iter()).filter_map(|i| if let Op::Var(Var::Param(k)) = i.op { Some(k) } else { None }).collect();
    let one = |v: Var, (lo, hi): (f64, f64)| -> Option<(f64, f64)> {
        let Some(Shown::Base(b)) = period(tapes, v, ctx) else { return None };
        let p = b.shift().hi;
        (p > 0.0 && hi - lo > p).then(|| (lo, (lo + p).next_up()))
    };
    let mut out = ctx.clone();
    let mut narrowed = false;
    for s in ctx.p_spans.iter().filter(|s| !s.derived) {
        let rests = read.iter().any(|&j| ctx.p_spans.is_derived(j) && ctx.p_spans.sources(j).contains(&s.k));
        if let Some(span) = (!rests).then(|| one(Var::Param(s.k), s.span)).flatten() {
            out.p_spans.narrow(s.k, span);
            out.params[s.k as usize] = span_mid(span);
            narrowed = true;
        }
    }
    if let Some(span) = ctx.t_span.and_then(|t| one(Var::T, t)) {
        (out.t_span, out.t) = (Some(span), span_mid(span));
        narrowed = true;
    }
    narrowed.then_some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Graph;

    /// The period in t of the expression `f` builds (params at `values`,
    /// those `ranged` over their bands), enclosed: None where none is
    /// shown, Some(None) where t is not read.
    fn shown_over(values: &[f64], ranged: &[(u32, (f64, f64))], f: impl Fn(&mut Graph, u32) -> u32) -> Option<Option<I>> {
        let mut g = Graph::new();
        let t = g.var(Var::T);
        let root = f(&mut g, t);
        let tape = Tape::new(&g, root);
        let ctx = Ctx { p_spans: crate::eval::ParamSpans::new(ranged.iter().copied()), ..Ctx::new(0.0, values.to_vec()) };
        period(&[&tape], Var::T, &ctx).map(|s| match s {
            Shown::Every => None,
            Shown::Base(b) => Some(b.shift()),
        })
    }

    fn shown(values: &[f64], f: impl Fn(&mut Graph, u32) -> u32) -> Option<Option<I>> {
        shown_over(values, &[], f)
    }

    fn trig(g: &mut Graph, op: UnOp, t: u32, k: f64, over: bool) -> u32 {
        let c = g.constant(k);
        let a = if over { g.bin(BinOp::Div, t, c) } else { g.mul(c, t) };
        g.un(op, a)
    }

    /// Each rule: whole multiples of one coefficient (written as a
    /// product with a float, or a quotient by a whole number) repeat
    /// together; two floats repeat together only over the period of their
    /// exact gcd; a form outside sin and cos, and a product of forms, do
    /// not repeat.
    #[test]
    fn periods_are_shown_from_how_a_tape_is_built() {
        let tau = std::f64::consts::TAU;
        let near = |got: Option<Option<I>>, want: f64| matches!(got, Some(Some(p)) if p.lo <= p.hi && (p.mid() - want).abs() < 1e-12 * want);
        let sum = |a: (UnOp, f64, bool), b: (UnOp, f64, bool)| {
            shown(&[], move |g, t| {
                let x = trig(g, a.0, t, a.1, a.2);
                let y = trig(g, b.0, t, b.1, b.2);
                g.add(x, y)
            })
        };
        assert!(near(sum((UnOp::Sin, 10.0, false), (UnOp::Cos, 20.0, false)), tau / 10.0));
        assert!(near(sum((UnOp::Sin, 3.0, true), (UnOp::Cos, 1.5, true)), 3.0 * tau));
        assert!(near(shown(&[], |g, t| trig(g, UnOp::Sin, t, 0.3, false)), tau / 0.3));
        // 0.3 and 0.7 are dyadic fractions with denominators 2^54 and 2^53.
        assert!(matches!(sum((UnOp::Sin, 0.3, false), (UnOp::Sin, 0.7, false)), Some(Some(p)) if p.lo > 1e15));
        assert!(
            shown(&[], |g, t| {
                let a = trig(g, UnOp::Sin, t, 1.0, false);
                g.add(a, t)
            })
            .is_none()
        );
        assert!(
            shown(&[], |g, t| {
                let s = g.mul(t, t);
                g.un(UnOp::Sin, s)
            })
            .is_none()
        );
        assert_eq!(shown(&[], |g, _| g.constant(1.0)), Some(None));
    }

    /// A coefficient that is arithmetic in real numbers is that real
    /// number, not its float: 1 + p read as a real (a warp's x moved to
    /// the param p = 0.1) is not 1.1's float, so sin((1 + p) t) and sin(1.1
    /// t) do not repeat together over 2 pi / 1.1 (at t = 1e16 they differ
    /// by 0.094), only over their exact gcd's period. Taken in floats
    /// (1 + p fixed), they are one coefficient.
    #[test]
    fn a_real_coefficient_is_its_real_number() {
        let pair = |real: bool| {
            shown(&[0.1], move |g, t| {
                let p = g.var(Var::Param(0));
                let p = if real { g.real(p) } else { p };
                let one = g.constant(1.0);
                let c = g.add(p, one);
                let a = g.mul(c, t);
                let a = g.un(UnOp::Sin, a);
                let b = trig(g, UnOp::Sin, t, 1.1, false);
                g.sub(a, b)
            })
        };
        assert!(matches!(pair(true), Some(Some(p)) if p.lo > 1e16), "{:?}", pair(true));
        assert!(matches!(pair(false), Some(Some(p)) if (p.mid() - std::f64::consts::TAU / 1.1).abs() < 1e-12), "{:?}", pair(false));
    }

    /// A select of a ranged param is real arithmetic whichever branch
    /// its condition picks, as the evaluator takes it (eval::fixed): 1 +
    /// (select c 0.1 q), c = -1, q ranged, is 1 + 0.1 in real numbers, not
    /// 1.1's float; with q a point it is the float the evaluator computes.
    #[test]
    fn a_select_of_a_ranged_param_is_real() {
        let pair = |ranged: &[(u32, (f64, f64))]| {
            shown_over(&[-1.0, 0.0], ranged, |g, t| {
                let (c, q) = (g.var(Var::Param(0)), g.var(Var::Param(1)));
                let tenth = g.constant(0.1);
                let s = g.select(c, tenth, q);
                let one = g.constant(1.0);
                let k = g.add(one, s);
                let a = g.mul(k, t);
                let a = g.un(UnOp::Sin, a);
                let b = trig(g, UnOp::Sin, t, 1.1, false);
                g.sub(a, b)
            })
        };
        let tau_11 = std::f64::consts::TAU / 1.1;
        assert!(!matches!(pair(&[(1, (-1.0, 1.0))]), Some(Some(p)) if p.lo < 1e16), "{:?}", pair(&[(1, (-1.0, 1.0))]));
        assert!(matches!(pair(&[]), Some(Some(p)) if (p.mid() - tau_11).abs() < 1e-12), "{:?}", pair(&[]));
    }

    /// A turn written in degrees (deg->rad, pi exact) has a period with no
    /// pi in it: sin(deg->rad(360 t)) repeats every 1 exactly, and with
    /// rad->deg back, every 2 pi.
    #[test]
    fn a_turn_in_degrees_has_an_exact_period() {
        let got = shown(&[], |g, t| {
            let c = g.constant(360.0);
            let a = g.mul(c, t);
            let a = g.un(UnOp::Rad, a);
            g.un(UnOp::Sin, a)
        });
        assert_eq!(got, Some(Some(I::point(1.0))));
        let back = shown(&[], |g, t| {
            let a = g.un(UnOp::Deg, t);
            let a = g.un(UnOp::Rad, a);
            g.un(UnOp::Cos, a)
        });
        assert!(matches!(back, Some(Some(p)) if p.lo <= std::f64::consts::TAU && std::f64::consts::TAU.next_up() <= p.hi), "{back:?}");
    }

    /// The period's enclosure holds it where its fraction's terms are not
    /// floats: sin(t a b / c / d) with a..d whole numbers about 2^30, the
    /// base 670773797348616147 / 1085022332563062331, its period
    /// 10.163480453275372698..., strictly between two floats.
    #[test]
    fn a_period_is_enclosed_outward() {
        let rates = [1128473121.0, 1188816614.0, 1250991818.0, 1734659359.0];
        let got = shown(&rates, |g, t| {
            let mut f = t;
            for k in 0..4 {
                let p = g.var(Var::Param(k));
                f = g.bin(if k < 2 { BinOp::Mul } else { BinOp::Div }, f, p);
            }
            g.un(UnOp::Sin, f)
        });
        // The float below it (10.16348045327537263915...) and the one above.
        let below = 10.163480453275373_f64;
        assert!(matches!(got, Some(Some(p)) if p.lo <= below && below.next_up() <= p.hi), "{got:?}");
    }
}
