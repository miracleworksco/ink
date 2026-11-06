//! ink's range (docs/CORE.md Answers): a claim outside it is not checked,
//! but unknown, saying why.

use crate::view::ModelView;
use crate::{Claim, Doc};
use ink_core::eval::Ctx;
use ink_core::graph::{Var, wild};
use ink_core::interval::I;
use ink_core::tape::Tape;

/// The farthest coordinate (mm): 1,000 km.
const REACH: f64 = 1e9;
/// The finest length asked (mm), and the finest relative to the farthest
/// coordinate: a nanometre, a micrometre at 1 km, a millimetre at 1,000 km.
const FINEST: f64 = 1e-6;
const RELATIVE: f64 = 1e-9;

/// Why `claim` (its fields `roots`, to accuracy `tol`, at `res`) is
/// outside ink's range, if it is.
pub(crate) fn outside(v: &mut ModelView, claim: &Claim, roots: &[u32], tol: f64, res: Option<f64>) -> Option<String> {
    let scope = scope(v);
    if let Some(why) = numbers(&v.doc, &scope, roots) {
        return Some(why);
    }
    let (mut far, mut finest) = (0.0f64, f64::INFINITY);
    match claim {
        Claim::Point { at, .. } => {
            for &n in at {
                let c = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, n), &scope, [I::point(0.0); 3], None);
                far = far.max(c.lo.abs()).max(c.hi.abs());
            }
        }
        Claim::Value(..) => {}
        _ => {
            for r in searched(v, claim, roots) {
                let reach = if v.motion.is_some() { v.at_point(scope.clone(), |v| v.reach(r)) } else { v.reach(r) };
                far = reach.lo.iter().chain(&reach.hi).fold(far, |m, c| m.max(c.abs()));
            }
            let length = v.doc.accuracy_quantity(claim) == Some(crate::units::Quantity::Length);
            finest = asked(if length { tol } else { f64::INFINITY }, res, v.finest());
            if let Claim::Clearance(_, _, d) = *claim
                && d > 0.0
            {
                finest = finest.min(d);
            }
        }
    }
    scaled(far, finest)
}

/// What the view asks over, its motion's span the moment's (a claim's,
/// interfere's): the guard judges the span asked, its ends included.
pub(crate) fn scope(v: &ModelView) -> Ctx {
    Ctx { t_span: v.motion.or(v.ctx.t_span), ..v.ctx.clone() }
}

/// The fields the claim's search looks in (claims::decide), so their reach
/// is the one it finds anyway.
fn searched(v: &mut ModelView, claim: &Claim, roots: &[u32]) -> Vec<u32> {
    let g = &mut v.doc.g;
    match *claim {
        Claim::Empty(f) | Claim::Range(f) => vec![f],
        Claim::Clearance(a, b, d) if d > 0.0 => vec![a, b],
        Claim::Disjoint(a, b) | Claim::Overlaps(a, b) | Claim::Clearance(a, b, _) => vec![g.max(a, b)],
        Claim::Subset(a, b) => {
            let nb = g.neg(b);
            vec![g.max(a, nb)]
        }
        Claim::Touches(a, b) => vec![g.max(a, b), a, b],
        Claim::Integral { ref parts, .. } => parts.iter().map(|p| p.0).collect(),
        Claim::Mean { ref parts, .. } => parts.iter().map(|p| p.0).collect(),
        Claim::Metric { ref metric, .. } => metric.shape.into_iter().collect(),
        _ => roots.iter().copied().filter(|&r| g.depends_on(r, |x| matches!(x, Var::X | Var::Y | Var::Z))).collect(),
    }
}

/// `outside` for one of `ink interfere`'s pairs: whether a field is empty
/// over `region`, moved by `moved` (a shared frame's placement).
pub(crate) fn outside_region(doc: &Doc, ctx: &Ctx, (lo, hi): (ink_core::V3, ink_core::V3), moved: Option<&ink_kernel::rigid::Rigid>, tol: f64, res: Option<f64>) -> Option<String> {
    if let Some(why) = numbers(doc, ctx, &[]) {
        return Some(why);
    }
    let mut far = lo.iter().chain(&hi).fold(0.0f64, |m, c| m.max(c.abs()));
    // A turn takes a coordinate no farther than sqrt(3) of the farthest.
    if let Some(at) = moved {
        far = 2.0 * far + at.t.iter().fold(0.0f64, |m, c| m.max(c.abs()));
    }
    scaled(far, asked(tol, res, doc.finest()))
}

/// The finest length a search to `tol` (0: exactly, to cells of `exact`)
/// at `res` asks about.
fn asked(tol: f64, res: Option<f64>, exact: f64) -> f64 {
    let a = if tol > 0.0 { tol } else { exact };
    res.filter(|&r| r > 0.0).map_or(a, |r| a.min(r))
}

/// A command's answer as a whole (`ink interfere`'s, `distance`'s,
/// `bounds`'): refused, saying why, where the model's numbers or those
/// `fields` compute are outside the range, or a box it reports or searches
/// (`boxes`), or the finest length it asks about.
pub fn refuse(doc: &Doc, ctx: &Ctx, fields: &[u32], boxes: impl IntoIterator<Item = (ink_core::V3, ink_core::V3)>, finest: f64) -> Result<(), crate::InkError> {
    let far = boxes.into_iter().flat_map(|(lo, hi)| lo.into_iter().chain(hi)).fold(0.0f64, |m, c| m.max(c.abs()));
    match numbers(doc, ctx, fields).or_else(|| scaled(far, finest)) {
        Some(why) => Err(crate::InkError::msg(why)),
        None => Ok(()),
    }
}

/// A number the model is built or asked with outside the range: one it
/// is built from (Graph::wild: written, read from a file, made or read at
/// load in Rust; a fold floats lost), a param's value or band, the moment
/// or span, or a value `fields` compute in floats from params that floats
/// lost (`lost`). Over a band or a span values are bounded outward (sound
/// at every scale): only its ends are judged. An unsolved param's NaN is
/// no number.
fn numbers(doc: &Doc, ctx: &Ctx, fields: &[u32]) -> Option<String> {
    const RANGE: &str = "outside ink's range (0, or 1e-100 to 1e100 in size)";
    let what = |x: f64| {
        if wild(x) {
            format!("the number {x:e} is {RANGE}")
        } else {
            format!("an {} (a value floats rounded to {:e}) is {RANGE}", if x.is_infinite() { "overflow" } else { "underflow" }, x.abs())
        }
    };
    if let Some(x) = doc.g.wild() {
        // Where it is written, if it is (not read from a file or made in Rust).
        return Some(match doc.wild_written.as_ref().filter(|w| w.0.to_bits() == x.to_bits()) {
            Some((_, at)) => format!("the number {x:e} written at {at} is {RANGE}"),
            None => what(x),
        });
    }
    let name = |k: usize| doc.params.get(k).map_or("?", |p| p.name.as_str());
    if let Some((k, x)) = ctx.params.iter().enumerate().find(|p| wild(*p.1)) {
        return Some(format!("param '{}' is {x:e}, {RANGE}", name(k)));
    }
    let band = ctx.p_spans.iter().find_map(|s| [s.span.0, s.span.1].into_iter().find(|&x| wild(x)).map(|x| (s.k, x)));
    if let Some((k, x)) = band {
        return Some(format!("param '{}' ranges to {x:e}, {RANGE}", name(k as usize)));
    }
    let (a, b) = ctx.t_span.unwrap_or((ctx.t, ctx.t));
    if let Some(t) = [ctx.t, a, b].into_iter().find(|&t| wild(t)) {
        return Some(format!("the moment {t:e} s is {RANGE}"));
    }
    fields.iter().find_map(|&f| lost_in(doc, &Tape::new(&doc.g, f), ctx)).map(what)
}

/// A value `tape` computes in floats from constants and params alone
/// (Tape::fixed; not a ranged one, bounded outward; not a solved one, whose
/// answers are measured) that floats lost.
fn lost_in(doc: &Doc, tape: &Tape, ctx: &Ctx) -> Option<f64> {
    use ink_core::graph::{BinOp, UnOp, bin_f64, lost, un_f64};
    use ink_core::tape::Op;
    let ranged = |i: usize| tape.pmask[i] != 0 && tape.params_of(i).iter().any(|&k| ctx.p_spans.get(k).is_some() || doc.solved_params.contains(&(k as usize)));
    let mut v = vec![f64::NAN; tape.ssa.len()];
    for (i, ins) in tape.ssa.iter().enumerate().filter(|&(i, _)| tape.fixed[i] && !ranged(i)) {
        let (a, b) = (v[ins.a as usize], v[ins.b as usize]);
        let (x, args, scales) = match ins.op {
            Op::Const(x) => (x, [x, x], false),
            Op::Var(Var::Param(p)) => (ctx.params.get(p as usize).copied().unwrap_or(f64::NAN), [0.0; 2], false),
            Op::Un(op) => (un_f64(op, a), [a, a], matches!(op, UnOp::Square | UnOp::Exp)),
            Op::Bin(op) => (bin_f64(op, a, b), [a, b], matches!(op, BinOp::Mul | BinOp::Div | BinOp::Pow | BinOp::Emax)),
            Op::Select => (if a < 0.0 { b } else { v[ins.c as usize] }, [a, b], false),
            _ => (a, [a, a], false),
        };
        if !matches!(ins.op, Op::Const(_) | Op::Var(_)) && lost(x, &args, scales) {
            return Some(x);
        }
        v[i] = x;
    }
    None
}

fn scaled(far: f64, finest: f64) -> Option<String> {
    if !(far <= REACH) {
        return Some(format!("a coordinate {far:e} mm from the origin is outside ink's range (about 1 nm to 1,000 km)"));
    }
    let floor = FINEST.max(RELATIVE * far);
    (finest < floor).then(|| format!("a length of {finest:e} mm, {far:e} mm from the origin, is outside ink's range (no finer than {floor:e} mm there)"))
}
