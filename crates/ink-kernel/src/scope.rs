//! The scope beyond space: ranged params (a context's `p_spans`). A proof
//! covers every param value in it; searches halve its axes where they keep
//! an answer undecided, as the octree halves space. A witness is a point of
//! it (`At`). A span of time is an axis of the scope too, where a field
//! reads it (`spans_t`).

use ink_core::eval::{Ctx, span_mid};
use ink_core::graph::Var;
use ink_core::tape::{Op, Tape};

/// Halvings of a ranged param's band one search path may make (down to
/// 1/256 of it).
pub const PARAM_SPLITS: u32 = 8;

/// A point of the scope beyond space: a moment, and each ranged param's
/// value (none when no param is ranged).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct At {
    pub t: f64,
    pub params: Vec<(u32, f64)>,
}

impl At {
    /// ctx's point: its moment and its ranged params' values.
    pub fn of(ctx: &Ctx) -> At {
        At { t: ctx.t, params: ctx.p_spans.values(&ctx.params) }
    }
}

/// The ranged params `tape` reads, as axes of the scope a search halves: a
/// derived one (a solve's unknown) stands for the params it is derived
/// from, which narrow it (ParamSpans).
pub fn read_params(tape: &Tape, ctx: &Ctx) -> Vec<u32> {
    if ctx.p_spans.is_empty() {
        return Vec::new();
    }
    let mut ks: Vec<u32> = tape
        .ssa
        .iter()
        .filter_map(|i| match i.op {
            Op::Var(Var::Param(k)) if ctx.p_spans.get(k).is_some() => Some(k),
            _ => None,
        })
        .flat_map(|k| ctx.p_spans.sources(k))
        .collect();
    ks.sort_unstable();
    ks.dedup();
    ks
}

/// Whether a span can be halved: wider than 2^-depth of `whole`, with its
/// middle strictly inside (a span of two adjacent floats cannot).
pub fn can_halve(span: (f64, f64), whole: (f64, f64), depth: u32) -> bool {
    let m = span_mid(span);
    span.0 < m && m < span.1 && (span.1 - span.0) * 2f64.powi(depth as i32) > whole.1 - whole.0
}

/// Whether `tape` reads time and ctx spans more than one moment of it, at
/// a middle strictly inside: time is then an axis of the scope as a ranged
/// param is (bounds first order in it, halved where it blurs them most).
pub fn spans_t(tape: &Tape, ctx: &Ctx) -> bool {
    ctx.t_span.is_some_and(|s| s.0 < span_mid(s) && span_mid(s) < s.1) && tape.ssa.iter().any(|i| matches!(i.op, Op::Var(Var::T)))
}

/// The middle of ctx's scope: each ranged param at the middle of its span
/// (Ctx::param_mids), and with `t`, time at the middle of its span, a
/// moment.
pub fn middle(ctx: &Ctx, t: bool) -> Ctx {
    let c = ctx.param_mids();
    match ctx.t_span.filter(|_| t) {
        Some(s) => Ctx { t: span_mid(s), t_span: None, ..c },
        None => c,
    }
}

/// The ranged params `tape` reads whose spans can still be halved (to
/// 2^-depth of their bands).
pub fn splittable(tape: &Tape, ctx: &Ctx, depth: u32) -> Vec<u32> {
    halvable(&read_params(tape, ctx), ctx, depth)
}

/// Those of params `ks` whose spans in ctx can still be halved (to 2^-depth
/// of their bands).
pub fn halvable(ks: &[u32], ctx: &Ctx, depth: u32) -> Vec<u32> {
    ks.iter().copied().filter(|&k| ctx.p_spans.iter().any(|s| s.k == k && can_halve(s.span, s.band, depth))).collect()
}

/// How many times param k's span in ctx has been halved from its band
/// (log2 of their widths' ratio).
pub fn halvings(ctx: &Ctx, k: u32) -> f64 {
    let s = ctx.p_spans.iter().find(|s| s.k == k).expect("a ranged param");
    ((s.band.1 - s.band.0) / (s.span.1 - s.span.0)).log2()
}

/// ctx with param k's span halved: the lower and the upper half, closed
/// (they share the middle), the param at each one's middle.
pub fn halves(ctx: &Ctx, k: u32) -> (Ctx, Ctx) {
    let (lo, hi) = ctx.p_spans.get(k).expect("a ranged param");
    let m = span_mid((lo, hi));
    (ctx.with_param(k, (lo, m)), ctx.with_param(k, (m, hi)))
}

/// ctx with each ranged param moved to
/// the end of its span that lowers the sum of `tapes` at `p` (by the sign
/// of its derivative there; spans unchanged), again from there while that
/// moves one (a max follows its winning branch, so the other's inputs show
/// only once it wins), at most once per input and once more: a witness
/// search, which at inputs entering the fields once and monotonically
/// reaches their worst values for that point, where the middle values show
/// nothing (a part moving with two params, at the end of both). A
/// witness found there is confirmed like any other. With each tape's value
/// at `p` there (in floats, as eval_point's). None where no derivative
/// moves an input.
pub fn pushed(tapes: &[&Tape], ctx: &Ctx, p: ink_core::V3) -> Option<(Ctx, Vec<f64>)> {
    let mut c = ctx.clone();
    let mut moved = false;
    // Derivatives along params a solve's unknowns are derived from count
    // what the unknowns pass on (ParamSpans::chain).
    let duals = |c: &Ctx| -> Vec<Vec<f64>> {
        tapes
            .iter()
            .map(|t| {
                let mut d = ink_core::eval::eval_dual(t, c, p);
                c.p_spans.chain(c, &mut d[5..]);
                d
            })
            .collect()
    };
    let mut ds = duals(&c);
    // To the end of `span` that the derivative `dk` lowers toward.
    let end = |span: (f64, f64), dk: f64| {
        if dk > 0.0 {
            Some(span.0)
        } else if dk < 0.0 {
            Some(span.1)
        } else {
            None
        }
    };
    for _ in 0..=ctx.p_spans.iter().count() {
        let mut again = false;
        for s in ctx.p_spans.iter().filter(|s| !s.derived) {
            let Some(v) = end(s.span, ds.iter().map(|d| d[5 + s.k as usize]).sum()) else { continue };
            again |= v != c.params[s.k as usize];
            c.params[s.k as usize] = v;
        }
        if !again {
            break;
        }
        moved = true;
        // Derived params follow (their values at the params' new ones).
        c.derive();
        ds = duals(&c);
    }
    moved.then(|| (c, ds.iter().map(|d| d[0]).collect()))
}
