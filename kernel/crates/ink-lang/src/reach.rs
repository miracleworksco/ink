//! Where analysis of a field looks (Reach): the window, a field's box
//! proven within it, and clipping to it.

use super::*;

/// Where analysis of a field looks (see `Doc::reach`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reach {
    pub lo: V3,
    pub hi: V3,
    /// The box holds all of the field's material. If not (the field is not
    /// proven bounded) it is the window, and a result that needs all of the
    /// shape (emptiness, an upper bound on volume, a box, a distance's lower
    /// bound) is unknown; one from a witness inside still holds.
    pub whole: bool,
}

impl Reach {
    /// Refuse a whole-shape answer if this box is only the fallback window.
    pub fn require_whole(&self, target: &str) -> Result<(), InkError> {
        if self.whole {
            return Ok(());
        }
        Err(InkError {
            hint: Some("intersect the field with a bounded shape that contains the entire intended part, e.g. (max f (box ...))".into()),
            ..InkError::msg(format!("{target} is not proven bounded: it may extend outside the 1 m window at the origin"))
        })
    }
}

impl Doc {
    /// The window: where fields that are not bounded are examined (and
    /// viewed). Bounded shapes are analysed whole wherever they lie.
    pub fn domain(&self) -> (V3, V3) {
        DEFAULT_DOMAIN
    }

    /// The finest cell a proof refines to unless given a resolution: half
    /// the accuracy (ACCURACY of the window), the cell that proves touching
    /// surfaces to it: the same for every model, so a part's claims look as
    /// fine in an assembly as the assembly's own.
    pub fn finest(&self) -> f64 {
        let (lo, hi) = self.domain();
        ACCURACY / 2.0 * crate::query::max_extent(lo, hi)
    }

    /// Where analysis of field `f` looks (over `ctx.t_span`, if set, at
    /// every moment of it):
    /// a box proven to hold all of its material (`query::extent`, a search
    /// over all of space), or, if the field is not proven bounded, the
    /// window (`domain`), and then results that need all of the shape say so.
    /// A shape inside the window gets the window itself.
    pub fn reach(&self, f: NodeId, ctx: &Ctx) -> Reach {
        reach_in(&self.g, self.domain(), f, ctx)
    }

    /// Field that is negative strictly inside an axis-aligned box.
    /// A side at infinity clips nothing, and is no term.
    pub fn box_field(&mut self, lo: V3, hi: V3) -> NodeId {
        let vars = [Var::X, Var::Y, Var::Z];
        let mut terms = Vec::new();
        for i in 0..3 {
            let v = self.g.var(vars[i]);
            if lo[i] > f64::NEG_INFINITY {
                let c = self.g.constant(lo[i]);
                terms.push(self.g.sub(c, v));
            }
            if hi[i] < f64::INFINITY {
                let c = self.g.constant(hi[i]);
                terms.push(self.g.sub(v, c));
            }
        }
        self.reduce(BinOp::Max, &terms)
    }

    /// `box_field` with corners that are fields (of constants and params).
    pub(super) fn box_nodes(&mut self, lo: [NodeId; 3], hi: [NodeId; 3]) -> NodeId {
        let vars = [Var::X, Var::Y, Var::Z];
        let mut terms = Vec::new();
        for i in 0..3 {
            let v = self.g.var(vars[i]);
            terms.push(self.g.sub(lo[i], v));
            terms.push(self.g.sub(v, hi[i]));
        }
        self.reduce(BinOp::Max, &terms)
    }

    /// Intersect a field with a box (keeps every query bounded and closed).
    pub fn clip(&mut self, f: NodeId, lo: V3, hi: V3) -> NodeId {
        let b = self.box_field(lo, hi);
        self.g.max(f, b)
    }

    /// `f` clipped to where analysis looks (see `reach`), and that region.
    /// Clipping a bounded field to its reach changes no sign: it has no
    /// material outside.
    pub fn clip_reach(&mut self, f: NodeId, ctx: &Ctx) -> (NodeId, Reach) {
        let r = self.reach(f, ctx);
        (self.clip(f, r.lo, r.hi), r)
    }
}

/// Doc::reach of `f` in graph `g` with window (`wlo`, `whi`).
pub fn reach_in(g: &ink_core::graph::Graph, window: (V3, V3), f: NodeId, ctx: &Ctx) -> Reach {
    reach_of(&std::sync::Arc::new(Tape::new(g, f)), window, ctx)
}

/// Where analysis of a field that reads no z (the same at every height)
/// looks in the plane: in x and y, the reach of `slab` (the field between
/// two heights); in z, a band about 0 as deep as its material is wide. The
/// field is clipped in x and y alone (ModelView::clipped), so no search
/// halves the band (nothing in it changes along z), and a point of it as
/// deep as any can be (no deeper than its material is wide) is as deep in
/// the band.
pub fn plane_reach(g: &ink_core::graph::Graph, window: (V3, V3), slab: NodeId, ctx: &Ctx) -> Reach {
    let tape = std::sync::Arc::new(Tape::new(g, slab));
    let r = reach_of(&tape, window, ctx);
    // Its material's box (none: empty, any band will do).
    let (lo, hi) = crate::query::bounds(&tape, ctx, r.lo, r.hi).map_or((r.lo, r.hi), |b| (b.0, b.1));
    let half = 0.5 * (hi[0] - lo[0]).max(hi[1] - lo[1]);
    Reach { lo: [r.lo[0], r.lo[1], -half], hi: [r.hi[0], r.hi[1], half], whole: r.whole }
}

/// Halvings of a ranged param's span `extent_in_pieces` may make.
const REACH_PIECES: u32 = 6;

/// The extent of `tape` over ctx's scope as the union of its extents over
/// pieces of it: each ranged param it reads halved (the widest first, its
/// halves each found again) where a piece's is unbounded, up to `depth`
/// halvings. Over a whole turn a turned part's far points bound as any
/// direction, so no direction is shown outside (unbounded); over pieces of
/// the turn each is. The union holds every piece's material, so all of
/// it. Unbounded where a piece still is at the last halving, or no param
/// is read.
fn extent_in_pieces(tape: &std::sync::Arc<Tape>, ctx: &Ctx, window: (V3, V3), depth: u32) -> crate::query::Extent {
    use crate::query::Extent;
    let ks: Vec<u32> = ctx
        .p_spans
        .iter()
        .filter(|s| !s.derived && s.span.0 < s.span.1 && tape.ssa.iter().any(|i| matches!(i.op, ink_core::tape::Op::Var(Var::Param(k)) if k == s.k)))
        .map(|s| s.k)
        .collect();
    let Some(&k) = ks.iter().max_by(|&&a, &&b| {
        let w = |k: u32| ctx.p_spans.get(k).map_or(0.0, |s| s.1 - s.0);
        w(a).total_cmp(&w(b))
    }) else {
        return Extent::Unbounded;
    };
    if depth == 0 {
        return Extent::Unbounded;
    }
    let (lo, hi) = ctx.p_spans.get(k).expect("a ranged param");
    let m = ink_core::eval::span_mid((lo, hi));
    if !(lo < m && m < hi) {
        return Extent::Unbounded;
    }
    let mut out = Extent::Empty;
    for half in [(lo, m), (m, hi)] {
        let c = ctx.with_param(k, half);
        let e = match crate::query::extent(tape, &c, window) {
            Extent::Unbounded => extent_in_pieces(tape, &c, window, depth - 1),
            e => e,
        };
        out = match (out, e) {
            (Extent::Unbounded, _) | (_, Extent::Unbounded) => return Extent::Unbounded,
            (Extent::Empty, e) => e,
            (o, Extent::Empty) => o,
            (Extent::Bounded(a, b), Extent::Bounded(c, d)) => {
                let (lo, hi) = ink_core::v3::box_union((a, b), (c, d));
                Extent::Bounded(lo, hi)
            }
        };
    }
    out
}

/// Whether the first of the finest pieces `extent_in_pieces` may reach (each
/// param it reads at the first 2^-depth of its span) has a bounded extent:
/// where even that one has none, no union of pieces does, and the search
/// for one (each extent's own up to its largest budget, down every halving)
/// costs as much as the claim (a belt's shape driven along its loop).
fn finest_piece_bounded(tape: &std::sync::Arc<Tape>, ctx: &Ctx, window: (V3, V3), depth: u32) -> bool {
    let mut c = ctx.clone();
    let mut any = false;
    for s in ctx.p_spans.iter().filter(|s| !s.derived && s.span.0 < s.span.1 && tape.ssa.iter().any(|i| matches!(i.op, ink_core::tape::Op::Var(Var::Param(k)) if k == s.k))) {
        let (lo, hi) = s.span;
        c = c.with_param(s.k, (lo, lo + (hi - lo) * 2f64.powi(-(depth as i32))));
        any = true;
    }
    any && !matches!(crate::query::extent(tape, &c, window), crate::query::Extent::Unbounded)
}

/// Doc::reach of the field `tape` computes, with window (`wlo`, `whi`).
pub fn reach_of(tape: &std::sync::Arc<Tape>, (wlo, whi): (V3, V3), ctx: &Ctx) -> Reach {
    let extent = match crate::query::extent(tape, ctx, (wlo, whi)) {
        crate::query::Extent::Unbounded if finest_piece_bounded(tape, ctx, (wlo, whi), REACH_PIECES) => extent_in_pieces(tape, ctx, (wlo, whi), REACH_PIECES),
        e => e,
    };
    match extent {
        crate::query::Extent::Bounded(lo, hi) if (0..3).any(|k| lo[k] < wlo[k] || hi[k] > whi[k]) => Reach { lo, hi, whole: true },
        // No material anywhere, or all of it inside the window.
        crate::query::Extent::Empty | crate::query::Extent::Bounded(..) => Reach { lo: wlo, hi: whi, whole: true },
        crate::query::Extent::Unbounded => Reach { lo: wlo, hi: whi, whole: false },
    }
}
