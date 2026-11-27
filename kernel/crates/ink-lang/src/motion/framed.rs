use super::*;

/// A pair claim searched in a moving frame (`moving_frame`), for its
/// answer in the world: the frame's point as fields of x, y, z and t's
/// param k (`world`), and the claim's shapes in the world (t read as k).
/// A witness found in the frame is a real point of the world claim's
/// violation (the frame is exact in real numbers), but the point reported
/// is its image in floats, which may be none: a frame 1e16 from the origin
/// rounds a witness off both shapes. So every witness reported is shown
/// one again, of the world's own fields at its moment, by bounds rounded
/// outward (`confirmed`); one that is not leaves the claim to the world's
/// search.
#[derive(Clone, Debug)]
pub(crate) struct Frame {
    world: [u32; 3],
    k: u32,
    span: (f64, f64),
    /// The world's fields: what an emptiness proof's witness must be inside
    /// (max(a, b); a subset's max(a, -b), whose witness is inside a and not
    /// inside b: `subset`), and the two shapes.
    empty: u32,
    subset: bool,
    overlap: u32,
    a: u32,
    b: u32,
    /// The distance a pair of points must be within: a clearance's, or a
    /// contact's (the claim's accuracy).
    d: f64,
    /// The scope the frame was searched in (t's param over the span, what
    /// solves on t set derived from it): a witness's point of it is its
    /// world's (`ctx`).
    scope: Ctx,
    /// How far the frame's shapes may be off the world's (`charge`).
    charge: f64,
}

/// Halvings of a witness's cube `Frame::inside` may make on each axis.
const CONFIRM_SPLITS: u32 = 3;

impl Frame {
    pub(super) fn new(v: &mut ModelView, claim: &Claim, world: [u32; 3], k: u32, span: (f64, f64), accuracy: Option<f64>, charge: f64) -> Frame {
        let t = v.doc.g.var(Var::Param(k));
        let world = world.map(|n| v.doc.g.subst(n, &[(Var::T, t)]));
        let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else {
            unreachable!("a pair claim has a frame")
        };
        let (a, b) = (v.doc.g.subst(a, &[(Var::T, t)]), v.doc.g.subst(b, &[(Var::T, t)]));
        let overlap = v.doc.g.max(a, b);
        let empty = match claim {
            Claim::Subset(..) => claims::subset_field(&mut v.doc.g, a, b),
            _ => overlap,
        };
        let d = match *claim {
            Claim::Clearance(_, _, d) => d,
            _ => claims::claim_accuracy(v, claim, accuracy),
        };
        let scope = if v.ctx.p_spans.get(k).is_some() { v.ctx.clone() } else { v.ctx.with_band(k, span) };
        Frame { world, k, span, empty, subset: matches!(claim, Claim::Subset(..)), overlap, a, b, d, scope, charge }
    }

    /// The witness's point of the question's scope (`scope`): the view's
    /// ranged params and t's, each at its witness value (`params`) as a
    /// span of that one value, the params derived from them following
    /// (Ctx::point): what a solve on t sets, at the witness's moment. A ranged param keeps its
    /// meaning there, its expressions' values in real numbers bounded
    /// outward: pinned as a float instead, t at 0.5 reads 1e16 + 0.001 t
    /// as its rounded float, and a pair 0.9 off both shapes shows inside.
    fn ctx(&self, params: &[(u32, f64)]) -> Ctx {
        let mut c = self.scope.clone();
        for &(k, x) in params {
            c.params[k as usize] = x;
        }
        c.point()
    }

    /// A witness's params with t's: where the frame took t out of the
    /// question, the witness holds at every moment of the span in the
    /// frame, and is given at its start.
    fn with_t(&self, params: &[(u32, f64)]) -> Vec<(u32, f64)> {
        let mut ps = params.to_vec();
        if !ps.iter().any(|p| p.0 == self.k) {
            ps.push((self.k, self.span.0));
        }
        ps
    }

    /// p of the frame moved into the world at ctx's moment (floats).
    fn point(&self, v: &ModelView, c: &Ctx, p: V3) -> V3 {
        self.world.map(|n| ink_core::eval::eval_point(&Tape::new(&v.doc.g, n), c, p))
    }

    /// Whether field f is below 0 all over box b at ctx's moment, by bounds
    /// rounded outward (halved up to CONFIRM_SPLITS times where undecided).
    fn inside(tape: &Tape, c: &Ctx, b: [interval::I; 3], splits: u32) -> bool {
        let r = ink_core::eval::eval_interval(tape, c, b, None);
        if r.all_inside() {
            return true;
        }
        if splits == 0 || r.none_inside() || b.iter().all(|x| x.lo == x.hi) {
            return false;
        }
        (0..8).all(|m: usize| {
            let h: [interval::I; 3] = std::array::from_fn(|j| {
                let x = b[j];
                let mid = x.mid();
                if m >> j & 1 == 0 { interval::I::new(x.lo, mid) } else { interval::I::new(mid, x.hi) }
            });
            Frame::inside(tape, c, h, splits - 1)
        })
    }

    /// An emptiness proof's witness (a point every point within `tol` of
    /// which is inside the solid: its cube) in the world, shown so of
    /// field f there; None where it is not. Any other outcome as it is (a
    /// proof holds in real numbers; an undecided box is not one of the
    /// world's: dropped).
    fn proof(&self, v: &ModelView, f: u32, p: ink_kernel::ranged::Proof) -> Option<ink_kernel::ranged::Proof> {
        use ink_kernel::{query::Verdict, ranged::Proof};
        let verdict = match p.verdict {
            Verdict::Disproven { witness, params, .. } => {
                let params = self.with_t(&params);
                let c = self.ctx(&params);
                let w = self.point(v, &c, witness);
                // Asked of shapes off the world's (`charge`), a witness is
                // shown one to the claim's own accuracy.
                let tol = if self.charge > 0.0 { p.tol.max(self.d) } else { p.tol };
                let cube = w.map(|x| if tol > 0.0 { interval::I::new((x - tol).next_down(), (x + tol).next_up()) } else { interval::I::point(x) });
                let tape = Tape::new(&v.doc.g, f);
                let shown = match self.subset && f == self.empty {
                    true => {
                        Frame::inside(&Tape::new(&v.doc.g, self.a), &c, cube, CONFIRM_SPLITS)
                            && ink_core::eval::eval_interval(&Tape::new(&v.doc.g, self.b), &c, cube, None).none_inside()
                    }
                    false => Frame::inside(&tape, &c, cube, CONFIRM_SPLITS),
                };
                if !shown {
                    return None;
                }
                let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                Verdict::Disproven { witness: w, t, value: ink_core::eval::eval_point(&tape, &c, w), params }
            }
            // The accuracy it holds to, off by the charge.
            Verdict::Unknown { cells, volume, ended, attained, hunt, .. } => {
                Verdict::Unknown { cells, volume, ended, cell: None, attained: attained.map(|a| if self.charge > 0.0 { interval::add_up(a, self.charge) } else { a }), hunt }
            }
            o => o,
        };
        Some(Proof { verdict, ..p })
    }

    /// How the claim held, its near point (found in the frame) moved into
    /// the world at its moment and confirmed there on the world's field,
    /// as a witness is (`proof`): outright in the frame, it holds outright
    /// in the world where the frame is exact in real numbers (a charged
    /// one's is not: Held::charged first); a near point the
    /// world does not show is not measured.
    pub(super) fn held(&self, v: &ModelView, h: claims::Held) -> claims::Held {
        let claims::Held::ToAccuracy { near: Some((_, Some((p, at)))) } = h else { return h };
        let params = self.with_t(&at.params);
        let c = self.ctx(&params);
        let w = self.point(v, &c, p);
        let at = ink_kernel::scope::At { t: params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1), params };
        match claims::Held::confirmed_at(false, [(w, at)], &Tape::new(&v.doc.g, self.empty), &c) {
            claims::Held::ToAccuracy { near: Some((_, None)) } => claims::Held::ToAccuracy { near: None },
            h => h,
        }
    }

    /// The evidence with its witnesses in the world, each shown one of the
    /// world's fields at its moment: an emptiness's (or an overlap's) point
    /// inside to its tolerance; a clearance's (or a contact's) pair, each
    /// point in its shape and nearer each other than its distance. None
    /// where one is not.
    pub(super) fn confirmed(&self, v: &ModelView, e: Evidence) -> Option<Evidence> {
        use ink_kernel::query::Verdict;
        Some(match e {
            Evidence::Empty(p) => Evidence::Empty(self.proof(v, self.empty, p)?),
            Evidence::Overlaps(p, x) => Evidence::Overlaps(self.proof(v, self.overlap, p)?, x),
            Evidence::ClearanceSearch { verdict, pair, res } => {
                let (verdict, pair) = match verdict {
                    Verdict::Disproven { params, .. } => {
                        let params = self.with_t(&params);
                        let c = self.ctx(&params);
                        let (p, q) = pair?;
                        // One in each shape (the claim's way round).
                        let (p, q) = (self.point(v, &c, p), self.point(v, &c, q));
                        let at = |f: u32, x: V3| Frame::inside(&Tape::new(&v.doc.g, f), &c, x.map(interval::I::point), 0);
                        let each = at(self.b, p) && at(self.a, q) || at(self.a, p) && at(self.b, q);
                        if !(each && ink_kernel::query::distance_up(p, q) < self.d) {
                            return None;
                        }
                        let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                        let value = ink_kernel::query::distance_up(p, q);
                        (Verdict::Disproven { witness: p, t, value, params }, Some((p, q)))
                    }
                    // The clearance it holds to, less the charge.
                    Verdict::Unknown { cells, volume, ended, attained, hunt, .. } => {
                        let attained = attained.map(|a| if self.charge > 0.0 { -interval::add_up(-a, self.charge) } else { a });
                        (Verdict::Unknown { cells, volume, ended, cell: None, attained, hunt }, None)
                    }
                    verdict => (verdict, pair),
                };
                Evidence::ClearanceSearch { verdict, pair, res }
            }
            Evidence::Touches { overlap, gap } => Evidence::Touches { overlap: Box::new(self.confirmed(v, *overlap)?), gap: Box::new(self.confirmed(v, *gap)?) },
            e => e,
        })
    }
}
