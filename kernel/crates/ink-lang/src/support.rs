use super::*;
use crate::claims::Evidence;
use crate::view::ModelView;
use ink_core::interval::{I, add_o, mul_o, sub_o};
use ink_kernel::query::{Proof, Verdict};
use std::sync::Arc;

const LAYER_START: f64 = 31.0 / 32.0;

/// Certify an air gap below the point. Confirmation must enclose the
/// whole square at all depths from LAYER_START of the chosen depth to the full
/// depth. The start is solid, so such a witness has positive thickness:
/// a single zero-valued face where two blocks touch cannot fail.
/// The enclosure includes the centre at full depth for every choice of
/// depth. An emptiness proof still has to exclude every vertical break.
pub(super) fn section(g: &mut Graph, f: NodeId, r: NodeId, depth: NodeId, name: &str) -> NodeId {
    let coordinate = |axis| graph::Fingerprint::default().str("overhang section coordinate").str(name).word(axis).value() as u64;
    let mut offset = |axis: Var, i| {
        let u = g.within(-1.0, 1.0, 0.0, coordinate(i));
        let d = g.bin_real(BinOp::Mul, r, u);
        let p = g.var(axis);
        (axis, g.add(p, d))
    };
    let xy = [offset(Var::X, 0), offset(Var::Y, 1)];
    let fraction = g.within(LAYER_START, 1.0, 1.0, coordinate(2));
    let down = g.bin_real(BinOp::Mul, depth, fraction);
    let z = g.var(Var::Z);
    let q = g.sub(z, down);
    let read = g.subst(f, &[xy[0], xy[1], (Var::Z, q)]);
    let no = g.constant(1.0);
    let yes = g.constant(-1.0);
    g.select(read, no, yes)
}

#[derive(Clone, Copy)]
pub(super) struct Overhang {
    pub(super) f: NodeId,
    pub(super) bed: NodeId,
    pub(super) h: NodeId,
    pub(super) slope: NodeId,
    pub(super) k: u32,
    pub(super) certain: NodeId,
    pub(super) possible: NodeId,
}

impl Overhang {
    /// Look for a witness on the underside nearest the uncovered point.
    /// Projection only suggests positions; outward intervals certify the
    /// material at the start and every point of an air layer below it.
    pub(super) fn witness(self, v: &mut ModelView, candidate: &Proof) -> Option<Evidence> {
        let centre = match candidate.verdict {
            Verdict::Disproven { witness, .. } => witness,
            Verdict::Unknown { cell: Some([lo, hi]), .. } => std::array::from_fn(|i| lo[i] + 0.5 * (hi[i] - lo[i])),
            _ => return None,
        };
        let scope = v.ctx.into_scope();
        let height = v.doc.scalar(self.h, &scope);
        if !(height > 0.0 && height.is_finite()) {
            return None;
        }
        let material = Arc::new(Tape::new(&v.doc.g, self.f));
        let mut gap = GapSearch {
            start: v.clipped(self.f),
            material: &material,
            bed: Tape::new(&v.doc.g, self.bed),
            slope: Tape::new(&v.doc.g, self.slope),
            height: Tape::new(&v.doc.g, self.h),
            slots: Vec::new(),
        };
        let span = scope.p_spans.get(self.k)?;
        let context = |s| {
            if !(span.0..=span.1).contains(&s) {
                return None;
            }
            let ctx = scope.with_param(self.k, (s, s)).into_scope().point();
            ctx.in_scope().then_some(ctx)
        };
        // The coverage search often already found a point inside the roof,
        // even when its thickness is far below the surface sampling grid.
        for i in 0..=32 {
            let s = 0.5_f64.powi(i);
            if let Some(ctx) = context(s)
                && let Some(proof) = gap.at(centre, &ctx, s)
            {
                return Some(proof);
            }
        }
        let reach = v.reach(self.f);
        let local = (centre.map(|x| x - 2.0 * height), centre.map(|x| x + 2.0 * height), 12.0);
        // Try the first uncovered neighbourhood, then other undersides.
        for (lo, hi, cells) in [local, (reach.lo, reach.hi, 24.0)] {
            if ink_kernel::stop::reached() || !(0..3).all(|i| hi[i] > lo[i] && (hi[i] - lo[i]).is_finite()) {
                continue;
            }
            let (points, _) = ink_kernel::feature::surface(&material, &scope, lo, hi, cells);
            let mut points: Vec<_> = points
                .into_iter()
                .filter(|(p, n, _)| {
                    let slope = eval::eval_point(&gap.slope, &scope, *p);
                    n[2] < -0.9 * slope / slope.hypot(1.0) && p[2] >= eval::eval_point(&gap.bed, &scope, *p)
                })
                .collect();
            points.sort_by(|a, b| a.1[2].total_cmp(&b.1[2]).then_with(|| a.0[0].total_cmp(&b.0[0])).then_with(|| a.0[1].total_cmp(&b.0[1])));
            let step = points.len().div_ceil(96).max(1);
            for i in [1, 3, 6, 10, 14] {
                let s = 0.5_f64.powi(i);
                let Some(ctx) = context(s) else { continue };
                for (p, _, _) in points.iter().step_by(step) {
                    let witness = [p[0], p[1], p[2] + height * s / 128.0];
                    if let Some(proof) = gap.at(witness, &ctx, s) {
                        return Some(proof);
                    }
                }
            }
        }
        None
    }
}

struct GapSearch<'a> {
    start: Arc<Tape>,
    material: &'a Tape,
    bed: Tape,
    slope: Tape,
    height: Tape,
    slots: Vec<I>,
}

impl GapSearch<'_> {
    fn at(&mut self, witness: V3, ctx: &Ctx, fraction: f64) -> Option<Evidence> {
        if ink_kernel::stop::reached() {
            return None;
        }
        let at = witness.map(I::point);
        let depth = mul_o(eval::eval_interval(&self.height, ctx, at, None), I::point(fraction));
        let z = sub_o(at[2], mul_o(depth, I::new(LAYER_START, 1.0)));
        if !(depth.lo > 0.0 && z.lo > eval::eval_interval(&self.bed, ctx, at, None).hi) {
            return None;
        }
        let width = mul_o(depth, eval::eval_interval(&self.slope, ctx, at, None));
        if !(width.lo >= 0.0 && width.hi.is_finite()) {
            return None;
        }
        // Cheap point readings can discard a candidate. Only the interval
        // checks below can confirm it (including air whose field equals zero).
        let mid = 0.5 * z.lo + 0.5 * z.hi;
        if eval::eval_point(self.material, ctx, [witness[0], witness[1], mid]) < 0.0 {
            return None;
        }
        if !eval::eval_interval_in(&self.start, ctx, at, None, &mut self.slots).all_inside() {
            return None;
        }
        let mut air = |x, y| {
            let plane = [add_o(at[0], mul_o(width, x)), add_o(at[1], mul_o(width, y)), z];
            eval::eval_interval_in(self.material, ctx, plane, None, &mut self.slots).none_inside()
        };
        clear_disk(&mut air, [-1.0; 2], [1.0; 2], 0).then(|| Evidence::OverhangGap { witness, at: ink_kernel::scope::At::of(ctx), layer_z: [z.lo, z.hi] })
    }
}

/// Cover the disk by boxes. The endpoints here are small exact dyadics:
/// the nearest corner's squared radius is exact, including equality on
/// the circle. A tile outside the disk is irrelevant to every allowed ray.
fn clear_disk(air: &mut impl FnMut(I, I) -> bool, lo: [f64; 2], hi: [f64; 2], level: u32) -> bool {
    let r = [0, 1].map(|i| lo[i].max(-hi[i]).max(0.0));
    if r[0] * r[0] + r[1] * r[1] > 1.0 || air(I::new(lo[0], hi[0]), I::new(lo[1], hi[1])) {
        return true;
    }
    if level == 8 {
        return false;
    }
    let k = level as usize % 2;
    let mid = 0.5 * (lo[k] + hi[k]);
    let (mut left, mut right) = (hi, lo);
    left[k] = mid;
    right[k] = mid;
    clear_disk(air, lo, left, level + 1) && clear_disk(air, right, hi, level + 1)
}
