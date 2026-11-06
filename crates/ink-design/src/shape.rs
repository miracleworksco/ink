//! Shape derivatives for the optimizers: how boundary-weighted integrals
//! (volume, mass, compliance) change with each param.

use ink_core::V3;
use ink_core::eval::{Ctx, eval_dual, eval_grads, eval_points};
use ink_core::tape::Tape;
use ink_kernel::octree::{descend, root_for};
use ink_kernel::query::samples;
use std::sync::Arc;

/// Shape derivatives of boundary-weighted integrals: for each weight w,
/// dQ/dp = -∮ w (∂f/∂p)/|∇f| dA over the part's surface, for every param p
/// (the boundary moves by -(∂f/∂p)/|∇f| along its normal). Weight 1 gives
/// the volume, a density the mass, and minus the strain energy density σ:ε
/// the compliance (Hadamard's formula, for boundaries free of loads and
/// supports). Evaluated as a band integral: each sample's partial occupancy
/// clamp(1/2 - f/w) moves with the field by -(∂f/∂p)/w, with exact field
/// derivatives (forward-mode duals). Returns the surface area and one
/// gradient per weight.
pub type Weight<'a> = &'a (dyn Fn(V3) -> f64 + Sync);

pub fn shape_derivatives(tape: &Arc<Tape>, weights: &[Weight], ctx: &Ctx, lo: V3, hi: V3, h: f64) -> (f64, Vec<Vec<f64>>) {
    let np = ctx.params.len();
    let nw = weights.len();
    let (grid, root) = root_for(lo, hi, h, 1);
    // (area, gradient per weight) per surface cell.
    let parts: Vec<(f64, Vec<Vec<f64>>)> = descend(
        tape,
        ctx,
        &grid,
        root,
        1,
        &|t: &Arc<Tape>, c| {
            let vol = (grid.h * c.size as f64).powi(3);
            // 4 samples per axis and a band as wide as their spacing: the band
            // then stays within cells the surface crosses (at most 0.87 of a
            // spacing wide, samples an eighth of a cell from the walls), so
            // no sample the moving surface sweeps is missed by the octree.
            let pts = samples(&grid, c, 4);
            let fg = eval_grads(t, ctx, &pts);
            let s = grid.h * c.size as f64 / 4.0;
            let per = vol / pts.len() as f64;
            let (mut area, mut dq) = (0.0, vec![vec![0.0; np]; nw]);
            // Band samples (half-open, so samples exactly at its edges count
            // once), with probes one band-width out and in along the normal.
            let mut band = Vec::new();
            let mut probes = Vec::new();
            for (p, f) in pts.iter().zip(&fg) {
                let w = s * (f[1].abs() + f[2].abs() + f[3].abs());
                let occ = 0.5 - f[0] / w;
                if !(w > 0.0 && (0.0..1.0).contains(&occ)) {
                    continue;
                }
                let gl = (f[1] * f[1] + f[2] * f[2] + f[3] * f[3]).sqrt();
                let step = w / gl;
                let n = [f[1] / gl, f[2] / gl, f[3] / gl];
                band.push((*p, *f));
                probes.push(std::array::from_fn(|k| p[k] + step * n[k]));
                probes.push(std::array::from_fn(|k| p[k] - step * n[k]));
            }
            let mut pv = Vec::new();
            eval_points(t, ctx, &probes, &mut pv);
            for (i, (p, f)) in band.iter().enumerate() {
                // A real boundary has outside on one side and inside on the
                // other; skip sheets of zero thickness (coincident faces of
                // touching parts in a union, which are interior) and slivers
                // thinner than the band.
                if !(pv[2 * i] > 0.0 && pv[2 * i + 1] < 0.0) {
                    continue;
                }
                let w = s * (f[1].abs() + f[2].abs() + f[3].abs());
                let d = eval_dual(t, ctx, *p);
                let gl = (f[1] * f[1] + f[2] * f[2] + f[3] * f[3]).sqrt();
                // The band of partial occupancy is w/|grad f| thick.
                area += per * gl / w;
                for (wf, g) in weights.iter().zip(dq.iter_mut()) {
                    let wv = wf(*p);
                    for k in 0..np {
                        g[k] -= per * d[5 + k] / w * wv;
                    }
                }
            }
            Some((area, dq))
        },
        &|_: &Arc<Tape>, _| None,
    );
    let mut area = 0.0;
    let mut total = vec![vec![0.0; np]; nw];
    for (a, dq) in parts {
        area += a;
        for (t, g) in total.iter_mut().zip(dq) {
            t.iter_mut().zip(g).for_each(|(t, v)| *t += v);
        }
    }
    (area, total)
}
