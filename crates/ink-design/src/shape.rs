//! Shape derivatives for the optimizers: how boundary-weighted integrals
//! (volume, mass, compliance) change with each param.

use ink_core::V3;
use ink_core::eval::{Ctx, eval_dual, eval_points};
use ink_core::sign::inside;
use ink_core::tape::Tape;
use ink_kernel::octree::{descend, root_for};
use std::sync::Arc;

/// Shape derivatives of boundary-weighted integrals: for each weight w,
/// dQ/dp = -∮ w (∂f/∂p)/|∇f| dA over the part's surface, for every param p
/// (the boundary moves by -(∂f/∂p)/|∇f| along its normal). Weight 1 gives
/// the volume, a density the mass, and minus the strain energy density σ:ε
/// the compliance (Hadamard's formula, for boundaries free of loads and
/// supports). Summed over the crossings between inside and outside on a
/// lattice of samples (`ink_core::sign`), with exact field derivatives
/// (forward-mode duals) just inside each. Returns the surface area and one
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
        &|_: &Arc<Tape>, c| {
            // A lattice 4 a side in the cell, one more each way: the
            // segments between neighbours on it where one end is inside and
            // the other not, each halved toward its crossing. A crossing
            // along axis k stands for s^2 of the surface's area projected
            // across k, |n_k| dA, so s^2 |n_k| per crossing sums to the
            // area over the three axes (n_k^2 adds to 1). Each is counted
            // in the cell holding it (half-open), read by the field's whole
            // tape (the lattice reaches past the cell), its derivatives
            // just inside it: a field reading 0 outside, with no slope
            // there, still moves with its inside.
            let s = grid.h * c.size as f64 / 4.0;
            let (clo, chi) = grid.cell_bounds(c);
            let at = |i: [usize; 3]| -> V3 { std::array::from_fn(|k| clo[k] + s * (i[k] as f64 - 0.5)) };
            let idx = |i: [usize; 3]| (i[2] * 6 + i[1]) * 6 + i[0];
            let mut pts = Vec::with_capacity(216);
            for z in 0..6 {
                for y in 0..6 {
                    for x in 0..6 {
                        pts.push(at([x, y, z]));
                    }
                }
            }
            let mut v = Vec::new();
            eval_points(tape, ctx, &pts, &mut v);
            // (inside end, outside end, axis)
            let mut br: Vec<(V3, V3, usize)> = Vec::new();
            for z in 0..6 {
                for y in 0..6 {
                    for x in 0..6 {
                        let i = [x, y, z];
                        for k in 0..3 {
                            if i[k] == 5 {
                                continue;
                            }
                            let mut j = i;
                            j[k] += 1;
                            let (a, b) = (idx(i), idx(j));
                            match (inside(v[a]), inside(v[b])) {
                                (true, false) => br.push((pts[a], pts[b], k)),
                                (false, true) => br.push((pts[b], pts[a], k)),
                                _ => {}
                            }
                        }
                    }
                }
            }
            if br.is_empty() {
                return None;
            }
            let mut mid = Vec::with_capacity(br.len());
            for _ in 0..12 {
                mid.clear();
                mid.extend(br.iter().map(|(p, q, _)| -> V3 { std::array::from_fn(|k| 0.5 * (p[k] + q[k])) }));
                eval_points(tape, ctx, &mid, &mut v);
                for ((b, m), &f) in br.iter_mut().zip(&mid).zip(&v) {
                    if inside(f) {
                        b.0 = *m;
                    } else {
                        b.1 = *m;
                    }
                }
            }
            let (mut area, mut dq) = (0.0, vec![vec![0.0; np]; nw]);
            for (p, q, k) in br {
                let c: V3 = std::array::from_fn(|a| 0.5 * (p[a] + q[a]));
                if !(0..3).all(|a| clo[a] <= c[a] && c[a] < chi[a]) {
                    continue;
                }
                let d = eval_dual(tape, ctx, p);
                let gl = (d[1] * d[1] + d[2] * d[2] + d[3] * d[3]).sqrt();
                if !(gl > 0.0 && gl.is_finite()) {
                    continue;
                }
                let da = s * s * d[1 + k].abs() / gl;
                area += da;
                for (wf, g) in weights.iter().zip(dq.iter_mut()) {
                    let wv = wf(c);
                    for j in 0..np {
                        g[j] -= da * wv * d[5 + j] / gl;
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
