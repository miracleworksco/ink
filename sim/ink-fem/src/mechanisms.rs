//! The rigid motions the supports leave free, and a study whose loads drive
//! one.

use crate::eigen::dense_sym_eigen;
use ink_core::{V3, v3};

/// The failure of a study whose loads drive a free rigid motion (what it
/// is, the loads' work along it per unit, and that work's unit): no
/// static solution.
pub(crate) fn driven_failure(what: &str, work: f64, unit: &str) -> String {
    format!(
        "the supports leave the part free to {what}, and the loads drive it ({} {unit}): a mechanism, \
         with no static solution. Hold that motion (another fix, or one without the free motion), or balance the loads",
        (work * 1e4).round() / 1e4,
    )
}

/// The six rigid motions of a piece at a point: translations along x, y, z,
/// and turns about x, y, z through its centroid, scaled by its size so
/// none moves a point of it more than 1.
pub(crate) fn rigid_mode(p: V3, cent: V3, size: f64, r: usize) -> V3 {
    let mut out = [0.0; 3];
    if r < 3 {
        out[r] = 1.0;
    } else {
        let q: V3 = std::array::from_fn(|k| (p[k] - cent[k]) / size);
        let (a, b) = ((r - 2) % 3, (r - 1) % 3);
        // e_(r-3) x q
        out[a] = -q[b];
        out[b] = q[a];
    }
    out
}

/// The rigid motions the supports leave free, as a description and their
/// displacement per unit (1 rad for a turn, 1 mm for a slide) over the
/// dofs. Per piece (`comp_of` per node, its centroid and size), `blocks`
/// is the supports' penalty energy over its six rigid motions (6 x 6) and
/// `scale` that of a translation were every support held on every axis;
/// the low-rank terms add theirs, coupling the pieces a support's free
/// motions span. The elastic stiffness gives a rigid motion no energy, so a
/// motion is free when the supports give it none: under 1e-9 of the scale
/// (a support's own weakest hold, a thin band turning about its axis, is
/// about its (thickness / size)^2, and rounding 1e-16).
#[allow(clippy::too_many_arguments)]
pub(crate) fn free_rigid_motions(
    pos: &[V3],
    comp_of: &[usize],
    cent: &[V3],
    size: &[f64],
    blocks: &[[f64; 36]],
    scale: &[f64],
    lowrank: &[(f64, Vec<f64>)],
) -> Vec<(String, Vec<f64>, bool)> {
    let nc = cent.len();
    // Each low-rank term over the pieces' rigid motions, and the pieces it
    // couples.
    let terms: Vec<(f64, Vec<f64>)> = lowrank
        .iter()
        .map(|(w, g)| {
            let mut t = vec![0.0; 6 * nc];
            for (n, p) in pos.iter().enumerate() {
                let gn = &g[n * 3..n * 3 + 3];
                if gn.iter().any(|v| *v != 0.0) {
                    let c = comp_of[n];
                    for r in 0..6 {
                        let m = rigid_mode(*p, cent[c], size[c], r);
                        t[c * 6 + r] += (0..3).map(|d| gn[d] * m[d]).sum::<f64>();
                    }
                }
            }
            (*w, t)
        })
        .collect();
    let mut group: Vec<usize> = (0..nc).collect();
    fn root(g: &mut [usize], mut x: usize) -> usize {
        while g[x] != x {
            g[x] = g[g[x]];
            x = g[x];
        }
        x
    }
    for (_, t) in &terms {
        let touched: Vec<usize> = (0..nc).filter(|&c| t[c * 6..c * 6 + 6].iter().any(|v| *v != 0.0)).collect();
        for w in touched.windows(2) {
            let (a, b) = (root(&mut group, w[0]), root(&mut group, w[1]));
            group[a] = b;
        }
    }
    let mut out = Vec::new();
    for g0 in 0..nc {
        let members: Vec<usize> = (0..nc).filter(|&c| root(&mut group, c) == g0).collect();
        if members.is_empty() {
            continue;
        }
        let dim = 6 * members.len();
        // The energy, scaled per piece by its translation scale.
        let mut a = vec![0.0; dim * dim];
        let at = |i: usize| members[i / 6] * 6 + i % 6;
        for i in 0..dim {
            for j in 0..dim {
                let (ci, cj) = (members[i / 6], members[j / 6]);
                let mut v = if ci == cj { blocks[ci][(i % 6) * 6 + j % 6] } else { 0.0 };
                for (w, t) in &terms {
                    v += w * t[at(i)] * t[at(j)];
                }
                a[i * dim + j] = v / (scale[ci] * scale[cj]).sqrt();
            }
        }
        let (vals, vecs) = dense_sym_eigen(&a, dim);
        let null: Vec<Vec<f64>> = (0..dim).filter(|&k| vals[k] < 1e-9).map(|k| (0..dim).map(|i| vecs[i * dim + k]).collect()).collect();
        if null.is_empty() {
            continue;
        }
        // A basis that separates slides from turns: combinations with no
        // turn first (the null space of the turn parts), then the turns,
        // with the slides' share taken out.
        let m = null.len();
        let turn_part = |x: &[f64]| -> Vec<f64> { (0..dim).filter(|i| i % 6 >= 3).map(|i| x[i]).collect() };
        let mut tt = vec![0.0; m * m];
        for i in 0..m {
            for j in 0..m {
                tt[i * m + j] = turn_part(&null[i]).iter().zip(turn_part(&null[j])).map(|(a, b)| a * b).sum();
            }
        }
        let (tv, tvec) = dense_sym_eigen(&tt, m);
        let combo = |k: usize| -> Vec<f64> { (0..dim).map(|i| (0..m).map(|j| tvec[j * m + k] * null[j][i]).sum()).collect() };
        let slides: Vec<Vec<f64>> = (0..m).filter(|&k| tv[k] < 1e-10).map(combo).collect();
        let mut basis = slides.clone();
        for k in (0..m).filter(|&k| tv[k] >= 1e-10) {
            let mut x = combo(k);
            for s in &slides {
                let c: f64 = x.iter().zip(s).map(|(a, b)| a * b).sum();
                x.iter_mut().zip(s).for_each(|(a, b)| *a -= c * b);
            }
            basis.push(x);
        }
        for (mode, x) in basis.into_iter().enumerate() {
            // Back to amplitudes of the rigid motions, then of the piece that
            // moves most: its translation of the centroid and its turn.
            let mut amp: Vec<f64> = (0..dim).map(|i| x[i] / scale[members[i / 6]].sqrt()).collect();
            let (k, _) = (0..members.len()).map(|k| (k, (0..6).map(|r| amp[k * 6 + r].powi(2)).sum::<f64>())).fold((0, -1.0), |a, b| if b.1 > a.1 { b } else { a });
            let c = members[k];
            // An eigenvector's sign is arbitrary: the motion is said along
            // its direction (the turn's axis, else the slide's) with the
            // largest component positive.
            let part = |o: usize| -> V3 { std::array::from_fn(|d| amp[k * 6 + o + d]) };
            let dir = if v3::norm(part(3)) > 1e-6 * v3::norm(part(0)) { part(3) } else { part(0) };
            if dir[(0..3).fold(0, |m, d| if dir[d].abs() > dir[m].abs() { d } else { m })] < 0.0 {
                amp.iter_mut().for_each(|v| *v = -*v);
            }
            let t: V3 = std::array::from_fn(|d| amp[k * 6 + d]);
            let w: V3 = std::array::from_fn(|d| amp[k * 6 + 3 + d] / size[c]);
            let (tl, wl) = (v3::norm(t), v3::norm(w));
            let r6 = |v: f64| (v * 1e4).round() / 1e4 + 0.0;
            let f3 = |v: V3| format!("({}, {}, {})", r6(v[0]), r6(v[1]), r6(v[2]));
            // Per unit: a turn per radian, a slide per mm.
            let unit = if wl * size[c] > 1e-6 * tl { wl } else { tl };
            let mut what = if wl * size[c] > 1e-6 * tl {
                let axis = w.map(|v| v / wl);
                // The axis's point nearest the origin (w x u0 / |w|^2, u0 the
                // motion there), and any pitch along it.
                let wc = v3::cross(w, cent[c]);
                let u0: V3 = std::array::from_fn(|d| t[d] - wc[d]);
                let p = v3::cross(w, u0).map(|x| x / (wl * wl));
                let pitch = (0..3).map(|d| w[d] * u0[d]).sum::<f64>() / (wl * wl);
                let mut s = format!("turn about the axis along {} through {}", f3(axis), f3(p));
                if pitch.abs() > 1e-6 * size[c] {
                    s += &format!(", advancing {} mm per radian along it", r6(pitch));
                }
                s
            } else {
                format!("slide along {}", f3(t.map(|v| v / tl)))
            };
            if nc > 1 {
                what += &format!(" (the piece about {})", f3(cent[c]));
            }
            let mut field = vec![0.0; pos.len() * 3];
            for (n, p) in pos.iter().enumerate() {
                if let Some(kk) = members.iter().position(|&cc| cc == comp_of[n]) {
                    let cc = members[kk];
                    for r in 0..6 {
                        let mo = rigid_mode(*p, cent[cc], size[cc], r);
                        for d in 0..3 {
                            field[n * 3 + d] += amp[kk * 6 + r] / unit * mo[d];
                        }
                    }
                }
            }
            out.push((what, field, mode >= slides.len()));
        }
    }
    out
}
