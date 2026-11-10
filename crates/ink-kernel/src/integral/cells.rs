//! Integrals estimated on cells of one size, halved until a halving
//! moves the first under the accuracy asked: where `integrate`'s bounds
//! are out of reach (a field that jumps, edges everywhere), at the cost
//! the estimate takes.

use super::{MOMENTS, Measured};
use crate::octree::{Cell, Grid, descend, root_for};
use crate::query::{max_extent, occupancy};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_points, span_mid};
use ink_core::interval::{I, add_o, mul_o, sub_o};
use ink_core::tape::Tape;
use rayon::prelude::*;
use std::sync::Arc;

/// Cells a halving may make: one that would make more is not made
/// (surfaces over a large area at fine cells: a few seconds' work).
pub const PASS_CELLS: usize = 1 << 21;

/// Integrals at one cell size, and the cells they took, full and surface.
struct Sums {
    sums: Vec<f64>,
    cells: [usize; 2],
}

/// `gs`' integrals over shape < 0 in `lo..hi` (as `integrate`'s), from
/// cells of a 100th of the region halved until the first moved less than
/// `a / 2` of it, at cells some of which are proven full of material
/// (cells none of which fits inside it can misread a wall alike at every
/// size); after the first halving, one that would pass PASS_CELLS is not
/// made, and it is unsettled. `known`: each integral's bounds, proven
/// before. The estimates (clamped into them); `sampled`, twice the last
/// move relative to it, and `moved`, twice each integral's.
#[allow(clippy::too_many_arguments)]
pub(super) fn estimate(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, a: f64, known: (Vec<f64>, Vec<f64>)) -> Measured {
    let mut h = max_extent(lo, hi) / 100.0;
    let mut last = sums(shape, gs, moments, ctx, lo, hi, h);
    let (mut spread, mut settled) = (None, last.sums.is_empty());
    let mut moved = vec![0.0; last.sums.len()];
    while !settled && last.sums[0].is_finite() && !(spread.is_some() && (last.cells[0] + last.cells[1]).saturating_mul(4) > PASS_CELLS) {
        let next = sums(shape, gs, moments, ctx, lo, hi, h / 2.0);
        h /= 2.0;
        let s = 2.0 * (next.sums[0] - last.sums[0]).abs() / next.sums[0].abs();
        settled = next.cells == [0; 2] || (next.cells[0] > 0 && s <= a);
        moved = next.sums.iter().zip(&last.sums).map(|(a, b)| 2.0 * (a - b).abs()).collect();
        (last, spread) = (next, Some(if s.is_finite() { s } else { 0.0 }));
    }
    let (lo, hi) = known;
    for (k, x) in last.sums.iter_mut().enumerate() {
        *x = x.clamp(lo[k], hi[k]);
    }
    let why = (!settled).then(|| match last.cells[0] {
        0 => "no cell is proven full of material (a wall thinner than every cell), and the next halving would take over 2097152 cells".to_string(),
        _ => format!("a halving moved it over {}% of it, and the next would take over 2097152 cells", a * 50.0),
    });
    Measured { sums: last.sums, lo, hi, res: h, boxes: last.cells[0] + last.cells[1], settled, why, sampled: spread, moved }
}

/// A cell's volume, as clipped (an estimate's; bounds use box_volume_out).
fn cell_volume(grid: &Grid, c: Cell) -> f64 {
    grid.cell_sides(c).iter().product()
}

/// `k^3` sample points at the centres of the sub-cells of a cell (as
/// clipped), `k` along each side.
pub fn samples(grid: &Grid, c: Cell, k: usize) -> Vec<V3> {
    let (lo, _) = grid.cell_bounds(c);
    let s = grid.cell_sides(c).map(|x| x / k as f64);
    let mut out = Vec::with_capacity(k * k * k);
    for a in 0..k {
        for b in 0..k {
            for d in 0..k {
                out.push([lo[0] + s[0] * (d as f64 + 0.5), lo[1] + s[1] * (b as f64 + 0.5), lo[2] + s[2] * (a as f64 + 0.5)]);
            }
        }
    }
    out
}

/// Bounds of the integral of `g` over a cell full of material (of every
/// moment and param value of ctx's spans), rounded outward, and whether
/// g's slopes are the same across it to 1e-9 of their size (their range's
/// width: g linear but for rounding, so splits would narrow little): the
/// tighter of its range times the cell's volume and the centred form. With
/// c a point of the cell and a_k(x) the mean of dg/dx_k on the segment from
/// c to x (in the range of that slope over the cell, which covers both
/// sides of a kink), g(x) = g(c) + sum_k a_k(x) (x_k - c_k), so the
/// integral is vol g(c) plus, per axis, at most the slope's range times
/// the integral of x_k - c_k where positive, less where negative: exact
/// for g linear in x, y, z (a centre of mass's moment). A g that jumps in
/// the cell (a select, a floor) has unbounded slopes there and keeps the
/// plain bound. Halved while the slopes vary over it, down to `leaf` grid
/// steps.
fn full_integral(g: &Tape, ctx: &Ctx, grid: &Grid, c: Cell, leaf: i64) -> (I, bool) {
    let b = grid.cell_box(c);
    let side: [I; 3] = b.map(|x| sub_o(I::point(x.hi), I::point(x.lo)));
    let vol = mul_o(mul_o(side[0], side[1]), side[2]);
    let plain = mul_o(vol, eval_interval(g, ctx, b, None));
    let gg = ink_core::eval::eval_interval_grad(g, ctx, b);
    let centre: V3 = b.map(|x| span_mid((x.lo, x.hi)));
    let mut centred = mul_o(vol, eval_interval(g, ctx, centre.map(I::point), None));
    for k in 0..3 {
        // The integral of x_k - c_k where positive and where negative: the
        // cross-section times half the squared distance to each face.
        let area = mul_o(side[(k + 1) % 3], side[(k + 2) % 3]);
        let half = |d: I| mul_o(mul_o(area, mul_o(d, d)), I::point(0.5));
        let up = half(sub_o(I::point(b[k].hi), I::point(centre[k])));
        let down = half(sub_o(I::point(centre[k]), I::point(b[k].lo)));
        let term = sub_o(mul_o(gg.d[k], up), mul_o(gg.d[k], down));
        centred = add_o(centred, term);
    }
    let both = I::new(plain.lo.max(centred.lo), plain.hi.min(centred.hi));
    let linear = |d: &I| d.lo.is_finite() && d.hi.is_finite() && d.hi - d.lo <= 1e-9 * ink_core::interval::mag(*d);
    let curved = !gg.d.iter().all(linear);
    if !curved || c.size <= leaf {
        return (both, !curved);
    }
    let halves = c.children().map(|k| if grid.meets(k) { full_integral(g, ctx, grid, k, leaf).0 } else { I::point(0.0) });
    let sum = halves.iter().fold(I::point(0.0), |s, x| add_o(s, *x));
    (I::new(sum.lo.max(both.lo), sum.hi.min(both.hi)), false)
}

/// Integrals of each of `gs` (None: 1) over the region shape < 0 at cells
/// of `h`; with `moments`, the first's times [1, x, y, z, xx, yy, zz, xy,
/// yz, zx] in its place (mass moments, rho the first). Cells full of
/// material use their 2x2x2 Gauss points, exact where an integrand is
/// linear across them, halved where the first is not and its integral
/// there is not bounded (`full_integral`) to ACCURACY of itself, down to
/// cells of one grid step: points alone miss a layer of density between
/// them at every cell size; bounds see it. Cells the surface may cross use
/// samples weighted by the material in them (`occupancy`). Every integrand
/// at the same points, in one descent.
fn sums(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Sums {
    let (grid, root) = root_for(lo, hi, h, 1);
    let n = gs.len() + if moments && !gs.is_empty() { 9 } else { 0 };
    // Each integral at points weighted `w`.
    let at = |pts: &[V3], w: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; n];
        let mut k = 0;
        for (i, g) in gs.iter().enumerate() {
            let mut v = Vec::new();
            match g {
                Some(t) => eval_points(t, ctx, pts, &mut v),
                None => v.resize(pts.len(), 1.0),
            }
            let terms: &[fn(V3) -> f64] = if i == 0 && moments { &MOMENTS } else { &MOMENTS[..1] };
            for term in terms {
                out[k] = v.iter().zip(pts).zip(w).map(|((a, p), b)| a * term(*p) * b).sum();
                k += 1;
            }
        }
        out
    };
    let gauss = |c: Cell| -> Vec<f64> {
        let (a, b) = grid.cell_bounds(c);
        let g = 0.5 / 3f64.sqrt();
        let pts: Vec<V3> = (0..8).map(|k| std::array::from_fn(|ax| 0.5 * (a[ax] + b[ax]) + if (k >> ax) & 1 == 1 { g } else { -g } * (b[ax] - a[ax]))).collect();
        at(&pts, &[(grid.h * c.size as f64).powi(3) / 8.0; 8])
    };
    fn full(first: Option<&Tape>, ctx: &Ctx, grid: &Grid, c: Cell, gauss: &(dyn Fn(Cell) -> Vec<f64> + Sync)) -> ([usize; 2], Vec<f64>) {
        if let Some(r) = first
            && c.size > 1
        {
            let (b, linear) = full_integral(r, ctx, grid, c, c.size);
            if !linear && !(b.hi - b.lo <= super::ACCURACY * ink_core::interval::mag(b)) {
                let halves: Vec<([usize; 2], Vec<f64>)> = c.children().par_iter().map(|&k| full(first, ctx, grid, k, gauss)).collect();
                return add_sums(halves);
            }
        }
        ([1, 0], gauss(c))
    }
    let first = gs.first().copied().flatten();
    let parts: Vec<([usize; 2], Vec<f64>)> = descend(
        shape,
        ctx,
        &grid,
        root,
        1,
        &|t: &Arc<Tape>, c| {
            let pts = samples(&grid, c, 2);
            let fg = eval_grads(t, ctx, &pts);
            let vol = cell_volume(&grid, c) / pts.len() as f64;
            let s = grid.cell_sides(c).map(|x| x / 2.0);
            let w: Vec<f64> = fg.iter().map(|f| occupancy(f, s) * vol).collect();
            Some(([0, 1], at(&pts, &w)))
        },
        &|_: &Arc<Tape>, c| Some(full(first, ctx, &grid, c, &gauss)),
    );
    let (cells, mut sums) = add_sums(parts);
    sums.resize(n, 0.0);
    Sums { sums, cells }
}

/// `sums`' parts, summed in order.
fn add_sums(parts: Vec<([usize; 2], Vec<f64>)>) -> ([usize; 2], Vec<f64>) {
    parts.into_iter().fold(([0; 2], Vec::new()), |(n, s), (m, x)| {
        let s = if s.is_empty() { x } else { s.iter().zip(&x).map(|(a, b)| a + b).collect() };
        ([n[0] + m[0], n[1] + m[1]], s)
    })
}
