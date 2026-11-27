use super::GridData;
use crate::V3;
use crate::interval::{I, add_up, mul_up};

fn orders(mask: u8) -> [usize; 3] {
    std::array::from_fn(|a| (mask >> a & 1) as usize + 2 * (mask >> (a + 3) & 1) as usize)
}

/// The cubic basis, or its derivative along this axis. A derivative
/// lowers the degree and takes a forward difference of the controls.
/// Apply that difference to the weights so each control is read once.
fn weights(t: f64, order: usize) -> [f64; 4] {
    let s = 1.0 - t;
    match order {
        0 => [s * s * s / 6.0, (4.0 - t * t * (6.0 - 3.0 * t)) / 6.0, (4.0 - s * s * (6.0 - 3.0 * s)) / 6.0, t * t * t / 6.0],
        1 => {
            let (a, b, c) = (0.5 * s * s, 0.5 + s * t, 0.5 * t * t);
            [-a, a - b, b - c, c]
        }
        2 => [s, t - 2.0 * s, s - 2.0 * t, t],
        3 => [-1.0, 3.0, -3.0, 1.0],
        _ => [0.0; 4],
    }
}

const DIFFERENCE: [[f64; 4]; 4] = [[1.0, 0.0, 0.0, 0.0], [-1.0, 1.0, 0.0, 0.0], [1.0, -2.0, 1.0, 0.0], [-1.0, 3.0, -3.0, 1.0]];

impl GridData {
    pub(super) fn spline_component(&self, u: V3, mask: u8) -> f64 {
        let d = orders(mask);
        let h = self.h;
        (0..d.iter().sum::<usize>()).fold(self.spline_local::<false>(u, d).0, |v, _| v / h)
    }

    pub(super) fn spline_grad(&self, u: V3, mask: u8) -> (f64, V3) {
        let d = orders(mask);
        let (mut v, mut g) = self.spline_local::<true>(u, d);
        for _ in 0..d.iter().sum::<usize>() {
            v /= self.h;
            g = g.map(|g| g / self.h);
        }
        (v, g.map(|g| g / self.h))
    }

    fn spline_local<const GRAD: bool>(&self, u: V3, d: [usize; 3]) -> (f64, V3) {
        let c = std::array::from_fn::<_, 3, _>(|a| (u[a] as usize).min(self.dims[a].saturating_sub(2)));
        let t = std::array::from_fn::<_, 3, _>(|a| u[a] - c[a] as f64);
        let [wx, wy, wz] = std::array::from_fn(|a| weights(t[a], d[a]));
        let [dx, dy, dz] = if GRAD { std::array::from_fn(|a| weights(t[a], d[a] + 1)) } else { [[0.0; 4]; 3] };
        let at = std::array::from_fn::<_, 3, _>(|a| std::array::from_fn::<_, 4, _>(|i| (c[a] + i).saturating_sub(1).min(self.dims[a] - 1)));
        let mut sum = 0.0;
        let mut grad = [0.0; 3];
        for (k, &z) in at[2].iter().enumerate() {
            let (mut plane, mut plane_x, mut plane_y) = (0.0, 0.0, 0.0);
            for (j, &y) in at[1].iter().enumerate() {
                let (mut row, mut row_x) = (0.0, 0.0);
                for (i, &x) in at[0].iter().enumerate() {
                    let v = self.at(x, y, z) as f64;
                    row += wx[i] * v;
                    if GRAD {
                        row_x += dx[i] * v;
                    }
                }
                plane += wy[j] * row;
                if GRAD {
                    plane_x += wy[j] * row_x;
                    plane_y += dy[j] * row;
                }
            }
            sum += wz[k] * plane;
            if GRAD {
                grad[0] += wz[k] * plane_x;
                grad[1] += wz[k] * plane_y;
                grad[2] += dz[k] * plane;
            }
        }
        (sum, grad)
    }

    pub(super) fn spline_bound(&self, b: [I; 3], mask: u8) -> I {
        let d = orders(mask);
        let order = d.iter().sum::<usize>();
        let u = std::array::from_fn::<_, 3, _>(|a| {
            let top = (self.dims[a] - 1) as f64;
            let (lo, hi) = self.local_range(b[a], a);
            (lo.min(top), hi.max(0.0))
        });
        let cells = std::array::from_fn(|a| {
            let last = self.dims[a].saturating_sub(2);
            ((u[a].0 as usize).min(last), (u[a].1 as usize).min(last))
        });
        let (lo, hi) = self.sample_bounds(std::array::from_fn(|a| (cells[a].0.saturating_sub(1), (cells[a].1 + 2).min(self.dims[a] - 1))));
        let spread = add_up(hi, -lo);
        // Absolute error for the point code and the differences below.
        // Each cubic weight takes at most seven rounded operations. Its
        // absolute error is under 16 eps times the l1 norm of its stencil;
        // the three weighted sums of four terms take under 8 eps per axis,
        // all measured against the largest absolute control. The l1 norm
        // of the difference stencils' tensor product is 2^order.
        let norm = (1 << order) as f64;
        let err = add_up(mul_up(256.0 * f64::EPSILON * norm, lo.abs().max(hi.abs())), f64::from_bits(1));
        let (mut low, mut high) = if order == 0 {
            (lo, hi)
        } else if (0..3).map(|a| cells[a].1 - cells[a].0 + 4 - d[a]).product::<usize>() <= 128 {
            self.spline_differences(cells, d)
        } else {
            // A nonzero difference sums to zero, with equal positive and
            // negative weight, so only half its norm acts on the spread.
            let mag = mul_up(norm * 0.5, spread);
            (-mag, mag)
        };
        let mid = u.map(|(lo, hi)| lo + 0.5 * (hi - lo));
        let radius = (0..3).fold(0.0, |s, a| add_up(s, add_up(mid[a], -u[a].0).max(add_up(u[a].1, -mid[a]))));
        if radius < 1.0 && !d.contains(&3) {
            // A tight bound as the box shrinks. Along each axis another
            // difference has half its norm times the spread as its bound.
            let v = self.spline_local::<false>(mid, d).0;
            let change = mul_up(radius, mul_up(norm, spread));
            low = low.max(-add_up(change, -v));
            high = high.min(add_up(v, change));
        }
        self.per_h(I::new(low - 2.0 * err, high + 2.0 * err), order as u32)
    }

    fn spline_differences(&self, cells: [(usize, usize); 3], d: [usize; 3]) -> (f64, f64) {
        let at = |i: usize, a: usize| i.saturating_sub(1).min(self.dims[a] - 1);
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        // Indices here are shifted up by one so the left ghost control is 0.
        for k in cells[2].0..=cells[2].1 + 3 - d[2] {
            for j in cells[1].0..=cells[1].1 + 3 - d[1] {
                for i in cells[0].0..=cells[0].1 + 3 - d[0] {
                    let mut v = 0.0;
                    for z in 0..=d[2] {
                        for y in 0..=d[1] {
                            for x in 0..=d[0] {
                                let w = DIFFERENCE[d[0]][x] * DIFFERENCE[d[1]][y] * DIFFERENCE[d[2]][z];
                                v += w * self.at(at(i + x, 0), at(j + y, 1), at(k + z, 2)) as f64;
                            }
                        }
                    }
                    lo = lo.min(v);
                    hi = hi.max(v);
                }
            }
        }
        (lo, hi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        eval::{Ctx, eval_point},
        graph::{Graph, Var},
        tape::Tape,
    };
    use std::sync::Arc;

    #[test]
    fn derivatives_join_at_a_knot() {
        let mut g = Graph::new();
        let p = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let values = vec![0.0, 0.0, 1.0, 0.0, 0.0];
        let grid = GridData::new([0.0; 3], 2.0, [5, 1, 1], values).cubic_bspline().holding();
        let f = g.grid(Arc::new(grid), p);
        let dx = g.diff(f, Var::X);
        let dxx = g.diff(dx, Var::X);
        let ctx = Ctx::new(0.0, vec![]);
        // Spike's right knot, one grid spacing from its centre. The two
        // cubic pieces agree through order two, even though their cubic
        // coefficients differ.
        for (f, want) in [(f, 1.0 / 6.0), (dx, -0.25), (dxx, 0.25)] {
            let tape = Tape::new(&g, f);
            for x in [6.0 - 1e-7, 6.0, 6.0 + 1e-7] {
                let got = eval_point(&tape, &ctx, [x, 0.0, 0.0]);
                assert!((got - want).abs() < 1e-7, "{x}: {got} vs {want}");
            }
        }
    }

    #[test]
    fn bounds_cover_cubic_components_on_both_sides_of_knots() {
        let dims = [5, 4, 6];
        let grid = GridData::new([0.0; 3], 0.5, dims, (0..dims.iter().product()).map(|i| ((i * 17) % 37) as f32 - 18.0).collect()).cubic_bspline();
        for b in [
            [I::new(0.51, 0.55), I::new(0.95, 0.97), I::new(1.37, 1.42)],
            [I::new(0.91, 1.01), I::new(0.95, 1.07), I::new(0.91, 1.08)],
            [I::new(-0.01, 0.01), I::new(1.49, 1.51), I::new(2.5, 3.0)],
            [I::new(-2.0, 4.0); 3],
        ] {
            for mask in 0..64 {
                let r = grid.component_bound(b, mask);
                for k in 0..27 {
                    let p = std::array::from_fn(|a| match k / 3_usize.pow(a as u32) % 3 {
                        0 => b[a].lo,
                        1 => (b[a].lo + b[a].hi) / 2.0,
                        _ => b[a].hi,
                    });
                    let v = grid.component(p, mask);
                    assert!(r.lo <= v && v <= r.hi, "mask {mask} at {p:?}: {v} outside {r:?}");
                    let (value, grad) = grid.component_grad(p, mask);
                    assert_eq!(value, v);
                    for (a, deriv) in grad.into_iter().enumerate() {
                        assert_eq!(deriv, grid.derivative(mask, a).map_or(0.0, |m| grid.component(p, m)));
                    }
                }
            }
        }
        let p = [0.75, 0.5, 1.0];
        let r = grid.component_bound(p.map(I::point), 0);
        assert!(r.hi - r.lo < 1e-10, "{r:?}");
    }
}
