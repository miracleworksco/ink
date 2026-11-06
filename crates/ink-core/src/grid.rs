//! Sampled fields: values on a regular grid, read by trilinear interpolation.
//!
//! Imported meshes (as signed distance), baked expensive fields, and
//! package results (study-field) become fields like any other. This module
//! is the interpolant alone: its components at a point clamped into the
//! grid, and bounds on them over a box (exact within a few cells, else the
//! min/max of the grid values the box covers, through a min/max pyramid).
//! Outside the grid a field continues as the boundary value plus the
//! distance to the grid box (a shape's field), or as the boundary value
//! alone where its maker says it holds (a quantity: `holding`); that
//! continuation is ordinary graph ops (`Graph::grid`), so the core rules
//! give its derivatives and ranges.
//!
//! A grid may say what it was sampled from (`Source`: a mesh, a field, a
//! package's result) and how far that may be from the samples (`Band`).
//! Answers read the samples (the interpolant: `component_bound`), and are
//! measured at the grid's spacing.

use crate::V3;
use crate::interval::I;
use std::sync::Arc;

/// Immutable once made: answers on it are cached by the hash of its
/// samples, taken when it is made (`fingerprint`), and its pyramid is made
/// from them. A grid that cannot change cannot go stale:
///
/// ```compile_fail
/// let mut g = ink_core::grid::GridData::new([0.0; 3], 1.0, [2; 3], vec![1.0; 8]);
/// g.values[0] = -1.0;
/// ```
pub struct GridData {
    origin: V3,
    h: f64,
    /// Points per axis.
    dims: [usize; 3],
    /// x-fastest.
    values: Vec<f32>,
    /// pyramid[k]: min/max over blocks of 2^(k+1) cells per axis (block b
    /// covers points b*2^(k+1) ..= (b+1)*2^(k+1)).
    pyramid: Vec<([usize; 3], Vec<(f32, f32)>)>,
    /// What the samples were taken from, if the grid is a sample of
    /// something else (None: the values are the model's own data).
    source: Option<Arc<Source>>,
    /// A content hash of the samples (origin, spacing, dims, values),
    /// taken once: `fingerprint`.
    samples: u128,
    /// Past the grid box it holds its boundary value (`holding`).
    holds: bool,
}

/// The position of point (or cell) `q` in an x-fastest array of a grid
/// `n` per axis. Every grid in ink (samples, a study's cells and nodes) is
/// laid out so.
pub fn index(q: [usize; 3], n: [usize; 3]) -> usize {
    (q[2] * n[1] + q[1]) * n[0] + q[0]
}

/// The point at position `l` of an x-fastest array of a grid `n` per axis
/// (`index`'s inverse).
pub fn ijk(l: usize, n: [usize; 3]) -> [usize; 3] {
    [l % n[0], (l / n[0]) % n[1], l / (n[0] * n[1])]
}

/// What a sampled field was sampled from, and how far it may be from it.
/// Answers that read the field name it (`Graph::sources`).
#[derive(Clone, Debug, PartialEq)]
pub struct Source {
    /// What was sampled: a file ("ref/base.stl"), or a field or result
    /// ("file.ink:12", "study pull: von-mises").
    pub what: String,
    /// Where it was sampled: the sampling form's file and line.
    pub at: String,
    /// How: "mesh" (import-stl: the signed distance to a triangle mesh),
    /// "bake" (grid-bake), "redistance", "result" (a package's result).
    pub kind: &'static str,
    /// The grid spacing.
    pub res: f64,
    /// A sound bound on how far the source may be from the samples, if one
    /// is known.
    pub band: Option<Band>,
    /// The largest |samples - source| measured at probe points, and how
    /// many points: an estimate, not a bound.
    pub measured: Option<(f64, usize)>,
    /// What the band rests on, or why there is none.
    pub note: Option<String>,
    /// What the samples measure, for clients' legends: its unit ("MPa",
    /// "mm"; "" where none is known), and the least and the most value the
    /// source has where it has one (a result's over its material, not the
    /// values carried past it).
    pub unit: &'static str,
    pub span: Option<(f64, f64)>,
    /// The vector it is a component of (what that is: "study vib:
    /// displacement, mode 2: 823.6 Hz") and the axis it is the component
    /// along (0 x, 1 y, 2 z): a result's displacement-x (clients stating
    /// how far a part warped by it is displaced).
    pub component: Option<(String, usize)>,
}

/// How far a source may be from its samples, soundly.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Band {
    /// Inside the grid box, |source - interpolant| <= err.
    pub err: f64,
    /// Beyond the grid box, if known: the source is 1-Lipschitz and at
    /// least `outside` plus the distance to the box there (a mesh's signed
    /// distance with the mesh `outside` inside the box). Then, with q the
    /// point clamped into the box and e the distance to it, the source is
    /// between outside + e and interpolant(q) + err + e. None: unknown
    /// there.
    pub outside: Option<f64>,
}

impl Source {
    /// The bound a sampled field gets from its spacing h and Lipschitz
    /// constant L alone: an L-Lipschitz d sampled exactly at the grid
    /// points differs from its trilinear interpolant I by at most
    /// L (sqrt 3 / 2) h. At a point p with weights w_c on the corners c of
    /// its cell, I(p) - d(p) = sum w_c (d(c) - d(p)), so |I(p) - d(p)| <=
    /// L sum w_c |c - p| <= L sqrt(sum w_c |c - p|^2) (Jensen) =
    /// L h sqrt(sum_a t_a (1 - t_a)) <= L h sqrt(3 / 4). Rounded up.
    pub fn interpolation_bound(h: f64, lipschitz: f64) -> f64 {
        crate::interval::mul_up(crate::interval::mul_up(h, lipschitz), 0.866_025_403_784_438_7_f64.next_up())
    }
}

impl std::fmt::Debug for GridData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GridData({:?} points, h {}, source {:?})", self.dims, self.h, self.source)
    }
}

impl GridData {
    pub fn new(origin: V3, h: f64, dims: [usize; 3], values: Vec<f32>) -> GridData {
        assert_eq!(values.len(), dims[0] * dims[1] * dims[2]);
        let mut fp = crate::graph::Fingerprint::default();
        fp.f64(h);
        origin.iter().for_each(|v| _ = fp.f64(*v));
        dims.iter().for_each(|v| _ = fp.word(*v as u64));
        fp.word(values.len() as u64);
        values.iter().for_each(|v| _ = fp.raw(&v.to_bits().to_le_bytes()));
        let mut g = GridData { origin, h, dims, values, pyramid: Vec::new(), source: None, samples: fp.value(), holds: false };
        g.build_pyramid();
        g
    }

    /// Hash what answers on the grid can read into `fp` (Tape::fingerprint):
    /// its samples (hashed when it was made; the pyramid is made from
    /// them), and of its source what is a number (the spacing, the error
    /// bound); not the words describing it (what was sampled and where),
    /// which clients add to answers.
    pub fn fingerprint(&self, fp: &mut crate::graph::Fingerprint) {
        fp.key(self.samples).word(self.holds as u64).word(self.source.is_some() as u64);
        if let Some(s) = &self.source {
            let Source { what: _, at: _, kind, res, band, measured, note: _, unit: _, span: _, component: _ } = &**s;
            fp.str(kind).f64(*res).str(&format!("{band:?} {measured:?}"));
        }
    }

    pub fn origin(&self) -> V3 {
        self.origin
    }

    /// The spacing.
    pub fn h(&self) -> f64 {
        self.h
    }

    /// Points per axis.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// The samples, x-fastest.
    pub fn values(&self) -> &[f32] {
        &self.values
    }

    pub fn source(&self) -> Option<&Arc<Source>> {
        self.source.as_ref()
    }

    /// The same samples, holding their boundary value past the grid box:
    /// a quantity (a stress, a displacement) is no distance, and grown by
    /// the distance it is wrong there (a displacement so grown warped
    /// empty space into material). Its maker says so: a package's result
    /// that is a quantity (`ResultSpec::shape` false). A shape's field (an
    /// imported mesh's signed distance, a baked field, a package's solid)
    /// continues as the boundary value plus the distance to the box, so it
    /// stays bounded.
    pub fn holding(mut self) -> GridData {
        self.holds = true;
        self
    }

    /// Whether it holds its boundary value past the grid box (`holding`).
    pub fn holds(&self) -> bool {
        self.holds
    }

    /// The min/max pyramid: level k's blocks are 2^(k+1) cells a side
    /// (block b covers points b 2^(k+1) ..= (b + 1) 2^(k+1)), its count of
    /// blocks per axis and each block's (min, max) of the samples, x-fastest.
    /// The last level has at most 2 blocks per axis; none where the grid has
    /// at most 2 points per axis.
    pub fn pyramid(&self) -> &[([usize; 3], Vec<(f32, f32)>)] {
        &self.pyramid
    }

    /// The same samples, saying what they were taken from.
    pub fn with_source(mut self, source: Source) -> GridData {
        self.source = Some(Arc::new(source));
        self
    }

    fn at(&self, i: usize, j: usize, k: usize) -> f32 {
        self.values[index([i, j, k], self.dims)]
    }

    fn build_pyramid(&mut self) {
        let mut size = 2usize;
        while (0..3).any(|a| self.dims[a] > size) {
            let bd: [usize; 3] = std::array::from_fn(|a| (self.dims[a] - 1).div_ceil(size).max(1));
            let mut blocks = vec![(f32::INFINITY, f32::NEG_INFINITY); bd[0] * bd[1] * bd[2]];
            for bk in 0..bd[2] {
                for bj in 0..bd[1] {
                    for bi in 0..bd[0] {
                        let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
                        for k in bk * size..=((bk + 1) * size).min(self.dims[2] - 1) {
                            for j in bj * size..=((bj + 1) * size).min(self.dims[1] - 1) {
                                for i in bi * size..=((bi + 1) * size).min(self.dims[0] - 1) {
                                    let v = self.at(i, j, k);
                                    lo = lo.min(v);
                                    hi = hi.max(v);
                                }
                            }
                        }
                        blocks[index([bi, bj, bk], bd)] = (lo, hi);
                    }
                }
            }
            self.pyramid.push((bd, blocks));
            size *= 2;
        }
    }

    /// The local coordinate u = (v - origin) / h of a point along `axis`,
    /// clamped into the grid: the one `component` computes.
    fn local(&self, v: f64, axis: usize) -> f64 {
        ((v - self.origin[axis]) / self.h).clamp(0.0, (self.dims[axis] - 1) as f64)
    }

    /// Bounds on the local coordinate u = (v - origin) / h of every point
    /// of `x` along `axis`, clamped into the grid: the formula at the ends
    /// (it is monotone), each step rounded outward (an ulp after each
    /// correctly rounded subtraction and division). So they hold for the
    /// real coordinates and for the rounded u `component` computes.
    fn local_range(&self, x: I, axis: usize) -> (f64, f64) {
        let top = (self.dims[axis] - 1) as f64;
        let o = self.origin[axis];
        (((x.lo - o).next_down() / self.h).next_down().max(0.0), ((x.hi - o).next_up() / self.h).next_up().min(top))
    }

    /// Whether some point of `x` may lie on or beyond either end of the
    /// grid along `axis`, where the coordinate clamps and the field stops
    /// changing along it.
    pub fn may_clamp(&self, x: I, axis: usize) -> bool {
        let (u0, u1) = self.local_range(x, axis);
        u0 <= 0.0 || u1 >= (self.dims[axis] - 1) as f64
    }

    /// The grid box's corners.
    pub fn corners(&self) -> (V3, V3) {
        (self.origin, std::array::from_fn(|a| self.origin[a] + self.h * (self.dims[a] - 1) as f64))
    }

    /// Component `mask` of the trilinear interpolant at `p` (clamped into
    /// the grid): the value for mask 0, else its derivative along the axes
    /// in the mask (bit 0 x, 1 y, 2 z). These 8 are all the derivatives a
    /// trilinear interpolant has (any axis twice gives 0 inside a cell), so
    /// the family is closed under differentiation. A point on a cell face
    /// belongs to the cell above it (the last cell is closed).
    pub fn component(&self, p: V3, mask: u8) -> f64 {
        let mut idx = [0usize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            let u = self.local(p[a], a);
            let i = (u.floor() as usize).min(self.dims[a].saturating_sub(2));
            idx[a] = i;
            t[a] = u - i as f64;
        }
        let v = |di: usize, dj: usize, dk: usize| {
            let i = (idx[0] + di).min(self.dims[0] - 1);
            let j = (idx[1] + dj).min(self.dims[1] - 1);
            let k = (idx[2] + dk).min(self.dims[2] - 1);
            self.at(i, j, k) as f64
        };
        // Reduce the cell's corners along x, then y, then z: a lerp
        // a + (b - a) t, or the difference b - a along a differentiated axis
        // (the same forms and order as the GPU).
        let mut c: Vec<f64> = (0..8).map(|n| v(n & 1, (n >> 1) & 1, n >> 2)).collect();
        for a in 0..3 {
            c = c.chunks(2).map(|ab| if mask >> a & 1 == 1 { ab[1] - ab[0] } else { ab[0] + (ab[1] - ab[0]) * t[a] }).collect();
        }
        // Divided by h once per axis of the mask, as `per_h` bounds it (and
        // the GPU computes it): h^k itself can underflow to 0 or overflow
        // where the quotient does not, which would make a finite
        // derivative infinite or 0, and 0 / 0 for zero data.
        (0..mask.count_ones()).fold(c[0], |v, _| v / self.h)
    }

    /// Component `mask` and its gradient: the components with one more axis
    /// (0 along an axis the mask already has, inside a cell).
    pub fn component_grad(&self, p: V3, mask: u8) -> (f64, V3) {
        let g = std::array::from_fn(|a| if mask >> a & 1 == 1 { 0.0 } else { self.component(p, mask | 1 << a) });
        (self.component(p, mask), g)
    }

    /// Sound bounds of component `mask` over a box (clamped into the grid).
    /// Mask 0 is `value_bound`. Otherwise, in a cell the component interpolates
    /// (linearly along the other axes) the cell's mixed differences along
    /// the mask's axes, so it lies between their min and max over the cells
    /// the box covers.
    pub fn component_bound(&self, b: [I; 3], mask: u8) -> I {
        if b.iter().any(|x| !(x.lo.is_finite() && x.hi.is_finite())) {
            return I::ENTIRE;
        }
        // Small boxes: exact. In each cell the component is multilinear in
        // the axes it does not differentiate (constant along the others),
        // so over the part of the box in that cell its extremes are at that
        // part's corners. Without this a box inside one cell got the range of
        // the whole cell (no proof within a cell of a sampled surface).
        if let Some(r) = self.corner_bound(b, mask) {
            return r;
        }
        if mask == 0 {
            return self.value_bound(b);
        }
        // Per axis: cells along mask axes (indices of their low corner),
        // points along the others.
        let range: [(usize, usize); 3] = std::array::from_fn(|a| {
            let n = self.dims[a];
            let (u0, u1) = self.local_range(b[a], a);
            if mask >> a & 1 == 1 {
                let last = n.saturating_sub(2);
                ((u0.floor().max(0.0) as usize).min(last), (u1.floor().max(0.0) as usize).min(last))
            } else {
                ((u0.floor().max(0.0) as usize).min(n - 1), (u1.ceil() as usize).min(n - 1))
            }
        });
        let count: usize = (0..3).map(|a| range[a].1 - range[a].0 + 1).product();
        if count > 1 << 20 {
            return I::ENTIRE;
        }
        let step: [usize; 3] = std::array::from_fn(|a| (mask >> a & 1) as usize);
        let (mut lo, mut hiv) = (f64::INFINITY, f64::NEG_INFINITY);
        // Rounding of the differences: up to a few ulps of the largest
        // corner involved (cancellation can lose small terms entirely).
        let mut big = 0.0f64;
        for k in range[2].0..=range[2].1 {
            for j in range[1].0..=range[1].1 {
                for i in range[0].0..=range[0].1 {
                    // The mixed difference over the corners i..i+step.
                    let mut d = 0.0;
                    for n in 0..8usize {
                        let o = [n & 1, (n >> 1) & 1, n >> 2];
                        if (0..3).any(|a| o[a] > step[a]) {
                            continue;
                        }
                        let sign = if (o[0] + o[1] + o[2]) % 2 == (step[0] + step[1] + step[2]) % 2 { 1.0 } else { -1.0 };
                        let c = self.at((i + o[0]).min(self.dims[0] - 1), (j + o[1]).min(self.dims[1] - 1), (k + o[2]).min(self.dims[2] - 1)) as f64;
                        big = big.max(c.abs());
                        d += sign * c;
                    }
                    lo = lo.min(d);
                    hiv = hiv.max(d);
                }
            }
        }
        // The samples are exact data (f32), and the f64 sums of up to 8 of
        // them round. With B the largest, the i-th partial sum is at most
        // i B (plus the errors so far, which are smaller), and rounds by at
        // most half an ulp of its binade: u 2^floor(log2(i B)), u = 2^-53.
        // Over i = 2..8 these come to at most 28.5 u B (the worst B is 4/3
        // of a power of 2), under the 32 u B allowed; the subtractions of e
        // round within an ulp of their result, which per_h allows for.
        let e = 16.0 * f64::EPSILON * big + f64::from_bits(1);
        self.per_h(I::new(lo - e, hiv + e), mask)
    }

    /// `component_bound` for a box over at most a few hundred cells: the
    /// extremes over each cell's part of the box, at its corners. None for
    /// larger boxes.
    fn corner_bound(&self, b: [I; 3], mask: u8) -> Option<I> {
        let n = self.dims;
        // Cells (by low index) the clamped box touches, and the box in local
        // coordinates u = (p - origin) / h.
        let mut cells = [(0usize, 0usize); 3];
        let mut u = [(0.0f64, 0.0f64); 3];
        for a in 0..3 {
            let last = n[a].saturating_sub(2);
            let (u0, u1) = self.local_range(b[a], a);
            u[a] = (u0, u1);
            cells[a] = ((u0.floor().max(0.0) as usize).min(last), (u1.floor().max(0.0) as usize).min(last));
        }
        if (0..3).map(|a| cells[a].1 - cells[a].0 + 1).product::<usize>() > 512 {
            return None;
        }
        let (mut lo, mut hiv, mut big) = (f64::INFINITY, f64::NEG_INFINITY, 0.0f64);
        for k in cells[2].0..=cells[2].1 {
            for j in cells[1].0..=cells[1].1 {
                for i in cells[0].0..=cells[0].1 {
                    let c = [i, j, k];
                    let corner = |di: usize, dj: usize, dk: usize| self.at((i + di).min(n[0] - 1), (j + dj).min(n[1] - 1), (k + dk).min(n[2] - 1)) as f64;
                    let vals: Vec<f64> = (0..8).map(|m| corner(m & 1, (m >> 1) & 1, m >> 2)).collect();
                    big = vals.iter().fold(big, |acc, v| acc.max(v.abs()));
                    // The box's part in this cell, in the cell's coordinates.
                    let t: [(f64, f64); 3] = std::array::from_fn(|a| ((u[a].0 - c[a] as f64).clamp(0.0, 1.0), (u[a].1 - c[a] as f64).clamp(0.0, 1.0)));
                    for m in 0..8usize {
                        let p: [f64; 3] = std::array::from_fn(|a| if m >> a & 1 == 1 { t[a].1 } else { t[a].0 });
                        // The cell's interpolant (its own, also on its upper
                        // faces), reduced as in `component`.
                        let mut r = vals.clone();
                        for a in 0..3 {
                            r = r.chunks(2).map(|ab| if mask >> a & 1 == 1 { ab[1] - ab[0] } else { ab[0] + (ab[1] - ab[0]) * p[a] }).collect();
                        }
                        lo = lo.min(r[0]);
                        hiv = hiv.max(r[0]);
                    }
                }
            }
        }
        // Rounding of the reductions: a few ulps of the largest corner. With
        // B the largest and u = 2^-53, a difference of values within E of
        // exact and at most M carries 2 E and rounds by 2 u M (magnitude
        // 2 M); a lerp a + (b - a) t, weights 1 - t and t, carries E and
        // rounds by 5 u M (the difference, the product, the sum). Over the
        // three axes, from E = 0 and M = B, at most 28 u B whatever the
        // mask (difference, difference, lerp: 2, 8, 28), under the 64 u B
        // allowed. Products that are subnormal round by half the least
        // subnormal, far below 64 u times any nonzero f32.
        let e = 32.0 * f64::EPSILON * big + f64::from_bits(1);
        Some(self.per_h(I::new(lo - e, hiv + e), mask))
    }

    /// Differences `r` (rounded from the ends by up to 2 ulps) divided by h
    /// once per axis of the mask, outward: h^k itself can underflow to 0 or
    /// overflow where the quotient does not, and an end that overflows stays
    /// a finite bound on its side.
    fn per_h(&self, r: I, mask: u8) -> I {
        let mut r = crate::interval::widen(r, 2);
        for _ in 0..mask.count_ones() {
            r = crate::interval::widen(I::new(r.lo / self.h, r.hi / self.h), 1);
        }
        r
    }

    /// Whether the box may cross a cell face along `axis` (where components
    /// differentiated along it jump). The ends' local coordinates are
    /// rounded outward (local_range): an end whose rounding could reach a
    /// face counts as crossing it, since a coordinate a hair below the face
    /// can round onto it.
    pub fn crosses_face(&self, b: I, axis: usize) -> bool {
        if !(b.lo.is_finite() && b.hi.is_finite()) {
            return true;
        }
        let (u0, u1) = self.local_range(b, axis);
        let cell = |u: f64| (u.floor() as usize).min(self.dims[axis].saturating_sub(2));
        cell(u0) != cell(u1)
    }

    /// Bounds of the value (component 0) over a box clamped into the grid:
    /// the min/max of the grid points it touches, each value a convex
    /// combination of them. Like `component`, no exterior distance: the
    /// graph adds that (`Graph::grid`).
    fn value_bound(&self, b: [I; 3]) -> I {
        // Grid points the clamped box touches (cells' corners included).
        let range: [(usize, usize); 3] = std::array::from_fn(|a| {
            let (u0, u1) = self.local_range(b[a], a);
            ((u0.floor() as usize).min(self.dims[a] - 1), (u1.ceil() as usize).min(self.dims[a] - 1))
        });
        let span = (0..3).map(|a| range[a].1 - range[a].0 + 1).max().unwrap();
        let (mut lo, mut hiv) = (f32::INFINITY, f32::NEG_INFINITY);
        // Coarsest detail needed: the finest pyramid level whose blocks cover
        // the range in at most ~8 per axis (or single points for tiny ranges).
        if span > 4 && !self.pyramid.is_empty() {
            let mut k = 0;
            while k + 1 < self.pyramid.len() && span / (2usize << k) > 8 {
                k += 1;
            }
            let size = 2usize << k;
            let (bd, blocks) = &self.pyramid[k];
            let br: [(usize, usize); 3] = std::array::from_fn(|a| ((range[a].0 / size).min(bd[a] - 1), (range[a].1 / size).min(bd[a] - 1)));
            for bk in br[2].0..=br[2].1 {
                for bj in br[1].0..=br[1].1 {
                    for bi in br[0].0..=br[0].1 {
                        let (a, c) = blocks[index([bi, bj, bk], *bd)];
                        lo = lo.min(a);
                        hiv = hiv.max(c);
                    }
                }
            }
        } else {
            for k in range[2].0..=range[2].1 {
                for j in range[1].0..=range[1].1 {
                    for i in range[0].0..=range[0].1 {
                        let v = self.at(i, j, k);
                        lo = lo.min(v);
                        hiv = hiv.max(v);
                    }
                }
            }
        }
        // The samples are exact (f32); a value is three lerps of them in f64,
        // which round by at most 28 u B (B the largest, u = 2^-53:
        // `corner_bound`'s reckoning with no axis differentiated), within
        // the 64 u B allowed; the subtractions of e round within an ulp.
        let big = lo.abs().max(hiv.abs()) as f64;
        let e = 32.0 * f64::EPSILON * big + f64::from_bits(1);
        crate::interval::widen(I::new(lo as f64 - e, hiv as f64 + e), 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::{Ctx, eval_dual, eval_grads, eval_interval, eval_interval_grad, eval_point};
    use crate::graph::{Graph, Var};
    use crate::tape::Tape;
    use std::sync::Arc;

    // The specification, written independently of `component`: each
    // component is the sum over the cell's corners of tensor-product weights
    // (1 - t, t), with -1/h, 1/h along a differentiated axis. The point is
    // clamped into the grid; a point on a cell face belongs to the cell
    // above it, the last cell is closed.
    fn spec(g: &GridData, p: V3, mask: u8) -> f64 {
        let mut cell = [0usize; 3];
        let mut t = [0.0; 3];
        for a in 0..3 {
            let top = (g.dims[a] - 1) as f64;
            let u = ((p[a] - g.origin[a]) / g.h).clamp(0.0, top);
            cell[a] = (u.floor() as usize).min(g.dims[a].saturating_sub(2));
            t[a] = u - cell[a] as f64;
        }
        let mut sum = 0.0;
        for n in 0..8usize {
            let o = [n & 1, (n >> 1) & 1, n >> 2];
            let w: f64 = (0..3)
                .map(|a| match (mask >> a & 1 == 1, o[a]) {
                    (true, 0) => -1.0 / g.h,
                    (true, _) => 1.0 / g.h,
                    (false, 0) => 1.0 - t[a],
                    (false, _) => t[a],
                })
                .product();
            let i: [usize; 3] = std::array::from_fn(|a| (cell[a] + o[a]).min(g.dims[a] - 1));
            sum += w * g.at(i[0], i[1], i[2]) as f64;
        }
        sum
    }

    // The continued field the graph specifies (`Graph::grid`): the
    // interpolant at the clamped point plus the distance to the grid box.
    // Its gradient: along an axis strictly inside, the interpolant's slope;
    // outside, the distance's; on a grid-box face, half the interpolant's
    // slope (the clamp's min/max tie) and none of the distance's.
    fn spec_field(g: &GridData, p: V3) -> (f64, V3) {
        let (lo, hi) = g.corners();
        let q: V3 = std::array::from_fn(|a| p[a].clamp(lo[a], hi[a]));
        let dist = (0..3).map(|a| (p[a] - q[a]).powi(2)).sum::<f64>().sqrt();
        let grad = std::array::from_fn(|a| {
            if p[a] < lo[a] || p[a] > hi[a] {
                (p[a] - q[a]) / dist
            } else if p[a] == lo[a] || p[a] == hi[a] {
                0.5 * spec(g, q, 1 << a)
            } else {
                spec(g, q, 1 << a)
            }
        });
        (spec(g, q, 0) + dist, grad)
    }

    fn grid() -> GridData {
        let dims = [4, 3, 3];
        let values = (0..36)
            .map(|n| {
                let (i, j, k) = (n % 4, (n / 4) % 3, n / 12);
                ((i * 7 + j * 13 + k * 5) % 11) as f32 * 0.375 - 2.0 + (i * j * k) as f32 * 0.25
            })
            .collect();
        GridData::new([-1.0, 0.5, 2.0], 0.5, dims, values)
    }

    // Local coordinates per axis: outside below, the grid-box faces, inside a
    // cell, every interior cell face and a hair to each side, outside above.
    // Their products cover exterior faces, edges and corners.
    fn coords(n: usize) -> Vec<f64> {
        let top = (n - 1) as f64;
        let hair = 2f64.powi(-20);
        let mut u = vec![-1.25, 0.0, 0.3, top - 0.4, top, top + 0.75];
        for k in 1..n - 1 {
            u.extend([k as f64 - hair, k as f64, k as f64 + hair]);
        }
        u
    }

    fn points(g: &GridData) -> Vec<V3> {
        let u: [Vec<f64>; 3] = std::array::from_fn(|a| coords(g.dims[a]));
        let mut pts = Vec::new();
        for &uz in &u[2] {
            for &uy in &u[1] {
                for &ux in &u[0] {
                    let uu = [ux, uy, uz];
                    pts.push(std::array::from_fn(|a| g.origin[a] + uu[a] * g.h));
                }
            }
        }
        pts
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + b.abs())
    }

    fn field_tape(data: &Arc<GridData>) -> (Graph, u32) {
        let mut g = Graph::new();
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let f = g.grid(data.clone(), xyz);
        (g, f)
    }

    const CTX: Ctx = Ctx::new(0.0, Vec::new());

    #[test]
    fn components_match_the_specification_everywhere() {
        // Every mask (mixed ones included), at cell faces and either side,
        // on and beyond the grid box: `component` is the specified
        // interpolant, and `component_grad` its next components.
        let g = grid();
        for p in points(&g) {
            for mask in 0..8u8 {
                let (v, d) = g.component_grad(p, mask);
                assert!(close(v, spec(&g, p, mask)), "mask {mask} at {p:?}: {v} vs {}", spec(&g, p, mask));
                for a in 0..3 {
                    let want = if mask >> a & 1 == 1 { 0.0 } else { spec(&g, p, mask | 1 << a) };
                    assert!(close(d[a], want), "mask {mask} d{a} at {p:?}: {} vs {want}", d[a]);
                }
            }
        }
    }

    #[test]
    fn the_graph_continues_the_field_as_specified() {
        // Value and gradient of `Graph::grid` from every point evaluator and
        // the symbolic derivative, inside, on the grid-box faces and outside
        // across its faces, edges and corners.
        let data = Arc::new(grid());
        let (mut g, f) = field_tape(&data);
        let d = [Var::X, Var::Y, Var::Z].map(|v| g.diff(f, v));
        let dxy = g.diff(d[0], Var::Y);
        let tape = Tape::new(&g, f);
        let dt = d.map(|n| Tape::new(&g, n));
        let dxy = Tape::new(&g, dxy);
        let pts = points(&data);
        let grads = eval_grads(&tape, &CTX, &pts);
        for (p, gr) in pts.iter().zip(&grads) {
            let (v, want) = spec_field(&data, *p);
            let dual = eval_dual(&tape, &CTX, *p);
            assert!(close(eval_point(&tape, &CTX, *p), v) && close(gr[0], v) && close(dual[0], v), "value at {p:?}");
            for a in 0..3 {
                let sym = eval_point(&dt[a], &CTX, *p);
                assert!(close(gr[1 + a], want[a]), "grads d{a} at {p:?}: {} vs {}", gr[1 + a], want[a]);
                assert!(close(dual[1 + a], want[a]), "dual d{a} at {p:?}: {} vs {}", dual[1 + a], want[a]);
                assert!(close(sym, want[a]), "symbolic d{a} at {p:?}: {sym} vs {}", want[a]);
            }
            // The mixed derivative strictly inside the grid box: the
            // interpolant's, on cell faces too (the cell above).
            let (lo, hi) = data.corners();
            if (0..3).all(|a| lo[a] < p[a] && p[a] < hi[a]) {
                assert!(close(eval_point(&dxy, &CTX, *p), spec(&data, *p, 3)), "d2/dxdy at {p:?}");
            }
        }
    }

    #[test]
    fn tiny_boxes_across_faces_bound_both_sides() {
        // A box a hair wide around each point (across a cell face or the grid
        // box face when the point is on one): the raw components' bounds hold
        // each corner's one-sided component, and the graph's value and
        // gradient bounds hold the continued field at the corners and center.
        let data = Arc::new(grid());
        let (g, f) = field_tape(&data);
        let tape = Tape::new(&g, f);
        let delta = data.h * 2f64.powi(-30);
        for p in points(&data) {
            let b = p.map(|v| I::new(v - delta, v + delta));
            let corners: Vec<V3> = (0..8).map(|n| std::array::from_fn(|a| if n >> a & 1 == 1 { b[a].hi } else { b[a].lo })).collect();
            for mask in 0..8u8 {
                let r = data.component_bound(b, mask);
                for c in &corners {
                    assert!(r.contains(spec(&data, *c, mask)), "mask {mask} over {b:?}: {r:?} misses {}", spec(&data, *c, mask));
                }
            }
            let v = eval_interval(&tape, &CTX, b, None);
            let gb = eval_interval_grad(&tape, &CTX, b);
            for c in corners.iter().chain([&p]) {
                let (val, grad) = spec_field(&data, *c);
                assert!(v.contains(val) && gb.v.contains(val), "value over {b:?}: {v:?} {:?} miss {val}", gb.v);
                for a in 0..3 {
                    assert!(gb.d[a].contains(grad[a]), "d{a} over {b:?}: {:?} misses {}", gb.d[a], grad[a]);
                }
            }
        }
    }

    #[test]
    fn raw_bounds_count_no_exterior_distance() {
        // A box beyond the grid over 23 x 23 cells (past the corner bound's
        // 512) takes the pyramid path. The raw value bound holds the clamped interpolant (the sample
        // range), and the graph adds the distance once: [5, ...] here, not
        // twice.
        let n = 24;
        let values = (0..n * n * n).map(|i| ((i * 37 % 101) as f32 - 50.0) * 0.1).collect();
        let data = Arc::new(GridData::new([0.0; 3], 1.0, [n; 3], values));
        let (lo_s, hi_s) = data.values.iter().fold((f64::MAX, f64::MIN), |(l, h), v| (l.min(*v as f64), h.max(*v as f64)));
        let top = (n - 1) as f64;
        let b = [I::new(top + 5.0, top + 10.0), I::new(0.0, top), I::new(0.0, top)];
        assert!(data.corner_bound(b, 0).is_none(), "the box must take the pyramid path");
        let r = data.component_bound(b, 0);
        assert!(r.lo <= lo_s && r.hi >= hi_s && r.lo > lo_s - 1e-3 && r.hi < hi_s + 1e-3, "{r:?} vs [{lo_s}, {hi_s}]");
        let (g, f) = field_tape(&data);
        let v = eval_interval(&Tape::new(&g, f), &CTX, b, None);
        assert!(v.lo <= lo_s + 5.0 && v.lo > lo_s + 5.0 - 1e-3 && v.hi >= hi_s + 10.0, "{v:?}");
    }
}
