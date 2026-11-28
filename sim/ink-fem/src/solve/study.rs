//! What a solve fixes before its cells (`Study`): the grid over the part,
//! the samples of a cell, and the material fields' reads over a box.

use crate::cells::{
    AIR, K, MATERIAL, MAX_DEPTH, PROPERTY, PROPERTY_ROUNDING, PROPERTY_SHRINK, PROPERTY_STEEP, PROPERTY_UNSEEN, THIN, pieces, sample_chords, thin_and_resolved, variation_resolved,
};
use crate::element::sample_xi;
use crate::eval::{Ctx, eval_interval, eval_points};
use crate::occupancy::{Classified, Cov, at_xi, occupied_spread, spread_moments};
use crate::octree::{Cell, Grid, root_for};
use crate::surface::{Facet, part_surface, surface_samples};
use crate::tape::Tape;
use crate::{Boundary, Fix, Held, Physics, RegionSel, Setup};
use ink_core::tape::Names;
use ink_core::{V3, eval, grid, interval};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use super::SolveError;

/// A solve's fixed context: its setup, the grid over the part, a cell's
/// samples, and the material fields it reads (each of their reads a method).
pub(super) struct Study<'a> {
    pub(super) s: &'a Setup,
    pub(super) ctx: &'a Ctx,
    pub(super) ndof: usize,
    /// A cell's degrees of freedom: 8 nodes of `ndof`.
    pub(super) ne: usize,
    pub(super) h: f64,
    /// The declared mirrors: (axis, coordinate).
    pub(super) mirrors: Vec<(usize, f64)>,
    pub(super) dims: [i64; 3],
    pub(super) lo: V3,
    pub(super) total_cells: i64,
    pub(super) grid: Grid,
    pub(super) root: Cell,
    /// A cell's samples (`sample_xi(K)`), each `w_q` of the reference cube.
    pub(super) xis: Vec<V3>,
    pub(super) w_q: f64,
    /// A full cell's samples (`sample_xi(2)`).
    pub(super) full_xi: Vec<V3>,
    pub(super) default_nu: f64,
    pub(super) props: Vec<&'a Tape>,
    pub(super) series: bool,
    pub(super) judge: f64,
    pub(super) tiles: usize,
    /// A fix with free motions (a pivot, a slide) holds its region's mean
    /// motion, as a distributing hold (RBE3) does, unless it says (rigid):
    /// its region's material is free material, deforming as it will, and its
    /// boundary's hold (its terms' kww, the weight) holds only the region's
    /// rigid motions its free ones leave (`freeing`). Held rigidly over a
    /// pinned column's end face, the bar could not spread there as it was
    /// squeezed, and that sideways stress, meeting the buckle's slope at the
    /// end that turns, read the column's factor 0.6% under Euler's (with
    /// shear) at nu 0.3; nu 0 read it right.
    pub(super) distributing: Vec<bool>,
    /// The fixes whose held value reads t (a transient's: `timed_holds`).
    pub(super) timed: Vec<bool>,
    /// The moment the held values are read at (`values_at`), as its
    /// bits: a transient's whose holds vary in time sets it to each step's
    /// end, and to its duration for the state it reports.
    held_t: AtomicU64,
}

impl<'a> Study<'a> {
    pub(super) fn new(s: &'a Setup, ctx: &'a Ctx) -> Result<Self, SolveError> {
        let ndof = s.physics.ndof();
        let ne = 8 * ndof;
        let h = s.h;

        let mirrors: Vec<_> = s
            .fixes
            .iter()
            .filter_map(|f| match f.boundary {
                Boundary::Mirror { axis, at } => Some((axis, at)),
                _ => None,
            })
            .collect();
        let mut seen_axes = [false; 3];
        for &(axis, at) in &mirrors {
            if axis >= 3 || !at.is_finite() || std::mem::replace(&mut seen_axes[axis], true) || ndof != 3 {
                return Err(SolveError::InvalidMirror("declare at most one finite mirror coordinate per axis in an elastic study".into()));
            }
            if s.body.as_ref().is_some_and(|(_, g)| g[axis] != 0.0) {
                return Err(SolveError::InvalidMirror("gravity acts across a declared mirror. Model the whole part with the same gravity vector on both sides".into()));
            }
            let mut plane = std::array::from_fn(|a| interval::I::new(s.lo[a], s.hi[a]));
            plane[axis] = interval::I::point(at);
            for f in &s.fixes {
                let u = match &f.value[axis] {
                    Held::Number(v) => interval::I::point(*v),
                    Held::Field(t) => eval_interval(t, ctx, plane, None),
                };
                if f.dofs[axis] && (u.lo != 0.0 || u.hi != 0.0) && f.region.may_touch(ctx, plane) {
                    // An inconclusive interval test cannot discard a legacy interpretation.
                    return Err(
                        format!("cannot establish zero normal displacement at the mirror under fix '{}'. Move the support off the cut or model the whole part", f.name).into()
                    );
                }
            }
        }

        // Grid covering the part, one cell of padding (which also absorbs the
        // division's rounding), its lines at multiples of h from the origin, not
        // from the bounds: those are only as tight as their search's cells, so a
        // grid placed from them sits off the part's faces by that slack, and
        // every cell along a face is cut. On the lattice a face at a multiple of
        // h (a box's, at round coordinates) lies on cell boundaries, and the grid
        // at twice the size (a claim's coarser grid) nests in this one.
        let first: [i64; 3] = std::array::from_fn(|a| (s.lo[a] / h).floor() as i64 - 1);
        let dims: [i64; 3] = std::array::from_fn(|a| (s.hi[a] / h).ceil() as i64 + 1 - first[a]);
        let lo: V3 = std::array::from_fn(|a| first[a] as f64 * h);
        let total_cells = dims.iter().product::<i64>();
        if total_cells > 400_000_000 {
            return Err(format!("grid of {dims:?} cells is too fine; increase res").into());
        }
        let (grid, root) = root_for(lo, [lo[0] + dims[0] as f64 * h, lo[1] + dims[1] as f64 * h, lo[2] + dims[2] as f64 * h], h, 1);
        let xis = sample_xi(K);
        let w_q = 8.0 / xis.len() as f64;

        let default_nu = match s.physics {
            Physics::Elastic { nu } => nu,
            Physics::Heat => 0.0,
        };
        let mass_read = s.modes > 0 || s.transient.is_some();
        let props: Vec<&Tape> =
            [s.coef_const.is_none().then_some(&*s.coef), s.nu.as_deref(), s.strength.as_deref(), s.body.as_ref().map(|b| &*b.0), s.density.as_deref().filter(|_| mass_read)]
                .into_iter()
                .flatten()
                .collect();
        let series = s.coef_const.is_none();
        // The grid cells' size thin features are judged by (`Setup::split_as`).
        let judge = s.split_as.unwrap_or(h);
        // A cell's material is judged as the grid it is split as would judge
        // it: over that grid's cells, as many as tile it, so a coarser grid
        // stops splitting at the level the claim's grid does, its cells there
        // twice or four times as large, and sees what any of them would.
        let tiles = (h / judge).round().max(1.0) as usize;
        let distributing: Vec<bool> =
            s.fixes.iter().map(|f| !f.free.is_empty() && !f.rigid && f.film.is_none() && !matches!(f.boundary, Boundary::Mirror { .. }) && s.physics.ndof() == 3).collect();
        let full_xi = sample_xi(2);
        // A held value that varies in time is held as it is at each step of
        // a transient (`over_time`); with no time to step, it holds nothing.
        let timed: Vec<bool> = s.fixes.iter().map(|f| f.value.iter().any(|v| matches!(v, Held::Field(t) if Names::of(&[&**t], ctx).reads_t))).collect();
        if let Some(f) = s.fixes.iter().zip(&timed).find(|x| *x.1).map(|x| x.0)
            && s.transient.is_none()
        {
            return Err(format!("fix '{}' holds a value that varies in time (t): only a transient study steps through time", f.name).into());
        }
        // (The state a transient reports is its last: its holds' then.)
        let held_t = AtomicU64::new(s.transient.filter(|_| timed.contains(&true)).map_or(ctx.t, |tr| tr.duration).to_bits());
        Ok(Study { s, ctx, ndof, ne, h, mirrors, dims, lo, total_cells, grid, root, xis, w_q, full_xi, default_nu, props, series, judge, tiles, distributing, timed, held_t })
    }

    pub(super) fn in_grid(&self, i: [i64; 3]) -> bool {
        let dims = self.dims;
        (0..3).all(|a| i[a] >= 0 && i[a] < dims[a])
    }

    pub(super) fn lin(&self, i: [i64; 3]) -> usize {
        let dims = self.dims;
        grid::index(i.map(|v| v as usize), dims.map(|v| v as usize))
    }

    /// By grid cell, then a child's place (on the finest lattice: children
    /// of a grid cell may be of several sizes).
    pub(super) fn order(&self, c: &Classified) -> (usize, u32, u32, u32) {
        let at = |a: usize| u32::from(c.sub[a]) << (MAX_DEPTH - c.level);
        (self.lin(c.idx), at(2), at(1), at(0))
    }

    pub(super) fn coef_at(&self, pts: &[V3]) -> Vec<f64> {
        let Study { s, ctx, .. } = *self;
        match s.coef_const {
            Some(c) => vec![c; pts.len()],
            None => {
                let mut v = Vec::new();
                eval_points(&s.coef, ctx, pts, &mut v);
                v
            }
        }
    }

    pub(super) fn nu_at(&self, pts: &[V3]) -> Vec<f64> {
        let Study { s, ctx, default_nu, .. } = *self;
        match &s.nu {
            None => vec![default_nu; pts.len()],
            Some(t) => {
                let mut v = Vec::new();
                eval_points(t, ctx, pts, &mut v);
                v
            }
        }
    }

    /// The material fields the solve reads, and the yield its claims read
    /// (`PROPERTY`): their ranges over a box of side `size` about `ctr`.
    /// (Density only where the solve reads it: a mass, a heat capacity.)
    /// The coefficient a solve averages in series (E, k), the first field
    /// where it varies: a cell's one value stands for a stack of what lies
    /// across it, whose sum is 1/k's (a series resistance, a bar's
    /// stretch), so its reciprocal is judged as a property too. Read as k,
    /// a dip toward 0 is smooth while 1/k is near singular at its bottom
    /// (k 0.01 + 100 |x|: 1/k 0.1 to 100 over a 0.2 mm cell, its one value
    /// 9 to 60% off); and a film of k 0.01 among k 100 to 200 moves k's
    /// range less than a slope beside it can, but 1/k's by 100 against
    /// thousandths. Unbounded where k reaches 0.
    pub(super) fn fields(&self, b: [interval::I; 3]) -> Vec<interval::I> {
        let Study { ctx, ref props, series, .. } = *self;
        let mut r: Vec<interval::I> = props.iter().map(|t| eval_interval(t, ctx, b, None)).collect();
        if series {
            let k = r[0];
            r.push(if k.lo > 0.0 && k.hi.is_finite() { interval::I::new(1.0 / k.hi, 1.0 / k.lo) } else { interval::I::new(f64::NEG_INFINITY, f64::INFINITY) });
        }
        r
    }

    pub(super) fn hulls(&self, ctr: V3, size: f64) -> Vec<interval::I> {
        let half = 0.5 * size * (1.0 - PROPERTY_SHRINK);
        self.fields(std::array::from_fn(|a| interval::I::new(ctr[a] - half, ctr[a] + half)))
    }

    /// Whether the properties' variation over a box (their ranges `r` over
    /// it) is resolved at its size (`variation_resolved`).
    pub(super) fn resolved_box(&self, ctr: V3, size: f64, r: &[interval::I]) -> bool {
        if r.iter().all(|i| width(i) == 0.0) {
            return true;
        }
        let (quarter, half) = (self.hulls(ctr, 0.25 * size), self.hulls(ctr, 0.5 * size));
        (0..r.len()).all(|k| variation_resolved(width(&quarter[k]), width(&half[k]), width(&r[k])))
    }

    /// Each property's slopes over a box of side `size` about `ctr`
    /// (`eval_interval_grad`), where it varies (its range `r`), and k's
    /// where 1/k does.
    pub(super) fn slopes(&self, ctr: V3, size: f64, r: &[interval::I]) -> Vec<Option<interval::G>> {
        let half = 0.5 * size * (1.0 - PROPERTY_SHRINK);
        let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - half, ctr[a] + half));
        let inv = self.series && width(&r[self.props.len()]) != 0.0;
        self.props.iter().enumerate().map(|(k, t)| (width(&r[k]) != 0.0 || k == 0 && inv).then(|| eval::eval_interval_grad(t, self.ctx, b))).collect()
    }

    /// Whether each property (and `1/k` where `series` averages in series) is
    /// no steeper over a box of side `size` (their ranges `r` over it, their
    /// slopes `g`) than its range allows (`PROPERTY_STEEP`): rejects jumps
    /// (`select`, steps) and features steeper than that ratio, including
    /// soft films in series where `1/k` is steep while `k`'s contrast is
    /// small next to `k`'s range.
    pub(super) fn gentle(&self, size: f64, r: &[interval::I], g: &[Option<interval::G>]) -> bool {
        let steep = |g: &Option<interval::G>, over: f64, w: f64| g.as_ref().is_some_and(|g| interval::lipschitz(g) / over * size <= PROPERTY_STEEP * w);
        let n = self.props.len();
        (0..n).all(|k| width(&r[k]) == 0.0 || steep(&g[k], 1.0, width(&r[k])))
            && (!self.series || width(&r[n]) == 0.0 || r[0].lo > 0.0 && r[0].hi.is_finite() && steep(&g[0], r[0].lo * r[0].lo, width(&r[n])))
    }

    /// Whether each property's range over a box of side `size` about `ctr`
    /// (`r`) reaches no further beyond what the centres of the box and of
    /// the 26 of its size round it take (`at(k)`: property k's value at the
    /// box's centre, and its least and greatest there) than PROPERTY_UNSEEN
    /// of that, on each side: a film between them is missed by the cells on
    /// every grid. Where its bounds reach further, the range is taken no
    /// wider than its slopes allow from the centre's value (`within`): a
    /// field written as a polynomial is bounded wider than it is (a
    /// smoothstep's flank, over a 0.25 mm box, 2.4 times its range), which
    /// reached past the centres on a monotone field. Where that still
    /// reaches too far, its eight halves are so taken: a slope bound's reach
    /// is first order in the box's side, and a smooth field's least or
    /// greatest beyond the centres second.
    pub(super) fn seen(&self, ctr: V3, size: f64, r: &[interval::I], at: impl Fn(usize) -> (f64, f64, f64)) -> bool {
        let Study { ctx, ref props, .. } = *self;
        let half = 0.5 * size * (1.0 - PROPERTY_SHRINK);
        let whole: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - half, ctr[a] + half));
        // Each half from the box's centre to its (shrunk) face: together the box.
        let halves: Vec<V3> = (0..8).map(|i| std::array::from_fn(|a| ctr[a] + 0.25 * size * if i >> a & 1 == 1 { 1.0 } else { -1.0 })).collect();
        let half_box =
            |c: V3| -> [interval::I; 3] { std::array::from_fn(|a| if c[a] < ctr[a] { interval::I::new(ctr[a] - half, ctr[a]) } else { interval::I::new(ctr[a], ctr[a] + half) }) };
        let mut v = Vec::new();
        props.iter().zip(r).enumerate().all(|(k, (t, r))| {
            if width(r) == 0.0 {
                return true;
            }
            let (centre, lo, hi) = at(k);
            if !centre.is_finite() || !lo.is_finite() || !hi.is_finite() {
                return false;
            }
            let inside = |b: interval::I| b.hi <= hi + (PROPERTY_UNSEEN + PROPERTY_ROUNDING) * hi.abs() && b.lo >= lo - (PROPERTY_UNSEEN + PROPERTY_ROUNDING) * lo.abs();
            inside(*r) || inside(within(*r, centre, &eval::eval_interval_grad(t, ctx, whole), half)) || {
                eval_points(t, ctx, &halves, &mut v);
                halves.iter().zip(&v).all(|(&c, &at)| {
                    let b = half_box(c);
                    inside(within(eval_interval(t, ctx, b, None), at, &eval::eval_interval_grad(t, ctx, b), 0.25 * size))
                })
            }
        })
    }

    /// `seen` over a box of side `size` about `ctr`, the centres round it
    /// inside the part read here (ignoring neighbours outside the part).
    pub(super) fn seen_round(&self, ctr: V3, size: f64, r: &[interval::I]) -> bool {
        let pts: Vec<V3> = (0..27).map(|i| std::array::from_fn(|a| ctr[a] + ((i / 3usize.pow(a as u32) % 3) as f64 - 1.0) * size)).collect();
        let mut sv = Vec::new();
        eval_points(&self.s.shape, self.ctx, &pts, &mut sv);
        let in_part: [bool; 27] = std::array::from_fn(|j| sv[j] <= 0.0);
        let values: Vec<Vec<f64>> = self
            .props
            .iter()
            .map(|t| {
                let mut v = Vec::new();
                eval_points(t, self.ctx, &pts, &mut v);
                v
            })
            .collect();
        self.seen(ctr, size, r, |k| {
            let round: [f64; 27] = std::array::from_fn(|j| values[k][j]);
            seen_bounds(&round, &in_part)
        })
    }

    /// A cell's material: (PROPERTY where not resolved, PROPERTY where it
    /// is): its (thin, resolved) bits. Resolved where the cell, as judged,
    /// and every box of the floor's size tiling it pass `variation_resolved`
    /// about their own centres. The cell's own ladder reads what lies in its
    /// outer part only in its whole box, so a film there passed beside a
    /// step or a slope nearer the centre that grew its ranges as a smooth
    /// field's grow (`< 30` held at 25.7 against 35.65). Its eight halves
    /// leave their own outer parts so, the cell's faces' eighth. Every point
    /// lies in a floor box, and from `variation_resolved`'s three conditions
    /// (r <= 1.25 d1, d1 < d2 <= 4 d1) a film whose contrast moves a box's
    /// range by J passes only beside other variation in that same box: in
    /// its quarter (r >= J), at least 0.8 J more over its half; in the ring
    /// to its half, more than J / 2.25 beyond it; past its half, at least
    /// 0.16 J within it (J <= r4 <= r + 5 d1 <= 6.25 d1). What a film can
    /// hide behind lies inside one floor cell, as an interface's does, and
    /// there only the field's slope tells it from what grows as a smooth
    /// field's does: the cell, as judged, is `gentle` too, which no jump in
    /// it is (judged so over every floor box too, a smooth field's solve
    /// took twice as long, for no jump the cell's own box does not hold).
    /// A film on steep finite slopes is as gentle as the box's range lets
    /// it be, which slivers beside it widen: every floor box is `seen` too,
    /// reaching no further than the centres of the floor's boxes round it
    /// take. (Judged by a coarser cell's own centres only, a film within
    /// what they take would pass, though no cell reads it at its strain: a
    /// stiff one at a bent bar's face, as stiff as a centre a cell deeper.)
    pub(super) fn property(&self, ctr: V3, size: f64, level: u8) -> (u8, u8) {
        let Study { ref props, tiles, .. } = *self;
        if props.is_empty() || self.hulls(ctr, size).iter().all(|i| width(i) == 0.0) {
            return (0, PROPERTY);
        }
        let side = size / tiles as f64;
        let judged = tiling(ctr, size, tiles).all(|c| {
            let r = self.hulls(c, side);
            self.resolved_box(c, side, &r) && self.gentle(side, &r, &self.slopes(c, side, &r))
        });
        // Every floor box, `seen` by the centres of the floor's boxes round
        // it: each property read once over them (n + 2 a side).
        let n = tiles << (MAX_DEPTH - level);
        let floor = size / n as f64;
        let floors = || {
            let m = n + 2;
            let lattice: Vec<V3> = (0..m * m * m).map(|i| std::array::from_fn(|a| ctr[a] - 0.5 * size + floor * ((i / m.pow(a as u32) % m) as f64 - 0.5))).collect();
            let mut shape_v = Vec::new();
            eval_points(&self.s.shape, self.ctx, &lattice, &mut shape_v);
            let values: Vec<Vec<f64>> = props
                .iter()
                .map(|t| {
                    let mut v = Vec::new();
                    eval_points(t, self.ctx, &lattice, &mut v);
                    v
                })
                .collect();
            tiling(ctr, size, n).enumerate().all(|(i, c)| {
                let at = |a: usize| i / n.pow(a as u32) % n + 1;
                let idx = |o: [usize; 3]| o[0] + m * (o[1] + m * o[2]);
                let in_part: [bool; 27] = std::array::from_fn(|j| shape_v[idx(std::array::from_fn(|a| at(a) + j / 3usize.pow(a as u32) % 3 - 1))] <= 0.0);
                // A floor box whose centre lies in air and whose entire box is
                // proven outside the part holds no material.
                if !in_part[13] {
                    let half = 0.5 * floor * (1.0 - PROPERTY_SHRINK);
                    let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(c[a] - half, c[a] + half));
                    if eval_interval(&self.s.shape, self.ctx, b, None).lo > 0.0 {
                        return true;
                    }
                }
                let r = self.hulls(c, floor);
                self.resolved_box(c, floor, &r)
                    && self.seen(c, floor, &r, |k| {
                        let v = &values[k];
                        let round: [f64; 27] = std::array::from_fn(|j| v[idx(std::array::from_fn(|a| at(a) + j / 3usize.pow(a as u32) % 3 - 1))]);
                        seen_bounds(&round, &in_part)
                    })
            })
        };
        if judged && floors() { (0, PROPERTY) } else { (PROPERTY, 0) }
    }

    /// A cell's samples (side `size`, at `ctr`): their points, the material
    /// each holds (`occupied`), and those holding any.
    pub(super) fn sampled(&self, t: &Tape, ctr: V3, size: f64) -> (Vec<V3>, Vec<(f64, V3, Cov)>, u64) {
        let Study { ctx, ref xis, .. } = *self;
        let pts: Vec<V3> = xis.iter().map(|&xi| at_xi(ctr, size, xi)).collect();
        let occ = occupied_spread(t, ctx, &pts, xis, size);
        let material = occ.iter().enumerate().fold(0u64, |m, (i, &(o, _, _))| if o > 0.0 { m | 1 << i } else { m });
        (pts, occ, material)
    }

    /// Full cells of side `size` (at `ctrs`): each one's
    /// material by the chord through its centre (`sample_chords`), as its
    /// thin and resolved phases and thinnest feature. Material filling a
    /// cell may still be thinner than two (a plate whose faces lie on the
    /// cells' planes).
    pub(super) fn full_cells(&self, ctrs: &[V3], size: f64) -> Vec<(u8, u8, Option<f64>)> {
        let Study { s, ctx, h, judge, .. } = *self;
        sample_chords(&s.shape, ctx, ctrs, size, THIN * size * (judge / h).max(1.0))
            .into_iter()
            .map(|c| match c[0] {
                Some(v) if v < THIN * size => (MATERIAL, 0, Some(v)),
                Some(v) => (0, MATERIAL, Some(v)),
                None => (0, 0, None),
            })
            .collect()
    }

    /// Of the samples `mask` of a cell (their coefficients `co`)
    /// holding two materials and no other: Some(the axis they lie in layers
    /// along, and which samples along it hold the larger), or Some(None)
    /// where they lie otherwise. None: one material, or more than two (a
    /// gradient).
    pub(super) fn layered(&self, co: &[f64], mask: u64) -> Option<Option<(usize, u8)>> {
        let xis = &self.xis;
        let same = |a: f64, b: f64| (a - b).abs() <= PROPERTY_ROUNDING * a.abs().max(b.abs());
        let of = || (0..co.len()).filter(|i| mask >> i & 1 == 1);
        let (lo, hi) = of().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), i| (l.min(co[i]), h.max(co[i])));
        if same(lo, hi) || !of().all(|i| same(co[i], lo) || same(co[i], hi)) {
            return None;
        }
        let along = |i: usize, a: usize| (((xis[i][a] + 1.0) * 0.5 * K as f64) as usize).min(K - 1);
        Some((0..3).find_map(|a| {
            let mut at = [None; K];
            of().all(|i| *at[along(i, a)].get_or_insert(same(co[i], hi)) == same(co[i], hi)).then(|| (a, (0..K).fold(0u8, |b, k| b | u8::from(at[k] == Some(true)) << k)))
        }))
    }

    /// A cut cell's pieces, each a cell of its own (a gap the samples
    /// resolve splits the cell).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn pieces_of(&self, t: &Arc<Tape>, idx: [i64; 3], level: u8, sub: [u16; 3], size: f64, pts: &[V3], occ: &[(f64, V3, Cov)], material: u64) -> Vec<Classified> {
        let Study { s, ctx, h, judge, w_q, .. } = *self;
        // (A piece holding two materials takes their laminate's element,
        // `laminate`, once split.)
        let co = self.coef_at(pts);
        let nus = self.nu_at(pts);
        // A uniform nu is kept exactly (its mean could round off it).
        let uniform = nus.iter().all(|&v| v == nus[0]).then_some(nus[0]);
        // Chords as long as the thinnest a cell here is split for
        // (`split_as`): two of its own cells, or of the grid it is split as.
        let mut chords = sample_chords(&s.shape, ctx, pts, size, THIN * size * (judge / h).max(1.0));
        // A sample holding a phase no chord's point fell in holds a feature
        // of it thinner than the points' spacing (a foil between them):
        // thin, unmeasured.
        for (c, &(o, _, _)) in chords.iter_mut().zip(occ) {
            if o > 0.0 && c[0].is_none() {
                c[0] = Some(0.0);
            }
            if o < 1.0 && c[1].is_none() {
                c[1] = Some(0.0);
            }
        }
        let [(thin_m, resolved_m), (thin_a, resolved_a)] = thin_and_resolved(occ, &chords, size);
        let air_narrowest = chords.iter().filter_map(|c| c[1]).min_by(f64::total_cmp);
        let mut kept = Vec::new();
        for mask in pieces(material) {
            let mut m = [0.0; 27];
            let (mut vol, mut csum, mut nusum) = (0.0, 0.0, 0.0);
            for (i, ((&(o, xi, cov), c), nu)) in occ.iter().zip(&co).zip(&nus).enumerate() {
                if mask >> i & 1 == 0 {
                    continue;
                }
                let w_q = w_q * o;
                vol += w_q;
                csum += w_q * c;
                nusum += w_q * nu;
                // The sample's material as it spreads about its centroid
                // (`occupied_spread`).
                for (m, v) in m.iter_mut().zip(spread_moments(xi, &cov)) {
                    *m += w_q * c * v;
                }
            }
            kept.push(Classified {
                idx,
                frac: vol / 8.0,
                coef: csum / vol,
                nu: uniform.unwrap_or(nusum / vol),
                moments: Some(m),
                tape: Some(t.clone()),
                mask,
                thin: if thin_m & mask != 0 { MATERIAL } else { 0 } | if thin_a != 0 { AIR } else { 0 },
                resolved: if resolved_m & mask != 0 { MATERIAL } else { 0 } | if resolved_a != 0 { AIR } else { 0 },
                level,
                sub,
                thinnest: (0..64).filter(|&i| mask >> i & 1 == 1).filter_map(|i| chords[i][0]).chain(air_narrowest).min_by(f64::total_cmp),
                lam: None,
            });
        }
        kept
    }

    /// Of a cell not of one material (`PROPERTY`, of side `size` about
    /// `ctr`, at `level`), whether it holds one interface between materials
    /// each at least THIN of its cells thick beside it: not split, as a
    /// thick body's face is not, its band of cells a cell wide on each of a
    /// claim's grids, its error the grid's, which their ladder measures (as a
    /// taper's edge within THIN floor cells of resolved material). Each
    /// material must be the one that reaches the step, joined to it: above
    /// the floor, along an axis, a column THIN cells deep beside the cell on
    /// each side, of one material (one value of each property), grown into
    /// the cell until the two leave less than a floor cell between them,
    /// where any other material lies under the floor; their two materials
    /// the ends of the cell's ranges. (A film of a material, then a sliver of
    /// another, then a block of the first: a block beside the cell reaches
    /// no further than the sliver, so the film is split.) At the floor, each
    /// end held by a box THIN cells across touching the cell, its variation
    /// resolved (a step between slopes, at any angle), its middle quarter
    /// material.
    /// (In air a field may hold a value only air takes: `by-body` continues
    /// a body's past its surface.)
    pub(super) fn one_material(&self, lo: V3, hi: V3) -> Option<Vec<f64>> {
        let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(lo[a], hi[a]));
        let r = self.fields(b);
        r.iter().all(|i| width(i) == 0.0).then(|| r.iter().map(|i| 0.5 * (i.lo + i.hi)).collect())
    }

    /// A box about `p` of side `size` whose middle quarter is all material
    /// (by the shape's bounds): one centred on a face is half air.
    pub(super) fn in_material(&self, p: V3, size: f64) -> bool {
        let Study { s, ctx, .. } = *self;
        let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(p[a] - 0.125 * size, p[a] + 0.125 * size));
        eval_interval(&s.shape, ctx, b, None).hi < 0.0
    }

    /// Whether the stiffness (E or k, or nu) varies over a box of side
    /// `size` about `ctr`.
    pub(super) fn stiffness_varies(&self, ctr: V3, size: f64) -> bool {
        let Study { s, ctx, .. } = *self;
        let half = 0.5 * size * (1.0 - PROPERTY_SHRINK);
        let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - half, ctr[a] + half));
        [s.coef_const.is_none().then_some(&*s.coef), s.nu.as_deref()].into_iter().flatten().any(|t| width(&eval_interval(t, ctx, b, None)) > 0.0)
    }

    /// Read the held values at moment `t` from now on (`held_t`).
    pub(super) fn hold_at(&self, t: f64) {
        self.held_t.store(t.to_bits(), Ordering::Relaxed);
    }

    pub(super) fn held_time(&self) -> f64 {
        f64::from_bits(self.held_t.load(Ordering::Relaxed))
    }

    // A fix's held value at points (each axis: a number or a field), at
    // the moment `held_t`.
    pub(super) fn values_at(&self, fx: &Fix, pts: &[V3]) -> Vec<[f64; 3]> {
        let (t, at);
        let ctx = if (t = self.held_time(), t.to_bits() == self.ctx.t.to_bits()).1 {
            self.ctx
        } else {
            at = Ctx { t, ..self.ctx.clone() };
            &at
        };
        let mut out = vec![[0.0; 3]; pts.len()];
        for (d, v) in fx.value.iter().enumerate() {
            match v {
                Held::Number(x) => out.iter_mut().for_each(|o| o[d] = *x),
                Held::Field(t) => {
                    let mut vals = Vec::new();
                    eval_points(t, ctx, pts, &mut vals);
                    out.iter_mut().zip(vals).for_each(|(o, x)| o[d] = x);
                }
            }
        }
        out
    }

    pub(super) fn region_overlaps(&self, r: &RegionSel, idx: [i64; 3]) -> bool {
        let Study { ctx, ref grid, .. } = *self;
        r.may_touch(ctx, grid.cell_box(Cell { i: idx, size: 1 }))
    }

    /// The part's surface, meshed at half the cells it lies in (`part_surface`).
    pub(super) fn whole_surface(&self, split: &[([i64; 3], u8)]) -> Vec<Facet> {
        let Study { s, ctx, ref grid, .. } = *self;
        part_surface(&s.shape, ctx, grid, s.lo, s.hi, split)
    }

    /// The mirror a facet lies on, if one.
    pub(super) fn cut_face(&self, fc: &Facet) -> Option<(usize, f64)> {
        let Study { h, ref mirrors, .. } = *self;
        mirrors.iter().copied().find(|&(axis, at)| fc.normal[axis].abs() > 1.0 - 1e-8 && fc.corners.iter().all(|p| (p[axis] - at).abs() <= 1e-8 * h))
    }

    /// Its facets' points on the surface (`surface_samples`), for regions.
    pub(super) fn sample_faces(&self, facets: &[Facet]) -> Vec<V3> {
        let Study { s, ctx, ref mirrors, .. } = *self;
        let mut points = surface_samples(&s.shape, ctx, facets);
        // The surface projection stops at a finite distance. On a known
        // planar cut use its exact coordinate: a sample nudged into
        // the solid can otherwise put bore pressure on an internal face.
        if !mirrors.is_empty() {
            for (fc, samples) in facets.iter().zip(points.as_chunks_mut::<4>().0) {
                if let Some((axis, at)) = self.cut_face(fc) {
                    for p in samples {
                        p[axis] = at;
                    }
                }
            }
        }
        points
    }

    pub(super) fn material_held(&self, ctr: V3, size: f64, level: u8) -> bool {
        let Study { s, ctx, .. } = *self;
        let r = self.hulls(ctr, size);
        let varies: Vec<usize> = (0..r.len()).filter(|&k| width(&r[k]) > 0.0).collect();
        if varies.is_empty() {
            return true;
        }
        let half = 0.5 * size * (1.0 - PROPERTY_SHRINK);
        if level < MAX_DEPTH {
            let n = 1usize << (MAX_DEPTH - level);
            let floor = size / n as f64;
            (0..3).any(|a| {
                // The column on side `sg` grown `k` floor cells into the
                // cell, if of one material: its values.
                let column = |sg: f64, k: usize| -> Option<Vec<f64>> {
                    let (mut lo, mut hi): (V3, V3) = (ctr.map(|c| c - half), ctr.map(|c| c + half));
                    let face = ctr[a] + sg * 0.5 * size;
                    let (near, far) = (face - sg * k as f64 * floor, face + sg * THIN * size);
                    (lo[a], hi[a]) = (near.min(far), near.max(far));
                    self.one_material(lo, hi)
                };
                // How far the column reaches into the cell, and its values.
                // Material all along it: each quarter of its depth holds some
                // (by the shape's bounds), so air beyond a face holds none
                // of it.
                let reach = |sg: f64| -> Option<(usize, Vec<f64>)> {
                    let deep = (0..4).all(|j| {
                        let (n, f) = (ctr[a] + sg * (0.5 + 0.5 * j as f64) * size, ctr[a] + sg * (0.5 + 0.5 * (j + 1) as f64) * size);
                        let b: [interval::I; 3] =
                            std::array::from_fn(|b| if b == a { interval::I::new(n.min(f), n.max(f)) } else { interval::I::new(ctr[b] - half, ctr[b] + half) });
                        eval_interval(&s.shape, ctx, b, None).lo < 0.0
                    });
                    let values = column(sg, 0).filter(|_| deep)?;
                    let (mut k, mut out) = (0, n);
                    while out - k > 1 {
                        let mid = (k + out) / 2;
                        if column(sg, mid).is_some() { k = mid } else { out = mid }
                    }
                    Some((k, values))
                };
                let (Some((kp, up)), Some((km, um))) = (reach(1.0), reach(-1.0)) else { return false };
                kp + km + 1 >= n
                    && varies.iter().all(|&k| {
                        let (lo, hi) = (r[k].lo, r[k].hi);
                        (same(up[k], lo) && same(um[k], hi)) || (same(up[k], hi) && same(um[k], lo))
                    })
            })
        } else {
            let mut ends: Vec<(usize, f64)> = varies.iter().flat_map(|&k| [(k, r[k].lo), (k, r[k].hi)]).collect();
            // Touching the cell: its middle 1.5 cells off along each axis
            // it is beside (a hair less, so its box shares the cell's face).
            for o in (0..27).map(|i| [i % 3, i / 3 % 3, i / 9].map(|v| v as f64 - 1.0)).filter(|o| *o != [0.0; 3]) {
                let p: V3 = std::array::from_fn(|a| ctr[a] + o[a] * (1.5 - 1e-5) * size);
                if !self.in_material(p, THIN * size) {
                    continue;
                }
                let r = self.hulls(p, THIN * size);
                let holds = |&(k, v): &(usize, f64)| r[k].lo <= v && v <= r[k].hi;
                if !ends.iter().any(holds) || !self.resolved_box(p, THIN * size, &r) || {
                    !self.gentle(THIN * size, &r, &self.slopes(p, THIN * size, &r)) || !self.seen_round(p, THIN * size, &r)
                } {
                    continue;
                }
                ends.retain(|e| !holds(e));
                if ends.is_empty() {
                    return true;
                }
            }
            false
        }
    }
}

/// Property bounds seen at a box's centre and its 26 same-size neighbours
/// inside the part (`in_part[j]` = `shape(p_j) <= 0.0`). Neighbours outside
/// the part (in air) and non-finite values are ignored, so neither a skin
/// clamped into air nor an inward dip extrapolated across the centre can
/// widen the allowed range beyond values actually taken inside the part.
fn seen_bounds(round: &[f64; 27], in_part: &[bool; 27]) -> (f64, f64, f64) {
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for j in 0..27 {
        if in_part[j] && round[j].is_finite() {
            lo = lo.min(round[j]);
            hi = hi.max(round[j]);
        }
    }
    (round[13], lo, hi)
}

// A property's range `r` over a box reaching `h` from its centre along
// each axis (its slopes `g`, its centre's value `at`), no wider than the
// slopes allow from that value: the mean value theorem, along the segment
// from the centre.
fn within(r: interval::I, at: f64, g: &interval::G, h: f64) -> interval::I {
    let reach = h * g.d.iter().map(|d| interval::mag(*d)).sum::<f64>().min(3f64.sqrt() * g.l);
    let (lo, hi) = (r.lo.max(at - reach), r.hi.min(at + reach));
    if lo <= hi { interval::I { lo, hi } } else { r }
}

// A range within rounding of its value (bounds round outward at each
// step of the tape: 200 x 1e-3 is not a point) is none.
pub(super) fn width(i: &interval::I) -> f64 {
    match i.hi - i.lo {
        w if !w.is_finite() => f64::INFINITY,
        w if w <= PROPERTY_ROUNDING * i.lo.abs().max(i.hi.abs()) => 0.0,
        w => w,
    }
}

pub(super) fn same(a: f64, b: f64) -> bool {
    (a - b).abs() <= PROPERTY_ROUNDING * a.abs().max(b.abs())
}

// The boxes of side `w`, `n` a side, tiling a cell of side `size` about
// `ctr`.
pub(super) fn tiling(ctr: V3, size: f64, n: usize) -> impl Iterator<Item = V3> {
    let w = size / n as f64;
    (0..n * n * n).map(move |i| -> V3 { std::array::from_fn(|a| ctr[a] - 0.5 * size + w * ((i / n.pow(a as u32) % n) as f64 + 0.5)) })
}
