//! What a cut cell holds: each sample's share of material and where it
//! lies, the moments they make, and the faces of the regions within it.

use crate::cells::{K, Placed, locate, place};
use crate::eval::{Ctx, eval_interval, eval_points};
use crate::hold;
use crate::laminate::Laminate;
use crate::octree::{Cell, Grid};
use crate::tape::Tape;
use ink_core::{V3, interval, v3};
use std::sync::Arc;

#[derive(Clone)]
pub(crate) struct Classified {
    pub(crate) idx: [i64; 3],
    pub(crate) frac: f64,
    pub(crate) coef: f64,
    /// Elastic: mean Poisson's ratio over the material.
    pub(crate) nu: f64,
    pub(crate) moments: Option<[f64; 27]>,
    /// For cut cells: the shape's tape simplified for this cell by the
    /// octree (much shorter than the whole shape for big unions).
    pub(crate) tape: Option<Arc<Tape>>,
    /// The samples of this piece of the cell's material (see SimCell).
    pub(crate) mask: u64,
    /// The phases (`MATERIAL`, `AIR`) of which its samples hold a thin
    /// feature, and those resolved in them (`thin_and_resolved`): its
    /// piece's material, the cell's air; a full cell's by the chord through
    /// its centre (material filling a cell may still be thinner than two:
    /// a fin's middle row of cells).
    pub(crate) thin: u8,
    pub(crate) resolved: u8,
    /// Refinement, as SimCell's.
    pub(crate) level: u8,
    pub(crate) sub: [u16; 3],
    /// The thinnest feature its samples lie in, of either phase (its
    /// piece's material, the cell's air: `sample_chords`).
    pub(crate) thinnest: Option<f64>,
    /// Two materials across a face (`laminate`).
    pub(crate) lam: Option<Box<Laminate>>,
}

/// Per sample (at xi in a cell of size h): the material fraction of its
/// sub-cell, and where that material's centroid sits (in xi). Integrals over
/// the material use the centroid: lumping a partly filled sub-cell at its
/// centre misplaces the boundary by up to a quarter sub-cell.
pub(crate) fn occupied(tape: &Tape, ctx: &Ctx, pts: &[V3], xis: &[V3], h: f64) -> Vec<(f64, V3)> {
    occupied_spread(tape, ctx, pts, xis, h).into_iter().map(|(o, at, _)| (o, at)).collect()
}

/// A sample's material's covariance about its centroid: xx, yy, zz, xy,
/// yz, xz.
pub(crate) type Cov = [f64; 6];

/// `occupied`, with the covariance of each sample's material about its
/// centroid (in xi^2). A cell's moments are each sample's own with it:
/// lumped at its centroid, a sheet one sample thick had none through its
/// thickness and bent freely (a 0.05 mm ligament 1000 times too soft), and
/// a whole sample was 6% short; spread along each axis alone, a sheet
/// lying across the axes was as thick as the sample (6 times too stiff).
pub(crate) fn occupied_spread(tape: &Tape, ctx: &Ctx, pts: &[V3], xis: &[V3], h: f64) -> Vec<(f64, V3, Cov)> {
    let to_xi = 2.0 / h;
    occupied_parts(tape, ctx, pts, h).into_iter().zip(xis).map(|((o, off, cov), xi)| (o, std::array::from_fn(|k| xi[k] + off[k] * to_xi), cov.map(|v| v * to_xi * to_xi))).collect()
}

/// The moments x^a y^b z^c (a, b, c each 0, 1 or 2; `m[a + 3 b + 9 c]`)
/// of material of unit volume with centroid `c` and covariance `cov`,
/// the central moments of four and six coordinates from the covariance
/// (Isserlis: exact for independent coordinates, a sample's own box).
pub(crate) fn spread_moments(c: V3, cov: &Cov) -> [f64; 27] {
    let s = |i: usize, j: usize| -> f64 {
        match (i.min(j), i.max(j)) {
            (0, 0) => cov[0],
            (1, 1) => cov[1],
            (2, 2) => cov[2],
            (0, 1) => cov[3],
            (1, 2) => cov[4],
            _ => cov[5],
        }
    };
    if cov[3] == 0.0 && cov[4] == 0.0 && cov[5] == 0.0 {
        let p: [[f64; 3]; 3] = std::array::from_fn(|a| [1.0, c[a], c[a] * c[a] + cov[a]]);
        return std::array::from_fn(|k| p[0][k % 3] * p[1][(k / 3) % 3] * p[2][k / 9]);
    }
    // Central moment of the coordinates listed (each axis at most twice).
    fn central(ax: &mut Vec<usize>, s: &dyn Fn(usize, usize) -> f64) -> f64 {
        if ax.is_empty() {
            return 1.0;
        }
        if ax.len() % 2 == 1 {
            return 0.0;
        }
        let first = ax.remove(0);
        let mut sum = 0.0;
        for j in 0..ax.len() {
            let other = ax.remove(j);
            sum += s(first, other) * central(ax, s);
            ax.insert(j, other);
        }
        ax.insert(0, first);
        sum
    }
    let mut cm = [0.0; 27];
    for (k, v) in cm.iter_mut().enumerate() {
        let mut ax: Vec<usize> = (0..3).flat_map(|a| std::iter::repeat_n(a, [k % 3, (k / 3) % 3, k / 9][a])).collect();
        *v = central(&mut ax, &s);
    }
    let binom = |n: usize, j: usize| if n == 2 && j == 1 { 2.0 } else { 1.0 };
    std::array::from_fn(|k| {
        let e = [k % 3, (k / 3) % 3, k / 9];
        let mut sum = 0.0;
        for j in 0..=e[2] {
            for i in 0..=e[1] {
                for h in 0..=e[0] {
                    let w = binom(e[0], h) * binom(e[1], i) * binom(e[2], j) * c[0].powi((e[0] - h) as i32) * c[1].powi((e[1] - i) as i32) * c[2].powi((e[2] - j) as i32);
                    sum += w * cm[h + 3 * i + 9 * j];
                }
            }
        }
        sum
    })
}

/// Per sample: the material fraction of its sub-cell (side h/K), its
/// centroid's offset from the sub-cell's centre (mm), and its covariance
/// about its centroid (mm^2).
///
/// The fraction is the sub-cell's part on the material side of the field's
/// tangent plane at its centre. Past a crease of the field (a max of two
/// branches: where a cut's face lies flush with the part's, the field in
/// the cut is |z - 5|, never negative) that plane puts material where there
/// is none: the sub-cell straddling the flush face read a third material,
/// a film bridging the cut's mouth that held a slit shut. So where the
/// material it estimates has its centroid in air, the part is the
/// sub-cell's on the material side of both tangent planes, the centre's
/// and the centroid's (exact for two planes); a union's crease, whose
/// estimate falls short, is left as it is. The covariance: of the
/// sample's part on the material side of the tangent plane (`clipped_cov`; s^2 / 12 on each axis across a whole
/// sample), or, past a crease, of the part measured at its points.
pub(crate) fn occupied_parts(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64) -> Vec<(f64, V3, Cov)> {
    let s = h / K as f64;
    let (tangent, g) = tangent_offsets(tape, ctx, pts, h);
    let whole = s * s / 12.0;
    let mut out: Vec<(f64, V3, Cov)> =
        tangent.into_iter().zip(&g).map(|((o, off), g)| (o, off, if o > 0.0 && o < 1.0 { clipped_cov(g, s) } else { [whole, whole, whole, 0.0, 0.0, 0.0] })).collect();
    // The tangent plane places the surface by the field's value over its slope, the depth only of a
    // distance. Where a sample's whole box is proven by the field's bounds to hold no material, or all
    // material, it is so, whatever the plane said: a bore written as the cube of a distance put
    // material in its air samples (at a third of their depth), and their chords, finding none, read it
    // as a foil thinner than they, which split the bore to the floor.
    for (i, o) in out.iter_mut().enumerate().filter(|(_, o)| o.0 > 0.0 && o.0 < 1.0) {
        let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(pts[i][a] - 0.5 * s, pts[i][a] + 0.5 * s));
        let r = eval_interval(tape, ctx, b, None);
        if r.lo > 0.0 {
            *o = (0.0, [0.0; 3], [0.0; 6]);
        } else if r.hi <= 0.0 {
            *o = (1.0, [0.0; 3], [whole, whole, whole, 0.0, 0.0, 0.0]);
        } else if r.lo == 0.0 {
            // A sample the part only touches, on its face or along its edge
            // (a face on the cells' plane, the sample beyond it), holds
            // none: its box a billionth smaller is all air. The plane at
            // its centre, a distance's past a box's edge, gave it a
            // fifteenth of itself, and a cell beyond a face on its plane a
            // sliver of material that is not there.
            let within: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(pts[i][a] - 0.5 * s * (1.0 - 1e-9), pts[i][a] + 0.5 * s * (1.0 - 1e-9)));
            if eval_interval(tape, ctx, within, None).lo > 0.0 {
                *o = (0.0, [0.0; 3], [0.0; 6]);
            }
        }
    }
    // A sample that several flat faces of the part cross (an edge's, a
    // corner's) is its part within them all, or outside them all at an
    // inside corner, exactly: where the field at sixteen points of the
    // sample is, to a rounding, the largest or the smallest of the planes
    // read at them (`region_faces`; past an outside edge a distance reads
    // the plane through the edge, which changes nothing of what lies
    // within). Its nearer face's plane alone counted the material past the
    // other: along every edge the cells held too much, a quarter of a
    // sample's section at most, by where the edge fell in the sample and
    // not by the cells' size, so three grids could agree on it. A 20 x 2 x
    // 2 bar with two sides 0.1 mm past the planes of 0.5 mm cells held
    // 0.16% too much (80.125 mm^3), and read its load's work 0.16% short
    // and its end's stretch with it.
    //
    // Not a wall thinner than the sample, two flat faces back to back
    // inside it: that stays as below, counted by 8 points across the
    // sample or not at all. Counted exactly, a foil under 3% of its cells
    // (`MIN_FRACTION`) is seen and then dropped with the slivers, and
    // nothing says so; as it is, its least count is an eighth of a sample,
    // 1/32 of the cell and over the floor, and a foil the points miss is
    // found by the lines between the samples and told (`slivers`). Exact
    // there needs the dropped told first.
    let crossed: Vec<usize> = (0..pts.len()).filter(|&i| out[i].0 > 0.0 && out[i].0 < 1.0).collect();
    let (_, faces) = region_faces(tape, ctx, &crossed.iter().map(|&i| pts[i]).collect::<Vec<V3>>(), h, true, true);
    let mut several = vec![false; pts.len()];
    for (&i, f) in crossed.iter().zip(&faces).filter(|x| x.1.planes.len() > 1) {
        out[i] = sets_part(None, &[f.set(true)], s);
        several[i] = true;
    }
    let partial: Vec<usize> = (0..pts.len()).filter(|&i| out[i].0 > 0.0 && out[i].0 < 1.0 && !several[i]).collect();
    if !partial.is_empty() {
        let at: Vec<V3> = partial.iter().map(|&i| v3::add(pts[i], out[i].1)).collect();
        let g2 = crate::eval::eval_grads(tape, ctx, &at);
        for ((&i, q), g2) in partial.iter().zip(&at).zip(&g2) {
            if g2[0] > 0.0 {
                out[i] = two_plane_occupancy(&g[i], pts[i], g2, *q, s);
            }
        }
    }
    out
}

/// Per sample: its sub-cell's part on the material side of the field's
/// tangent plane at its centre, and that part's centroid's offset (mm); and
/// the field's values and gradients there. A region is this: one drawn
/// thinner than a sample (a symmetry plane's layer, an edge's strip) holds
/// what it lies on, as it is meant to.
pub(crate) fn tangent_offsets(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64) -> (Vec<(f64, V3)>, Vec<[f64; 4]>) {
    let (out, g, _) = tangent_near(tape, ctx, pts, h);
    (out, g)
}

/// `tangent_offsets`, and the samples whose centres the field reads within
/// a sample of its zero.
fn tangent_near(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64) -> (Vec<(f64, V3)>, Vec<[f64; 4]>, Vec<usize>) {
    let s = h / K as f64;
    let mut g = crate::eval::eval_grads(tape, ctx, pts);
    // The plane is the field's tangent at the sample, placed where the field's zero lies along its
    // gradient: value over slope places it there only for a distance (a cube of one put it at a third of
    // its depth). Within a sample of the plane's place, the zero is found along the gradient by the
    // secant on the field itself, and the plane moved to it; a distance's lands at once, unmoved.
    let slope = |g: &[f64; 4]| (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    let near: Vec<usize> = (0..pts.len()).filter(|&i| slope(&g[i]) > 0.0 && (g[i][0] / slope(&g[i])).abs() < s).collect();
    if !near.is_empty() {
        let dir: Vec<V3> = near
            .iter()
            .map(|&i| {
                let k = slope(&g[i]);
                [g[i][1] / k, g[i][2] / k, g[i][3] / k]
            })
            .collect();
        // t: the zero's place along the gradient (mm), from the plane's.
        let mut t: Vec<f64> = near.iter().map(|&i| -g[i][0] / slope(&g[i])).collect();
        let (mut t0, mut f0): (Vec<f64>, Vec<f64>) = (vec![0.0; near.len()], near.iter().map(|&i| g[i][0]).collect());
        for _ in 0..SECANT_STEPS {
            let at: Vec<V3> = near.iter().zip(&dir).zip(&t).map(|((&i, d), t)| v3::add(pts[i], v3::scale(*d, *t))).collect();
            let mut f = Vec::new();
            eval_points(tape, ctx, &at, &mut f);
            let mut moved = false;
            for k in 0..near.len() {
                let (df, dt) = (f[k] - f0[k], t[k] - t0[k]);
                if f[k] == 0.0 || df == 0.0 || !(df * dt > 0.0) {
                    continue;
                }
                let next = t[k] - f[k] * dt / df;
                if !(next.abs() < 2.0 * s) {
                    continue;
                }
                (t0[k], f0[k]) = (t[k], f[k]);
                moved |= (next - t[k]).abs() > 1e-9 * s;
                t[k] = next;
            }
            if !moved {
                break;
            }
        }
        for (k, &i) in near.iter().enumerate() {
            g[i][0] = -slope(&g[i]) * t[k];
        }
    }
    let out = g
        .iter()
        .map(|g| {
            let o = crate::query::occupancy(g, [s; 3]);
            (o, if o > 0.0 && o < 1.0 { occupied_centroid(g, s) } else { [0.0; 3] })
        })
        .collect();
    (out, g, near)
}

/// A region's faces through a sample: planes (the field's value and
/// gradient, as from the sample's centre), the region where all of them
/// are negative, or (`any`) where any one is.
#[derive(Clone, Debug)]
pub(crate) struct Faces {
    pub(crate) planes: Vec<[f64; 4]>,
    pub(crate) any: bool,
}

impl Faces {
    /// The region's part of a sample, or the rest of it (`sets_part`'s set).
    pub(crate) fn set(&self, inside: bool) -> (Vec<[f64; 4]>, bool) {
        // (the part all the planes bound; any of them, what all leave out)
        if self.any { (self.planes.iter().map(|g| g.map(|v| -v)).collect(), !inside) } else { (self.planes.clone(), inside) }
    }

    /// The region's field at `d` from the sample's centre as its planes
    /// give it, each over its slope (mm): negative inside.
    pub(crate) fn at(&self, d: V3) -> f64 {
        let v = self.planes.iter().map(|g| (g[0] + g[1] * d[0] + g[2] * d[1] + g[3] * d[2]) / (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt().max(1e-300));
        if self.any { v.fold(f64::INFINITY, f64::min) } else { v.fold(f64::NEG_INFINITY, f64::max) }
    }
}

/// How far along the line between two of a sample's points the field's
/// gradient is read to tell a corner from a curve (`region_faces`): 1 /
/// phi^2, a share no face written at a round place in a sample's quarters
/// crosses the line at.
const BETWEEN: f64 = 0.381_966_011_250_105;

/// `tangent_offsets`, and each sample's faces of the field's region: one
/// tangent plane where one face crosses it, flat or curved, as before; and
/// where several meet in it (a corner or an edge of the region, a slab
/// thinner than the sample), with `several`, each one's plane, and whether
/// the region is inside all or any (a max of fields or a min). The field's
/// plane is read at sixteen points of the sample (its eight corners and
/// the eight a quarter of it from its centre, each a millionth of it off,
/// that none lie on a tie between faces at whole slopes). Two points
/// beside each other read one face where their planes agree, or where the
/// gradient `BETWEEN` them is nearer the middle of theirs than either: it
/// turns with a curved face, and at a corner it is one face's or the
/// other's. The faces so found are kept where, together, they read the
/// field at the sixteen points nearer than the one plane does; else the
/// one plane stays, the nearest of them and it. (One plane a sample, a
/// block's face inside a sample
/// ran on past its other faces and the second of two was not traced: 4%
/// of its boundary on 0.3 mm cells, its reading 1.8% apart across a cell's
/// plane; and at a tie between three faces the plane was their mean,
/// which cut a quarter of a bracket's corner sample off as free.)
///
/// `flat`, for the part's own material (`occupied_parts`): only the flat
/// faces that cross a sample, and only where they read the field at the
/// sixteen points to a rounding; a curved face's sample, and a wall's
/// thinner than the sample, keep their one plane.
pub(crate) fn region_faces(tape: &Tape, ctx: &Ctx, pts: &[V3], h: f64, several: bool, flat: bool) -> (Vec<(f64, V3)>, Vec<Faces>) {
    let s = h / K as f64;
    let (mut share, g, near) = tangent_near(tape, ctx, pts, h);
    let mut faces: Vec<Faces> = g.iter().map(|g| Faces { planes: vec![*g], any: false }).collect();
    let slope = |g: &[f64; 4]| (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    if near.is_empty() {
        return (share, faces);
    }
    // (Each point a millionth of a sample off its place, by three unlike
    // amounts: a field written as a max of faces has no gradient on their
    // tie, and a point exactly on one reads either face's. The answers
    // move with them by their size.)
    let offset = |q: usize| -> V3 {
        let r = if q < 8 { 0.25 * s } else { 0.5 * s };
        std::array::from_fn(|a| if q >> a & 1 == 1 { r } else { -r } + [1e-6, 1.618e-6, 2.718e-6][a] * s)
    };
    let value = |g: &[f64; 4], d: V3| g[0] + g[1] * d[0] + g[2] * d[1] + g[3] * d[2];
    let unit = |g: &[f64; 4]| -> V3 { [g[1], g[2], g[3]].map(|v| v / slope(g).max(1e-300)) };
    let read = crate::eval::eval_grads(tape, ctx, &near.iter().flat_map(|&i| (0..16).map(move |q| (i, q))).map(|(i, q)| v3::add(pts[i], offset(q))).collect::<Vec<V3>>());
    // Each point's plane, as from its sample's centre.
    let plane = |n: usize, q: usize| -> [f64; 4] {
        let (r, d) = (read[16 * n + q], offset(q));
        [r[0] - (r[1] * d[0] + r[2] * d[1] + r[3] * d[2]), r[1], r[2], r[3]]
    };
    // The points beside each other: the inner eight along the axes, and each
    // corner with the inner one before it.
    let beside: Vec<(usize, usize)> = (0..8).flat_map(|p| (0..3).filter(move |a| p >> a & 1 == 0).map(move |a| (p, p | 1 << a))).chain((0..8).map(|p| (p, 8 + p))).collect();
    let apart = |n: usize, p: usize, q: usize| v3::norm(v3::sub(unit(&read[16 * n + p]), unit(&read[16 * n + q])));
    // A sample all of whose points read its own plane holds one flat face.
    let one_plane = |n: usize| {
        let gi = &g[near[n]];
        (0..16).all(|q| {
            let l = plane(n, q);
            v3::norm(v3::sub(unit(&l), unit(gi))) <= 1e-9 && (l[0] - gi[0]).abs() <= 1e-9 * s * slope(gi)
        })
    };
    let open: Vec<usize> = (0..near.len()).filter(|&n| !one_plane(n)).collect();
    // (the part's flat faces alone, `flat`, ask nothing between the points)
    let asked: Vec<(usize, usize, usize)> =
        open.iter().filter(|_| !flat).flat_map(|&n| beside.iter().map(move |&(p, q)| (n, p, q))).filter(|&(n, p, q)| apart(n, p, q) > 1e-9).collect();
    let between = crate::eval::eval_grads(
        tape,
        ctx,
        &asked
            .iter()
            .map(|&(n, p, q)| {
                let (a, b) = (offset(p), offset(q));
                v3::add(pts[near[n]], std::array::from_fn(|k| a[k] + BETWEEN * (b[k] - a[k])))
            })
            .collect::<Vec<V3>>(),
    );
    let mut turn = asked.iter().zip(&between).peekable();
    // (`flat`: the samples whose planes read the field at the sixteen points
    // to a rounding, each with its planes, each plane's point and how many
    // points read it, and whether the region is inside any of them)
    let mut exact_at: Vec<(usize, Vec<[f64; 4]>, Vec<(usize, usize)>, bool)> = Vec::new();
    for &n in &open {
        let i = near[n];
        // The points reading one face, joined.
        let mut of: Vec<usize> = (0..16).collect();
        let root = |of: &mut Vec<usize>, mut p: usize| {
            while of[p] != p {
                p = of[p];
            }
            p
        };
        for &(p, q) in beside.iter().filter(|_| !flat) {
            let delta = apart(n, p, q);
            let one = if delta <= 1e-9 {
                true
            } else {
                let (_, m) = turn.next().expect("a gradient asked for");
                let to = |e: usize| v3::norm(v3::sub(unit(m), unit(&read[16 * n + e])));
                to(p).min(to(q)) >= 0.5 * BETWEEN * delta
            };
            if one {
                let (a, b) = (root(&mut of, p), root(&mut of, q));
                of[a.max(b)] = a.min(b);
            }
        }
        // (and any two whose planes are one, wherever in the sample: a flat
        // face's points lie either side of another's)
        for p in 0..16 {
            for q in p + 1..16 {
                let (a, b) = (plane(n, p), plane(n, q));
                if apart(n, p, q) <= 1e-9 && (a[0] - b[0]).abs() <= 1e-9 * s * slope(&a) {
                    let (a, b) = (root(&mut of, p), root(&mut of, q));
                    of[a.max(b)] = a.min(b);
                }
            }
        }
        // Each face's plane: its point nearest the field's zero.
        let mut planes: Vec<[f64; 4]> = Vec::new();
        let mut read_at: Vec<(usize, usize)> = Vec::new();
        for r in (0..16).filter(|&r| root(&mut of.clone(), r) == r) {
            let its: Vec<usize> = (0..16).filter(|&q| root(&mut of.clone(), q) == r).collect();
            let q = its
                .iter()
                .copied()
                .min_by(|&a, &b| (read[16 * n + a][0] / slope(&read[16 * n + a])).abs().total_cmp(&(read[16 * n + b][0] / slope(&read[16 * n + b])).abs()))
                .expect("a face's point");
            planes.push(plane(n, q));
            read_at.push((q, its.len()));
        }
        if planes.len() == 1 {
            continue;
        }
        // How far each reading of the field is from it at the sixteen points.
        let off = |model: &dyn Fn(V3) -> f64| (0..16).map(|q| (model(offset(q)) - read[16 * n + q][0]).abs()).fold(0.0, f64::max);
        let all = off(&|d| planes.iter().map(|l| value(l, d)).fold(f64::NEG_INFINITY, f64::max));
        let any = off(&|d| planes.iter().map(|l| value(l, d)).fold(f64::INFINITY, f64::min));
        let one = off(&|d| value(&g[i], d));
        if flat {
            if several && all.min(any) <= 1e-9 * s * slope(&g[i]) && read_at.iter().any(|r| r.1 > 1) {
                exact_at.push((n, planes, read_at, any < all));
            }
            continue;
        }
        if several && all.min(any) < one {
            let f = Faces { planes, any: any < all };
            let (v, c, _) = sets_part(None, &[f.set(true)], s);
            share[i] = (v, if v > 0.0 && v < 1.0 { c } else { [0.0; 3] });
            faces[i] = f;
        }
    }
    // `flat`: a sample's flat faces. The largest of the tangent planes of
    // any field that bulges one way (a distance to a ball or a round bar,
    // or past a box's edge, where it is the distance to the edge) reads it
    // at their own points to a rounding, so that alone tells nothing of
    // flat. A plane is a flat face's where the field's plane a
    // ten-millionth of a sample from its point is still it (asked only in
    // a sample where two points or more read one plane: a ball's have
    // none); the others must cut nothing from the sample's part within the
    // flat ones (a box's edge's do not: each passes through the edge), or
    // the sample stays as its one plane reads it.
    let step: V3 = [0.618e-7 * s, -0.414e-7 * s, 0.732e-7 * s];
    let asked: Vec<(usize, usize)> = exact_at.iter().enumerate().flat_map(|(k, e)| (0..e.2.len()).map(move |j| (k, j))).collect();
    let again = crate::eval::eval_grads(tape, ctx, &asked.iter().map(|&(k, j)| v3::add(v3::add(pts[near[exact_at[k].0]], offset(exact_at[k].2[j].0)), step)).collect::<Vec<V3>>());
    let mut is_flat: Vec<Vec<bool>> = exact_at.iter().map(|e| vec![false; e.1.len()]).collect();
    for (&(k, j), r) in asked.iter().zip(&again) {
        let (l, d) = (&exact_at[k].1[j], v3::add(offset(exact_at[k].2[j].0), step));
        is_flat[k][j] = v3::norm(v3::sub(unit(r), unit(l))) <= 1e-12 && (r[0] - (r[1] * d[0] + r[2] * d[1] + r[3] * d[2]) - l[0]).abs() <= 1e-9 * s * slope(l);
    }
    for ((n, planes, _, any), is_flat) in exact_at.into_iter().zip(is_flat) {
        let f = Faces { planes: planes.iter().zip(&is_flat).filter(|x| *x.1).map(|x| *x.0).collect(), any };
        // (Two of them back to back, within a tenth of a radian: a wall or
        // a gap thinner than the sample, which stays as it is read today.
        // `occupied_parts` says why.)
        let back_to_back = f.planes.iter().enumerate().any(|(j, a)| f.planes[..j].iter().any(|b| v3::dot(unit(a), unit(b)) < -0.995));
        if f.planes.len() < 2 || back_to_back {
            continue;
        }
        let (v, c, _) = sets_part(None, &[f.set(true)], s);
        if (sets_part(None, &[Faces { planes, any }.set(true)], s).0 - v).abs() <= 1e-9 {
            share[near[n]] = (v, if v > 0.0 && v < 1.0 { c } else { [0.0; 3] });
            faces[near[n]] = f;
        }
    }
    (share, faces)
}

/// Secant steps `tangent_offsets` takes to the field's zero along its gradient: one for a distance (it
/// lands), a few for a power of one.
pub(crate) const SECANT_STEPS: usize = 8;

/// The part of the cube of side `s` centred at `p` on the material side of
/// two tangent planes (field value and gradient `a` at `p`, `b` at `q`):
/// its share, its centroid's offset and its covariance about it, by 8^3
/// points. (The material's own; a hold's part of a sample is exact,
/// `sets_part`. Exact here too moves what is read of every curved part: a
/// change of its own.)
fn two_plane_occupancy(a: &[f64; 4], p: V3, b: &[f64; 4], q: V3, s: f64) -> (f64, V3, Cov) {
    const N: usize = 8;
    let b = [b[0] + (0..3).map(|k| b[1 + k] * (p[k] - q[k])).sum::<f64>(), b[1], b[2], b[3]];
    let (mut n, mut sum, mut sq) = (0usize, [0.0; 3], [0.0; 6]);
    for i in 0..N * N * N {
        let d: V3 = std::array::from_fn(|k| ((([i % N, (i / N) % N, i / (N * N)][k]) as f64 + 0.5) / N as f64 - 0.5) * s);
        if [a, &b].iter().any(|g| g[0] + g[1] * d[0] + g[2] * d[1] + g[3] * d[2] > 0.0) {
            continue;
        }
        n += 1;
        sum = v3::add(sum, d);
        let pairs = [d[0] * d[0], d[1] * d[1], d[2] * d[2], d[0] * d[1], d[1] * d[2], d[0] * d[2]];
        sq.iter_mut().zip(pairs).for_each(|(a, b)| *a += b);
    }
    if n == 0 {
        return (0.0, [0.0; 3], [0.0; 6]);
    }
    let (c, cell): (V3, f64) = (sum.map(|v| v / n as f64), (s / N as f64).powi(2) / 12.0);
    let centred = [c[0] * c[0], c[1] * c[1], c[2] * c[2], c[0] * c[1], c[1] * c[2], c[0] * c[2]];
    (n as f64 / (N * N * N) as f64, c, std::array::from_fn(|k| sq[k] / n as f64 - centred[k] + if k < 3 { cell } else { 0.0 }))
}

/// The part of the cube of side `s` (about its centre) on the negative side
/// of each plane (value and gradient at the centre), `first` too if any:
/// its share, its centroid's offset and its covariance about it, exactly.
pub(crate) fn planes_part(first: Option<&[f64; 4]>, planes: &[[f64; 4]], s: f64) -> (f64, V3, Cov) {
    sets_part(first, &planes.iter().map(|g| (vec![*g], true)).collect::<Vec<_>>(), s)
}

/// The part of the cube of side `s` (about its centre) inside `first`, if
/// any, and inside or outside each set of planes (inside: on the negative
/// side of them all), as its flag says: its share, its centroid's offset
/// and its covariance about it, exactly. Inside them all is a polytope
/// (`hold::within`); outside one is the rest, by its moments' difference.
/// (By 8^3 points, a strip of a hold 0.01 mm square along a part's edge
/// held no material.)
pub(crate) fn sets_part(first: Option<&[f64; 4]>, sets: &[(Vec<[f64; 4]>, bool)], s: f64) -> (f64, V3, Cov) {
    let inner: Vec<[f64; 4]> = first.into_iter().copied().chain(sets.iter().filter(|x| x.1).flat_map(|x| x.0.iter().copied())).collect();
    let outer: Vec<&Vec<[f64; 4]>> = sets.iter().filter(|x| !x.1).map(|x| &x.0).collect();
    let (mut v, mut m1, mut m2) = (0.0, [0.0; 3], [0.0; 6]);
    for t in 0..1u32 << outer.len() {
        let planes: Vec<[f64; 4]> = inner.iter().copied().chain(outer.iter().enumerate().filter(|(q, _)| t >> q & 1 == 1).flat_map(|(_, o)| o.iter().copied())).collect();
        let (pv, p1, p2) = hold::within(&planes, s);
        let sign = if t.count_ones() % 2 == 0 { 1.0 } else { -1.0 };
        v += sign * pv;
        (0..3).for_each(|k| m1[k] += sign * p1[k]);
        (0..6).for_each(|k| m2[k] += sign * p2[k]);
    }
    let cube = s * s * s;
    if !(v > 1e-12 * cube) {
        return (0.0, [0.0; 3], [0.0; 6]);
    }
    let c = m1.map(|x| x / v);
    let centred = [c[0] * c[0], c[1] * c[1], c[2] * c[2], c[0] * c[1], c[1] * c[2], c[0] * c[2]];
    (v / cube, c, std::array::from_fn(|k| m2[k] / v - centred[k]))
}

/// The covariance (mm^2) of the cube of side `s` cut by the tangent plane
/// (field value and gradient `g` at its centre), material where the field
/// is negative: along the axis the plane is steepest across, each of 8 x 8
/// columns holds an interval exactly, its moments integrated as it lies
/// (so a sheet thinner than the columns' spacing keeps its thickness);
/// across them, at the columns' centres, each a column's width wide.
pub(crate) fn clipped_cov(g: &[f64; 4], s: f64) -> Cov {
    const M: usize = 8;
    let k = (0..3).max_by(|&a, &b| g[1 + a].abs().total_cmp(&g[1 + b].abs())).unwrap_or(2);
    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
    let (gk, half, w) = (g[1 + k], 0.5 * s, s / M as f64);
    // Sums over the material of d, and d_a d_b (as Cov's pairs).
    let (mut vol, mut m1, mut m2) = (0.0, [0.0; 3], [0.0; 6]);
    for c in 0..M * M {
        let mut d = [0.0; 3];
        d[i] = ((c % M) as f64 + 0.5) * w - half;
        d[j] = ((c / M) as f64 + 0.5) * w - half;
        // Along k, a + gk t <= 0 within the cube.
        let a = g[0] + g[1 + i] * d[i] + g[1 + j] * d[j];
        let (lo, hi) = if gk == 0.0 {
            if a <= 0.0 { (-half, half) } else { (0.0, 0.0) }
        } else if gk > 0.0 {
            (-half, (-a / gk).clamp(-half, half))
        } else {
            ((-a / gk).clamp(-half, half), half)
        };
        let len = hi - lo;
        if len <= 0.0 {
            continue;
        }
        let (t1, t2) = ((hi * hi - lo * lo) / 2.0, (hi * hi * hi - lo * lo * lo) / 3.0);
        vol += len;
        m1[i] += len * d[i];
        m1[j] += len * d[j];
        m1[k] += t1;
        let mut sq = [0.0; 3];
        sq[i] = len * (d[i] * d[i] + w * w / 12.0);
        sq[j] = len * (d[j] * d[j] + w * w / 12.0);
        sq[k] = t2;
        let mut cross = [[0.0; 3]; 3];
        cross[i][j] = len * d[i] * d[j];
        cross[i][k] = d[i] * t1;
        cross[j][k] = d[j] * t1;
        let at = |a: usize, b: usize| cross[a.min(b)][a.max(b)] + cross[a.max(b)][a.min(b)];
        let add = [sq[0], sq[1], sq[2], at(0, 1), at(1, 2), at(0, 2)];
        m2.iter_mut().zip(add).for_each(|(m, v)| *m += v);
    }
    if vol <= 0.0 {
        // Under the columns' spacing: a corner's speck, its spread its own.
        return [0.0; 6];
    }
    let c = m1.map(|v| v / vol);
    let centred = [c[0] * c[0], c[1] * c[1], c[2] * c[2], c[0] * c[1], c[1] * c[2], c[0] * c[2]];
    std::array::from_fn(|n| m2[n] / vol - centred[n])
}

/// Where a cut cell's samples join material across air the part does not
/// join it round: a gap glued shut, its width and its middle. The cells'
/// nodes tie together what the cells round them hold, so the
/// discretization joins material within the block of 3^3 cells about a
/// cell; the part joins it only through material. A sample counts as
/// material where any of its sub-cell is (`occupied`), so air under about
/// a sample wide is never parted by the samples. From each partly full
/// sample, one sample along the field's gradient (away from its own wall,
/// into the air): where that reaches material again with the field
/// falling on, the material on the two sides is checked for one piece
/// within the block (`block_pieces`). A slit or a clearance is not; a
/// corner's sides, a groove's under a cell deep, a wedge's near its apex
/// are, and their error is the grid's, which shrinks with the cells (a
/// claim's coarser grids measure it). `t` is the cell's tape, `whole` the
/// shape's (the step may reach the next cell).
pub(crate) fn glued_gap(t: &Tape, whole: &Arc<Tape>, ctx: &Ctx, pts: &[V3], occ: &[(f64, V3, Cov)], h: f64, centre: V3) -> Option<(f64, V3)> {
    let s = h / K as f64;
    let partial: Vec<V3> = pts.iter().zip(occ).filter(|(_, o)| o.0 > 0.0 && o.0 < 1.0).map(|(p, _)| *p).collect();
    if partial.is_empty() {
        return None;
    }
    let steps: Vec<(V3, V3)> = partial.iter().zip(crate::eval::eval_grads(t, ctx, &partial)).filter_map(|(p, g)| v3::direction([g[1], g[2], g[3]]).map(|n| (*p, n))).collect();
    let next: Vec<V3> = steps.iter().map(|(p, n)| std::array::from_fn(|a| p[a] + s * n[a])).collect();
    let ahead: Vec<(V3, V3, V3)> = steps
        .iter()
        .zip(&next)
        .zip(crate::eval::eval_grads(whole, ctx, &next))
        .filter(|((pn, _), g)| {
            let (n, gq) = (pn.1, [g[1], g[2], g[3]]);
            -v3::dot(gq, n) > 0.5 * v3::norm(gq) && crate::query::occupancy(g, [s; 3]) > 0.0
        })
        .map(|((&(p, n), &q), _)| (p, q, n))
        .collect();
    if ahead.is_empty() {
        return None;
    }
    let lo: V3 = std::array::from_fn(|a| centre[a] - 1.5 * h);
    let hi: V3 = std::array::from_fn(|a| centre[a] + 1.5 * h);
    let block = crate::query::focus(whole, ctx, lo, hi);
    let piece = block_pieces(&block, ctx, lo, h);
    let b = 3.0 * h / FINE as f64;
    let label = |q: V3| -> u32 {
        let i: [i64; 3] = std::array::from_fn(|a| ((q[a] - lo[a]) / b).floor() as i64);
        if i.iter().all(|&v| (0..FINE as i64).contains(&v)) { piece[i[0] as usize + FINE * (i[1] as usize + FINE * i[2] as usize)] } else { 0 }
    };
    // The piece of material first met from `from` along `d`, within a sample.
    let seed = |from: V3, d: V3| (0..=(2.0 * s / b).ceil() as usize).map(|k| label(std::array::from_fn(|a| from[a] + 0.5 * b * k as f64 * d[a]))).find(|&l| l != 0);
    let mut out: Option<(f64, V3)> = None;
    for (p, q, n) in ahead {
        // Joined only where bounds show the two sides one piece of solid
        // material (`block_pieces`, flush seams solid by the kernel's
        // rule); a side too thin for them to show solid joins nothing.
        if let (Some(a), Some(c)) = (seed(p, n.map(|v| -v)), seed(q, n))
            && a == c
        {
            continue;
        }
        // The air between: the longest run of it along the step, from a
        // sample behind `p` to one past `q`.
        const ALONG: usize = 192;
        let line: Vec<V3> = (0..ALONG).map(|k| std::array::from_fn(|a| p[a] + s * (3.0 * (k as f64 + 0.5) / ALONG as f64 - 1.0) * n[a])).collect();
        let mut f = Vec::new();
        eval_points(whole, ctx, &line, &mut f);
        let (mut best, mut run) = ((0, 0), 0);
        for (k, v) in f.iter().enumerate() {
            run = if *v > 0.0 { run + 1 } else { 0 };
            if run > best.1 {
                best = (k + 1 - run, run);
            }
        }
        // Material met again along the step by the field's sign, past the air: else the next sample's
        // material was the tangent plane's alone (a field that is no distance there: a bore written as
        // (1 - r)(1 + 0.95 sin 25x) read a 0.215 mm gap it does not have).
        if best.1 > 0 && !f[best.0 + best.1..].iter().any(|v| *v <= 0.0) {
            continue;
        }
        // Not shown joined: air found is a gap; none found proves nothing
        // (a hairline passes between any points), so it is a gap of
        // unknown width, 0 here, under the points' spacing.
        let width = 3.0 * s * best.1 as f64 / ALONG as f64;
        let mid = if best.1 == 0 { std::array::from_fn(|a| 0.5 * (p[a] + q[a])) } else { line[best.0 + best.1 / 2] };
        if out.is_none_or(|o| width < o.0) {
            out = Some((width, mid));
        }
    }
    out
}

/// The part's surface normal at points on or near it: its gradient at
/// each, after Newton's steps onto its zero.
pub(crate) fn surface_normals(part: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<Option<V3>> {
    let mut at = pts.to_vec();
    for _ in 0..3 {
        let g = crate::eval::eval_grads(part, ctx, &at);
        for (p, g) in at.iter_mut().zip(&g) {
            let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
            if n2 > 0.0 {
                *p = std::array::from_fn(|a| p[a] - g[0] * g[1 + a] / n2);
            }
        }
    }
    crate::eval::eval_grads(part, ctx, &at).iter().map(|g| v3::direction([g[1], g[2], g[3]])).collect()
}

/// Boxes per axis of the block of 3^3 cells round a cell (`block_pieces`).
const FINE: usize = 32;

/// Over the block of 3^3 cells from `lo` (cells of h), on a grid of
/// FINE^3 boxes, the piece of material each box lies in: boxes bounds show
/// solid by the kernel's rule (`sign::solid`: inside, or on a 0 where
/// pieces meet flush, which is no surface), joined across their faces,
/// numbered from 1; 0 where a box is not shown so. Air of any width parts
/// the pieces it separates; pieces meeting flush are one; a wall thinner
/// than about two boxes is no piece.
fn block_pieces(tape: &Tape, ctx: &Ctx, lo: V3, h: f64) -> Vec<u32> {
    let b = 3.0 * h / FINE as f64;
    let mut full = vec![false; FINE * FINE * FINE];
    let mut stack = vec![([0usize; 3], FINE)];
    while let Some((i, n)) = stack.pop() {
        let bx: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(lo[a] + i[a] as f64 * b, lo[a] + (i[a] + n) as f64 * b));
        let r = eval_interval(tape, ctx, bx, None);
        // Solid by the kernel's one rule (`sign::solid`): inside by
        // bounds, or inside but on a 0 of no volume, pieces meeting flush
        // (a lattice set on its skin, a union of boxes), read once a box
        // is as small as it gets.
        if r.all_inside() || n == 1 && r.may_be_inside() && ink_core::sign::solid(tape, ctx, bx) {
            for z in i[2]..i[2] + n {
                for y in i[1]..i[1] + n {
                    full[i[0] + FINE * (y + FINE * z)..i[0] + n + FINE * (y + FINE * z)].fill(true);
                }
            }
        } else if r.may_be_inside() && n > 1 {
            let m = n / 2;
            for k in 0..8 {
                stack.push(([i[0] + (k & 1) * m, i[1] + (k >> 1 & 1) * m, i[2] + (k >> 2 & 1) * m], m));
            }
        }
    }
    let mut piece = vec![0u32; full.len()];
    let mut next = 0;
    for start in 0..full.len() {
        if !full[start] || piece[start] != 0 {
            continue;
        }
        next += 1;
        piece[start] = next;
        let mut todo = vec![start];
        while let Some(k) = todo.pop() {
            let at = [k % FINE, (k / FINE) % FINE, k / (FINE * FINE)];
            for a in 0..3 {
                for up in [false, true] {
                    if (up && at[a] + 1 == FINE) || (!up && at[a] == 0) {
                        continue;
                    }
                    let step = [1, FINE, FINE * FINE][a];
                    let j = if up { k + step } else { k - step };
                    if full[j] && piece[j] == 0 {
                        piece[j] = next;
                        todo.push(j);
                    }
                }
            }
        }
    }
    piece
}

/// The cell taking a point `q` of the surface (normal `n`), of `cells`
/// (`at`: a grid cell's; `finest`, their deepest level), and `q` in its
/// coordinates: the cell it lies in, where `takes` it; else the one half
/// that cell's size inside; or, where no cell holds `q` (a sliver dropped,
/// the cell beyond a face on its plane), the one half of each size the
/// cells have inside, finest first. Every read and integral over the
/// surface (stress reads, pressures, films, contacts) finds its cell so: a
/// step of half a grid cell from a split thin fin's face lands in air.
#[allow(clippy::too_many_arguments)]
pub(crate) fn surface_cell<T: Placed>(
    cells: &[T],
    at: &dyn Fn([i64; 3]) -> std::ops::Range<usize>,
    grid: &Grid,
    finest: u8,
    q: V3,
    n: V3,
    takes: &dyn Fn(usize) -> bool,
) -> Option<(usize, V3)> {
    let holding = |p: V3| -> Option<(usize, f64, V3)> {
        let idx: [i64; 3] = std::array::from_fn(|a| ((p[a] - grid.origin[a]) / grid.h).floor() as i64);
        let (ci, xi) = locate(cells, at(idx), grid, p)?;
        let size = place(grid, &cells[ci]).1;
        Some((ci, size, std::array::from_fn(|a| xi[a] + (q[a] - p[a]) * 2.0 / size)))
    };
    let taken = |c: Option<(usize, f64, V3)>| c.filter(|c| takes(c.0)).map(|c| (c.0, c.2));
    let here = holding(q);
    if let Some(c) = taken(here) {
        return Some(c);
    }
    let sizes: Vec<f64> = match here {
        Some(c) => vec![c.1],
        None => (0..=finest).rev().map(|l| grid.h / f64::from(1u32 << l)).collect(),
    };
    sizes.iter().find_map(|&size| taken(holding(std::array::from_fn(|a| q[a] - 0.5 * size * n[a]))))
}

pub(crate) fn center(grid: &Grid, idx: [i64; 3]) -> V3 {
    let (a, b) = grid.cell_bounds(Cell { i: idx, size: 1 });
    [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0]
}

pub(crate) fn at_xi(c: V3, h: f64, xi: V3) -> V3 {
    [c[0] + 0.5 * h * xi[0], c[1] + 0.5 * h * xi[1], c[2] + 0.5 * h * xi[2]]
}

/// With query::occupancy's slab model, the offset from the cube's centre to the
/// centroid of its material part: half the empty part's depth, against the
/// gradient. Lumping a partly filled sub-cell at its centre misplaces the
/// boundary by up to a quarter of the sub-cell (first order in stiffness).
pub(crate) fn occupied_centroid(f: &[f64; 4], s: f64) -> V3 {
    let g2 = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
    let w = s * (f[1].abs() + f[2].abs() + f[3].abs());
    if !(g2 > 0.0 && w.is_finite() && f[0].is_finite()) {
        return [0.0; 3];
    }
    // Along the gradient (field units), material spans -w/2 .. min(-f0, w/2).
    // The slab is deeper than the cube along a diagonal: keep the centroid
    // inside the cube.
    let top = (-f[0]).clamp(-0.5 * w, 0.5 * w);
    let c = 0.5 * (top - 0.5 * w);
    [f[1], f[2], f[3]].map(|g| (g * c / g2).clamp(-0.5 * s, 0.5 * s))
}
