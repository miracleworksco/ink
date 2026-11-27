use super::analysis::find;
use super::brick::B;
use super::report::Ran;
use crate::octree::{Cell, Grid};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_points};
use ink_core::interval::I;
use ink_core::tape::{Pass, Tape};
use std::collections::HashMap;
use std::sync::Arc;

/// Whether what a mesh draws besides its shape (its paint) is the same
/// throughout a box, so a block of it may be drawn whole (`Flat`).
pub type Even<'a> = dyn Fn([I; 3]) -> bool + Sync + 'a;

/// The largest block drawn whole (`Flat`), in cells: the least is a brick.
pub(super) const FLAT_MOST: i64 = 64 * B;

/// A slab of a block of cells drawn whole: the block's cells `lo..hi` along
/// `axis` (grid indices). Over it, grown by half a cell across `axis`
/// (where its rim's vertices may lie), the field's rate along `axis` is
/// shown never 0 (at least `rate` across: it changes sign at most once
/// along each line that way), and on its two faces square to `axis` one
/// sign each, apart: the surface there is a height field over the square
/// across it, one sheet crossing each line along `axis` once, between the
/// faces, and nothing else lies there. `slopes`: the rates along the next
/// two axes over the rate along `axis` there, bounding the height field's
/// slopes. `tape`: the field there; `way`: the sign of the rate along
/// `axis`; `span`: where its faces lie along `axis` (a cell thick, the
/// field may keep one sign between one and its cells' side). The block's other slabs
/// along `axis` are drawn whole too, hold no surface, or are a cell thick
/// and meshed as the grid is (`plain`: a face of a wall too thin for a
/// plane of the grid between it and the wall's middle) (`slabs`). Slabs beside each
/// other with one axis and way are joined in a region (`regions`), meshed
/// only round its rim (`rim_cell`) and triangulated inside it (`disc`).
pub(super) struct Flat {
    pub(super) cell: Cell,
    pub(super) axis: usize,
    pub(super) lo: i64,
    pub(super) hi: i64,
    span: (f64, f64),
    pub(super) plain: bool,
    slopes: [I; 2],
    rate: f64,
    way: f64,
    pub(super) tape: Arc<Tape>,
}

/// Whether cell `c` within block `q` touches one of its sides along other
/// axes than `a`: holds a cell of its rim.
pub(super) fn rim(c: Cell, q: Cell, a: usize) -> bool {
    (0..3).filter(|&k| k != a).any(|k| c.i[k] == q.i[k] || c.i[k] + c.size == q.i[k] + q.size)
}

/// Block `c` as one drawn whole, where it is one: along the axis the
/// field's rate may be largest along over it, its slabs (`Flat`) each
/// drawn whole, holding no surface, or a cell thick and meshed as the grid
/// is (no more than PLAIN_SLABS of those), halving each that is none of
/// these while its field is not monotone along the axis; and one drawn
/// whole. Not tried where the rates across the axis spread over the block
/// more than a slab drawn whole's slopes may, at the largest rate along it:
/// an edge's other face, a hole's wall, a curve, where a slab's spread as
/// wide; but also a wall's face whose block reaches into a hollow behind
/// it as far as where another wall is nearer, which a slab could draw:
/// trying those too cost a swerve base's meshes 7% more than it saved.
pub(super) fn flat(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, c: Cell, within: f64, even: &Even<'_>, ran: &Ran) -> Option<Vec<Flat>> {
    let h = grid.h;
    let (lo, hi) = (grid.pos(c.i), grid.pos(c.i.map(|x| x + c.size)));
    ran.boxes.fetch_add(4 * tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
    let g = ink_core::eval::eval_interval_grad(tape, ctx, std::array::from_fn(|k| I::new(lo[k] - 0.5 * h, hi[k] + 0.5 * h))).d;
    let a = (0..3).max_by(|&x, &y| ink_core::interval::mag(g[x]).total_cmp(&ink_core::interval::mag(g[y])))?;
    let across = [1, 2].map(|k| ink_core::interval::div_o(g[(a + k) % 3], I::point(ink_core::interval::mag(g[a]))));
    if !(spread(across) * 2f64.sqrt() * (hi[a] - lo[a] + h) <= 0.5 * within) {
        return None;
    }
    let mut out = Vec::new();
    (slabs(tape, ctx, grid, c, a, c.i[a], c.i[a] + c.size, within, even, ran, &mut out) && out.iter().any(|f| !f.plain)).then_some(out)
}

/// The most slabs of a block drawn whole that are meshed as the grid is:
/// a face of a wall too thin for a plane of the grid between it and the
/// wall's middle, each side.
const PLAIN_SLABS: usize = 2;

/// Slab `lo..hi` along `a` of block `c`, added to `out` where it is drawn
/// whole (`slab`; a cell thick, on one side of a plane through it where
/// the field keeps one sign on the other), or is a cell thick (meshed as
/// the grid is, while no more than PLAIN_SLABS are); else true where the
/// field over it (grown by half a cell across `a`) keeps one sign, or its
/// halves along `a` are each one of these.
#[allow(clippy::too_many_arguments)]
fn slabs(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, c: Cell, a: usize, lo: i64, hi: i64, within: f64, even: &Even<'_>, ran: &Ran, out: &mut Vec<Flat>) -> bool {
    let h = grid.h;
    let (clo, chi) = (grid.pos(c.i), grid.pos(c.i.map(|x| x + c.size)));
    let mut b: [I; 3] = std::array::from_fn(|k| I::new(clo[k] - 0.5 * h, chi[k] + 0.5 * h));
    b[a] = I::new(grid.pos([lo; 3])[a], grid.pos([hi; 3])[a]);
    ran.boxes.fetch_add(tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
    let r = eval_interval(tape, ctx, b, None);
    if r.clear_of_surface() || r.all_inside() {
        return true;
    }
    match slab(tape, ctx, c, a, (lo, hi), b, None, within, even, ran) {
        Slab::Flat(f) => {
            out.push(f);
            return true;
        }
        Slab::Empty => return true,
        Slab::Stop => return false,
        Slab::Halve => {}
    }
    if hi - lo == 1 {
        // A cell thick: drawn whole on one side of a plane through it (at
        // a half, a quarter or three quarters) where the field keeps one
        // sign on the other (a wall's middle, where its faces' pieces meet),
        // its proof reaching to the cell's side, which may read 0 there
        // where the quarter cell past it holds no material (a face on a
        // plane of the grid: its samples read 0, drawn as air, so it is the
        // cells' on its material's side). Never past the side: the cells
        // there are meshed apart, and a sheet in them drawn twice.
        let (s0, s1) = (b[a].lo, b[a].hi);
        let q = 0.25 * (s1 - s0);
        for t in [0.5, 0.25, 0.75] {
            let at = s0 + t * (s1 - s0);
            // Past the cells' side above, air is there where the field rises
            // along the axis; below, where it falls.
            for (one, rest, beyond, past) in [((s0, at), (at, s1), (s1, s1 + q), 1.0), ((at, s1), (s0, at), (s0 - q, s0), -1.0)] {
                let part = |r: (f64, f64)| {
                    let mut p = b;
                    p[a] = I::new(r.0, r.1);
                    p
                };
                ran.boxes.fetch_add(tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
                let r = eval_interval(tape, ctx, part(one), None);
                if !(r.all_inside() || r.clear_of_surface()) {
                    continue;
                }
                if let Slab::Flat(f) = slab(tape, ctx, c, a, (lo, hi), part(rest), None, within, even, ran) {
                    out.push(f);
                    return true;
                }
                ran.boxes.fetch_add(tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
                if eval_interval(tape, ctx, part(beyond), None).none_inside()
                    && let Slab::Flat(f) = slab(tape, ctx, c, a, (lo, hi), part(rest), Some(past), within, even, ran)
                    && f.way == past
                {
                    out.push(f);
                    return true;
                }
            }
        }
        let plain = Flat { cell: c, axis: a, lo, hi, span: (b[a].lo, b[a].hi), plain: true, slopes: [I::point(0.0); 2], rate: 0.0, way: 0.0, tape: tape.clone() };
        out.push(plain);
        return out.iter().filter(|f| f.plain).count() <= PLAIN_SLABS;
    }
    let mid = lo + (hi - lo) / 2;
    slabs(tape, ctx, grid, c, a, lo, mid, within, even, ran, out) && slabs(tape, ctx, grid, c, a, mid, hi, within, even, ran, out)
}

/// Slab `lo..hi` along `a` of block `c` (`b` its box, grown by half a cell
/// across `a`) as one drawn whole (`Flat`), where it is one: the surface
/// flat enough that a triangle across it (no longer than its square's
/// diagonal, grown) strays no more than half of `within`, and `even`
/// holding over it. Its halves are worth trying where its field is not
/// shown monotone along `a`, or its faces one sign each. `zero`: its face
/// at that end (1, above; -1, below) may read 0 as air. `ran` counts the
/// instructions run.
#[allow(clippy::too_many_arguments)]
fn slab(tape: &Arc<Tape>, ctx: &Ctx, c: Cell, a: usize, (lo, hi): (i64, i64), b: [I; 3], zero: Option<f64>, within: f64, even: &Even<'_>, ran: &Ran) -> Slab {
    ran.boxes.fetch_add(4 * tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
    let g = ink_core::eval::eval_interval_grad(tape, ctx, b).d;
    let rate = if g[a].lo > 0.0 {
        g[a].lo
    } else if g[a].hi < 0.0 {
        -g[a].hi
    } else {
        0.0
    };
    if !(rate > 0.0 && rate.is_finite()) {
        return Slab::Halve;
    }
    // Its slopes range over its whole square: a thinner slab's do too.
    let slopes = [1, 2].map(|k| ink_core::interval::div_o(g[(a + k) % 3], g[a]));
    let side = b[(a + 1) % 3].hi - b[(a + 1) % 3].lo;
    if !(spread(slopes) * 2f64.sqrt() * side <= 0.5 * within) {
        return Slab::Stop;
    }
    let face = |x: f64| {
        let mut f = b;
        f[a] = I::point(x);
        eval_interval(tape, ctx, f, None)
    };
    ran.boxes.fetch_add(2 * tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed);
    let (below, above) = (face(b[a].lo), face(b[a].hi));
    let side = |r: I, end: f64| {
        if r.all_inside() {
            -1
        } else if r.clear_of_surface() || zero == Some(end) && r.none_inside() {
            1
        } else {
            0
        }
    };
    match (side(below, -1.0), side(above, 1.0)) {
        (0, _) | (_, 0) => Slab::Halve,
        // Monotone, one sign at both ends: no sign change across it.
        (x, y) if x == y => Slab::Empty,
        _ if !even(b) => Slab::Stop,
        _ => Slab::Flat(Flat { cell: c, axis: a, lo, hi, span: (b[a].lo, b[a].hi), plain: false, slopes, rate, way: g[a].lo.signum(), tape: tape.clone() }),
    }
}

/// What a slab is (`slab`): drawn whole; holding no surface; not shown
/// either, but its halves may be; or not, nor its halves.
enum Slab {
    Flat(Flat),
    Empty,
    Halve,
    Stop,
}

/// How far a height field with slopes in `s` strays from a triangle with
/// its corners on it, per unit of the triangle's longest side (seen along
/// its axis): half the diagonal of the slopes' box. A point of the
/// triangle is the corners' weighted mean; the field there differs from
/// theirs by its slopes along the way to each, less one slope's (which
/// the weights cancel), by at most this times the way's length.
fn spread(s: [I; 2]) -> f64 {
    0.5 * (s[0].hi - s[0].lo).hypot(s[1].hi - s[1].lo)
}

/// The blocks drawn whole, each's slabs by its cell (`at`, indices in the
/// slabs), the sizes they come in (`sizes`), with each slab's region
/// (`region`, none for one sampled whole again).
pub(super) struct Blocks {
    pub(super) at: HashMap<([i64; 3], i64), Vec<u32>>,
    sizes: Vec<i64>,
    pub(super) region: Vec<Option<u32>>,
}

impl Blocks {
    pub(super) fn new(flats: &[Flat]) -> Blocks {
        let mut at: HashMap<([i64; 3], i64), Vec<u32>> = HashMap::new();
        for (k, f) in flats.iter().enumerate() {
            at.entry((f.cell.i, f.cell.size)).or_default().push(k as u32);
        }
        let mut sizes: Vec<i64> = flats.iter().map(|f| f.cell.size).collect();
        sizes.sort_unstable();
        sizes.dedup();
        Blocks { at, sizes, region: vec![None; flats.len()] }
    }

    /// The slabs of the block drawn whole holding cell `c`.
    pub(super) fn block(&self, c: [i64; 3]) -> &[u32] {
        self.sizes.iter().find_map(|&s| self.at.get(&(c.map(|x| x.div_euclid(s) * s), s))).map_or(&[], Vec::as_slice)
    }

    /// The slab holding cell `c`.
    pub(super) fn of(&self, flats: &[Flat], c: [i64; 3]) -> Option<u32> {
        self.block(c).iter().copied().find(|&k| {
            let f = &flats[k as usize];
            (f.lo..f.hi).contains(&c[f.axis])
        })
    }

    pub(super) fn region(&self, flats: &[Flat], c: [i64; 3]) -> Option<u32> {
        self.of(flats, c).and_then(|k| self.region[k as usize])
    }
}

/// The slabs drawn whole joined into regions: each joined to those
/// beside it on its sides (along the other axes than its own) with its
/// axis and way, overlapping it along the axis, while the region's slopes (Flat::slopes, their hull)
/// over the diagonal of its box seen along the axis (as long as a
/// triangle across it may be) stray no more than half of `within`
/// (`spread`). Each region's slopes.
pub(super) fn regions(flats: &[Flat], blocks: &mut Blocks, h: f64, within: f64) -> Vec<[I; 2]> {
    let n = flats.len();
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for (k, f) in flats.iter().enumerate().filter(|f| !f.1.plain) {
        for side in 0..4 {
            let axis = (f.axis + 1 + side / 2) % 3;
            let mut c = f.cell.i;
            c[axis] += if side % 2 == 0 { -1 } else { f.cell.size };
            for &j in blocks.block(c) {
                let g = &flats[j as usize];
                if !g.plain && g.axis == f.axis && g.way == f.way && g.lo < f.hi && f.lo < g.hi {
                    pairs.push(((k as u32).min(j), (k as u32).max(j)));
                }
            }
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    // Each region's slopes and box across its axis (cells), by its root.
    let uv = |f: &Flat| ((f.axis + 1) % 3, (f.axis + 2) % 3);
    let mut slopes: Vec<[I; 2]> = flats.iter().map(|f| f.slopes).collect();
    let mut span: Vec<[i64; 4]> = flats
        .iter()
        .map(|f| {
            let (u, v) = uv(f);
            [f.cell.i[u], f.cell.i[v], f.cell.i[u] + f.cell.size, f.cell.i[v] + f.cell.size]
        })
        .collect();
    for (a, b) in pairs {
        let (ra, rb) = (find(&mut parent, a) as usize, find(&mut parent, b) as usize);
        if ra == rb {
            continue;
        }
        let s: [I; 2] = std::array::from_fn(|k| I::new(slopes[ra][k].lo.min(slopes[rb][k].lo), slopes[ra][k].hi.max(slopes[rb][k].hi)));
        let (x, y) = (span[ra], span[rb]);
        let box_: [i64; 4] = [x[0].min(y[0]), x[1].min(y[1]), x[2].max(y[2]), x[3].max(y[3])];
        let diagonal = ((box_[2] - box_[0] + 1) as f64 * h).hypot((box_[3] - box_[1] + 1) as f64 * h);
        if spread(s) * diagonal <= 0.5 * within {
            let (keep, gone) = (ra.min(rb), ra.max(rb));
            parent[gone] = keep as u32;
            (slopes[keep], span[keep]) = (s, box_);
        }
    }
    let mut index: HashMap<u32, u32> = HashMap::new();
    let mut out = Vec::new();
    for k in (0..n as u32).filter(|&k| !flats[k as usize].plain) {
        let r = find(&mut parent, k);
        let next = index.len() as u32;
        let i = *index.entry(r).or_insert_with(|| {
            out.push(slopes[r as usize]);
            next
        });
        blocks.region[k as usize] = Some(i);
    }
    out
}

/// Whether cell `c` of block `q` (axis `a`) is meshed: in a slab's region,
/// on the block's rim, beside a cell (along the other axes, sides or corners) not in
/// the region. Of a quad round an edge of the grid, its cells lie at one
/// height along the axis (no edge near the block's faces square to it
/// changes sign), so a quad holding a cell outside the region holds only
/// cells meshed, or outside.
pub(super) fn rim_cell(flats: &[Flat], blocks: &Blocks, q: Cell, a: usize, c: [i64; 3]) -> bool {
    let (u, v) = ((a + 1) % 3, (a + 2) % 3);
    if !rim(Cell { i: c, size: 1 }, q, a) {
        return false;
    }
    let Some(r) = blocks.region(flats, c) else { return false };
    // A cell beside it in the block lies at its height, in its slab.
    let held = |n: [i64; 3]| (0..3).all(|k| q.i[k] <= n[k] && n[k] < q.i[k] + q.size);
    (0..9).any(|k| {
        let mut n = c;
        n[u] += k % 3 - 1;
        n[v] += k / 3 - 1;
        k != 4 && !held(n) && blocks.region(flats, n) != Some(r)
    })
}

/// The triangles drawing a region whole (`regions`), from its rim's edges
/// (`kept`, each as the quads outside it run it): those the quads do not
/// also run back make the loops its cells inside close, one round it and
/// one round each hole, triangulated (`bridge`, `ear_clip`) as seen along
/// its axis; and how far the triangles may stray from the surface: its
/// slopes' `spread` times the longest side seen so, and the farthest the
/// loops' vertices lie from the surface along the axis (each's value over
/// its block's rate: the line along the axis through it crosses the
/// surface in the block, its value changing at least that fast on the
/// way). None where the loops are not simple and apart, or more than one
/// turns the way the surface faces, or the triangles may stray more than
/// half of `within`. `cell`: each vertex's cell; `ran` counts the
/// instructions run.
#[allow(clippy::too_many_arguments)]
pub(super) fn disc(
    flats: &[Flat],
    blocks: &Blocks,
    slopes: [I; 2],
    ctx: &Ctx,
    verts: &[V3],
    cell: &[[i64; 3]],
    kept: &[[u32; 2]],
    within: f64,
    ran: &Ran,
) -> Option<(Vec<[u32; 3]>, f64)> {
    let have: std::collections::HashSet<[u32; 2]> = kept.iter().copied().collect();
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &[a, b] in kept {
        if !have.contains(&[b, a]) && next.insert(b, a).is_some() {
            return None;
        }
    }
    let first = &flats[blocks.of(flats, cell[*next.keys().next()? as usize])? as usize];
    let (axis, way) = (first.axis, first.way);
    // The loops, each from its least vertex.
    let mut starts: Vec<u32> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut loops: Vec<Vec<u32>> = Vec::new();
    for s in starts {
        if seen.contains(&s) {
            continue;
        }
        let mut ring = vec![s];
        seen.insert(s);
        let mut x = *next.get(&s)?;
        while x != s {
            if !seen.insert(x) {
                return None;
            }
            ring.push(x);
            x = *next.get(&x)?;
        }
        if ring.len() < 3 {
            return None;
        }
        loops.push(ring);
    }
    let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
    let ids: Vec<u32> = loops.iter().flatten().copied().collect();
    let place: HashMap<u32, usize> = ids.iter().enumerate().map(|(k, &i)| (i, k)).collect();
    let p: Vec<[f64; 2]> = ids.iter().map(|&i| [verts[i as usize][u], verts[i as usize][v]]).collect();
    let rings: Vec<Vec<usize>> = loops.iter().map(|l| l.iter().map(|i| place[i]).collect()).collect();
    if !simple(&p, &rings) {
        return None;
    }
    let area = |r: &[usize]| (0..r.len()).map(|k| cross2(p[r[k]], p[r[(k + 1) % r.len()]])).sum::<f64>() * way;
    let (outer, holes): (Vec<Vec<usize>>, Vec<Vec<usize>>) = rings.into_iter().partition(|r| area(r) > 0.0);
    let [outer] = <[Vec<usize>; 1]>::try_from(outer).ok()?;
    let joined = bridge(&p, outer, holes)?;
    let tris = ear_clip(&p, &joined, way)?;
    // How far each vertex lies from the surface along the axis.
    let mut off = 0.0f64;
    let mut by: HashMap<u32, Vec<usize>> = HashMap::new();
    for (k, &i) in ids.iter().enumerate() {
        by.entry(blocks.of(flats, cell[i as usize])?).or_default().push(k);
    }
    for (b, ks) in by {
        let f = &flats[b as usize];
        let at: Vec<V3> = ks.iter().map(|&k| verts[ids[k] as usize]).collect();
        // Each in its slab, between its faces.
        let (lo, hi) = f.span;
        if at.iter().any(|p| !(lo <= p[axis] && p[axis] <= hi)) {
            return None;
        }
        let mut vals = Vec::new();
        eval_points(&f.tape, ctx, &at, &mut vals);
        ran.points.fetch_add(at.len() * f.tape.work(Pass::Point), std::sync::atomic::Ordering::Relaxed);
        off = vals.iter().fold(off, |m, x| m.max(x.abs() / f.rate));
    }
    let side = |a: usize, b: usize| (p[a][0] - p[b][0]).hypot(p[a][1] - p[b][1]);
    let longest = tris.iter().map(|t| side(t[0], t[1]).max(side(t[1], t[2])).max(side(t[2], t[0]))).fold(0.0f64, f64::max);
    let stray = spread(slopes) * longest + off;
    (stray <= 0.5 * within).then(|| (tris.iter().map(|t| t.map(|k| ids[k])).collect(), stray))
}

fn cross2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}

/// Twice the signed area of triangle a b c.
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Whether closed polygons `rings` (indices into `p`) are simple and apart:
/// no two of their sides meet but neighbours at their shared corner (sides
/// found near each other by squares of the longest side's length).
fn simple(p: &[[f64; 2]], rings: &[Vec<usize>]) -> bool {
    // Each side: its ring, its place in it, and its ends.
    let sides: Vec<(usize, usize, [f64; 2], [f64; 2])> =
        rings.iter().enumerate().flat_map(|(r, ring)| (0..ring.len()).map(move |k| (r, k, p[ring[k]], p[ring[(k + 1) % ring.len()]]))).collect();
    let s = sides.iter().map(|(_, _, a, b)| (a[0] - b[0]).abs().max((a[1] - b[1]).abs())).fold(0.0f64, f64::max);
    if !(s > 0.0 && s.is_finite()) {
        return false;
    }
    let key = |q: [f64; 2]| [(q[0] / s).floor() as i64, (q[1] / s).floor() as i64];
    let mut near: HashMap<[i64; 2], Vec<usize>> = HashMap::new();
    for (k, (_, _, a, b)) in sides.iter().enumerate() {
        let (a, b) = (key(*a), key(*b));
        for x in a[0].min(b[0])..=a[0].max(b[0]) {
            for y in a[1].min(b[1])..=a[1].max(b[1]) {
                near.entry([x, y]).or_default().push(k);
            }
        }
    }
    let meet = |(a, b): ([f64; 2], [f64; 2]), (c, d): ([f64; 2], [f64; 2])| {
        let (o1, o2, o3, o4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
        if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
            return true;
        }
        // Touching or overlapping in a line: a corner on the other side.
        let on = |o: f64, a: [f64; 2], b: [f64; 2], q: [f64; 2]| o == 0.0 && q[0] >= a[0].min(b[0]) && q[0] <= a[0].max(b[0]) && q[1] >= a[1].min(b[1]) && q[1] <= a[1].max(b[1]);
        on(o1, a, b, c) || on(o2, a, b, d) || on(o3, c, d, a) || on(o4, c, d, b)
    };
    near.values().all(|ks| {
        ks.iter().enumerate().all(|(x, &i)| {
            ks[x + 1..].iter().all(|&j| {
                let ((ri, ki, a, b), (rj, kj, c, d)) = (sides[i], sides[j]);
                let n = rings[ri].len();
                if ri == rj && (ki + 1) % n == kj {
                    // Neighbours meet at their shared corner only: not
                    // folded back along one line.
                    orient(a, b, d) != 0.0 || (b[0] - a[0]) * (d[0] - b[0]) + (b[1] - a[1]) * (d[1] - b[1]) > 0.0
                } else if ri == rj && (kj + 1) % n == ki {
                    orient(c, d, b) != 0.0 || (d[0] - c[0]) * (b[0] - d[0]) + (d[1] - c[1]) * (b[1] - d[1]) > 0.0
                } else {
                    !meet((a, b), (c, d))
                }
            })
        })
    })
}

/// Polygon `outer` with its `holes` (indices into `p`; the holes turning
/// the other way, inside it, all apart) made one polygon: each hole, the
/// farthest along x first, joined by a bridge there and back from its
/// corner farthest along x to a corner of the polygon it sees (Eberly,
/// Triangulation by Ear Clipping, 2002): of the side the ray from it along
/// x first meets, the end farther along x, or of the polygon's corners in
/// the triangle that corner, the ray's point and the hole's make, the one
/// at the least angle to the ray. None where the ray meets no side.
fn bridge(p: &[[f64; 2]], mut poly: Vec<usize>, mut holes: Vec<Vec<usize>>) -> Option<Vec<usize>> {
    let far = |r: &[usize]| (0..r.len()).max_by(|&a, &b| p[r[a]][0].total_cmp(&p[r[b]][0]).then(p[r[b]][1].total_cmp(&p[r[a]][1]))).unwrap();
    holes.sort_by(|a, b| p[b[far(b)]][0].total_cmp(&p[a[far(a)]][0]));
    for hole in holes {
        let m = far(&hole);
        let mp = p[hole[m]];
        let n = poly.len();
        let mut best: Option<(f64, usize)> = None;
        for i in 0..n {
            let (a, b) = (p[poly[i]], p[poly[(i + 1) % n]]);
            if (a[1] - mp[1]) * (b[1] - mp[1]) > 0.0 || a[1] == b[1] {
                continue;
            }
            let x = a[0] + (mp[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if x >= mp[0] && best.is_none_or(|(bx, _)| x < bx) {
                best = Some((x, i));
            }
        }
        let (x, i) = best?;
        let ray = [x, mp[1]];
        let mut pick = if p[poly[i]][0] >= p[poly[(i + 1) % n]][0] { i } else { (i + 1) % n };
        if p[poly[pick]] != ray {
            let tri = [mp, ray, p[poly[pick]]];
            let sign = orient(tri[0], tri[1], tri[2]).signum();
            let angle = |q: [f64; 2]| ((q[1] - mp[1]).abs()).atan2(q[0] - mp[0]);
            let mut least = angle(p[poly[pick]]);
            for (k, &q) in poly.iter().enumerate() {
                let q = p[q];
                let inside = (0..3).all(|e| orient(tri[e], tri[(e + 1) % 3], q) * sign >= 0.0);
                if k != pick && inside && q != mp && angle(q) < least {
                    (pick, least) = (k, angle(q));
                }
            }
        }
        let mut joined = Vec::with_capacity(n + hole.len() + 2);
        joined.extend_from_slice(&poly[..=pick]);
        joined.extend(hole[m..].iter().chain(&hole[..=m]));
        joined.extend_from_slice(&poly[pick..]);
        poly = joined;
    }
    Some(poly)
}

/// Polygon `poly` (indices into `p`, turning the way of `way`'s sign;
/// simple, but where a bridge to a hole runs there and back) cut into
/// triangles by its ears: corners turning its way whose triangle holds no
/// other corner (but one where its own corners are). Each triangle's
/// corners, as indices in `p`, in its order. None where no ear is left
/// before the end.
fn ear_clip(p: &[[f64; 2]], poly: &[usize], way: f64) -> Option<Vec<[usize; 3]>> {
    let n = poly.len();
    let at = |k: usize| p[poly[k]];
    let (mut prev, mut next): (Vec<usize>, Vec<usize>) = ((0..n).map(|k| (k + n - 1) % n).collect(), (0..n).map(|k| (k + 1) % n).collect());
    let convex = |a: usize, b: usize, c: usize| orient(at(a), at(b), at(c)) * way > 0.0;
    let mut reflex: std::collections::BTreeSet<usize> = (0..n).filter(|&k| !convex(prev[k], k, next[k])).collect();
    let ear = |prev: &[usize], next: &[usize], reflex: &std::collections::BTreeSet<usize>, k: usize| {
        let (a, c) = (prev[k], next[k]);
        let corners = [at(a), at(k), at(c)];
        convex(a, k, c)
            && !reflex.iter().any(|&q| {
                let pq = at(q);
                !corners.contains(&pq) && {
                    let (o1, o2, o3) = (orient(corners[0], corners[1], pq) * way, orient(corners[1], corners[2], pq) * way, orient(corners[2], corners[0], pq) * way);
                    o1 >= 0.0 && o2 >= 0.0 && o3 >= 0.0
                }
            })
    };
    let mut out = Vec::with_capacity(n.saturating_sub(2));
    let (mut left, mut k, mut tried) = (n, 0, 0);
    while left > 3 {
        if tried > left {
            return None;
        }
        if ear(&prev, &next, &reflex, k) {
            let (a, c) = (prev[k], next[k]);
            out.push([poly[a], poly[k], poly[c]]);
            next[a] = c;
            prev[c] = a;
            reflex.remove(&k);
            left -= 1;
            for q in [a, c] {
                if convex(prev[q], q, next[q]) {
                    reflex.remove(&q);
                } else {
                    reflex.insert(q);
                }
            }
            tried = 0;
            k = a;
        } else {
            tried += 1;
            k = next[k];
        }
    }
    if !convex(prev[k], k, next[k]) {
        return None;
    }
    out.push([poly[prev[k]], poly[k], poly[next[k]]]);
    Some(out)
}
