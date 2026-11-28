//! Material between a cell's samples.
//!
//! A cell knows the part by its 4^3 samples, so a foil, web or ligament
//! thinner than their spacing can lie between them unseen. Each plane
//! through the samples, across each axis, is walked between them (and from
//! the outer samples to the cell's faces): from a point outside reading f,
//! no material lies within f / L of it, L a bound on the field's slope
//! there (its interval gradient: a field need not be a distance, the cube
//! of one reads 8 at 2 mm, the root of one 2 at 4). A plane, not the lines
//! between the samples: a rod a quarter of a cell thick between those
//! lines crossed none (a sample sees material in its box only where the
//! field's tangent plane at its centre reaches it, as a distance's does),
//! and its stiffness was lost. At the floor (a 32nd of the cell, the
//! smallest it splits to) the field is read at 4 points a side, and what
//! the bounds still cannot settle is taken as found. (At an eighth of the
//! cell, a field whose bounds are loose, a lofted blade's, left hundreds
//! of cells along its edges taken as found, by their bounds and not by any
//! material.)
//!
//! The planes are walked through every sample's box that sees no
//! material, in a cell whose samples see some as in one whose samples see
//! none. Beside a box that sees material, the floor is an eighth of the
//! cell, its points an eighth of a sample apart, and material read at one
//! is that material's own face (a corner's tip, a face grazed by a plane)
//! where it lies within an eighth of a sample past the tangent plane of
//! every seen sample beside it (the plane the cells place their material
//! by), or where a straight line from it to a seen sample's centre inside
//! the part crosses no air; else it is material the samples miss beside
//! material they see. What no point reads there and the bounds cannot
//! settle is the face's. (Taken as the face's own a sample's spacing out,
//! rods a sample under a face were never looked for; nor, in a cell whose
//! samples saw material, a rod in the same cell as the face.) Material in
//! a box that sees some is that sample's own, and is not looked for.
//!
//! Material found is sorted by what it joins (`changes_connection`): a tip
//! of seen material, touching it at one place, changes nothing; a foil
//! joining nothing seen, two pieces apart, one piece at two places (a
//! ligament beside it), or running beside it over more than a tip's few
//! cells makes every claim unknown, naming the place.

use crate::eval::{Ctx, eval_interval, eval_points};
use crate::occupancy::tangent_offsets;
use crate::tape::Tape;
use ink_core::{V3, eval, interval, v3};
use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

/// The samples' lattice over the whole grid: sample `g` (each cell's 4 a
/// side, `g = 4 cell + j`) at `o + (g + 1/2) s`, its box the lattice cube
/// round it; which see material (read off the cells' own samples, no
/// evaluation); and the whole shape, whose field a seen sample's tangent
/// plane is read from where a walk needs it.
pub(crate) struct Lattice<'a> {
    pub o: V3,
    pub s: f64,
    /// A cell's samples that see material, as 64 bits (bit a + 4 b + 16 c):
    /// all of them for a full cell, none for one the samples see nothing in.
    pub cell: &'a (dyn Fn([i64; 3]) -> u64 + Sync),
    pub shape: &'a Tape,
    /// The seen samples' tangent planes read so far, shared by every
    /// cell's walk.
    pub planes: Planes,
}

/// Seen samples' tangent planes (`Near::plane`), in shards locked one at
/// a time.
pub(crate) struct Planes(Vec<Mutex<HashMap<[i64; 3], Option<[f64; 4]>>>>);

impl Planes {
    pub(crate) fn new() -> Planes {
        Planes((0..64).map(|_| Mutex::new(HashMap::new())).collect())
    }
    fn shard(&self, k: [i64; 3]) -> &Mutex<HashMap<[i64; 3], Option<[f64; 4]>>> {
        let h = (k[0].wrapping_mul(73_856_093) ^ k[1].wrapping_mul(19_349_663) ^ k[2].wrapping_mul(83_492_791)) as u64;
        &self.0[(h % 64) as usize]
    }
}

/// How far past the tangent planes of the seen samples beside it material
/// in a box that sees none is still their face's own, in samples: an
/// eighth (a 32nd of the cell, the least share of a sample the cells
/// count).
const OWN: f64 = 0.125;

/// Seen samples within a cell of cell `i` (its own and its 26 neighbours'),
/// by lattice offset from the cell's first sample less 4 (-4 .. 8 each
/// way): read once per cell walked; and by the same offset, the planes and
/// the boxes of air this walk has read.
struct Near {
    base: [i64; 3],
    bits: Vec<bool>,
    planes: RefCell<HashMap<usize, Option<[f64; 4]>>>,
    air: RefCell<HashMap<usize, bool>>,
}

impl Near {
    fn of(lat: &Lattice, i: [i64; 3]) -> Near {
        let mut bits = vec![false; 12 * 12 * 12];
        for n in 0..27 {
            let d = [n % 3 - 1, n / 3 % 3 - 1, n / 9 - 1];
            let m = (lat.cell)([i[0] + d[0], i[1] + d[1], i[2] + d[2]]);
            if m == 0 {
                continue;
            }
            for k in 0..64 {
                if m >> k & 1 == 1 {
                    let g = [(d[0] + 1) * 4 + (k % 4) as i64, (d[1] + 1) * 4 + (k / 4 % 4) as i64, (d[2] + 1) * 4 + (k / 16) as i64];
                    bits[(g[0] + 12 * (g[1] + 12 * g[2])) as usize] = true;
                }
            }
        }
        Near { base: i.map(|v| 4 * v - 4), bits, planes: RefCell::new(HashMap::new()), air: RefCell::new(HashMap::new()) }
    }
    /// Lattice sample `g`'s place in this walk's tables, if within them.
    fn index(&self, g: [i64; 3]) -> Option<usize> {
        let l: [i64; 3] = std::array::from_fn(|a| g[a] - self.base[a]);
        l.iter().all(|&v| (0..12).contains(&v)).then(|| (l[0] + 12 * (l[1] + 12 * l[2])) as usize)
    }
    /// Whether lattice sample `g` (global) sees material (outside the 27
    /// cells: taken as false).
    fn seen(&self, g: [i64; 3]) -> bool {
        self.index(g).is_some_and(|i| self.bits[i])
    }
    /// The seen samples beside lattice sample `g` (its 26 nearest).
    fn beside(&self, g: [i64; 3]) -> impl Iterator<Item = [i64; 3]> + '_ {
        (0..27).map(move |n| [g[0] + n % 3 - 1, g[1] + n / 3 % 3 - 1, g[2] + n / 9 - 1]).filter(move |&k| k != g && self.seen(k))
    }
    /// Whether box `g` (of the walked cell, `tape` its field) holds no
    /// material by one interval over it, or none but on its face (a
    /// billionth inside it, none): no plane through it meets any.
    fn air(&self, lat: &Lattice, tape: &Tape, ctx: &Ctx, g: [i64; 3]) -> bool {
        let i = self.index(g).expect("a walked cell's box");
        if let Some(&a) = self.air.borrow().get(&i) {
            return a;
        }
        let bx = |e: f64| -> [interval::I; 3] { std::array::from_fn(|a| interval::I::new(lat.o[a] + (g[a] as f64 + e) * lat.s, lat.o[a] + (g[a] as f64 + 1.0 - e) * lat.s)) };
        let r = eval_interval(tape, ctx, bx(0.0), None).lo;
        let a = r > 0.0 || r == 0.0 && eval_interval(tape, ctx, bx(1e-9), None).lo > 0.0;
        self.air.borrow_mut().insert(i, a);
        a
    }
    /// Seen sample `k`'s tangent plane (`plane_at`).
    fn plane(&self, lat: &Lattice, ctx: &Ctx, k: [i64; 3]) -> Option<[f64; 4]> {
        let i = self.index(k).expect("a sample beside a walked cell's box");
        if let Some(&p) = self.planes.borrow().get(&i) {
            return p;
        }
        let cached = lat.planes.shard(k).lock().unwrap().get(&k).copied();
        let p = cached.unwrap_or_else(|| {
            // (read outside the lock)
            let p = plane_at(lat, ctx, k);
            lat.planes.shard(k).lock().unwrap().insert(k, p);
            p
        });
        self.planes.borrow_mut().insert(i, p);
        p
    }
    /// Whether material at points `xs` of box `g` (one that sees nothing)
    /// lies within `OWN` of a sample past the tangent plane of every seen
    /// sample beside it AND within `OWN` of at least one such plane on
    /// either side (so a thin plate's far face pointing away from the box
    /// does not licence material past its near face). A sample with no
    /// plane holds nothing of it.
    fn own(&self, lat: &Lattice, ctx: &Ctx, g: [i64; 3], xs: &[V3]) -> bool {
        let tol = OWN * lat.s;
        let mut planes = Vec::new();
        for k in self.beside(g) {
            let Some(p) = self.plane(lat, ctx, k) else { return false };
            planes.push(p);
        }
        if planes.is_empty() {
            return false;
        }
        xs.iter().all(|x| {
            let mut near_one = false;
            for p in &planes {
                let d = v3::dot([p[0], p[1], p[2]], *x) - p[3];
                if d > tol {
                    return false;
                }
                if d >= -tol {
                    near_one = true;
                }
            }
            near_one
        })
    }
}

/// Seen sample `k`'s tangent plane, as the cells place its material
/// (`tangent_offsets`: through the field's zero along its gradient),
/// turned to the field's gradient at that zero where it lies within a
/// sample: a field that is no distance turns its level sets, and the
/// sample's own gradient tilted the plane off the face it places (a
/// lofted blade's skin read a sample past its neighbours' planes). A unit
/// normal out of the material and the offset: n.x - d positive past it.
/// None where the field is flat there (a field held at 0 inside, as a
/// root of a distance's).
fn plane_at(lat: &Lattice, ctx: &Ctx, k: [i64; 3]) -> Option<[f64; 4]> {
    let p: V3 = std::array::from_fn(|a| lat.o[a] + (k[a] as f64 + 0.5) * lat.s);
    let g = tangent_offsets(lat.shape, ctx, &[p], 4.0 * lat.s).1[0];
    let slope = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    if !(slope > 0.0 && slope.is_finite() && g[0].is_finite()) {
        return None;
    }
    let t = g[0] / slope;
    let z: V3 = std::array::from_fn(|a| p[a] - t * g[a + 1] / slope);
    let mut n = [g[1] / slope, g[2] / slope, g[3] / slope];
    if t.abs() < lat.s {
        let gz = eval::eval_grads(lat.shape, ctx, &[z])[0];
        let k = (gz[1] * gz[1] + gz[2] * gz[2] + gz[3] * gz[3]).sqrt();
        if k > 0.0 && k.is_finite() {
            n = [gz[1] / k, gz[2] / k, gz[3] / k];
        }
    }
    Some([n[0], n[1], n[2], v3::dot(n, z)])
}

/// Whether every point on the straight segment from `a` to `b` (each read
/// in material) is shown to lie in material, but for air thinner along it
/// than `eps`: by interval bounds over sub-segments' boxes, halved, each
/// midpoint read in material, to `eps` long, where its ends, both read in
/// material, are taken as joined. (Where pieces meet flush, a leaf set on
/// a block's face, the field is 0 across the seam, and the bounds of a box
/// across it reach past 0 by as far as it reaches past the seam, however
/// small; the kernel cuts such a box at a constant's level,
/// `sign::holds_by_cuts`, but a seam placed by a param has none. A length,
/// not a field value: a field need not be a distance, and one scaled by
/// 1e-12 read a sample of air within a billionth of one of 0.) Fails
/// closed (`false`) where a point reads air, or the budget runs out.
fn segment_in_material(shape: &Tape, ctx: &Ctx, a: V3, b: V3, eps: f64) -> bool {
    let mut stack = vec![(a, b)];
    let mut evals = 0usize;
    while let Some((p, q)) = stack.pop() {
        evals += 1;
        if evals > 256 {
            return false;
        }
        let bx: [interval::I; 3] = std::array::from_fn(|k| interval::I::new(p[k].min(q[k]), p[k].max(q[k])));
        let r = eval_interval(shape, ctx, bx, None);
        if r.hi <= 0.0 {
            continue;
        }
        if r.lo > 0.0 {
            return false;
        }
        if v3::norm(v3::sub(q, p)) <= eps {
            continue;
        }
        let mid: V3 = std::array::from_fn(|k| 0.5 * (p[k] + q[k]));
        let mut v = Vec::new();
        eval_points(shape, ctx, &[mid], &mut v);
        if v[0] > 0.0 {
            return false;
        }
        stack.push((p, mid));
        stack.push((mid, q));
    }
    true
}

/// Whether material at `x` (in box `g`, beside seen samples) is proven to
/// be joined to the material they see along a straight line: interval
/// bounds over pieces of the segment from `x` to the centre of a seen
/// sample beside it inside the part (the nearest three) hold no air.
fn joined(lat: &Lattice, ctx: &Ctx, near: &Near, g: [i64; 3], x: V3) -> bool {
    let seen: Vec<V3> = near.beside(g).map(|k| std::array::from_fn(|a| lat.o[a] + (k[a] as f64 + 0.5) * lat.s)).collect();
    let mut f = Vec::new();
    eval_points(lat.shape, ctx, &seen, &mut f);
    let mut inside: Vec<(f64, V3)> = seen.iter().zip(&f).filter(|(_, f)| **f <= 0.0).map(|(&c, _)| (v3::norm(v3::sub(c, x)), c)).collect();
    inside.sort_by(|a, b| a.0.total_cmp(&b.0));
    inside.iter().take(3).any(|&(_, t)| segment_in_material(lat.shape, ctx, x, t, 1e-9 * lat.s))
}

/// The lattice box rectangle `b` lies in, if one (its midline splits put
/// each quarter in one).
fn box_of(lat: &Lattice, b: &[interval::I; 3]) -> Option<[i64; 3]> {
    let eps = 1e-9 * lat.s;
    let g: [i64; 3] = std::array::from_fn(|a| ((b[a].mid() - lat.o[a]) / lat.s).floor() as i64);
    (0..3).all(|a| b[a].lo >= lat.o[a] + g[a] as f64 * lat.s - eps && b[a].hi <= lat.o[a] + (g[a] + 1) as f64 * lat.s + eps).then_some(g)
}

/// The walk over a rectangle `b` of a plane through samples (a box flat
/// along its plane's axis): where it finds material, or a part of it the
/// bounds cannot settle. A rectangle over several samples' boxes is split
/// into them; a seen sample's box is its own, and not walked. Settled at
/// once where the field's bounds over the box or the rectangle hold no
/// material, where the field's value at its centre over the bound on its
/// slope there (its interval gradient: a field need not be a distance)
/// reaches its corners, or, beside seen samples, where all of it is their
/// face's; else quartered, to the floor (`floor` beside seen samples, a
/// quarter of it elsewhere), and read there at 4 points a side: material a
/// quarter of the floor thick lies over one.
fn walk(tape: &Tape, ctx: &Ctx, b: [interval::I; 3], floor: f64, lat: &Lattice, near: &Near) -> Option<V3> {
    let w: V3 = std::array::from_fn(|a| b[a].hi - b[a].lo);
    let c: V3 = std::array::from_fn(|a| b[a].mid());
    let split = |fl: f64| -> Option<V3> {
        let halves = |a: usize| -> Vec<interval::I> { if w[a] <= fl * (1.0 + 1e-9) { vec![b[a]] } else { vec![interval::I::new(b[a].lo, c[a]), interval::I::new(c[a], b[a].hi)] } };
        let (h0, h1, h2) = (halves(0), halves(1), halves(2));
        h0.iter().find_map(|x| h1.iter().find_map(|y| h2.iter().find_map(|z| walk(tape, ctx, [*x, *y, *z], floor, lat, near))))
    };
    let Some(g) = box_of(lat, &b) else {
        // (one interval first where none of its boxes sees material)
        let k = |a: usize, x: f64| ((x - lat.o[a]) / lat.s).floor() as i64;
        let eps = 1e-9 * lat.s;
        let seen = (0..8).any(|n| near.seen(std::array::from_fn(|a| k(a, if n >> a & 1 == 1 { b[a].hi - eps } else { b[a].lo + eps }))));
        if !seen && eval_interval(tape, ctx, b, None).lo > 0.0 {
            return None;
        }
        return split(floor);
    };
    if near.seen(g) || near.air(lat, tape, ctx, g) {
        return None;
    }
    let r = eval_interval(tape, ctx, b, None);
    if r.lo > 0.0 {
        return None;
    }
    // (touched only, on its edge: a face on the samples' plane, the box
    // beside it; a billionth inside, all air)
    if r.lo == 0.0 {
        let within: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(b[a].lo + 1e-9 * w[a], b[a].hi - 1e-9 * w[a]));
        if eval_interval(tape, ctx, within, None).lo > 0.0 {
            return None;
        }
    }
    let beside = near.beside(g).next().is_some();
    if beside {
        let corners: Vec<V3> = (0..8).map(|n| std::array::from_fn(|a| if n >> a & 1 == 1 { b[a].hi } else { b[a].lo })).collect();
        if near.own(lat, ctx, g, &corners) {
            return None;
        }
    }
    // (material read at x: lost, unless beside seen samples it is their
    // face's)
    let lost = |x: V3| !beside || !near.own(lat, ctx, g, &[x]) && !joined(lat, ctx, near, g, x);
    let fl = if beside { floor } else { floor / 4.0 };
    if w.iter().all(|&w| w <= fl * (1.0 + 1e-9)) {
        // What no point finds and the bounds cannot settle is taken as
        // found, but beside seen samples, where it is their face's.
        let n: [usize; 3] = std::array::from_fn(|a| if w[a] > 0.0 { 4 } else { 1 });
        let pts: Vec<V3> = (0..n[0] * n[1] * n[2])
            .map(|i| {
                let j = [i % n[0], i / n[0] % n[1], i / (n[0] * n[1])];
                std::array::from_fn(|a| b[a].lo + (j[a] as f64 + 0.5) / n[a] as f64 * w[a])
            })
            .collect();
        let mut v = Vec::new();
        eval_points(tape, ctx, &pts, &mut v);
        return pts.iter().zip(&v).find(|&(&x, &v)| v <= 0.0 && lost(x)).map(|(x, _)| *x).or((!beside).then_some(c));
    }
    let mut v = Vec::new();
    eval_points(tape, ctx, &[c], &mut v);
    if v[0] <= 0.0 && lost(c) {
        return Some(c);
    }
    if v[0] > 0.0 && v[0] >= interval::lipschitz(&eval::eval_interval_grad(tape, ctx, b)) * 0.5 * interval::norm(w) {
        return None;
    }
    split(fl)
}

/// A cell (index `i`, centre `ctr`, side `h`; `tape` the field over it)
/// with a sample that sees no material: where its planes through the
/// samples, across each axis and face to face, meet material the samples
/// do not see. Material reaching a spacing along any
/// axis crosses one of them, whichever way it lies (a rod between the
/// lines through the samples, along an axis or not, meets none of them):
/// what crosses none lies within the box between eight samples, so joins
/// nothing. A plane's rectangle between samples that all see material is
/// theirs, and not walked: a part of one material with no feature thinner
/// than a quarter of a cell walks only the boxes beside its face, each
/// settled by its interval or the seen planes beside it.
pub(crate) fn missed_in_cell(tape: &Tape, ctx: &Ctx, i: [i64; 3], ctr: V3, h: f64, lat: &Lattice) -> Option<V3> {
    // the cell's 4 samples a side at -0.75 .. 0.75, and its faces at +-1
    let xs = [-1.0, -0.75, -0.25, 0.25, 0.75, 1.0];
    let at = |a: usize, j: usize| ctr[a] + 0.5 * h * xs[j];
    let floor = h / 8.0;
    let nb = Near::of(lat, i);
    // (no material in the cell's box at all, by one interval: then no plane
    // of it meets any)
    let cb: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * h, ctr[a] + 0.5 * h));
    if eval_interval(tape, ctx, cb, None).lo > 0.0 {
        return None;
    }
    // The samples (`j` 1..=4 along an axis) a point between xs[k] and
    // xs[k + 1] is nearest: a face's rectangle its outer sample's.
    const NEAREST: [&[usize]; 5] = [&[1], &[1, 2], &[2, 3], &[3, 4], &[4]];
    for ax in 0..3 {
        let (u, w) = ((ax + 1) % 3, (ax + 2) % 3);
        for n in 1..5 {
            let rect = |k: [usize; 2], l: [usize; 2]| -> [interval::I; 3] {
                let mut b = [interval::I::new(at(ax, n), at(ax, n)); 3];
                b[u] = interval::I::new(at(u, k[0]), at(u, k[1]));
                b[w] = interval::I::new(at(w, l[0]), at(w, l[1]));
                b
            };
            // (one interval over the whole plane first)
            if eval_interval(tape, ctx, rect([0, 5], [0, 5]), None).lo > 0.0 {
                continue;
            }
            for k in 0..5 {
                for l in 0..5 {
                    let seen = NEAREST[k].iter().all(|&ju| {
                        NEAREST[l].iter().all(|&jw| {
                            let mut j = [0usize; 3];
                            (j[ax], j[u], j[w]) = (n, ju, jw);
                            nb.seen(std::array::from_fn(|a| 4 * i[a] + j[a] as i64 - 1))
                        })
                    });
                    if seen {
                        continue;
                    }
                    if let Some(x) = walk(tape, ctx, rect([k, k + 1], [l, l + 1]), floor, lat, &nb) {
                        return Some(x);
                    }
                }
            }
        }
    }
    None
}

/// The most cells a tip of seen material, touching it at one place, spans
/// along an axis among the seen cells it touches: a cell with its
/// neighbours, and one more.
const TIP: i64 = 4;

/// Of the cells material was found in between the samples (each with where),
/// those joined (26 neighbours) into clusters, each sorted by the cells
/// the samples see material in (`seen`) among and round it: the first
/// cluster that touches none (material the cells lost), seen material in
/// two places apart (a web joining two pieces, or a ligament beside one),
/// or one place longer than a tip's (`TIP`: a rod or a web beside it,
/// joined to it however far apart, or all along), its cell, where, and
/// what it joins. A cluster touching seen material in one place within a
/// tip's cells is its tip, and changes nothing. (Sorted by the seen cells
/// touched alone, rods within a cell of a face touched its cells all along
/// and were one place, however far apart they joined it.)
pub(crate) fn changes_connection(found: Vec<([i64; 3], V3)>, seen: &dyn Fn([i64; 3]) -> bool) -> Option<([i64; 3], V3, String)> {
    joins(found, seen, &|_| false)
}

/// As `changes_connection`, of material in the cells the mesh leaves out
/// (`found`) beside the cells it holds (`seen`), but a cluster touching one
/// connected group of held cells whose every cell is `skin` (each of its
/// material within a sample of held material) is that material's own face,
/// however far it runs along it: a bar's face a hair past a cell's plane,
/// held by a fix, leaves a sheet of such cells all along it. A cluster
/// joining two separate groups (or none) is still refused even if every
/// sliver touches some held cell.
pub(crate) fn joins(found: Vec<([i64; 3], V3)>, seen: &dyn Fn([i64; 3]) -> bool, skin: &dyn Fn([i64; 3]) -> bool) -> Option<([i64; 3], V3, String)> {
    let at: HashMap<[i64; 3], V3> = found.into_iter().collect();
    let nb = |i: [i64; 3]| -> Vec<[i64; 3]> {
        let mut out = Vec::with_capacity(26);
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if (dx, dy, dz) != (0, 0, 0) {
                        out.push([i[0] + dx, i[1] + dy, i[2] + dz]);
                    }
                }
            }
        }
        out
    };
    let mut done: HashSet<[i64; 3]> = HashSet::new();
    let mut keys: Vec<[i64; 3]> = at.keys().copied().collect();
    keys.sort();
    for f in keys {
        if !done.insert(f) {
            continue;
        }
        let (mut cluster, mut stack) = (vec![f], vec![f]);
        while let Some(c) = stack.pop() {
            for n in nb(c) {
                if at.contains_key(&n) && done.insert(n) {
                    cluster.push(n);
                    stack.push(n);
                }
            }
        }
        // (a cell of the cluster may see material itself: what was found
        // lies beside it there)
        let touch: HashSet<[i64; 3]> = cluster.iter().flat_map(|&c| nb(c).into_iter().filter(|n| !at.contains_key(n)).chain([c])).filter(|&n| seen(n)).collect();
        let mut grouped: HashSet<[i64; 3]> = HashSet::new();
        let mut groups = 0;
        for &t in &touch {
            if !grouped.insert(t) {
                continue;
            }
            groups += 1;
            let mut st = vec![t];
            while let Some(c) = st.pop() {
                for n in nb(c) {
                    if touch.contains(&n) && grouped.insert(n) {
                        st.push(n);
                    }
                }
            }
        }
        let span = (0..3).map(|a| touch.iter().map(|t| t[a]).max().unwrap_or(0) - touch.iter().map(|t| t[a]).min().unwrap_or(0) + 1).max().unwrap_or(0);
        if groups == 1 && (span <= TIP || cluster.iter().all(|&c| skin(c))) {
            continue;
        }
        let what = match groups {
            0 => "joins no material the cells hold".to_string(),
            1 => format!("runs beside the material the cells hold over {span} cells, more than a tip's"),
            _ => format!("joins the material the cells hold in {groups} places apart"),
        };
        return Some((f, at[&f], what));
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_line_across_a_flush_seam_is_in_material_and_one_across_a_slit_is_not() {
        // Two slabs, x in [L - 1, L] and [L + gap, L + 1], |y|, |z| < 1, at
        // L a param (no constant level for the kernel to cut at) and at L 0
        // (a constant one), as a leaf set on a block's face. A line from one
        // to the other a sample long crosses the seam where the field is 0:
        // material. With the second slab a hundredth of a sample off, the
        // line crosses air, and is not.
        use ink_core::graph::{Graph, UnOp, Var};
        let s = 0.0625;
        for (param, at) in [(true, 40.0), (false, 0.0)] {
            for gap in [0.0, 0.01 * s] {
                let mut g = Graph::new();
                let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
                let l = if param { g.var(Var::Param(0)) } else { g.constant(at) };
                let slab = |g: &mut Graph, lo: f64, hi: f64| {
                    let (lo, hi) = (g.constant(lo), g.constant(hi));
                    let (lo, hi) = (g.add(l, lo), g.add(l, hi));
                    let (a, b) = (g.sub(x, hi), g.sub(lo, x));
                    let one = g.constant(1.0);
                    let (ay, az) = (g.un(UnOp::Abs, y), g.un(UnOp::Abs, z));
                    let (cy, cz) = (g.sub(ay, one), g.sub(az, one));
                    let m = g.max(a, b);
                    let m = g.max(m, cy);
                    g.max(m, cz)
                };
                let (left, right) = (slab(&mut g, -1.0, 0.0), slab(&mut g, gap, 1.0));
                let part = g.min(left, right);
                let tape = Tape::new(&g, part);
                let ctx = Ctx::new(0.0, if param { vec![at] } else { Vec::new() });
                let (a, b) = ([at - 0.43 * s, 0.17 * s, 0.0], [at + 0.61 * s, -0.29 * s, 0.13 * s]);
                assert_eq!(segment_in_material(&tape, &ctx, a, b, 1e-9 * s), gap == 0.0, "param {param}, gap {gap}");
                // (the field scaled by 1e-12 reads the slit within a
                // billionth of a sample of 0: still air)
                let k = g.constant(1e-12);
                let scaled = g.mul(k, part);
                let tape = Tape::new(&g, scaled);
                assert_eq!(segment_in_material(&tape, &ctx, a, b, 1e-9 * s), gap == 0.0, "scaled: param {param}, gap {gap}");
            }
        }
    }

    #[test]
    fn changes_connection_span_and_empty_touch() {
        // 5-cell run beside seen cells -> refused (span > TIP = 4)
        let run5: Vec<([i64; 3], V3)> = (0..5).map(|x| ([x, 1, 0], [x as f64, 1.0, 0.0])).collect();
        let res = changes_connection(run5, &|c| c[1] == 0 && (-1..=5).contains(&c[0]) && c[2] == 0);
        assert!(res.as_ref().is_some_and(|(_, _, w)| w.contains("more than a tip's")), "{res:?}");

        // 2-cell run beside seen cells -> accepted as a tip (1 group, span <= 4)
        let run2: Vec<([i64; 3], V3)> = (0..2).map(|x| ([x, 1, 0], [x as f64, 1.0, 0.0])).collect();
        let res = changes_connection(run2, &|c| c[1] == 0 && (0..=1).contains(&c[0]) && c[2] == 0);
        assert!(res.is_none(), "{res:?}");

        // Cluster touching no seen cells -> refused ("joins no material")
        let isolated = vec![([0, 0, 0], [0.0, 0.0, 0.0])];
        let res = changes_connection(isolated, &|_| false);
        assert!(res.as_ref().is_some_and(|(_, _, w)| w.contains("joins no material")), "{res:?}");

        // A skin cluster that bridges two separate seen groups is still refused.
        let bridge = vec![([0, 1, 0], [0.0, 1.0, 0.0]), ([0, 2, 0], [0.0, 2.0, 0.0])];
        let res = joins(bridge, &|c| c == [0, 0, 0] || c == [0, 3, 0], &|_| true);
        assert!(res.as_ref().is_some_and(|(_, _, w)| w.contains("2 places apart")), "{res:?}");
    }

    #[test]
    fn a_skin_along_held_material_joins_nothing_and_a_foil_beside_it_does() {
        // A sheet of 20 cells along seen cells: refused as a run, but a skin
        // (each cell's material against them) changes nothing; a sheet with
        // one cell off them is sorted as the run it is.
        let sheet: Vec<([i64; 3], V3)> = (0..20).map(|x| ([x, 1, 0], [x as f64, 1.0, 0.0])).collect();
        let seen = |c: [i64; 3]| c[1] == 0 && (0..20).contains(&c[0]) && c[2] == 0;
        assert!(joins(sheet.clone(), &seen, &|_| false).is_some());
        assert!(joins(sheet.clone(), &seen, &|_| true).is_none());
        let res = joins(sheet, &seen, &|c| c[0] != 10);
        assert!(res.as_ref().is_some_and(|(_, _, w)| w.contains("more than a tip's")), "{res:?}");
    }
}
