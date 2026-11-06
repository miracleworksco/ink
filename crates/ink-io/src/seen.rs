//! Meshes for seeing and for files (a GLB, a render, STL and 3MF): a field
//! meshed and simplified, with the field's own normals at each triangle's corners for
//! shading, and its paint (colour, metallic, roughness, transmission) cut where it
//! changes: no triangle spans two paints, so a label or a seam between
//! colours is a sharp line, not the zigzag of the triangles it crosses. Cached on disk by the question (ink_kernel::cache): the tape,
//! the box, the context it reads and every option. A part's mesh depends on the part alone, so an unchanged part
//! is not meshed again, in any assembly.

use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::tape::Tape;
use ink_kernel::mesh::{Mesh, Unresolved, mesh_checked};
use ink_kernel::octree::Local;
use std::sync::Arc;

/// How a viewing mesh is made: sampled at `res`, simplified to `tol` (else
/// `tol_per_res` of res), keeping triangles short along x for a strip bent
/// round arcs as tight as `bend`.
#[derive(Clone, Copy, Debug)]
pub struct How {
    pub res: f64,
    pub tol: Option<f64>,
    pub tol_per_res: f64,
    pub bend: Option<f64>,
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Viewed {
    pub verts: Vec<V3>,
    pub tris: Vec<[u32; 3]>,
    /// The field's unit normal at each triangle's corners, three per
    /// triangle (`corner_normals`).
    pub normals: Vec<[f32; 3]>,
    /// Where each vertex's paint is read: `BENEATH_PER_RES` of the res
    /// under the surface by it, in the material there (`beneath`).
    pub paint_at: Vec<V3>,
    /// The res it was sampled at and the tolerance it was simplified to.
    pub res: f64,
    pub tol: f64,
    /// Triangles at res, before simplification.
    pub uniform_triangles: usize,
    /// The instructions reading its paint ran (`Paint::reads`): beside its
    /// triangles, what making it cost.
    pub paint_work: usize,
    /// The farthest a vertex at res moved from the result (<= tol).
    pub moved: f64,
    /// Sign changes between samples, which the mesh leaves out.
    pub unresolved: Unresolved,
    /// Found in the cache, not made.
    #[serde(skip)]
    pub cached: bool,
}

impl Viewed {
    pub fn mesh(&self) -> Mesh {
        Mesh { verts: self.verts.clone(), tris: self.tris.clone() }
    }
}

/// Meshes that were found in the cache and made, in this process.
pub static HITS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
pub static MADE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// A change in a paint value (a colour channel, metallic or roughness, 0
/// to 1) that shows.
const PAINT_STEP: f64 = 1.0 / 64.0;

/// Paint is read this fraction of the res beneath the surface (`beneath`):
/// a surface has the paint of the material under it, so a paint chosen by
/// a body's own field (a belt black where the belt is, `select belt`) does
/// not flip from point to point on the surface, where that field is zero.
/// Small, so material that thin, sampled at res, still has its own paint.
const BENEATH_PER_RES: f64 = 1e-3;

/// The mesh of `tape` over `lo..hi`, from the cache if it was made before.
/// `paint`: the body's paint values that vary over space; the triangles
/// across which one steps from one paint to another are kept as sampled
/// (no simplified triangle spans two paints) and cut along the step
/// (`cut_paint`).
pub fn mesh(tape: &Arc<Tape>, paint: &[Arc<Tape>], ctx: &Ctx, lo: V3, hi: V3, how: How) -> Viewed {
    let (key, ctx) = key(tape, paint, ctx, lo, hi, how);
    if let Some(v) = key.and_then(cached) {
        return v;
    }
    MADE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let v = make(tape, paint, &ctx, lo, hi, how);
    ink_kernel::cache::computed();
    if let Some(k) = key {
        ink_kernel::cache::store(KIND, k, &v);
    }
    v
}

/// About as many points as `surface` takes of a surface.
const AREA_SAMPLES: f64 = 1000.0;

/// The finest `surface` samples one: cells across its box's longest side.
const AREA_CELLS: f64 = 2048.0;

/// A shape's surface, as a view mesh of it is judged by (`surface`).
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct Surface {
    /// Its area, all of it: what a mesh of it samples, about two triangles
    /// per res squared.
    pub area: f64,
    /// Points where it bounds material (the field changes sign across
    /// them), each with its unit outward normal and the area it stands
    /// for. A face of no thickness (where a subtract leaves one) bounds
    /// none: no res draws it.
    pub bounds: Vec<(V3, V3, f64)>,
}

/// How far a mesh lies from a surface's boundary points (`Surface::off`).
#[derive(Clone, Copy, Debug)]
pub struct Off {
    /// Their area-weighted mean distance from it: infinite where it has no
    /// triangles, 0 where the surface bounds no material.
    pub mean: f64,
    /// The thickest the material is at those farther than its res from it
    /// (0: none): a res under that draws them.
    pub thickest: f64,
}

impl Surface {
    /// How far `mesh`, sampled at res `h`, lies from the surface's
    /// boundary points, and how thick the material it leaves out is (the
    /// field `tape`'s).
    pub fn off(&self, mesh: &Mesh, h: f64, tape: &Tape, ctx: &Ctx) -> Off {
        use rayon::prelude::*;
        let area: f64 = self.bounds.iter().map(|b| b.2).sum();
        if area <= 0.0 {
            return Off { mean: 0.0, thickest: 0.0 };
        }
        let tris: Vec<[V3; 3]> = mesh.tris.iter().map(|t| t.map(|i| mesh.verts[i as usize])).collect();
        let near = crate::stl::MeshDistance::new(&tris);
        let d: Vec<f64> = match near.triangles() {
            0 => vec![f64::INFINITY; self.bounds.len()],
            _ => self.bounds.par_iter().map(|b| near.distance(b.0).abs()).collect(),
        };
        let mean = d.iter().zip(&self.bounds).map(|(d, b)| d * b.2).sum::<f64>() / area;
        let far: Vec<(V3, V3)> = d.iter().zip(&self.bounds).filter(|(d, _)| **d > h).map(|(_, b)| (b.0, b.1)).collect();
        Off { mean, thickest: thickest(tape, ctx, &far, h) }
    }
}

/// Samples a cell across, either way, on a line through the material.
const ACROSS: i32 = 16;

/// The thickest the material is behind boundary points (point, outward
/// normal): at each, its shortest crossing along the normal and the three
/// axes from a point just inside, from where the field crosses zero either
/// side (between samples a sixteenth of `h` apart, up to `h` away), as the
/// mesher measures material it leaves out (ink_kernel::mesh::Unresolved).
fn thickest(tape: &Tape, ctx: &Ctx, pts: &[(V3, V3)], h: f64) -> f64 {
    let step = h / ACROSS as f64;
    let n = 2 * ACROSS as usize + 1;
    let mut line = Vec::with_capacity(4 * n * pts.len());
    for (p, u) in pts {
        let q: V3 = std::array::from_fn(|a| p[a] - 1e-3 * h * u[a]);
        for u in [*u, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]] {
            line.extend((-ACROSS..=ACROSS).map(|i| std::array::from_fn(|a| q[a] + step * i as f64 * u[a])));
        }
    }
    let mut v = Vec::new();
    ink_core::eval::eval_points(tape, ctx, &line, &mut v);
    let inside = |f: f64| f < 0.0;
    let across = |f: &[f64]| {
        // The run of samples inside through the middle, and how far past
        // each end of it the field crosses zero.
        let (mut lo, mut hi) = (ACROSS as usize, ACROSS as usize);
        if !inside(f[lo]) {
            return 0.0;
        }
        while lo > 0 && inside(f[lo - 1]) {
            lo -= 1;
        }
        while hi + 1 < f.len() && inside(f[hi + 1]) {
            hi += 1;
        }
        let past = |a: usize, b: Option<usize>| b.map_or(0.0, |b| f[a] / (f[a] - f[b]));
        step * ((hi - lo) as f64 + past(lo, lo.checked_sub(1)) + past(hi, (hi + 1 < f.len()).then_some(hi + 1)))
    };
    v.chunks(4 * n).map(|f| f.chunks(n).map(across).fold(f64::INFINITY, f64::min)).fold(0.0, f64::max)
}

/// A shape's `Surface` in `lo..hi` (which holds its material): one point
/// per cell it crosses (ink_kernel::feature::surface), on it, with the area
/// of its tangent plane in the cell, from 16 cells across, finer by the
/// square root of the points wanted over those found (they go as the
/// square of the cells across where the surface is resolved, more slowly
/// where it is thinner than a cell: a long extrusion's profile, a belt),
/// until about AREA_SAMPLES, up to AREA_CELLS across. Its area is a
/// sample: a belt's runs read 44 to 82 k mm2 frame to frame against 50 k
/// drawn, and a surface folded within a cell (a slotted extrusion's
/// profile, a sheet's two sides) counts once there, a half to a third of
/// it. Empty where none is found. Cached.
pub fn surface(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Surface {
    ink_kernel::cache::memo_region(
        "seen surface",
        &[&**tape],
        ctx,
        (lo, hi),
        |_| {},
        |ctx| {
            let mut cells: f64 = 16.0;
            loop {
                let (found, h) = ink_kernel::feature::surface(tape, ctx, lo, hi, cells);
                let n = found.len() as f64;
                if n >= AREA_SAMPLES / 2.0 || cells >= AREA_CELLS {
                    // At most twice AREA_SAMPLES kept, each standing for those
                    // skipped; of them, those with material on one side only.
                    let every = found.len().div_ceil(2 * AREA_SAMPLES as usize).max(1);
                    let kept: Vec<(V3, V3, f64)> = found.chunks(every).map(|c| (c[0].0, c[0].1, c.iter().map(|p| p.2).sum::<f64>())).collect();
                    let e = 1e-3 * h;
                    let probes: Vec<V3> = kept.iter().flat_map(|(p, u, _)| [-e, e].map(|s| std::array::from_fn(|a| p[a] + s * u[a]))).collect();
                    let mut v = Vec::new();
                    ink_core::eval::eval_points(tape, ctx, &probes, &mut v);
                    let bounds = kept.into_iter().zip(v.chunks(2)).filter(|(_, f)| (f[0] < 0.0) != (f[1] < 0.0)).map(|(b, _)| b).collect();
                    return Surface { area: found.iter().map(|p| p.2).sum::<f64>(), bounds };
                }
                cells = (cells * (AREA_SAMPLES / n.max(1.0)).sqrt()).min(AREA_CELLS);
            }
        },
    )
}

fn key<'a>(tape: &Arc<Tape>, paint: &[Arc<Tape>], ctx: &'a Ctx, lo: V3, hi: V3, how: How) -> (Option<u128>, std::borrow::Cow<'a, Ctx>) {
    let tapes: Vec<&Tape> = std::iter::once(&**tape).chain(paint.iter().map(|t| &**t)).collect();
    ink_kernel::cache::ask_region("seen mesh", &tapes, ctx, (lo, hi), |fp| {
        fp.f64(how.res).f64(how.tol.unwrap_or(-1.0)).f64(how.tol_per_res);
        fp.f64(how.bend.unwrap_or(-1.0));
    })
}

fn cached(key: u128) -> Option<Viewed> {
    let v = ink_kernel::cache::load::<Viewed>(KIND, key)?;
    HITS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    ink_kernel::cache::served();
    Some(Viewed { cached: true, ..v })
}

const KIND: &str = "meshes";

fn make(tape: &Arc<Tape>, paint: &[Arc<Tape>], ctx: &Ctx, lo: V3, hi: V3, how: How) -> Viewed {
    made(tape, paint, ctx, lo, hi, how).0
}

/// `make`, and the points its paint was read at and the instructions that
/// cost (`Paint::reads`; the second is the mesh's `paint_work`).
fn made(tape: &Arc<Tape>, paint: &[Arc<Tape>], ctx: &Ctx, lo: V3, hi: V3, how: How) -> (Viewed, [usize; 2]) {
    let h = how.res;
    let (fine, unresolved) = mesh_checked(tape, ctx, lo, hi, h);
    let tol = how.tol.unwrap_or(how.tol_per_res * h);
    let uniform_triangles = fine.tris.len();
    // The field near the surface, a tape per block of 8 cells.
    let local = Local::new(tape, ctx, lo, 8.0 * h);
    let field = |p: &[V3]| local.eval(p);
    let under = |p: &[V3]| beneath(&local, p, BENEATH_PER_RES * h);
    let paint = Paint { tapes: paint, ctx, under: &under, origin: lo, size: 8.0 * h, reads: Default::default() };
    let fine = split_paint(fine, &paint, &local, h);
    let (m, moved) = if tol > 0.0 {
        let keep = paint_steps(&fine, &paint);
        let simple = crate::simplify::simplify(&fine, tol, Some(&field), how.bend.map(|r| (8.0 * r * tol).sqrt()), &keep);
        (simple.mesh, simple.moved.iter().fold(0.0f64, |a, &d| a.max(d)))
    } else {
        (fine, 0.0)
    };
    let m = cut_paint(m, &paint, h);
    let normals = corner_normals(&m, &|p: &[V3]| local.grads(p));
    let paint_at = under(&m.verts);
    let reads = paint.reads.get();
    let v = Viewed { verts: m.verts, tris: m.tris, normals, paint_at, res: h, tol, uniform_triangles, paint_work: reads[1], moved, unresolved, cached: false };
    (v, reads)
}

/// The vertices of the triangles across which a paint value steps: along
/// one of their edges it changes by more than `PAINT_STEP`, and the change
/// is a jump, not a gradient or a pattern finer than the sampling (a
/// motor's laminations, which simplifying blends as the samples do): it
/// lies between one pair of the edge's quarter points, then of that
/// quarter's, and is still over half of it there.
fn paint_steps(m: &Mesh, paint: &Paint) -> Vec<bool> {
    let mut keep = vec![false; if paint.tapes.is_empty() { 0 } else { m.verts.len() }];
    if paint.tapes.is_empty() {
        return keep;
    }
    for (f, s) in stepping(m, &paint.at(&m.verts), paint).into_iter().enumerate() {
        if s {
            m.tris[f].iter().for_each(|&i| keep[i as usize] = true);
        }
    }
    keep
}

/// Whether a paint value steps across each triangle (`paint_steps`), given
/// the paints at its vertices (`Paint::at`).
fn stepping(m: &Mesh, at: &[f64], paint: &Paint) -> Vec<bool> {
    let mut steps = vec![false; m.tris.len()];
    let lerp = |a: V3, b: V3, s: f64| -> V3 { std::array::from_fn(|k| a[k] + s * (b[k] - a[k])) };
    let n = paint.tapes.len();
    for j in 0..n {
        let eval = |pts: &[V3]| paint.read(j..j + 1, pts);
        let v: Vec<f64> = at.iter().skip(j).step_by(n).copied().collect();
        // Each edge it changes along: its triangle, and the stretch the
        // change lies in (ends and their values).
        let mut open: Vec<(usize, [V3; 2], [f64; 2], f64)> = Vec::new();
        for (f, t) in m.tris.iter().enumerate() {
            for e in 0..3 {
                let (a, b) = (t[e] as usize, t[(e + 1) % 3] as usize);
                let d = (v[a] - v[b]).abs();
                if d > PAINT_STEP {
                    open.push((f, [m.verts[a], m.verts[b]], [v[a], v[b]], d));
                }
            }
        }
        for _ in 0..2 {
            let pts: Vec<V3> = open.iter().flat_map(|(_, [a, b], _, _)| [0.25, 0.5, 0.75].map(|s| lerp(*a, *b, s))).collect();
            let q = eval(&pts);
            open = open
                .into_iter()
                .enumerate()
                .filter_map(|(k, (f, [a, b], [va, vb], d))| {
                    let vals = [va, q[3 * k], q[3 * k + 1], q[3 * k + 2], vb];
                    let mut jumps = (0..4).filter(|&i| (vals[i + 1] - vals[i]).abs() > PAINT_STEP);
                    let (i, None) = (jumps.next()?, jumps.next()) else { return None };
                    let (s0, s1) = (i as f64 / 4.0, (i + 1) as f64 / 4.0);
                    Some((f, [lerp(a, b, s0), lerp(a, b, s1)], [vals[i], vals[i + 1]], d))
                })
                .collect();
        }
        for (f, _, [va, vb], d) in open {
            if (vb - va).abs() > 0.5 * d {
                steps[f] = true;
            }
        }
    }
    steps
}

/// Paint finer than the sampling, which a mesh drawing each triangle in its
/// corners' paint would leave out: a label's strokes narrower than res lie
/// between vertices, so whole letters (a T's stem) vanish or break, as the
/// sample points happen to fall. Each edge along which a paint value leaves
/// the range of its ends' by more than a step, at one of its quarter points
/// (a stroke or line crossing it, not a gradient, which stays within its
/// ends'), is halved: its middle moved onto the surface, the triangles round
/// it split so the mesh stays closed; twice, so a stroke down to a quarter
/// of res across has vertices in it, and `cut_paint` cuts along its sides.
/// Only edges in blocks of 8 cells over which some paint is not proven
/// constant to a step (by its range there) are looked at.
fn split_paint(mut m: Mesh, paint: &Paint, local: &Local, h: f64) -> Mesh {
    use ink_core::interval::I;
    use std::collections::HashMap;
    let n = paint.tapes.len();
    if n == 0 {
        return m;
    }
    let size = 8.0 * h;
    let block = |p: V3| -> [i64; 3] { std::array::from_fn(|a| (p[a] / size).floor() as i64) };
    let mut varies: HashMap<[i64; 3], bool> = HashMap::new();
    let mid = |a: V3, b: V3, s: f64| -> V3 { std::array::from_fn(|k| a[k] + s * (b[k] - a[k])) };
    // Each paint at each vertex (read once; the middles added are read in
    // the next pass), and the first vertex the last pass added: an edge
    // without one was looked at then and not split, so is not again.
    let mut at_vert: Vec<f64> = Vec::new();
    let mut fresh = 0u32;
    for _ in 0..2 {
        at_vert.extend(paint.at(&m.verts[at_vert.len() / n..]));
        let mut edges: Vec<(u32, u32)> =
            m.tris.iter().flat_map(|t| (0..3).map(move |e| (t[e].min(t[(e + 1) % 3]), t[e].max(t[(e + 1) % 3])))).filter(|&(_, b)| b >= fresh).collect();
        edges.sort_unstable();
        edges.dedup();
        // Those in blocks where a paint may change: a block's box, grown by
        // a cell (an edge's reach past its middle's block), each paint's range.
        let keys: Vec<[i64; 3]> = edges.iter().map(|&(a, b)| block(mid(m.verts[a as usize], m.verts[b as usize], 0.5))).collect();
        let mut new_keys: Vec<[i64; 3]> = keys.iter().filter(|k| !varies.contains_key(*k)).copied().collect();
        new_keys.sort_unstable();
        new_keys.dedup();
        let found: Vec<bool> = {
            use rayon::prelude::*;
            new_keys
                .par_iter()
                .map(|k| {
                    let b: [I; 3] = std::array::from_fn(|a| I { lo: k[a] as f64 * size - h, hi: (k[a] + 1) as f64 * size + h });
                    paint.tapes.iter().any(|p| {
                        let r = ink_core::eval::eval_interval(p, paint.ctx, b, None);
                        !(r.hi.clamp(0.0, 1.0) - r.lo.clamp(0.0, 1.0) <= PAINT_STEP)
                    })
                })
                .collect()
        };
        varies.extend(new_keys.into_iter().zip(found));
        let edges: Vec<(u32, u32)> = edges.into_iter().zip(&keys).filter(|(_, k)| varies[*k]).map(|(e, _)| e).collect();
        // Each paint at the quarter points of those edges.
        let pts: Vec<V3> = edges
            .iter()
            .flat_map(|&(a, b)| {
                let (a, b) = (m.verts[a as usize], m.verts[b as usize]);
                [0.25, 0.5, 0.75].map(|s| mid(a, b, s))
            })
            .collect();
        let vals = paint.at(&pts);
        let crossed = |k: usize| {
            let (a, b) = (edges[k].0 as usize, edges[k].1 as usize);
            (0..n).any(|j| {
                let (va, vb) = (at_vert[a * n + j], at_vert[b * n + j]);
                (0..3).map(|s| vals[(3 * k + s) * n + j]).any(|q| q < va.min(vb) - PAINT_STEP || q > va.max(vb) + PAINT_STEP)
            })
        };
        let split: Vec<(u32, u32)> = (0..edges.len()).filter(|&k| crossed(k)).map(|k| edges[k]).collect();
        if split.is_empty() {
            break;
        }
        // Each split edge's middle, on the surface (a Newton step from the
        // chord's middle, kept within a quarter of the edge of it).
        let chords: Vec<V3> = split.iter().map(|&(a, b)| mid(m.verts[a as usize], m.verts[b as usize], 0.5)).collect();
        let g = local.grads(&chords);
        let base = m.verts.len() as u32;
        fresh = base;
        let mut middle: HashMap<(u32, u32), u32> = HashMap::with_capacity(split.len());
        for (k, (&(a, b), (c, g))) in split.iter().zip(chords.iter().zip(&g)).enumerate() {
            let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
            let step: V3 = if n2 > 0.0 && n2.is_finite() { [1, 2, 3].map(|i| g[0] * g[i] / n2) } else { [0.0; 3] };
            let reach = 0.25 * ink_core::v3::norm(ink_core::v3::sub(m.verts[a as usize], m.verts[b as usize]));
            let p = if ink_core::v3::norm(step) <= reach { ink_core::v3::sub(*c, step) } else { *c };
            m.verts.push(p);
            middle.insert((a, b), base + k as u32);
        }
        // Each triangle split by its halved edges, its corners' order kept.
        let half = |a: u32, b: u32| middle.get(&(a.min(b), a.max(b))).copied();
        let mut tris = Vec::with_capacity(m.tris.len() + 3 * split.len());
        for t in &m.tris {
            let mids = [0, 1, 2].map(|e| half(t[e], t[(e + 1) % 3]));
            let r = |i: usize| t[i % 3];
            match mids.iter().filter(|x| x.is_some()).count() {
                0 => tris.push(*t),
                1 => {
                    let e = mids.iter().position(Option::is_some).unwrap();
                    let (p, q, s, pq) = (r(e), r(e + 1), r(e + 2), mids[e].unwrap());
                    tris.extend([[p, pq, s], [pq, q, s]]);
                }
                2 => {
                    // The whole edge is s-p: p-q and q-s are halved.
                    let u = mids.iter().position(Option::is_none).unwrap();
                    let (s, p, q) = (r(u), r(u + 1), r(u + 2));
                    let (pq, qs) = (mids[(u + 1) % 3].unwrap(), mids[(u + 2) % 3].unwrap());
                    tris.extend([[pq, q, qs], [p, pq, qs], [p, qs, s]]);
                }
                _ => {
                    let [ab, bc, ca] = mids.map(Option::unwrap);
                    tris.extend([[t[0], ab, ca], [ab, t[1], bc], [ca, bc, t[2]], [ab, bc, ca]]);
                }
            }
        }
        m.tris = tris;
    }
    m
}

/// Points of the surface moved to where their paint is read (`beneath`).
type Under<'a> = dyn Fn(&[V3]) -> Vec<V3> + 'a;

/// Points by the surface moved into the material under them: to `d` under
/// the surface nearest them, as the field's value and gradient place it
/// (Newton steps: one on a face, one per face meeting at an edge or
/// corner), kept once at least half that deep. So a point is read beside
/// its own stretch of surface, not across a gap in material beyond (a gap
/// can be stepped over only where it is narrower than the point's distance
/// from the surface, plus `d`), and a sheet down to `d` thick is read
/// inside. Points it cannot take there (no gradient, material thinner than
/// `d`) are left where they are. (Not at the mesher's own samples in the
/// material: those are on its cells' edges, up to a cell sideways from a
/// vertex placed within the cell, which would move a label's edge.)
fn beneath(local: &Local, pts: &[V3], d: f64) -> Vec<V3> {
    let mut at = pts.to_vec();
    let mut todo: Vec<usize> = (0..pts.len()).collect();
    for step in 0..4 {
        let g = local.grads(&todo.iter().map(|&i| at[i]).collect::<Vec<_>>());
        todo = todo
            .into_iter()
            .zip(g)
            .filter_map(|(i, g)| {
                let n = ink_core::v3::norm([g[1], g[2], g[3]]);
                if g[0] <= -0.5 * d * n {
                    return None; // deep enough
                }
                if !(n > 0.0 && n.is_finite()) || step == 3 {
                    at[i] = pts[i];
                    return None;
                }
                let t = (g[0] + d * n) / (n * n);
                at[i] = std::array::from_fn(|k| at[i][k] - t * g[k + 1]);
                Some(i)
            })
            .collect();
    }
    at
}

/// A body's paint values that vary over space, read at points of its
/// surface: each beneath it (`beneath`), on the paint's tape simplified to
/// each block of 8 cells the points are in, as the field's is (`Local`),
/// so a paint chosen among many parts (a wheel's rollers, each its own
/// colour) reads in a block only the parts near it. Points are read in
/// runs, each of the points in `RUN_BLOCKS` blocks together in space
/// (`runs`), the run's simplified tapes dropped after it: the blocks held
/// are a run's for each paint, not the whole surface's.
struct Paint<'a> {
    tapes: &'a [Arc<Tape>],
    ctx: &'a Ctx,
    under: &'a Under<'a>,
    origin: V3,
    size: f64,
    /// Paint values read so far, and the instructions they cost.
    reads: std::cell::Cell<[usize; 2]>,
}

/// The blocks of 8 cells a run of a `Paint`'s points lies in.
const RUN_BLOCKS: usize = 2048;

impl Paint<'_> {
    /// Every paint value at each point (0 to 1), the point's values together.
    fn at(&self, pts: &[V3]) -> Vec<f64> {
        self.read(0..self.tapes.len(), pts)
    }

    /// The paint values `which` at each point, the point's values together.
    fn read(&self, which: std::ops::Range<usize>, pts: &[V3]) -> Vec<f64> {
        let n = which.len();
        let mut out = vec![0.0; n * pts.len()];
        for run in self.runs(pts) {
            let v = self.run().read(which.clone(), &run.iter().map(|&i| pts[i as usize]).collect::<Vec<_>>());
            for (&i, v) in run.iter().zip(v.chunks(n.max(1))) {
                out[i as usize * n..(i as usize + 1) * n].copy_from_slice(v);
            }
        }
        out
    }

    /// The points in runs (their indices), each the points in at most
    /// `RUN_BLOCKS` blocks, taken in Morton order (so a run's blocks are
    /// together in space and share their larger blocks).
    fn runs(&self, pts: &[V3]) -> Vec<Vec<u32>> {
        use rayon::prelude::*;
        // Each block coordinate's bits, every third bit, to 21 bits.
        let spread = |v: u64| (0..21).fold(0u64, |m, b| m | ((v >> b) & 1) << (3 * b));
        let (origin, size) = (self.origin, self.size);
        let key = |p: &V3| -> u64 {
            (0..3).fold(0, |m, a| {
                let k = ((p[a] - origin[a]) / size).floor() + 1.0;
                m | spread(k.clamp(0.0, (1u64 << 21) as f64 - 1.0) as u64) << a
            })
        };
        let mut order: Vec<(u64, u32)> = pts.par_iter().enumerate().map(|(i, p)| (key(p), i as u32)).collect();
        order.par_sort_unstable();
        let mut runs: Vec<Vec<u32>> = vec![Vec::new()];
        let (mut blocks, mut last) = (0, None);
        for (k, i) in order {
            if last != Some(k) {
                if blocks == RUN_BLOCKS {
                    runs.push(Vec::new());
                    blocks = 0;
                }
                (blocks, last) = (blocks + 1, Some(k));
            }
            runs.last_mut().unwrap().push(i);
        }
        runs
    }

    /// A reader for one run of points (which keeps the blocks it reads).
    fn run(&self) -> Run<'_> {
        Run { paint: self, local: self.tapes.iter().map(|t| Local::new(t, self.ctx, self.origin, self.size)).collect() }
    }
}

/// A `Paint` reading one run of points: each paint's blocks, kept while
/// it is.
struct Run<'a> {
    paint: &'a Paint<'a>,
    local: Vec<Local<'a>>,
}

impl Run<'_> {
    fn at(&self, pts: &[V3]) -> Vec<f64> {
        self.read(0..self.local.len(), pts)
    }

    fn read(&self, which: std::ops::Range<usize>, pts: &[V3]) -> Vec<f64> {
        let n = which.len();
        let mut out = vec![0.0; n * pts.len()];
        let at = (self.paint.under)(pts);
        for (j, local) in self.local[which].iter().enumerate() {
            let work = local.work();
            local.eval(&at).into_iter().enumerate().for_each(|(i, v)| out[i * n + j] = v.clamp(0.0, 1.0));
            let [r, w] = self.paint.reads.get();
            self.paint.reads.set([r + pts.len(), w + local.work() - work]);
        }
        out
    }
}

/// Whether two points have the same paint: no value differs by a step.
fn same(a: &[f64], b: &[f64]) -> bool {
    apart(a, b) <= PAINT_STEP
}

/// How far apart two paints are: the most any value differs.
fn apart(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |d, (x, y)| d.max((x - y).abs()))
}

/// The mesh with its paint cut where it steps: each triangle whose corners
/// have two paints, one corner against two, is cut where the paint changes
/// along its two edges from that corner (bisected, each middle going to the
/// end whose paint it is nearer, to a millionth of the edge): the corner's
/// triangle and the other two's quad, each with its own vertices at the cut
/// (a millionth of an edge apart), so each has one paint. Its triangles are
/// as sampled at `res` (`paint_steps` kept them from simplifying), so a
/// curve of the paint's boundary within one is cut as a chord of a cell. A
/// triangle with three paints at its corners is left whole, and one whose
/// paint changes gradually (no step) is not cut.
fn cut_paint(mut m: Mesh, paint: &Paint, res: f64) -> Mesh {
    use std::collections::HashMap;
    let n = paint.tapes.len();
    if n == 0 {
        return m;
    }
    let at = paint.at(&m.verts);
    let steps = stepping(&m, &at, paint);
    let at = |i: u32| &at[i as usize * n..(i as usize + 1) * n];
    // Each cut triangle: its lone corner, then the other two, in its order.
    let lone: Vec<Option<[u32; 3]>> = m
        .tris
        .iter()
        .zip(&steps)
        .map(|(&t, &s)| {
            if !s {
                return None;
            }
            let p = t.map(at);
            let (ab, bc, ca) = (same(p[0], p[1]), same(p[1], p[2]), same(p[2], p[0]));
            match (ab, bc, ca) {
                (false, true, false) => Some(t),
                (false, false, true) => Some([t[1], t[2], t[0]]),
                (true, false, false) => Some([t[2], t[0], t[1]]),
                _ => None,
            }
        })
        .collect();
    // The edges cut, from a lone corner to another, and where: bisected
    // together, `ends[k]` keeping one end with the lone corner's paint.
    let mut cuts: HashMap<(u32, u32), usize> = HashMap::new();
    let mut ends: Vec<(u32, u32, V3, V3)> = Vec::new();
    for c in lone.iter().flatten() {
        for o in [c[1], c[2]] {
            cuts.entry((c[0].min(o), c[0].max(o))).or_insert_with(|| {
                ends.push((c[0], o, m.verts[c[0] as usize], m.verts[o as usize]));
                ends.len() - 1
            });
        }
    }
    // In runs of edges together in space, each reading its blocks once.
    let firsts: Vec<V3> = ends.iter().map(|e| e.2).collect();
    for edges in paint.runs(&firsts) {
        let run = paint.run();
        for _ in 0..20 {
            let mids: Vec<V3> = edges.iter().map(|&k| std::array::from_fn(|a| 0.5 * (ends[k as usize].2[a] + ends[k as usize].3[a]))).collect();
            let pm = run.at(&mids);
            for (&k, (pm, x)) in edges.iter().zip(pm.chunks(n).zip(mids)) {
                let e = &mut ends[k as usize];
                if apart(pm, at(e.0)) <= apart(pm, at(e.1)) {
                    e.2 = x;
                } else {
                    e.3 = x;
                }
            }
        }
    }
    // Each cut edge's two new vertices: the one by its end `a` first.
    let base = m.verts.len() as u32;
    for &(_, _, p, q) in &ends {
        m.verts.extend([p, q]);
    }
    let near = |a: u32, b: u32| -> u32 {
        let k = cuts[&(a.min(b), a.max(b))];
        base + 2 * k as u32 + if ends[k].0 == a { 0 } else { 1 }
    };
    let mut tris = Vec::with_capacity(m.tris.len() + 2 * cuts.len());
    for (t, c) in m.tris.iter().zip(&lone) {
        match c {
            None => tris.push(*t),
            Some([l, a, b]) => {
                let (la, lb) = (near(*l, *a), near(*l, *b));
                let (al, bl) = (near(*a, *l), near(*b, *l));
                tris.extend([[*l, la, lb], [al, *a, *b], [al, *b, bl]]);
            }
        }
    }
    // Less those of no area: a cut through a corner on the paint's edge.
    let area = |t: &[u32; 3]| {
        ink_core::v3::norm(ink_core::v3::cross(
            ink_core::v3::sub(m.verts[t[1] as usize], m.verts[t[0] as usize]),
            ink_core::v3::sub(m.verts[t[2] as usize], m.verts[t[0] as usize]),
        ))
    };
    tris.retain(|t| area(t) > 1e-12 * res * res);
    m.tris = tris;
    m
}

/// Corners within this angle of their vertex's normal are on a smooth
/// surface there (cos 20 degrees): the triangles around a vertex of a mesh
/// simplified to tol differ by far less on a smooth surface, and the
/// sides of an edge by more.
const SMOOTH_COS: f64 = 0.94;

/// The field's unit normal at each triangle's corners: at the vertex where
/// the surface is smooth, so shading is continuous across triangles; else
/// (the vertex is on an edge or corner, where the gradient takes one side)
/// a quarter of the way in from the corner to the triangle's centre, the
/// side the triangle lies on. Where the field has no gradient, the face's.
fn corner_normals(m: &Mesh, grads: &dyn Fn(&[V3]) -> Vec<[f64; 4]>) -> Vec<[f32; 3]> {
    let inner: Vec<V3> = m
        .tris
        .iter()
        .flat_map(|t| {
            let p = t.map(|i| m.verts[i as usize]);
            let c: V3 = std::array::from_fn(|k| (p[0][k] + p[1][k] + p[2][k]) / 3.0);
            p.map(|v| std::array::from_fn(|k| v[k] + 0.25 * (c[k] - v[k])))
        })
        .collect();
    let (at_corner, at_vertex) = (grads(&inner), grads(&m.verts));
    let unit = |g: [f64; 4]| {
        let v = [g[1], g[2], g[3]];
        let l = ink_core::v3::norm(v);
        (l > 0.0 && l.is_finite()).then(|| v.map(|x| x / l))
    };
    let vertex: Vec<Option<V3>> = at_vertex.into_iter().map(unit).collect();
    m.tris
        .iter()
        .enumerate()
        .flat_map(|(f, t)| {
            let face = ink_core::v3::normal(t.map(|i| m.verts[i as usize]));
            let face = unit([0.0, face[0], face[1], face[2]]);
            let (at_corner, vertex) = (&at_corner, &vertex);
            (0..3).map(move |c| {
                let own = unit(at_corner[3 * f + c]).or(face);
                let n = match (vertex[t[c] as usize], own) {
                    (Some(v), Some(o)) if ink_core::v3::dot(v, o) >= SMOOTH_COS => Some(v),
                    (v, o) => o.or(v),
                };
                n.map_or([0.0, 0.0, 1.0], |n| n.map(|x| x as f32))
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph::{Graph, Var};

    /// A unit cube, |p|_inf - 1, as a tape.
    fn cube() -> Arc<Tape> {
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let [ax, ay, az] = [x, y, z].map(|v| {
            let n = g.neg(v);
            g.max(v, n)
        });
        let m = g.max(ax, ay);
        let m = g.max(m, az);
        let one = g.constant(1.0);
        let f = g.sub(m, one);
        Arc::new(Tape::new(&g, f))
    }

    #[test]
    fn a_cube_is_shaded_flat_to_its_edges() {
        // Every corner normal is its face's axis exactly, also at the
        // vertices on the cube's edges and corners (where the field's
        // gradient there would take one side).
        let how = How { res: 0.2, tol: None, tol_per_res: 0.25, bend: None };
        let v = make(&cube(), &[], &Ctx::new(0.0, vec![]), [-1.2; 3], [1.2; 3], how);
        assert_eq!(v.normals.len(), 3 * v.tris.len());
        for (f, t) in v.tris.iter().enumerate() {
            let c: V3 = std::array::from_fn(|k| t.iter().map(|&i| v.verts[i as usize][k]).sum::<f64>() / 3.0);
            let axis = (0..3).max_by(|&a, &b| c[a].abs().total_cmp(&c[b].abs())).unwrap();
            for n in &v.normals[3 * f..3 * f + 3] {
                let want = c[axis].signum() as f32;
                assert!((n[axis] - want).abs() < 1e-5 && n.iter().map(|x| x * x).sum::<f32>() > 0.999, "{n:?} at {c:?}");
            }
        }
    }

    #[test]
    fn a_change_of_paint_is_cut_where_it_steps() {
        // A cube with a round label on its x faces: paint 1 inside y^2 + z^2
        // < 0.25, 0 outside. No triangle spans the two paints (its corners
        // and centre have one; the cut is a chord of the circle, within its
        // sagitta of it), yet away from the label the faces are still
        // simplified to large triangles.
        let tape = cube();
        let mut g = Graph::new();
        let (y, z) = (g.var(Var::Y), g.var(Var::Z));
        let (yy, zz) = (g.mul(y, y), g.mul(z, z));
        let r2 = g.add(yy, zz);
        let (quarter, k, one, zero) = (g.constant(0.25), g.constant(1e12), g.constant(1.0), g.constant(0.0));
        let d = g.sub(quarter, r2);
        let s = g.mul(d, k);
        let s = g.max(s, zero);
        let s = g.min(s, one);
        let paint = vec![Arc::new(Tape::new(&g, s))];
        let how = How { res: 0.1, tol: None, tol_per_res: 0.25, bend: None };
        let ctx = Ctx::new(0.0, vec![]);
        let v = make(&tape, &paint, &ctx, [-1.2; 3], [1.2; 3], how);
        let p = |q: V3| ink_core::eval::eval_point(&paint[0], &ctx, q);
        let mut cut = 0;
        for t in &v.tris {
            let c = t.map(|i| v.verts[i as usize]);
            let centre: V3 = std::array::from_fn(|k| (c[0][k] + c[1][k] + c[2][k]) / 3.0);
            let at = [c[0], c[1], c[2], centre].map(p);
            assert!(at.iter().all(|&q| (q - at[0]).abs() <= PAINT_STEP), "a triangle spans the label's edge: {at:?} at {c:?}");
            cut += c.iter().any(|q| (q[1] * q[1] + q[2] * q[2] - 0.25).abs() < 1e-5) as usize;
        }
        assert!(cut > 100, "{cut} triangles meet at the label's edge");
        let long = |t: &[u32; 3]| (0..3).map(|e| ink_core::v3::norm(ink_core::v3::sub(v.verts[t[e] as usize], v.verts[t[(e + 1) % 3] as usize]))).fold(0.0, f64::max);
        assert!(v.tris.iter().any(|t| long(t) > 1.0), "simplified");
        // and without paint, the cube's faces are a few triangles each
        assert!(make(&tape, &[], &ctx, [-1.2; 3], [1.2; 3], how).tris.len() < v.tris.len() / 4);
    }

    #[test]
    fn a_stroke_finer_than_res_is_drawn() {
        // A cube at res 0.2 with a stripe 0.07 wide (a third of res) painted
        // across its faces at y in [0.03, 0.10], lying along the cells' rows
        // and between their samples: without splitting the edges it crosses
        // no triangle has the stripe's paint and it vanishes. Its area on
        // the +x face (2 x 0.07) is drawn, sharp: no triangle spans two paints.
        let tape = cube();
        let mut g = Graph::new();
        let y = g.var(Var::Y);
        let (c, w, k, one, zero) = (g.constant(0.065), g.constant(0.035), g.constant(1e12), g.constant(1.0), g.constant(0.0));
        let d = g.sub(y, c);
        let n = g.neg(d);
        let d = g.max(d, n);
        let d = g.sub(w, d);
        let s = g.mul(d, k);
        let s = g.max(s, zero);
        let s = g.min(s, one);
        let paint = vec![Arc::new(Tape::new(&g, s))];
        let how = How { res: 0.2, tol: None, tol_per_res: 0.25, bend: None };
        let ctx = Ctx::new(0.0, vec![]);
        let v = make(&tape, &paint, &ctx, [-1.2; 3], [1.2; 3], how);
        let p = |q: V3| ink_core::eval::eval_point(&paint[0], &ctx, q);
        let mut area = 0.0;
        for t in &v.tris {
            let c = t.map(|i| v.verts[i as usize]);
            let centre: V3 = std::array::from_fn(|k| (c[0][k] + c[1][k] + c[2][k]) / 3.0);
            let at = [c[0], c[1], c[2], centre].map(p);
            assert!(at.iter().all(|&q| (q - at[0]).abs() <= PAINT_STEP), "a triangle spans the stripe's edge: {at:?} at {c:?}");
            if centre[0] > 0.99 && at[3] > 0.5 {
                area += 0.5 * ink_core::v3::norm(ink_core::v3::normal(c));
            }
        }
        assert!((area - 0.14).abs() < 0.01, "the stripe's area on +x {area}, not 0.14");
    }

    #[test]
    fn paint_among_many_parts_is_read_near_each_part_only() {
        // A ring of 16 balls 1.15 apart, each its own colour, painted as
        // catalog parts are (each channel the nearest part's, a chain of
        // selects over every part's field: the paint's tapes read all 16).
        // Meshing reads
        // each paint on its tape simplified to the parts near each block,
        // and at a few points per triangle sampled; a mesh reading every
        // part at every point (paint on whole tapes), or each edge's ends
        // again per edge and pass, took a mecanum wheel from 37 to 80 s.
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let ball = |g: &mut Graph, k: usize| {
            let a = k as f64 * std::f64::consts::TAU / 16.0;
            let d = [(x, 6.0 * a.cos()), (y, 6.0 * a.sin()), (z, 0.0)].map(|(v, c)| {
                let c = g.constant(c);
                let d = g.sub(v, c);
                g.mul(d, d)
            });
            let r2 = g.add(d[0], d[1]);
            let r2 = g.add(r2, d[2]);
            let r = g.un(ink_core::graph::UnOp::Sqrt, r2);
            let radius = g.constant(0.6);
            g.sub(r, radius)
        };
        let balls: Vec<_> = (0..16).map(|k| ball(&mut g, k)).collect();
        let colour = |k: usize, c: usize| ((k * 7 + c * 3) % 16) as f64 / 15.0;
        let mut channels: Vec<_> = (0..5).map(|c| g.constant(colour(0, c))).collect();
        let mut near = balls[0];
        for (k, &f) in balls.iter().enumerate().skip(1) {
            let d = g.sub(f, near);
            for (c, ch) in channels.iter_mut().enumerate() {
                let v = g.constant(colour(k, c));
                *ch = g.select(d, v, *ch);
            }
            near = g.select(d, f, near);
        }
        let shape = balls.iter().copied().reduce(|a, b| g.min(a, b)).unwrap();
        let tape = Arc::new(Tape::new(&g, shape));
        let paint: Vec<Arc<Tape>> = channels.iter().map(|&c| Arc::new(Tape::new(&g, c))).collect();
        let how = How { res: 0.1, tol: None, tol_per_res: 0.25, bend: None };
        let (v, [read, work]) = made(&tape, &paint, &Ctx::new(0.0, vec![]), [-6.8, -6.8, -0.8], [6.8, 6.8, 0.8], how);
        let per_value = work as f64 / read as f64;
        assert!(per_value < paint[0].len() as f64 / 4.0, "{per_value} instructions a paint value, of the paint's {}", paint[0].len());
        let per_triangle = read as f64 / (paint.len() * v.uniform_triangles) as f64;
        assert!(per_triangle < 6.0, "{per_triangle} points read a triangle sampled");
    }

    #[test]
    fn a_surface_has_the_paint_beneath_it() {
        // A block standing out of a cube's +x face, one body, painted 1
        // inside the block and 0 outside: the block's own surface, where its
        // field is zero, is all 1 and the cube's all 0 (read at the surface
        // the paint would flip from vertex to vertex on the block, whose
        // vertices lie a little either side of it). Not simplified, so
        // every vertex as sampled.
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let abs = |g: &mut Graph, v| {
            let n = g.neg(v);
            g.max(v, n)
        };
        let (one, half, k, zero) = (g.constant(1.0), g.constant(0.5), g.constant(1e12), g.constant(0.0));
        let [ax, ay, az] = [x, y, z].map(|v| abs(&mut g, v));
        let m = g.max(ax, ay);
        let m = g.max(m, az);
        let cube = g.sub(m, one);
        let bx = g.sub(x, one);
        let bx = abs(&mut g, bx);
        let m = g.max(bx, ay);
        let m = g.max(m, az);
        let block = g.sub(m, half);
        let body = g.min(cube, block);
        let s = g.neg(block);
        let s = g.mul(s, k);
        let s = g.max(s, zero);
        let s = g.min(s, one);
        let (tape, block, paint) = (Arc::new(Tape::new(&g, body)), Tape::new(&g, block), vec![Arc::new(Tape::new(&g, s))]);
        let how = How { res: 0.1, tol: None, tol_per_res: 0.0, bend: None };
        let ctx = Ctx::new(0.0, vec![]);
        let v = make(&tape, &paint, &ctx, [-1.2, -1.2, -1.2], [1.7, 1.2, 1.2], how);
        let on_block = v.verts.iter().zip(&v.paint_at).filter(|(q, _)| ink_core::eval::eval_point(&block, &ctx, **q).abs() < 0.02 && q[0] > 1.05);
        let mut n = 0;
        for (q, at) in on_block {
            assert_eq!(ink_core::eval::eval_point(&paint[0], &ctx, *at), 1.0, "the block's surface at {q:?}");
            n += 1;
        }
        assert!(n > 50, "{n} vertices on the block");
        // every vertex read in the material, at edges and corners too
        for (q, at) in v.verts.iter().zip(&v.paint_at) {
            assert!(ink_core::eval::eval_point(&tape, &ctx, *at) < 0.0, "{q:?} read outside, at {at:?}");
        }
        for (q, at) in v.verts.iter().zip(&v.paint_at).filter(|(q, _)| ink_core::eval::eval_point(&block, &ctx, **q) > 0.02) {
            assert_eq!(ink_core::eval::eval_point(&paint[0], &ctx, *at), 0.0, "the cube's surface at {q:?}");
        }
    }

    /// 40 x 40 plates at res 2, z from `lo` to `hi` each, as one body, the
    /// first painted 1 and the rest 0: its vertices on the faces (within
    /// 19 of the axis), by side (below z = 0, above), each checked to read 1.
    fn a_plate_reads_its_own_paint(plates: &[(f64, f64)]) -> [usize; 2] {
        let mut g = Graph::new();
        let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
        let abs = |g: &mut Graph, v| {
            let n = g.neg(v);
            g.max(v, n)
        };
        let [ax, ay] = [x, y].map(|v| abs(&mut g, v));
        let (twenty, k, zero, one) = (g.constant(20.0), g.constant(1e12), g.constant(0.0), g.constant(1.0));
        let m = g.max(ax, ay);
        let side = g.sub(m, twenty);
        let plates: Vec<_> = plates
            .iter()
            .map(|&(lo, hi)| {
                let (c, r) = (g.constant(0.5 * (lo + hi)), g.constant(0.5 * (hi - lo)));
                let d = g.sub(z, c);
                let d = abs(&mut g, d);
                let d = g.sub(d, r);
                g.max(side, d)
            })
            .collect();
        let body = plates[1..].iter().fold(plates[0], |b, &p| g.min(b, p));
        let s = g.neg(plates[0]);
        let s = g.mul(s, k);
        let s = g.max(s, zero);
        let s = g.min(s, one);
        let (tape, paint) = (Arc::new(Tape::new(&g, body)), vec![Arc::new(Tape::new(&g, s))]);
        let how = How { res: 2.0, tol: None, tol_per_res: 0.0, bend: None };
        let ctx = Ctx::new(0.0, vec![]);
        let v = make(&tape, &paint, &ctx, [-22.0, -22.0, -4.0], [22.0, 22.0, 4.0], how);
        let mut faces = [0, 0];
        for (q, at) in v.verts.iter().zip(&v.paint_at).filter(|(q, _)| q[0].abs().max(q[1].abs()) < 19.0 && q[2].abs() < 0.1) {
            assert_eq!(ink_core::eval::eval_point(&paint[0], &ctx, *at), 1.0, "the plate's surface at {q:?} read at {at:?}");
            faces[(q[2] > 0.0) as usize] += 1;
        }
        faces
    }

    #[test]
    fn a_sheet_thinner_than_the_step_or_beside_another_solid_has_its_own_paint() {
        // 0.1 thick: both faces read inside it, not past the far face.
        let faces = a_plate_reads_its_own_paint(&[(-0.05, 0.05)]);
        assert!(faces.iter().all(|&n| n > 100), "{faces:?} vertices on the faces");
        // 0.01 thick, and another plate from 0.16 up (the gap between them
        // finer than the cells, so the sheet's top face is not meshed): its
        // bottom face reads inside it, not in the other plate.
        let faces = a_plate_reads_its_own_paint(&[(-0.005, 0.005), (0.16, 2.0)]);
        assert!(faces[0] > 100, "{faces:?} vertices on the faces");
    }
}
