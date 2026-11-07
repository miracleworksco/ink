//! Simplification of a closed triangle mesh to a distance tolerance: large
//! triangles where the surface is flat, small ones where it curves.
//!
//! Half-edge collapses: a vertex moves onto a neighbour and its triangles
//! follow, so the output's vertices are some of the input's (on the surface
//! wherever the input's are). Collapses are tried in order of quadric error
//! (Garland and Heckbert 1997) and kept only if:
//!
//! - the edge's link condition holds, and both triangles' far vertices keep
//!   three neighbours or more: the surface keeps its topology (a closed
//!   2-manifold stays one, its shells and holes as they were);
//! - no triangle turns over or collapses to a sliver of no area;
//! - every input vertex removed so far stays within `tol` of the output.
//!   Each removed vertex is tracked to one output triangle, and its distance
//!   to that triangle (an upper bound on its distance to the output) is
//!   rechecked whenever the triangle changes;
//! - given a field (negative inside), the surface lies within `tol` of each
//!   new triangle's centroid and each new edge's midpoint, along the
//!   normal: the field changes sign there, from inside behind to outside in
//!   front. This is the other direction, from the output to the surface:
//!   between input vertices no triangle strays (along a long thin
//!   triangle's edge, or bridging a concave corner, where every input
//!   vertex is near some triangle but the triangle is not near the
//!   surface).
//!
//! Sharp edges and corners stay: moving a vertex off a crease takes the
//! crease's input vertices out of tolerance.
//!
//! Vertices on an edge with other than two triangles, or whose triangles do
//! not form one closed fan, or that the caller marks (`keep`), are kept as
//! they are; so is a vertex whose
//! every collapse the field refuses (around it the input already strays
//! from the surface by more than tol: a curve too tight for the sampling,
//! a feature under it), though a neighbour may still collapse onto it.
//!
//! Passes run on all cores. Every vertex proposes its cheapest collapse
//! that keeps the topology and turns no triangle over, against the mesh as
//! it stands. A proposal is chosen when its cost is the least among those
//! whose neighbourhoods (the triangles of both ends) overlap its own (then
//! again among the rest, a few rounds), so the chosen ones are disjoint and
//! still valid when applied together. The chosen ones are then checked for
//! distance: tracked vertices, and the field, evaluated for all at once. A
//! refused collapse stays refused (the field's verdict while the vertex's
//! triangles stay, a tracked vertex's while its neighbourhood does), and
//! meanwhile its vertex proposes every collapse that passes its tracked
//! vertices, the field choosing the cheapest it passes. The next pass
//! proposes again only next to a change.

use ink_core::V3;
use ink_core::sign::inside;
use ink_core::v3::{direction, dot, norm, normal, sub};
use ink_kernel::mesh::Mesh;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

const NONE: u32 = u32::MAX;
const DEAD: [u32; 3] = [NONE; 3];

pub struct Simplified {
    pub mesh: Mesh,
    /// For each input vertex, its distance to the output: 0 where it was
    /// kept, else to the output triangle it is tracked to (an upper bound).
    /// Each is at most the tolerance.
    pub moved: Vec<f64>,
}

/// Values of a field at points: negative inside, positive outside.
pub type Field<'a> = &'a (dyn Fn(&[V3]) -> Vec<f64> + Sync);

/// Simplify `m` so that every one of its vertices lies within `tol` of the
/// result and, given `field`, each new triangle's centroid and new edge's
/// midpoint within `tol` of the field's surface (see the module notes).
/// `tol <= 0` returns the mesh unchanged. Given `along_x`, no collapse
/// makes an edge longer than that along x (a strip that will be bent
/// along x: a carried body's). The vertices `keep` marks (by index; it may
/// be empty) stay as they are, like those on open edges: none moves, and
/// none is collapsed onto (where a body's colour changes, so the change
/// stays where it was sampled).
pub fn simplify(m: &Mesh, tol: f64, field: Option<Field>, along_x: Option<f64>, keep: &[bool]) -> Simplified {
    if !(tol > 0.0) {
        return Simplified { mesh: m.clone(), moved: vec![0.0; m.verts.len()] };
    }
    let mut s = State::new(m);
    s.fixed.iter_mut().zip(keep).for_each(|(f, k)| *f |= k);
    s.along_x = along_x.unwrap_or(f64::INFINITY);
    let grid = s.grid();
    for shift in [0.0, 0.5] {
        s.by_blocks(tol, field, grid, shift);
    }
    s.grain = GRAIN;
    s.run(tol, field);
    s.finish()
}

/// A block of space simplified alone is about this many edges across.
const BLOCK_EDGES: f64 = 48.0;

/// The fewest items a pass splits off to another thread.
const GRAIN: usize = 1024;

/// A block simplified alone: its vertices and triangles (the whole mesh's
/// numbers), and its state after, numbered as those lists are.
struct Block {
    verts: Vec<u32>,
    tris: Vec<u32>,
    out: Vec<[u32; 3]>,
    quad: Vec<Quadric>,
    stays: Vec<bool>,
    owner: Vec<u32>,
    moved: Vec<f64>,
}

/// A symmetric 4x4 quadric: p^T A p + 2 b.p + c, as
/// [a00 a01 a02 a11 a12 a22 b0 b1 b2 c].
#[derive(Clone, Copy, Default)]
struct Quadric([f64; 10]);

impl Quadric {
    /// Squared distance to the plane n.p + d = 0 (|n| = 1), times w.
    fn plane(n: V3, d: f64, w: f64) -> Quadric {
        Quadric([w * n[0] * n[0], w * n[0] * n[1], w * n[0] * n[2], w * n[1] * n[1], w * n[1] * n[2], w * n[2] * n[2], w * d * n[0], w * d * n[1], w * d * n[2], w * d * d])
    }

    fn add(&mut self, o: &Quadric) {
        for k in 0..10 {
            self.0[k] += o.0[k];
        }
    }

    fn eval(&self, p: V3) -> f64 {
        let q = &self.0;
        let [x, y, z] = p;
        q[0] * x * x + q[3] * y * y + q[5] * z * z + 2.0 * (q[1] * x * y + q[2] * x * z + q[4] * y * z) + 2.0 * (q[6] * x + q[7] * y + q[8] * z) + q[9]
    }
}

/// Whether a triangle has area for its size: |n| > 1e-8 of its longest edge
/// squared (a sliver's normal means nothing).
pub fn has_area(t: [V3; 3], n: V3) -> bool {
    let l2 = (0..3).map(|k| dot(sub(t[(k + 1) % 3], t[k]), sub(t[(k + 1) % 3], t[k]))).fold(0.0, f64::max);
    dot(n, n) > 1e-16 * l2 * l2
}

/// A triangle set up for many distance queries: what depends only on the
/// triangle, computed once.
struct Near {
    a: V3,
    ab: V3,
    ac: V3,
    bc: V3,
    /// ab.ab, ab.ac, ac.ac, bc.bc, and 1 / (ab.ab ac.ac - (ab.ac)^2).
    d: [f64; 5],
}

impl Near {
    fn new([a, b, c]: [V3; 3]) -> Near {
        let (ab, ac, bc) = (sub(b, a), sub(c, a), sub(c, b));
        let (aa, ax, cc) = (dot(ab, ab), dot(ab, ac), dot(ac, ac));
        Near { a, ab, ac, bc, d: [aa, ax, cc, dot(bc, bc), 1.0 / (aa * cc - ax * ax)] }
    }

    /// Squared distance from p: the plane's foot if it falls inside, else
    /// the nearest edge. A sliver (no inverse) gives NaN, out of tolerance.
    fn dist2(&self, p: V3) -> f64 {
        let ap = sub(p, self.a);
        let (e, f) = (dot(ap, self.ab), dot(ap, self.ac));
        let [aa, ax, cc, bb, inv] = self.d;
        let v = (cc * e - ax * f) * inv;
        let w = (aa * f - ax * e) * inv;
        if v >= 0.0 && w >= 0.0 && v + w <= 1.0 {
            let q = [ap[0] - v * self.ab[0] - w * self.ac[0], ap[1] - v * self.ab[1] - w * self.ac[1], ap[2] - v * self.ab[2] - w * self.ac[2]];
            return dot(q, q);
        }
        let seg = |o: V3, d: V3, dd: f64| {
            let t = (dot(o, d) / dd).clamp(0.0, 1.0);
            let q = [o[0] - t * d[0], o[1] - t * d[1], o[2] - t * d[2]];
            dot(q, q)
        };
        let bp = sub(ap, self.ab);
        seg(ap, self.ab, aa).min(seg(ap, self.ac, cc)).min(seg(bp, self.bc, bb))
    }
}

/// For each list of probes, whether the field passes them all (`probes`).
fn passes(field: Field, probes: &[Vec<[V3; 2]>]) -> Vec<bool> {
    let f = field(&probes.iter().flatten().flatten().copied().collect::<Vec<_>>());
    let mut at = 0;
    probes
        .iter()
        .map(|p| {
            let ok = (0..p.len()).all(|k| inside(f[at + 2 * k]) && !inside(f[at + 2 * k + 1]));
            at += 2 * p.len();
            ok
        })
        .collect()
}

/// A proposed collapse of u onto v; `checked` if its distances (tracked
/// vertices and field) have been checked already.
#[derive(Clone, Copy)]
struct Candidate {
    cost: f64,
    u: u32,
    v: u32,
    checked: bool,
}

/// A vertex's proposal: its cheapest collapse that keeps the topology and
/// turns no triangle over; or, once one of its collapses has been refused,
/// every such collapse that also keeps its tracked vertices within tol, in
/// cost order, with the field's probes for each (the field decides among
/// them, for all such vertices at once).
enum Proposal {
    Shape(Option<Candidate>),
    Checked(Vec<(Candidate, Vec<[V3; 2]>)>),
}

impl Candidate {
    /// Orders candidates by cost, then in a scrambled order of vertices
    /// (a bijection, so keys stay distinct): vertex numbers follow space,
    /// and in that order, among equal costs, few would be least locally.
    fn key(&self) -> u64 {
        let mut h = self.u.wrapping_mul(0x9E37_79B1);
        h ^= h >> 15;
        h = h.wrapping_mul(0x85EB_CA77);
        h ^= h >> 13;
        ((self.cost as f32).to_bits() as u64) << 32 | h as u64
    }
}

/// What collapsing u onto v changes: triangles rewritten (or DEAD), and the
/// tracked vertices it reassigns, with their distances.
#[derive(Default)]
struct Plan {
    u: u32,
    v: u32,
    rewrite: Vec<(u32, [u32; 3])>,
    assign: Vec<(u32, u32, f64)>,
}

struct State<'a> {
    pos: &'a [V3],
    tris: Vec<[u32; 3]>,
    quad: Vec<Quadric>,
    fixed: Vec<bool>,
    /// Vertices whose every collapse the field refused.
    stays: Vec<bool>,
    /// Each vertex's triangles.
    vt: Vec<Vec<u32>>,
    /// Removed vertices tracked to each triangle: a list through `next`.
    head: Vec<u32>,
    next: Vec<u32>,
    /// For a removed vertex: the triangle it is tracked to, and how far.
    owner: Vec<u32>,
    moved: Vec<f64>,
    /// The longest an edge may be along x.
    along_x: f64,
    /// The fewest items a pass hands to another thread (usize::MAX: none).
    grain: usize,
}

impl<'a> State<'a> {
    fn new(m: &'a Mesh) -> State<'a> {
        let nv = m.verts.len();
        let tris = m.tris.clone();
        let mut quad = vec![Quadric::default(); nv];
        let mut vt = vec![Vec::new(); nv];
        for (ti, t) in tris.iter().enumerate() {
            let p = t.map(|i| m.verts[i as usize]);
            let n = normal(p);
            let l = norm(n);
            let q = (l > 0.0).then(|| {
                let nn = n.map(|c| c / l);
                Quadric::plane(nn, -dot(nn, p[0]), 0.5 * l)
            });
            for &i in t {
                if let Some(q) = &q {
                    quad[i as usize].add(q);
                }
                vt[i as usize].push(ti as u32);
            }
        }
        let fixed = (0..nv).into_par_iter().map(|u| !closed_fan(&tris, &vt[u], u as u32)).collect();
        let nt = tris.len();
        State {
            pos: &m.verts,
            tris,
            quad,
            fixed,
            stays: vec![false; nv],
            vt,
            head: vec![NONE; nt],
            next: vec![NONE; nv],
            owner: vec![NONE; nv],
            moved: vec![0.0; nv],
            along_x: f64::INFINITY,
            grain: usize::MAX,
        }
    }

    /// The blocks of space the mesh is simplified in first: a corner, and
    /// a side about BLOCK_EDGES of its edges across.
    fn grid(&self) -> (V3, f64) {
        let mut lo = [f64::INFINITY; 3];
        let mut edges = 0.0;
        for t in &self.tris {
            let [a, b] = [t[0], t[1]].map(|i| self.pos[i as usize]);
            (0..3).for_each(|c| lo[c] = lo[c].min(a[c]));
            edges += norm(sub(a, b));
        }
        (lo, BLOCK_EDGES * edges / self.tris.len().max(1) as f64)
    }

    /// The mesh cut into blocks of space (`grid`, moved `shift` of a block
    /// along each axis), each simplified alone on one thread, all at once,
    /// then taken back as it left them. A vertex whose triangles lie in
    /// more than one block has an open fan in each and stays; the blocks
    /// moved half a block take those between the first ones, and the whole
    /// mesh's pass (`run`) those at their corners. A pass over the whole
    /// mesh waits on every thread at each of its steps, a hundred passes
    /// for a fine mesh: slower on many threads than on one, and slower
    /// still on a loaded machine.
    fn by_blocks(&mut self, tol: f64, field: Option<Field>, (lo, side): (V3, f64), shift: f64) {
        if !(side > 0.0) {
            return;
        }
        let key = |t: u32| -> [i64; 3] {
            let p = self.tri(self.tris[t as usize]);
            std::array::from_fn(|c| (((p[0][c] + p[1][c] + p[2][c]) / 3.0 - lo[c]) / side + shift).floor() as i64)
        };
        let mut keyed: Vec<([i64; 3], u32)> = (0..self.tris.len() as u32).into_par_iter().filter(|&t| self.tris[t as usize] != DEAD).map(|t| (key(t), t)).collect();
        keyed.par_sort_unstable();
        let blocks: Vec<&[([i64; 3], u32)]> = keyed.chunk_by(|a, b| a.0 == b.0).collect();
        if blocks.len() < 2 {
            return;
        }
        let done: Vec<Block> = blocks.par_iter().with_min_len(1).map(|b| self.alone(&b.iter().map(|x| x.1).collect::<Vec<_>>(), tol, field)).collect();
        for b in done {
            for (k, &t) in b.tris.iter().enumerate() {
                self.tris[t as usize] = match b.out[k] {
                    DEAD => DEAD,
                    o => o.map(|w| b.verts[w as usize]),
                };
            }
            for (k, &w) in b.verts.iter().enumerate() {
                let w = w as usize;
                (self.quad[w], self.stays[w], self.moved[w]) = (b.quad[k], b.stays[k], b.moved[k]);
                if b.owner[k] != NONE {
                    self.owner[w] = b.tris[b.owner[k] as usize];
                }
            }
        }
        self.vt.iter_mut().for_each(Vec::clear);
        for (t, tri) in self.tris.iter().enumerate().filter(|t| *t.1 != DEAD) {
            tri.iter().for_each(|&w| self.vt[w as usize].push(t as u32));
        }
        self.head.fill(NONE);
        for w in 0..self.pos.len() {
            let t = self.owner[w];
            if t != NONE {
                self.next[w] = self.head[t as usize];
                self.head[t as usize] = w as u32;
            }
        }
    }

    /// The triangles `tris` simplified alone, on this thread, with what is
    /// tracked to them: a state of their own, as this one has them.
    fn alone(&self, tris: &[u32], tol: f64, field: Option<Field>) -> Block {
        let tracked = |t: u32| std::iter::successors(Some(self.head[t as usize]).filter(|&w| w != NONE), |&w| Some(self.next[w as usize]).filter(|&w| w != NONE));
        let mut verts: Vec<u32> = tris.iter().flat_map(|&t| self.tris[t as usize].into_iter().chain(tracked(t))).collect();
        verts.sort_unstable();
        verts.dedup();
        let local = |w: u32| verts.binary_search(&w).unwrap() as u32;
        let m = Mesh { verts: verts.iter().map(|&w| self.pos[w as usize]).collect(), tris: tris.iter().map(|&t| self.tris[t as usize].map(local)).collect() };
        let mut s = State::new(&m);
        for (k, &w) in verts.iter().enumerate() {
            let w = w as usize;
            s.fixed[k] |= self.fixed[w];
            (s.quad[k], s.stays[k], s.moved[k]) = (self.quad[w], self.stays[w], self.moved[w]);
        }
        for (k, &t) in tris.iter().enumerate() {
            for w in tracked(t).map(local) {
                (s.owner[w as usize], s.next[w as usize], s.head[k]) = (k as u32, s.head[k], w);
            }
        }
        s.along_x = self.along_x;
        s.run(tol, field);
        let State { tris: out, quad, stays, owner, moved, .. } = s;
        Block { verts, tris: tris.to_vec(), out, quad, stays, owner, moved }
    }

    fn run(&mut self, tol: f64, field: Option<Field>) {
        let (nv, g) = (self.pos.len(), self.grain);
        let mut cand: Vec<Option<Candidate>> = vec![None; nv];
        let mut dirty: Vec<u32> = (0..nv as u32).filter(|&u| !self.vt[u as usize].is_empty()).collect();
        let mut active: Vec<u32> = Vec::new();
        let mut is_active = vec![false; nv];
        // Collapses refused, by vertex. The field's verdict reads only the
        // positions of the vertex's neighbours and its target, so it stands
        // while the vertex's triangles do; a tracked vertex's, also what is
        // tracked to them and the target's triangles, so it stands until
        // its neighbourhood changes.
        let mut off_surface: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut off_track: HashMap<u32, Vec<u32>> = HashMap::new();
        let slot: Vec<AtomicU64> = (0..self.tris.len()).map(|_| AtomicU64::new(u64::MAX)).collect();
        let mark: Vec<AtomicBool> = (0..nv).map(|_| AtomicBool::new(false)).collect();
        let touched: Vec<AtomicBool> = (0..nv).map(|_| AtomicBool::new(false)).collect();
        while !dirty.is_empty() {
            let st = &*self;
            let fresh: Vec<Proposal> = dirty
                .par_iter()
                .with_min_len(g)
                .map_init(Plan::default, |scratch, &u| {
                    let refused = [&off_surface, &off_track].map(|m| m.get(&u).map_or(&[][..], |b| &b[..]));
                    st.propose(u, tol, refused, scratch)
                })
                .collect();
            // The field decides among each refused vertex's collapses: the
            // cheapest it passes; those it fails stay refused.
            let lists: Vec<Vec<[V3; 2]>> = fresh
                .iter()
                .flat_map(|p| match p {
                    Proposal::Checked(l) => l.iter().map(|c| c.1.clone()).collect(),
                    Proposal::Shape(_) => Vec::new(),
                })
                .collect();
            let mut verdicts = match field {
                Some(f) => passes(f, &lists),
                None => vec![true; lists.len()],
            }
            .into_iter();
            for (&u, p) in dirty.iter().zip(fresh) {
                let c = match p {
                    Proposal::Shape(c) => c,
                    Proposal::Checked(l) => {
                        let mut first = None;
                        let all = !l.is_empty();
                        for (c, _) in l {
                            match verdicts.next().unwrap() {
                                true => first = first.or(Some(c)),
                                false => off_surface.entry(u).or_default().push(c.v),
                            }
                        }
                        // The field refuses every collapse of u: the surface
                        // strays from u's triangles more than tol already,
                        // and coarser ones only stray further. u stays (a
                        // neighbour may still collapse onto it).
                        if all && first.is_none() {
                            self.stays[u as usize] = true;
                        }
                        first
                    }
                };
                cand[u as usize] = c;
                if c.is_some() && !is_active[u as usize] {
                    is_active[u as usize] = true;
                    active.push(u);
                }
            }
            active.retain(|&u| {
                let keep = cand[u as usize].is_some();
                is_active[u as usize] = keep;
                keep
            });

            // Apply each collapse whose key is least on every triangle of
            // its neighbourhood; then again among those whose neighbourhood
            // no chosen one touches, a few rounds.
            let st = &*self;
            let region = |u: u32| {
                let c = cand[u as usize].unwrap();
                st.vt[c.u as usize].iter().chain(&st.vt[c.v as usize]).copied()
            };
            const TAKEN: u64 = 0;
            let mut chosen: Vec<Candidate> = Vec::new();
            let mut open = active.clone();
            for _ in 0..4 {
                open.par_iter().with_min_len(g).for_each(|&u| {
                    let k = cand[u as usize].unwrap().key();
                    region(u).for_each(|t| {
                        slot[t as usize].fetch_min(k, Ordering::Relaxed);
                    });
                });
                let won: Vec<u32> =
                    open.par_iter().with_min_len(g).copied().filter(|&u| region(u).all(|t| slot[t as usize].load(Ordering::Relaxed) == cand[u as usize].unwrap().key())).collect();
                won.par_iter().with_min_len(g).for_each(|&u| region(u).for_each(|t| slot[t as usize].store(TAKEN, Ordering::Relaxed)));
                open.par_iter().with_min_len(g).for_each(|&u| {
                    region(u).for_each(|t| {
                        let _ = slot[t as usize].try_update(Ordering::Relaxed, Ordering::Relaxed, |k| (k != TAKEN).then_some(u64::MAX));
                    });
                });
                chosen.extend(won.iter().map(|&u| cand[u as usize].unwrap()));
                open.retain(|&u| region(u).all(|t| slot[t as usize].load(Ordering::Relaxed) != TAKEN));
                if open.is_empty() {
                    break;
                }
            }
            chosen.par_iter().with_min_len(g).for_each(|c| region(c.u).for_each(|t| slot[t as usize].store(u64::MAX, Ordering::Relaxed)));
            if chosen.is_empty() {
                break;
            }
            // Plan the chosen ones; refuse those that would take a tracked
            // vertex, or (given the field) the surface, out of tolerance.
            let planned: Vec<(Plan, (bool, bool))> = chosen
                .par_iter()
                .with_min_len(g)
                .filter_map(|c| {
                    let mut p = Plan::default();
                    let valid = st.plan(c.u, c.v, &mut p);
                    debug_assert!(valid, "a proposed collapse stays valid while its neighbourhood does");
                    valid.then(|| {
                        let within = st.track(tol, &mut p);
                        (p, (c.checked, within))
                    })
                })
                .collect();
            let tracked: Vec<(u32, u32)> = planned.iter().filter(|p| !p.1.1).map(|p| (p.0.u, p.0.v)).collect();
            let mut plans: Vec<(Plan, bool)> = planned.into_iter().filter(|p| p.1.1).map(|p| (p.0, p.1.0)).collect();
            let surfaced = match field {
                Some(f) => self.refuse(&mut plans, tol, f),
                None => Vec::new(),
            };

            // What changed: the vertices of each collapsed vertex's triangles.
            let st = &*self;
            let changed: Vec<u32> = plans
                .par_iter()
                .with_min_len(g)
                .flat_map_iter(|(p, _)| st.vt[p.u as usize].iter().flat_map(|&t| st.tris[t as usize]))
                .filter(|&x| !touched[x as usize].swap(true, Ordering::Relaxed))
                .collect();
            for (p, _) in &plans {
                self.apply(p);
                cand[p.u as usize] = None;
            }
            // Propose again next to a change: a proposal reads its vertex's
            // triangles and its neighbours'. A refused collapse stays refused
            // until then.
            let st = &*self;
            let fresh = |x: u32| !mark[x as usize].swap(true, Ordering::Relaxed);
            let mut next: Vec<u32> =
                changed.par_iter().with_min_len(g).flat_map_iter(|&w| st.vt[w as usize].iter().flat_map(|&t| st.tris[t as usize])).filter(|&x| fresh(x)).collect();
            for &w in &changed {
                off_surface.remove(&w);
            }
            for &w in &next {
                off_track.remove(&w);
            }
            for (map, refused) in [(&mut off_surface, &surfaced), (&mut off_track, &tracked)] {
                for &(u, v) in refused {
                    map.entry(u).or_default().push(v);
                    if fresh(u) {
                        next.push(u);
                    }
                }
            }
            next.par_iter().with_min_len(g).for_each(|&w| mark[w as usize].store(false, Ordering::Relaxed));
            changed.par_iter().with_min_len(g).for_each(|&w| touched[w as usize].store(false, Ordering::Relaxed));
            dirty = next;
        }
    }

    /// The field's probes for a collapse: at each new triangle's centroid
    /// and each new edge's midpoint (where a long thin triangle strays
    /// most), the points tol behind and in front along the normal. It
    /// passes where the field is inside (f < 0, as for meshing) behind and
    /// outside in front: the surface lies within tol, the right way round.
    fn probes(&self, p: &Plan, tol: f64) -> Vec<[V3; 2]> {
        let pv = self.pos[p.v as usize];
        let mut at: Vec<(V3, V3)> = Vec::with_capacity(2 * p.rewrite.len());
        let mut edges: Vec<(u32, V3)> = Vec::new();
        for t in p.rewrite.iter().filter(|r| r.1 != DEAD).map(|r| r.1) {
            let c = self.tri(t);
            let n = normal(c);
            at.push((std::array::from_fn(|a| (c[0][a] + c[1][a] + c[2][a]) / 3.0), n));
            for &x in t.iter().filter(|&&x| x != p.v) {
                match edges.iter_mut().find(|e| e.0 == x) {
                    Some(e) => e.1 = std::array::from_fn(|a| e.1[a] + n[a]),
                    None => edges.push((x, n)),
                }
            }
        }
        at.extend(edges.iter().map(|&(x, n)| (std::array::from_fn(|a| 0.5 * (pv[a] + self.pos[x as usize][a])), n)));
        at.into_iter()
            .map(|(m, n)| {
                let s = tol / norm(n);
                [std::array::from_fn(|a| m[a] - s * n[a]), std::array::from_fn(|a| m[a] + s * n[a])]
            })
            .collect()
    }

    /// Drop the plans not yet checked whose probes the field fails.
    /// Returns the collapses dropped.
    fn refuse(&self, plans: &mut Vec<(Plan, bool)>, tol: f64, field: Field) -> Vec<(u32, u32)> {
        let g = self.grain;
        let probes: Vec<Vec<[V3; 2]>> = plans.par_iter().with_min_len(g).map(|(p, checked)| if *checked { Vec::new() } else { self.probes(p, tol) }).collect();
        let ok = passes(field, &probes);
        let mut refused = Vec::new();
        let mut keep = ok.into_iter();
        plans.retain(|(p, _)| {
            let k = keep.next().unwrap();
            if !k {
                refused.push((p.u, p.v));
            }
            k
        });
        refused
    }

    fn apply(&mut self, p: &Plan) {
        let (u, v) = (p.u as usize, p.v as usize);
        for &t in &self.vt[u] {
            self.head[t as usize] = NONE;
        }
        for &(t, tri) in &p.rewrite {
            if tri == DEAD {
                for w in self.tris[t as usize] {
                    if w != p.u {
                        self.vt[w as usize].retain(|&x| x != t);
                    }
                }
            } else {
                self.vt[v].push(t);
            }
            self.tris[t as usize] = tri;
        }
        self.vt[u].clear();
        for &(w, t, d) in &p.assign {
            self.owner[w as usize] = t;
            self.moved[w as usize] = d;
            self.next[w as usize] = self.head[t as usize];
            self.head[t as usize] = w;
        }
        let qu = self.quad[u];
        self.quad[v].add(&qu);
    }

    fn tri(&self, t: [u32; 3]) -> [V3; 3] {
        t.map(|i| self.pos[i as usize])
    }

    /// u's proposal (see `Proposal`), among its neighbours not refused
    /// (by the field, by a tracked vertex).
    fn propose(&self, u: u32, tol: f64, refused: [&[u32]; 2], scratch: &mut Plan) -> Proposal {
        if self.fixed[u as usize] || self.stays[u as usize] {
            return Proposal::Shape(None);
        }
        let q = self.quad[u as usize];
        let mut costs: Vec<(f64, u32)> = self.vt[u as usize]
            .iter()
            .flat_map(|&t| self.tris[t as usize])
            .filter(|&v| v != u && !self.fixed[v as usize] && !refused.iter().any(|r| r.contains(&v)))
            .map(|v| {
                let mut qv = q;
                qv.add(&self.quad[v as usize]);
                (qv.eval(self.pos[v as usize]).max(0.0), v)
            })
            .collect();
        costs.sort_unstable_by(|p, q| p.0.total_cmp(&q.0).then(p.1.cmp(&q.1)));
        costs.dedup_by_key(|c| c.1);
        if refused.iter().all(|r| r.is_empty()) {
            let c = costs.into_iter().find(|&(_, v)| self.plan(u, v, scratch));
            return Proposal::Shape(c.map(|(cost, v)| Candidate { cost, u, v, checked: false }));
        }
        let mut list = Vec::new();
        for (cost, v) in costs {
            if self.plan(u, v, scratch) && self.track(tol, scratch) {
                list.push((Candidate { cost, u, v, checked: true }, self.probes(scratch, tol)));
            }
        }
        Proposal::Checked(list)
    }

    /// The far vertices of edge uv's two triangles, if it has two.
    fn apexes(&self, u: u32, v: u32) -> Option<(u32, u32)> {
        let mut far = self.vt[u as usize].iter().map(|&t| self.tris[t as usize]).filter(|t| t.contains(&v)).map(|t| t.into_iter().find(|&w| w != u && w != v).unwrap());
        let (a, b) = (far.next()?, far.next()?);
        far.next().is_none().then_some((a, b))
    }

    /// Fill `p` with what collapsing u onto v changes to the triangles;
    /// false if that changes the topology or turns a triangle over (`track`
    /// checks the distances). Far vertices keep three neighbours or more even when
    /// several collapses apply at once: each takes 3 consecutive triangles
    /// of its far vertex's fan, which no other applied collapse touches.
    fn plan(&self, u: u32, v: u32, p: &mut Plan) -> bool {
        p.u = u;
        p.v = v;
        p.rewrite.clear();
        p.assign.clear();
        let (tu, tv) = (&self.vt[u as usize], &self.vt[v as usize]);
        let Some((a, b)) = self.apexes(u, v) else { return false };
        // Link condition: u and v share no neighbour but a and b.
        let shared = tu.iter().flat_map(|&t| self.tris[t as usize]).any(|x| x != u && x != v && x != a && x != b && tv.iter().any(|&t| self.tris[t as usize].contains(&x)));
        if shared || self.vt[a as usize].len() <= 3 || self.vt[b as usize].len() <= 3 {
            return false;
        }
        // No triangle turns over or loses its area. A sliver among u's input
        // triangles is compared with u's fan as a whole.
        let fan = tu.iter().fold([0.0; 3], |s, &t| {
            let n = normal(self.tri(self.tris[t as usize]));
            [s[0] + n[0], s[1] + n[1], s[2] + n[2]]
        });
        for &t in tu {
            let old = self.tris[t as usize];
            if old.contains(&v) {
                p.rewrite.push((t, DEAD));
                continue;
            }
            let new = old.map(|w| if w == u { v } else { w });
            let (po, pn) = (self.tri(old), self.tri(new));
            let (no, nn) = (normal(po), normal(pn));
            let reference = if has_area(po, no) { no } else { fan };
            if !has_area(pn, nn) || dot(reference, nn) <= 0.0 {
                return false;
            }
            if (0..3).any(|k| (pn[k][0] - pn[(k + 1) % 3][0]).abs() > self.along_x) {
                return false;
            }
            p.rewrite.push((t, new));
        }
        true
    }

    /// Track u, and every vertex tracked to u's triangles, for the collapse
    /// `p` (planned): false if one would lie farther than tol from v's fan.
    fn track(&self, tol: f64, p: &mut Plan) -> bool {
        let (u, v) = (p.u, p.v);
        let (tu, tv) = (&self.vt[u as usize], &self.vt[v as usize]);
        p.assign.clear();
        // u, and every vertex tracked to u's triangles, to a triangle of v's
        // fan as it will be within tol: its own triangle as rewritten if it
        // is; for one on a removed triangle (u, v, x), first those at x.
        let after: Vec<(u32, Near)> = tv
            .iter()
            .filter(|&&t| !self.tris[t as usize].contains(&u))
            .map(|&t| (t, Near::new(self.tri(self.tris[t as usize]))))
            .chain(p.rewrite.iter().filter(|r| r.1 != DEAD).map(|&(t, tri)| (t, Near::new(self.tri(tri)))))
            .collect();
        let tol2 = tol * tol;
        let track = |w: u32, first: &[usize], assign: &mut Vec<(u32, u32, f64)>| {
            let q = self.pos[w as usize];
            for k in first.iter().copied().chain(0..after.len()) {
                let d = after[k].1.dist2(q);
                if d <= tol2 {
                    assign.push((w, after[k].0, d.sqrt()));
                    return true;
                }
            }
            false
        };
        if !track(u, &[], &mut p.assign) {
            return false;
        }
        let mut first: Vec<usize> = Vec::with_capacity(4);
        for &t in tu {
            let mut w = self.head[t as usize];
            if w == NONE {
                continue;
            }
            first.clear();
            match after.iter().position(|a| a.0 == t) {
                Some(k) => first.push(k),
                None => {
                    let x = self.tris[t as usize].into_iter().find(|&x| x != u && x != v).unwrap();
                    let has_x = |k: usize| {
                        let t = after[k].0;
                        p.rewrite.iter().find(|r| r.0 == t).map_or(self.tris[t as usize], |r| r.1).contains(&x)
                    };
                    first.extend((0..after.len()).filter(|&k| has_x(k)));
                }
            }
            while w != NONE {
                if !track(w, &first, &mut p.assign) {
                    return false;
                }
                w = self.next[w as usize];
            }
        }
        true
    }

    fn finish(self) -> Simplified {
        let nv = self.pos.len();
        let mut index = vec![NONE; nv];
        let mut verts = Vec::new();
        let tris = self
            .tris
            .iter()
            .filter(|t| **t != DEAD)
            .map(|t| {
                t.map(|w| {
                    if index[w as usize] == NONE {
                        index[w as usize] = verts.len() as u32;
                        verts.push(self.pos[w as usize]);
                    }
                    index[w as usize]
                })
            })
            .collect();
        Simplified { mesh: Mesh { verts, tris }, moved: self.moved }
    }
}

/// Whether u's triangles form one closed fan, each of its edges shared by
/// exactly two of them, in opposite directions.
fn closed_fan(tris: &[[u32; 3]], tu: &[u32], u: u32) -> bool {
    if tu.len() < 3 {
        return false;
    }
    // Each triangle (u, p, q), in its winding, steps p -> q around u.
    let step: Vec<(u32, u32)> = tu
        .iter()
        .map(|&t| {
            let tri = tris[t as usize];
            let k = tri.iter().position(|&w| w == u).unwrap();
            (tri[(k + 1) % 3], tri[(k + 2) % 3])
        })
        .collect();
    let mut from: Vec<u32> = step.iter().map(|s| s.0).collect();
    from.sort_unstable();
    if from.windows(2).any(|w| w[0] == w[1]) {
        return false;
    }
    let (start, mut at, mut n) = (step[0].0, step[0].1, 1);
    while at != start {
        match step.iter().find(|s| s.0 == at) {
            Some(s) if n < step.len() => {
                at = s.1;
                n += 1;
            }
            _ => return false,
        }
    }
    n == step.len()
}

// ---------------------------------------------------------------- deviation

/// How far a mesh lies from a field's surface, by samples (`deviation`).
#[derive(Clone, Copy, Debug)]
pub struct Deviation {
    /// The largest bound over the samples, and the sample it is at.
    pub max: f64,
    pub at: V3,
    /// Samples with no sign change within `reach` along their normal (not
    /// in `max`): the mesh may lie farther than that from the surface there.
    pub beyond: usize,
    pub first_beyond: Option<V3>,
}

/// How far a mesh lies from the zero set of `field` (negative inside), at
/// samples at most `spacing` apart: its vertices, points along its edges,
/// and a grid over each triangle (rows parallel to its longest edge). Each
/// sample's distance is bounded above by a change of sign of the field
/// along a line through it: the mesh's normal there, and where that finds
/// none within a few `step`s, also the field's gradient (`grads`: value
/// and gradient), the nearer kept. Searched both ways in steps doubling from
/// `step` up to `reach`, then halved 8 times at the samples deciding the
/// largest (the rest stay at their bracket, within `step` or a doubling).
/// A sample, not a bound on the whole mesh: a mesh can stray between
/// samples, and surface far from the mesh (a missed feature) is not seen
/// from it.
pub fn deviation(field: Field, grads: &(dyn Fn(&[V3]) -> Vec<[f64; 4]> + Sync), m: &Mesh, spacing: f64, step: f64, reach: f64) -> Deviation {
    let area: Vec<V3> = m.tris.par_iter().map(|t| normal(t.map(|i| m.verts[i as usize]))).collect();
    // Vertex normals weighted by their triangles' angles there, so a vertex
    // where large and small triangles meet leans to neither.
    let mut vn = vec![[0.0; 3]; m.verts.len()];
    for (t, n) in m.tris.iter().zip(&area) {
        let l = norm(*n);
        if l == 0.0 {
            continue;
        }
        for k in 0..3 {
            let [a, b, c] = [t[k], t[(k + 1) % 3], t[(k + 2) % 3]].map(|i| m.verts[i as usize]);
            let (u, v) = (sub(b, a), sub(c, a));
            let cos = dot(u, v) / (dot(u, u) * dot(v, v)).sqrt();
            let w = libm::acos(cos.clamp(-1.0, 1.0)) / l;
            if w.is_finite() {
                vn[t[k] as usize] = std::array::from_fn(|j| vn[t[k] as usize][j] + w * n[j]);
            }
        }
    }
    let mix = |a: V3, b: V3, t: f64| -> V3 { std::array::from_fn(|c| a[c] + t * (b[c] - a[c])) };
    let len = |a: V3, b: V3| norm(sub(a, b));
    let parts = |l: f64| ((l / spacing).ceil() as usize).clamp(1, 1 << 16);
    // Directions: a vertex's normal is its triangles' areas summed; along an
    // edge, its ends' mixed; inside a triangle, its normal (its corners'
    // mixed, where it is a sliver). An edge is sampled by the triangle
    // holding it in increasing vertex order: each edge of a closed mesh once.
    let sample = |k: usize, out: &mut Vec<(V3, V3)>| {
        if k < m.verts.len() {
            out.push((m.verts[k], vn[k]));
            return;
        }
        let t = m.tris[k - m.verts.len()];
        let p = t.map(|i| m.verts[i as usize]);
        for e in 0..3 {
            let (i, j) = (t[e] as usize, t[(e + 1) % 3] as usize);
            if i < j {
                let n = 2 * parts(len(m.verts[i], m.verts[j]));
                for s in 1..n {
                    let f = s as f64 / n as f64;
                    out.push((mix(m.verts[i], m.verts[j], f), mix(vn[i], vn[j], f)));
                }
            }
        }
        // Rows across the triangle from its longest edge ab to c, each cut
        // into parts, sampled at the parts' middles.
        let e = (0..3).max_by(|&x, &y| len(p[x], p[(x + 1) % 3]).total_cmp(&len(p[y], p[(y + 1) % 3]))).unwrap();
        let (a, b, c) = (p[e], p[(e + 1) % 3], p[(e + 2) % 3]);
        let n = area[k - m.verts.len()];
        let corners = [vn[t[e] as usize], vn[t[(e + 1) % 3] as usize], vn[t[(e + 2) % 3] as usize]];
        let flat = has_area(p, n);
        let foot = {
            let ab = len(a, b);
            let u = ((c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1]) * (b[1] - a[1]) + (c[2] - a[2]) * (b[2] - a[2])) / (ab * ab).max(f64::MIN_POSITIVE);
            mix(a, b, u)
        };
        let rows = parts(len(foot, c));
        for r in 0..rows {
            let h = (r as f64 + 0.5) / rows as f64;
            let (l, rr) = (mix(a, c, h), mix(b, c, h));
            let cols = parts(len(l, rr));
            for q in 0..cols {
                let w = (q as f64 + 0.5) / cols as f64;
                let x = mix(l, rr, w);
                let d = if flat { n } else { mix(mix(corners[0], corners[2], h), mix(corners[1], corners[2], h), w) };
                out.push((x, d));
            }
        }
    };
    let mut d = Deviation { max: 0.0, at: [0.0; 3], beyond: 0, first_beyond: None };
    let total = m.verts.len() + m.tris.len();
    // Chunks of 2^18 vertices or triangles, sampled in parallel.
    for lo in (0..total).step_by(1 << 18) {
        let chunk: Vec<(V3, V3)> = (lo..(lo + (1 << 18)).min(total))
            .into_par_iter()
            .flat_map_iter(|k| {
                let mut v = Vec::new();
                sample(k, &mut v);
                v.into_iter().filter_map(|(p, n)| {
                    let l = norm(n);
                    (l > 0.0).then(|| (p, n.map(|c| c / l)))
                })
            })
            .collect();
        let c = deviation_of(field, grads, &chunk, step, reach);
        d.beyond += c.beyond;
        d.first_beyond = d.first_beyond.or(c.first_beyond);
        if c.max > d.max {
            (d.max, d.at) = (c.max, c.at);
        }
    }
    d
}

/// `deviation` over samples (point, unit normal).
fn deviation_of(field: Field, grads: &(dyn Fn(&[V3]) -> Vec<[f64; 4]> + Sync), samples: &[(V3, V3)], step: f64, reach: f64) -> Deviation {
    let f0 = field(&samples.par_iter().map(|s| s.0).collect::<Vec<_>>());
    // Each sample's bracket along its line: a sign change within (lo, hi]
    // along d * dir.
    type Bracket = Option<(f64, f64, f64)>;
    let search = |dir: &[V3], which: &[usize], out: &mut [Bracket]| {
        let point = |i: usize, s: f64| -> V3 { std::array::from_fn(|c| samples[i].0[c] + s * dir[i][c]) };
        let mut open: Vec<usize> = which.to_vec();
        let (mut lo, mut s) = (0.0, step);
        while !open.is_empty() && lo < reach {
            let f = field(&open.par_iter().flat_map_iter(|&i| [point(i, s), point(i, -s)]).collect::<Vec<_>>());
            let found: Vec<Option<f64>> =
                open.par_iter().enumerate().map(|(k, &i)| [(1.0, f[2 * k]), (-1.0, f[2 * k + 1])].into_iter().find(|&(_, v)| inside(v) != inside(f0[i])).map(|x| x.0)).collect();
            for (&i, d) in open.iter().zip(&found) {
                if let Some(d) = d {
                    out[i] = Some((lo, s, *d));
                }
            }
            open = open.iter().zip(&found).filter(|x| x.1.is_none()).map(|x| *x.0).collect();
            (lo, s) = (s, 2.0 * s);
        }
    };
    let normals: Vec<V3> = samples.iter().map(|s| s.1).collect();
    let mut along = vec![None; samples.len()];
    let all: Vec<usize> = (0..samples.len()).filter(|&i| f0[i] != 0.0).collect();
    search(&normals, &all, &mut along);
    // Where the normal finds no sign change within a few steps, the
    // gradient too: the normal can run along a face the sample lies just
    // off (where small and large triangles meet).
    let far: Vec<usize> = all.iter().copied().filter(|&i| along[i].is_none_or(|b| b.1 > 4.0 * step)).collect();
    let g = grads(&far.iter().map(|&i| samples[i].0).collect::<Vec<_>>());
    let mut dir = normals;
    for (&i, g) in far.iter().zip(&g) {
        dir[i] = direction([g[1], g[2], g[3]]).unwrap_or([0.0; 3]);
    }
    let mut across = vec![None; samples.len()];
    search(&dir, &far, &mut across);
    // Keep the nearer bracket, and the line it lies on.
    let mut bracket: Vec<Bracket> = vec![None; samples.len()];
    for i in 0..samples.len() {
        bracket[i] = match (along[i], across[i]) {
            (Some(a), Some(b)) if b.1 < a.1 => Some(b),
            (Some(a), _) => Some(a),
            (None, b) => b,
        };
        if across[i].is_none() || bracket[i] != across[i] {
            dir[i] = samples[i].1;
        }
        if f0[i] == 0.0 {
            bracket[i] = Some((0.0, 0.0, 1.0));
        }
    }
    let point = |i: usize, s: f64| -> V3 { std::array::from_fn(|c| samples[i].0[c] + s * dir[i][c]) };
    let top = bracket.iter().flatten().fold(0.0f64, |a, b| a.max(b.1));
    let worst: Vec<usize> = (0..samples.len()).filter(|&i| bracket[i].is_some_and(|b| b.1 == top && b.0 > 0.0)).collect();
    for _ in 0..8 {
        if worst.is_empty() {
            break;
        }
        let f = field(
            &worst
                .iter()
                .map(|&i| {
                    let (lo, hi, d) = bracket[i].unwrap();
                    point(i, d * 0.5 * (lo + hi))
                })
                .collect::<Vec<_>>(),
        );
        for (&i, v) in worst.iter().zip(f) {
            let (lo, hi, d) = bracket[i].unwrap();
            let mid = 0.5 * (lo + hi);
            bracket[i] = Some(if inside(v) != inside(f0[i]) { (lo, mid, d) } else { (mid, hi, d) });
        }
    }
    let open: Vec<usize> = (0..samples.len()).filter(|&i| bracket[i].is_none()).collect();
    let mut d = Deviation { max: 0.0, at: [0.0; 3], beyond: open.len(), first_beyond: open.first().map(|&i| samples[i].0) };
    for (i, b) in bracket.iter().enumerate() {
        if let Some((_, hi, _)) = b
            && *hi > d.max
        {
            (d.max, d.at) = (*hi, samples[i].0);
        }
    }
    d
}
