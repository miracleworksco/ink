//! Meshes to fields: STL reading and the signed distance to a triangle
//! mesh, in floats.

use ink_core::V3;
use ink_core::v3::{EMPTY_BOX, box_union, dot, norm, normal, sub};

/// A map keyed by small integers, hashed fast.
type Map<K, V> = std::collections::HashMap<K, V, rustc_hash::FxBuildHasher>;

/// Read a binary or ASCII STL into triangles.
pub fn read_stl(bytes: &[u8]) -> Result<Vec<[V3; 3]>, String> {
    // Binary headers can begin with "solid". An exact binary length wins;
    // otherwise recognize ASCII even with whitespace or no complete facet.
    let binary = bytes.len() >= 84 && 84 + u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as u64 * 50 == bytes.len() as u64;
    let ascii = !binary && std::str::from_utf8(&bytes[..bytes.len().min(512)]).is_ok_and(|s| s.split_whitespace().next().is_some_and(|w| w.eq_ignore_ascii_case("solid")));
    if ascii {
        let text = std::str::from_utf8(bytes).map_err(|_| "STL is not valid UTF-8")?;
        return read_ascii(text);
    }
    if bytes.len() < 84 {
        return Err("STL too short".into());
    }
    let n = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    if bytes.len() < 84 + 50 * n {
        return Err(format!("binary STL says {n} triangles but is truncated"));
    }
    Ok((0..n)
        .map(|t| {
            let o = 84 + 50 * t + 12;
            let f = |k: usize| f32::from_le_bytes(bytes[o + 4 * k..o + 4 * k + 4].try_into().unwrap()) as f64;
            [[f(0), f(1), f(2)], [f(3), f(4), f(5)], [f(6), f(7), f(8)]]
        })
        .collect())
}

/// An ASCII STL, strictly: solids of facets, each `facet ...`, `outer
/// loop`, three `vertex x y z`, `endloop`, `endfacet`, and `endsolid`. What
/// is not that is an error at its line, never skipped: a dropped vertex
/// would shift every triangle after it.
fn read_ascii(text: &str) -> Result<Vec<[V3; 3]>, String> {
    // What comes next: 0 solid, 1 facet or endsolid, 2 outer loop, 3..=5
    // a vertex, 6 endloop, 7 endfacet.
    let mut at = 0;
    let mut tri = [[0.0; 3]; 3];
    let mut tris = Vec::new();
    let mut last = 0;
    for (i, line) in text.lines().enumerate() {
        let words: Vec<&str> = line.split_whitespace().collect();
        let Some(word) = words.first() else { continue };
        last = i + 1;
        let err = |want: &str| format!("STL line {}: expected {want}, got '{}'", i + 1, line.trim());
        at = match (at, word.to_ascii_lowercase().as_str()) {
            (0, "solid") => 1,
            (1, "facet") => {
                if words.len() != 5 || !words[1].eq_ignore_ascii_case("normal") || words[2..].iter().any(|n| n.parse::<f64>().is_err()) {
                    return Err(err("'facet normal nx ny nz', three numbers"));
                }
                2
            }
            (1, "endsolid") => 0,
            (2, "outer") if words.get(1).is_some_and(|w| w.eq_ignore_ascii_case("loop")) && words.len() == 2 => 3,
            (3..=5, "vertex") => {
                let v: Vec<f64> = words[1..].iter().filter_map(|x| x.parse().ok()).collect();
                if words.len() != 4 || v.len() != 3 {
                    return Err(err("'vertex x y z', three numbers"));
                }
                tri[at - 3] = [v[0], v[1], v[2]];
                at + 1
            }
            (6, "endloop") if words.len() == 1 => 7,
            (7, "endfacet") if words.len() == 1 => {
                tris.push(tri);
                1
            }
            (0, _) => return Err(err("'solid'")),
            (1, _) => return Err(err("'facet' or 'endsolid'")),
            (2, _) => return Err(err("'outer loop'")),
            (3..=5, _) => return Err(err(&format!("vertex {} of 3", at - 2))),
            (6, _) => return Err(err("'endloop' (a facet has three vertices)")),
            _ => return Err(err("'endfacet'")),
        };
    }
    if at != 0 {
        return Err(format!("STL ends inside a {} at line {last}: the file is truncated", if at == 1 { "solid (no endsolid)" } else { "facet" }));
    }
    Ok(tris)
}

// The triangles `closest_on_triangle` measured on this thread: the tests'
// measure of a search's work.
#[cfg(test)]
thread_local!(static MEASURED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) });

/// A lower bound on the distance from a point to a triangle, 1-Lipschitz
/// in the point (so at least its value at c less the distance from c):
/// the distance to its box, and to the capsule about its longest edge
/// that holds it (radius the third corner's distance from that edge),
/// each found by well-conditioned float steps, less a generous allowance
/// for their rounding. Unlike `closest_on_triangle`'s distance, it never
/// exceeds the true one on a thin facet.
#[derive(Clone, Copy)]
struct Bound {
    lo: V3,
    hi: V3,
    a: V3,
    b: V3,
    radius: f64,
    slack: f64,
}

impl Bound {
    fn of(v: [V3; 3], lo: V3, hi: V3) -> Bound {
        let len = |k: usize| norm(sub(v[(k + 1) % 3], v[k]));
        let k = (0..3).max_by(|&x, &y| len(x).total_cmp(&len(y))).unwrap();
        let (a, b, c) = (v[k], v[(k + 1) % 3], v[(k + 2) % 3]);
        let slack = 1e-12 * (1.0 + v.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs())));
        Bound { lo, hi, a, b, radius: segment_distance(c, a, b) + slack, slack }
    }

    fn at(&self, p: V3) -> f64 {
        let boxed = (0..3).map(|k| (self.lo[k] - p[k]).max(0.0).max(p[k] - self.hi[k]).powi(2)).sum::<f64>().sqrt();
        let far = 1e-12 * p.iter().fold(0.0f64, |m, x| m.max(x.abs()));
        boxed.max(segment_distance(p, self.a, self.b) - self.radius) - self.slack - far
    }

    /// Whether `at(p)` is below d, by squares (`far`: as `at` finds it).
    fn below(&self, p: V3, far: f64, d: f64) -> bool {
        let d = d + self.slack + far;
        (0..3).map(|k| (self.lo[k] - p[k]).max(0.0).max(p[k] - self.hi[k]).powi(2)).sum::<f64>() < d * d && segment_d2(p, self.a, self.b) < (d + self.radius).powi(2)
    }
}

/// The distance from p to the segment ab, in floats.
fn segment_distance(p: V3, a: V3, b: V3) -> f64 {
    segment_d2(p, a, b).sqrt()
}

/// Its square.
fn segment_d2(p: V3, a: V3, b: V3) -> f64 {
    let ab = sub(b, a);
    let l = dot(ab, ab);
    let t = if l > 0.0 { (dot(sub(p, a), ab) / l).clamp(0.0, 1.0) } else { 0.0 };
    let e = sub(p, [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]]);
    dot(e, e)
}

/// Closest point on triangle (Ericson), and which feature: 0 face, 1..=3
/// vertex a/b/c, 4..=6 edge ab/bc/ca.
fn closest_on_triangle(p: V3, a: V3, b: V3, c: V3) -> (V3, u8) {
    #[cfg(test)]
    MEASURED.with(|m| m.set(m.get() + 1));
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, 1);
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return (b, 2);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return ([a[0] + v * ab[0], a[1] + v * ab[1], a[2] + v * ab[2]], 4);
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return (c, 3);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return ([a[0] + w * ac[0], a[1] + w * ac[1], a[2] + w * ac[2]], 6);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return ([b[0] + w * (c[0] - b[0]), b[1] + w * (c[1] - b[1]), b[2] + w * (c[2] - b[2])], 5);
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    ([a[0] + ab[0] * v + ac[0] * w, a[1] + ab[1] * v + ac[1] * w, a[2] + ab[2] * v + ac[2] * w], 0)
}

/// A node of MeshDistance's tree: its box (lo, hi), then its children's
/// node ids, or for a leaf its first triangle in the order and how many.
type Node = (V3, V3, u32, u32, bool);
const EMPTY: V3 = [0.0; 3];

/// A box turned to axes (orthonormal, in floats): what it holds lies
/// between lo and hi along each. Across a slanted sliver it is as thin as
/// the sliver, where its box is as wide as it is long.
struct Oriented {
    axes: [V3; 3],
    lo: V3,
    hi: V3,
}

/// Distance to a triangle mesh, in floats: signed by angle-weighted
/// pseudo-normals (`distance`), and for a closed mesh by parity, the side
/// of the solid it encloses (`signed`, `sample`).
pub struct MeshDistance {
    verts: Vec<V3>,
    tris: Vec<[u32; 3]>,
    face_n: Vec<V3>,
    vert_n: Vec<V3>,
    /// Each edge's pseudo-normal (its triangles' normals summed) and how
    /// many triangles it is in.
    edge_n: Map<(u32, u32), (V3, u32)>,
    // AABB tree (`build`)
    nodes: Vec<Node>,
    order: Vec<u32>,
    /// Each node's box turned to its largest triangle (`oriented`).
    turned: Vec<Oriented>,
    /// Every edge in an even number of triangles (`odd_edges`).
    closed: bool,
    /// The edges in an odd number of triangles (`odd_edges`).
    odd: usize,
}

impl MeshDistance {
    pub fn new(tris_in: &[[V3; 3]]) -> MeshDistance {
        Self::made(tris_in, false).expect("not timed")
    }

    /// `new`, or None where the call's time limit (ink_kernel::stop) passes
    /// first: each step checks it.
    pub fn in_time(tris_in: &[[V3; 3]]) -> Option<MeshDistance> {
        Self::made(tris_in, true)
    }

    fn made(tris_in: &[[V3; 3]], timed: bool) -> Option<MeshDistance> {
        let stopped = || timed && ink_kernel::stop::reached();
        // Weld identical positions into shared vertices.
        let mut index: Map<[u64; 3], u32> = Map::default();
        let mut verts = Vec::new();
        let mut tris = Vec::with_capacity(tris_in.len());
        for (i, t) in tris_in.iter().enumerate() {
            if i % (1 << 16) == 0 && stopped() {
                return None;
            }
            let ids = t.map(|v| {
                *index.entry(v.map(|c| (c + 0.0).to_bits())).or_insert_with(|| {
                    verts.push(v);
                    (verts.len() - 1) as u32
                })
            });
            if ids[0] != ids[1] && ids[1] != ids[2] && ids[0] != ids[2] {
                tris.push(ids);
            }
        }
        let face_n: Vec<V3> = tris
            .iter()
            .map(|t| {
                let n = normal(t.map(|i| verts[i as usize]));
                let l = norm(n).max(1e-300);
                n.map(|c| c / l)
            })
            .collect();
        let mut vert_n = vec![[0.0; 3]; verts.len()];
        let mut edge_n: Map<(u32, u32), (V3, u32)> = Map::default();
        for (i, (t, n)) in tris.iter().zip(&face_n).enumerate() {
            if i % (1 << 16) == 0 && stopped() {
                return None;
            }
            for k in 0..3 {
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let (e1, e2) = (sub(verts[b as usize], verts[a as usize]), sub(verts[c as usize], verts[a as usize]));
                let ang = libm::acos((dot(e1, e2) / (norm(e1) * norm(e2)).max(1e-300)).clamp(-1.0, 1.0));
                for d in 0..3 {
                    vert_n[a as usize][d] += ang * n[d];
                }
                let e = edge_n.entry((a.min(b), a.max(b))).or_insert(([0.0; 3], 0));
                for d in 0..3 {
                    e.0[d] += n[d];
                }
                e.1 += 1;
            }
        }
        let odd = edge_n.values().filter(|e| e.1 % 2 == 1).count();
        let mut md = MeshDistance { verts, tris, face_n, vert_n, edge_n, nodes: Vec::new(), order: Vec::new(), turned: Vec::new(), closed: odd == 0, odd };
        let mut order: Vec<u32> = (0..md.tris.len() as u32).collect();
        let boxes: Vec<(V3, V3)> = order.iter().map(|&t| md.tri_box(t)).collect();
        let mut sizes = std::collections::HashMap::new();
        let mut nodes = vec![(EMPTY, EMPTY, 0, 0, true); Self::tree_size(order.len(), &mut sizes)];
        Self::build(&boxes, &mut order, 0, &mut nodes, 0, &sizes, &stopped);
        if stopped() {
            return None;
        }
        (md.nodes, md.order) = (nodes, order);
        md.turned = md.oriented();
        if stopped() {
            return None;
        }
        Some(md)
    }

    /// How many nodes `build` makes over n triangles (each count it splits
    /// into kept in `sizes`).
    fn tree_size(n: usize, sizes: &mut std::collections::HashMap<usize, usize>) -> usize {
        if n <= 4 {
            return 1;
        }
        if let Some(&k) = sizes.get(&n) {
            return k;
        }
        let k = 1 + Self::tree_size(n / 2, sizes) + Self::tree_size(n - n / 2, sizes);
        sizes.insert(n, k);
        k
    }

    fn tri_box(&self, t: u32) -> (V3, V3) {
        let v = self.tris[t as usize].map(|i| self.verts[i as usize]);
        (std::array::from_fn(|a| v[0][a].min(v[1][a]).min(v[2][a])), std::array::from_fn(|a| v[0][a].max(v[1][a]).max(v[2][a])))
    }

    /// The tree over the triangles `order` holds (from `start` in the whole
    /// order), in `nodes` (its own, from node `id`) in preorder: each node
    /// the box of its triangles, split at the middle of them sorted (stably)
    /// along its box's longest axis by their boxes' centres. Its halves are
    /// built side by side, each into the nodes its size takes. A large part
    /// is not begun once `stopped` (the tree is then not made).
    fn build(
        boxes: &[(V3, V3)],
        order: &mut [u32],
        start: usize,
        nodes: &mut [Node],
        id: u32,
        sizes: &std::collections::HashMap<usize, usize>,
        stopped: &(dyn Fn() -> bool + Sync),
    ) {
        use rayon::prelude::*;
        let n = order.len();
        if n > 1 << 12 && stopped() {
            return;
        }
        let (lo, hi) = order.iter().map(|&t| boxes[t as usize]).fold(EMPTY_BOX, box_union);
        if n <= 4 {
            nodes[0] = (lo, hi, start as u32, n as u32, true);
            return;
        }
        let axis = (0..3).max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b]))).unwrap();
        let key = |t: &u32| boxes[*t as usize].0[axis] + boxes[*t as usize].1[axis];
        let by = |x: &u32, y: &u32| key(x).total_cmp(&key(y));
        if n > 1 << 14 {
            order.par_sort_by(by)
        } else {
            order.sort_by(by)
        }
        let mid = n / 2;
        let left = if mid <= 4 { 1 } else { sizes[&mid] };
        nodes[0] = (lo, hi, id + 1, id + 1 + left as u32, false);
        let (lo_o, hi_o) = order.split_at_mut(mid);
        let (ln, rn) = nodes[1..].split_at_mut(left);
        let rid = id + 1 + left as u32;
        if n > 1 << 12 {
            rayon::join(|| Self::build(boxes, lo_o, start, ln, id + 1, sizes, stopped), || Self::build(boxes, hi_o, start + mid, rn, rid, sizes, stopped));
        } else {
            Self::build(boxes, lo_o, start, ln, id + 1, sizes, stopped);
            Self::build(boxes, hi_o, start + mid, rn, rid, sizes, stopped);
        }
    }

    /// Each node's box turned to its largest triangle: along its longest
    /// edge, its normal and their cross, each a little larger for the
    /// rounding of its corners' coordinates there.
    fn oriented(&self) -> Vec<Oriented> {
        use rayon::prelude::*;
        // Each node's run of `order`.
        let mut runs = vec![(0u32, 0u32); self.nodes.len()];
        fn run(nodes: &[Node], n: usize, runs: &mut [(u32, u32)]) -> (u32, u32) {
            let (_, _, a, b, leaf) = nodes[n];
            runs[n] = if leaf {
                (a, a + b)
            } else {
                let (l, r) = (run(nodes, a as usize, runs), run(nodes, b as usize, runs));
                (l.0, r.1)
            };
            runs[n]
        }
        if !self.nodes.is_empty() {
            run(&self.nodes, 0, &mut runs);
        }
        let slack = 1e-12 * (1.0 + self.verts.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs())));
        let corners = |t: u32| self.tris[t as usize].map(|i| self.verts[i as usize]);
        let area: Vec<f64> = (0..self.tris.len() as u32).into_par_iter().map(|t| norm(normal(corners(t)))).collect();
        runs.par_iter()
            .map(|&(s, e)| {
                let ts = &self.order[s as usize..e as usize];
                let Some(&big) = ts.iter().max_by(|&&x, &&y| area[x as usize].total_cmp(&area[y as usize])) else {
                    return Oriented { axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], lo: [f64::INFINITY; 3], hi: [f64::NEG_INFINITY; 3] };
                };
                let v = corners(big);
                let len = |k: usize| norm(sub(v[(k + 1) % 3], v[k]));
                let k = (0..3).max_by(|&x, &y| len(x).total_cmp(&len(y))).unwrap();
                let unit = |x: V3| -> Option<V3> {
                    let l = norm(x);
                    (l > 0.0 && l.is_finite()).then(|| x.map(|c| c / l))
                };
                let axes = unit(sub(v[(k + 1) % 3], v[k]))
                    .and_then(|u| {
                        let n = normal(v);
                        let along = dot(n, u);
                        unit(std::array::from_fn(|a| n[a] - along * u[a])).map(|n| (u, n))
                    })
                    .map(|(u, n)| [u, n, [n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]]])
                    .unwrap_or([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
                let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
                for &t in ts {
                    for c in corners(t) {
                        for a in 0..3 {
                            let x = dot(axes[a], c);
                            lo[a] = lo[a].min(x);
                            hi[a] = hi[a].max(x);
                        }
                    }
                }
                Oriented { axes, lo: lo.map(|x| x - slack), hi: hi.map(|x| x + slack) }
            })
            .collect()
    }

    /// The squared distance from p to node n's boxes (the larger), no more
    /// than to what it holds (`slack`: for the rounding of p's coordinates
    /// along the turned axes, generously).
    fn node_d2(&self, n: u32, p: V3, slack: f64) -> f64 {
        let (lo, hi, ..) = self.nodes[n as usize];
        let o = &self.turned[n as usize];
        let (mut aa, mut oo) = (0.0, 0.0);
        for a in 0..3 {
            aa += (lo[a] - p[a]).max(0.0).max(p[a] - hi[a]).powi(2);
            let x = dot(o.axes[a], p);
            oo += ((o.lo[a] - x).max(x - o.hi[a]) - slack).max(0.0).powi(2);
        }
        aa.max(oo)
    }

    /// The triangle nearest p, in floats: (squared distance, nearest point,
    /// triangle, feature as `closest_on_triangle`).
    fn closest(&self, p: V3) -> (f64, V3, u32, u8) {
        self.closest_within(p, f64::INFINITY)
    }

    /// `closest`, given that a triangle lies within `reach` of p (only
    /// nearer ones are looked at).
    fn closest_within(&self, p: V3, reach: f64) -> (f64, V3, u32, u8) {
        self.closest_from(p, (reach * reach * (1.0 + 1e-9), [0.0; 3], 0u32, 0u8))
    }

    /// `closest`, from `best`: what a triangle gives (or a squared reach a
    /// triangle lies within). Only nearer ones are looked at.
    fn closest_from(&self, p: V3, mut best: (f64, V3, u32, u8)) -> (f64, V3, u32, u8) {
        let mut stack = Vec::with_capacity(64);
        stack.push(0u32);
        let slack = 1e-12 * (1.0 + p.iter().fold(0.0f64, |m, x| m.max(x.abs())));
        let box_d2 = |n: u32| self.node_d2(n, p, slack);
        while let Some(n) = stack.pop() {
            let (_, _, a, b, leaf) = self.nodes[n as usize];
            if box_d2(n) >= best.0 {
                continue;
            }
            if leaf {
                for &t in &self.order[a as usize..(a + b) as usize] {
                    let v = self.tris[t as usize].map(|i| self.verts[i as usize]);
                    let (q, feat) = closest_on_triangle(p, v[0], v[1], v[2]);
                    let d2 = dot(sub(p, q), sub(p, q));
                    if d2 < best.0 {
                        best = (d2, q, t, feat);
                    }
                }
            } else {
                let (dl, dr) = (box_d2(a), box_d2(b));
                if dl < dr {
                    stack.push(b);
                    stack.push(a);
                } else {
                    stack.push(a);
                    stack.push(b);
                }
            }
        }
        best
    }

    /// The triangles that may lie within `reach` of a point of box lo..hi
    /// (centre c, half-diagonal at most r): their boxes within it of the
    /// box, and their `Bound` within reach + r of c. Each with its bound
    /// and that bound at c, nearest c first; None where more than `most`
    /// are.
    fn near(&self, bounds: &[Bound], lo: V3, hi: V3, c: V3, r: f64, reach: f64, most: usize) -> Option<Vec<(f64, u32, Bound)>> {
        let apart = |a: V3, b: V3| -> f64 { (0..3).map(|k| (a[k] - hi[k]).max(lo[k] - b[k]).max(0.0).powi(2)).sum() };
        let r2 = reach * reach;
        let mut out = Vec::new();
        let mut stack = vec![0u32];
        let slack = 1e-12 * (1.0 + c.iter().fold(0.0f64, |m, x| m.max(x.abs())));
        while let Some(n) = stack.pop() {
            let (nlo, nhi, a, b, leaf) = self.nodes[n as usize];
            if apart(nlo, nhi) > r2 || self.node_d2(n, c, slack) > (reach + r) * (reach + r) {
                continue;
            }
            if !leaf {
                stack.extend([a, b]);
                continue;
            }
            for &t in &self.order[a as usize..(a + b) as usize] {
                let bound = bounds[t as usize];
                if apart(bound.lo, bound.hi) > r2 {
                    continue;
                }
                let at = bound.at(c);
                if at <= reach + r {
                    if out.len() == most {
                        return None;
                    }
                    out.push((at, t, bound));
                }
            }
        }
        out.sort_by(|x, y| x.0.total_cmp(&y.0));
        Some(out)
    }

    /// Signed distance at p (negative inside), in floats, signed by the
    /// nearest feature's pseudo-normal: exact in real arithmetic for a
    /// closed, consistently oriented mesh facing out that does not cross
    /// itself.
    pub fn distance(&self, p: V3) -> f64 {
        self.distance_within(p, f64::INFINITY)
    }

    /// `distance`, given that the mesh lies within `reach` of p (a
    /// neighbour's distance plus how far it is: only nearer triangles are
    /// looked at, and the answer is the same).
    pub fn distance_within(&self, p: V3, reach: f64) -> f64 {
        let (d2, q, t, feat) = self.closest_within(p, reach);
        let tri = self.tris[t as usize];
        let n = match feat {
            0 => self.face_n[t as usize],
            1..=3 => self.vert_n[tri[(feat - 1) as usize] as usize],
            _ => {
                let (a, b) = match feat {
                    4 => (tri[0], tri[1]),
                    5 => (tri[1], tri[2]),
                    _ => (tri[2], tri[0]),
                };
                self.edge_n[&(a.min(b), a.max(b))].0
            }
        };
        let s = if dot(sub(p, q), n) < 0.0 { -1.0 } else { 1.0 };
        s * d2.sqrt()
    }

    /// How many triangles it has: those of nonzero area (`new` drops the
    /// rest).
    pub fn triangles(&self) -> usize {
        self.tris.len()
    }

    pub fn bounds(&self) -> (V3, V3) {
        (self.nodes[0].0, self.nodes[0].1)
    }

    /// Points on the mesh (its distance is 0 there): every vertex, where
    /// sharp features are, and every triangle's centroid.
    pub fn surface_points(&self) -> Vec<V3> {
        let centroid = |t: &[u32; 3]| -> V3 { std::array::from_fn(|a| t.iter().map(|&i| self.verts[i as usize][a]).sum::<f64>() / 3.0) };
        self.verts.iter().copied().chain(self.tris.iter().map(centroid)).collect()
    }

    /// The number of edges (of the welded mesh) with an odd number of
    /// triangles. With none, the mesh bounds a solid (a mod-2 cycle: a path
    /// crosses it an even number of times between any two points off it, so
    /// "inside" is well defined: an odd number of crossings from outside),
    /// whichever way its triangles face and even where it crosses itself.
    pub fn odd_edges(&self) -> usize {
        self.odd
    }

    /// Where the line through (y, z) along x crosses the mesh, in order:
    /// None where that is unclear (the line meets an edge or a vertex of
    /// it, in its yz projection, or crosses it an odd number of times,
    /// which no line through a closed mesh does).
    fn crossings(&self, y: f64, z: f64) -> Option<Vec<f64>> {
        let mut xs = Vec::new();
        let mut stack = vec![0u32];
        while let Some(n) = stack.pop() {
            let (lo, hi, a, b, leaf) = self.nodes[n as usize];
            if y < lo[1] || y > hi[1] || z < lo[2] || z > hi[2] {
                continue;
            }
            if !leaf {
                stack.extend([a, b]);
                continue;
            }
            for &t in &self.order[a as usize..(a + b) as usize] {
                let v = self.tris[t as usize].map(|i| self.verts[i as usize]);
                // Each edge's side of (y, z): the weight of the corner
                // across from it.
                let side = |p: V3, q: V3| (q[1] - p[1]) * (z - p[2]) - (q[2] - p[2]) * (y - p[1]);
                let w = [side(v[1], v[2]), side(v[2], v[0]), side(v[0], v[1])];
                if w.iter().all(|&x| x > 0.0) || w.iter().all(|&x| x < 0.0) {
                    let total = w[0] + w[1] + w[2];
                    xs.push((w[0] * v[0][0] + w[1] * v[1][0] + w[2] * v[2][0]) / total);
                } else if w.contains(&0.0) && !(w.iter().any(|&x| x > 0.0) && w.iter().any(|&x| x < 0.0)) && w.iter().any(|&x| x != 0.0) {
                    return None;
                }
            }
        }
        xs.sort_by(f64::total_cmp);
        (xs.len() % 2 == 0).then_some(xs)
    }

    /// The distance at p, negative inside: for a closed mesh by parity
    /// (inside where the line from p along +x crosses it an odd number of
    /// times: `xs`, the crossings of the line through p), else (or where
    /// those are unclear) as `distance` signs it.
    /// `reach`: a triangle lies within it of p.
    fn signed_on(&self, xs: Option<&[f64]>, p: V3, reach: f64) -> f64 {
        match xs.filter(|_| self.closed) {
            Some(xs) => {
                let inside = xs.iter().filter(|&&c| c > p[0]).count() % 2 == 1;
                (if inside { -1.0 } else { 1.0 }) * self.closest_within(p, reach).0.sqrt()
            }
            None => self.distance(p),
        }
    }

    /// `signed_on` at p.
    pub fn signed(&self, p: V3) -> f64 {
        self.signed_on(self.crossings(p[1], p[2]).as_deref(), p, f64::INFINITY)
    }

    /// `signed` at the points origin + h (i, j, k) of a grid (i fastest),
    /// and how many lines of it (along x) were signed by `distance` for
    /// being unclear (none for a mesh that is not closed: it is signed so
    /// throughout). A block of BLOCK^3 points gathers once the triangles
    /// that may lie within d + r of a point of it (d its centre's
    /// distance, r its half-diagonal: no point is farther from the centre's
    /// nearest), by a lower bound on their distance (`Bound`). Each point,
    /// from the triangle nearest the point before it, reads those nearest
    /// the centre first by that bound, passing over one bounded no nearer
    /// than its nearest so far and stopping where the rest are, and keeps
    /// what it finds where that is within the reach gathered: every
    /// triangle not gathered is bounded farther from it. Else, or where
    /// the block is near more than MOST triangles (far from the surface,
    /// where a ball reaches much of it), it asks the tree itself, from that
    /// triangle too.
    pub fn sample(&self, origin: V3, h: f64, dims: [usize; 3]) -> (Vec<f64>, usize) {
        use rayon::prelude::*;
        const BLOCK: usize = 4;
        // Triangles a block reads at most, past which its points ask the tree.
        const MOST: usize = 2048;
        let lines: Vec<Option<Vec<f64>>> =
            (0..dims[1] * dims[2]).into_par_iter().map(|r| self.crossings(origin[1] + h * (r % dims[1]) as f64, origin[2] + h * (r / dims[1]) as f64)).collect();
        let unclear = lines.iter().filter(|xs| self.closed && xs.is_none()).count();
        if self.tris.is_empty() {
            return (vec![f64::INFINITY; dims[0] * dims[1] * dims[2]], unclear);
        }
        let bounds: Vec<Bound> = (0..self.tris.len() as u32)
            .into_par_iter()
            .map(|t| {
                let (lo, hi) = self.tri_box(t);
                Bound::of(self.tris[t as usize].map(|i| self.verts[i as usize]), lo, hi)
            })
            .collect();
        let at = |i: usize, j: usize, k: usize| -> V3 { [origin[0] + h * i as f64, origin[1] + h * j as f64, origin[2] + h * k as f64] };
        let layer = dims[0] * dims[1];
        let mut values = vec![0.0; layer * dims[2]];
        values.par_chunks_mut(layer * BLOCK).enumerate().for_each(|(kb, slab)| {
            let k0 = kb * BLOCK;
            let nk = slab.len() / layer;
            for j0 in (0..dims[1]).step_by(BLOCK) {
                let nj = BLOCK.min(dims[1] - j0);
                for i0 in (0..dims[0]).step_by(BLOCK) {
                    let ni = BLOCK.min(dims[0] - i0);
                    let half = |n: usize| 0.5 * h * (n - 1) as f64;
                    let c: V3 = [origin[0] + h * i0 as f64 + half(ni), origin[1] + h * j0 as f64 + half(nj), origin[2] + h * k0 as f64 + half(nk)];
                    let r = (half(ni).powi(2) + half(nj).powi(2) + half(nk).powi(2)).sqrt();
                    let (blo, bhi) = (at(i0, j0, k0), at(i0 + ni - 1, j0 + nj - 1, k0 + nk - 1));
                    // The half-diagonal, rounded up generously.
                    let r = r * (1.0 + 1e-9);
                    let (dc2, _, mut seed, _) = self.closest(c);
                    let reach = (dc2.sqrt() + r) * (1.0 + 1e-9);
                    let near = self.near(&bounds, blo, bhi, c, r, reach, MOST);
                    // Rounding in the distances to boxes from c and p's
                    // distance from c, generously.
                    let slack = 1e-12 * (1.0 + c.iter().fold(0.0f64, |m, x| m.max(x.abs())));

                    for k in 0..nk {
                        for j in j0..j0 + nj {
                            let xs = lines[j + dims[1] * (k0 + k)].as_deref().filter(|_| self.closed);
                            for i in i0..i0 + ni {
                                let p = at(i, j, k0 + k);
                                let pc = norm(sub(p, c));
                                let far = 1e-12 * p.iter().fold(0.0f64, |m, x| m.max(x.abs()));
                                let found = match (xs, &near) {
                                    (Some(xs), Some(near)) => {
                                        // From the triangle nearest the point before.
                                        let v = self.tris[seed as usize].map(|i| self.verts[i as usize]);
                                        let q = closest_on_triangle(p, v[0], v[1], v[2]).0;
                                        let mut best = dot(sub(p, q), sub(p, q));
                                        let mut root = best.sqrt();
                                        for (dc, t, bound) in near {
                                            // It, and those after it, at least this far from p.
                                            if dc - pc - slack > root {
                                                break;
                                            }
                                            if !bound.below(p, far, root) {
                                                continue;
                                            }
                                            let v = self.tris[*t as usize].map(|i| self.verts[i as usize]);
                                            let q = closest_on_triangle(p, v[0], v[1], v[2]).0;
                                            let d2 = dot(sub(p, q), sub(p, q));
                                            if d2 < best {
                                                (best, root, seed) = (d2, d2.sqrt(), *t);
                                            }
                                        }
                                        let inside = xs.iter().filter(|&&x| x > p[0]).count() % 2 == 1;
                                        (best <= reach * reach).then(|| (if inside { -1.0 } else { 1.0 }) * best.sqrt())
                                    }
                                    _ => None,
                                };
                                slab[i + dims[0] * j + layer * k] = match (found, xs) {
                                    (Some(v), _) => v,
                                    // The tree, from the triangle nearest the point before.
                                    (None, Some(xs)) => {
                                        let v = self.tris[seed as usize].map(|i| self.verts[i as usize]);
                                        let (q, feat) = closest_on_triangle(p, v[0], v[1], v[2]);
                                        let best = self.closest_from(p, (dot(sub(p, q), sub(p, q)), q, seed, feat));
                                        seed = best.2;
                                        let inside = xs.iter().filter(|&&x| x > p[0]).count() % 2 == 1;
                                        (if inside { -1.0 } else { 1.0 }) * best.0.sqrt()
                                    }
                                    (None, None) => self.distance(p),
                                };
                            }
                        }
                    }
                }
            }
        });
        (values, unclear)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A coin (radius 10, 1 thick) whose faces are fans of n slivers
    /// about its centre: each sliver's box overlaps hundreds of others'.
    fn coin_of(n: usize) -> Vec<[V3; 3]> {
        let rim = |k: usize, z: f64| -> V3 {
            let a = 2.0 * std::f64::consts::PI * (k % n) as f64 / n as f64;
            [10.0 * a.cos(), 10.0 * a.sin(), z]
        };
        let mut t = Vec::new();
        for k in 0..n {
            t.push([[0.0, 0.0, 1.0], rim(k, 1.0), rim(k + 1, 1.0)]);
            t.push([[0.0, 0.0, 0.0], rim(k + 1, 0.0), rim(k, 0.0)]);
            t.push([rim(k, 0.0), rim(k + 1, 0.0), rim(k + 1, 1.0)]);
            t.push([rim(k, 0.0), rim(k + 1, 1.0), rim(k, 1.0)]);
        }
        t
    }

    fn coin() -> Vec<[V3; 3]> {
        coin_of(720)
    }

    /// A block keeps a thin facet's nearness where `closest_on_triangle`
    /// overstates it: a tetrahedron of four long facets 1e-5 thin sampled
    /// at 1e-5 reads 7.1e-6 mm at a point, as `signed` there, not the 1.6e-3
    /// of a block that gathered triangles by that function's distances.
    #[test]
    fn a_blocks_candidates_keep_the_nearest_thin_triangle() {
        let a = [-2.8977450831680853, -3.149304231523823, 4.1636086179063545];
        let b = [-7.726046376475313, -4.309005578162985, -5.142638580228725];
        let c = [-7.378670127498295, -4.225570277481827, -4.4730936361585885];
        let d = [-7.378670229554812, -4.225570145124801, -4.473093193836915];
        let md = MeshDistance::new(&[[a, b, c], [a, d, b], [a, c, d], [b, d, c]]);
        assert_eq!(md.odd_edges(), 0);
        let (origin, h, dims) = ([-7.2046590605025855, -4.183785961964447, -4.137683642101432], 9.999999999999999e-6, [4, 4, 4]);
        let (grid, _) = md.sample(origin, h, dims);
        for (n, v) in grid.iter().enumerate() {
            let p = [origin[0] + h * (n % 4) as f64, origin[1] + h * (n / 4 % 4) as f64, origin[2] + h * (n / 16) as f64];
            assert_eq!(v.to_bits(), md.signed(p).to_bits(), "at {p:?}: {v}, signed {}", md.signed(p));
        }
    }

    /// A grid's samples are each point's `signed`, to the bit, and near a
    /// fan of slivers a point measures 16 triangles (each point's own
    /// search of the tree measured 260: the slivers' boxes overlap it).
    #[test]
    fn a_grid_near_slivers_measures_few_triangles_a_point() {
        let md = MeshDistance::new(&coin());
        assert!(md.closed);
        let (origin, h, dims) = ([2.0191, 2.0137, -0.3191], 0.05, [24, 24, 32]);
        let pool = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let (values, measured) = pool.install(|| {
            MEASURED.with(|m| m.set(0));
            let (v, unclear) = md.sample(origin, h, dims);
            assert_eq!(unclear, 0);
            (v, MEASURED.with(|m| m.get()))
        });
        for (n, v) in values.iter().enumerate() {
            let [i, j, k] = [n % dims[0], n / dims[0] % dims[1], n / (dims[0] * dims[1])];
            let p = [origin[0] + h * i as f64, origin[1] + h * j as f64, origin[2] + h * k as f64];
            assert_eq!(v.to_bits(), md.signed(p).to_bits(), "at {p:?}");
        }
        let points = values.len();
        assert!(measured < 60 * points, "{measured} triangles measured for {points} points");
    }

    /// Where a block is near more triangles than it gathers (the centre of
    /// a coin of 4800 triangles, 2400 face slivers meeting at its centres),
    /// its points ask the tree, and read as `signed` does, to the bit.
    #[test]
    fn a_block_near_too_many_triangles_reads_as_signed() {
        let md = MeshDistance::new(&coin_of(1200));
        let (origin, h, dims) = ([-0.4019, -0.3987, 0.1013], 0.05, [16, 16, 16]);
        let (values, unclear) = md.sample(origin, h, dims);
        assert_eq!(unclear, 0);
        for (n, v) in values.iter().enumerate() {
            let [i, j, k] = [n % dims[0], n / dims[0] % dims[1], n / (dims[0] * dims[1])];
            let p = [origin[0] + h * i as f64, origin[1] + h * j as f64, origin[2] + h * k as f64];
            assert_eq!(v.to_bits(), md.signed(p).to_bits(), "at {p:?}");
        }
    }

    /// When a block is tiny relative to its distance from the mesh (`r << d`),
    /// scaling the reach rounding allowance with `d + r` keeps tied triangles
    /// from dropping out of the block's gather by rounding.
    #[test]
    fn a_blocks_reach_keeps_a_tied_triangle_on_a_tiny_grid() {
        let tet = |apex: V3, sign: f64| -> [[V3; 3]; 4] {
            let b = [apex[0] + sign * 1.0, apex[1] + sign * 0.2, apex[2] + sign * 0.2];
            let c = [apex[0] + sign * 0.2, apex[1] + sign * 1.0, apex[2] + sign * 0.2];
            let d = [apex[0] + sign * 0.2, apex[1] + sign * 0.2, apex[2] + sign * 1.0];
            if sign < 0.0 { [[apex, b, c], [apex, c, d], [apex, d, b], [b, d, c]] } else { [[apex, c, b], [apex, d, c], [apex, b, d], [b, c, d]] }
        };
        let mut tris = Vec::new();
        tris.extend_from_slice(&tet([1.5658009414360923, -2.827669168683515, -2.0668825626934093], -1.0));
        tris.extend_from_slice(&tet([3.6435262642795645, -0.7499438458400433, 0.010842760150062785], 1.0));
        let md = MeshDistance::new(&tris);
        let (origin, h, dims) = ([2.604663576492978, -1.7888065336266292, -1.0280199276365236], 8.788283415904935e-9, [4, 4, 4]);
        let (grid, _) = md.sample(origin, h, dims);
        for (n, v) in grid.iter().enumerate() {
            let p = [origin[0] + h * (n % 4) as f64, origin[1] + h * (n / 4 % 4) as f64, origin[2] + h * (n / 16) as f64];
            assert_eq!(v.to_bits(), md.signed(p).to_bits(), "at {p:?}: {v}, signed {}", md.signed(p));
        }
    }

    /// On a coin turned off all three axes, the oriented node boxes never
    /// prune a closer triangle: `closest` matches a brute-force scan over
    /// every triangle to the bit.
    #[test]
    fn a_rotated_coins_closest_matches_a_linear_scan() {
        let (ca, sa) = (0.73_f64.cos(), 0.73_f64.sin());
        let (cb, sb) = (1.11_f64.cos(), 1.11_f64.sin());
        let rot = |p: V3| -> V3 {
            let [x, y, z] = [p[0] * ca - p[1] * sa, p[0] * sa + p[1] * ca, p[2]];
            [x, y * cb - z * sb, y * sb + z * cb]
        };
        let tris: Vec<[V3; 3]> = coin_of(120).into_iter().map(|t| t.map(rot)).collect();
        let md = MeshDistance::new(&tris);
        for k in 0..64 {
            let a = 2.0 * std::f64::consts::PI * (k as f64 + 0.37) / 64.0;
            let p = rot([6.3 * a.cos(), 6.3 * a.sin(), if k % 2 == 0 { 1.21 } else { -0.35 }]);
            let got = md.closest(p).0;
            let brute = tris
                .iter()
                .map(|v| {
                    let q = closest_on_triangle(p, v[0], v[1], v[2]).0;
                    dot(sub(p, q), sub(p, q))
                })
                .fold(f64::INFINITY, f64::min);
            assert_eq!(got.to_bits(), brute.to_bits(), "at {p:?}: tree {got}, brute {brute}");
        }
    }

    /// Above a coin's face, away from its centre, a search of the tree
    /// measures 41 triangles a point: its nodes' turned boxes are as thin
    /// as the slivers across the axes, whose own boxes hold the point
    /// below it (by those alone it measured 168).
    #[test]
    fn a_search_beside_slanted_slivers_measures_few_triangles() {
        let md = MeshDistance::new(&coin());
        let points: Vec<V3> = (0..64)
            .map(|k| {
                let a = 2.0 * std::f64::consts::PI * (k as f64 + 0.37) / 64.0;
                [6.3 * a.cos(), 6.3 * a.sin(), 1.21]
            })
            .collect();
        MEASURED.with(|m| m.set(0));
        for p in &points {
            assert!((md.closest(*p).0.sqrt() - 0.21).abs() < 1e-9);
        }
        let measured = MEASURED.with(|m| m.get());
        assert!(measured < 60 * points.len(), "{measured} triangles measured for {} points", points.len());
    }
}
