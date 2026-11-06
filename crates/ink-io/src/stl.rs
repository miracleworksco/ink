//! Meshes to fields: STL reading and the signed distance to a triangle
//! mesh, in floats.

use ink_core::V3;
use ink_core::v3::{dot, norm, normal, sub};

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

/// Closest point on triangle (Ericson), and which feature: 0 face, 1..=3
/// vertex a/b/c, 4..=6 edge ab/bc/ca.
fn closest_on_triangle(p: V3, a: V3, b: V3, c: V3) -> (V3, u8) {
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

/// Distance to a triangle mesh, in floats: signed by angle-weighted
/// pseudo-normals (`distance`), and for a closed mesh by parity, the side
/// of the solid it encloses (`signed`, `sample`).
pub struct MeshDistance {
    verts: Vec<V3>,
    tris: Vec<[u32; 3]>,
    face_n: Vec<V3>,
    vert_n: Vec<V3>,
    edge_n: std::collections::HashMap<(u32, u32), V3>,
    // AABB tree: nodes (lo, hi, left/first, right/count, is_leaf)
    nodes: Vec<(V3, V3, u32, u32, bool)>,
    order: Vec<u32>,
    /// Every edge in an even number of triangles (`odd_edges`).
    closed: bool,
}

impl MeshDistance {
    pub fn new(tris_in: &[[V3; 3]]) -> MeshDistance {
        // Weld identical positions into shared vertices.
        let mut index: std::collections::HashMap<[u64; 3], u32> = std::collections::HashMap::new();
        let mut verts = Vec::new();
        let mut tris = Vec::with_capacity(tris_in.len());
        for t in tris_in {
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
        let mut edge_n: std::collections::HashMap<(u32, u32), V3> = std::collections::HashMap::new();
        for (t, n) in tris.iter().zip(&face_n) {
            for k in 0..3 {
                let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let (e1, e2) = (sub(verts[b as usize], verts[a as usize]), sub(verts[c as usize], verts[a as usize]));
                let ang = libm::acos((dot(e1, e2) / (norm(e1) * norm(e2)).max(1e-300)).clamp(-1.0, 1.0));
                for d in 0..3 {
                    vert_n[a as usize][d] += ang * n[d];
                }
                let e = edge_n.entry((a.min(b), a.max(b))).or_insert([0.0; 3]);
                for d in 0..3 {
                    e[d] += n[d];
                }
            }
        }
        let mut md = MeshDistance { verts, tris, face_n, vert_n, edge_n, nodes: Vec::new(), order: Vec::new(), closed: false };
        md.closed = md.odd_edges() == 0;
        md.order = (0..md.tris.len() as u32).collect();
        let n = md.order.len();
        md.build(0, n);
        md
    }

    fn tri_box(&self, t: u32) -> (V3, V3) {
        let v = self.tris[t as usize].map(|i| self.verts[i as usize]);
        (std::array::from_fn(|a| v[0][a].min(v[1][a]).min(v[2][a])), std::array::from_fn(|a| v[0][a].max(v[1][a]).max(v[2][a])))
    }

    fn build(&mut self, start: usize, end: usize) -> u32 {
        let (lo, hi) = self.order[start..end].iter().map(|&t| self.tri_box(t)).fold(ink_core::v3::EMPTY_BOX, ink_core::v3::box_union);
        let id = self.nodes.len() as u32;
        if end - start <= 4 {
            self.nodes.push((lo, hi, start as u32, (end - start) as u32, true));
            return id;
        }
        self.nodes.push((lo, hi, 0, 0, false));
        let axis = (0..3).max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b]))).unwrap();
        let mid = (start + end) / 2;
        let key = |t: &u32, md: &MeshDistance| {
            let (a, b) = md.tri_box(*t);
            a[axis] + b[axis]
        };
        let mut slice: Vec<u32> = self.order[start..end].to_vec();
        slice.sort_by(|x, y| key(x, self).total_cmp(&key(y, self)));
        self.order[start..end].copy_from_slice(&slice);
        let l = self.build(start, mid);
        let r = self.build(mid, end);
        self.nodes[id as usize].2 = l;
        self.nodes[id as usize].3 = r;
        id
    }

    /// The triangle nearest p, in floats: (squared distance, nearest point,
    /// triangle, feature as `closest_on_triangle`).
    fn closest(&self, p: V3) -> (f64, V3, u32, u8) {
        self.closest_within(p, f64::INFINITY)
    }

    /// `closest`, given that a triangle lies within `reach` of p (only
    /// nearer ones are looked at).
    fn closest_within(&self, p: V3, reach: f64) -> (f64, V3, u32, u8) {
        let mut best = (reach * reach * (1.0 + 1e-9), [0.0; 3], 0u32, 0u8);
        let mut stack = Vec::with_capacity(64);
        stack.push(0u32);
        let box_d2 = |lo: V3, hi: V3| -> f64 { (0..3).map(|a| (lo[a] - p[a]).max(0.0).max(p[a] - hi[a]).powi(2)).sum() };
        while let Some(n) = stack.pop() {
            let (lo, hi, a, b, leaf) = self.nodes[n as usize];
            if box_d2(lo, hi) >= best.0 {
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
                let (dl, dr) = (box_d2(self.nodes[a as usize].0, self.nodes[a as usize].1), box_d2(self.nodes[b as usize].0, self.nodes[b as usize].1));
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

    /// Signed distance at p (negative inside), in floats, signed by the
    /// nearest feature's pseudo-normal: exact in real arithmetic for a
    /// closed, consistently oriented mesh facing out that does not cross
    /// itself.
    pub fn distance(&self, p: V3) -> f64 {
        let (d2, q, t, feat) = self.closest(p);
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
                self.edge_n[&(a.min(b), a.max(b))]
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
        let mut edges: std::collections::HashMap<(u32, u32), u32> = std::collections::HashMap::with_capacity(3 * self.tris.len());
        for t in &self.tris {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        edges.values().filter(|&&n| n % 2 == 1).count()
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
    /// a line of it at a time, and how many lines were signed by
    /// `distance` for being unclear (none for a mesh that is not closed:
    /// it is signed so throughout).
    pub fn sample(&self, origin: V3, h: f64, dims: [usize; 3]) -> (Vec<f64>, usize) {
        use rayon::prelude::*;
        let mut values = vec![0.0; dims[0] * dims[1] * dims[2]];
        let unclear = values
            .par_chunks_mut(dims[0])
            .enumerate()
            .map(|(r, line)| {
                let (y, z) = (origin[1] + h * (r % dims[1]) as f64, origin[2] + h * (r / dims[1]) as f64);
                let xs = self.crossings(y, z);
                // Each point is within h of the last one, so the mesh within
                // h more than that was.
                let mut reach = f64::INFINITY;
                for (i, v) in line.iter_mut().enumerate() {
                    *v = self.signed_on(xs.as_deref(), [origin[0] + h * i as f64, y, z], reach);
                    reach = v.abs() + h;
                }
                (self.closed && xs.is_none()) as usize
            })
            .sum();
        (values, unclear)
    }
}
