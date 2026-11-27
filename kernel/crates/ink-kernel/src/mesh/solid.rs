use ink_core::V3;
use ink_core::v3::{cross, dot, norm};
use rayon::prelude::*;
use std::collections::HashMap;

/// A closed mesh's triangles by the squares of side `s` their extents
/// across y and z meet, to find what a ray along +x crosses.
pub(super) struct Solid {
    s: f64,
    squares: HashMap<[i64; 2], Vec<u32>>,
    cubes: HashMap<[i64; 3], Vec<u32>>,
}

impl Solid {
    /// None unless the mesh is closed: each edge two triangles'.
    pub(super) fn new(verts: &[V3], tris: &[[u32; 3]], s: f64) -> Option<Solid> {
        let mut edges: HashMap<(u32, u32), u32> = HashMap::with_capacity(3 * tris.len());
        for t in tris {
            for e in 0..3 {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
        if tris.is_empty() || edges.values().any(|&n| n != 2) {
            return None;
        }
        let (mut squares, mut cubes): (HashMap<[i64; 2], Vec<u32>>, HashMap<[i64; 3], Vec<u32>>) = (HashMap::new(), HashMap::new());
        for (k, t) in tris.iter().enumerate() {
            let p = t.map(|i| verts[i as usize]);
            let span = |a: usize| {
                let (lo, hi) = (p.iter().map(|v| v[a]).fold(f64::INFINITY, f64::min), p.iter().map(|v| v[a]).fold(f64::NEG_INFINITY, f64::max));
                ((lo / s).floor() as i64)..=((hi / s).floor() as i64)
            };
            for y in span(1) {
                for z in span(2) {
                    squares.entry([y, z]).or_default().push(k as u32);
                    for x in span(0) {
                        cubes.entry([x, y, z]).or_default().push(k as u32);
                    }
                }
            }
        }
        Some(Solid { s, squares, cubes })
    }

    /// Whether `p` lies inside the solid (a ray from it along +x crosses
    /// its triangles an odd number of times); None where the ray meets an
    /// edge or vertex, or `p` lies on a triangle, as rounding can tell.
    pub(super) fn inside(&self, verts: &[V3], tris: &[[u32; 3]], p: V3) -> Option<bool> {
        let key = [(p[1] / self.s).floor() as i64, (p[2] / self.s).floor() as i64];
        let mut odd = false;
        for &k in self.squares.get(&key).map_or(&[][..], Vec::as_slice) {
            let [a, b, c] = tris[k as usize].map(|i| verts[i as usize]);
            // Twice the signed areas, across y and z, of p with each edge.
            let side = |u: V3, v: V3| (u[1] - p[1]) * (v[2] - p[2]) - (u[2] - p[2]) * (v[1] - p[1]);
            let (d0, d1, d2) = (side(b, c), side(c, a), side(a, b));
            let scale = [a, b, c].iter().flat_map(|v| [(v[1] - p[1]).abs(), (v[2] - p[2]).abs()]).fold(0.0f64, f64::max);
            let eps = 1e-9 * scale * scale;
            let all = d0 + d1 + d2;
            if (d0 > eps && d1 > eps && d2 > eps) || (d0 < -eps && d1 < -eps && d2 < -eps) {
                let x = (d0 * a[0] + d1 * b[0] + d2 * c[0]) / all;
                if (x - p[0]).abs() <= 1e-9 * (x.abs() + p[0].abs() + self.s) {
                    return None;
                }
                odd ^= x > p[0];
            } else if !(d0 < -eps || d1 < -eps || d2 < -eps) || !(d0 > eps || d1 > eps || d2 > eps) {
                // On an edge or vertex across y and z (or the triangle is
                // edge-on there): only if p is not plainly outside it.
                return None;
            }
        }
        Some(odd)
    }

    /// Whether a triangle lies within `r` of `p` (r taken no more than `s`).
    pub(super) fn near(&self, verts: &[V3], tris: &[[u32; 3]], p: V3, r: f64) -> bool {
        let r = r.min(self.s);
        if r <= 0.0 {
            return false;
        }
        let c = p.map(|x| (x / self.s).floor() as i64);
        (0..27).any(|q| {
            let key = [c[0] + q % 3 - 1, c[1] + q / 3 % 3 - 1, c[2] + q / 9 - 1];
            self.cubes.get(&key).is_some_and(|ks| ks.iter().any(|&k| norm(sub3(closest_on_triangle(p, tris[k as usize].map(|i| verts[i as usize])), p)) <= r))
        })
    }
}

/// The cubes `s` across that triangle `p` may pass through, each once: in
/// each column along the axis its normal is nearest, the cubes between its
/// plane's heights over the column, of the columns its shadow meets. As
/// many as its area takes, not its box.
fn cubes_of(p: [V3; 3], s: f64, mut f: impl FnMut([i64; 3])) {
    let at = |x: f64| (x / s).floor() as i64;
    let lo: V3 = std::array::from_fn(|a| p[0][a].min(p[1][a]).min(p[2][a]));
    let hi: V3 = std::array::from_fn(|a| p[0][a].max(p[1][a]).max(p[2][a]));
    let n = cross(sub3(p[1], p[0]), sub3(p[2], p[0]));
    let w = (0..3).max_by(|&a, &b| n[a].abs().total_cmp(&n[b].abs())).unwrap();
    let (u, v) = ((w + 1) % 3, (w + 2) % 3);
    let cell = |iu: i64, iv: i64, iw: i64| {
        let mut c = [0; 3];
        (c[u], c[v], c[w]) = (iu, iv, iw);
        c
    };
    if !(n[w].abs() > 0.0) {
        for iu in at(lo[u])..=at(hi[u]) {
            for iv in at(lo[v])..=at(hi[v]) {
                (at(lo[w])..=at(hi[w])).for_each(|iw| f(cell(iu, iv, iw)));
            }
        }
        return;
    }
    // The plane's height over (x, y) along w.
    let height = |x: f64, y: f64| p[0][w] - (n[u] * (x - p[0][u]) + n[v] * (y - p[0][v])) / n[w];
    for iu in at(lo[u])..=at(hi[u]) {
        let (a, b) = (iu as f64 * s, (iu + 1) as f64 * s);
        // The shadow's v over this strip: its corners in it and where its
        // edges cross the strip's sides.
        let (mut v0, mut v1) = (f64::INFINITY, f64::NEG_INFINITY);
        for k in 0..3 {
            let (q, r) = (p[k], p[(k + 1) % 3]);
            if (a..=b).contains(&q[u]) {
                (v0, v1) = (v0.min(q[v]), v1.max(q[v]));
            }
            for side in [a, b] {
                if (q[u] - side) * (r[u] - side) < 0.0 {
                    let y = q[v] + (r[v] - q[v]) * (side - q[u]) / (r[u] - q[u]);
                    (v0, v1) = (v0.min(y), v1.max(y));
                }
            }
        }
        if v0 > v1 {
            continue;
        }
        for iv in at(v0)..=at(v1) {
            let (c, d) = (iv as f64 * s, (iv + 1) as f64 * s);
            let hs = [height(a, c), height(a, d), height(b, c), height(b, d)];
            let (h0, h1) = (hs.iter().fold(f64::INFINITY, |m, &x| m.min(x)).max(lo[w]), hs.iter().fold(f64::NEG_INFINITY, |m, &x| m.max(x)).min(hi[w]));
            (at(h0)..=at(h1)).for_each(|iw| f(cell(iu, iv, iw)));
        }
    }
}

fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// The point of triangle `t` nearest `p` (Ericson, Real-Time Collision
/// Detection, 5.1.5).
fn closest_on_triangle(p: V3, t: [V3; 3]) -> V3 {
    let [a, b, c] = t;
    let (ab, ac, ap) = (sub3(b, a), sub3(c, a), sub3(p, a));
    let at = |u: f64, v: f64| -> V3 { std::array::from_fn(|k| a[k] + u * ab[k] + v * ac[k]) };
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub3(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return at(d1 / (d1 - d3), 0.0);
    }
    let cp = sub3(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return at(0.0, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return std::array::from_fn(|k| b[k] + w * (c[k] - b[k]));
    }
    let denom = 1.0 / (va + vb + vc);
    at(vb * denom, vc * denom)
}

/// Vertices placed by `qef` may make triangles that pass through others
/// where sheets come within a cell of each other: the slack past the cell,
/// or two patches of one cell placed across each other. Each vertex of a
/// triangle that crosses another goes to its crossings' mean inside its
/// cell (`safe`), round after round while that changes one; sharp edges
/// and corners elsewhere keep their place.
pub(super) fn untangle(verts: &mut [V3], safe: &[V3], tris: &[[u32; 3]], h: f64) {
    for _ in 0..4 {
        let mut moved = false;
        for t in crossing(verts, tris, 2.0 * h) {
            for v in tris[t as usize] {
                let v = v as usize;
                if v < safe.len() && verts[v] != safe[v] {
                    (verts[v], moved) = (safe[v], true);
                }
            }
        }
        if !moved {
            return;
        }
    }
}

/// The triangles that pass through another sharing no vertex with them:
/// an edge of one through the other's inside. Found in a grid of cubes `s`
/// across, each triangle in those it passes through (a pair in several
/// tested in each).
pub fn crossing(verts: &[V3], tris: &[[u32; 3]], s: f64) -> Vec<u32> {
    let bounds = |t: &[u32; 3]| {
        let p = t.map(|i| verts[i as usize]);
        let lo: V3 = std::array::from_fn(|a| p[0][a].min(p[1][a]).min(p[2][a]));
        let hi: V3 = std::array::from_fn(|a| p[0][a].max(p[1][a]).max(p[2][a]));
        (lo, hi)
    };
    let mut placed: Vec<([i64; 3], u32)> = tris
        .par_iter()
        .enumerate()
        .flat_map_iter(|(k, t)| {
            let mut v = Vec::new();
            cubes_of(t.map(|i| verts[i as usize]), s, |c| v.push((c, k as u32)));
            v
        })
        .collect();
    placed.par_sort_unstable();
    let through = |p: V3, q: V3, t: [V3; 3]| {
        let n = cross(sub3(t[1], t[0]), sub3(t[2], t[0]));
        let (s1, s2) = (dot(sub3(p, t[0]), n), dot(sub3(q, t[0]), n));
        if s1 * s2 >= 0.0 {
            return false;
        }
        let f = s1 / (s1 - s2);
        let x: V3 = std::array::from_fn(|a| p[a] + (q[a] - p[a]) * f);
        let side = |a: V3, b: V3| dot(cross(sub3(b, a), sub3(x, a)), n);
        let (u, v, w) = (side(t[0], t[1]), side(t[1], t[2]), side(t[2], t[0]));
        (u > 0.0 && v > 0.0 && w > 0.0) || (u < 0.0 && v < 0.0 && w < 0.0)
    };
    let cells: Vec<Vec<u32>> = placed.chunk_by(|a, b| a.0 == b.0).map(|r| r.iter().map(|p| p.1).collect()).collect();
    let mut out: Vec<u32> = cells
        .par_iter()
        .flat_map_iter(|ks| {
            let mut found = Vec::new();
            for (x, &a) in ks.iter().enumerate() {
                let (ta, (la, ha)) = (tris[a as usize], bounds(&tris[a as usize]));
                for &b in &ks[x + 1..] {
                    let tb = tris[b as usize];
                    if ta.iter().any(|v| tb.contains(v)) {
                        continue;
                    }
                    let (lb, hb) = bounds(&tb);
                    if (0..3).any(|k| ha[k] < lb[k] || hb[k] < la[k]) {
                        continue;
                    }
                    let (pa, pb) = (ta.map(|i| verts[i as usize]), tb.map(|i| verts[i as usize]));
                    if (0..3).any(|k| through(pa[k], pa[(k + 1) % 3], pb)) || (0..3).any(|k| through(pb[k], pb[(k + 1) % 3], pa)) {
                        found.extend([a, b]);
                    }
                }
            }
            found
        })
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}
