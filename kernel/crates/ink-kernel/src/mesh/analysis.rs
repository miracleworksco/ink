use ink_core::V3;
use ink_core::v3::{EMPTY_BOX, box_union, cross, dot, norm, normal};
use rayon::prelude::*;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct Mesh {
    pub verts: Vec<V3>,
    pub tris: Vec<[u32; 3]>,
}

impl Mesh {
    /// The sum of its triangles' areas.
    pub fn area(&self) -> f64 {
        self.tris.iter().map(|t| 0.5 * norm(normal(t.map(|i| self.verts[i as usize])))).sum()
    }
}

#[derive(Clone, Debug)]
pub struct Shell {
    pub triangles: usize,
    pub volume: f64,
    pub area: f64,
    pub lo: V3,
    pub hi: V3,
}

#[derive(Clone, Debug)]
pub struct Analysis {
    pub area: f64,
    pub volume: f64,
    pub boundary_edges: usize,
    pub nonmanifold_edges: usize,
    pub shells: Vec<Shell>,
}

fn tri_terms(m: &Mesh, t: &[u32; 3]) -> (f64, f64) {
    let (a, b, c) = (m.verts[t[0] as usize], m.verts[t[1] as usize], m.verts[t[2] as usize]);
    let area = 0.5 * norm(normal([a, b, c]));
    let vol = dot(a, cross(b, c)) / 6.0;
    (area, vol)
}

/// The root of `x` in a union-find forest of parent links, halving the
/// path to it on the way.
pub fn find(p: &mut [u32], mut x: u32) -> u32 {
    while p[x as usize] != x {
        p[x as usize] = p[p[x as usize] as usize];
        x = p[x as usize];
    }
    x
}

pub fn analyze(m: &Mesh) -> Analysis {
    let mut keys: Vec<u64> = m
        .tris
        .par_iter()
        .flat_map_iter(|t| {
            (0..3).map(move |e| {
                let (a, b) = (t[e], t[(e + 1) % 3]);
                ((a.min(b) as u64) << 32) | a.max(b) as u64
            })
        })
        .collect();
    keys.par_sort_unstable();
    let (mut boundary_edges, mut nonmanifold_edges) = (0, 0);
    for run in keys.chunk_by(|a, b| a == b) {
        match run.len() {
            1 => boundary_edges += 1,
            2 => {}
            _ => nonmanifold_edges += 1,
        }
    }

    // Union-find over vertices to split into connected shells.
    let mut parent: Vec<u32> = (0..m.verts.len() as u32).collect();
    for t in &m.tris {
        let r0 = find(&mut parent, t[0]);
        for &v in &t[1..] {
            let r = find(&mut parent, v);
            parent[r as usize] = r0;
        }
    }
    let mut shells: HashMap<u32, Shell> = HashMap::new();
    let (mut area, mut volume) = (0.0, 0.0);
    for t in &m.tris {
        let (a, v) = tri_terms(m, t);
        area += a;
        volume += v;
        let r = find(&mut parent, t[0]);
        let (lo, hi) = EMPTY_BOX;
        let s = shells.entry(r).or_insert(Shell { triangles: 0, volume: 0.0, area: 0.0, lo, hi });
        s.triangles += 1;
        s.area += a;
        s.volume += v;
        for &i in t {
            let p = m.verts[i as usize];
            (s.lo, s.hi) = box_union((s.lo, s.hi), (p, p));
        }
    }
    let mut shells: Vec<Shell> = shells.into_values().collect();
    shells.sort_by(|a, b| b.volume.abs().total_cmp(&a.volume.abs()));
    Analysis { area, volume, boundary_edges, nonmanifold_edges, shells }
}
