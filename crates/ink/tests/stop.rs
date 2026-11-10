//! The call's time limit (ink::stop) holds for the whole process, so its
//! tests are a binary of their own: a timed mesher and a mesh's distance
//! tree check it in every step, so neither answers past it.

use ink::eval::Ctx;
use ink::graph::{BinOp, Graph, UnOp, Var};
use ink::tape::Tape;
use std::sync::Arc;
use std::time::Instant;

/// `run` stopped at small shares of its quickest of three unstopped runs
/// gives nothing each time: no run finishes several times faster than its
/// quickest, so each stop falls before its end (no clock is compared to
/// the stop but the limit's own).
fn stops_every_time<T>(what: &str, run: impl Fn() -> Option<T>) {
    let now = || std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_secs_f64();
    // A limit far off, so stop::within sets each run's own.
    const FAR: f64 = 1e6;
    ink::stop::at(now() + FAR);
    let best = (0..3)
        .map(|_| {
            let t0 = Instant::now();
            assert!(run().is_some(), "{what}: no stop, an answer");
            t0.elapsed().as_secs_f64()
        })
        .fold(f64::INFINITY, f64::min);
    for share in [0.15, 0.25] {
        let (r, _) = ink::stop::within(share * best / FAR, &run);
        assert!(r.is_none(), "{what}: an answer, stopped at {share} of {best:.3} s");
    }
}

/// A sphere of radius 10, meshed at 0.05: about a million triangles.
fn sphere() -> (Arc<Tape>, Ctx) {
    let mut g = Graph::new();
    let mut s = None;
    for v in [Var::X, Var::Y, Var::Z] {
        let a = g.var(v);
        let q = g.un(UnOp::Square, a);
        s = Some(s.map_or(q, |s| g.bin(BinOp::Add, s, q)));
    }
    let r = g.un(UnOp::Sqrt, s.unwrap());
    let ten = g.constant(10.0);
    let f = g.bin(BinOp::Sub, r, ten);
    (Arc::new(Tape::new(&g, f)), Ctx::new(0.0, vec![]))
}

#[test]
fn a_timed_mesh_and_its_distance_tree_stop_at_the_limit() {
    let (tape, ctx) = sphere();
    let mesh = || ink::mesh::mesh_in_time(&tape, &ctx, [-10.5; 3], [10.5; 3], 0.05);
    stops_every_time("the mesh", mesh);
    let m = ink::mesh::mesh(&tape, &ctx, [-10.5; 3], [10.5; 3], 0.05);
    let tris: Vec<[ink::V3; 3]> = m.tris.iter().map(|t| t.map(|i| m.verts[i as usize])).collect();
    assert!(tris.len() > 500_000, "{} triangles", tris.len());
    stops_every_time("the distance tree", || ink::grid::MeshDistance::in_time(&tris));
}
