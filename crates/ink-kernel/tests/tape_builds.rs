//! Tapes a query builds, counted: a process of its own, so no other test's
//! builds are in the count (Tape::built is the process's).

use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, Graph, NodeId, UnOp, Var};
use ink_core::tape::Tape;
use ink_kernel::integral::{ACCURACY, Until, integrate};
use ink_kernel::octree::{descend, root_for};
use std::sync::Arc;

/// A plate 40 x 40 x 6 with four holes of radius 3: the max of three
/// slabs and the holes' complement, as a part's model writes it.
fn plate(g: &mut Graph) -> NodeId {
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let slab = |g: &mut Graph, v: NodeId, h: f64| {
        let a = g.un(UnOp::Abs, v);
        let h = g.constant(h);
        g.bin(BinOp::Sub, a, h)
    };
    let (sx, sy, sz) = (slab(g, x, 20.0), slab(g, y, 20.0), slab(g, z, 3.0));
    let plate = g.bin(BinOp::Max, sx, sy);
    let plate = g.bin(BinOp::Max, plate, sz);
    let mut holes = None;
    for (cx, cy) in [(-10.0, -10.0), (10.0, -10.0), (-10.0, 10.0), (10.0, 10.0)] {
        let (cx, cy) = (g.constant(cx), g.constant(cy));
        let (dx, dy) = (g.bin(BinOp::Sub, x, cx), g.bin(BinOp::Sub, y, cy));
        let (dx, dy) = (g.un(UnOp::Square, dx), g.un(UnOp::Square, dy));
        let r = g.bin(BinOp::Add, dx, dy);
        let r = g.un(UnOp::Sqrt, r);
        let three = g.constant(3.0);
        let hole = g.bin(BinOp::Sub, r, three);
        holes = Some(holes.map_or(hole, |h| g.bin(BinOp::Min, h, hole)));
    }
    let out = g.un(UnOp::Neg, holes.unwrap());
    g.bin(BinOp::Max, plate, out)
}

/// A mass, measured to the accuracy a published mass is, and the plate's
/// cells down to a quarter of a millimetre (as integrals' estimates on
/// cells take them): the tapes its boxes and cells are simplified to, and
/// their operands' (integral::shares), are each built once, not once per
/// box or cell that reads them. One test, so the count is one query's.
#[test]
fn a_query_builds_each_tape_once() {
    let mut g = Graph::new();
    let f = plate(&mut g);
    let shape = Arc::new(Tape::new(&g, f));
    let ctx = Ctx::new(0.0, vec![]);
    let before = Tape::built();
    let m = integrate(&shape, &[None], true, &ctx, [-21.0; 3], [21.0; 3], Until::Accuracy(ACCURACY));
    let built = Tape::built() - before;
    assert!(m.settled, "{m:?}");
    // Each built once: a few hundred tapes, against one or more a box.
    assert!(built * 10 < m.boxes, "{built} tapes built for {} boxes", m.boxes);

    let (grid, root) = root_for([-21.0; 3], [21.0; 3], 0.25, 1);
    let before = Tape::built();
    let cells = descend(&shape, &ctx, &grid, root, 1, &|_, c| Some(c), &|_, c| Some(c)).len();
    let built = Tape::built() - before;
    // A few dozen, against one for every dozen cells.
    assert!(built * 100 < cells, "{built} tapes built for {cells} cells");
}
