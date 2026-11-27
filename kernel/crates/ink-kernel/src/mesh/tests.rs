use super::brick::{B, Brick, brick};
use super::contour::split_quads;
use super::proof::{sheet, split};
use super::report::Ran;
use super::{Limit, mesh_checked};
use crate::octree::{Cell, Grid, Local};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_points};
use ink_core::graph::{Graph, UnOp, Var};
use ink_core::interval::I;
use ink_core::tape::{Pass, Sheet, Tape};
use std::sync::Arc;

/// The brick at the origin's cell 0 (`[0, 1]^3` at h = 1), and whether
/// it is shown drawn (not left open).
fn cell_0(g: &Graph, f: ink_core::graph::NodeId) -> (Brick, bool) {
    let (tape, ctx) = (Tape::new(g, f), Ctx::new(0.0, vec![]));
    let b = brick(&Arc::new(tape.clone()), &tape, &ctx, &Grid::new([0.0; 3], 1.0), Cell { i: [0; 3], size: B }, None, (&Default::default(), usize::MAX), Some(&Ran::default()))
        .unwrap();
    let drawn = !b.open.iter().any(|o| o.0 == 0);
    (b, drawn)
}

#[test]
fn a_cell_is_shown_drawn_over_all_its_vertex_may_reach() {
    // z + max(|x - .5|, |y - .5|) - 1.2, a cone: one sheet round cell 0,
    // whose vertex is placed above it (z 1.19): it lies in the box
    // shown, the cell grown by half a cell.
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let (half, top) = (g.constant(0.5), g.constant(1.2));
    let (dx, dy) = (g.sub(x, half), g.sub(y, half));
    let (ax, ay) = (g.un(UnOp::Abs, dx), g.un(UnOp::Abs, dy));
    let m = g.max(ax, ay);
    let s = g.add(z, m);
    let f = g.sub(s, top);
    let (b, drawn) = cell_0(&g, f);
    let p = b.patches.iter().find(|q| q.0 == [0; 3]).unwrap().2;
    assert!(drawn && p.iter().all(|&v| (-0.5..=1.5).contains(&v)), "{p:?}");
    // (z - .5) (1 - 8 max(x - 1.05, 0)): a plane in cell 0, turning over
    // past x = 1.175, within half a cell of it. One sheet in the cell,
    // not round it: not shown drawn.
    let mut g = Graph::new();
    let [x, z] = [Var::X, Var::Z].map(|v| g.var(v));
    let (half, at, zero, eight, one) = (g.constant(0.5), g.constant(1.05), g.constant(0.0), g.constant(8.0), g.constant(1.0));
    let (dz, dx) = (g.sub(z, half), g.sub(x, at));
    let past = g.max(dx, zero);
    let k = g.mul(eight, past);
    let k = g.sub(one, k);
    let f = g.mul(dz, k);
    assert!(!cell_0(&g, f).1);
}

#[test]
fn a_crossed_cells_box_holding_a_piece_besides_its_sheet_is_not_shown_drawn() {
    // A plane, z = 0.25, crossing cell 0, and in its box a double cone
    // from x 0.3 to 0.7 at y 0.5, z 1.25, widest (0.05) at x 0.5: a
    // piece of its own beside the sheet. Not shown drawn. (A guard: a
    // rule accepting the box part by part is argued unsound here, a
    // closed piece across parts each one sheet along its own way, but
    // this cone is not accepted by such a rule either.)
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let c = |g: &mut Graph, v: f64| g.constant(v);
    let (q, h) = (c(&mut g, 0.25), c(&mut g, 0.5));
    let plane = g.sub(z, q);
    let (dy, zc) = (g.sub(y, h), c(&mut g, 1.25));
    let dz = g.sub(z, zc);
    let (yy, zz) = (g.un(UnOp::Square, dy), g.un(UnOp::Square, dz));
    let s = g.add(yy, zz);
    let radial = g.un(UnOp::Sqrt, s);
    let dx = g.sub(x, h);
    let ax = g.un(UnOp::Abs, dx);
    let (k, w) = (c(&mut g, 0.25), c(&mut g, 0.05));
    let taper = g.mul(k, ax);
    let r = g.sub(w, taper);
    let cone = g.sub(radial, r);
    let f = g.min(plane, cone);
    assert!(!cell_0(&g, f).1);
}

#[test]
fn a_comparison_that_jumps_past_equality_is_a_seam() {
    // select(a - y, x + 0.2 a^2 + z^2 - 0.05, x + 0.2 y^2 + z^2 - 0.05),
    // a = 2 floor(x) + 1: the sides agree where a = y, but a jumps from
    // -1 to 1 at x = 0, never meeting y, and the field jumps there,
    // a closed piece 0 <= x < 0.05 beside the sheet near x = -0.15.
    // Cell [-0.25, 0.25]^3 at h 0.5 is not shown drawn.
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let (offset, two, one, k) = (g.constant(0.05), g.constant(2.0), g.constant(1.0), g.constant(0.2));
    let fl = g.un(UnOp::Floor, x);
    let twice = g.mul(two, fl);
    let a = g.add(twice, one);
    let cond = g.sub(a, y);
    let (aa, bb, zz) = (g.un(UnOp::Square, a), g.un(UnOp::Square, y), g.un(UnOp::Square, z));
    let (aa, bb) = (g.mul(k, aa), g.mul(k, bb));
    let (xa, xb) = (g.add(x, aa), g.add(x, bb));
    let (xa, xb) = (g.add(xa, zz), g.add(xb, zz));
    let (left, right) = (g.sub(xa, offset), g.sub(xb, offset));
    let f = g.select(cond, left, right);
    let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
    let b = brick(&Arc::new(tape.clone()), &tape, &ctx, &Grid::new([-0.25; 3], 0.5), Cell { i: [0; 3], size: B }, None, (&Default::default(), usize::MAX), Some(&Ran::default()))
        .unwrap();
    let mut v = Vec::new();
    eval_points(&tape, &ctx, &[[-0.2, 0.0, 0.0], [-0.1, 0.0, 0.0], [0.025, 0.0, 0.0], [0.1, 0.0, 0.0]], &mut v);
    assert_eq!(v.iter().map(|v| *v < 0.0).collect::<Vec<_>>(), [true, false, true, false]);
    assert!(b.open.iter().any(|o| o.0 == 0));
}

#[test]
fn a_seams_reach_reads_the_boxs_bounds_once() {
    // select(x + y + z, y - 0.1, y - 0.2): it may switch along every
    // axis, so no axis shows one sheet, and its reach (0.1 over a rate
    // of 1) is found along the direction tried from the box's bounds
    // and slopes, found once, and the slopes at its centre: 11 passes
    // of the tape (bounds 1, slopes 3 and 3, four walks), not the
    // box's bounds and slopes again (12, counting the centre's as 1).
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let (lower, upper) = (g.constant(0.1), g.constant(0.2));
    let xy = g.add(x, y);
    let cond = g.add(xy, z);
    let (a, b) = (g.sub(y, lower), g.sub(y, upper));
    let f = g.select(cond, a, b);
    let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
    let mut read = 0;
    let (held, seam) = sheet(&tape, &ctx, [-0.5; 3], [0.5; 3], &mut read);
    assert!(seam && matches!(held, Sheet::Thin(t) if (0.1..0.11).contains(&t)), "{held:?}");
    assert!(read <= 11 * tape.work(Pass::Box), "{read} instructions, {} a pass", tape.work(Pass::Box));
}

#[test]
fn a_split_box_is_bounded_once_for_every_direction() {
    // max(-x, y - 0.25) over [-1, 1]^3, a face over the half x > 0,
    // split where y meets 0.25: asked along five directions at once
    // its halves are bounded once, each direction getting what it gets
    // alone, for under 3/5 of the instructions of asking each alone
    // (about half: 204 against 402).
    let mut g = Graph::new();
    let [x, y] = [Var::X, Var::Y].map(|v| g.var(v));
    let quarter = g.constant(0.25);
    let (nx, top) = (g.un(UnOp::Neg, x), g.sub(y, quarter));
    let f = g.max(nx, top);
    let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
    let b = [I::new(-1.0, 1.0); 3];
    let ls = ink_core::sign::levels(&tape).to_vec();
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let dirs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [s, s, 0.0], [s, -s, 0.0]];
    let (mut all, mut alone) = (0, 0);
    let together = split(&tape, &ctx, b, &dirs, ink_core::sign::SPLITS, &ls, &mut all);
    let apart: Vec<i8> = dirs.iter().map(|&d| split(&tape, &ctx, b, &[d], ink_core::sign::SPLITS, &ls, &mut alone)[0]).collect();
    assert_eq!(together, apart);
    assert!(together.iter().any(|&w| w != 0), "{together:?}");
    assert!(5 * all < 3 * alone, "{all} instructions at once, {alone} apart");
}

#[test]
fn a_seams_bound_is_kept_whole() {
    // At h = 1e6 a seam's reach is about 800000.03; a ball 800000.019
    // across fits in the slab it may hide. Its bound is no less.
    let h = 1e6;
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let c = |g: &mut Graph, v: f64| g.constant(v);
    let (lower, switch, upper, half, radius) = (c(&mut g, 0.1 * h), c(&mut g, 0.10000001 * h), c(&mut g, 0.90000003 * h), c(&mut g, 0.5 * h), c(&mut g, 0.45 * h));
    let (b, top) = (g.sub(y, lower), g.sub(y, upper));
    let (dx, dz) = (g.sub(x, half), g.sub(z, half));
    let (ax, az) = (g.un(UnOp::Abs, dx), g.un(UnOp::Abs, dz));
    let radial = g.max(ax, az);
    let radial = g.sub(radial, radius);
    let above = g.sub(switch, y);
    let cond = g.max(radial, above);
    let f = g.select(cond, top, b);
    let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
    let br =
        brick(&Arc::new(tape.clone()), &tape, &ctx, &Grid::new([0.0; 3], h), Cell { i: [0; 3], size: B }, None, (&Default::default(), usize::MAX), Some(&Ran::default())).unwrap();
    let k = br.open.iter().position(|o| o.0 == 0).unwrap();
    let (_, u) = mesh_checked(&Arc::new(tape), &ctx, [0.0; 3], [h; 3], h, Limit::GRID, 0.0).unwrap();
    assert!(br.thin[k] >= 0.800000019 * h && u.thinner_than >= 0.800000019 * h, "{} {}", br.thin[k], u.thinner_than);
}

#[test]
fn a_bricks_samples_run_the_field_of_its_own_box() {
    // z - 4.5 (off the grid's points: none on 0) and a branch, min'd
    // with it, that loses throughout the brick at the origin ([0, 8]^3
    // at h = 1) but may win in the half cell round it its tape holds:
    // the brick's samples run the field as simplified over its own box,
    // as handed that tape.
    let mut g = Graph::new();
    let (x, z) = (g.var(Var::X), g.var(Var::Z));
    let (four, k, eight) = (g.constant(4.0), g.constant(8.25), g.constant(8.0));
    let mid = g.constant(4.5);
    let f = g.sub(z, mid);
    let d = g.sub(k, x);
    let b = g.mul(d, eight);
    let b = g.add(b, four);
    let f = g.min(f, b);
    let (tape, ctx) = (Tape::new(&g, f), Ctx::new(0.0, vec![]));
    let (grid, c) = (Grid::new([0.0; 3], 1.0), Cell { i: [0; 3], size: B });
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(&tape, &ctx, grid.cell_box(c), Some(&mut ch));
    let own = tape.simplify(&ch);
    assert!(own.len() < tape.len(), "the branch is dropped over the brick");
    let ran = |t: &Tape| {
        let r = Ran::default();
        let b = brick(&Arc::new(t.clone()), &tape, &ctx, &grid, c, None, (&Default::default(), usize::MAX), Some(&r)).unwrap();
        (b.patches.len(), r.points.load(std::sync::atomic::Ordering::Relaxed))
    };
    assert_eq!(ran(&tape), ran(&own));
}

#[test]
fn a_walls_face_is_one_sheet_in_a_box_reaching_its_middle() {
    // |z| - w: cell 0's box, z from -0.5 to 1.5, reaches the wall's
    // middle (z = 0). With w = 0.6 it holds one face only: drawn. With
    // w = 0.3 it holds both: not shown drawn.
    let wall = |w: f64| {
        let mut g = Graph::new();
        let z = g.var(Var::Z);
        let a = g.un(UnOp::Abs, z);
        let w = g.constant(w);
        let f = g.sub(a, w);
        cell_0(&g, f).1
    };
    assert!(wall(0.6));
    assert!(!wall(0.3));
}

#[test]
fn split_quads_breaks_roundoff_ties_by_grid_coords_and_keeps_creases() {
    let ctx = Ctx::new(0.0, vec![]);
    let run = |g: &Graph, f: ink_core::graph::NodeId, verts: &[V3], cells: &[[i64; 3]], quads: Vec<[u32; 4]>| {
        let tape = Arc::new(Tape::new(g, f));
        let local = Local::new(&tape, &ctx, [0.0; 3], 8.0);
        split_quads(&local, verts, cells, quads, 1.0, &|| false).unwrap()
    };

    // Flat face z = 0.5 perturbed by +/- 1e-14 * x: the two diagonal midpoints
    // differ only by floating-point roundoff, so both perturbations must pick
    // the same diagonal (ordered grid coords: ([0,0,0],[1,1,0]) < ([0,1,0],[1,0,0])),
    // regardless of cyclic vertex order.
    let flat_verts = [[0.1, 0.1, 0.5], [0.9, 0.2, 0.5], [0.9, 0.9, 0.5], [0.2, 0.9, 0.5]];
    let flat_cells = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]];
    let tilted = |eps: f64| {
        let mut g = Graph::new();
        let (x, z) = (g.var(Var::X), g.var(Var::Z));
        let (half, k) = (g.constant(0.5), g.constant(eps));
        let dz = g.sub(z, half);
        let tilt = g.mul(k, x);
        let f = g.add(dz, tilt);
        (g, f)
    };
    let (gp, fp) = tilted(1e-14);
    let (gn, fn_) = tilted(-1e-14);
    let tp = run(&gp, fp, &flat_verts, &flat_cells, vec![[0, 1, 2, 3]]);
    let tn = run(&gn, fn_, &flat_verts, &flat_cells, vec![[0, 1, 2, 3]]);
    assert_eq!(tp, vec![[0, 1, 2], [0, 2, 3]]);
    assert_eq!(tp, tn);
    // Rotating the quad keeps the same geometric diagonal (0-2).
    let rot = run(&gp, fp, &flat_verts, &flat_cells, vec![[1, 2, 3, 0]]);
    assert_eq!(rot, vec![[1, 2, 0], [2, 3, 0]]);

    // Sharp V-crease z = |x - y|: diagonal 1-3 runs along the crease (f = 0)
    // while diagonal 0-2 cuts across (f = 1). Even though 0-2 has smaller
    // ordered grid coordinates, the crease diagonal 1-3 must win.
    let mut g = Graph::new();
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|v| g.var(v));
    let d = g.sub(x, y);
    let a = g.un(UnOp::Abs, d);
    let f = g.sub(z, a);
    let crease_verts = [[0.0, 1.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]];
    let crease_cells = [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]];
    let tc = run(&g, f, &crease_verts, &crease_cells, vec![[0, 1, 2, 3]]);
    assert_eq!(tc, vec![[0, 1, 3], [1, 2, 3]]);
}
