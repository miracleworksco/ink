use super::contour::{EDGES, REFINE_STEPS, face_center, face_key, faces, patches, qef, refine_crossings};
use super::proof::{CROSSED, Proofs, unshown};
use super::report::Ran;
use crate::octree::{Cell, Grid};
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_points};
use ink_core::interval::I;
use ink_core::sign::{drawn, inside};
use ink_core::tape::{Pass, Tape};
use std::collections::HashMap;
use std::sync::Arc;

pub(super) const B: i64 = 8;
pub(super) const N: usize = B as usize + 1;

pub(super) struct Brick {
    pub(super) i0: [i64; 3],
    /// Its corners' values, where a cell of it changes sign (else none).
    pub(super) vals: Vec<f64>,
    /// (cell, mask of the 12 cell edges this patch owns, vertex position,
    /// its crossings' mean inside the cell: where it goes if its
    /// triangles cross others, `untangle`)
    pub(super) patches: Vec<([i64; 3], u16, V3, V3)>,
    /// Its cells not shown to hold nothing the mesh does not draw
    /// (`unshown`): each cell's index in the brick, and its state
    /// (`INSIDE`, `PROVEN`, `FOUND`, `CROSSED`, `SEAM`).
    pub(super) open: Vec<(u16, u8)>,
    /// Its samples on 0 (corners, saddle-face centres, crossings' steps)
    /// nothing near decides (`ink_core::sign::drawn`): drawn as air,
    /// reported (a wall may be drawn where pieces meet flush).
    pub(super) zeros: Vec<V3>,
    /// Beside each open cell the surface crosses, how thick what it may
    /// hide is (Tape::sheet: a seam's reach, else its box, 2 h); 0 for one
    /// it does not.
    pub(super) thin: Vec<f64>,
    /// Its tape (the field simplified over it and half a cell round it),
    /// where it has open cells the surface does not cross (what the search
    /// reads there) or cells left untried.
    pub(super) tape: Option<Arc<Tape>>,
    /// Its crossed cells left untried, and whether their blocks were tried
    /// (`Proofs::untried`).
    pub(super) untried: Vec<(u16, bool)>,
}

#[inline]
pub(super) fn vidx(i: usize, j: usize, k: usize) -> usize {
    ink_core::grid::index([i, j, k], [N; 3])
}

/// Cell `k` (its index in the brick) of the brick at `i0`, in the grid.
pub(super) fn cell_of(i0: [i64; 3], k: u16) -> [i64; 3] {
    let (b, k) = (B as u16, k);
    [i0[0] + (k % b) as i64, i0[1] + (k / b % b) as i64, i0[2] + (k / (b * b)) as i64]
}

/// A brick's cells: `tape` simplified over it (and half a cell round it),
/// `whole` the field anywhere; with `check` (counting the instructions its
/// proofs run), those to search too, and its samples on 0 nothing near
/// decides. Where `meshed` is given, only the cells it takes (by their
/// index in the grid) are meshed and may be open: the others lie in a
/// region drawn whole (`Flat`), which shows them drawn. None, unsampled,
/// where its samples would take those `sampled` past `most`.
#[allow(clippy::too_many_arguments)]
pub(super) fn brick(
    tape: &Arc<Tape>,
    whole: &Tape,
    ctx: &Ctx,
    grid: &Grid,
    c: Cell,
    meshed: Option<&(dyn Fn([i64; 3]) -> bool + Sync)>,
    (sampled, most): (&std::sync::atomic::AtomicUsize, usize),
    ran: Option<&Ran>,
) -> Option<Brick> {
    let check = ran.is_some();
    // Its samples and crossings lie in its own box: they run `tape`
    // simplified over that, the branches only the half cell round it needs
    // dropped (its proofs and samples on 0, which read there, run `tape`).
    // Bounding and simplifying it is counted as a proof's (`focus`).
    let own = grid.cell_box(c);
    let near = {
        let mut ch = vec![ink_core::tape::BOTH; tape.len()];
        eval_interval(tape, ctx, own, Some(&mut ch));
        let near = if ch.iter().any(|&x| x != ink_core::tape::BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
        let simplified = if Arc::ptr_eq(&near, tape) { 0 } else { tape.simplify_work(&near) };
        _ = ran.map(|r| r.boxes.fetch_add(tape.work(Pass::Box) + simplified, std::sync::atomic::Ordering::Relaxed));
        near
    };
    let count = |points: usize, grads: usize, t: &Tape| {
        _ = ran.map(|r| r.points.fetch_add(points * t.work(Pass::Point) + grads * t.work(Pass::Gradient), std::sync::atomic::Ordering::Relaxed))
    };
    // What deciding samples on 0 runs (`drawn`), for each: its gradient and
    // two samples beside it, and at least its cube's bounds, solid's and
    // air's (more where they split it).
    let zeros_ran = |vals: &[f64]| {
        let z = vals.iter().filter(|v| **v == 0.0).count();
        count(2 * z, z, tape);
        _ = ran.map(|r| r.boxes.fetch_add(4 * z * tape.work(Pass::Box), std::sync::atomic::Ordering::Relaxed));
    };
    let i0 = c.i;
    let mut pts = Vec::with_capacity(N * N * N);
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                pts.push(grid.pos([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64]));
            }
        }
    }
    // The cells meshed (`meshed`; all where none is given), and the
    // corners sampled: theirs (the others NaN).
    let meshed_given = meshed.is_some();
    let cells_meshed: Vec<bool> =
        (0..(B * B * B) as usize).map(|q| meshed.is_none_or(|m| m([i0[0] + q as i64 % B, i0[1] + q as i64 / B % B, i0[2] + q as i64 / (B * B)]))).collect();
    let meshed = |c: [usize; 3]| cells_meshed[c[0] + B as usize * (c[1] + B as usize * c[2])];
    let mut corner = vec![false; N * N * N];
    for (q, _) in cells_meshed.iter().enumerate().filter(|m| *m.1) {
        let (i, j, k) = (q % B as usize, q / B as usize % B as usize, q / (B * B) as usize);
        for c in 0..8 {
            corner[vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))] = true;
        }
    }
    let at: Vec<usize> = (0..N * N * N).filter(|&q| corner[q]).collect();
    if at.is_empty() || sampled.fetch_add(at.len(), std::sync::atomic::Ordering::Relaxed) + at.len() > most {
        return None;
    }
    let taken: Vec<V3> = at.iter().map(|&q| pts[q]).collect();
    let mut got = Vec::new();
    eval_points(&near, ctx, &taken, &mut got);
    count(taken.len(), 0, &near);
    // Where `tape` holds: the brick and half a cell round it.
    let held: [I; 3] = std::array::from_fn(|a| I::new(pts[0][a] - 0.5 * grid.h, pts[pts.len() - 1][a] + 0.5 * grid.h));
    let within = Some((held, whole));
    // Samples on 0 as drawn (`drawn`): those nothing near decides are
    // reported.
    zeros_ran(&got);
    let mut zeros = drawn(tape, ctx, &taken, &mut got, within);
    if !check {
        zeros.clear();
    }
    let mut vals = vec![f64::NAN; N * N * N];
    for (&q, v) in at.iter().zip(got) {
        vals[q] = v;
    }
    // Pass 1: sign-changing cells, and the centres of their saddle faces.
    let mut mixed: Vec<([usize; 3], [bool; 8])> = Vec::new();
    let mut face_keys: HashMap<([i64; 3], usize), usize> = HashMap::new();
    let mut centers: Vec<V3> = Vec::new();
    for k in 0..B as usize {
        for j in 0..B as usize {
            for i in 0..B as usize {
                if !meshed([i, j, k]) {
                    continue;
                }
                let corner = |c: usize| vals[vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1))];
                let ins: [bool; 8] = std::array::from_fn(|c| inside(corner(c)));
                if ins.iter().all(|&v| v) || ins.iter().all(|&v| !v) {
                    continue;
                }
                for f in faces() {
                    let r = f.ring.map(|c| ins[c]);
                    if r[0] == r[2] && r[1] == r[3] && r[0] != r[1] {
                        let key = face_key([i0[0] + i as i64, i0[1] + j as i64, i0[2] + k as i64], f);
                        face_keys.entry(key).or_insert_with(|| {
                            centers.push(face_center(grid, key));
                            centers.len() - 1
                        });
                    }
                }
                mixed.push(([i, j, k], ins));
            }
        }
    }
    let (mut open, mut thin, mut searched, mut untried) = (Vec::new(), Vec::new(), None, Vec::new());
    if let Some(ran) = ran {
        let mut crossed = vec![false; (B * B * B) as usize];
        for (c, _) in &mixed {
            crossed[c[0] + B as usize * (c[1] + B as usize * c[2])] = true;
        }
        // A sample on 0 nothing near decides is reported 2 h thick.
        let mut p = Proofs { settled: !zeros.is_empty(), meshed: meshed_given.then_some(&cells_meshed[..]), ..Default::default() };
        unshown(tape, ctx, grid, &pts, &vals, &crossed, [0; 3], B as usize, false, false, &mut p);
        ran.tapes.fetch_add(p.tapes, std::sync::atomic::Ordering::Relaxed);
        ran.boxes.fetch_add(p.read, std::sync::atomic::Ordering::Relaxed);
        p.open.sort_unstable_by_key(|o| o.0);
        (open, thin) = p.open.into_iter().map(|(k, state, t)| ((k, state), t)).unzip();
        untried = p.untried;
        searched = (open.iter().any(|o: &(u16, u8)| o.1 & CROSSED == 0) || !untried.is_empty()).then(|| tape.clone());
    }
    if mixed.is_empty() {
        thin.resize(open.len(), 0.0);
        return (!open.is_empty() || !zeros.is_empty()).then(|| Brick { i0, vals: Vec::new(), patches: Vec::new(), open, zeros, thin, tape: searched, untried });
    }
    let mut center_vals = Vec::new();
    eval_points(&near, ctx, &centers, &mut center_vals);
    count(centers.len(), 0, &near);
    zeros_ran(&center_vals);
    let open_centres = drawn(tape, ctx, &centers, &mut center_vals, within);
    if check {
        zeros.extend(open_centres);
    }

    // Pass 2: patches and their edge crossings.
    let mut owners: Vec<([i64; 3], u16, usize, usize)> = Vec::new(); // cell, mask, crossing start, count
    let mut brackets: Vec<(V3, f64, V3, f64)> = Vec::new();
    for ([i, j, k], ins) in &mixed {
        let cell = [i0[0] + *i as i64, i0[1] + *j as i64, i0[2] + *k as i64];
        let corner = |c: usize| vidx(i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
        let fc = |fi: usize| inside(center_vals[face_keys[&face_key(cell, &faces()[fi])]]);
        for mask in patches(ins, fc) {
            let start = brackets.len();
            for (e, &(a, b)) in EDGES.iter().enumerate() {
                if mask & (1 << e) == 0 {
                    continue;
                }
                let (a, b) = (corner(a), corner(b));
                let (pa, va, pb, vb) = (pts[a], vals[a], pts[b], vals[b]);
                brackets.push(if inside(va) { (pa, va, pb, vb) } else { (pb, vb, pa, va) });
            }
            owners.push((cell, mask, start, brackets.len() - start));
        }
    }
    let (cross, steps) = refine_crossings(&near, Some((own, whole)), ctx, &mut brackets);
    if check {
        zeros.extend(steps);
    }
    let grads = eval_grads(&near, ctx, &cross);
    count(brackets.len() * REFINE_STEPS, cross.len(), &near);
    let placed: Vec<(V3, V3)> = owners.iter().map(|&(cell, _, s, n)| qef(&cross[s..s + n], &grads[s..s + n], grid.pos(cell), grid.h)).collect();
    let patches = owners.iter().zip(placed).map(|(&(cell, mask, _, _), (v, q))| (cell, mask, v, q)).collect();
    thin.resize(open.len(), 0.0);
    Some(Brick { i0, vals, patches, open, zeros, thin, tape: searched, untried })
}
