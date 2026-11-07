//! Geometric queries with proof-carrying answers: bounds, integrals (volume,
//! mass), emptiness proofs, and planar slices.

use crate::octree::{Cell, Grid, Local, descend, root_for, try_root_for};
use crate::scope::{self, At};
use crate::search;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval, eval_interval_in, eval_interval_scope, eval_point, eval_points, span_mid};
use ink_core::interval::{I, add_o, mul_o, sub_o};
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

pub fn max_extent(lo: V3, hi: V3) -> f64 {
    (0..3).map(|a| hi[a] - lo[a]).fold(0.0, f64::max)
}

/// A cell's volume, as clipped (an estimate's; bounds use box_volume_out).
fn cell_volume(grid: &Grid, c: Cell) -> f64 {
    grid.cell_sides(c).iter().product()
}

/// `k^3` sample points at the centres of the sub-cells of a cell (as
/// clipped), `k` along each side.
pub fn samples(grid: &Grid, c: Cell, k: usize) -> Vec<V3> {
    let (lo, _) = grid.cell_bounds(c);
    let s = grid.cell_sides(c).map(|x| x / k as f64);
    let mut out = Vec::with_capacity(k * k * k);
    for a in 0..k {
        for b in 0..k {
            for d in 0..k {
                out.push([lo[0] + s[0] * (d as f64 + 0.5), lo[1] + s[1] * (b as f64 + 0.5), lo[2] + s[2] * (a as f64 + 0.5)]);
            }
        }
    }
    out
}

/// Fraction of a box of sides `s` centred on a sample that lies inside,
/// the field taken as linear there: the slab fraction along the gradient.
pub fn occupancy(f: &[f64; 4], s: V3) -> f64 {
    let w = s[0] * f[1].abs() + s[1] * f[2].abs() + s[2] * f[3].abs(); // field range across the box
    if !(w > 0.0 && w.is_finite() && f[0].is_finite()) {
        return if f[0] < 0.0 { 1.0 } else { 0.0 };
    }
    (0.5 - f[0] / w).clamp(0.0, 1.0)
}

/// The tape simplified for the region `lo..hi`, the field there and only
/// there: one interval pass settles which min/max branches can matter in
/// it, so a union of thousands of parts costs only the parts nearby in
/// everything that follows. Searches with it examine boxes inside the
/// region (their cells clipped to it: octree::try_root_for).
pub fn focus(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Arc<Tape> {
    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let mut ch = vec![BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    if ch.iter().any(|&c| c != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() }
}

/// Proven outer bounding box of the region where the field is < 0, found in
/// two passes (coarse, then refined). Returns `(lo, hi, cell_size)`; `None`
/// means the field is proven empty.
pub fn bounds(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Option<(V3, V3, f64)> {
    crate::cache::memo_region("bounds", &[tape], ctx, (lo, hi), |_| {}, |c| bounds_passes(tape, c, lo, hi, &[64.0, 256.0]))
}

/// Where a field can have material (be negative), over all of space. Its
/// zero set without interior (a surface of a field never negative) is not
/// material: a field 0 everywhere is empty.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Extent {
    /// Interval bounds prove the field positive everywhere.
    Empty,
    /// The field is proven positive outside this box.
    Bounded(V3, V3),
    /// Not proven bounded: material may reach arbitrarily far (a half-space,
    /// an unclipped lattice), or interval bounds cannot show that it does not.
    Unbounded,
}

/// Cells one search of `extent` may evaluate before it gives up and
/// reports `Unbounded`.
const EXTENT_BUDGET: usize = 1 << 13;

/// 2^40: farther than any model reaches, in window sizes (`extent_sectors`).
const TB: f64 = 1_099_511_627_776.0;

/// Where the field can be negative, over all of space, proven by covering all of space with pieces on which
/// interval bounds (`eval_interval_linear`, sound with infinite ends, whose
/// reals are finite) show the field positive, except pieces that join the
/// box. Two covers, as each settles shapes the other cannot: far sectors
/// (`extent_sectors`: narrow directions out to infinity, so a part rotated
/// about a slanted axis is settled), then, if that fails, boxes in
/// compactified coordinates (`extent_boxes`: slabs out to infinity, so a
/// twist's angle narrows with its axis). So the answer cannot miss material
/// however far it lies, and a field neither cover can bound (a half-space,
/// an unclipped lattice, a twist turning without limit) is reported as
/// such. A shape inside the window gets the window as its box.
///
/// The covers search the field's members, not just the whole
/// (`extent_within`): a union's material is its members' material, and an
/// intersection's lies in any member's. So a union of parts that need
/// different covers is bounded, and a search runs only on small tapes.
/// Every search starts with a small budget, doubled up to `EXTENT_BUDGET`
/// while undecided, so no member spends more than its cheapest proof needs
/// (an intersection whose members are unbounded, as an extruded profile,
/// is settled whole at once); the doubling at most doubles the cost of a
/// search that fails, and only the members that fail are searched again.
/// Members are searched cheapest first (`members`), so an intersection's
/// box is its cheapest member's proven bounded at the least budget.
pub fn extent(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3)) -> Extent {
    crate::cache::memo_region("extent", &[tape], ctx, window, |_| {}, |c| extent_now(tape, c, window))
}

fn extent_now(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3)) -> Extent {
    let mut found = HashMap::new();
    let mut budget = EXTENT_BUDGET >> 10;
    loop {
        if let Some(e) = extent_within(tape, tape.len() - 1, ctx, window, budget, &mut found) {
            return e;
        }
        if budget >= EXTENT_BUDGET {
            return Extent::Unbounded;
        }
        budget *= 2;
    }
}

/// The extent of instruction `i` of `tape` within `budget` evaluations per
/// search (None: undecided at that budget). `found` holds each instruction's
/// extent, or the budget it was undecided at.
fn extent_within(tape: &Arc<Tape>, i: usize, ctx: &Ctx, window: (V3, V3), budget: usize, found: &mut HashMap<usize, Result<Extent, usize>>) -> Option<Extent> {
    use ink_core::graph::BinOp;
    use ink_core::tape::Op;
    let i = value_of(tape, i);
    match found.get(&i) {
        Some(Ok(e)) => return Some(*e),
        Some(Err(b)) if *b >= budget => return None,
        _ => {}
    }
    let first = !found.contains_key(&i);
    let whole = |i: usize| if i == tape.len() - 1 { tape.clone() } else { Arc::new(tape.cone(i)) };
    let e = match tape.ssa[i].op {
        // Material in any member: the union of their boxes.
        Op::Bin(BinOp::Min) => {
            let mut u = Extent::Empty;
            for m in members(tape, i, BinOp::Min) {
                u = match (u, extent_within(tape, m, ctx, window, budget, found)) {
                    (_, None) => {
                        found.insert(i, Err(budget));
                        return None;
                    }
                    (Extent::Bounded(l0, h0), Some(Extent::Bounded(l1, h1))) => {
                        let (lo, hi) = ink_core::v3::box_union((l0, h0), (l1, h1));
                        Extent::Bounded(lo, hi)
                    }
                    (Extent::Empty, Some(e)) => e,
                    (e, Some(_)) => e,
                };
            }
            Some(u)
        }
        // Material in every member: inside any member's box, or found whole
        // (slabs and profiles are unbounded alone); or none anywhere, by the
        // first-order bound of the max over all of space, which the covers'
        // plain bounds cannot see: members that never meet, as two parallel
        // slabs, whose exactly opposite slopes cancel however far they reach.
        // That bound does not depend on the budget: it is tried once, where
        // the covers fail at the first (smallest) budget, so an extrusion,
        // which they settle there, never pays for it.
        Op::Bin(BinOp::Max) => members(tape, i, BinOp::Max).into_iter().find_map(|m| extent_within(tape, m, ctx, window, budget, found)).or_else(|| {
            let t = whole(i);
            let never = || crate::centred::centred_at_least(&[(&t, 0.0)], ctx, [I::ENTIRE; 3], 0.0);
            extent_search(&t, ctx, window, budget).or_else(|| (first && never()).then_some(Extent::Empty))
        }),
        _ => extent_search(&whole(i), ctx, window, budget),
    };
    found.insert(i, e.ok_or(budget));
    e
}

/// The instruction whose value instruction `i` passes on (through tags,
/// Lipschitz promises and derivatives' marks).
fn value_of(tape: &Tape, mut i: usize) -> usize {
    use ink_core::tape::Op;
    while let Op::Tag(_) | Op::Lipschitz | Op::Deriv(..) | Op::Real = tape.ssa[i].op {
        i = tape.ssa[i].a as usize;
    }
    i
}

/// The members of the `op` (min or max) at instruction `i`: its arguments,
/// with nested `op`s flattened, cheapest first: by the instructions each
/// reads (its cone: what evaluating it costs), then in tape order, which
/// is the hash order of commutative arguments (Graph::bin) and put a
/// pair's long belt before its small part. Each instruction is taken once
/// (min and max are idempotent), so a shared subgraph does not unfold.
fn members(tape: &Tape, i: usize, op: ink_core::graph::BinOp) -> Vec<usize> {
    use ink_core::tape::Op;
    // Each instruction marked with the last walk that reached it: the
    // flattening's (i), then each member's own (m).
    let (mut out, mut stack, mut mark) = (Vec::new(), vec![i], vec![usize::MAX; tape.ssa.len()]);
    while let Some(j) = stack.pop() {
        let j = value_of(tape, j);
        if std::mem::replace(&mut mark[j], i) == i {
            continue;
        }
        match tape.ssa[j].op {
            Op::Bin(o) if o == op => stack.extend([tape.ssa[j].b as usize, tape.ssa[j].a as usize]),
            _ => out.push(j),
        }
    }
    out.sort_by_cached_key(|&m| {
        let (mut n, mut stack) = (0, vec![m]);
        mark[m] = m;
        while let Some(j) = stack.pop() {
            n += 1;
            let (args, k) = tape.ssa[j].args();
            for &a in &args[..k] {
                if std::mem::replace(&mut mark[a as usize], m) != m {
                    stack.push(a as usize);
                }
            }
        }
        (n, m)
    });
    out
}

/// One field's extent by the two covers, each within `budget` evaluations.
fn extent_search(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3), budget: usize) -> Option<Extent> {
    match extent_sectors(tape, ctx, window, budget) {
        Extent::Unbounded => match extent_boxes(tape, ctx, window, budget) {
            Extent::Unbounded => None,
            e => Some(e),
        },
        e => Some(e),
    }
}

/// A window's centre c and half-sides r, as the covers take them: their
/// pieces are c + r (s q), whatever window they came from, the central box
/// the one they start from (s = 1), so together they cover all of space.
/// (r must not be 0, or the sectors would collapse onto c; both are halved
/// before the sum, which overflows for a window's ends near 1e308.)
fn centre_radius(window: (V3, V3)) -> (V3, V3) {
    (std::array::from_fn(|k| 0.5 * window.0[k] + 0.5 * window.1[k]), std::array::from_fn(|k| (0.5 * window.1[k] - 0.5 * window.0[k]).max(f64::MIN_POSITIVE)))
}

/// `extent` by far sectors: around the window, six cones (one per face,
/// scaled out from its centre), each split into sectors: a range of
/// directions and a range of distance, the farthest reaching infinity. Over
/// a sector a linear term's distance factor is taken once, so a rotated
/// coordinate stays one-signed on a narrow sector however far it reaches.
/// Undecided sectors split by direction, then by distance (doubling); a
/// finite one small for its distance joins the box, a far one never does.
fn extent_sectors(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3), budget: usize) -> Extent {
    use ink_core::eval::{Region, eval_interval_linear};
    let (c, r) = centre_radius(window);
    /// A sector: the face (axis, side), its direction ranges along the
    /// other two axes (in [-1, 1]), and its range of distance (scale of the
    /// window, from 1; hi infinite for the far end).
    #[derive(Clone)]
    struct Sector {
        axis: usize,
        side: f64,
        v: [(f64, f64); 2],
        s: (f64, f64),
        tape: Arc<Tape>,
    }
    let others = |axis: usize| -> [usize; 2] { [(axis + 1) % 3, (axis + 2) % 3] };
    let region = |x: &Sector| {
        let mut q = [I::point(0.0); 3];
        q[x.axis] = I::point(x.side);
        for (j, &k) in others(x.axis).iter().enumerate() {
            q[k] = I::new(x.v[j].0, x.v[j].1);
        }
        Region::Sector { c, r, q, s: I::new(x.s.0, x.s.1) }
    };
    let (mut blo, mut bhi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut join = |b: [I; 3]| {
        for k in 0..3 {
            blo[k] = blo[k].min(b[k].lo);
            bhi[k] = bhi[k].max(b[k].hi);
        }
    };
    // The window: material there joins the box as it is.
    // Rounded outward, so it holds the real box c +- r.
    let wbox: [I; 3] = std::array::from_fn(|k| I::new((c[k] - r[k]).next_down(), (c[k] + r[k]).next_up()));
    if eval_interval_linear(tape, ctx, &Region::Box(wbox), None).lo < 0.0 {
        join(wbox);
    }
    let mut stack: Vec<Sector> =
        (0..6).map(|f| Sector { axis: f / 2, side: if f % 2 == 0 { -1.0 } else { 1.0 }, v: [(-1.0, 1.0); 2], s: (1.0, f64::INFINITY), tape: tape.clone() }).collect();
    let mut evals = 0;
    while let Some(x) = stack.pop() {
        // Doubled past the largest float: no real point is that far.
        if x.s.0.is_infinite() {
            continue;
        }
        let reg = region(&x);
        let mut ch = vec![BOTH; x.tape.len()];
        let v = eval_interval_linear(&x.tape, ctx, &reg, Some(&mut ch));
        evals += 1;
        if v.lo >= 0.0 {
            continue;
        }
        let far = x.s.1.is_infinite();
        let width = (x.v[0].1 - x.v[0].0).max(x.v[1].1 - x.v[1].0);
        // A finite sector small for its distance (directions an eighth of a
        // face, or half of one twice the window's size out: the box is
        // tightened later, and a big shape far out stays cheap) is material,
        // or too fine to settle: it joins the box.
        let fine = if x.s.0 >= 2.0 { 0.5 } else { 0.125 };
        if !far && width <= fine && x.s.1 <= 2.0 * x.s.0 {
            join(reg.bounds());
            continue;
        }
        if evals >= budget {
            return Extent::Unbounded;
        }
        let sub = if ch.iter().any(|&k| k != BOTH) { Arc::new(x.tape.simplify(&ch)) } else { x.tape.clone() };
        let mut kid = |f: &dyn Fn(&mut Sector)| {
            let mut y = Sector { tape: sub.clone(), ..x.clone() };
            f(&mut y);
            stack.push(y);
        };
        // A finite sector splits its directions. A far one splits them only
        // if they do not settle its far end out where no model reaches (2^40
        // window sizes: there only narrower directions can help); else its
        // distance: the far end doubles away (material in between stays in
        // the finite pieces).
        let narrower = !far
            || width > 1.0 / 64.0 && {
                evals += 1;
                let mut y = region(&x);
                if let Region::Sector { s, .. } = &mut y {
                    *s = I::new(x.s.0.max(TB), f64::INFINITY);
                }
                eval_interval_linear(&sub, ctx, &y, None).lo < 0.0
            };
        if narrower {
            let j = if x.v[0].1 - x.v[0].0 >= x.v[1].1 - x.v[1].0 { 0 } else { 1 };
            let m = 0.5 * (x.v[j].0 + x.v[j].1);
            kid(&|y| y.v[j].1 = m);
            kid(&|y| y.v[j].0 = m);
        } else {
            let m = 2.0 * x.s.0;
            kid(&|y| y.s.1 = m);
            kid(&|y| y.s.0 = m);
        }
    }
    if blo[0] > bhi[0] { Extent::Empty } else { Extent::Bounded(blo, bhi) }
}

/// The compactified coordinate of `extent_boxes` on one axis: u in (-1, 1)
/// maps to x = c + r u / (1 - |u|), so u = +-1/2 is the window's faces c +- r,
/// and u = +-1 is infinity.
fn from_u(u: f64, c: f64, r: f64) -> f64 {
    if u >= 1.0 {
        f64::INFINITY
    } else if u <= -1.0 {
        f64::NEG_INFINITY
    } else if u == 0.5 {
        c + r
    } else if u == -0.5 {
        c - r
    } else {
        c + r * u / (1.0 - u.abs())
    }
}

/// `extent` by boxes: an octree over compactified coordinates centred on
/// the window (its root is all of space, the cell at the middle the window,
/// cells at the rim reach infinity). Only cells reaching infinity are
/// refined, widest axis first (so a slab's finite axes narrow too); finite
/// cells not proven empty join the box.
fn extent_boxes(tape: &Arc<Tape>, ctx: &Ctx, window: (V3, V3), budget: usize) -> Extent {
    use ink_core::eval::{Region, eval_interval_linear};
    let (c, r) = centre_radius(window);
    struct C {
        lo: V3,
        hi: V3,
        tape: Arc<Tape>,
    }
    let mut stack = vec![C { lo: [-1.0; 3], hi: [1.0; 3], tape: tape.clone() }];
    let (mut blo, mut bhi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    let mut evals = 0;
    while let Some(cell) = stack.pop() {
        let b: [I; 3] = std::array::from_fn(|k| I::new(from_u(cell.lo[k], c[k], r[k]), from_u(cell.hi[k], c[k], r[k])));
        let mut ch = vec![BOTH; cell.tape.len()];
        let v = eval_interval_linear(&cell.tape, ctx, &Region::Box(b), Some(&mut ch));
        evals += 1;
        if v.lo >= 0.0 {
            continue;
        }
        // A cell at u = +-1 alone holds no real point.
        if (0..3).any(|k| cell.lo[k] >= 1.0 || cell.hi[k] <= -1.0) {
            continue;
        }
        let infinite: Vec<usize> = (0..3).filter(|&k| cell.lo[k] <= -1.0 || cell.hi[k] >= 1.0).collect();
        if infinite.is_empty() {
            for k in 0..3 {
                blo[k] = blo[k].min(b[k].lo);
                bhi[k] = bhi[k].max(b[k].hi);
            }
            continue;
        }
        if evals >= budget {
            return Extent::Unbounded;
        }
        let sub = if ch.iter().any(|&x| x != BOTH) { Arc::new(cell.tape.simplify(&ch)) } else { cell.tape.clone() };
        // Halve the widest axis (one that reaches infinity on a tie): finite
        // axes split too, so a twist's angle over a slab narrows.
        let width = |k: usize| (cell.hi[k] - cell.lo[k], infinite.contains(&k));
        let k = (0..3).max_by(|&i, &j| width(i).partial_cmp(&width(j)).unwrap()).unwrap();
        let m = 0.5 * (cell.lo[k] + cell.hi[k]);
        let (mut h1, mut l2) = (cell.hi, cell.lo);
        (h1[k], l2[k]) = (m, m);
        stack.push(C { lo: cell.lo, hi: h1, tape: sub.clone() });
        stack.push(C { lo: l2, hi: cell.hi, tape: sub });
    }
    if blo[0] > bhi[0] { Extent::Empty } else { Extent::Bounded(blo, bhi) }
}

/// Two quick passes (1/32, then 1/64 of the narrowed region): enough to focus a proof's search.
pub fn bounds_quick(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> Option<(V3, V3, f64)> {
    crate::cache::memo_region("bounds-quick", &[tape], ctx, (lo, hi), |_| {}, |c| bounds_passes(tape, c, lo, hi, &[32.0, 64.0]))
}

/// How far the material of a shape (where the field is < 0) inside `lo..hi`
/// lies from the z axis, at every value of ctx's scope (t a step, ranged
/// params): a lower bound (a point some material reaches: the farthest
/// corner of a box wholly material, or a box's middle proven material) and
/// an upper, within `tol` of each other where `budget` boxes allow. The one
/// search (crate::search), farthest first (a box's bound is its farthest
/// corner from the axis, negated: none of its points is farther; the upper
/// bound is the farthest of the boxes left and those kept), halving
/// across x and y, and across z where `split_z` (a field that reads z) and
/// z is the box's longest side; a box proven empty is settled, one wholly
/// material, within `tol` across (from the axis or the lower bound) or no
/// longer halved keeps its farthest corner in the upper bound. None:
/// proven empty.
pub fn radial_extent(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, split_z: bool, tol: f64, budget: usize) -> Option<(f64, f64)> {
    use std::cell::Cell;
    let far = |b: &[I; 3]| radius(b[0].lo.abs().max(b[0].hi.abs()), b[1].lo.abs().max(b[1].hi.abs())).hi;
    // The farthest corner itself, rounded down: a point of the box.
    let corner = |b: &[I; 3]| {
        let pick = |x: I| if x.lo.abs() > x.hi.abs() { x.lo } else { x.hi };
        radius(pick(b[0]), pick(b[1])).lo
    };
    // `kept`: the farthest of the boxes that may hold material and are
    // halved no further; `seen`: whether any may.
    let (lower, kept, seen) = (Cell::new(0.0f64), Cell::new(0.0f64), Cell::new(false));
    let test = |n: &search::Node<f64>| -> search::Test<f64, ()> {
        let v = eval_interval(tape, &n.ctx, n.b, None);
        if v.lo > 0.0 {
            return search::Test::Settled;
        }
        seen.set(true);
        let gap = |x: I| {
            if x.lo > 0.0 {
                x.lo
            } else if x.hi < 0.0 {
                -x.hi
            } else {
                0.0
            }
        };
        if v.hi < 0.0 || far(&n.b) - radius(gap(n.b[0]), gap(n.b[1])).lo.max(lower.get()) <= tol {
            if v.hi < 0.0 {
                lower.set(lower.get().max(corner(&n.b)));
            }
            kept.set(kept.get().max(far(&n.b)));
            return search::Test::Settled;
        }
        let mid = n.b.map(|x| x.mid());
        if eval_interval(tape, &n.ctx, mid.map(I::point), None).hi < 0.0 {
            lower.set(lower.get().max(radius(mid[0], mid[1]).lo));
        }
        let w = n.b.map(|x| x.hi - x.lo);
        let space = [w[0], w[1], if split_z && w[2] > w[0].max(w[1]) { w[2] } else { 0.0 }];
        search::Test::Open { bound: -far(&n.b), smear: search::Smear { space, ..Default::default() }, state: far(&n.b) }
    };
    let root = search::Node { b: std::array::from_fn(|k| I::new(lo[k], hi[k])), ctx: Arc::new(ctx.clone()), state: 0.0 };
    let first = test(&root);
    let mut s = search::Search::start(vec![(root, first)], search::Precision { finest: 0.0, splits: 0, maps: true });
    s.serially(budget, 8, test);
    let upper = s.left().map(|n| n.state).fold(kept.get(), f64::max);
    (seen.get() || s.left().next().is_some()).then(|| (lower.get(), upper.max(lower.get())))
}

/// sqrt(x^2 + y^2) bounded: lo <= the real radius <= hi. Scaled by a
/// power of two near the larger term so the squares neither underflow nor
/// overflow (a radius of 1e-200 is not 0), and each operation (scale,
/// square, sum, sqrt, rescale) rounded outward: one ulp after a float norm
/// bounds nothing.
fn radius(x: f64, y: f64) -> I {
    use ink_core::graph::UnOp;
    use ink_core::interval::un;
    let m = x.abs().max(y.abs());
    if m == 0.0 || !m.is_finite() {
        return I::point(m);
    }
    let e = (m.log2().floor() as i32).clamp(-1000, 1000);
    let sq = |v: f64| un(UnOp::Square, mul_o(I::point(v), I::point(2f64.powi(-e))));
    mul_o(un(UnOp::Sqrt, add_o(sq(x), sq(y))), I::point(2f64.powi(e)))
}

/// A box narrowed to under 3 cells of the last pass (a small part far
/// inside a big region) gets more passes, while they narrow it: each such
/// pass shrinks it at least 85 times, and a box of floats cannot shrink for
/// ever (a point-like set reaches the precision of floats in a few).
fn bounds_passes(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, passes: &[f64]) -> Option<(V3, V3, f64)> {
    let (mut lo, mut hi) = (lo, hi);
    let mut h = 0.0;
    let last = passes[passes.len() - 1];
    let extra = passes.iter().copied().chain(std::iter::repeat(last));
    let mut n = 0;
    let mut before = None;
    for divisions in extra {
        n += 1;
        if n > passes.len() && (max_extent(lo, hi) >= 3.0 * h || before == Some((lo, hi))) {
            break;
        }
        before = Some((lo, hi));
        h = max_extent(lo, hi) / divisions;
        // No root covers the region: no bounds tighter than it.
        let Some((grid, root)) = try_root_for(lo, hi, h, 1) else { return Some((lo, hi, h)) };
        let ((nlo, nhi), _) = bounds_cells(tape, ctx, &grid, root);
        if nlo[0] > nhi[0] {
            return None;
        }
        for k in 0..3 {
            lo[k] = lo[k].max(nlo[k]);
            hi[k] = hi[k].min(nhi[k]);
        }
    }
    Some((lo, hi, h))
}

/// Cells already enclosed by the box found so far cannot enlarge it.
/// For a cell crossing that box, evaluate only its outside slabs first:
/// its inside material must not force subdivision of an entire face just
/// to exclude the empty strip beyond the current bound. The result is
/// an enclosure of the material, as in a full descent.
/// With the cells it evaluated: a work count a test can bound, whatever the timing.
fn bounds_cells(tape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, root: Cell) -> ((V3, V3), usize) {
    let (mut lo, mut hi) = ink_core::v3::EMPTY_BOX;
    let mut pending = vec![(tape.clone(), root)];
    let mut visited = 0;
    while let Some((t, c)) = pending.pop() {
        if !grid.meets(c) {
            continue;
        }
        let bx = grid.cell_box(c);
        if lo[0] <= hi[0] {
            let mut empty = true;
            for k in 0..3 {
                for upper in [false, true] {
                    if if upper { bx[k].hi <= hi[k] } else { bx[k].lo >= lo[k] } {
                        continue;
                    }
                    let mut slab = bx;
                    slab[k] = if upper { I::new(bx[k].lo.max(hi[k]), bx[k].hi) } else { I::new(bx[k].lo, bx[k].hi.min(lo[k])) };
                    visited += 1;
                    if eval_interval(&t, ctx, slab, None).lo <= 0.0 {
                        empty = false;
                        break;
                    }
                }
                if !empty {
                    break;
                }
            }
            if empty {
                continue;
            }
        }
        visited += 1;
        let mut ch = vec![BOTH; t.len()];
        let r = eval_interval(&t, ctx, bx, Some(&mut ch));
        if r.lo > 0.0 {
            continue;
        }
        if r.hi < 0.0 || c.size == 1 {
            (lo, hi) = ink_core::v3::box_union((lo, hi), grid.cell_bounds(c));
        } else {
            let sub = if ch.iter().any(|&c| c != BOTH) { Arc::new(t.simplify(&ch)) } else { t };
            pending.extend(c.children().into_iter().map(|c| (sub.clone(), c)));
        }
    }
    ((lo, hi), visited)
}

/// The relative spread an estimated integral settles to (`refine`): a
/// mass, a volume, a published mass's density, a claim's integrals at
/// the latest.
pub const ACCURACY: f64 = 0.02;

/// Cells a halving of `refine` may make: one that would make more is not
/// made (surfaces over a large area at fine cells: a few seconds' work).
pub const PASS_CELLS: usize = 1 << 21;

/// Integrals at one cell size (`sums`), and the cells they took, full and
/// surface.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Sums {
    pub sums: Vec<f64>,
    pub cells: [usize; 2],
}

/// Integrals as `refine` left them, at cells of `res`: each with its
/// spread (none before a halving), the size before (`coarse`: it, and how far the first integral
/// moved from there, relative to it), the cells it took (full of
/// material, and the surface's), and whether it settled.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Measured {
    pub sums: Vec<f64>,
    pub spread: Vec<f64>,
    pub res: f64,
    pub coarse: Option<(f64, f64)>,
    pub cells: [usize; 2],
    pub settled: bool,
}

impl Measured {
    /// `at`'s integrals at cells of `h`, before any halving: no spread.
    /// Settled where there are none.
    pub fn first(h: f64, at: impl FnOnce(f64) -> Sums) -> Measured {
        let s = at(h);
        Measured { settled: s.sums.is_empty(), sums: s.sums, spread: Vec::new(), res: h, coarse: None, cells: s.cells }
    }

    /// Halved once (`refine`).
    pub fn halved(&self, at: impl FnOnce(f64) -> Sums) -> Measured {
        let next = at(self.res / 2.0);
        let spread: Vec<f64> = next.sums.iter().zip(&self.sums).map(|(a, b)| 2.0 * (a - b).abs()).collect();
        let (first, change) = (next.sums[0], (next.sums[0] - self.sums[0]).abs());
        let settled = next.cells == [0; 2] || (next.cells[0] > 0 && spread[0] <= ACCURACY * first.abs());
        let coarse = Some((self.res, if change == 0.0 { 0.0 } else { change / first.abs() }));
        Measured { sums: next.sums, spread, res: self.res / 2.0, coarse, cells: next.cells, settled }
    }

    /// Whether `refine` halves it again: unsettled, finite, and (past the
    /// first halving) the next within PASS_CELLS.
    pub fn halves(&self) -> bool {
        let cells = self.cells[0] + self.cells[1];
        !self.settled && self.sums[0].is_finite() && !(self.coarse.is_some() && cells.saturating_mul(4) > PASS_CELLS)
    }

    /// Of `measured` with moments: rho * [1, x, y, z, xx, yy, zz, xy, yz,
    /// zx], then (the second integrand's, where it is 1) the volume.
    pub fn moments(&self) -> [f64; 10] {
        std::array::from_fn(|k| self.sums[k])
    }
    pub fn volume(&self) -> f64 {
        self.sums[10]
    }

    /// The first integral's spread, relative to it (None before a halving).
    pub fn relative_spread(&self) -> Option<f64> {
        self.spread.first().map(|&s| if s == 0.0 { 0.0 } else { s / self.sums[0].abs() })
    }

    /// Why it did not settle (None: it did).
    pub fn unsettled(&self) -> Option<String> {
        let why = match self.cells[0] {
            0 => "no cell is proven full of material (a wall thinner than every cell)".to_string(),
            _ => format!("its spread is over {}% of it", ACCURACY * 100.0),
        };
        (!self.settled).then(|| format!("{why}, and the next halving would take over {PASS_CELLS} cells"))
    }
}

/// Integrals estimated, the one way they are: `at` (integrals at a cell
/// size) at cells of `h`, halved until the first's spread is within
/// ACCURACY of it, at cells some of which are proven full of material. A
/// spread is a sampled estimate of the change left: twice what the last
/// halving moved it. So at least one halving is made, and the answer is an estimate, not a
/// bound. Cells none of which fits inside the material resolve none of
/// its walls, and can read them alike at every size that leaves them so
/// (a 0.25 mm plate 200 mm across is 48% heavy at cells of 2 and of 1
/// mm). A halving takes about four times the cells (they cover a
/// surface): after the first, one that would pass PASS_CELLS is not made,
/// and the answer is unsettled. A claim's parts are halved this way each,
/// and stop sooner where the claim's status survives their spreads
/// (claims::estimate_parts).
pub fn refine(h: f64, mut at: impl FnMut(f64) -> Sums) -> Measured {
    let mut m = Measured::first(h, &mut at);
    while m.halves() {
        m = m.halved(&mut at);
    }
    m
}

/// `sums` refined (`refine`) from cells of a 100th of the box, to
/// ACCURACY; at `res`, once, where given. Cached, as the kernel's other
/// answers are. The cells follow walls, not length (a 1200 mm tube's
/// 0.71 mm walls measure 59% heavy at a 200th of its length).
pub fn measured(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>) -> Measured {
    let tapes: Vec<&Tape> = std::iter::once(&**shape).chain(gs.iter().flatten().copied()).collect();
    crate::cache::memo_region(
        "measured",
        &tapes,
        ctx,
        (lo, hi),
        |fp| {
            fp.word(res.map_or(0, |r| 1 + r.to_bits())).word(gs.len() as u64).word(moments as u64);
            gs.iter().for_each(|g| _ = fp.word(g.is_some() as u64));
        },
        |ctx| {
            let at = |h| sums(shape, gs, moments, ctx, lo, hi, h);
            match res {
                Some(h) => Measured { settled: true, ..Measured::first(h, at) },
                None => refine(max_extent(lo, hi) / 100.0, at),
            }
        },
    )
}

/// Bounds of the integral of `g` over a cell full of material (of every
/// moment and param value of ctx's spans), rounded outward, and whether
/// g's slopes are the same across it to 1e-9 of their size (their range's
/// width: g linear but for rounding, so splits would narrow little): the
/// tighter of its range times the cell's volume and the centred form. With
/// c a point of the cell and a_k(x) the mean of dg/dx_k on the segment from
/// c to x (in the range of that slope over the cell, which covers both
/// sides of a kink), g(x) = g(c) + sum_k a_k(x) (x_k - c_k), so the
/// integral is vol g(c) plus, per axis, at most the slope's range times
/// the integral of x_k - c_k where positive, less where negative: exact
/// for g linear in x, y, z (a centre of mass's moment). A g that jumps in
/// the cell (a select, a floor) has unbounded slopes there and keeps the
/// plain bound. Halved while the slopes vary over it, down to `leaf` grid
/// steps.
fn full_integral(g: &Tape, ctx: &Ctx, grid: &Grid, c: Cell, leaf: i64) -> (I, bool) {
    let b = grid.cell_box(c);
    let side: [I; 3] = b.map(|x| sub_o(I::point(x.hi), I::point(x.lo)));
    let vol = mul_o(mul_o(side[0], side[1]), side[2]);
    let plain = mul_o(vol, eval_interval(g, ctx, b, None));
    let gg = ink_core::eval::eval_interval_grad(g, ctx, b);
    let centre: V3 = b.map(|x| span_mid((x.lo, x.hi)));
    let mut centred = mul_o(vol, eval_interval(g, ctx, centre.map(I::point), None));
    for k in 0..3 {
        // The integral of x_k - c_k where positive and where negative: the
        // cross-section times half the squared distance to each face.
        let area = mul_o(side[(k + 1) % 3], side[(k + 2) % 3]);
        let half = |d: I| mul_o(mul_o(area, mul_o(d, d)), I::point(0.5));
        let up = half(sub_o(I::point(b[k].hi), I::point(centre[k])));
        let down = half(sub_o(I::point(centre[k]), I::point(b[k].lo)));
        let term = sub_o(mul_o(gg.d[k], up), mul_o(gg.d[k], down));
        centred = add_o(centred, term);
    }
    let both = I::new(plain.lo.max(centred.lo), plain.hi.min(centred.hi));
    let linear = |d: &I| d.lo.is_finite() && d.hi.is_finite() && d.hi - d.lo <= 1e-9 * ink_core::interval::mag(*d);
    let curved = !gg.d.iter().all(linear);
    if !curved || c.size <= leaf {
        return (both, !curved);
    }
    let halves = c.children().map(|k| if grid.meets(k) { full_integral(g, ctx, grid, k, leaf).0 } else { I::point(0.0) });
    let sum = halves.iter().fold(I::point(0.0), |s, x| add_o(s, *x));
    (I::new(sum.lo.max(both.lo), sum.hi.min(both.hi)), false)
}

#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Verdict {
    Proven,
    /// A point inside, at time t and the ranged params' values (none when
    /// no param is ranged), with the field's value there.
    Disproven {
        witness: V3,
        t: f64,
        value: f64,
        params: Vec<(u32, f64)>,
    },
    /// Boxes left undecided and their volume, and one of them, (lo, hi):
    /// the one that ended the proof (stuck at the finest cells), else the
    /// first left in its order (none where it had none); why the proof
    /// ended (`ended`); the hunt for a violation after it, none found, and
    /// why it stopped (`hunt`: search::Search::hunt); and the accuracy it holds to
    /// (`attained`: nothing violates it deeper than that, from the boxes
    /// left: no point of an emptiness's solid deeper; a clearance d, no
    /// pair closer than d less it; None: none proven, or the time limit
    /// stopped it).
    Unknown {
        cells: usize,
        volume: f64,
        ended: Ended,
        cell: Option<[V3; 2]>,
        attained: Option<f64>,
        hunt: Option<search::Hunt>,
    },
}

/// Why a search ended undecided.
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Ended {
    /// A box undecided with every side at its finest cells (`cell`), the
    /// first the search reached: the claim cannot be proven at them.
    Finest,
    /// It tested the boxes it may (BOXES, or the budget asked), or would
    /// have tested over four times as many (search::Search::round).
    Limit,
    /// The call's time limit passed (crate::stop): the boxes are those it
    /// left.
    Stopped,
}

/// An emptiness search's answer: its verdict, the cell size it reached,
/// the tolerance it was asked to (mm; 0: exactly), and whether a proof
/// needed it (`exact` false: some cell was shown only to hold no point
/// deeper than it, so the solid may not be empty; `Proof::exact`), with
/// up to LOOSE such cells, the lowest bounds first, where a point of the
/// solid may be (`loose`).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Proof {
    pub verdict: Verdict,
    pub res: f64,
    pub tol: f64,
    pub exact: bool,
    pub loose: Vec<Loose>,
}

/// A cell a proof settled only by its tolerance: its box (lo, hi) and the
/// piece of the scope it was settled over, each ranged param's span
/// (`params`, those derived from them left out) and time's over a step
/// (`t`).
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Loose {
    pub cell: [V3; 2],
    pub params: Vec<(u32, (f64, f64))>,
    pub t: Option<(f64, f64)>,
}

impl Loose {
    fn of(b: [I; 3], ctx: &Ctx) -> Loose {
        let params = ctx.p_spans.iter().filter(|s| !s.derived).map(|s| (s.k, s.span)).collect();
        Loose { cell: [b.map(|x| x.lo), b.map(|x| x.hi)], params, t: ctx.t_span }
    }

    /// ctx narrowed to the cell's scope.
    pub fn scope(&self, ctx: &Ctx) -> Ctx {
        let mut c = self.params.iter().fold(ctx.clone(), |c, &(k, span)| c.with_param(k, span));
        if self.t.is_some() {
            c.t_span = self.t;
        }
        c.into_scope()
    }
}

/// How many cells a proof that needed its tolerance keeps (Proof::loose).
pub const LOOSE: usize = 64;

/// Confirms candidate witnesses: whether a tape is below a level at a point
/// (at time ctx.t) by bounds rounded outward. Searches find witnesses with
/// the f64 evaluator, whose value can be a rounding below the level where
/// the real value is above it (x x - (2 x - 1) rounds to 0 at x = 1 + 2^-27,
/// really 2^-54); a witness counts only once this confirms it. ctx.t must
/// lie in the span a simplified tape was simplified over. Made once per leaf
/// or search: the point context and the evaluator's storage are reused.
struct Confirm {
    at: Ctx,
    in_span: bool,
    slots: Vec<I>,
}

impl Confirm {
    /// At ctx's point (ctx.t, ctx.params), which must lie in its spans; a
    /// ranged param keeps the meaning it has in the proof (Ctx::point).
    fn new(ctx: &Ctx) -> Confirm {
        Confirm { at: ctx.point(), in_span: ctx.in_scope(), slots: Vec::new() }
    }
    fn below(&mut self, tape: &Tape, p: V3, level: f64) -> bool {
        self.in_span && ink_core::eval::eval_interval_in(tape, &self.at, p.map(I::point), None, &mut self.slots).hi < level
    }
    /// Every point of the box is inside (below 0).
    fn inside(&mut self, tape: &Tape, b: [I; 3]) -> bool {
        self.in_span && ink_core::eval::eval_interval_in(tape, &self.at, b, None, &mut self.slots).hi < 0.0
    }
}

/// An upper bound on the distance between two points (rounded up). A
/// difference of 0 is exact: rounded up, the distance between a point and
/// itself would be 1.5e-323, reported as an overlapping clearance's value.
pub fn distance_up(p: V3, q: V3) -> f64 {
    ink_core::interval::norm(std::array::from_fn(|k| match (p[k] - q[k]).abs() {
        0.0 => 0.0,
        d => d.next_up(),
    }))
}

/// A point of the solid {tape < 0} at p or just under its surface there,
/// confirmed by bounds rounded outward: p itself, else p stepped along the
/// field's descent to a depth of 1e-12 `scale` (a length), growing tenfold
/// to 1e-6 `scale`. None: no point is confirmed (no gradient at p, or the
/// field is not below 0 just past it).
pub fn confirm_inside(tape: &Tape, ctx: &Ctx, p: V3, scale: f64) -> Option<V3> {
    let mut confirm = Confirm::new(ctx);
    let ok = |c: &mut Confirm, x: V3| x.iter().all(|v| v.is_finite()) && c.below(tape, x, 0.0);
    if ok(&mut confirm, p) {
        return Some(p);
    }
    let g = eval_grads(tape, ctx, &[p])[0];
    let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    if !(n > 0.0 && n.is_finite()) {
        return None;
    }
    let mut depth = 1e-12 * scale;
    for _ in 0..7 {
        // Newton's step onto the surface, then depth further along -grad.
        let k = (g[0] / n + depth) / n;
        let x = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        if ok(&mut confirm, x) {
            return Some(x);
        }
        depth *= 10.0;
    }
    None
}

/// Newton steps onto a field's zero set from p, until a step is shorter
/// than `tol` (a length: the field's own scale does not matter). None: the
/// steps did not settle (no gradient, or a field that is not a distance
/// sends them back and forth across a crease).
pub fn to_surface(tape: &Tape, ctx: &Ctx, mut p: V3, tol: f64) -> Option<V3> {
    for _ in 0..30 {
        let g = eval_grads(tape, ctx, &[p])[0];
        let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
        if g[0] == 0.0 {
            return Some(p);
        }
        let k = g[0] / n2;
        if !(n2 > 0.0 && k.is_finite()) {
            return None;
        }
        p = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        if !(k.abs() * n2.sqrt() > tol) {
            return p.iter().all(|x| x.is_finite()).then_some(p);
        }
    }
    None
}

/// A pair of points confirmed inside A and B, for a distance search that
/// found none (`min_distance` walks from B along A's slope, which a field
/// flat away from its solid, a union with a constant, does not have): a
/// point proven inside each (the emptiness search's witness, `finest` its
/// cells, over the boxes `ab` and `bb`), the segment between them cut where
/// it enters each solid, then moved closer (`closer_pair`). Its distance,
/// rounded up, is an upper bound on the solids' minimum distance: the
/// search goes on from it. None: a solid has no point found inside.
pub fn inside_pair(ta: &Arc<Tape>, tb: &Arc<Tape>, ctx: &Ctx, ab: (V3, V3), bb: (V3, V3), finest: f64) -> Option<(f64, V3, V3)> {
    let scale = max_extent(ab.0, ab.1).max(max_extent(bb.0, bb.1));
    let inside = |t: &Arc<Tape>, (lo, hi): (V3, V3)| match prove_empty(t, ctx, lo, hi, None, finest, 0.0).verdict {
        Verdict::Disproven { witness, .. } => confirm_inside(t, ctx, witness, scale),
        _ => None,
    };
    let (p, q) = (inside(ta, ab)?, inside(tb, bb)?);
    // Where the segment from x (outside t's solid) to y (inside) first
    // enters it, confirmed; else y.
    let enter = |t: &Tape, x: V3, y: V3| {
        let at = |s: f64| std::array::from_fn(|k| x[k] + s * (y[k] - x[k]));
        let f = |s: f64| eval_point(t, ctx, at(s));
        match f(0.0) {
            fx if fx >= 0.0 => confirm_inside(t, ctx, at(crossing(&f, Side::Inside, (0.0, fx), (1.0, f(1.0)))), scale).unwrap_or(y),
            _ => y,
        }
    };
    let p = enter(ta, q, p);
    let q = enter(tb, p, q);
    // Each in turn slid along its surface toward the other (`slide`),
    // while that brings them closer.
    let (mut p, mut q) = (p, q);
    let mut dist = distance_up(p, q);
    for _ in 0..64 {
        let p2 = slide(ta, ctx, p, q, scale);
        let q2 = slide(tb, ctx, q, p2, scale);
        let d2 = distance_up(p2, q2);
        if !(d2 < dist) {
            break;
        }
        (dist, p, q) = (d2, p2, q2);
    }
    let (d2, p2, q2) = closer_pair(ta, tb, ctx, p, q, scale);
    Some(if d2 < dist { (d2, p2, q2) } else { (dist, p, q) })
}

/// p, confirmed inside t's solid at its surface, moved along the surface
/// toward x: the part of x - p across t's slope at p, a step back onto the
/// surface (to_surface) and confirmed inside, halved until it comes closer
/// to x; p where none does. Only t's slope at its own surface is read, so
/// a field flat away from its solid moves it.
fn slide(t: &Tape, ctx: &Ctx, p: V3, x: V3, scale: f64) -> V3 {
    let g = eval_grads(t, ctx, &[p])[0];
    let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
    if !(n2 > 0.0 && n2.is_finite()) {
        return p;
    }
    let to: V3 = std::array::from_fn(|k| x[k] - p[k]);
    let along = (to[0] * g[1] + to[1] * g[2] + to[2] * g[3]) / n2;
    let across: V3 = std::array::from_fn(|k| to[k] - along * g[k + 1]);
    let gap = distance_up(p, x);
    let mut h = 1.0;
    for _ in 0..20 {
        let y = std::array::from_fn(|k| p[k] + h * across[k]);
        if let Some(y) = to_surface(t, ctx, y, 1e-12 * scale).and_then(|y| confirm_inside(t, ctx, y, scale))
            && distance_up(y, x) < gap
        {
            return y;
        }
        h *= 0.5;
    }
    p
}

/// A pair of points confirmed inside A and B (p, q) moved closer by
/// alternating projections (q onto A, then that onto B), a move kept only
/// while both points stay confirmed inside (`confirm_inside`, `scale` a
/// length) and the pair comes closer: the pair and its distance rounded
/// up, an upper bound on the solids' minimum distance.
pub fn closer_pair(ta: &Tape, tb: &Tape, ctx: &Ctx, mut p: V3, mut q: V3, scale: f64) -> (f64, V3, V3) {
    let tol = 1e-12 * scale;
    let mut dist = distance_up(p, q);
    for _ in 0..20 {
        let Some(p2) = to_surface(ta, ctx, q, tol).and_then(|x| confirm_inside(ta, ctx, x, scale)) else { break };
        let Some(q2) = to_surface(tb, ctx, p2, tol).and_then(|x| confirm_inside(tb, ctx, x, scale)) else { break };
        let d2 = distance_up(p2, q2);
        if !(d2 < dist) {
            break;
        }
        (dist, p, q) = (d2, p2, q2);
    }
    (dist, p, q)
}

/// Prove that a field's solid S in lo..hi (where it is below 0, as clipped
/// to lo..hi: a point past it is outside) holds no point deeper than `tol`
/// (every point of S within tol of a point outside it; 0: no point at all),
/// or find a witness (a point whose cube of half-side tol is inside: so is
/// the ball of radius tol about it, and it is deeper): the one search
/// (crate::search) of boxes of the material's box (a proven bound of S in
/// lo..hi), halved down to cells of the floor (positive: the claim's
/// precision, half its tolerance, or for an exact question a floor set by
/// the whole model's size), or of a `res` given where that is finer. Its
/// boxes lie on a lattice whose cells are the first cells (`first_cells`:
/// 1/32 of the material's box, coarsened to the floor times a power of two)
/// halved to the floor, so a claim over a whole part is decided as finely
/// as one over a feature, and a region under 32 floors keeps its own 1/32,
/// finer than the floor. Best first, the box
/// where the field can be lowest is refined first, so a witness is
/// usually found in the first few hundred boxes; a box undecided with
/// every side at the finest cells ends the proof (Ended::Finest), with the
/// accuracy it holds to (`attained`: no point of S deeper than that, from
/// the boxes left; at or under `tol`, it is proven to it). Else the search
/// hunts a point deeper than `tol` in the rest (search::Search::hunt), and
/// is unknown where it finds none. At most BOXES boxes are tested
/// (Ended::Limit).
/// Returns the answer, the cells it reached (the smallest side of a box
/// tested) and whether a proof needed the tolerance (`exact` false).
pub fn prove_empty(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64) -> Proof {
    empty_now(tape, ctx, lo, hi, res, floor, tol).0
}

/// `prove_empty`, and the boxes it tested.
fn empty_now(tape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, res: Option<f64>, floor: f64, tol: f64) -> (Proof, usize) {
    assert!(floor > 0.0, "an emptiness search needs a positive floor to stop at");
    let ctx = Arc::new(in_span(ctx));
    let proven = Proof { verdict: Verdict::Proven, res: max_extent(lo, hi) / 32.0, tol, exact: true, loose: Vec::new() };
    let Some((blo, bhi, _)) = bounds_quick(tape, &ctx, lo, hi) else { return (proven, 0) };
    let (first, floor) = first_cells(max_extent(blo, bhi), res, floor);
    // Narrowed once more, not again inside the narrowed box.
    let Some((blo, bhi, _)) = bounds_quick(tape, &ctx, blo, bhi) else { return (proven, 0) };
    let mut finest = first;
    while finest > floor {
        finest *= 0.5;
    }
    let Some(root) = lattice(blo, bhi, finest) else {
        return (
            Proof {
                verdict: Verdict::Unknown { cells: 1, volume: f64::INFINITY, ended: Ended::Limit, cell: None, attained: None, hunt: None },
                res: finest,
                tol,
                exact: false,
                loose: Vec::new(),
            },
            0,
        );
    };
    let q =
        Empty { tol, region: std::array::from_fn(|k| I::new(blo[k], bhi[k])), asked: std::array::from_fn(|k| I::new(lo[k], hi[k])), loose: Default::default(), t_band: ctx.t_span };
    let (searched, s) = empty_run(&q, vec![search::Node { b: root, ctx: ctx.clone(), state: (Kept::of(tape), Kept::of(tape)) }], finest);
    (searched, s.tested())
}

/// The emptiness search `q` from `roots` (boxes of a lattice of cells
/// `finest`), and its answer.
fn empty_run(q: &Empty, roots: Vec<search::Node<Pieces>>, finest: f64) -> (Proof, search::Search<Pieces, Witness>) {
    // The lattice's sides are finest times powers of two: a side up to half
    // as wide again is at the finest cells (its ends rounded).
    let mut s = search::Search::new(q, roots, search::Precision::holds(1.5 * finest));
    let mut outcome = s.finish(q, BOXES);
    // Not proven, and not to the tolerance by the boxes left (`attained`):
    // a point of S deeper than it may be there, so the search hunts one.
    let mut attained = None;
    if let search::Outcome::Unknown(u) = &outcome
        && !u.stopped()
    {
        let a = q.attained(&s.left().collect::<Vec<_>>());
        attained = Some(a);
        if !(q.tol > 0.0 && a <= q.tol) {
            s.hunt();
            outcome = s.finish(q, BOXES);
            attained = None;
        }
    }
    // The cells it reached: the lattice's, the smallest side tested (as
    // rounded) taken to the power of two of the finest it is.
    let res = finest * (s.res() / finest).log2().round().exp2();
    let loose: Vec<Loose> = q.loose.lock().unwrap().iter().map(|c| c.1.clone()).collect();
    let verdict = match outcome {
        search::Outcome::Proven => Verdict::Proven,
        search::Outcome::Found((witness, at, value)) => Verdict::Disproven { witness, t: at.t, value, params: at.params },
        search::Outcome::Unknown(u) => {
            let left: Vec<&search::Node<Pieces>> = s.left().collect();
            let attained = attained.or_else(|| (!u.stopped()).then(|| q.attained(&left)));
            if attained.is_some_and(|a| q.tol > 0.0 && a <= q.tol) {
                // Every box left is held to the tolerance only.
                let loose = left.iter().filter_map(|n| Some(Loose::of(q.clip(n.b)?, &n.ctx))).take(LOOSE).collect();
                return (Proof { verdict: Verdict::Proven, res, tol: q.tol, exact: false, loose }, s);
            }
            let volume = left.iter().filter_map(|n| q.clip(n.b)).map(|b| b.iter().map(|x| x.hi - x.lo).product::<f64>()).sum();
            let cell = u.cell.and_then(|[a, b]| q.clip(std::array::from_fn(|k| I::new(a[k], b[k])))).map(|b| [b.map(|x| x.lo), b.map(|x| x.hi)]);
            Verdict::Unknown { cells: u.cells, volume, ended: u.ended, cell, attained, hunt: u.hunt }
        }
    };
    (Proof { verdict, res, tol: q.tol, exact: loose.is_empty(), loose }, s)
}

/// The cube of a lattice of cells of side h anchored at lo that holds
/// lo..hi: its side h times the least power of two that reaches hi on
/// every axis (a cube, so halving every axis at once keeps it one: an
/// octree's cells). None where a side is not finite or needs more than
/// 2^62 cells.
fn lattice(lo: V3, hi: V3, h: f64) -> Option<[I; 3]> {
    let mut side = h;
    for _ in 0..=62 {
        let ends: V3 = std::array::from_fn(|k| lo[k] + side);
        if !(0..3).all(|k| lo[k].is_finite() && ends[k].is_finite()) {
            return None;
        }
        if (0..3).all(|k| ends[k] >= hi[k]) {
            return Some(std::array::from_fn(|k| I::new(lo[k], ends[k])));
        }
        side *= 2.0;
    }
    None
}

/// Boxes a search may test (an emptiness's or a clearance's).
pub const BOXES: usize = 1 << 22;

/// Witness samples a box confirms against the tolerance (the deepest
/// first: each costs a bound over a cube).
const DEEP_TRIES: usize = 4;

/// A point inside, its moment and ranged params, and the field's value there.
type Witness = (V3, At, f64);

/// What a box of an emptiness search keeps for its halves: the field
/// simplified over it, and over it grown by the tolerance (where a
/// witness's cube and the points that cover it are bounded).
type Pieces = (Kept, Kept);

/// A tape a box keeps for its halves: the tape its test read and the
/// choices made there, a step at a time, simplified when one of its halves
/// is first tested (cloned for them, they share it). A box waiting holds
/// choices, not a tape: the tapes alive are those of boxes whose halves
/// wait (an emptiness's eight, depth first; a clearance's boxes, best
/// first, wait apart from their siblings and keep their own) and of those
/// being tested.
#[derive(Clone)]
struct Kept(Arc<std::sync::Mutex<(Arc<Tape>, Vec<Vec<u8>>)>>);

impl Kept {
    fn of(t: &Arc<Tape>) -> Kept {
        Kept::after(t, Vec::new())
    }

    /// `from` simplified by each of `steps` in turn.
    fn after(from: &Arc<Tape>, steps: Vec<Vec<u8>>) -> Kept {
        let steps = steps.into_iter().filter(|ch| ch.iter().any(|&c| c != BOTH)).collect();
        Kept(Arc::new(std::sync::Mutex::new((from.clone(), steps))))
    }

    /// The tape, simplified once and kept.
    fn tape(&self) -> Arc<Tape> {
        let mut k = self.0.lock().unwrap();
        for ch in std::mem::take(&mut k.1) {
            k.0 = simple(&k.0, &ch);
        }
        k.0.clone()
    }

    /// The tape, not kept (a box read once, as the search ends).
    fn read(&self) -> Arc<Tape> {
        let k = self.0.lock().unwrap();
        k.1.iter().fold(k.0.clone(), |t, ch| simple(&t, ch))
    }

    /// The tape it holds: its own once made, else the one it is made from.
    #[cfg(test)]
    fn held(&self) -> *const Tape {
        Arc::as_ptr(&self.0.lock().unwrap().0)
    }
}

/// Emptiness as a question of the one search: no point of the solid S in
/// the material's box `region` (a proven bound of S in the region asked,
/// `asked`, where the tape holds) deeper than `tol`. Boxes are a lattice's
/// (prove_empty), each tested as clipped to `region`. A box is settled
/// where a test passes, each sound for any field:
/// - the field's bounds over it are at least 0 (over the scope's axes,
///   ranged params and time over a step of a motion, also first order in
///   them, and the facing form, whose pairs cancel the rates of two faces
///   one param or a motion moves together: a shaft and its bore of one
///   diameter);
/// - its first-order bound is, where no sample of it is inside (across a
///   gap thinner than the box, or at a branch cut);
/// - every point of it lies within the tolerance of a point shown outside
///   S: by a pair of faces facing each other along an axis, at any box
///   size (centred::centred_settles: faces that touch, or overlap thinner
///   than the tolerance), or by points shown outside near each of its
///   samples (`covered`, at boxes of a few tolerances): none of its solid
///   is deeper (`loose`: a proof that needed the tolerance, which keeps
///   the LOOSE such boxes of lowest bound).
///
/// A witness is its centre where the field is below 0 all over it and the
/// cube about it fits (or is shown inside), or one of its samples, the
/// deepest first, whose cube is shown inside: at the box's moment and
/// ranged params' middle values (its own: time halved inside a step finds
/// a witness at a moment of its half), or each param pushed to the end of
/// its span that lowers the field there (scope::pushed: a stack's worst
/// corner).
///
/// Open, a box proves no depth short of the tolerance (its depth is
/// bounded only for the boxes left: `attained`), so its bound is 0, and
/// open boxes tie: the smaller goes first, so the search is depth first,
/// as interval solvers prove a set empty (IBEX, RealPaver), its memory as
/// its depth, not as a contact's area (best first on the field's lower
/// bound refines a contact breadth first, every box near it bounded
/// alike), and a box no size decides is reached in as many rounds as
/// halvings. Witnesses are sought in every box tested. Its smear on each space
/// axis is its half-width on the lattice (the field's rate taken alike
/// along each: a cube is halved as an octree's cell is); where it is
/// decided at the middle of the scope's axes but not over their spans,
/// what the spans add is the cause, and space's smear is the field's
/// spread over the box (each axis's rate bound times its half-width),
/// shared out by width, each param's and time's its blur of the bound
/// (`scope_shares`): a step's moments are halved where they blur the box
/// more than its size does. Covered at the middle but not over the spans,
/// the axis halved least so far is halved first.
struct Empty {
    tol: f64,
    region: [I; 3],
    asked: [I; 3],
    loose: std::sync::Mutex<Vec<(f64, Loose)>>,
    /// Time's span over the search (a step of a motion), which a box's
    /// share of it is of.
    t_band: Option<(f64, f64)>,
}

impl search::Question for Empty {
    type State = Pieces;
    type Witness = Witness;

    fn test(&self, n: &search::Node<Pieces>) -> search::Test<Pieces, Witness> {
        use search::Test::{Found, Open, Settled};
        let Some(b) = self.clip(n.b) else { return Settled };
        let ctx = &*n.ctx;
        let (from, given) = (n.state.0.tape(), n.state.1.tape());
        let mut ch = vec![BOTH; from.len()];
        let r = eval_interval(&from, ctx, b, Some(&mut ch));
        if r.lo >= 0.0 {
            return Settled;
        }
        let mut sub = simple(&from, &ch);
        // `wide` simplified over the box grown by the tolerance where the
        // box's own tape has halved the one it was tested with (narrowing
        // costs a bound over the grown box).
        let mut chw = vec![BOTH; given.len()];
        let wide = match self.tol > 0.0 {
            false => sub.clone(),
            true if 2 * sub.len() <= given.len() => {
                eval_interval(&given, ctx, self.grown(b), Some(&mut chw));
                simple(&given, &chw)
            }
            true => given.clone(),
        };
        let mut chp = Vec::new();
        let c: V3 = std::array::from_fn(|k| b[k].mid());
        // Inside throughout: its centre is a witness where the box holds the
        // cube about it, else where that cube is inside.
        if r.hi < 0.0 && ((0..3).all(|k| c[k] - b[k].lo >= self.tol && b[k].hi - c[k] >= self.tol) || self.deep(&sub, &wide, &mut Confirm::new(ctx), c)) {
            return Found((c, At::of(ctx), eval_point(&sub, ctx, c)));
        }
        // The scope's axes the field reads: ranged params, and time over a
        // step of a motion.
        let (ps, ts) = (scope::splittable(&sub, ctx, scope::PARAM_SPLITS), scope::spans_t(&sub, ctx));
        let mid = (!ps.is_empty() || ts).then(|| scope::middle(ctx, ts));
        let (mut params, mut t, mut rates, mut centred) = (Vec::new(), 0.0, None, false);
        if let Some(mid) = &mid {
            let m = eval_interval(&sub, mid, b, None);
            if m.hi < 0.0 && self.deep(&sub, &wide, &mut Confirm::new(mid), c) {
                return Found((c, At::of(mid), eval_point(&sub, mid, c)));
            }
            if m.lo >= 0.0 {
                chp = vec![BOTH; sub.len()];
                let (r1, blur) = ink_core::eval::eval_interval_scope_blur(&sub, ctx, b, Some(&mut chp));
                if r1.lo >= 0.0 {
                    return Settled;
                }
                sub = simple(&sub, &chp);
                if crate::centred::centred_at_least(&[(&sub, 0.0)], ctx, b, 0.0) {
                    return Settled;
                }
                centred = true;
                let (shares, tb) = scope_shares(ctx, &[(&blur, 1.0)]);
                params = ps.iter().map(|&k| (k, share(&shares, k))).collect();
                t = if ts { tb } else { 0.0 };
                rates = Some(ink_core::eval::eval_interval_grad(&sub, ctx, b));
            }
        }
        let pts = box_samples(b, 4);
        let mut v = Vec::new();
        eval_points(&sub, ctx, &pts, &mut v);
        let within = |k: usize| (0..3).all(|a| b[a].contains(pts[k][a]));
        // The deepest samples below 0.
        let mut below: Vec<usize> = (0..pts.len()).filter(|&k| v[k] < 0.0 && within(k)).collect();
        below.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        let tries = if self.tol == 0.0 { below.len() } else { DEEP_TRIES };
        let mut confirm = (!below.is_empty()).then(|| Confirm::new(ctx));
        if let Some(&k) = below.iter().take(tries).find(|&&k| confirm.as_mut().is_some_and(|c| self.deep(&sub, &wide, c, pts[k]))) {
            return Found((pts[k], At::of(ctx), v[k]));
        }
        // A sample inside leaves no lower bound over the box to try (a
        // speed rule: skipping a bound never makes a wrong answer); with a
        // tolerance, a pair of faces that touch or overlap thinner than it
        // still covers the box, whatever its size.
        use crate::centred::Settles;
        let bounded = below.is_empty() && !centred;
        let settles = match self.tol > 0.0 {
            true => crate::centred::centred_settles(&sub, ctx, b, self.tol, bounded),
            false if bounded && crate::centred::centred_at_least(&[(&sub, 0.0)], ctx, b, 0.0) => Settles::AtLeast,
            false => Settles::No,
        };
        match settles {
            Settles::AtLeast => return Settled,
            Settles::Covered => {
                self.keep_loose(r.lo, b, ctx);
                return Settled;
            }
            Settles::No => {}
        }
        let sides = b.map(|x| x.hi - x.lo);
        if self.covered(&sub, &wide, ctx, b, sides, &pts, &v) {
            self.keep_loose(r.lo, b, ctx);
            return Settled;
        }
        // Covered at the scope's middle but not over it: the axis of it
        // halved least so far is halved first (time's halvings of the
        // search's span).
        if let Some(mid) = mid.as_ref().filter(|_| params.is_empty() && t == 0.0)
            && self.covered(&sub, &wide, mid, b, sides, &pts, &v)
        {
            let halved = ps.iter().map(|&k| (Some(k), scope::halvings(ctx, k))).chain(ts.then(|| (None, self.t_halvings(ctx))));
            match halved.min_by(|x, y| x.1.total_cmp(&y.1)).map(|x| x.0) {
                Some(Some(k)) => params.push((k, f64::INFINITY)),
                Some(None) => t = f64::INFINITY,
                None => {}
            }
        }
        if !ps.is_empty() || !scope::read_params(&sub, ctx).is_empty() {
            let mut order: Vec<usize> = (0..pts.len()).filter(|&k| within(k)).collect();
            order.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
            for &k in order.iter().take(4) {
                let Some((pushed, x)) = scope::pushed(&[&sub], ctx, pts[k]) else { continue };
                if x[0] < 0.0 && self.deep(&sub, &wide, &mut Confirm::new(&pushed), pts[k]) {
                    return Found((pts[k], At::of(&pushed), x[0]));
                }
            }
        }
        // Space is halved as an octree halves its cells (on the lattice's
        // widths: cubes, which the tolerance's cover settles), its spread
        // shared out by width where params compete with it.
        let hw = n.b.map(|x| 0.5 * (x.hi - x.lo));
        let space = match &rates {
            Some(g) => {
                let spread = (0..3).map(|k| ink_core::interval::mag(g.d[k]) * 0.5 * sides[k]).sum::<f64>();
                hw.map(|h| h * spread / hw.iter().sum::<f64>())
            }
            None => hw,
        };
        // The halves start from the field simplified over the box, and over
        // it grown by the tolerance only where that has halved the one it
        // was tested with (no tolerance reads no `wide`).
        let wide = if self.tol > 0.0 && 2 * wide.len() <= given.len() { Kept::after(&given, vec![chw]) } else { n.state.1.clone() };
        Open { bound: 0.0, smear: search::Smear { space, params, t }, state: (Kept::after(&from, vec![ch, chp]), wide) }
    }

    /// From the deepest samples of the box, a step down the field's slope:
    /// asked exactly, to just under its surface (`confirm_inside`: an
    /// overlap thinner than the samples' spacing, which none lies in);
    /// with a tolerance, to where the cube about it would fit were the
    /// field a distance (its half-diagonal below 0), confirmed (`deep`).
    fn seek(&self, n: &search::Node<Pieces>) -> Option<Witness> {
        let b = self.clip(n.b)?;
        let ctx = &*n.ctx;
        let (sub, wide) = (&n.state.0.tape(), &n.state.1.tape());
        let pts = box_samples(b, 4);
        let mut v = Vec::new();
        eval_points(sub, ctx, &pts, &mut v);
        let mut order: Vec<usize> = (0..pts.len()).collect();
        order.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        let inside = |p: &V3| (0..3).all(|a| b[a].contains(p[a]));
        let mut confirm = Confirm::new(ctx);
        order.into_iter().take(DEEP_TRIES).find_map(|k| {
            let p = match self.tol {
                0.0 => confirm_inside(sub, ctx, pts[k], max_extent(b.map(|x| x.lo), b.map(|x| x.hi))).filter(inside)?,
                tol => {
                    let [f, gx, gy, gz] = eval_grads(sub, ctx, &[pts[k]])[0];
                    let g2 = gx * gx + gy * gy + gz * gz;
                    let s = (f + 1.01 * 3f64.sqrt() * tol * g2.sqrt()) / g2;
                    let p = [pts[k][0] - s * gx, pts[k][1] - s * gy, pts[k][2] - s * gz];
                    (s.is_finite() && inside(&p) && self.deep(sub, wide, &mut confirm, p)).then_some(p)?
                }
            };
            Some((p, At::of(ctx), eval_point(sub, ctx, p)))
        })
    }
}

impl Empty {
    /// Keep box `b` over ctx's scope, settled only by the tolerance with
    /// `lo` its field's lower bound, among the LOOSE of lowest bound (ties
    /// by place, then scope, so the cells kept do not depend on the order
    /// threads settle them in).
    fn keep_loose(&self, lo: f64, b: [I; 3], ctx: &Ctx) {
        let key = |c: &(f64, Loose)| -> Vec<f64> {
            let [l, h] = c.1.cell;
            let scope = c.1.t.iter().chain(c.1.params.iter().map(|p| &p.1)).flat_map(|s| [s.0, s.1]);
            [c.0, l[0], l[1], l[2], h[0], h[1], h[2]].into_iter().chain(scope).collect()
        };
        let cell = (lo, Loose::of(b, ctx));
        let new = key(&cell);
        let mut kept = self.loose.lock().unwrap();
        let at = kept.partition_point(|c| key(c).iter().zip(&new).map(|(x, y)| x.total_cmp(y)).find(|o| o.is_ne()).is_none_or(|o| o.is_lt()));
        if at < LOOSE {
            kept.insert(at, cell);
            kept.truncate(LOOSE);
        }
    }

    /// How many times time's span has been halved in ctx from the search's
    /// (log2 of their widths' ratio).
    fn t_halvings(&self, ctx: &Ctx) -> f64 {
        let w = |s: Option<(f64, f64)>| s.map_or(1.0, |(a, b)| b - a);
        (w(self.t_band) / w(ctx.t_span)).log2()
    }

    /// A box of the lattice as clipped to the material's box (None: it
    /// holds none of it).
    fn clip(&self, b: [I; 3]) -> Option<[I; 3]> {
        let c: [I; 3] = std::array::from_fn(|k| I { lo: b[k].lo.max(self.region[k].lo), hi: b[k].hi.min(self.region[k].hi) });
        c.iter().all(|x| x.lo <= x.hi).then_some(c)
    }

    /// Box `b` grown by the tolerance, within the region asked (where the
    /// tape holds).
    fn grown(&self, b: [I; 3]) -> [I; 3] {
        std::array::from_fn(|k| I::new((b[k].lo - self.tol).next_down().max(self.asked[k].lo), (b[k].hi + self.tol).next_up().min(self.asked[k].hi)))
    }

    /// Whether p (in box b) is deeper than the tolerance: the cube of
    /// half-side tol about it is inside, at ctx's point of the scope
    /// (`wide`: the field over b grown by tol); with no tolerance, p itself
    /// is inside. The cube must lie in the material's box (outside it there
    /// is none).
    fn deep(&self, sub: &Tape, wide: &Tape, confirm: &mut Confirm, p: V3) -> bool {
        if self.tol == 0.0 {
            return confirm.below(sub, p, 0.0);
        }
        let cube: [I; 3] = std::array::from_fn(|k| I::new((p[k] - self.tol).next_down(), (p[k] + self.tol).next_up()));
        (0..3).all(|k| self.region[k].lo <= cube[k].lo && cube[k].hi <= self.region[k].hi) && confirm.inside(wide, cube)
    }

    /// Whether every point of box b (sides `side`) lies within the
    /// tolerance of a point shown outside the solid (the field at least 0
    /// there by bounds rounded outward, over the whole scope), so none of
    /// the box's solid is deeper than it. Needs no slope: it holds for any
    /// field. The points may lie beyond the box, within the tolerance of it
    /// (`wide` holds there): a box inside a contact thinner than twice the
    /// tolerance is settled from outside it. `pts` are the box's samples,
    /// the centres of its 4 x 4 x 4 sub-boxes (`box_samples`), `v` the field
    /// there (`sub` holds only in the box).
    fn covered(&self, sub: &Tape, wide: &Tape, ctx: &Ctx, b: [I; 3], side: V3, pts: &[V3], v: &[f64]) -> bool {
        if self.tol == 0.0 {
            return false;
        }
        let tol = self.tol;
        let grown = self.grown(b);
        let within = |p: &V3, r: &[I; 3]| (0..3).all(|k| r[k].contains(p[k]));
        // Shown outside: in the box by its own tape; past the region asked
        // or the material's box (a proven bound), by them; else in the
        // grown box by `wide`.
        let past = |q: &V3| !within(q, &self.asked) || !within(q, &self.region);
        let shown = |q: &V3, slots: &mut Vec<I>| -> bool {
            let t = if within(q, &b) {
                sub
            } else if past(q) {
                return true;
            } else if within(q, &grown) {
                wide
            } else {
                return false;
            };
            ink_core::eval::eval_interval_in(t, ctx, q.map(I::point), None, slots).lo >= 0.0
        };
        // A point shown outside within the tolerance of every point of the
        // box about c (field x there): c itself, else a step from c up the
        // field's slope to where it would reach 0 were it a distance, else
        // the farthest the box allows up and down the slope and along each
        // axis (the nearest outside may not be up the slope: the side of an
        // overlap whose field is the lower, scaled). The steps are guesses;
        // the check is not.
        let near = |c: V3, x: f64, half: V3, slots: &mut Vec<I>| -> Option<V3> {
            if x >= 0.0 && far(c, c, half, b) <= tol && shown(&c, slots) {
                return Some(c);
            }
            let reach = tol - ink_core::interval::norm(half);
            if !(reach > 0.0) {
                return None;
            }
            let [_, gx, gy, gz] = ink_core::eval::eval_grads(sub, ctx, &[c])[0];
            // (A slope that is 0 or not finite gives no step up it: its
            // points are NaN, which no test below may take as near or
            // outside, so they are dropped here.)
            let n = (gx * gx + gy * gy + gz * gz).sqrt();
            let step = ((-x / n).max(0.0) * 1.01 + 0.01 * reach).min(reach);
            let up = [gx / n, gy / n, gz / n];
            let axis = |k: usize, s: f64| -> V3 { std::array::from_fn(|a| if a == k { s } else { 0.0 }) };
            let ways = [
                (up, step),
                (up, reach),
                (up.map(|u| -u), reach),
                (axis(0, 1.0), reach),
                (axis(0, -1.0), reach),
                (axis(1, 1.0), reach),
                (axis(1, -1.0), reach),
                (axis(2, 1.0), reach),
                (axis(2, -1.0), reach),
            ];
            ways.into_iter()
                .map(|(u, s)| [c[0] + s * u[0], c[1] + s * u[1], c[2] + s * u[2]])
                .filter(|q| q.iter().all(|x| x.is_finite()))
                .find(|q| far(*q, c, half, b) <= tol && (past(q) || within(q, &grown) && eval_point(wide, ctx, *q) >= 0.0) && shown(q, slots))
        };
        let mut slots = Vec::new();
        // First the fewest boxes: the n x n x n whose corners are no more
        // than the tolerance from their centres (one for a box under about
        // the tolerance across, 8 for one twice that).
        let n = (side.iter().fold(0.0f64, |m, &x| m.max(x)) * 3f64.sqrt() / (2.0 * tol)).ceil().max(1.0);
        if n < 4.0 {
            let n = n as usize;
            let part: V3 = side.map(|x| x / n as f64);
            let mut grid = Vec::with_capacity(n * n * n);
            for a in 0..n {
                for c in 0..n {
                    for d in 0..n {
                        grid.push([b[0].lo + part[0] * (d as f64 + 0.5), b[1].lo + part[1] * (c as f64 + 0.5), b[2].lo + part[2] * (a as f64 + 0.5)]);
                    }
                }
            }
            let mut g = Vec::new();
            eval_points(sub, ctx, &grid, &mut g);
            let half = part.map(|x| 0.5 * x);
            if grid.iter().zip(&g).all(|(q, &x)| near(*q, x, half, &mut slots).is_some()) {
                return true;
            }
        }
        // Else each sub-box by a point already shown outside where one is
        // near enough (one covers its neighbours), else a sample outside
        // (its own centre first, then the nearest the box's centre), else
        // a step from its centre. The deepest first: a box that cannot be
        // covered fails at once.
        let centre: V3 = std::array::from_fn(|k| 0.5 * b[k].lo + 0.5 * b[k].hi);
        let mut order: Vec<usize> = (0..pts.len()).filter(|&k| v[k] >= 0.0).collect();
        let d2 = |k: usize| (0..3).map(|a| (pts[k][a] - centre[a]).powi(2)).sum::<f64>();
        order.sort_by(|&i, &j| d2(i).total_cmp(&d2(j)));
        let half: V3 = side.map(|x| 0.125 * x);
        // Samples outside, confirmed (lazily: None not yet asked).
        let mut out: Vec<Option<bool>> = vec![None; pts.len()];
        let mut shown_pts: Vec<V3> = Vec::new();
        let mut deepest: Vec<usize> = (0..pts.len()).collect();
        deepest.sort_by(|&i, &j| v[i].total_cmp(&v[j]));
        deepest.into_iter().all(|j| {
            if shown_pts.iter().any(|&q| far(q, pts[j], half, b) <= tol) {
                return true;
            }
            let sample = std::iter::once(j)
                .chain(order.iter().copied().filter(|&k| k != j))
                .find(|&k| v[k] >= 0.0 && far(pts[k], pts[j], half, b) <= tol && *out[k].get_or_insert_with(|| shown(&pts[k], &mut slots)));
            match sample.map(|k| pts[k]).or_else(|| near(pts[j], v[j], half, &mut slots)) {
                Some(q) => {
                    shown_pts.push(q);
                    true
                }
                None => false,
            }
        })
    }

    /// How far a point of the material's box in box b (as clipped) can lie
    /// from the material box's faces, beyond which S has no point: per
    /// axis, the most the nearer face can be away, and of the axes the
    /// least (rounded up). Per axis that is the most of min(p - lo, hi -
    /// p) over the box's p: half the material box's width where its real
    /// middle is in the box, else the far side's distance to the face
    /// beyond it, so the least of the three, never a sample at a float
    /// middle (which may lie on a face: 1e16 + 1 is no float).
    fn to_faces(&self, b: [I; 3]) -> f64 {
        use ink_core::interval::{add_up, mul_up};
        (0..3)
            .map(|k| {
                let (r, x) = (self.region[k], b[k]);
                let half = mul_up(add_up(r.hi, -r.lo), 0.5);
                add_up(x.hi, -r.lo).min(add_up(r.hi, -x.lo)).min(half).max(0.0)
            })
            .fold(f64::INFINITY, f64::min)
    }

    /// How deep a point of S in box n can be (its distance from the
    /// outside): for each of the box's 4 x 4 x 4 sub-boxes, the farthest
    /// it reaches from the nearest of its samples shown outside (the field
    /// at least 0 there over the whole scope), or from the material box's
    /// faces where nearer; the most of those.
    fn depth(&self, n: &search::Node<Pieces>) -> f64 {
        let Some(b) = self.clip(n.b) else { return 0.0 };
        let (ctx, sub) = (&*n.ctx, &n.state.0.read());
        let pts = box_samples(b, 4);
        let half: V3 = b.map(|x| 0.125 * (x.hi - x.lo));
        let mut v = Vec::new();
        eval_points(sub, ctx, &pts, &mut v);
        let mut slots = Vec::new();
        let within = |p: &V3| (0..3).all(|k| b[k].contains(p[k]));
        let out: Vec<V3> = (0..pts.len())
            .filter(|&k| v[k] >= 0.0 && within(&pts[k]) && ink_core::eval::eval_interval_in(sub, ctx, pts[k].map(I::point), None, &mut slots).lo >= 0.0)
            .map(|k| pts[k])
            .collect();
        pts.iter()
            .map(|&c| {
                let around: [I; 3] = std::array::from_fn(|k| I::new((c[k] - half[k]).next_down().max(b[k].lo), (c[k] + half[k]).next_up().min(b[k].hi)));
                out.iter().map(|&q| far(q, c, half, b)).fold(self.to_faces(around), f64::min)
            })
            .fold(0.0, f64::max)
    }

    /// The accuracy boxes `left` hold to: the deepest a point of S in them
    /// can be (`depth`), each bounded first by the material box's faces
    /// alone, and further only where that bound is more than the most found
    /// so far (it cannot raise the most otherwise).
    fn attained(&self, left: &[&search::Node<Pieces>]) -> f64 {
        let mut order: Vec<(f64, usize)> = left.iter().enumerate().map(|(i, n)| (self.clip(n.b).map_or(0.0, |b| self.to_faces(b)), i)).collect();
        order.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
        let mut most = 0.0f64;
        for chunk in order.chunks(search::BATCH) {
            if chunk[0].0 <= most {
                break;
            }
            let ds: Vec<f64> = chunk.par_iter().map(|&(f, i)| if f <= most { f } else { self.depth(left[i]) }).collect();
            most = ds.into_iter().fold(most, f64::max);
        }
        most
    }
}

/// The farthest a point of the box about c (half-sides `half`) is from q,
/// rounded up, with room for the coordinates' rounding (sub-boxes tile box
/// b to within a few ulps of its coordinates).
fn far(q: V3, c: V3, half: V3, b: [I; 3]) -> f64 {
    ink_core::interval::norm(std::array::from_fn(|k| {
        let ulps = 16.0 * f64::EPSILON * c[k].abs().max(q[k].abs()).max(b[k].lo.abs()).max(b[k].hi.abs());
        ((q[k] - c[k]).abs().next_up() + half[k] + ulps).next_up().next_up()
    }))
}

/// The first cells and the floor of an emptiness search (prove_empty)
/// over a region `extent` across that halves down to `floor`, or to a `res`
/// given where that is finer: 1/32 of the region, coarsened (under twice)
/// to the floor times a power of two, so the halvings land on the floor
/// itself, not on cells just above it and then again just below; doubled
/// while no coarser than a res given. A region under 32 floors keeps its
/// own 1/32.
fn first_cells(extent: f64, res: Option<f64>, floor: f64) -> (f64, f64) {
    let floor = res.map_or(floor, |r| r.min(floor));
    let mut h = extent / 32.0;
    if floor > 0.0 && h > floor {
        h = floor * (h / floor).log2().ceil().exp2();
    }
    if let Some(r) = res {
        while 2.0 * h <= r {
            h *= 2.0;
        }
    }
    (h, floor)
}

/// `ctx` with its moment and ranged params inside their spans: a witness
/// carries ctx.t and ctx.params, which must be in what the proof covered.
fn in_span(ctx: &Ctx) -> Ctx {
    ctx.into_scope()
}

// ---------------------------------------------------------------- slices

#[derive(Clone, Debug)]
pub struct Loop {
    pub pts: Vec<[f64; 2]>,
    pub closed: bool,
    /// Signed area: positive for material boundaries (counter-clockwise),
    /// negative for holes.
    pub area: f64,
    pub perimeter: f64,
}

/// Plane axes `(u, v)` for a slice normal to `axis`.
pub fn slice_axes(axis: usize) -> (usize, usize) {
    match axis {
        0 => (1, 2),
        1 => (0, 2),
        _ => (0, 1),
    }
}

type EdgeKey = (i64, i64, u8);

struct Seg {
    a: EdgeKey,
    b: EdgeKey,
    pa: [f64; 2],
    pb: [f64; 2],
}

/// Marching squares on an interval-pruned quadtree in the plane `axis = at`.
pub fn slice(tape: &Arc<Tape>, ctx: &Ctx, axis: usize, at: f64, lo: V3, hi: V3, h: f64) -> Vec<Loop> {
    let (u, v) = slice_axes(axis);
    let o = [lo[u] - h, lo[v] - h];
    let ext = (hi[u] - lo[u]).max(hi[v] - lo[v]) + 2.0 * h;
    let mut size = 1i64;
    while (size as f64) * h < ext {
        size *= 2;
    }
    let to3 = |p: [f64; 2]| {
        let mut q = [0.0; 3];
        q[axis] = at;
        q[u] = p[0];
        q[v] = p[1];
        q
    };
    let segs = quad_descend(tape, tape, ctx, &to3, o, h, axis, at, (u, v), [0, 0], size);
    join(segs)
}

/// A mark's slice: where field `g` crosses 0 in the plane `axis = at`,
/// inside `region` (r < 0) and the box lo..hi. Marching squares on `g` (as
/// `slice`), each chain then cut where it leaves them, the cut found by
/// bisection along its segment. A closed chain wholly inside stays closed;
/// the rest are open, with area 0. A mark is where g changes sign, as an
/// outline is where f does: a g touching 0 draws nothing there.
#[allow(clippy::too_many_arguments)]
pub fn slice_mark(g: &Arc<Tape>, region: &Tape, ctx: &Ctx, axis: usize, at: f64, lo: V3, hi: V3, h: f64) -> Vec<Loop> {
    let (u, v) = slice_axes(axis);
    let to3 = |p: [f64; 2]| {
        let mut q = [at; 3];
        q[u] = p[0];
        q[v] = p[1];
        q
    };
    // Inside where negative: the region within the box.
    let r = |p: [f64; 2]| eval_point(region, ctx, to3(p)).max(lo[u] - p[0]).max(p[0] - hi[u]).max(lo[v] - p[1]).max(p[1] - hi[v]);
    // The crossing on a -> b, r(a) and r(b) of opposite signs: well below
    // the cell (slice's own error), 40 halvings at most.
    let cross = |a: [f64; 2], b: [f64; 2], ra: f64| {
        let (mut a, mut b) = (a, b);
        for _ in 0..40 {
            if ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt() < h * 1e-6 {
                break;
            }
            let m = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
            if (r(m) < 0.0) == (ra < 0.0) { a = m } else { b = m }
        }
        [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]
    };
    let length = |pts: &[[f64; 2]]| pts.windows(2).map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt()).sum::<f64>();
    let mut out = Vec::new();
    for l in slice(g, ctx, axis, at, lo, hi, h) {
        let mut pts = l.pts;
        let mut rs: Vec<f64> = pts.iter().map(|&p| r(p)).collect();
        if l.closed {
            // Start a closed chain outside the region, if it leaves it.
            match rs.iter().position(|&x| x >= 0.0) {
                None => {
                    let mut p = pts.clone();
                    p.push(pts[0]);
                    out.push(Loop { perimeter: length(&p), pts, closed: true, area: 0.0 });
                    continue;
                }
                Some(k) => {
                    pts.rotate_left(k);
                    rs.rotate_left(k);
                    pts.push(pts[0]);
                    rs.push(rs[0]);
                }
            }
        }
        let mut cur: Vec<[f64; 2]> = Vec::new();
        for k in 0..pts.len() {
            let inside = rs[k] < 0.0;
            if k > 0 && inside != (rs[k - 1] < 0.0) {
                cur.push(cross(pts[k - 1], pts[k], rs[k - 1]));
                if !inside {
                    out.push(Loop { perimeter: length(&cur), pts: std::mem::take(&mut cur), closed: false, area: 0.0 });
                }
            }
            if inside {
                cur.push(pts[k]);
            }
        }
        if cur.len() >= 2 {
            out.push(Loop { perimeter: length(&cur), pts: cur, closed: false, area: 0.0 });
        }
    }
    out.retain(|l| l.pts.len() >= 2);
    out
}

/// `full`: the field unsimplified, for samples outside a cell (`march`).
#[allow(clippy::too_many_arguments)]
fn quad_descend(
    full: &Tape,
    tape: &Arc<Tape>,
    ctx: &Ctx,
    to3: &(dyn Fn([f64; 2]) -> V3 + Sync),
    o: [f64; 2],
    h: f64,
    axis: usize,
    at: f64,
    uv: (usize, usize),
    i: [i64; 2],
    size: i64,
) -> Vec<Seg> {
    // Same formula as the corner positions in `march`, so bounds agree bit-for-bit.
    let lo = [o[0] + h * i[0] as f64, o[1] + h * i[1] as f64];
    let hi = [o[0] + h * (i[0] + size) as f64, o[1] + h * (i[1] + size) as f64];
    let mut b = [I::point(0.0); 3];
    b[axis] = I::point(at);
    b[uv.0] = I::new(lo[0], hi[0]);
    b[uv.1] = I::new(lo[1], hi[1]);
    let mut ch = vec![BOTH; tape.len()];
    let r = eval_interval(tape, ctx, b, Some(&mut ch));
    if r.lo > 0.0 || r.hi < 0.0 {
        return Vec::new();
    }
    let sub = if ch.iter().any(|&c| c != BOTH) { Arc::new(tape.simplify(&ch)) } else { tape.clone() };
    if size == 1 {
        return march(full, &sub, ctx, to3, o, h, i);
    }
    let hs = size / 2;
    let kids = [[0, 0], [hs, 0], [0, hs], [hs, hs]];
    let parts: Vec<Vec<Seg>> = kids.par_iter().map(|k| quad_descend(full, &sub, ctx, to3, o, h, axis, at, uv, [i[0] + k[0], i[1] + k[1]], hs)).collect();
    parts.into_iter().flatten().collect()
}

fn march(full: &Tape, tape: &Tape, ctx: &Ctx, to3: &(dyn Fn([f64; 2]) -> V3 + Sync), o: [f64; 2], h: f64, i: [i64; 2]) -> Vec<Seg> {
    let p = |a: i64, b: i64| [o[0] + h * a as f64, o[1] + h * b as f64];
    let (x, y) = (i[0], i[1]);
    let corner = [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)];
    let c = corner.map(|(a, b)| p(a, b));
    let center = [o[0] + h * (x as f64 + 0.5), o[1] + h * (y as f64 + 0.5)];
    let mut v = Vec::new();
    eval_points(tape, ctx, &[to3(c[0]), to3(c[1]), to3(c[2]), to3(c[3]), to3(center)], &mut v);
    // The curve is where the field changes sign, between f < 0 and f > 0.
    // A corner at exactly 0 is inside unless a grid neighbour is positive:
    // so the curve runs through it where the field crosses there, and a
    // field touching 0 (negative or 0 on every side) draws nothing, no
    // sliver of zero width. Every cell sharing the corner reads the same
    // neighbours (unsimplified: they lie outside the cell).
    let inside = |k: usize| {
        v[k] < 0.0
            || v[k] == 0.0 && {
                let (a, b) = corner[k];
                let mut n = Vec::new();
                let near: Vec<V3> = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)].iter().map(|&(da, db)| to3(p(a + da, b + db))).collect();
                eval_points(full, ctx, &near, &mut n);
                n.iter().all(|&w| w <= 0.0)
            }
    };
    let ins = [inside(0), inside(1), inside(2), inside(3)];
    // Edges (from lower to higher grid corner): e0 c0-c1, e1 c1-c2, e2 c3-c2, e3 c0-c3.
    let edges: [(usize, usize, EdgeKey); 4] = [(0, 1, (x, y, 0)), (1, 2, (x + 1, y, 1)), (3, 2, (x, y + 1, 0)), (0, 3, (x, y, 1))];
    let cross = |e: usize| {
        let (a, b, _) = edges[e];
        // Both corners 0 (one with a positive neighbour): the middle.
        let d = v[a] - v[b];
        let t = if d == 0.0 { 0.5 } else { (v[a] / d).clamp(0.0, 1.0) };
        [c[a][0] + (c[b][0] - c[a][0]) * t, c[a][1] + (c[b][1] - c[a][1]) * t]
    };
    let crossing: Vec<usize> = (0..4).filter(|&e| ins[edges[e].0] != ins[edges[e].1]).collect();
    // Pair crossing edges into segments.
    let pairs: Vec<(usize, usize)> = match crossing.len() {
        2 => vec![(crossing[0], crossing[1])],
        4 => {
            if (v[4] < 0.0) == ins[0] {
                vec![(0, 1), (2, 3)]
            } else {
                vec![(3, 0), (1, 2)]
            }
        }
        _ => vec![],
    };
    // Orientation is topological, never geometric (segments can have zero
    // length when grid corners lie exactly on the surface). Walking the cell
    // counter-clockwise, a segment runs from the edge where the walk leaves
    // material to the edge where it enters, so material stays on the left.
    let ccw = [(0, 1), (1, 2), (2, 3), (3, 0)];
    let exits = |e: usize| ins[ccw[e].0] && !ins[ccw[e].1];
    let mut out = Vec::new();
    for (ea, eb) in pairs {
        let (from, to) = if exits(ea) { (ea, eb) } else { (eb, ea) };
        out.push(Seg { a: edges[from].2, b: edges[to].2, pa: cross(from), pb: cross(to) });
    }
    out
}

fn join(segs: Vec<Seg>) -> Vec<Loop> {
    let starts: HashMap<EdgeKey, usize> = segs.iter().enumerate().map(|(i, s)| (s.a, i)).collect();
    let mut used = vec![false; segs.len()];
    let mut loops = Vec::new();
    // Walk open chains from their heads first, then the closed loops.
    let ends: std::collections::HashSet<EdgeKey> = segs.iter().map(|s| s.b).collect();
    let heads = (0..segs.len()).filter(|&i| !ends.contains(&segs[i].a));
    let order: Vec<usize> = heads.chain(0..segs.len()).collect();
    for s0 in order {
        if used[s0] {
            continue;
        }
        used[s0] = true;
        let mut pts = vec![segs[s0].pa, segs[s0].pb];
        let mut end = segs[s0].b;
        let mut closed = false;
        while let Some(&n) = starts.get(&end) {
            if n == s0 {
                closed = true;
                pts.pop();
                break;
            }
            if used[n] {
                break;
            }
            used[n] = true;
            pts.push(segs[n].pb);
            end = segs[n].b;
        }
        let m = pts.len();
        let mut area = 0.0;
        let mut perim = 0.0;
        for k in 0..m {
            if !closed && k == m - 1 {
                break;
            }
            let (a, b) = (pts[k], pts[(k + 1) % m]);
            area += a[0] * b[1] - b[0] * a[1];
            perim += ((b[0] - a[0]).powi(2) + (b[1] - a[1]).powi(2)).sqrt();
        }
        loops.push(Loop { pts, closed, area: area / 2.0, perimeter: perim });
    }
    loops.sort_by(|a, b| b.area.abs().total_cmp(&a.area.abs()));
    loops
}

// ---------------------------------------------------------------- clearance

/// Samples a box seeks a witness pair from (those nearest both shapes first).
const WITNESS_STARTS: usize = 8;

/// Every point within d of the box, rounded outward (rounded to nearest,
/// lo - d may land above the true value and drop points within d: a bound
/// over it then certifies nothing).
fn grow(bx: [I; 3], d: f64) -> [I; 3] {
    bx.map(|x| I::new((x.lo - d).next_down(), (x.hi + d).next_up()))
}

/// A tape simplified by the choices an evaluation over a box made (the
/// same where none was made).
fn simple(t: &Arc<Tape>, ch: &[u8]) -> Arc<Tape> {
    if ch.iter().any(|&x| x != BOTH) { Arc::new(t.simplify(ch)) } else { t.clone() }
}

/// A box's centre (halved first, so it cannot overflow; on an unbounded
/// side a finite end, else 0) and the distance from it to the box's far
/// side on each axis, rounded up: the centre may round to one end of a
/// narrow side.
fn centre_reach(bx: [I; 3]) -> (V3, V3) {
    let c: V3 = bx.map(|x| match x.mid() {
        m if m.is_finite() => m,
        _ if x.lo.is_finite() => x.lo,
        _ if x.hi.is_finite() => x.hi,
        _ => 0.0,
    });
    (c, std::array::from_fn(|k| (c[k] - bx[k].lo).max(bx[k].hi - c[k]).next_up()))
}

/// The share of a field's slope bound along axis k (0 to 1; 1 where the
/// slope is unbounded): times the box's half-width there, how far along k
/// it spreads the field's bounds, as a length.
fn slope_share(g: &ink_core::interval::G, l: f64, k: usize) -> f64 {
    let m = ink_core::interval::mag(g.d[k]);
    if !(m > 0.0) {
        0.0
    } else if !(m.is_finite() && l.is_finite() && l > 0.0) {
        1.0
    } else {
        m / l
    }
}

/// A clearance d as a question of the one search: no point of B (the
/// covered shape, `b`) in a box within d of A. `a` is A's field near the
/// region, where witness pairs are sought with `b`; each box keeps A
/// simplified over it grown by d and B over it.
///
/// A box is settled where any test passes, each sound for any fields:
/// - b >= 0 over it (no point of B), by its bounds;
/// - a > 0 over it grown by d, by its bounds;
/// - a >= d L over it, L bounding |grad a| over it grown by d, or the
///   mean-value form a(centre) >= L (r + d), r its half-diagonal;
/// - max(b, a - L d) >= 0 by its first-order bound, whose pairs see B's
///   surface and A's facing each other across the gap, so the boxes need
///   only be about sqrt(margin x radius), where the tests of a alone need
///   boxes as small as the margin;
/// - over the scope's axes (ranged params, time over a step of a motion),
///   a's or b's bound first order in them.
///
/// Open, it proves every point of B in it at least a / L from A (its
/// bound, below d), and its smear on each axis is a's share of its slope
/// bound along it times the half-width: the spread a's bounds take along
/// it, as a length (B's slope would cost a bound per box for no fewer
/// boxes). A ranged param's, and time's, is its blur of a's bound over L,
/// or of b's over b's slope bound where more: a step's moments are halved
/// where the motion blurs a box more than its size does (a face sliding
/// along itself, little). A witness pair is at the box's moment.
struct Clear<'a> {
    a: &'a Arc<Tape>,
    b: &'a Arc<Tape>,
    d: f64,
    /// The region the boxes cover, and its size (the length witness points
    /// are confirmed at).
    region: [I; 3],
    scale: f64,
}

/// A point of B, a point of A within d of it, and their moment and ranged
/// params (the same for both).
type Near = (V3, V3, At);

impl search::Question for Clear<'_> {
    type State = (Arc<Tape>, Arc<Tape>);
    type Witness = Near;

    fn test(&self, n: &search::Node<Self::State>) -> search::Test<Self::State, Near> {
        use ink_core::interval::{add_up, lipschitz, mul_up, norm};
        use search::Test::{Open, Settled};
        let (bx, ctx, d) = (n.b, &*n.ctx, self.d);
        let (sa, sb) = &n.state;
        let mut chb = vec![BOTH; sb.len()];
        if eval_interval(sb, ctx, bx, Some(&mut chb)).lo >= 0.0 {
            return Settled;
        }
        // A field never negative within d of the box has no material there
        // (its zero set has no interior: a polygon's hollow turned about an
        // edge on the axis).
        let grown = grow(bx, d);
        let mut cha = vec![BOTH; sa.len()];
        if eval_interval(sa, ctx, grown, Some(&mut cha)).lo >= 0.0 {
            return Settled;
        }
        let (sa, sb) = (simple(sa, &cha), simple(sb, &chb));
        let (c, hw) = centre_reach(bx);
        let r = norm(hw);
        // The centre values as point intervals (outward-rounded, so sound),
        // first order in the ranged params: the mean-value forms stay tight
        // where a formula has dependencies that inflate plain bounds
        // (helices, rotations).
        let ga = ink_core::eval::eval_interval_grad(&sa, ctx, grown);
        let la = lipschitz(&ga);
        let ac = eval_interval_scope(&sa, ctx, c.map(I::point), None).lo;
        let alo = eval_interval(&sa, ctx, bx, None).lo;
        if la.is_finite() && ac > 0.0 && (d > 0.0 && alo >= mul_up(d, la) || ac >= mul_up(la, add_up(r, d))) {
            return Settled;
        }
        if la.is_finite() && crate::centred::centred_at_least(&[(&sb, 0.0), (&sa, mul_up(la, d))], ctx, bx, 0.0) {
            return Settled;
        }
        // Over the scope's axes (ranged params; time over a step of a
        // motion): bounds first order in them, and how much each blurs
        // them.
        let (mut params, mut t) = (Vec::new(), 0.0);
        let mut ps = scope::splittable(&sa, ctx, scope::PARAM_SPLITS);
        ps.extend(scope::splittable(&sb, ctx, scope::PARAM_SPLITS));
        let ts = scope::spans_t(&sa, ctx) || scope::spans_t(&sb, ctx);
        if !ps.is_empty() || ts {
            let (ra, blur_a) = ink_core::eval::eval_interval_scope_blur(&sa, ctx, grown, None);
            let (rb, blur_b) = ink_core::eval::eval_interval_scope_blur(&sb, ctx, bx, None);
            if ra.lo >= 0.0 || rb.lo >= 0.0 {
                return Settled;
            }
            // Each axis's blur of A's and of B's bounds as a length (over
            // each one's slope bound), the larger.
            let lb = lipschitz(&ink_core::eval::eval_interval_grad(&sb, ctx, bx));
            let (shares, tb) = scope_shares(ctx, &[(&blur_a, la), (&blur_b, lb)]);
            ps.sort_unstable();
            ps.dedup();
            params = ps.iter().map(|&k| (k, share(&shares, k))).collect();
            t = if ts { tb } else { 0.0 };
        }
        // Every point of B here is at least a / L from A, where a > 0 there:
        // a(p) >= alo, or a(centre) - L r, and A's field rises no faster
        // than L within d of the box.
        let bound = if la.is_finite() && la > 0.0 { (alo / la).next_down().max(((ac / la).next_down() - r.next_up()).next_down()) } else { 0.0 };
        let space = std::array::from_fn(|k| match slope_share(&ga, la, k) {
            0.0 => 0.0,
            s => s * hw[k],
        });
        Open { bound: bound.clamp(0.0, d), smear: search::Smear { space, params, t }, state: (sa, sb) }
    }

    /// From samples of B in the box (its 4 x 4 x 4), nearest A first,
    /// Newton steps along grad a toward A's surface, keeping within d;
    /// else, from the samples nearest both shapes, a point of B at or
    /// under its surface and one of A beside it, moved closer by
    /// alternating projections (closer_pair): a gap under d that no sample
    /// lies within d of. Over ranged params, the same with each param
    /// pushed to the end of its span that lowers both fields at the sample
    /// (scope::pushed).
    fn seek(&self, n: &search::Node<Self::State>) -> Option<Near> {
        let (bx, ctx) = (n.b, &*n.ctx);
        let (sa, sb) = &n.state;
        if !bx.iter().all(|x| x.lo.is_finite() && x.hi.is_finite()) {
            return None;
        }
        let pts = box_samples(bx, 4);
        let (mut vb, mut va) = (Vec::new(), Vec::new());
        eval_points(sb, ctx, &pts, &mut vb);
        eval_points(sa, ctx, &pts, &mut va);
        let mut order: Vec<usize> = (0..pts.len()).filter(|&k| pts[k].iter().all(|x| x.is_finite()) && (0..3).all(|j| bx[j].contains(pts[k][j]))).collect();
        // A sample where a >= L d has no point of A within d (L bounds
        // |grad a| on the box grown by d): no walk starts there.
        let reach = ink_core::interval::mul_up(ink_core::interval::lipschitz(&ink_core::eval::eval_interval_grad(sa, ctx, grow(bx, self.d))), self.d);
        order.sort_by(|&i, &j| va[i].total_cmp(&va[j]));
        let walks = order.iter().filter(|&&k| vb[k] < 0.0 && !(va[k] >= reach)).take(WITNESS_STARTS);
        if let Some(w) = walks.copied().find_map(|k| self.walk(&n.state, ctx, pts[k])) {
            return Some(w);
        }
        order.sort_by(|&i, &j| va[i].max(vb[i]).total_cmp(&va[j].max(vb[j])));
        if let Some(w) = order.iter().take(WITNESS_STARTS).find_map(|&k| self.pair(&n.state, ctx, pts[k])) {
            return Some(w);
        }
        if scope::read_params(sa, ctx).is_empty() && scope::read_params(sb, ctx).is_empty() {
            return None;
        }
        order.into_iter().take(4).find_map(|k| {
            let pushed = scope::pushed(&[sa, sb], ctx, pts[k])?.0;
            self.walk(&n.state, &pushed, pts[k]).or_else(|| self.pair(&n.state, &pushed, pts[k]))
        })
    }
}

impl Clear<'_> {
    /// From q in the box, Newton steps toward A's surface along the box's
    /// tape `sa` (`t`: A over the box grown by d, B over the box),
    /// overshooting by a sliver of d so a step lands inside: a point of A
    /// within d of q (so where `sa` is A), confirmed with B's and A's
    /// whole fields by bounds rounded outward?
    fn walk(&self, (sa, _): &(Arc<Tape>, Arc<Tape>), ctx: &Ctx, q: V3) -> Option<Near> {
        let mut p = q;
        for _ in 0..12 {
            let g = ink_core::eval::eval_grads(sa, ctx, &[p])[0];
            if g[0] < 0.0 {
                return self.confirmed(ctx, q, p);
            }
            let n2 = g[1] * g[1] + g[2] * g[2] + g[3] * g[3];
            if !(n2 > 0.0) {
                return None;
            }
            let k = (g[0] + 1e-3 * self.d) / n2;
            p = [p[0] - k * g[1], p[1] - k * g[2], p[2] - k * g[3]];
        }
        None
    }

    /// A pair of points closer than d near x, along the box's tapes: x
    /// onto B, that onto A, moved closer in turn (closer_pair); then each
    /// a quarter of the room left (d less their gap) deeper into its shape
    /// along its slope, so the pair stays within d and inside where the
    /// scope moves the shapes a little (a touches claim over ranged
    /// params). Confirmed with the whole fields (`confirmed`), the deeper
    /// pair first.
    fn pair(&self, (sa, sb): &(Arc<Tape>, Arc<Tape>), ctx: &Ctx, x: V3) -> Option<Near> {
        let onto = |t: &Tape, p: V3| match Confirm::new(ctx).below(t, p, 0.0) {
            true => Some(p),
            false => to_surface(t, ctx, p, 1e-12 * self.scale).and_then(|y| confirm_inside(t, ctx, y, self.scale)),
        };
        let q = onto(sb, x)?;
        let p = onto(sa, q)?;
        let (gap, p, q) = closer_pair(sa, sb, ctx, p, q, self.scale);
        if !(gap < self.d) {
            return None;
        }
        let room = 0.25 * (self.d - gap);
        let deeper = |t: &Tape, x: V3| {
            let g = eval_grads(t, ctx, &[x])[0];
            let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
            [x[0] - room * g[1] / n, x[1] - room * g[2] / n, x[2] - room * g[3] / n]
        };
        self.confirmed(ctx, deeper(sb, q), deeper(sa, p)).or_else(|| self.confirmed(ctx, q, p))
    }

    /// q in the region and B, p in A and within d of q, by B's and A's
    /// whole fields (A's holds within d of the region) and bounds rounded
    /// outward: a witness pair.
    fn confirmed(&self, ctx: &Ctx, q: V3, p: V3) -> Option<Near> {
        let mut confirm = Confirm::new(ctx);
        let ok = q.iter().chain(&p).all(|v| v.is_finite()) && (0..3).all(|k| self.region[k].contains(q[k])) && distance_up(p, q) < self.d;
        (ok && confirm.below(self.b, q, 0.0) && confirm.below(self.a, p, 0.0)).then(|| (q, p, At::of(ctx)))
    }
}

/// Each input's blur of the bounds `blurs` gives (eval_interval_scope_blur:
/// every input of the scope a bound reads), each over its scale (a slope
/// bound, to make it a length; 1, as it is), the largest: each ranged
/// param's, sorted by param, and time's. A derived param's (a solve's
/// unknown) counts for those it is derived from, which are what halve.
fn scope_shares(ctx: &Ctx, blurs: &[(&Vec<(ink_core::graph::Var, f64)>, f64)]) -> (Vec<(u32, f64)>, f64) {
    use ink_core::graph::Var;
    let (mut shares, mut t): (Vec<(u32, f64)>, f64) = (Vec::new(), 0.0);
    for &(blur, l) in blurs {
        for &(v, x) in blur.iter() {
            let x = if x == 0.0 {
                0.0
            } else if l.is_finite() && l > 0.0 {
                x / l
            } else {
                f64::INFINITY
            };
            match v {
                // NaN: as large as any.
                Var::T => t = if x.is_nan() || t.is_nan() { f64::NAN } else { t.max(x) },
                Var::Param(j) if ctx.p_spans.is_derived(j) => shares.extend(ctx.p_spans.sources(j).into_iter().map(|k| (k, x))),
                Var::Param(j) => shares.push((j, x)),
                _ => {}
            }
        }
    }
    shares.sort_unstable_by(|p, q| p.0.cmp(&q.0).then(q.1.total_cmp(&p.1)));
    shares.dedup_by_key(|x| x.0);
    (shares, t)
}

/// Param k's share (`scope_shares`): none where no bound reads it.
fn share(shares: &[(u32, f64)], k: u32) -> f64 {
    shares.binary_search_by_key(&k, |x| x.0).map_or(0.0, |i| shares[i].1)
}

/// `k^3` sample points at the centres of a box's sub-boxes, `k` along each
/// side.
fn box_samples(bx: [I; 3], k: usize) -> Vec<V3> {
    let s: V3 = std::array::from_fn(|j| (bx[j].hi - bx[j].lo) / k as f64);
    let mut out = Vec::with_capacity(k * k * k);
    for a in 0..k {
        for b in 0..k {
            for c in 0..k {
                out.push([bx[0].lo + s[0] * (c as f64 + 0.5), bx[1].lo + s[1] * (b as f64 + 0.5), bx[2].lo + s[2] * (a as f64 + 0.5)]);
            }
        }
    }
    out
}

/// The moment a cached answer's witness is at: the one asked (in_span)
/// where no tape of the question reads t, its cache key leaving the moment
/// out (the answer is the same at any); None where one does: the key holds
/// the moment and its span, and the witness keeps its own (a search over a
/// span finds it at a moment of it).
pub fn asked_moment(ctx: &Ctx, tapes: &[&Tape]) -> Option<f64> {
    (!ink_core::tape::Names::of(tapes, ctx).reads_t).then(|| in_span(ctx).t)
}

/// A verdict's witness at the moment `asked_moment` gives.
pub(crate) fn at_moment(mut v: Verdict, ctx: &Ctx, tapes: &[&Tape]) -> Verdict {
    if let Verdict::Disproven { t, .. } = &mut v
        && let Some(m) = asked_moment(ctx, tapes)
    {
        *t = m;
    }
    v
}

/// One way round for `clearance`: boxes cover `covered` in lo..hi (all of
/// it within d of the other shape; infinite ends: all of space, tested
/// whole), and `grown` is grown by d around each. `swapped`: the covered
/// shape is the claim's first (a pair found is reported the claim's way
/// round).
pub struct Way {
    pub grown: Arc<Tape>,
    pub covered: Arc<Tape>,
    pub lo: V3,
    pub hi: V3,
    pub swapped: bool,
}

/// No point of the second shape within d of the first, exactly: the one
/// search (crate::search) of each way round, its boxes halved down to
/// `finest`, testing at most `budget` boxes between them.
///
/// Either way decides it, and their costs can differ by orders of
/// magnitude: where one shape's field is far from a distance within d of
/// the other (a part cut by a half-space reads 0.1 at points 0.3 away),
/// the boxes covering the other stay undecided at any size. So each round
/// goes to the way with the least volume undecided (the first on a tie),
/// and a way that ends undecided leaves the rest to the other. A way's box
/// stuck at `finest` (the gap there within about it of d, or a field that
/// reads well under the distance near the other shape: a box of slabs
/// beside an edge) ends the way's proof, as best first it is the box
/// nearest a violation. Undecided both ways, both hunt a pair closer than d
/// in the rest (search::Search::hunt), in turn as before, each to its cap
/// or the budget (`hunt`); none found, it gives the way that proved more:
/// the accuracy it holds to (`attained`: every point of the covered
/// shape at least d less that from the other), its boxes left and the box
/// that ended it. The answer is a function of the question; it is cached
/// by it. Returns the verdict, a pair closer than d (the claim's way
/// round), and how finely it looked (search::Search::res).
pub fn clearance(ways: &[Way], ctx: &Ctx, d: f64, finest: f64, budget: usize) -> (Verdict, Option<(V3, V3)>, f64) {
    let tapes: Vec<&Tape> = ways.iter().flat_map(|w| [&*w.grown, &*w.covered]).collect();
    let region = ways.first().map_or(([0.0; 3], [0.0; 3]), |w| (w.lo, w.hi));
    let inputs = |fp: &mut ink_core::graph::Fingerprint| {
        for w in ways {
            w.lo.iter().chain(&w.hi).for_each(|v| _ = fp.f64(*v));
            fp.word(w.swapped as u64);
        }
        fp.f64(d).f64(finest).word(budget as u64);
    };
    let (v, pair, h) = crate::cache::memo_region("clearance", &tapes, ctx, region, inputs, |c| {
        let (v, pair, h, _) = clearance_now(ways, c, d, finest, budget);
        (v, pair, h)
    });
    (at_moment(v, ctx, &tapes), pair, h)
}

/// `clearance`, uncached, and the boxes it tested.
fn clearance_now(ways: &[Way], ctx: &Ctx, d: f64, finest: f64, budget: usize) -> (Verdict, Option<(V3, V3)>, f64, usize) {
    let ctx = Arc::new(in_span(ctx));
    // The grown shape within d of the region (its cells clipped to it,
    // each grown by d): all of it that matters.
    let near: Vec<Arc<Tape>> = ways.iter().map(|w| focus(&w.grown, &ctx, w.lo.map(|x| (x - d).next_down()), w.hi.map(|x| (x + d).next_up()))).collect();
    let qs: Vec<Clear> = ways
        .iter()
        .zip(&near)
        .map(|(w, a)| Clear { a, b: &w.covered, d, region: std::array::from_fn(|k| I::new(w.lo[k], w.hi[k])), scale: max_extent(w.lo, w.hi).min(f64::MAX) })
        .collect();
    let mut searches: Vec<search::Search<(Arc<Tape>, Arc<Tape>), Near>> = ways
        .iter()
        .zip(&qs)
        .map(|(w, q)| {
            let b = std::array::from_fn(|k| I::new(w.lo[k], w.hi[k]));
            search::Search::new(q, vec![search::Node { b, ctx: ctx.clone(), state: (q.a.clone(), w.covered.clone()) }], search::Precision::holds(finest))
        })
        .collect();
    let mut ended: Vec<Option<search::Undecided>> = vec![None; ways.len()];
    let swap = |w: &Way, (q, p): (V3, V3)| if w.swapped { (p, q) } else { (q, p) };
    let tested = |s: &[search::Search<(Arc<Tape>, Arc<Tape>), Near>]| s.iter().map(|s| s.tested()).sum::<usize>();
    for hunt in [false, true] {
        // Neither way proves it: unless stopped, both hunt a pair closer
        // than d (search::Search::hunt), in turn as before.
        if hunt {
            if ended.iter().flatten().any(|u| u.stopped()) {
                break;
            }
            searches.iter_mut().for_each(|s| s.hunt());
            ended.iter_mut().for_each(|u| *u = None);
        }
        while let Some(k) = (0..ways.len()).filter(|&k| ended[k].is_none()).min_by(|&i, &j| searches[i].volume().total_cmp(&searches[j].volume())) {
            let others: usize = searches.iter().enumerate().filter(|x| x.0 != k).map(|x| x.1.tested()).sum();
            match searches[k].round(&qs[k], budget.saturating_sub(others)) {
                None => {}
                Some(search::Outcome::Proven) => return (Verdict::Proven, None, searches[k].res(), tested(&searches)),
                Some(search::Outcome::Found((q, p, at))) => {
                    let v = Verdict::Disproven { witness: q, t: at.t, value: distance_up(p, q), params: at.params };
                    return (v, Some(swap(&ways[k], (q, p))), searches[k].res(), tested(&searches));
                }
                Some(search::Outcome::Unknown(u)) => ended[k] = Some(u),
            }
        }
    }
    // Both ways' hunts: their cells, stopped by the time limit if either
    // was, by the refinement limit if both were, else by a cap or the budget.
    let hunts: Vec<search::Hunt> = ended.iter().flatten().filter_map(|u| u.hunt).collect();
    let hunt = (!hunts.is_empty()).then(|| {
        use search::HuntEnd::{Budget, Cap, Stuck, Time};
        let by = |e| hunts.iter().any(|h| h.stopped == e);
        let stopped = if by(Time) {
            Time
        } else if hunts.iter().all(|h| h.stopped == Stuck) {
            Stuck
        } else if by(Cap) {
            Cap
        } else {
            Budget
        };
        search::Hunt { cells: hunts.iter().map(|h| h.cells).sum(), stopped }
    });
    // The way that proved the larger clearance (a stop ends both).
    let Some((k, u)) =
        ended.into_iter().enumerate().filter_map(|(k, u)| Some((k, u?))).max_by(|x, y| x.1.stopped().cmp(&y.1.stopped()).then(x.1.bound.total_cmp(&y.1.bound)).then(y.0.cmp(&x.0)))
    else {
        return (Verdict::Unknown { cells: 0, volume: 0.0, ended: Ended::Finest, cell: None, attained: None, hunt: None }, None, finest, tested(&searches));
    };
    let attained = (u.bound > 0.0).then(|| (d - u.bound).next_up());
    (Verdict::Unknown { cells: u.cells, volume: u.volume, ended: u.ended, cell: u.cell, attained, hunt }, None, searches[k].res(), tested(&searches))
}

/// The minimum distance between two shapes, bracketed (`min_distance`):
/// every point of each at least `lower` from the other, proven; a pair
/// confirmed inside them (on the first, on the second, the claim's way
/// round) and its distance rounded up, `upper`; why the search ended short
/// of the accuracy (None: it reached it); and how finely it looked.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct MinDistance {
    pub lower: f64,
    pub upper: Option<(f64, V3, V3)>,
    pub ended: Option<Ended>,
    pub res: f64,
    /// Boxes tested.
    pub tested: usize,
}

/// The minimum distance between the shapes of `ways` (`clearance`'s, each
/// way round), to `tol`: one search, clearance's question asked of a d
/// that shrinks. It starts at `far` (no pair is farther apart); each pair
/// a box's seek finds closer than d, moved closer by alternating
/// projections (closer_pair), is the upper bound, and d becomes it less
/// `tol`; it starts from `seed`'s where given (a pair confirmed inside
/// both, the claim's way round: (distance, on the first, on the second)). The boxes already tested keep what they proved (a box settled at
/// a d is settled at any smaller one, and one whose bound reaches the new d
/// is settled by it: Search::settle_from), so the search goes on where it
/// was. Proven at d, the distance is within [d, upper]; ended undecided,
/// `lower` is the bound the way that proved more attained. Boxes halved to
/// `finest`, at most `budget` tested. Cached by the question.
pub fn min_distance(ways: &[Way], ctx: &Ctx, far: f64, seed: Option<(f64, V3, V3)>, tol: f64, finest: f64, budget: usize) -> MinDistance {
    let tapes: Vec<&Tape> = ways.iter().flat_map(|w| [&*w.grown, &*w.covered]).collect();
    let region = ways.first().map_or(([0.0; 3], [0.0; 3]), |w| (w.lo, w.hi));
    let inputs = |fp: &mut ink_core::graph::Fingerprint| {
        for w in ways {
            w.lo.iter().chain(&w.hi).for_each(|v| _ = fp.f64(*v));
            fp.word(w.swapped as u64);
        }
        fp.f64(far).f64(tol).f64(finest).word(budget as u64);
        if let Some((g, p, q)) = seed {
            [g].iter().chain(&p).chain(&q).for_each(|v| _ = fp.f64(*v));
        }
    };
    crate::cache::memo_region("min distance", &tapes, ctx, region, inputs, |c| min_distance_now(ways, c, far, seed, tol, finest, budget))
}

fn min_distance_now(ways: &[Way], ctx: &Ctx, far: f64, seed: Option<(f64, V3, V3)>, tol: f64, finest: f64, budget: usize) -> MinDistance {
    let ctx = Arc::new(in_span(ctx));
    let mut upper = seed;
    let mut d = seed.map_or(far, |u| u.0 - tol);
    if !(d > 0.0) {
        return MinDistance { lower: 0.0, upper, ended: None, res: f64::INFINITY, tested: 0 };
    }
    let near: Vec<Arc<Tape>> = ways.iter().map(|w| focus(&w.grown, &ctx, w.lo.map(|x| (x - d).next_down()), w.hi.map(|x| (x + d).next_up()))).collect();
    let region = |w: &Way| std::array::from_fn(|k| I::new(w.lo[k], w.hi[k]));
    let scale = |w: &Way| max_extent(w.lo, w.hi).min(f64::MAX);
    let ask = |k: usize, d: f64| Clear { a: &near[k], b: &ways[k].covered, d, region: region(&ways[k]), scale: scale(&ways[k]) };
    let mut searches: Vec<search::Search<(Arc<Tape>, Arc<Tape>), Near>> = ways
        .iter()
        .enumerate()
        .map(|(k, w)| {
            search::Search::new(&ask(k, d), vec![search::Node { b: region(w), ctx: ctx.clone(), state: (near[k].clone(), w.covered.clone()) }], search::Precision::holds(finest))
        })
        .collect();
    let mut ended: Vec<Option<search::Undecided>> = vec![None; ways.len()];
    let tested = |s: &[search::Search<(Arc<Tape>, Arc<Tape>), Near>]| s.iter().map(|s| s.tested()).sum::<usize>();
    let res = |s: &[search::Search<(Arc<Tape>, Arc<Tape>), Near>]| s.iter().map(|s| s.res()).fold(f64::INFINITY, f64::min);
    while let Some(k) = (0..ways.len()).filter(|&k| ended[k].is_none()).min_by(|&i, &j| searches[i].volume().total_cmp(&searches[j].volume())) {
        let others: usize = searches.iter().enumerate().filter(|x| x.0 != k).map(|x| x.1.tested()).sum();
        match searches[k].round(&ask(k, d), budget.saturating_sub(others)) {
            None => {}
            Some(search::Outcome::Proven) => return MinDistance { lower: d.max(0.0), upper, ended: None, res: res(&searches), tested: tested(&searches) },
            Some(search::Outcome::Found((q, p, _))) => {
                let w = &ways[k];
                let (gap, p, q) = closer_pair(&w.grown, &w.covered, &ctx, p, q, scale(w));
                // The claim's way round: on the first, on the second.
                let pair = if w.swapped { (gap, q, p) } else { (gap, p, q) };
                if upper.is_none_or(|u| gap < u.0) {
                    upper = Some(pair);
                }
                d = gap - tol;
                if !(d > 0.0) {
                    return MinDistance { lower: 0.0, upper, ended: None, res: res(&searches), tested: tested(&searches) };
                }
                searches.iter_mut().for_each(|s| s.settle_from(d));
            }
            Some(search::Outcome::Unknown(u)) => ended[k] = Some(u),
        }
    }
    // The way that proved more (a stop ends both).
    let u = ended.into_iter().flatten().max_by(|x, y| x.stopped().cmp(&y.stopped()).then(x.bound.total_cmp(&y.bound)));
    let lower = u.as_ref().map_or(0.0, |u| u.bound.clamp(0.0, d.max(0.0)));
    MinDistance { lower, upper, ended: Some(u.map_or(Ended::Finest, |u| u.ended)), res: res(&searches), tested: tested(&searches) }
}

// ---------------------------------------------------------------- rays

/// A side of a field's zero set: inside (the field < 0) or outside (>= 0).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    Inside,
    Outside,
}

impl Side {
    pub fn holds(self, v: f64) -> bool {
        match self {
            Side::Inside => v < 0.0,
            Side::Outside => v >= 0.0,
        }
    }
}

/// Where a ray reached a side of a field (ray_to): at t, and, where
/// bounds rounded outward do not show a point within its eps after t on
/// that side, how far from 0 the field's real value may be there
/// (`rounding`: the ray runs in the surface, so only rounding decides the
/// side, and it may be on either).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RayHit {
    pub t: f64,
    pub rounding: Option<f64>,
}

/// The first t in [t0, t1] where o + t dir is on `side` of the field.
/// Interval arithmetic over segments of the ray, nearest first, each with
/// the tape `local` keeps for the smallest block holding it, rules out the
/// stretches wholly on the other side, as does the field's value at a
/// segment's middle (bounded outward) and its rate along the ray (a mean
/// value bound, finite only where the field is continuous: no undecided
/// select, floor or atan2 across a jump); a segment neither rules out is
/// halved down to `eps`, where the values at its ends decide. So a stretch
/// of `side` longer than eps is never stepped over, and a touch that does
/// not cross is not a crossing. A segment whose middle's bound holds 0 and
/// whose rate times its half-length is under that bound's width is a ray
/// lying in the surface: halving it decides nothing, so it ends the search,
/// a hit where it crosses if bounds show that one sure, else from its
/// start, within rounding of 0 there (RayHit::rounding). The crossing found
/// is bisected to the precision of floats: t is on `side`, in floats, and
/// a hit is sure where bounds show a point of it within eps after t.
/// None: no point of the ray is on that side, but for stretches shorter
/// than eps.
pub fn ray_to(local: &Local, o: V3, dir: V3, (t0, t1): (f64, f64), side: Side, eps: f64) -> Option<RayHit> {
    let ctx = local.ctx();
    let at = |t: f64| [o[0] + t * dir[0], o[1] + t * dir[1], o[2] + t * dir[2]];
    let f = |t: f64| {
        let p = at(t);
        eval_point(&local.tape_over(&p.map(I::point)), ctx, p)
    };
    let eps = eps.max(1e-9 * t1.abs());
    // Bounds `r` show: no point on `side`, every point on it.
    let decide = |r: I| match side {
        Side::Inside => (r.lo >= 0.0, r.hi < 0.0),
        Side::Outside => (r.hi < 0.0, r.lo >= 0.0),
    };
    // The bounds at t: of the ray's real point there, and its float.
    let point = |t: f64| {
        let (p, b) = (at(t), segment_box(o, dir, t, t));
        let b = std::array::from_fn(|k| I::new(b[k].lo.min(p[k]), b[k].hi.max(p[k])));
        eval_interval(&local.tape_over(&b), ctx, b, None)
    };
    // t, sure where bounds show it, or a point up to eps after it and
    // within the ray, on `side` (2^-20 eps first).
    let hit = |t: f64| {
        let r = point(t);
        let sure = decide(r).1 || (0..=20).rev().map(|i| t + eps * 0.5f64.powi(i)).filter(|&s| s <= t1).any(|s| decide(point(s)).1);
        RayHit { t, rounding: (!sure).then(|| r.lo.abs().max(r.hi.abs())) }
    };
    let start = f(t0);
    if side.holds(start) {
        return Some(hit(t0));
    }
    // Segments from t0, each twice as long as the one before, so a ray that
    // crosses soon is searched in the small blocks near where it starts;
    // each depth-first, nearest half first. The first is as long as the
    // field's value at t0 (a distance field's gap to the surface; a field
    // scaled by k costs log2 k more segments), at least eps.
    let (mut from, mut len) = (t0, if start.is_finite() { start.abs().max(eps) } else { eps });
    let (mut stack, mut slots) = (Vec::new(), Vec::new());
    loop {
        let (a, b) = match stack.pop() {
            Some(s) => s,
            None if from < t1 => {
                let s = (from, (from + len).min(t1));
                (from, len) = (s.1, 2.0 * len);
                s
            }
            None => return None,
        };
        let bx = segment_box(o, dir, a, b);
        let tape = local.tape_over(&bx);
        let r = eval_interval_in(&tape, ctx, bx, None, &mut slots);
        let (mut none, mut all) = decide(r);
        if !none && !all {
            // Its middle's bound and its rate along the ray: the real field
            // over the segment lies within their mean value bound.
            let g = ink_core::eval::eval_interval_grad(&tape, ctx, bx);
            let rate = (0..3).fold(I::point(0.0), |s, k| add_o(s, mul_o(g.d[k], I::point(dir[k]))));
            let m = span_mid((a, b));
            let mid = eval_interval_in(&tape, ctx, segment_box(o, dir, m, m), None, &mut slots);
            if rate.lo.is_finite() && rate.hi.is_finite() && mid.lo.is_finite() && mid.hi.is_finite() {
                let e = add_o(mid, mul_o(rate, sub_o(I::new(a, b), I::point(m))));
                (none, all) = decide(e);
                let undecided = !(decide(mid).0 || decide(mid).1);
                if !none && !all && undecided && ink_core::interval::mag(rate) * (b - a) * 0.5 <= mid.hi - mid.lo {
                    // Where it crosses in floats, if it does, and sure there
                    // where bounds show it (`hit`); else from a, to rounding.
                    let (fa, fb) = (f(a), f(b));
                    let t = match (side.holds(fa), side.holds(fb)) {
                        (false, true) => crossing(&f, side, (a, fa), (b, fb)),
                        _ => a,
                    };
                    return Some(match hit(t) {
                        h @ RayHit { rounding: None, .. } => h,
                        _ => RayHit { t: a, rounding: Some(e.lo.abs().max(e.hi.abs())) },
                    });
                }
            }
        }
        if none {
            continue;
        }
        if b - a <= eps.max(b.abs() * 1e-12) || all {
            // On that side somewhere in [a, b]: find the first crossing.
            let (fa, fb) = (f(a), f(b));
            if side.holds(fa) {
                return Some(hit(a));
            }
            if !side.holds(fb) {
                // No crossing at the ends; refine this stretch further.
                if b - a > eps {
                    let m = span_mid((a, b));
                    stack.push((m, b));
                    stack.push((a, m));
                }
                continue;
            }
            return Some(hit(crossing(&f, side, (a, fa), (b, fb))));
        }
        let m = span_mid((a, b));
        stack.push((m, b));
        stack.push((a, m));
    }
}

/// The crossing between a (off `side`) and b (on it), to the precision of
/// floats: the end on `side`. False position, Illinois's way (an end kept
/// twice has its value halved), and a bisection after any step that leaves
/// more than half the bracket, so it takes at most about twice the steps of
/// bisection and, where the field is smooth, a handful.
fn crossing(f: &dyn Fn(f64) -> f64, side: Side, (mut a, mut fa): (f64, f64), (mut b, mut fb): (f64, f64)) -> f64 {
    let (mut kept, mut halve) = (0i8, false);
    loop {
        let (w, mid) = (b - a, span_mid((a, b)));
        if mid <= a || mid >= b {
            break;
        }
        let m = if halve { mid } else { a + w * (fa / (fa - fb)) };
        let m = if m > a && m < b { m } else { mid };
        let fm = f(m);
        if side.holds(fm) {
            (b, fb) = (m, fm);
            if kept == -1 {
                fa *= 0.5;
            }
            kept = -1;
        } else {
            (a, fa) = (m, fm);
            if kept == 1 {
                fb *= 0.5;
            }
            kept = 1;
        }
        halve = b - a > 0.5 * w;
    }
    b
}

/// A box holding the points o + t dir for t in [t0, t1], rounded outward
/// (nearest ends can miss a sliver of the segment, and the crossing in it).
fn segment_box(o: V3, dir: V3, t0: f64, t1: f64) -> [I; 3] {
    std::array::from_fn(|k| add_o(I::point(o[k]), mul_o(I::new(t0, t1), I::point(dir[k]))))
}

// ---------------------------------------------------------------- integrals

/// Integrals of each of `gs` (None: 1) over the region shape < 0 at cells
/// of `h`; with `moments`, the first's times [1, x, y, z, xx, yy, zz, xy,
/// yz, zx] in its place (mass moments, rho the first). Cells full of
/// material use their 2x2x2 Gauss points, exact where an integrand is
/// linear across them, halved where the first is not and its integral
/// there is not bounded (`full_integral`) to ACCURACY of itself, down to
/// cells of one grid step: points alone miss a layer of density between
/// them at every cell size; bounds see it. Cells the surface may cross use
/// samples weighted by the material in them (`occupancy`). Every integrand
/// at the same points, in one descent.
pub fn sums(shape: &Arc<Tape>, gs: &[Option<&Tape>], moments: bool, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Sums {
    let (grid, root) = root_for(lo, hi, h, 1);
    let n = gs.len() + if moments && !gs.is_empty() { 9 } else { 0 };
    // Each integral at points weighted `w`.
    let at = |pts: &[V3], w: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; n];
        let mut k = 0;
        for (i, g) in gs.iter().enumerate() {
            let mut v = Vec::new();
            match g {
                Some(t) => eval_points(t, ctx, pts, &mut v),
                None => v.resize(pts.len(), 1.0),
            }
            let terms: &[fn(V3) -> f64] = if i == 0 && moments { &MOMENTS } else { &MOMENTS[..1] };
            for term in terms {
                out[k] = v.iter().zip(pts).zip(w).map(|((a, p), b)| a * term(*p) * b).sum();
                k += 1;
            }
        }
        out
    };
    let gauss = |c: Cell| -> Vec<f64> {
        let (a, b) = grid.cell_bounds(c);
        let g = 0.5 / 3f64.sqrt();
        let pts: Vec<V3> = (0..8).map(|k| std::array::from_fn(|ax| 0.5 * (a[ax] + b[ax]) + if (k >> ax) & 1 == 1 { g } else { -g } * (b[ax] - a[ax]))).collect();
        at(&pts, &[(grid.h * c.size as f64).powi(3) / 8.0; 8])
    };
    fn full(first: Option<&Tape>, ctx: &Ctx, grid: &Grid, c: Cell, gauss: &(dyn Fn(Cell) -> Vec<f64> + Sync)) -> ([usize; 2], Vec<f64>) {
        if let Some(r) = first
            && c.size > 1
        {
            let (b, linear) = full_integral(r, ctx, grid, c, c.size);
            if !linear && !(b.hi - b.lo <= ACCURACY * ink_core::interval::mag(b)) {
                let halves: Vec<([usize; 2], Vec<f64>)> = c.children().par_iter().map(|&k| full(first, ctx, grid, k, gauss)).collect();
                return add_sums(halves);
            }
        }
        ([1, 0], gauss(c))
    }
    let first = gs.first().copied().flatten();
    let parts: Vec<([usize; 2], Vec<f64>)> = descend(
        shape,
        ctx,
        &grid,
        root,
        1,
        &|t: &Arc<Tape>, c| {
            let pts = samples(&grid, c, 2);
            let fg = eval_grads(t, ctx, &pts);
            let vol = cell_volume(&grid, c) / pts.len() as f64;
            let s = grid.cell_sides(c).map(|x| x / 2.0);
            let w: Vec<f64> = fg.iter().map(|f| occupancy(f, s) * vol).collect();
            Some(([0, 1], at(&pts, &w)))
        },
        &|_: &Arc<Tape>, c| Some(full(first, ctx, &grid, c, &gauss)),
    );
    let (cells, mut sums) = add_sums(parts);
    sums.resize(n, 0.0);
    Sums { sums, cells }
}

/// The terms of mass moments: [1, x, y, z, xx, yy, zz, xy, yz, zx].
const MOMENTS: [fn(V3) -> f64; 10] = [|_| 1.0, |p| p[0], |p| p[1], |p| p[2], |p| p[0] * p[0], |p| p[1] * p[1], |p| p[2] * p[2], |p| p[0] * p[1], |p| p[1] * p[2], |p| p[2] * p[0]];

/// `sums`' parts, summed in order.
fn add_sums(parts: Vec<([usize; 2], Vec<f64>)>) -> ([usize; 2], Vec<f64>) {
    parts.into_iter().fold(([0; 2], Vec::new()), |(n, s), (m, x)| {
        let s = if s.is_empty() { x } else { s.iter().zip(&x).map(|(a, b)| a + b).collect() };
        ([n[0] + m[0], n[1] + m[1]], s)
    })
}

/// Eigenvalues and eigenvectors (columns) of a symmetric 3x3 matrix (Jacobi).
pub fn sym_eigen(m: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let (mut p, mut q, mut big) = (0, 1, 0.0);
        for i in 0..3 {
            for j in i + 1..3 {
                if a[i][j].abs() > big {
                    (p, q, big) = (i, j, a[i][j].abs());
                }
            }
        }
        if big < 1e-15 * (a[0][0].abs() + a[1][1].abs() + a[2][2].abs()).max(1e-300) {
            break;
        }
        let theta = 0.5 * libm::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
        let (c, s) = (libm::cos(theta), libm::sin(theta));
        for k in 0..3 {
            let (akp, akq) = (a[k][p], a[k][q]);
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
        }
        for k in 0..3 {
            let (apk, aqk) = (a[p][k], a[q][k]);
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
        }
        for k in 0..3 {
            let (vkp, vkq) = (v[k][p], v[k][q]);
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::graph::{Graph, Var};

    /// A field rewritten (x to p, y to
    /// a name of p) and then untagged is p - p = 0 and 1e16 + 0.5 above it
    /// in real numbers: +0.25 everywhere, the probe outside. Untagging is
    /// `Graph::rebuild`, like every rewrite: float folds of its own would
    /// fold the marked arithmetic to -0.25 and prove it inside (and its
    /// negation empty), exactly.
    #[test]
    fn an_untagged_rewrite_keeps_spatial_arithmetic_real() {
        let ctx = Ctx::new(0.0, vec![0.5]);
        for negated in [false, true] {
            let mut g = Graph::new();
            let (x, y, p) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Param(0)));
            let (half, big, quarter) = (g.constant(0.5), g.constant(1e16), g.constant(0.25));
            let d = g.sub(x, y);
            let a = g.add(d, half);
            let b = g.add(big, a);
            let c = g.sub(b, big);
            let f = g.sub(c, quarter);
            let f = if negated { g.neg(f) } else { f };
            let named = g.tag("same-input", p);
            let rewritten = g.subst(f, &[(Var::X, p), (Var::Y, named)]);
            let untagged = g.untagged(rewritten);
            let tape = Arc::new(Tape::new(&g, untagged));
            let v = eval_interval(&tape, &ctx, [I::new(-5.0, 5.0); 3], None);
            let real = if negated { -0.25 } else { 0.25 };
            assert!(v.lo <= real && real <= v.hi, "untagged, {} is {real}, bounded [{}, {}]", if negated { "its negation" } else { "the field" }, v.lo, v.hi);
            let proof = prove_empty(&tape, &ctx, [-5.0; 3], [5.0; 3], Some(0.1), 0.05, 0.05);
            let wrong = if negated { matches!(proof.verdict, Verdict::Proven) } else { matches!(proof.verdict, Verdict::Disproven { .. }) };
            assert!(!wrong, "untagged, {} is decided wrong: {:?}", if negated { "its negation" } else { "the field" }, proof.verdict);
        }
    }

    /// The same question with its param
    /// numbered 0 or 71 (a part imported after others) has one cache key,
    /// marked real or not.
    #[test]
    fn a_real_question_keeps_its_key_whatever_its_params_number() {
        let opts = crate::ranged::ProofOptions { res: Some(0.1), finest: 0.05, depth: 12, tol: 0.05 };
        let key = |k: u32, real: bool| {
            let mut g = Graph::new();
            let (x, p) = (g.var(Var::X), g.var(Var::Param(k)));
            let (big, quarter) = (g.constant(1e16), g.constant(0.25));
            let on = if real { x } else { p };
            let s = g.add(big, on);
            let s = g.sub(s, big);
            let f = g.sub(s, quarter);
            let f = if real { g.subst(f, &[(Var::X, p)]) } else { f };
            assert_eq!(g.reachable(&[f]).iter().any(|&n| matches!(g.nodes[n as usize], ink_core::graph::Node::Real(_))), real);
            let mut params = vec![0.0; k as usize + 1];
            params[k as usize] = 0.5;
            crate::ranged::prove_key(&Tape::new(&g, f), &Ctx::new(0.0, params), [-5.0; 3], [5.0; 3], &opts)
        };
        for real in [true, false] {
            let (a, b) = (key(0, real), key(71, real));
            assert!(a.is_some() && a == b, "real {real}: keys {a:?}, {b:?}");
        }
        assert!(key(0, true) != key(0, false), "a real question and its float are different questions");
    }

    #[test]
    fn a_density_varying_inside_full_cells_is_measured() {
        // A 10 mm cube of 0.001 g/mm^3 with a 0.2 mm layer of 0.011 across
        // x = 0: 1.2 g. Its full cells' 8 points each miss the layer at
        // every cell size (1.084 g, settled); their bounds see it.
        use crate::test_fields::*;
        use ink_core::graph::{BinOp, UnOp};
        let ctx = at(0.0);
        let mut b = B::new();
        let cube = b.cuboid([0.0; 3], [5.0; 3]);
        let tape = b.clip(cube);
        let x = b.v(Var::X);
        let ax = b.un(UnOp::Abs, x);
        let edge = b.c(0.1);
        let side = b.bin(BinOp::Sub, ax, edge);
        let (dense, none) = (b.c(0.01), b.c(0.0));
        let layer = b.0.select(side, dense, none);
        let base = b.c(0.001);
        let rho = b.bin(BinOp::Add, base, layer);
        let rho = b.tape(rho);
        let (lo, hi, _) = bounds(&tape, &ctx, LO, HI).unwrap();
        let m = measured(&tape, &[Some(&rho), None], true, &ctx, lo, hi, None);
        assert!((m.sums[0] / 1.2 - 1.0).abs() <= ACCURACY, "{m:?}");
        // At a cell size asked for, the layer is resolved to it.
        let m = measured(&tape, &[Some(&rho), None], true, &ctx, lo, hi, Some(0.025));
        assert!((m.sums[0] / 1.2 - 1.0).abs() <= ACCURACY, "{m:?}");
    }

    #[test]
    fn sums_with_an_integrand_in_another_slot_are_another_question() {
        // 65 integrands, one a tape and the rest None (1), the tape at slot
        // 0 or at slot 64: the presence of each slot once hashed into one
        // word, bit k % 64, so the two shared a key and one's cached sums
        // answered the other.
        let mut g = Graph::new();
        let (inside, two) = (g.constant(-1.0), g.constant(2.0));
        let shape = Arc::new(Tape::new(&g, inside));
        let two = Tape::new(&g, two);
        let ctx = Ctx::new(0.0, vec![]);
        let sums = |at: usize| {
            let mut gs = vec![None; 65];
            gs[at] = Some(&two);
            measured(&shape, &gs, false, &ctx, [-1.0; 3], [1.0; 3], Some(0.5)).sums
        };
        let (a, b) = (sums(0), sums(64));
        let near = |x: f64, y: f64| (x - y).abs() <= 1e-9 * y;
        assert!(near(a[0], 16.0) && near(a[64], 8.0), "{a:?}");
        assert!(near(b[0], 8.0) && near(b[64], 16.0), "{b:?}");
    }

    #[test]
    fn bounds_visit_extrema_without_walking_the_whole_surface() {
        use ink_core::graph::BinOp;
        let mut b = crate::test_fields::B::new();
        let cube = b.cuboid([0.0; 3], [10.0; 3]);
        // A disconnected component must still enlarge the box, even
        // after the first component has supplied all six extrema.
        let remote = b.cuboid([30.0, 2.0, -1.0], [1.0; 3]);
        let both = b.bin(BinOp::Min, cube, remote);
        let ctx = Ctx::new(0.0, vec![]);
        for f in [cube, both] {
            let tape = Arc::new(b.tape(f));
            let (grid, root) = try_root_for([-40.0; 3], [40.0; 3], 0.25, 1).unwrap();
            let (found, work) = bounds_cells(&tape, &ctx, &grid, root);
            let leaves: Vec<(V3, V3)> = descend(&tape, &ctx, &grid, root, 1, &|_, c| Some(grid.cell_bounds(c)), &|_, c| Some(grid.cell_bounds(c)));
            assert_eq!(found, leaves.into_iter().fold(ink_core::v3::EMPTY_BOX, ink_core::v3::box_union));
            assert!(work < 10_000, "bounding a flat face must not visit all its leaves: {work} cells");
        }
    }

    #[test]
    fn a_witness_over_a_span_is_at_a_moment_of_its_own() {
        // Two cubes of side 2, the second at x = 10 t: over t in [0.1, 0.4]
        // they overlap until 0.2 and are apart after, at the span's middle
        // (0.25) too. The search over the span halves its moments and finds
        // the overlap at one before 0.2, which the witness carries (inside
        // both there), through a proof's cached answer too. With time not
        // halved it would be undecided over the span; one reported at the
        // span's middle would be outside them.
        let mut b = crate::test_fields::B::new();
        let a = b.cuboid([0.0; 3], [1.0; 3]);
        let (x, y, z, t) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z), b.v(Var::T));
        let ten = b.c(10.0);
        let moved = b.bin(ink_core::graph::BinOp::Mul, ten, t);
        let x = b.bin(ink_core::graph::BinOp::Sub, x, moved);
        let (sx, sy, sz) = (b.slab(x, 0.0, 1.0), b.slab(y, 0.0, 1.0), b.slab(z, 0.0, 1.0));
        let m = b.max3(sx, sy, sz);
        let both = b.bin(ink_core::graph::BinOp::Max, a, m);
        let tape = Arc::new(b.tape(both));
        let (lo, hi) = ([-2.0; 3], [6.0, 2.0, 2.0]);
        let ctx = Ctx { t_span: Some((0.1, 0.4)), ..Ctx::new(0.25, vec![]) };
        let inside = |t: f64, w: V3| (0.1..0.2).contains(&t) && eval_point(&tape, &Ctx::new(t, vec![]), w) < 0.0;
        let s = prove_empty(&tape, &ctx, lo, hi, None, 0.05, 0.0);
        assert!(matches!(s.verdict, Verdict::Disproven { witness, t, .. } if inside(t, witness)), "{s:?}");
        let opts = crate::ranged::ProofOptions { res: None, finest: 0.05, depth: 8, tol: 0.0 };
        let p = crate::ranged::prove(&tape, &ctx, lo, hi, &opts);
        assert!(matches!(p.verdict, crate::query::Verdict::Disproven { witness, t, .. } if inside(t, witness)), "{p:?}");
    }

    /// A witness at a corner of the params is found in a few thousand
    /// boxes: a block at the end of a stack of three lengths a, b, c (each
    /// 10 +-0.1) against a wall at g (30.2 +-0.1), which overlap by up to
    /// 0.2 only where all four reach their worst ends (a cube of side 0.1,
    /// twice the tolerance, fits only within 0.1 of that corner). Each box
    /// tries its samples with every param pushed to the end of its span
    /// that lowers the field there (scope::pushed), so the corner is found
    /// once a box's samples fall in the overlap, whatever order the
    /// params' halves are searched in. The passes this replaced halved the
    /// lower halves first and searched each to its finest splits before a
    /// witness after it could count (benchmark 35: 22 s, millions of
    /// cells). Nominal (no ranges), it holds.
    #[test]
    fn a_witness_at_a_corner_of_the_params_is_found_in_few_boxes() {
        use crate::test_fields::B;
        use ink_core::graph::BinOp;
        let mut b = B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (pa, pb, pc, pg) = (b.v(Var::Param(0)), b.v(Var::Param(1)), b.v(Var::Param(2)), b.v(Var::Param(3)));
        let ab = b.bin(BinOp::Add, pa, pb);
        let end = b.bin(BinOp::Add, ab, pc);
        let past = b.bin(BinOp::Sub, x, end);
        let before = b.bin(BinOp::Sub, pg, x);
        let (sy, sz) = (b.slab(y, 0.0, 5.0), b.slab(z, 0.0, 5.0));
        let f = b.max3(past, before, sy);
        let f = b.bin(BinOp::Max, f, sz);
        let tape = Arc::new(b.tape(f));
        let band = |v: f64| (v - 0.1, v + 0.1);
        let ranged =
            Ctx { p_spans: ink_core::eval::ParamSpans::new([(0, band(10.0)), (1, band(10.0)), (2, band(10.0)), (3, band(30.2))]), ..Ctx::new(0.0, vec![10.0, 10.0, 10.0, 30.2]) };
        let (lo, hi) = ([29.0, -6.0, -6.0], [31.0, 6.0, 6.0]);
        let (s, n) = empty_now(&tape, &ranged, lo, hi, None, 0.025, 0.05);
        let Verdict::Disproven { witness, params, .. } = &s.verdict else { panic!("{s:?}") };
        let at = |k: u32| params.iter().find(|p| p.0 == k).unwrap().1;
        assert!(at(0) + at(1) + at(2) - at(3) >= 0.1 && (at(3)..=at(0) + at(1) + at(2)).contains(&witness[0]) && n < 10_000, "{s:?} after {n} boxes");
        let nominal = Ctx::new(0.0, vec![10.0, 10.0, 10.0, 30.2]);
        assert!(matches!(prove_empty(&tape, &nominal, lo, hi, None, 0.025, 0.05).verdict, Verdict::Proven));
    }

    /// A search that would test far more boxes than it may ends in a few
    /// thousand: a ball exactly inside an envelope of its size written
    /// otherwise ((r^2 - 400) / 40, the same sphere), asked to 0.0001,
    /// whose curved faces a pair's first-order bound cancels only to the
    /// square of a box's size, so every box along them settles only at
    /// about 0.02 across: millions of boxes. The passes this replaced
    /// stopped at their pass limit in a second; the one search tested its
    /// four million boxes first (25 s), until it ends where the boxes it
    /// would test, at the rate each depth has opened them, pass four times
    /// its budget.
    #[test]
    fn a_search_that_will_not_fit_ends_early() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        let mut b = B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let sq = [x, y, z].map(|v| b.un(UnOp::Square, v));
        let rr = b.bin(BinOp::Add, sq[0], sq[1]);
        let rr = b.bin(BinOp::Add, rr, sq[2]);
        let (r, size, area, scale) = (b.un(UnOp::Sqrt, rr), b.c(20.0), b.c(400.0), b.c(40.0));
        let ball = b.bin(BinOp::Sub, r, size);
        let envelope = b.bin(BinOp::Sub, rr, area);
        let envelope = b.bin(BinOp::Div, envelope, scale);
        let outside = b.un(UnOp::Neg, envelope);
        let f = b.bin(BinOp::Max, ball, outside);
        let tape = Arc::new(b.tape(f));
        let (s, n) = empty_now(&tape, &Ctx::new(0.0, vec![]), [-21.0; 3], [21.0; 3], None, 0.00005, 0.0001);
        assert!(matches!(s.verdict, Verdict::Unknown { ended: Ended::Limit, .. }) && n < 100_000, "{s:?} after {n} boxes");
    }

    /// Faces that touch are settled whatever their area, to the tolerance:
    /// a plate on its envelope's floor, whose bounds overlap the floor's by
    /// a few ulps across the face,
    /// settled box by box only by the tolerance's cover of its samples, at
    /// about four tolerances across (it ended unknown after 4.2 million
    /// boxes), is proven to 0.01 in one: a pair of faces facing each other
    /// along an axis leaves from any point of both within the tolerance,
    /// whatever the box's size. Moved 0.1 through the floor, it is a
    /// witness. Where two
    /// such faces meet (the plate as wide as its envelope: its edges on the
    /// envelope's), the envelope's walls are a choice, through which no
    /// exit is shown, and the edges are refined to the tolerance: boxes as
    /// their length, not as the faces' area (127 thousand here).
    #[test]
    fn faces_that_touch_are_covered_whatever_their_area() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        for (side, lift, holds, most) in [(49.0, 2.0, true, 5_000), (49.0, 1.9, false, 5_000), (50.0, 2.0, true, 200_000)] {
            let mut b = B::new();
            let plate = b.cuboid([0.0, 0.0, lift], [side, side, 2.0]);
            let envelope = b.cuboid([0.0, 0.0, 15.0], [50.0, 50.0, 15.0]);
            let outside = b.un(UnOp::Neg, envelope);
            let f = b.bin(BinOp::Max, plate, outside);
            let tape = Arc::new(b.tape(f));
            let (s, n) = empty_now(&tape, &Ctx::new(0.0, vec![]), [-51.0, -51.0, -1.0], [51.0, 51.0, 31.0], None, 0.005, 0.01);
            let ok = match s.verdict {
                Verdict::Proven => holds && !s.exact,
                Verdict::Disproven { witness, .. } => !holds && witness[2] < 0.0,
                _ => false,
            };
            assert!(ok && n < most, "plate {side} lifted {lift}: {s:?} after {n} boxes");
        }
    }

    /// A block 1.5 x 4 x 4 about x = 1e16 + 1, which is no float: its
    /// narrowed box is [1e16, 1e16 + 2], whose float middle is a face. How
    /// deep a point can lie is the real half-width, 1, not the depth at
    /// that middle (0), so no search shows it empty to 0.1 (and none can
    /// show it full: no float point is inside).
    #[test]
    fn depth_is_half_the_width_not_the_float_middles() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        let mut b = B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (big, one, w) = (b.c(1e16), b.c(1.0), b.c(0.75));
        let d = b.bin(BinOp::Sub, x, big);
        let d = b.bin(BinOp::Sub, d, one);
        let d = b.un(UnOp::Abs, d);
        let fx = b.bin(BinOp::Sub, d, w);
        let (fy, fz) = (b.slab(y, 0.0, 2.0), b.slab(z, 0.0, 2.0));
        let f = b.max3(fx, fy, fz);
        let tape = Arc::new(b.tape(f));
        let ctx = Ctx::new(0.0, vec![]);
        let Extent::Bounded(lo, hi) = extent(&tape, &ctx, ([-500.0; 3], [500.0; 3])) else { panic!("unbounded") };
        let (lo, hi, _) = bounds(&tape, &ctx, lo, hi).unwrap();
        for tol in [0.0, 0.1, 0.5] {
            let s = prove_empty(&tape, &ctx, lo, hi, None, 0.05, tol);
            match s.verdict {
                Verdict::Unknown { attained, .. } => assert!(attained.is_none_or(|a| a >= 0.75), "{tol}: {s:?}"),
                _ => panic!("{tol}: {s:?}"),
            }
        }
    }

    /// The radii bound the real ones: squares that underflow (1e-200) or
    /// overflow (1e154) and a norm whose float rounds below the real one.
    #[test]
    fn radial_extent_bounds_the_real_radii() {
        let mut g = Graph::new();
        let n = g.constant(-1.0);
        let tape = Tape::new(&g, n);
        let ctx = Ctx::new(0.0, vec![]);
        let at = |lo: V3, hi: V3| radial_extent(&tape, &ctx, lo, hi, false, 0.0, 100).unwrap();
        let (lower, upper) = at([1e-200, 0.0, 0.0], [2e-200, 0.0, 0.0]);
        assert!(lower > 0.0 && lower <= 2e-200 && upper >= 2e-200, "{lower} {upper}");
        let (lower, upper) = at([1e154, 1e154, 0.0], [1e154, 1e154, 0.0]);
        // sqrt(2) 1e154 = 1.41421356237...e154, not (f64::MAX, inf).
        assert!(lower < 1.4142135624e154 && upper > 1.4142135623e154 && upper - lower < 1e140, "{lower} {upper}");
        // sqrt(x^2 + y^2) in floats rounds to 1.5950777814546493, below the
        // real radius (by an exact rational check), and its next float down.
        let p = [1.0506403964780062, 1.2001782726657002, 0.0];
        let (lower, upper) = at(p, p);
        assert!(lower < 1.5950777814546493 && upper > 1.5950777814546493 && upper - lower < 1e-14, "{lower} {upper}");
    }

    /// Faces whose slopes are so small that weight times slope underflows
    /// to 0 bound no distance: the slab 0 < x < 2 scaled by 1e-323 is not
    /// covered by its faces (their reach once came out negative).
    #[test]
    fn faces_of_a_field_scaled_to_subnormals_cover_nothing() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp, Var};
        let mut b = B::new();
        let (k, x) = (b.c(1e-323), b.v(Var::X));
        let kx = b.bin(BinOp::Mul, k, x);
        let k2 = b.c(2e-323);
        let right = b.bin(BinOp::Sub, kx, k2);
        let left = b.un(UnOp::Neg, kx);
        let f = b.bin(BinOp::Max, right, left);
        let tape = Arc::new(b.tape(f));
        let (s, _) = empty_now(&tape, &Ctx::new(0.0, vec![]), [-10.0; 3], [10.0; 3], None, 0.01, 0.01);
        assert!(!matches!(s.verdict, Verdict::Proven), "{s:?}");
    }

    /// A part outside a zone that is its own feature (a screw's thread, an
    /// overlap a contact leaves out: max(part, other, -feature)) is settled
    /// at any box size: the max holds the feature's value and its negation,
    /// whose mean is 0 at every point however the feature curves: one box
    /// (bounds of the two taken apart lose it, and the sphere's surface
    /// inside the block settled only by the tolerance, at about four
    /// tolerances across: 273 thousand boxes). An overlap outside the
    /// feature (the head) is still found.
    #[test]
    fn a_part_outside_its_own_feature_is_settled_at_any_size() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        for (head, holds) in [(9.0, true), (6.5, false)] {
            let mut b = B::new();
            let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
            let sq = [x, y, z].map(|v| b.un(UnOp::Square, v));
            let rr = b.bin(BinOp::Add, sq[0], sq[1]);
            let rr = b.bin(BinOp::Add, rr, sq[2]);
            let (r, five) = (b.un(UnOp::Sqrt, rr), b.c(5.0));
            let feature = b.bin(BinOp::Sub, r, five);
            let cap = b.cuboid([0.0, 0.0, head], [1.0, 1.0, 1.0]);
            let part = b.bin(BinOp::Min, feature, cap);
            let other = b.cuboid([0.0, 0.0, 4.0], [4.5, 4.5, 2.5]);
            let out = b.un(UnOp::Neg, feature);
            let f = b.bin(BinOp::Max, part, other);
            let f = b.bin(BinOp::Max, f, out);
            let tape = Arc::new(b.tape(f));
            let (s, n) = empty_now(&tape, &Ctx::new(0.0, vec![]), [-6.0; 3], [6.0, 6.0, 11.0], None, 0.005, 0.01);
            let ok = match s.verdict {
                Verdict::Proven => holds,
                Verdict::Disproven { witness, .. } => !holds && witness[2] > 5.5,
                _ => false,
            };
            assert!(ok && n < 5_000, "head at {head}: {s:?} after {n} boxes");
        }
    }

    /// A box waiting keeps the choices its test made, not a tape: the
    /// tapes the boxes left hold are their parents', each shared by its
    /// halves, so fewer than the boxes (681 boxes hold 304; kept each its
    /// own, 623: benchmark 02's mesh-sweep peaked at 386 MB, 135 MB so).
    #[test]
    fn boxes_waiting_share_their_parents_tapes() {
        use crate::test_fields::B;
        // Sixty-four cubes in rows, each touching the next, under a plate
        // pressed 0.015 into their tops (thinner than twice the tolerance,
        // too thick for the faces' own cover): the boxes along the contact
        // keep narrowing the union to the cubes near them.
        let mut b = B::new();
        let cubes: Vec<_> = (0..64).map(|k| b.cuboid([(k % 8) as f64, (k / 8) as f64, 0.0], [0.5, 0.5, 0.5])).collect();
        let union = cubes.into_iter().reduce(|x, y| b.bin(ink_core::graph::BinOp::Min, x, y)).unwrap();
        let plate = b.cuboid([3.5, 3.5, 0.985], [4.0, 4.0, 0.5]);
        let f = b.bin(ink_core::graph::BinOp::Max, union, plate);
        let tape = Arc::new(b.tape(f));
        let (lo, hi) = ([-0.5, -0.5, -0.5], [7.5, 7.5, 1.5]);
        let region = std::array::from_fn(|k| I::new(lo[k], hi[k]));
        let q = Empty { tol: 0.01, region, asked: region, loose: Default::default(), t_band: None };
        let root = search::Node { b: lattice(lo, hi, 0.01).unwrap(), ctx: Arc::new(Ctx::new(0.0, vec![])), state: (Kept::of(&tape), Kept::of(&tape)) };
        let mut s = search::Search::new(&q, vec![root], search::Precision::holds(0.015));
        for _ in 0..20 {
            if s.round(&q, BOXES).is_some() {
                break;
            }
        }
        let left: Vec<_> = s.left().collect();
        let held: std::collections::HashSet<*const Tape> = left.iter().map(|n| n.state.0.held()).collect();
        assert!(left.len() > 500 && 3 * held.len() < 2 * left.len(), "{} boxes hold {} tapes", left.len(), held.len());
    }

    /// A box's tape is made when it is read, not when it is cloned.
    #[test]
    fn a_kept_tape_is_made_when_read() {
        use crate::test_fields::B;
        let mut b = B::new();
        let (u, v) = (b.cuboid([0.0; 3], [1.0; 3]), b.cuboid([3.0, 0.0, 0.0], [1.0; 3]));
        let f = b.bin(ink_core::graph::BinOp::Min, u, v);
        let tape = Arc::new(b.tape(f));
        let mut ch = vec![BOTH; tape.len()];
        eval_interval(&tape, &Ctx::new(0.0, vec![]), [I::new(-1.0, 1.0); 3], Some(&mut ch));
        let k = Kept::after(&tape, vec![ch]);
        let c = k.clone();
        assert!(std::ptr::eq(c.held(), Arc::as_ptr(&tape)), "made as cloned");
        assert!(c.tape().len() < tape.len() && std::ptr::eq(k.held(), Arc::as_ptr(&c.tape())), "made once, shared");
    }

    /// A ball about (cx, 0, 0) of radius r: its distance, or (squared) the
    /// square of it less r^2, a field with no slope bound of 1.
    fn ball(b: &mut crate::test_fields::B, cx: f64, r: f64, squared: bool) -> ink_core::graph::NodeId {
        use ink_core::graph::{BinOp, UnOp};
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let c = b.c(cx);
        let dx = b.bin(BinOp::Sub, x, c);
        let (xx, yy, zz) = (b.bin(BinOp::Mul, dx, dx), b.bin(BinOp::Mul, y, y), b.bin(BinOp::Mul, z, z));
        let s = b.bin(BinOp::Add, xx, yy);
        let s = b.bin(BinOp::Add, s, zz);
        match squared {
            true => {
                let r2 = b.c(r * r);
                b.bin(BinOp::Sub, s, r2)
            }
            false => {
                let l = b.un(UnOp::Sqrt, s);
                let r = b.c(r);
                b.bin(BinOp::Sub, l, r)
            }
        }
    }

    /// A field x + 10 y cos p - 10 y cos(p + 0.001) (plus a millionth of
    /// each of eight other params, all over [-2, 2]) moves with p only
    /// through its slope: a millionth's change in value, but over p's span
    /// its slope bound reaches 40, so a box 0.02 across beside it at 1 to 2
    /// along x is cleared by 0.5 only once p is halved. With p the ninth
    /// param, past the eight a first-order bound takes at once, its blur
    /// was not measured (taken as 0, then as what its span adds to the
    /// plain bound, which y's two reads swamp) and the search ended unknown.
    #[test]
    fn a_param_past_the_first_eight_is_halved_like_them() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        for special in [0, 8] {
            let mut b = B::new();
            let (x, y) = (b.v(Var::X), b.v(Var::Y));
            let ten = b.c(10.0);
            let yy = b.bin(BinOp::Mul, y, ten);
            let p = b.v(Var::Param(special));
            let eps = b.c(0.001);
            let shifted = b.bin(BinOp::Add, p, eps);
            let (cp, cq) = (b.un(UnOp::Cos, p), b.un(UnOp::Cos, shifted));
            let (xp, xq) = (b.bin(BinOp::Mul, yy, cp), b.bin(BinOp::Mul, yy, cq));
            let mut f = b.bin(BinOp::Sub, xp, xq);
            for k in (0..9).filter(|&k| k != special) {
                let (pk, tiny) = (b.v(Var::Param(k)), b.c(1e-6));
                let t = b.bin(BinOp::Mul, pk, tiny);
                f = b.bin(BinOp::Add, f, t);
            }
            let f = b.bin(BinOp::Add, f, x);
            let bar = b.cuboid([1.5, 0.0, 0.0], [0.5, 0.01, 0.01]);
            let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..9).map(|k| (k, (-2.0, 2.0)))), ..Ctx::new(0.0, vec![0.0; 9]) };
            let way = Way { grown: Arc::new(b.tape(f)), covered: Arc::new(b.tape(bar)), lo: [1.0, -0.01, -0.01], hi: [2.0, 0.01, 0.01], swapped: false };
            let (v, ..) = clearance_now(&[way], &ctx, 0.5, 0.025, 20_000);
            assert!(matches!(v, Verdict::Proven), "p{special}: {v:?}");
        }
    }

    /// A box reaching infinity is tested whole, never halved, so its
    /// infinite spread does not stop a param from being halved: B the
    /// half-space x > 1 + sin p - sin(p + 0.001) (plus a millionth of each
    /// of eight other params over [-0.5, 0.5]), A x < 0, asked over all of
    /// space: proven once p is halved, first or ninth. Ranked with them,
    /// space's infinite smear would outrank every param, and the box would
    /// be stuck at once.
    #[test]
    fn a_box_over_all_of_space_halves_its_params() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        for special in [0, 8] {
            let mut b = B::new();
            let x = b.v(Var::X);
            let p = b.v(Var::Param(special));
            let eps = b.c(0.001);
            let shifted = b.bin(BinOp::Add, p, eps);
            let (s0, s1) = (b.un(UnOp::Sin, p), b.un(UnOp::Sin, shifted));
            let mut f = b.bin(BinOp::Sub, s0, s1);
            for k in (0..9).filter(|&k| k != special) {
                let (pk, tiny) = (b.v(Var::Param(k)), b.c(1e-6));
                let t = b.bin(BinOp::Mul, pk, tiny);
                f = b.bin(BinOp::Add, f, t);
            }
            let one = b.c(1.0);
            let f = b.bin(BinOp::Add, one, f);
            let covered = b.bin(BinOp::Sub, f, x);
            let ctx = Ctx { p_spans: ink_core::eval::ParamSpans::new((0..9).map(|k| (k, (-0.5, 0.5)))), ..Ctx::new(0.0, vec![0.0; 9]) };
            let way = Way { grown: Arc::new(b.tape(x)), covered: Arc::new(b.tape(covered)), lo: [f64::NEG_INFINITY; 3], hi: [f64::INFINITY; 3], swapped: false };
            let (v, ..) = clearance_now(&[way], &ctx, 0.5, 0.1, 20_000);
            assert!(matches!(v, Verdict::Proven), "p{special}: {v:?}");
        }
    }

    /// The one search's cost, counted, both ways round. A pin of radius 5
    /// in a bore of 5.5 cleared by 0.49 (a margin of 0.01) is proven by the
    /// facing surfaces' first-order bound at boxes of about sqrt(margin x
    /// radius), not of the margin. A box of slabs reads 30 beside its edge
    /// where the other box is 42.4 away: claimed clear by 31, the box
    /// nearest a violation is refined first and, stuck at the finest
    /// cells, ends the search, which holds to the 30 its field proves. It
    /// tests 45 and 12 thousand boxes: counts far above these mean it
    /// refines where it need not.
    #[test]
    fn the_one_search_decides_in_few_boxes() {
        use crate::test_fields::B;
        use ink_core::graph::BinOp;
        let mut b = B::new();
        let (x, y) = (b.v(Var::X), b.v(Var::Y));
        let (xx, yy) = (b.bin(BinOp::Mul, x, x), b.bin(BinOp::Mul, y, y));
        let r = b.bin(BinOp::Add, xx, yy);
        let r = b.un(ink_core::graph::UnOp::Sqrt, r);
        let (five, bore) = (b.c(5.0), b.c(5.5));
        let pin = b.bin(BinOp::Sub, r, five);
        let hole = b.bin(BinOp::Sub, bore, r);
        let plate = b.cuboid([0.0; 3], [20.0, 20.0, 5.0]);
        let plate = b.bin(BinOp::Max, plate, hole);
        let (plate, pin) = (Arc::new(b.tape(plate)), Arc::new(b.tape(pin)));
        let ways = [(plate.clone(), pin.clone(), false), (pin, plate, true)].map(|(grown, covered, swapped)| Way { grown, covered, lo: [-6.0; 3], hi: [6.0; 3], swapped });
        let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 0.49, 0.0125, BOXES);
        assert!(matches!(v, Verdict::Proven) && n < 100_000, "{v:?} after {n} boxes");
        let mut b = B::new();
        let (a, other) = (b.cuboid([0.0; 3], [50.0; 3]), b.cuboid([130.0, 130.0, 0.0], [50.0; 3]));
        let (a, other) = (Arc::new(b.tape(a)), Arc::new(b.tape(other)));
        let ways = [
            Way { grown: a.clone(), covered: other.clone(), lo: [80.0, 80.0, -50.0], hi: [180.0, 180.0, 50.0], swapped: false },
            Way { grown: other, covered: a, lo: [-50.0; 3], hi: [50.0; 3], swapped: true },
        ];
        let (v, _, _, n) = clearance_now(&ways, &Ctx::new(0.0, vec![]), 31.0, 0.025, BOXES);
        let Verdict::Unknown { ended: Ended::Finest, attained: Some(e), .. } = v else { panic!("{v:?}") };
        assert!((1.0..1.01).contains(&e) && n < 25_000, "holds to {e} after {n} boxes");
    }

    /// An unknown clearance holds to the accuracy it attained: two balls
    /// 0.5 apart, their fields squared (no slope bound of 1, so a box's
    /// bound proves less than its distance), asked just under and just
    /// over the gap with coarse cells. Over it, never proven; where
    /// unknown, the clearance less its attained accuracy is no more than
    /// the gap.
    #[test]
    fn an_unknown_clearance_holds_to_the_accuracy_it_attained() {
        let mut b = crate::test_fields::B::new();
        let (p, q) = (ball(&mut b, 0.0, 2.0, true), ball(&mut b, 4.5, 2.0, true));
        let (p, q) = (Arc::new(b.tape(p)), Arc::new(b.tape(q)));
        let mut attained = 0;
        for d in [0.45, 0.499, 0.501, 0.55] {
            let way = Way { grown: p.clone(), covered: q.clone(), lo: [2.0, -2.0, -2.0], hi: [6.5, 2.0, 2.0], swapped: false };
            match clearance_now(&[way], &Ctx::new(0.0, vec![]), d, 0.05, BOXES).0 {
                Verdict::Proven => assert!(d <= 0.5, "{d} proven past the gap"),
                Verdict::Unknown { attained: Some(e), .. } => {
                    assert!(d - e <= 0.5, "{d} holds to {e}: past the gap");
                    attained += 1;
                }
                v => assert!(d > 0.5 || matches!(v, Verdict::Unknown { .. }), "{d}: {v:?}"),
            }
        }
        assert!(attained > 0, "no unknown attained an accuracy");
    }

    #[test]
    fn a_clearance_that_does_not_read_t_has_its_witness_at_the_moment_asked() {
        // Two balls 1 apart, neither moving: one answer at every moment of
        // the motion (its key holds no t), its witness at the moment asked
        // (a cached one carried the moment it was made at).
        let mut g = Graph::new();
        let ball = |g: &mut Graph, cx: f64| {
            let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
            let c = g.constant(cx);
            let dx = g.sub(x, c);
            let (xx, yy, zz) = (g.mul(dx, dx), g.mul(y, y), g.mul(z, z));
            let s = g.add(xx, yy);
            let s = g.add(s, zz);
            let r = g.un(ink_core::graph::UnOp::Sqrt, s);
            let one = g.constant(1.0);
            g.sub(r, one)
        };
        let (a, b) = (ball(&mut g, -1.5), ball(&mut g, 1.5));
        let (a, b) = (Arc::new(Tape::new(&g, a)), Arc::new(Tape::new(&g, b)));
        for t in [0.25, 0.75] {
            let ctx = Ctx { t_span: Some((0.0, 1.0)), ..Ctx::new(t, vec![]) };
            let way = Way { grown: a.clone(), covered: b.clone(), lo: [0.0, -1.5, -1.5], hi: [3.0, 1.5, 1.5], swapped: false };
            match clearance(&[way], &ctx, 2.0, 0.25, BOXES).0 {
                Verdict::Disproven { t: at, .. } => assert_eq!(at, t),
                v => panic!("{v:?}"),
            }
        }
    }

    #[test]
    fn witnesses_carry_finite_in_span_moments() {
        // 0.5 (t0 + t1) would overflow for t0 = 1e308, t1 = 1.1e308.
        let m = span_mid((1e308, 1.1e308));
        assert!((1e308..=1.1e308).contains(&m), "{m}");
        // 0.5 - t is below 0 over [1, 2]: a whole-cell witness, which must
        // carry a moment of the span, not ctx.t = 0 (where it is 0.5).
        let mut g = Graph::new();
        let (h, t) = (g.constant(0.5), g.var(Var::T));
        let f = g.sub(h, t);
        let tape = Arc::new(Tape::new(&g, f));
        let ctx = Ctx { t_span: Some((1.0, 2.0)), ..Ctx::new(0.0, vec![]) };
        match prove_empty(&tape, &ctx, [-1.0; 3], [1.0; 3], Some(0.5), 0.5, 0.0).verdict {
            Verdict::Disproven { t, witness, .. } => assert!((1.0..=2.0).contains(&t) && witness.iter().all(|x| x.is_finite()), "{t} {witness:?}"),
            v => panic!("{v:?}"),
        }
    }

    #[test]
    fn an_inside_pair_with_a_target_without_interior_ends() {
        // A field 0 over its box (an agent's `(move 0 0 5 (box 4 4 1))`,
        // arguments swapped): a mesh pass finds no surface, and each finer
        // one would test every cell, to 32 768 along the box (past 30 s on
        // 28 cores). The search proves it holds no point: no pair.
        let mut b = crate::test_fields::B::new();
        let slab = b.cuboid([0.0; 3], [5.0, 5.0, 1.0]);
        let zero = b.c(0.0);
        let (ta, tb) = (Arc::new(b.tape(zero)), Arc::new(b.tape(slab)));
        let started = std::time::Instant::now();
        let c = inside_pair(&ta, &tb, &Ctx::new(0.0, vec![]), ([-2.0, -2.0, 4.0], [2.0, 2.0, 6.0]), ([-5.0, -5.0, -1.0], [5.0, 5.0, 1.0]), 1e-3);
        assert!(c.is_none(), "{c:?}");
        assert!(started.elapsed().as_secs_f64() < 10.0, "{:?}", started.elapsed());
    }

    #[test]
    fn a_proof_refines_past_a_pass_with_many_undecided_cells() {
        // A rod (radius 20, 100 long) in a tube 0.002 clear of it: cells of
        // 0.256 leave the gap undecided all round (over 200 000 of them),
        // cells of 0.128 separate the faces. Refining stopped after a pass
        // of 100 000 undecided cells, unknown; the search tests up to BOXES
        // boxes, which this is well within. The tube's radius is written
        // apart (its squares): of one radius, the gap is exact at once (Same).
        use ink_core::graph::{BinOp, UnOp};
        let mut b = crate::test_fields::B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (xx, yy) = (b.bin(BinOp::Mul, x, x), b.bin(BinOp::Mul, y, y));
        let s = b.bin(BinOp::Add, xx, yy);
        let r = b.un(UnOp::Sqrt, s);
        let (sx, sy) = (b.un(UnOp::Square, x), b.un(UnOp::Square, y));
        let s2 = b.bin(BinOp::Add, sx, sy);
        let r2 = b.un(UnOp::Sqrt, s2);
        let (inner, outer) = (b.c(20.0), b.c(20.002));
        let rod = b.bin(BinOp::Sub, r, inner);
        let tube = b.bin(BinOp::Sub, outer, r2);
        let (span, around) = (b.slab(z, 0.0, 50.0), b.cuboid([0.0; 3], [30.0, 30.0, 50.0]));
        let rod = b.bin(BinOp::Max, rod, span);
        let tube = b.bin(BinOp::Max, tube, around);
        let both = b.bin(BinOp::Max, rod, tube);
        let tape = Arc::new(b.tape(both));
        let s = prove_empty(&tape, &Ctx::new(0.0, vec![]), [-31.0, -31.0, -51.0], [31.0, 31.0, 51.0], None, 0.001, 0.0);
        assert!(matches!(s.verdict, Verdict::Proven) && s.res < 0.2 && s.exact, "{s:?}");
    }

    #[test]
    fn a_witness_lies_in_the_region_asked() {
        // max(x - 0.007, -19 x - 0.0075), in a unit cube's slabs: material
        // 0.007 deep at the region's face x = 0. The first term dominates
        // over the region (x in [0, 1]), so a tape focused on it is x -
        // 0.007 alone, which below x = 0 reads ever deeper material where
        // the field is above 0 (from x = -0.0004). A root cell reaching
        // below the region, which that tape reads full, would give the
        // deepest witness (-0.044, 0.956, 0.056), outside the region and the
        // shape. Cells are clipped to the region.
        use ink_core::graph::BinOp;
        let mut b = crate::test_fields::B::new();
        let x = b.v(Var::X);
        let (k3, k19, k31) = (b.c(0.007), b.c(-19.0), b.c(0.0075));
        let a = b.bin(BinOp::Sub, x, k3);
        let steep = b.bin(BinOp::Mul, k19, x);
        let steep = b.bin(BinOp::Sub, steep, k31);
        let f = b.bin(BinOp::Max, a, steep);
        let (y, z) = (b.v(Var::Y), b.v(Var::Z));
        let (sy, sz) = (b.slab(y, 0.5, 0.5), b.slab(z, 0.5, 0.5));
        let f = b.max3(f, sy, sz);
        let tape = Arc::new(b.tape(f));
        let ctx = Ctx::new(0.0, vec![]);
        let opts = crate::ranged::ProofOptions { res: Some(0.05), finest: 0.05, depth: 0, tol: 0.0 };
        let (lo, hi) = ([0.0; 3], [1.0; 3]);
        match crate::ranged::prove(&tape, &ctx, lo, hi, &opts).verdict {
            crate::query::Verdict::Disproven { witness, .. } => {
                assert!((0..3).all(|k| (lo[k]..=hi[k]).contains(&witness[k])), "{witness:?}");
                assert!(eval_point(&tape, &ctx, witness) < 0.0, "{witness:?}");
            }
            o => panic!("{o:?}"),
        }
    }

    /// A ray finds the material behind a select whose condition straddles 0
    /// by rounding (cos of pi/2 y at y = 1): its bounds there are 10 wide,
    /// the branches' gap, not rounding (a rule taking them for rounding
    /// would drop it), and the material (x in 0.5
    /// .. 1.25 at y = 1, f -0.25 at x = 1) is found where it starts,
    /// unsure: only rounding decides the select there, so bounds show no
    /// point of it inside.
    #[test]
    fn a_ray_finds_material_behind_an_undecided_select() {
        use crate::test_fields::B;
        use ink_core::graph::{BinOp, UnOp};
        let mut b = B::new();
        let (x, y, z) = (b.v(Var::X), b.v(Var::Y), b.v(Var::Z));
        let (half, k) = (b.c(0.5), b.c(std::f64::consts::FRAC_PI_2));
        let ky = b.bin(BinOp::Mul, k, y);
        let cond = b.un(UnOp::Cos, ky);
        let (ten, zero) = (b.c(10.0), b.c(0.0));
        let sel = b.0.select(cond, ten, zero);
        let ramp = b.bin(BinOp::Sub, half, x);
        let a = b.bin(BinOp::Add, ramp, sel);
        let (sx, sy, sz) = (b.slab(x, 0.75, 0.5), b.slab(y, 1.0, 1.0), b.slab(z, 0.0, 1.0));
        let m = b.max3(sx, sy, sz);
        let f = b.bin(BinOp::Max, a, m);
        let tape = Arc::new(b.tape(f));
        let ctx = Ctx::new(0.0, vec![]);
        let t = ray_to(&Local::whole(&tape, &ctx), [0.0, 1.0, 0.0], [1.0, 0.0, 0.0], (0.0, 2.0), Side::Inside, 2e-9);
        assert!(t.is_some_and(|h| (h.t - 0.5).abs() < 1e-6 && h.rounding.is_some()), "{t:?}");
    }

    #[test]
    fn ray_segment_boxes_hold_the_real_segment() {
        // o + t dir correctly rounded (fma) lies in the box; nearest ends
        // (a rounded product, then a rounded sum) often miss it.
        let mut r = 0x9e3779b97f4a7c15u64;
        let mut f = || {
            r ^= r << 13;
            r ^= r >> 7;
            r ^= r << 17;
            (r >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        for _ in 0..10000 {
            let (o, dir, t) = ([f(), f(), f()], [f(), f(), f()], f().abs() * 10.0);
            let b = segment_box(o, dir, t, t);
            for k in 0..3 {
                assert!(b[k].contains(t.mul_add(dir[k], o[k])), "{o:?} {dir:?} {t}");
            }
        }
    }
}
