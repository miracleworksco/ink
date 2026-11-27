//! Planar material boundaries and marked sections.

use ink_core::V3;
use ink_core::eval::{Ctx, eval_interval, eval_point, eval_points};
use ink_core::interval::I;
use ink_core::sign::{drawn, inside};
use ink_core::tape::{BOTH, Tape};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

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
    slice_checked(tape, ctx, axis, at, lo, hi, h).0
}

/// `slice`, and the samples on 0 its loops read that no bound near them
/// decides (`ink_core::sign::drawn`): drawn as air, unresolved.
pub fn slice_checked(tape: &Arc<Tape>, ctx: &Ctx, axis: usize, at: f64, lo: V3, hi: V3, h: f64) -> (Vec<Loop>, Vec<V3>) {
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
    let open = std::sync::Mutex::new(Vec::new());
    let segs = quad_descend(tape, tape, ctx, &to3, o, h, axis, at, (u, v), [0, 0], size, &open);
    let mut open = open.into_inner().unwrap();
    open.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    open.dedup();
    (join(segs), open)
}

/// A mark's slice: where field `g` crosses 0 in the plane `axis = at`,
/// inside `region` (r < 0) and the box lo..hi. Marching squares on `g` (as
/// `slice`), each chain then cut where it leaves them, the cut found by
/// bisection along its segment. A closed chain wholly inside stays closed;
/// the rest are open, with area 0. A mark is where g changes sign, as an
/// outline is where f does: a g touching 0 draws nothing there. With the
/// samples on 0 in the region that no bound decides (`slice_checked`).
#[allow(clippy::too_many_arguments)]
pub fn slice_mark(g: &Arc<Tape>, region: &Tape, ctx: &Ctx, axis: usize, at: f64, lo: V3, hi: V3, h: f64) -> (Vec<Loop>, Vec<V3>) {
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
            if inside(r(m)) == inside(ra) { a = m } else { b = m }
        }
        [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])]
    };
    let length = |pts: &[[f64; 2]]| pts.windows(2).map(|w| ((w[1][0] - w[0][0]).powi(2) + (w[1][1] - w[0][1]).powi(2)).sqrt()).sum::<f64>();
    let mut out = Vec::new();
    let (loops, open) = slice_checked(g, ctx, axis, at, lo, hi, h);
    // Those whose cells (a cell each way about it, in the plane) the region
    // may reach within the box, by its bounds there: where it decides what
    // is drawn. Never by the region's value at it, which is no distance.
    let open: Vec<V3> = open
        .into_iter()
        .filter(|q| {
            let c: [I; 3] = std::array::from_fn(|k| if k == axis { I::point(at) } else { I::new((q[k] - h).max(lo[k]), (q[k] + h).min(hi[k])) });
            c.iter().all(|x| x.lo <= x.hi) && eval_interval(region, ctx, c, None).may_be_inside()
        })
        .collect();
    for l in loops {
        let mut pts = l.pts;
        let mut rs: Vec<f64> = pts.iter().map(|&p| r(p)).collect();
        if l.closed {
            // Start a closed chain outside the region, if it leaves it.
            match rs.iter().position(|&x| !inside(x)) {
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
            let within = inside(rs[k]);
            if k > 0 && within != inside(rs[k - 1]) {
                cur.push(cross(pts[k - 1], pts[k], rs[k - 1]));
                if !within {
                    out.push(Loop { perimeter: length(&cur), pts: std::mem::take(&mut cur), closed: false, area: 0.0 });
                }
            }
            if within {
                cur.push(pts[k]);
            }
        }
        if cur.len() >= 2 {
            out.push(Loop { perimeter: length(&cur), pts: cur, closed: false, area: 0.0 });
        }
    }
    out.retain(|l| l.pts.len() >= 2);
    (out, open)
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
    open: &std::sync::Mutex<Vec<V3>>,
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
    if r.none_inside() || r.all_inside() {
        return Vec::new();
    }
    let sub = tape.narrowed(&ch);
    if size == 1 {
        return march(full, &sub, ctx, to3, o, h, i, open);
    }
    let hs = size / 2;
    let kids = [[0, 0], [hs, 0], [0, hs], [hs, hs]];
    let parts: Vec<Vec<Seg>> = kids.par_iter().map(|k| quad_descend(full, &sub, ctx, to3, o, h, axis, at, uv, [i[0] + k[0], i[1] + k[1]], hs, open)).collect();
    parts.into_iter().flatten().collect()
}

#[allow(clippy::too_many_arguments)]
fn march(full: &Tape, tape: &Tape, ctx: &Ctx, to3: &(dyn Fn([f64; 2]) -> V3 + Sync), o: [f64; 2], h: f64, i: [i64; 2], open: &std::sync::Mutex<Vec<V3>>) -> Vec<Seg> {
    let p = |a: i64, b: i64| [o[0] + h * a as f64, o[1] + h * b as f64];
    let (x, y) = (i[0], i[1]);
    let c = [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)].map(|(a, b)| p(a, b));
    let center = [o[0] + h * (x as f64 + 0.5), o[1] + h * (y as f64 + 0.5)];
    let mut v = Vec::new();
    let pts = [to3(c[0]), to3(c[1]), to3(c[2]), to3(c[3]), to3(center)];
    eval_points(tape, ctx, &pts, &mut v);
    // The curve is where the field changes sign (ink_core::sign): a corner
    // on 0 inside the material reads with it, so a 0 with one sign round it
    // draws nothing. Every cell sharing the corner reads the same
    // neighbours (unsimplified: they lie outside the cell).
    let undecided = drawn(full, ctx, &pts, &mut v, None);
    if !undecided.is_empty() {
        open.lock().unwrap().extend(undecided);
    }
    let ins = [0, 1, 2, 3].map(|k| inside(v[k]));
    // Edges (from lower to higher grid corner): e0 c0-c1, e1 c1-c2, e2 c3-c2, e3 c0-c3.
    let edges: [(usize, usize, EdgeKey); 4] = [(0, 1, (x, y, 0)), (1, 2, (x + 1, y, 1)), (3, 2, (x, y + 1, 0)), (0, 3, (x, y, 1))];
    let cross = |e: usize| {
        let (a, b, _) = edges[e];
        // Both corners on 0 (one in air the field reads as 0): the middle.
        let d = v[a] - v[b];
        let t = if d == 0.0 { 0.5 } else { (v[a] / d).clamp(0.0, 1.0) };
        [c[a][0] + (c[b][0] - c[a][0]) * t, c[a][1] + (c[b][1] - c[a][1]) * t]
    };
    let crossing: Vec<usize> = (0..4).filter(|&e| ins[edges[e].0] != ins[edges[e].1]).collect();
    // Pair crossing edges into segments.
    let pairs: Vec<(usize, usize)> = match crossing.len() {
        2 => vec![(crossing[0], crossing[1])],
        4 => {
            if inside(v[4]) == ins[0] {
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
