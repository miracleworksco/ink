//! Material between a cell's samples.
//!
//! A cell knows the part by its 4^3 samples, so a foil, web or ligament
//! thinner than their spacing can lie between them unseen. Each line
//! between two neighbouring samples outside the part (and from the outer
//! samples to the cell's faces) is walked: from a point outside reading f,
//! no material lies within f / L of it along the line, L a bound on the
//! field's slope over the line (its interval gradient: a field need not be a
//! distance, the cube of one reads 8 at 2 mm). Where that reach falls under
//! the floor (an eighth of the cell, the smallest it splits to), the stretch
//! ahead is settled by the field's bounds over it, and where they cannot
//! settle it, what lies there is taken as found. A line within a sample's
//! spacing of a sample inside the part is that sample's own surface (a
//! corner's tip, a line grazing a face) and is not walked.
//! Material found in a cell whose samples see none is sorted by what it
//! joins (`changes_connection`): a tip of seen material changes nothing; a
//! foil joining nothing seen, two pieces apart, or one piece at two places
//! (a ligament beside it), makes every claim unknown, naming the place.

use crate::eval::{Ctx, eval_interval, eval_points};
use crate::tape::Tape;
use ink_core::{V3, eval, interval};
use std::collections::{HashMap, HashSet};

/// Whether the segment p..q holds no material, by the field's bounds over
/// its box: where they cannot say, what lies there is taken as found and
/// sorted by what it joins (`changes_connection`), a tip of the part's own
/// face (a fillet grazed) joining what the cells hold in one place.
fn clear(tape: &Tape, ctx: &Ctx, p: V3, q: V3) -> bool {
    let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(p[a].min(q[a]), p[a].max(q[a])));
    eval_interval(tape, ctx, b, None).lo > 0.0
}

/// The samples' lattice over the whole grid: sample `g` (each cell's 4 a
/// side, `g = 4 cell + j`) at `o + (g + 1/2) s`, and whether it is inside
/// the part (`inside`: read off the cells' own samples, no evaluation).
pub(crate) struct Lattice<'a> {
    pub o: V3,
    pub s: f64,
    /// A cell's samples inside the part, as its 64 bits (bit a + 4 b + 16 c),
    /// all of them for a full cell, none for one the samples see nothing in.
    pub cell: &'a (dyn Fn([i64; 3]) -> u64 + Sync),
}

/// The samples within a cell of cell `i` (its own and its 26 neighbours'),
/// inside the part or not, by lattice offset from the cell's first sample
/// less 4 (-4 .. 8 each way): read once per cell walked.
struct Near {
    base: [i64; 3],
    bits: Vec<bool>,
}

impl Near {
    fn of(lat: &Lattice, i: [i64; 3]) -> Near {
        let mut bits = vec![false; 12 * 12 * 12];
        for n in 0..27 {
            let d = [n % 3 - 1, n / 3 % 3 - 1, n / 9 - 1];
            let m = (lat.cell)([i[0] + d[0], i[1] + d[1], i[2] + d[2]]);
            if m == 0 {
                continue;
            }
            for k in 0..64 {
                if m >> k & 1 == 1 {
                    let g = [(d[0] + 1) * 4 + (k % 4) as i64, (d[1] + 1) * 4 + (k / 4 % 4) as i64, (d[2] + 1) * 4 + (k / 16) as i64];
                    bits[(g[0] + 12 * (g[1] + 12 * g[2])) as usize] = true;
                }
            }
        }
        Near { base: i.map(|v| 4 * v - 4), bits }
    }
    /// Whether lattice sample `g` (global) is inside the part (outside the
    /// 27 cells: taken as not).
    fn at(&self, g: [i64; 3]) -> bool {
        let l: [i64; 3] = std::array::from_fn(|a| g[a] - self.base[a]);
        l.iter().all(|&v| (0..12).contains(&v)) && self.bits[(l[0] + 12 * (l[1] + 12 * l[2])) as usize]
    }
    /// Whether a sample inside the part lies within a spacing of lattice
    /// sample `g` (its 27 nearest).
    fn around(&self, g: [i64; 3]) -> bool {
        (0..27).any(|n| self.at([g[0] + n % 3 - 1, g[1] + n / 3 % 3 - 1, g[2] + n / 9 - 1]))
    }
    /// The same of a point `x`: there, material is that sample's own surface
    /// (a corner's tip, the part's face beside a cell's face), not a feature
    /// the samples lost.
    fn near_inside(&self, lat: &Lattice, x: V3) -> bool {
        self.around(std::array::from_fn(|a| ((x[a] - lat.o[a]) / lat.s - 0.5).round() as i64))
    }
}

/// The walk along p..q: where it finds material (or cannot settle a stretch).
fn walk(tape: &Tape, ctx: &Ctx, p: V3, q: V3, floor: f64, lat: &Lattice, near: &Near) -> Option<V3> {
    let b: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(p[a].min(q[a]), p[a].max(q[a])));
    // (no material anywhere on it: one interval)
    if eval_interval(tape, ctx, b, None).lo > 0.0 {
        return None;
    }
    let l = interval::lipschitz(&eval::eval_interval_grad(tape, ctx, b));
    let d: V3 = std::array::from_fn(|a| q[a] - p[a]);
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let at = |t: f64| -> V3 { std::array::from_fn(|a| p[a] + d[a] * (t / len).min(1.0)) };
    let mut t = 0.0;
    let mut v = Vec::new();
    loop {
        let x = at(t);
        v.clear();
        eval_points(tape, ctx, &[x], &mut v);
        if v[0] <= 0.0 && !near.near_inside(lat, x) {
            return Some(x);
        }
        if t >= len || l == 0.0 {
            return None;
        }
        let reach = if l.is_finite() && v[0] > 0.0 { v[0] / l } else { 0.0 };
        if reach >= floor {
            t += reach;
            continue;
        }
        let end = (t + floor).min(len);
        if !near.near_inside(lat, x) && !clear(tape, ctx, x, at(end)) {
            return Some(x);
        }
        t = end;
    }
}

/// A cell (index `i`, centre `ctr`, side `h`; `tape` the field over it)
/// whose samples all lie outside the part: where a line between its
/// samples, or from them to its faces, meets material the samples do not
/// see. A line within a spacing of a sample inside the part (its own or a
/// neighbour cell's) is that sample's surface and is not walked: a part of
/// one material with no feature thinner than a quarter of a cell walks only
/// the lines a spacing or more from its face, each settled by one interval.
pub(crate) fn missed_in_cell(tape: &Tape, ctx: &Ctx, i: [i64; 3], ctr: V3, h: f64, lat: &Lattice) -> Option<V3> {
    // the cell's 4 samples a side at -0.75 .. 0.75, and its faces at +-1
    let xs = [-1.0, -0.75, -0.25, 0.25, 0.75, 1.0];
    let pt = |j: [usize; 3]| -> V3 { std::array::from_fn(|a| ctr[a] + 0.5 * h * xs[j[a]]) };
    let floor = h / 8.0;
    // its samples' lattice indices (`j` 1..4 a sample; 0 and 5 its faces,
    // read as the neighbour's sample beyond)
    let g = |j: [usize; 3]| -> [i64; 3] { std::array::from_fn(|a| 4 * i[a] + j[a] as i64 - 1) };
    let nb = Near::of(lat, i);
    let near = |j: [usize; 3]| -> bool { nb.around(g(j)) };
    // (no material in the cell's box at all, by one interval: then no line
    // of it meets any)
    let cb: [interval::I; 3] = std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * h, ctr[a] + 0.5 * h));
    if eval_interval(tape, ctx, cb, None).lo > 0.0 {
        return None;
    }
    for ax in 0..3 {
        for a in 1..5 {
            for b in 1..5 {
                let idx = |n: usize| -> [usize; 3] {
                    let mut r = [0; 3];
                    r[ax] = n;
                    r[(ax + 1) % 3] = a;
                    r[(ax + 2) % 3] = b;
                    r
                };
                // the line's stretches not beside a sample inside the part,
                // each a run of its segments: one interval over a run first
                let mut k = 0;
                while k < 5 {
                    if near(idx(k)) || near(idx(k + 1)) {
                        k += 1;
                        continue;
                    }
                    let mut e = k + 1;
                    while e < 5 && !(near(idx(e)) || near(idx(e + 1))) {
                        e += 1;
                    }
                    let (p, q) = (pt(idx(k)), pt(idx(e)));
                    let rb: [interval::I; 3] = std::array::from_fn(|c| interval::I::new(p[c].min(q[c]), p[c].max(q[c])));
                    if eval_interval(tape, ctx, rb, None).lo <= 0.0 {
                        for m in k..e {
                            if let Some(x) = walk(tape, ctx, pt(idx(m)), pt(idx(m + 1)), floor, lat, &nb) {
                                return Some(x);
                            }
                        }
                    }
                    k = e;
                }
            }
        }
    }
    None
}

/// Of the cells material was found in between the samples (each with where),
/// those joined (26 neighbours) into clusters, each sorted by the cells
/// round it the samples see material in (`seen`): the first cluster that
/// touches none (material the cells lost), or seen material in two places
/// apart (a web joining two pieces, or a ligament beside one), its cell,
/// where, and what it joins. A cluster touching seen material in one place
/// is its tip, and changes nothing.
pub(crate) fn changes_connection(found: Vec<([i64; 3], V3)>, seen: &dyn Fn([i64; 3]) -> bool) -> Option<([i64; 3], V3, String)> {
    let at: HashMap<[i64; 3], V3> = found.into_iter().collect();
    let nb = |i: [i64; 3]| -> Vec<[i64; 3]> {
        let mut out = Vec::with_capacity(26);
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if (dx, dy, dz) != (0, 0, 0) {
                        out.push([i[0] + dx, i[1] + dy, i[2] + dz]);
                    }
                }
            }
        }
        out
    };
    let mut done: HashSet<[i64; 3]> = HashSet::new();
    let mut keys: Vec<[i64; 3]> = at.keys().copied().collect();
    keys.sort();
    for f in keys {
        if !done.insert(f) {
            continue;
        }
        let (mut cluster, mut stack) = (vec![f], vec![f]);
        while let Some(c) = stack.pop() {
            for n in nb(c) {
                if at.contains_key(&n) && done.insert(n) {
                    cluster.push(n);
                    stack.push(n);
                }
            }
        }
        let touch: HashSet<[i64; 3]> = cluster.iter().flat_map(|&c| nb(c)).filter(|&n| !at.contains_key(&n) && seen(n)).collect();
        let mut grouped: HashSet<[i64; 3]> = HashSet::new();
        let mut groups = 0;
        for &t in &touch {
            if !grouped.insert(t) {
                continue;
            }
            groups += 1;
            let mut st = vec![t];
            while let Some(c) = st.pop() {
                for n in nb(c) {
                    if touch.contains(&n) && grouped.insert(n) {
                        st.push(n);
                    }
                }
            }
        }
        if groups != 1 {
            let what = if groups == 0 { "joins no material the cells hold".to_string() } else { format!("joins the material the cells hold in {groups} places apart") };
            return Some((f, at[&f], what));
        }
    }
    None
}
