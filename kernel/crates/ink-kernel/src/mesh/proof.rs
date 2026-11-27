use super::brick::{B, Brick, N, cell_of, vidx};
use super::report::Ran;
use crate::octree::Grid;
use ink_core::V3;
use ink_core::eval::{Ctx, eval_grads, eval_interval};
use ink_core::interval::I;
use ink_core::sign::inside;
use ink_core::tape::{Pass, Sheet, Tape};
use ink_core::v3::norm;
use rayon::prelude::*;

/// An open cell's corners are inside; the search showed it holds nothing
/// the mesh does not (each part of it keeps one sign, or lies in the box of
/// a crossed cell shown drawn); a sample in it has the other sign; its
/// corners do not agree (the surface crosses it, drawn, but may hold more).
/// Else it is undecided.
pub(super) const INSIDE: u8 = 1;
pub(super) const PROVEN: u8 = 2;
pub(super) const FOUND: u8 = 4;
pub(super) const CROSSED: u8 = 8;
/// A crossed cell whose box holds a seam (Tape::switches).
pub(super) const SEAM: u8 = 16;

/// Whether the field's range over `lo..hi` (within the grid's region)
/// proves it drawn on side `inside_` throughout: solid, or air
/// (ink_core::sign: split where pieces meet flush; a region outside the
/// grid's is: nothing there is meshed).
pub(super) fn proven(tape: &Tape, ctx: &Ctx, grid: &Grid, lo: V3, hi: V3, inside_: bool) -> bool {
    let b: [I; 3] = std::array::from_fn(|a| I::new(lo[a].max(grid.clip.0[a]), hi[a].min(grid.clip.1[a])));
    if b.iter().any(|i| !(i.lo <= i.hi)) {
        return true;
    }
    if inside_ { ink_core::sign::solid(tape, ctx, b) } else { ink_core::sign::air(tape, ctx, b) }
}

/// The cells of the block of n^3 at `c` of a brick (`pts` its corners,
/// `vals` the field there, `meshed` the cells meshed, the others' corners
/// not sampled, `crossed` those meshed whose corners do not agree), of
/// those meshed, not shown to hold nothing the mesh does not draw, added to `open` with
/// how thick what a crossed one may hide is (0 for the others):
///
/// - A cell whose corners agree, where the field's range does not prove it
///   keeps their sign throughout, nor split where pieces meet flush
///   (`proven`'s solid or air).
/// - A cell the surface crosses (CROSSED, and SEAM where a select may
///   switch there), where the surface is not one sheet over the box its
///   vertex may be placed in, the cell grown by half a cell (`sheet`). Not
///   where only each part of it is: that would not rule out a closed piece
///   across parts, each part one sheet along its own way (argued; no such
///   piece was found to pass).
///
/// Each is shown over a block at once (its box; for crossed cells their
/// boxes together: one sheet over it is one over each), else each eighth
/// is tried, down to a cell (`drawn`: the block's crossed cells are shown
/// already). Each block's tape is its parent's focused on its box, so a
/// proof reads only what can reach it (a sweep's nearest segments, not all
/// of them); a cell reads its block's (2 cells across). `read` counts the
/// instructions run.
/// What a brick's proofs found and ran (`unshown`): its open cells (each
/// cell's index in the brick, state and how thick what a crossed one may
/// hide is), the instructions run (`read`) and tapes made; whether a cell
/// of it is reported with what may hide there as thick as a crossed
/// cell's box (2 h, the most any cell reports): then its region's
/// `thinner_than` is settled (`settled`), and a crossed cell is left
/// untried (`untried`, and whether its blocks were tried), as no proof of
/// it can change the report but its counts; so is every crossed cell of a
/// block, untried as a block. `only`: trying again
/// just these cells (beside a cell the search reads, whose box their
/// proofs may cover), each block of theirs tried that was not.
#[derive(Default)]
pub(super) struct Proofs<'a> {
    pub(super) open: Vec<(u16, u8, f64)>,
    pub(super) read: usize,
    pub(super) tapes: usize,
    pub(super) settled: bool,
    pub(super) untried: Vec<(u16, bool)>,
    pub(super) only: Option<(&'a [bool], &'a [bool])>,
    /// The cells meshed, the others' corners not sampled (all: none).
    pub(super) meshed: Option<&'a [bool]>,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn unshown(tape: &Tape, ctx: &Ctx, grid: &Grid, pts: &[V3], vals: &[f64], crossed: &[bool], c: [usize; 3], n: usize, drawn: bool, untry: bool, p: &mut Proofs) {
    let (b, h) = (B as usize, grid.h);
    let cell = (c[0] + b * (c[1] + b * c[2])) as u16;
    let (lo, hi) = (pts[vidx(c[0], c[1], c[2])], pts[vidx(c[0] + n, c[1] + n, c[2] + n)]);
    let block = |q: usize| (c[0] + q % n) + b * ((c[1] + q / n % n) + b * (c[2] + q / (n * n)));
    if p.only.is_some_and(|(only, _)| !(0..n * n * n).any(|q| only[block(q)])) {
        return;
    }
    // Only the cells meshed are shown (Proofs::meshed): none here, nothing
    // to show.
    let Some(first) = (0..n * n * n).map(block).find(|&k| p.meshed.is_none_or(|m| m[k])) else { return };
    let any_crossed = (0..n * n * n).any(|q| crossed[block(q)]);
    let (focused, drawn) = match (any_crossed, drawn) {
        (false, _) => {
            // The side of a corner sampled.
            let side = inside(vals[vidx(first % b, first / b % b, first / (b * b))]);
            // Within the grid's region: past it nothing is meshed.
            let bx: [I; 3] = std::array::from_fn(|a| I::new(lo[a].max(grid.clip.0[a]), hi[a].min(grid.clip.1[a])));
            if bx.iter().any(|i| !(i.lo <= i.hi)) {
                return;
            }
            // A cell's proof reads its block's tape: a tape of its own costs
            // more to make than reading it saves.
            let (r, focused) = match n {
                1 => {
                    p.read += tape.work(Pass::Box);
                    (eval_interval(tape, ctx, bx, None), None)
                }
                _ => focus(tape, ctx, bx, &mut p.read),
            };
            p.tapes += focused.is_some() as usize;
            if if side { r.all_inside() } else { r.none_inside() } {
                return;
            }
            // Else split where pieces meet flush (ink_core::sign): solid
            // (a join inside the material) or air (a flush mouth).
            let t = focused.as_ref().unwrap_or(tape);
            p.read += 4 * t.work(Pass::Box);
            if if side { ink_core::sign::solid(t, ctx, bx) } else { ink_core::sign::air(t, ctx, bx) } {
                return;
            }
            if n == 1 {
                p.open.push((cell, if side { INSIDE } else { 0 }, 0.0));
                return;
            }
            (focused, true)
        }
        (true, true) => (None, true),
        (true, false) => {
            let (lo, hi) = (lo.map(|x| x - 0.5 * h), hi.map(|x| x + 0.5 * h));
            let bx: [I; 3] = std::array::from_fn(|a| I::new(lo[a], hi[a]));
            let focused = (n > 1).then(|| focus(tape, ctx, bx, &mut p.read).1).flatten();
            p.tapes += focused.is_some() as usize;
            let t = focused.as_ref().unwrap_or(tape);
            // Trying cells again: a block tried before was not shown drawn.
            let blocks_untried = |q: usize| p.only.is_some_and(|(only, tried)| only[block(q)] && !tried[block(q)]);
            if untry && n == 1 {
                p.untried.push((cell, false));
                return;
            }
            if untry || n > 1 && p.only.is_some() && !(0..n * n * n).any(blocks_untried) {
                (focused, false)
            } else {
                // Settled, a cell or block is left untried.
                if p.settled && p.only.is_none() {
                    if n == 1 {
                        p.untried.push((cell, true));
                        return;
                    }
                    let (t, m) = (focused.as_ref().unwrap_or(tape), n / 2);
                    for k in 0..8 {
                        unshown(t, ctx, grid, pts, vals, crossed, [c[0] + (k & 1) * m, c[1] + (k >> 1 & 1) * m, c[2] + (k >> 2 & 1) * m], m, false, true, p);
                    }
                    return;
                }
                let (held, seam) = sheet(t, ctx, lo, hi, &mut p.read);
                if n == 1 && held != Sheet::One {
                    let thin = match held {
                        Sheet::Thin(t) => t.min(2.0 * h),
                        _ => 2.0 * h,
                    };
                    p.settled |= thin >= 2.0 * h;
                    p.open.push((cell, CROSSED | if seam { SEAM } else { 0 }, thin));
                }
                (focused, held == Sheet::One)
            }
        }
    };
    if n == 1 {
        return;
    }
    let (t, m) = (focused.as_ref().unwrap_or(tape), n / 2);
    for k in 0..8 {
        unshown(t, ctx, grid, pts, vals, crossed, [c[0] + (k & 1) * m, c[1] + (k >> 1 & 1) * m, c[2] + (k >> 2 & 1) * m], m, drawn, untry, p);
    }
}

/// `tape`'s range over the box `b`, and `tape` with each min, max and
/// select whose side it decides there taken (None: it decides none): the
/// same field within `b` (`read` counts the instructions run, its
/// simplifying's too).
fn focus(tape: &Tape, ctx: &Ctx, b: [I; 3], read: &mut usize) -> (I, Option<Tape>) {
    *read += tape.work(Pass::Box);
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    let r = eval_interval(tape, ctx, b, Some(&mut ch));
    let simple = ch.iter().any(|&x| x != ink_core::tape::BOTH).then(|| tape.simplify(&ch));
    *read += simple.as_ref().map_or(0, |s| tape.simplify_work(s));
    (r, simple)
}

/// What the field over `lo..hi`, round a cell the surface crosses (the box
/// its vertex may be placed in, the cell grown by half a cell), holds besides what the mesh draws
/// there (Tape::sheet). One sheet: it changes sign at most once along one
/// direction, so no piece lies wholly inside, and that sheet and the
/// cell's vertex both lie in the box, within its diagonal (2 sqrt 3 h) of
/// each other: drawn. Else how thick what it may hide is, where a seam's
/// reach bounds it. In a cell the surface does not cross no vertex is
/// placed: a sheet there (a spike between samples) is left out, so this is
/// no proof there. Checked along each axis (a face square to one is
/// constant along the others), then along the middle of the field's
/// gradient's bounds (between the faces at an edge) and along its
/// gradient at the box's middle (any direction would do): where a select
/// may switch (a seam), as its reach; elsewhere over the box, and else
/// over it split where pieces meet (`split`). The best found, and whether
/// there is a seam. `read` counts the instructions run.
pub(super) fn sheet(tape: &Tape, ctx: &Ctx, lo: V3, hi: V3, read: &mut usize) -> (Sheet, bool) {
    *read += 4 * tape.work(Pass::Box);
    let b: [I; 3] = std::array::from_fn(|a| I::new(lo[a], hi[a]));
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    let g = ink_core::eval::eval_interval_grad_every(tape, ctx, b);
    let mid: V3 = std::array::from_fn(|a| {
        let i = g[g.len() - 1].d[a];
        0.5 * (i.lo + i.hi)
    });
    let centre: V3 = std::array::from_fn(|a| 0.5 * (lo[a] + hi[a]));
    let at = eval_grads(tape, ctx, &[centre])[0];
    let unit = |v: V3| {
        let n = norm(v);
        (n > 0.0 && n.is_finite()).then(|| v.map(|x| x / n))
    };
    let sheets = tape.sheets(&ch, &g);
    let seam = sheets.switches();
    let dirs = [unit(mid), unit([at[1], at[2], at[3]])];
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    *read += 3 * sheets.reads();
    if axes.iter().any(|&d| sheets.sheet(None, d) == Sheet::One) {
        return (Sheet::One, seam);
    }
    if seam {
        // A seam's reach over the box, from its bounds there and the
        // field's at its centre. (Over its eighths too, one direction for
        // all, it proved a third of a duct's seam cells drawn, 3% of a
        // sample's cells, for a fifth more time.)
        *read += 3 * tape.work(Pass::Box);
        let point = ink_core::eval::eval_interval_grad_every(tape, ctx, centre.map(I::point));
        let half: V3 = std::array::from_fn(|a| 0.5 * (hi[a] - lo[a]));
        let mut best = Sheet::Unknown;
        for d in dirs.into_iter().flatten() {
            *read += sheets.reads();
            let held = sheets.sheet(Some((&point, half)), d);
            if held == Sheet::One {
                return (Sheet::One, seam);
            }
            best = least(best, held);
        }
        return (best, seam);
    }
    let dirs: Vec<V3> = dirs.into_iter().flatten().chain(axes).collect();
    *read += 2 * sheets.reads();
    if dirs[..dirs.len() - 3].iter().any(|&d| sheets.sheet(None, d) == Sheet::One) || {
        // Split where pieces meet as the field simplified over the box
        // would be (its levels, not those of pieces decided away there).
        let ls = ink_core::sign::levels_in(tape, &ch);
        split(tape, ctx, b, &dirs, ink_core::sign::SPLITS, &ls, read).iter().any(|&w| w != 0)
    } {
        return (Sheet::One, seam);
    }
    (Sheet::Unknown, seam)
}

/// Of two sound accounts of a region (Tape::sheet), the stronger.
fn least(a: Sheet, b: Sheet) -> Sheet {
    match (a, b) {
        (Sheet::One, _) | (_, Sheet::One) => Sheet::One,
        (Sheet::Thin(x), Sheet::Thin(y)) => Sheet::Thin(x.min(y)),
        (Sheet::Unknown, x) | (x, Sheet::Unknown) => x,
    }
}

/// The field's way of changing sign along each unit d of `dirs` over box b
/// where it is one sheet there (Tape::sheet_way: 1, -1, 2; else 0,
/// unknown), or where that is unknown, over b split by a plane where
/// pieces meet (`split`), up to `depth` splits. The box's bounds are found
/// once for every direction.
fn ways(tape: &Tape, ctx: &Ctx, b: [I; 3], dirs: &[V3], depth: usize, ls: &[(u32, f64)], read: &mut usize) -> Vec<i8> {
    *read += 4 * tape.work(Pass::Box);
    let mut ch = vec![ink_core::tape::BOTH; tape.len()];
    eval_interval(tape, ctx, b, Some(&mut ch));
    let g = ink_core::eval::eval_interval_grad_every(tape, ctx, b);
    let sheets = tape.sheets(&ch, &g);
    let mut w: Vec<i8> = dirs
        .iter()
        .map(|&d| {
            *read += sheets.reads();
            match sheets.way(None, d) {
                (Sheet::One, w) => w,
                _ => 0,
            }
        })
        .collect();
    let open: Vec<usize> = (0..dirs.len()).filter(|&k| w[k] == 0).collect();
    if depth > 0 && !open.is_empty() {
        let found = split(tape, ctx, b, &open.iter().map(|&k| dirs[k]).collect::<Vec<_>>(), depth, ls, read);
        open.iter().zip(found).for_each(|(&k, f)| w[k] = f);
    }
    w
}

/// Along each unit d of `dirs`, the field's way over box b split by a
/// plane where pieces meet (`ls`: ink_core::sign::levels of the field
/// there, the first a coordinate's that crosses b; 0 where none does): one
/// half air (ink_core::sign::air), the other changing sign once along d,
/// where each line along d leaves it into the air after going from inside
/// to not, or enters it from the air before going from not to inside (the
/// air adds no change of sign); or, d along the plane (no line crosses
/// it), the halves changing sign the same way (or one keeping one sign).
/// Each half's bounds are found once, for every direction.
pub(super) fn split(tape: &Tape, ctx: &Ctx, b: [I; 3], dirs: &[V3], depth: usize, ls: &[(u32, f64)], read: &mut usize) -> Vec<i8> {
    let plane = ls.iter().find_map(|&(i, c)| ink_core::sign::axis(tape, i as usize).filter(|&k| b[k].lo < c && c < b[k].hi).map(|k| (k, c)));
    let Some((k, c)) = plane else { return vec![0; dirs.len()] };
    let (mut lower, mut upper) = (b, b);
    lower[k] = I::new(b[k].lo, c);
    upper[k] = I::new(c, b[k].hi);
    // Each half's ways and whether it is air, found when first asked.
    let mut halves: [Option<Vec<i8>>; 2] = [None, None];
    let mut airs: [Option<bool>; 2] = [None, None];
    let mut half = |h: usize, read: &mut usize| -> Vec<i8> { halves[h].get_or_insert_with(|| ways(tape, ctx, [lower, upper][h], dirs, depth - 1, ls, read)).clone() };
    let mut out = vec![0; dirs.len()];
    for (q, d) in dirs.iter().enumerate() {
        // d along the plane: each line along d stays in one half, so the
        // halves need only agree (one changing sign throughout one way,
        // the other too or keeping one sign: the plane between them is
        // then a face of the sheet, as at an edge).
        if d[k] == 0.0 {
            let (wl, wu) = (half(0, read)[q], half(1, read)[q]);
            let agreed = match (wl, wu) {
                (0, _) | (_, 0) => 0,
                (2, w) | (w, 2) => w,
                (a, b) if a == b => a,
                _ => 0,
            };
            if agreed != 0 {
                out[q] = agreed;
                continue;
            }
        }
        // s: the side of the plane the air is on.
        for (air, other, s) in [(0, 1, -1.0), (1, 0, 1.0)] {
            let is_air = *airs[air].get_or_insert_with(|| {
                *read += 2 * tape.work(Pass::Box);
                ink_core::sign::air(tape, ctx, [lower, upper][air])
            });
            if !is_air {
                continue;
            }
            let toward = d[k] * s;
            out[q] = match half(other, read)[q] {
                w @ (1 | -1) if toward == 0.0 || f64::from(w) * toward > 0.0 => w,
                // Inside or air throughout: its surface, if any, is the plane.
                2 if toward == 0.0 => 2,
                2 => {
                    if toward > 0.0 {
                        1
                    } else {
                        -1
                    }
                }
                _ => 0,
            };
            break;
        }
    }
    out
}

/// The bricks' crossed cells left untried (`Proofs::untried`) beside an
/// open cell the surface does not cross that the search reads (the
/// `first` in brick order, `searched_cells`), tried: a box of theirs shown
/// drawn may hold what it samples. Each brick's again from its tape, its
/// blocks focused but not tried again.
pub(super) fn try_untried(bricks: &mut [Brick], ctx: &Ctx, grid: &Grid, ran: &Ran, first: usize) {
    let searched: std::collections::HashSet<[i64; 3]> =
        bricks.iter().flat_map(|b| b.open.iter().filter(|o| o.1 & CROSSED == 0).map(move |o| cell_of(b.i0, o.0))).take(first).collect();
    let read = |c: [i64; 3]| (0..27).any(|q| q != 13 && searched.contains(&[c[0] + q % 3 - 1, c[1] + q / 3 % 3 - 1, c[2] + q / 9 - 1]));
    bricks.par_iter_mut().filter(|b| !b.untried.is_empty()).for_each(|b| {
        let (mut only, mut blocks) = (vec![false; (B * B * B) as usize], vec![false; (B * B * B) as usize]);
        let mut any = false;
        for &(k, tried) in &b.untried {
            if read(cell_of(b.i0, k)) {
                (only[k as usize], blocks[k as usize], any) = (true, tried, true);
            }
        }
        let Some(tape) = b.tape.clone().filter(|_| any) else { return };
        let pts: Vec<V3> = (0..N * N * N).map(|q| grid.pos([b.i0[0] + (q % N) as i64, b.i0[1] + (q / N % N) as i64, b.i0[2] + (q / (N * N)) as i64])).collect();
        let crossed: Vec<bool> = (0..(B * B * B) as usize)
            .map(|k| {
                let (i, j, l) = (k % B as usize, k / B as usize % B as usize, k / (B * B) as usize);
                let at = |c: usize| b.vals[vidx(i + (c & 1), j + (c >> 1 & 1), l + (c >> 2 & 1))];
                // A cell not meshed (in a region drawn whole) has corners not sampled.
                (0..8).all(|c| !at(c).is_nan()) && (1..8).any(|c| inside(at(c)) != inside(at(0)))
            })
            .collect();
        let mut p = Proofs { only: Some((&only, &blocks)), ..Default::default() };
        unshown(&tape, ctx, grid, &pts, &b.vals, &crossed, [0; 3], B as usize, false, false, &mut p);
        ran.tapes.fetch_add(p.tapes, std::sync::atomic::Ordering::Relaxed);
        ran.boxes.fetch_add(p.read, std::sync::atomic::Ordering::Relaxed);
        b.untried.retain(|&(k, _)| !only[k as usize]);
        let mut all: Vec<((u16, u8), f64)> = b.open.drain(..).zip(b.thin.drain(..)).chain(p.open.into_iter().map(|(k, state, t)| ((k, state), t))).collect();
        all.sort_unstable_by_key(|o| o.0.0);
        (b.open, b.thin) = all.into_iter().unzip();
    });
}
