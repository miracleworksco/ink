//! Which side of a shape a point is on: the one rule every part of ink
//! reads (docs/CORE.md). A point is in a shape where the shape's field is
//! below 0 there; 0 is not inside. A surface lies between a point inside
//! and one that is not, so a 0 the field keeps one sign around (faces
//! meeting flush inside the material, a cut ending flush with a face) is
//! no surface, and a region where the field is 0 is air, however thin.

use crate::V3;
use crate::eval::{Ctx, eval_points};
use crate::interval::I;
use crate::tape::{Op, Tape};

/// Whether a field's value puts its point inside.
#[inline]
pub fn inside(f: f64) -> bool {
    f < 0.0
}

impl I {
    /// Bounds over a region show every point of it inside.
    #[inline]
    pub fn all_inside(self) -> bool {
        self.hi < 0.0
    }

    /// Bounds over a region show no point of it inside.
    #[inline]
    pub fn none_inside(self) -> bool {
        self.lo >= 0.0
    }

    /// Bounds over a region leave room for a point inside.
    #[inline]
    pub fn may_be_inside(self) -> bool {
        !self.none_inside()
    }

    /// Bounds over a closed region show no point inside nor on the surface
    /// (a surface point, the limit of points inside but not inside itself,
    /// reads 0): no surface meets it, even on its faces.
    #[inline]
    pub fn clear_of_surface(self) -> bool {
        self.lo > 0.0
    }
}

/// The least distance a point's neighbours are read at: a billionth of its
/// size (nearer, the field's rounding may give a 0 either sign).
fn least(p: V3) -> f64 {
    1e-9 * p.iter().fold(1.0f64, |m, x| m.max(x.abs()))
}

/// How many times a sample on 0's cube is halved looking for air.
const AIR_DEPTH: u32 = 2;

/// Whether some part of box b (halved up to `depth` times) holds no point
/// inside: air there, shown by bounds (`air`).
fn air_in(tape: &Tape, ctx: &Ctx, b: [I; 3], depth: u32) -> bool {
    if air(tape, ctx, b) {
        return true;
    }
    if depth == 0 || crate::eval::eval_interval(tape, ctx, b, None).all_inside() {
        return false;
    }
    (0..8).any(|o| air_in(tape, ctx, halves(b, o), depth - 1))
}

/// Box b's half `o` (bit k: the upper half along axis k).
fn halves(b: [I; 3], o: usize) -> [I; 3] {
    std::array::from_fn(|k| if o >> k & 1 == 0 { I::new(b[k].lo, b[k].mid()) } else { I::new(b[k].mid(), b[k].hi) })
}

/// The value a sample on 0 is drawn with where it is material (`drawn`):
/// the least below 0.
pub const SOLID: f64 = -f64::MIN_POSITIVE;

/// Samples on 0 (`vals`, the field at `pts`) as drawn, and those nothing
/// near decides, returned for the caller to report. A sample on 0 is not
/// inside; a surface is where the field changes sign, so where the 0 is
/// one the field keeps one sign around it is none. Each is read a few
/// billionths of its size about it, the finest length ink speaks of: a
/// surface crossing it (a point inside one way along the field's slope and
/// one not the other, each decided by its own value), or a part of the
/// cube about it with no point inside by bounds (`air`, halved twice: air
/// read as 0), leave it air; the cube solid by bounds (`solid`: inside but
/// on a 0 of no volume, pieces meeting flush) draws it as material (its
/// value set to SOLID). Else it is drawn as air and returned: neither
/// samples nor bounds show which. `within`: the box `tape` holds over (a
/// brick's simplified tape) and the whole field: a side sample past the
/// box is read by the whole field, the cube clipped to the box.
pub fn drawn(tape: &Tape, ctx: &Ctx, pts: &[V3], vals: &mut [f64], within: Option<([I; 3], &Tape)>) -> Vec<V3> {
    let zeros: Vec<usize> = (0..vals.len()).filter(|&i| vals[i] == 0.0).collect();
    if zeros.is_empty() {
        return Vec::new();
    }
    let at: Vec<V3> = zeros.iter().map(|&i| pts[i]).collect();
    let grads = crate::eval::eval_grads(tape, ctx, &at);
    let mut sides = Vec::with_capacity(2 * zeros.len());
    for (p, g) in at.iter().zip(&grads) {
        let n = (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
        let d = 4.0 * least(*p) / if n > 0.0 && n.is_finite() { n } else { f64::INFINITY };
        sides.extend([-d, d].map(|t| -> V3 { std::array::from_fn(|k| p[k] + t * g[k + 1]) }));
    }
    let mut sv = Vec::new();
    eval_points(tape, ctx, &sides, &mut sv);
    if let Some((w, whole)) = within {
        let past: Vec<usize> = (0..sides.len()).filter(|&i| (0..3).any(|k| !w[k].contains(sides[i][k]))).collect();
        if !past.is_empty() {
            let mut pv = Vec::new();
            eval_points(whole, ctx, &past.iter().map(|&i| sides[i]).collect::<Vec<_>>(), &mut pv);
            past.iter().zip(pv).for_each(|(&i, v)| sv[i] = v);
        }
    }
    let mut open = Vec::new();
    for (z, (&i, p)) in zeros.iter().zip(&at).enumerate() {
        if inside(sv[2 * z]) != inside(sv[2 * z + 1]) {
            continue;
        }
        let near = 4.0 * least(*p);
        let mut cube: [I; 3] = std::array::from_fn(|k| I::new(p[k] - near, p[k] + near));
        if let Some((w, _)) = within {
            cube = std::array::from_fn(|k| I::new(cube[k].lo.max(w[k].lo), cube[k].hi.min(w[k].hi)));
        }
        if !cube.iter().all(|x| x.lo <= x.hi) {
            open.push(*p);
        } else if solid(tape, ctx, cube) {
            vals[i] = SOLID;
        } else if !air_in(tape, ctx, cube, AIR_DEPTH) {
            open.push(*p);
        }
    }
    open
}

/// Whether bounds show every point of box b inside, or inside but on a 0
/// of no volume (pieces meeting flush: the field at most 0 throughout, by
/// `holds_by_cuts`, and its 0 holding no ball, Tape::zeros_thin): no
/// surface meets the box, since a surface needs air beside material, and
/// there the material is dense. Drawn, it is material throughout.
pub fn solid(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> bool {
    let r = crate::eval::eval_interval(tape, ctx, b, None);
    r.all_inside() || r.may_be_inside() && holds_by_cuts(tape, ctx, b, &|r: I| r.hi <= 0.0) && tape.zeros_thin(&crate::eval::eval_interval_grad_every(tape, ctx, b))
}

/// Whether bounds show no point of box b inside, by `holds_by_cuts`: a 0
/// between pieces that meet flush in air (a cut ending flush with a face)
/// is decided across it.
pub fn air(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> bool {
    holds_by_cuts(tape, ctx, b, &|r: I| r.none_inside())
}

/// Where pieces of a tape may meet flush, as the levels a value is split
/// at: (instruction, c), a value the tape reads more than once (or a
/// coordinate) and a constant it is compared with, v - c (c - v, v + -c),
/// or |v - c1| - c2 (a slab, as std's box writes its faces: c1 +- c2 where
/// that sum is exact), or 0 (pieces meeting on a face each read it). An
/// instruction is the same in the SSA form and the slot code.
pub fn levels(tape: &Tape) -> &[(u32, f64)] {
    tape.levels.get_or_init(|| found_levels(tape, None))
}

/// `levels` of the tape simplified by `choices` (Tape::simplify), in the
/// tape's own numbering: what a region's proofs split at, read from the
/// unsimplified tape.
pub fn levels_in(tape: &Tape, choices: &[u8]) -> Vec<(u32, f64)> {
    found_levels(tape, Some(choices))
}

/// `levels`, through the branches `choices` pick (each value read as the
/// simplified tape reads it, `Tape::picked`; the values it no longer
/// reads left out).
fn found_levels(tape: &Tape, choices: Option<&[u8]>) -> Vec<(u32, f64)> {
    use crate::graph::{BinOp, UnOp, Var};
    let ssa = &tape.ssa;
    let live = choices.map(|c| tape.live(c));
    let kept = |i: usize| live.as_ref().is_none_or(|l| l[i]) && choices.is_none_or(|c| tape.picked(c, i) == i);
    let at = |j: u32| choices.map_or(j, |c| tape.picked(c, j as usize) as u32);
    let mut reads = vec![0usize; ssa.len()];
    let constant = |j: u32| match ssa[at(j) as usize].op {
        Op::Const(c) => Some(c),
        _ => None,
    };
    let mut found: Vec<(u32, f64)> = Vec::new();
    for (i, ins) in ssa.iter().enumerate() {
        if !kept(i) {
            continue;
        }
        let (args, n) = Tape::reads(ssa, i);
        for &x in &args[..n] {
            reads[at(x) as usize] += 1;
        }
        match ins.op {
            Op::Bin(BinOp::Sub) => {
                if let Some(c) = constant(ins.b) {
                    found.push((at(ins.a), c));
                    // |v - c1| - c2.
                    let abs = &ssa[at(ins.a) as usize];
                    if let Op::Un(UnOp::Abs) = abs.op
                        && let inner = &ssa[at(abs.a) as usize]
                        && let (Op::Bin(BinOp::Sub), Some(c1)) = (inner.op, constant(inner.b))
                    {
                        for v in [c1 - c, c1 + c] {
                            // Exact: the face is that float.
                            if (v - c1).abs() == c.abs() {
                                found.push((at(inner.a), v));
                            }
                        }
                    }
                } else if let Some(c) = constant(ins.a) {
                    found.push((at(ins.b), c));
                }
            }
            Op::Bin(BinOp::Add) => {
                for (v, k) in [(ins.a, ins.b), (ins.b, ins.a)] {
                    if let Some(c) = constant(k) {
                        found.push((at(v), -c));
                    }
                }
            }
            _ => {}
        }
    }
    let shared = |i: usize| match ssa[i].op {
        Op::Var(Var::X | Var::Y | Var::Z) => true,
        Op::Const(_) | Op::Var(_) => false,
        _ => reads[i] >= 2,
    };
    found.extend((0..ssa.len()).filter(|&i| kept(i) && shared(i)).map(|i| (i as u32, 0.0)));
    let mut out: Vec<(u32, f64)> = Vec::new();
    for (i, c) in found {
        if shared(i as usize) && c.is_finite() && !out.contains(&(i, c)) {
            out.push((i, c));
        }
    }
    out
}

/// The axis of a coordinate instruction `i` reads.
pub fn axis(tape: &Tape, i: usize) -> Option<usize> {
    match tape.ssa[i].op {
        Op::Var(crate::graph::Var::X) => Some(0),
        Op::Var(crate::graph::Var::Y) => Some(1),
        Op::Var(crate::graph::Var::Z) => Some(2),
        _ => None,
    }
}

/// How many times a box is split where pieces meet (`holds_by_cuts`).
pub const SPLITS: usize = 3;

/// Whether `test` holds of the tape's bounds over every point of box b,
/// shown by bounds over the box, else over parts of it split where a value
/// that pieces each read crosses a level (`levels`): each part bounded with
/// that value on one side of it wherever it is read (a coordinate's: the
/// box cut there; else `eval_interval_clamped`), up to SPLITS deep. Where
/// pieces meet on a face, or a cut ends flush with one, the field is
/// exactly 0 across it, and the bounds of the box alone (each read taken
/// apart) can never show which side of 0 it keeps; split by the face, each
/// part is within one piece or the other.
pub fn holds_by_cuts(tape: &Tape, ctx: &Ctx, b: [I; 3], test: &dyn Fn(I) -> bool) -> bool {
    #[allow(clippy::too_many_arguments)]
    fn all(tape: &Tape, ctx: &Ctx, b: [I; 3], ls: &[(u32, f64)], clamps: &mut Vec<(usize, I)>, depth: usize, test: &dyn Fn(I) -> bool, vals: &mut Vec<I>) -> bool {
        let Some(r) = crate::eval::eval_interval_clamped(tape, ctx, b, clamps, vals) else { return true };
        if test(r) {
            return true;
        }
        if depth == SPLITS {
            return false;
        }
        let Some(&(i, c)) = ls.iter().find(|&&(i, c)| vals[i as usize].lo < c && c < vals[i as usize].hi) else { return false };
        let i = i as usize;
        match axis(tape, i) {
            Some(k) => {
                let (mut lo, mut hi) = (b, b);
                lo[k] = I::new(b[k].lo, c);
                hi[k] = I::new(c, b[k].hi);
                all(tape, ctx, lo, ls, clamps, depth + 1, test, vals) && all(tape, ctx, hi, ls, clamps, depth + 1, test, vals)
            }
            None => {
                let v = vals[i];
                for part in [I { lo: v.lo, hi: c }, I { lo: c, hi: v.hi }] {
                    clamps.push((i, part));
                    let ok = all(tape, ctx, b, ls, clamps, depth + 1, test, vals);
                    clamps.pop();
                    if !ok {
                        return false;
                    }
                }
                true
            }
        }
    }
    let mut vals = Vec::new();
    if test(crate::eval::eval_interval_clamped(tape, ctx, b, &[], &mut vals).unwrap_or(I::point(0.0))) {
        return true;
    }
    let ls = levels(tape);
    !ls.is_empty() && all(tape, ctx, b, ls, &mut Vec::new(), 0, test, &mut vals)
}

/// Where a surface crosses the line through `p` along unit `n` within
/// `tol` of p: a point inside and one not, either side of it, `tol` or
/// less apart, found from pairs at distances from `d` down to a billionth
/// of p's size, each pair read as samples that far apart
/// (`drawn`) and halved toward its crossing. None where the only
/// crossings are farther than `tol` (a 0 with one sign round it, a real
/// surface beside it, is no surface).
pub fn crossing(tape: &Tape, ctx: &Ctx, p: V3, n: V3, d: f64, tol: f64) -> Option<V3> {
    let at = |t: f64| -> V3 { std::array::from_fn(|k| p[k] + t * n[k]) };
    // A side's membership as drawn; None where a sample on 0 there is
    // undecided (`drawn`): no crossing is shown across it.
    let side = |t: f64| {
        let q = at(t);
        let mut v = Vec::new();
        eval_points(tape, ctx, &[q], &mut v);
        drawn(tape, ctx, &[q], &mut v, None).is_empty().then(|| inside(v[0]))
    };
    let least = least(p);
    let mut s = d;
    loop {
        let (mut a, mut b) = (-s, s);
        if let (Some(ia), Some(ib)) = (side(a), side(b))
            && ia != ib
        {
            // Halve toward a crossing; it counts where its bracket comes
            // within tol of p.
            let mut shown = true;
            while b - a > tol {
                let m = 0.5 * (a + b);
                if m == a || m == b {
                    break;
                }
                match side(m) {
                    Some(i) if i == ia => a = m,
                    Some(_) => b = m,
                    None => {
                        shown = false;
                        break;
                    }
                }
            }
            if shown && a <= tol && b >= -tol {
                return Some(at(0.5 * (a + b)));
            }
        }
        if s <= least {
            return None;
        }
        s = (0.25 * s).max(least);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{Graph, UnOp, Var};

    #[test]
    fn a_0_inside_the_material_is_drawn_as_material_and_one_beside_air_as_air() {
        // -|x|: 0 on x = 0 inside the material, of no volume; a slit 1.5
        // wide read as 0: air.
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let ax = g.un(UnOp::Abs, x);
        let nax = g.neg(ax);
        let w = g.constant(0.75);
        let cut = g.sub(w, ax);
        let zero = g.constant(0.0);
        let slit = g.min(zero, cut);
        let ctx = Ctx::new(0.0, vec![]);
        for (f, drawn_as) in [(nax, SOLID), (slit, 0.0)] {
            let t = Tape::new(&g, f);
            let mut v = [0.0];
            assert!(drawn(&t, &ctx, &[[0.0; 3]], &mut v, None).is_empty());
            assert_eq!(v[0], drawn_as);
        }
    }

    #[test]
    fn faces_meeting_flush_are_cut_where_they_meet_off_the_axes() {
        // A slab and a hole that ends flush with its face, both read along
        // u = (x + y) / sqrt 2: 0 across the mouth, air either side; the
        // halves of a slab meeting flush on u = 1: 0 inside the material.
        let mut g = Graph::new();
        let (x, y) = (g.var(Var::X), g.var(Var::Y));
        let s = g.add(x, y);
        let k = g.constant(std::f64::consts::FRAC_1_SQRT_2);
        let u = g.mul(s, k);
        let (one, two) = (g.constant(1.0), g.constant(2.0));
        let slab = g.sub(u, one);
        let hole = g.sub(one, u);
        let mouth = g.max(slab, hole);
        let (below, above) = (g.sub(u, one), g.sub(u, two));
        let above = g.max(above, hole);
        let joined = g.min(below, above);
        let ctx = Ctx::new(0.0, vec![]);
        let at = std::f64::consts::FRAC_1_SQRT_2;
        let b = [I::new(at - 0.1, at + 0.1), I::new(at - 0.1, at + 0.1), I::new(-1.0, 1.0)];
        let (m, j) = (Tape::new(&g, mouth), Tape::new(&g, joined));
        assert!(!crate::eval::eval_interval(&m, &ctx, b, None).none_inside(), "the box alone leaves it open");
        assert!(air(&m, &ctx, b));
        assert!(!crate::eval::eval_interval(&j, &ctx, b, None).all_inside());
        assert!(solid(&j, &ctx, b));
    }

    #[test]
    fn a_zero_beside_a_surface_is_no_crossing() {
        // 0 at x = 0, positive either side; a surface at x = -0.5.
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let ax = g.un(UnOp::Abs, x);
        let half = g.constant(0.5);
        let wall = g.add(x, half);
        let f = g.min(ax, wall);
        let t = Tape::new(&g, f);
        let ctx = Ctx::new(0.0, vec![]);
        assert_eq!(crossing(&t, &ctx, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 1e-8), None);
        let c = crossing(&t, &ctx, [-0.5, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1e-8).expect("the surface");
        assert!((c[0] + 0.5).abs() <= 1e-8, "{c:?}");
    }
}
