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

/// How many times `undecided` halves the cube about a sample on 0 looking
/// for air.
const AIR_DEPTH: u32 = 2;

/// Whether some part of box b (halved up to `depth` times) holds no point
/// inside: air there, shown by bounds.
fn air_in(tape: &Tape, ctx: &Ctx, b: [I; 3], depth: u32) -> bool {
    let r = crate::eval::eval_interval(tape, ctx, b, None);
    if r.none_inside() {
        return true;
    }
    if r.all_inside() || depth == 0 {
        return false;
    }
    (0..8).any(|o| air_in(tape, ctx, halves(b, o), depth - 1))
}

/// Box b's half `o` (bit k: the upper half along axis k).
fn halves(b: [I; 3], o: usize) -> [I; 3] {
    std::array::from_fn(|k| if o >> k & 1 == 0 { I::new(b[k].lo, b[k].mid()) } else { I::new(b[k].mid(), b[k].hi) })
}

/// The samples on 0 (`vals`, the field at `pts`) that nothing near them
/// decides. A sample on 0 is not inside, whatever is round it: drawn as
/// air, never filled. Neither bounds nor samples can show that a 0 has no
/// thickness (faces meeting flush, or air the field reads as 0 however
/// thin, read alike), so none is read as material. It is decided where a
/// few billionths of its size away the field shows air beside it: a
/// surface crossing it (a point inside one way along the field's slope and
/// one not the other, each decided by its own value), or a part of the
/// cube about it with no point inside by bounds (halved twice: air read as
/// 0). Else it is returned, for the caller to report: a wall may be drawn
/// there where pieces meet flush. `within`: the box `tape` holds over (a
/// brick's simplified tape) and the whole field: a side sample past the
/// box is read by the whole field, the cube clipped to the box.
pub fn undecided(tape: &Tape, ctx: &Ctx, pts: &[V3], vals: &[f64], within: Option<([I; 3], &Tape)>) -> Vec<V3> {
    let zeros: Vec<V3> = pts.iter().zip(vals).filter(|(_, v)| **v == 0.0).map(|(p, _)| *p).collect();
    if zeros.is_empty() {
        return Vec::new();
    }
    let grads = crate::eval::eval_grads(tape, ctx, &zeros);
    let mut sides = Vec::with_capacity(2 * zeros.len());
    for (p, g) in zeros.iter().zip(&grads) {
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
    for (i, p) in zeros.iter().enumerate() {
        if inside(sv[2 * i]) != inside(sv[2 * i + 1]) {
            continue;
        }
        let near = 4.0 * least(*p);
        let mut cube: [I; 3] = std::array::from_fn(|k| I::new(p[k] - near, p[k] + near));
        if let Some((w, _)) = within {
            cube = std::array::from_fn(|k| I::new(cube[k].lo.max(w[k].lo), cube[k].hi.min(w[k].hi)));
        }
        if !(cube.iter().all(|x| x.lo <= x.hi) && air_in(tape, ctx, cube, AIR_DEPTH)) {
            open.push(*p);
        }
    }
    open
}

/// Where a tape's value may be split by sign: a coordinate less a
/// constant (a plane square to an axis: the box is cut there), and a value
/// it reads more than once (pieces meeting on a curved face each read it:
/// its sign is fixed for every read). Coordinates and constants aside.
enum Cut {
    Plane(usize, f64),
    Value,
}

fn cuts(tape: &Tape) -> Vec<(usize, Cut)> {
    let mut writer = vec![usize::MAX; tape.nslots];
    let mut reads = vec![0usize; tape.code.len()];
    let mut out = Vec::new();
    let axis = |op: &Op| match op {
        Op::Var(crate::graph::Var::X) => Some(0),
        Op::Var(crate::graph::Var::Y) => Some(1),
        Op::Var(crate::graph::Var::Z) => Some(2),
        _ => None,
    };
    for (i, ins) in tape.code.iter().enumerate() {
        let (args, n) = ins.args();
        for &a in &args[..n] {
            let w = writer[a as usize];
            if w != usize::MAX {
                reads[w] += 1;
            }
        }
        // v - c (or c - v), and |v - c1| - c2 (a slab, as std's box writes
        // its faces): planes v = c, v = c1 +- c2 where that sum is exact.
        if let Op::Bin(crate::graph::BinOp::Sub) = ins.op {
            let (wa, wb) = (writer[ins.a as usize], writer[ins.b as usize]);
            if wa != usize::MAX && wb != usize::MAX {
                let (oa, ob) = (&tape.code[wa].op, &tape.code[wb].op);
                match (axis(oa), ob, axis(ob), oa) {
                    (Some(k), Op::Const(c), _, _) | (_, _, Some(k), Op::Const(c)) => out.push((i, Cut::Plane(k, *c))),
                    _ => {}
                }
                if let (Op::Un(crate::graph::UnOp::Abs), Op::Const(c2)) = (oa, ob) {
                    let inner = writer[tape.code[wa].a as usize];
                    if inner != usize::MAX && matches!(tape.code[inner].op, Op::Bin(crate::graph::BinOp::Sub)) {
                        let (wv, wc) = (writer[tape.code[inner].a as usize], writer[tape.code[inner].b as usize]);
                        if wv != usize::MAX
                            && wc != usize::MAX
                            && let (Some(k), Op::Const(c1)) = (axis(&tape.code[wv].op), &tape.code[wc].op)
                        {
                            for at in [c1 - c2, c1 + c2] {
                                // Exact: the face is that float.
                                if (at - c1).abs() == c2.abs() {
                                    out.push((i, Cut::Plane(k, at)));
                                }
                            }
                        }
                    }
                }
            }
        }
        writer[ins.out as usize] = i;
    }
    let planes: Vec<usize> = out.iter().map(|c| c.0).collect();
    out.extend((0..tape.code.len()).filter(|&i| reads[i] >= 2 && !planes.contains(&i) && !matches!(tape.code[i].op, Op::Const(_) | Op::Var(_))).map(|i| (i, Cut::Value)));
    out
}

/// Whether `test` holds of the tape's bounds over every point of box b,
/// shown by bounds over the box, else over parts of it split where a plane
/// square to an axis crosses it (the box cut there) or where a value it
/// reads more than once changes sign (each part bounded with that value's
/// sign known wherever it is read, `eval_interval_clamped`), up to SPLITS
/// deep. Where pieces meet on a face, or a cut ends flush with one, the
/// field is exactly 0 across it, and the bounds of the box alone (each read
/// taken apart) can never show which side of 0 it keeps; split by the face,
/// each part is within one piece or the other.
pub fn holds_by_cuts(tape: &Tape, ctx: &Ctx, b: [I; 3], test: &dyn Fn(I) -> bool) -> bool {
    const SPLITS: usize = 3;
    #[allow(clippy::too_many_arguments)]
    fn all(tape: &Tape, ctx: &Ctx, b: [I; 3], cs: &[(usize, Cut)], clamps: &mut Vec<(usize, I)>, depth: usize, test: &dyn Fn(I) -> bool, vals: &mut Vec<I>) -> bool {
        let Some(r) = crate::eval::eval_interval_clamped(tape, ctx, b, clamps, vals) else { return true };
        if test(r) {
            return true;
        }
        if depth == SPLITS {
            return false;
        }
        let straddles = |i: usize, vals: &[I]| vals[i].lo < 0.0 && vals[i].hi > 0.0;
        let Some((i, cut)) =
            cs.iter().find(|(i, c)| straddles(*i, vals) && !clamps.iter().any(|k| k.0 == *i) && !matches!(c, Cut::Plane(k, at) if !(b[*k].lo < *at && *at < b[*k].hi)))
        else {
            return false;
        };
        match *cut {
            Cut::Plane(k, at) => {
                let (mut lo, mut hi) = (b, b);
                lo[k] = I::new(b[k].lo, at);
                hi[k] = I::new(at, b[k].hi);
                all(tape, ctx, lo, cs, clamps, depth + 1, test, vals) && all(tape, ctx, hi, cs, clamps, depth + 1, test, vals)
            }
            Cut::Value => {
                let v = vals[*i];
                for part in [I { lo: v.lo, hi: 0.0 }, I { lo: 0.0, hi: v.hi }] {
                    clamps.push((*i, part));
                    let ok = all(tape, ctx, b, cs, clamps, depth + 1, test, vals);
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
    let cs = cuts(tape);
    !cs.is_empty() && all(tape, ctx, b, &cs, &mut Vec::new(), 0, test, &mut vals)
}

/// Where a surface crosses the line through `p` along unit `n` within
/// `tol` of p: a point inside and one not, either side of it, `tol` or
/// less apart, found from pairs at distances from `d` down to a billionth
/// of p's size, each pair read as samples that far apart
/// (`undecided`) and halved toward its crossing. None where the only
/// crossings are farther than `tol` (a 0 with one sign round it, a real
/// surface beside it, is no surface).
pub fn crossing(tape: &Tape, ctx: &Ctx, p: V3, n: V3, d: f64, tol: f64) -> Option<V3> {
    let at = |t: f64| -> V3 { std::array::from_fn(|k| p[k] + t * n[k]) };
    // A side's membership; None where a sample on 0 there is undecided
    // (`undecided`): no crossing is shown across it.
    let side = |t: f64| {
        let q = at(t);
        let mut v = Vec::new();
        eval_points(tape, ctx, &[q], &mut v);
        undecided(tape, ctx, &[q], &v, None).is_empty().then(|| inside(v[0]))
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
    fn a_0_is_never_material_and_is_decided_only_beside_air() {
        // 0 on x = 0 inside the material: undecided; in a slit 1.5 wide
        // read as 0: air, decided.
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let ax = g.un(UnOp::Abs, x);
        let nax = g.neg(ax);
        let w = g.constant(0.75);
        let cut = g.sub(w, ax);
        let zero = g.constant(0.0);
        let slit = g.min(zero, cut);
        let ctx = Ctx::new(0.0, vec![]);
        for (f, decided) in [(nax, false), (slit, true)] {
            let t = Tape::new(&g, f);
            assert_eq!(undecided(&t, &ctx, &[[0.0; 3]], &[0.0], None).is_empty(), decided);
        }
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
