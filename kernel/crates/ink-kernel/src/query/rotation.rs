//! An enclosure for a field read in rotating coordinates. The warp has
//! already been substituted into the field. Find its two actual coordinate
//! expressions, so the proof also sees moves and fields combined with it.

use ink_core::V3;
use ink_core::eval::{Ctx, Region, eval_interval_clamped, eval_interval_linear};
use ink_core::graph::{BinOp, UnOp, Var};
use ink_core::interval::{I, add_o, mag, mul_o, sub_o};
use ink_core::tape::{Op, Tape};
use std::collections::HashMap;

struct Rotation {
    angle: usize,
    input: [usize; 2],
    output: [usize; 2],
}

/// cos(a) u + sin(a) v and cos(a) v - sin(a) u, with exactly the same
/// angle and two inputs. Look through commuted sums and products only.
/// Fixed operations use float arithmetic, so do not use the real-number
/// rotation identity for them.
fn rotations(t: &Tape) -> Vec<Rotation> {
    let product = |i: u32, op| {
        let p = t.ssa[i as usize];
        if t.fixed[i as usize] || !matches!(p.op, Op::Bin(BinOp::Mul)) {
            return None;
        }
        [(p.a, p.b), (p.b, p.a)].into_iter().find_map(|(trig, coord)| {
            let f = t.ssa[trig as usize];
            (matches!(f.op, Op::Un(o) if o == op) && !t.fixed[trig as usize]).then_some((f.a, coord))
        })
    };
    let mut sums = HashMap::new();
    let mut differences = Vec::new();
    for (i, ins) in t.ssa.iter().enumerate() {
        if t.fixed[i] || !matches!(ins.op, Op::Bin(BinOp::Add | BinOp::Sub)) {
            continue;
        }
        let terms = |a, b| {
            let (angle, u) = product(a, UnOp::Cos)?;
            let (other, v) = product(b, UnOp::Sin)?;
            (angle == other).then_some((angle, u, v))
        };
        if matches!(ins.op, Op::Bin(BinOp::Add)) {
            if let Some(key) = terms(ins.a, ins.b).or_else(|| terms(ins.b, ins.a)) {
                sums.insert(key, i);
            }
        } else if let Some((a, v, u)) = terms(ins.a, ins.b) {
            differences.push(((a, u, v), i));
        }
    }
    differences.into_iter().filter_map(|(key @ (a, u, v), j)| Some(Rotation { angle: a as usize, input: [u as usize, v as usize], output: [*sums.get(&key)?, j] })).collect()
}

/// Finite at every finite point, even on a slab extending to infinity.
/// An actual infinite angle (e.g. ln 0) cannot be treated as a rotation.
fn finite_values(t: &Tape, values: &[I]) -> Vec<bool> {
    let mut finite = Vec::with_capacity(values.len());
    for (i, ins) in t.ssa.iter().enumerate() {
        finite.push(
            values[i].lo.is_finite() && values[i].hi.is_finite()
                || !t.fixed[i]
                    && match ins.op {
                        Op::Var(Var::X | Var::Y | Var::Z) => true,
                        Op::Bin(BinOp::Add | BinOp::Sub | BinOp::Mul) => finite[ins.a as usize] && finite[ins.b as usize],
                        Op::Un(UnOp::Rad | UnOp::Deg | UnOp::Neg | UnOp::Abs) | Op::Tag(_) | Op::Real => finite[ins.a as usize],
                        _ => false,
                    },
        );
    }
    finite
}

/// Every material point must have |coordinate| < cap. Prove each side
/// over *all* space, including any disconnected material. Start with the
/// largest coordinate in the proposed window, then tighten the enclosure.
fn half_width(t: &Tape, ctx: &Ctx, coordinate: usize, limit: f64) -> Option<f64> {
    let mut scratch = Vec::new();
    let mut fits = |cap: f64| {
        [I::new(f64::NEG_INFINITY, -cap), I::new(cap, f64::INFINITY)]
            .into_iter()
            .all(|side| eval_interval_clamped(t, ctx, [I::ENTIRE; 3], &[(coordinate, side)], &mut scratch).is_none_or(|v| v.none_inside()))
    };
    let (mut lo, mut hi) = (0.0, limit.next_up());
    if !hi.is_finite() || !fits(hi) {
        return None;
    }
    for _ in 0..16 {
        let mid = I::new(lo, hi).mid();
        if fits(mid) { hi = mid } else { lo = mid }
    }
    Some(hi)
}

/// In the input plane, material lies inside the same cylinder as in the
/// rotated plane, since (cu + sv)² + (cv - su)² = u² + v². The field
/// itself supplies the axial extent (a twist about z leaves z unchanged).
struct Cylinder {
    xy: [Tape; 2],
    radius_squared: f64,
}

impl Cylinder {
    fn excludes(&self, ctx: &Ctx, region: &Region) -> bool {
        let [u, v] = self.xy.each_ref().map(|t| {
            let r = eval_interval_linear(t, ctx, region, None);
            I::point(r.lo.max(-r.hi).max(0.0))
        });
        add_o(mul_o(u, u), mul_o(v, v)).lo >= self.radius_squared
    }
}

/// Prove an entire outward cone empty. Divide its directions where it
/// straddles both the cylinder and the field's end planes (a tilted part).
/// Giving up leaves the original refusal. No finite distance cutoff.
fn exterior_empty(t: &Tape, ctx: &Ctx, cylinders: &[Cylinder], region: Region, depth: u8) -> bool {
    if cylinders.iter().any(|c| c.excludes(ctx, &region)) || eval_interval_linear(t, ctx, &region, None).none_inside() {
        return true;
    }
    let Region::Sector { c, r, mut q, s } = region else { unreachable!() };
    if depth == 0 {
        return false;
    }
    let k = (0..3).max_by(|&a, &b| (q[a].hi - q[a].lo).total_cmp(&(q[b].hi - q[b].lo))).unwrap();
    let (span, mid) = (q[k], q[k].mid());
    [I::new(span.lo, mid), I::new(mid, span.hi)].into_iter().all(|part| {
        q[k] = part;
        exterior_empty(t, ctx, cylinders, Region::Sector { c, r, q, s }, depth - 1)
    })
}

/// Grant the working window only if every point outside is proven air.
pub(super) fn inside_window(t: &Tape, ctx: &Ctx, (lo, hi): (V3, V3)) -> bool {
    let pairs = rotations(t);
    if pairs.is_empty() {
        return false;
    }
    let window: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let mut values = Vec::new();
    eval_interval_clamped(t, ctx, [I::ENTIRE; 3], &[], &mut values);
    let finite = finite_values(t, &values);
    eval_interval_clamped(t, ctx, window, &[], &mut values);
    let cylinders: Vec<_> = pairs
        .into_iter()
        .filter_map(|p| {
            if !finite[p.angle] || !p.input.iter().all(|&i| finite[i]) {
                return None;
            }
            let [a, b] = p.output;
            let w = I::point(half_width(t, ctx, a, mag(values[a]))?);
            let h = I::point(half_width(t, ctx, b, mag(values[b]))?);
            let radius_squared = add_o(mul_o(w, w), mul_o(h, h)).hi;
            radius_squared.is_finite().then(|| Cylinder { xy: p.input.map(|i| t.cone(i)), radius_squared })
        })
        .collect();
    if cylinders.is_empty() {
        return false;
    }
    // An inner box, even if the midpoint rounded. Its exterior covers
    // the window's entire exterior.
    let c = window.map(|w| w.mid());
    let r: V3 = std::array::from_fn(|k| sub_o(I::point(c[k]), I::point(lo[k])).lo.min(sub_o(I::point(hi[k]), I::point(c[k])).lo));
    r.into_iter().all(|h| h > 0.0 && h.is_finite())
        && (0..6).all(|face| {
            let mut q = [I::new(-1.0, 1.0); 3];
            q[face / 2] = I::point(if face % 2 == 0 { -1.0 } else { 1.0 });
            exterior_empty(t, ctx, &cylinders, Region::Sector { c, r, q, s: I::new(1.0, f64::INFINITY) }, 6)
        })
}
