//! Batched point evaluation and point derivatives.

use super::context::Ctx;
use super::corner::corner_reach_point;
use crate::V3;
use crate::graph::{BinOp, UnOp, Var, bin_f64, dbin, defined, div_f64, dun, un_f64};
use crate::tape::{Op, Tape};

/// Lane width for batched evaluation (8 AVX-512 vectors of f64).
pub const L: usize = 64;

pub fn eval_points(tape: &Tape, ctx: &Ctx, pts: &[V3], out: &mut Vec<f64>) {
    // A lone point in one lane, a few (a cell's samples) in narrow
    // lanes, as eval_grads.
    match pts.len() {
        1 => points_in::<1>(tape, ctx, pts, out),
        ..=8 => points_in::<8>(tape, ctx, pts, out),
        _ => points_in::<L>(tape, ctx, pts, out),
    }
}

/// `eval_points` in chunks of N points.
fn points_in<const N: usize>(tape: &Tape, ctx: &Ctx, pts: &[V3], out: &mut Vec<f64>) {
    out.clear();
    out.reserve(pts.len());
    let mut s = vec![[0.0f64; N]; tape.nslots];
    let (mut x, mut y, mut z) = ([0.0; N], [0.0; N], [0.0; N]);
    for chunk in pts.chunks(N) {
        for k in 0..N {
            let p = chunk.get(k).copied().unwrap_or([0.0; 3]);
            x[k] = p[0];
            y[k] = p[1];
            z[k] = p[2];
        }
        run_lanes(tape, ctx, &x, &y, &z, &mut s);
        out.extend_from_slice(&s[tape.out_slot()][..chunk.len()]);
    }
}

/// `eval_points` split across all cores (for large point sets).
pub fn eval_points_par(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<f64> {
    use rayon::prelude::*;
    pts.par_chunks(4096)
        .flat_map_iter(|c| {
            let mut out = Vec::with_capacity(c.len());
            eval_points(tape, ctx, c, &mut out);
            out
        })
        .collect()
}

pub fn eval_point(tape: &Tape, ctx: &Ctx, p: V3) -> f64 {
    let mut s = vec![[0.0f64; 1]; tape.nslots];
    run_lanes(tape, ctx, &[p[0]], &[p[1]], &[p[2]], &mut s);
    s[tape.out_slot()][0]
}

/// The tape over N points at once, one lane each (L for batches, 1 for a
/// lone point: the same code, so the same values).
fn run_lanes<const N: usize>(tape: &Tape, ctx: &Ctx, x: &[f64; N], y: &[f64; N], z: &[f64; N], s: &mut [[f64; N]]) {
    for (i, ins) in tape.code.iter().enumerate() {
        let o = ins.out as usize;
        match ins.op {
            Op::Const(v) => s[o] = [v; N],
            Op::Var(Var::X) => s[o] = *x,
            Op::Var(Var::Y) => s[o] = *y,
            Op::Var(Var::Z) => s[o] = *z,
            Op::Var(Var::T) => s[o] = [ctx.t; N],
            Op::Var(Var::Param(p)) => s[o] = [ctx.params[p as usize]; N],
            Op::Un(op) => {
                let a = ins.a as usize;
                match op {
                    UnOp::Neg => lanes(s, o, [a; 3], |x| -x[0]),
                    UnOp::Abs => lanes(s, o, [a; 3], |x| x[0].abs()),
                    UnOp::Sqrt => lanes(s, o, [a; 3], |x| x[0].max(0.0).sqrt()),
                    UnOp::Square => lanes(s, o, [a; 3], |x| x[0] * x[0]),
                    UnOp::Floor => lanes(s, o, [a; 3], |x| x[0].floor()),
                    _ => lanes(s, o, [a; 3], |x| un_f64(op, x[0])),
                }
            }
            Op::Bin(op) => {
                let ab = [ins.a as usize, ins.b as usize, ins.b as usize];
                match op {
                    BinOp::Add => lanes(s, o, ab, |x| defined(x[0] + x[1])),
                    BinOp::Sub => lanes(s, o, ab, |x| defined(x[0] - x[1])),
                    BinOp::Mul => lanes(s, o, ab, |x| defined(x[0] * x[1])),
                    BinOp::Div => lanes(s, o, ab, |x| div_f64(x[0], x[1])),
                    BinOp::Min => lanes(s, o, ab, |x| x[0].min(x[1])),
                    BinOp::Max => lanes(s, o, ab, |x| x[0].max(x[1])),
                    _ => lanes(s, o, ab, |x| bin_f64(op, x[0], x[1])),
                }
            }
            Op::Select => lanes(s, o, [ins.a as usize, ins.b as usize, ins.c as usize], |x| if x[0] < 0.0 { x[1] } else { x[2] }),
            Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real => s[o] = s[ins.a as usize],
            Op::Within(_) => s[o] = s[ins.c as usize],
            Op::Frame => s[o] = [0.0; N],
            Op::Grid(gi, m) => {
                let (x, y, z) = (s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
                let g = &tape.grids[gi as usize];
                let r = &mut s[o];
                (0..N).for_each(|k| r[k] = g.component([x[k], y[k], z[k]], m));
            }
            Op::CornerReach(fa, fb) => {
                let rv = s[ins.a as usize];
                let gv = s[ins.c as usize];
                let f = &tape.code[tape.ssa[i].b as usize];
                let (cx, cy, cz) = (s[f.a as usize], s[f.b as usize], s[f.c as usize]);
                let (ta, tb) = (&tape.forms[fa as usize], &tape.forms[fb as usize]);
                let mut out = [0.0; N];
                for k in 0..N {
                    if gv[k] < 0.0 {
                        out[k] = corner_reach_point(ta, tb, ctx, rv[k], [cx[k], cy[k], cz[k]]);
                    }
                }
                s[o] = out;
            }
        }
    }
}

/// Slot `o` of each lane as `f` of slots `args` there, read where they are:
/// copying a slot of L lanes per operand cost more than the arithmetic. An
/// output that reuses an operand's slot is computed whole first and then
/// stored: each lane reads only its own lane, so the values are the same
/// as writing lane by lane, and the loop is the same straight one.
#[inline(always)]
fn lanes<const N: usize>(s: &mut [[f64; N]], o: usize, args: [usize; 3], f: impl Fn([f64; 3]) -> f64) {
    if args.contains(&o) {
        let x = args.map(|a| &s[a]);
        let mut r = [0.0; N];
        for k in 0..N {
            r[k] = f([x[0][k], x[1][k], x[2][k]]);
        }
        s[o] = r;
        return;
    }
    let (below, rest) = s.split_at_mut(o);
    let (r, above) = rest.split_first_mut().unwrap();
    let at = |a: usize| if a < o { &below[a] } else { &above[a - o - 1] };
    let x = args.map(at);
    for k in 0..N {
        r[k] = f([x[0][k], x[1][k], x[2][k]]);
    }
}

/// Value and spatial gradient `[f, df/dx, df/dy, df/dz]` at each point.
pub fn eval_grads(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<[f64; 4]> {
    // A few points (a cell's samples) in narrow lanes: the slots are
    // allocated, and each op run, per lane.
    match pts.len() {
        1 => grads_in::<1>(tape, ctx, pts),
        0 | 2..=8 => grads_in::<8>(tape, ctx, pts),
        9..=16 => grads_in::<16>(tape, ctx, pts),
        17..=32 => grads_in::<32>(tape, ctx, pts),
        _ => grads_in::<L>(tape, ctx, pts),
    }
}

/// Gradients at intermediate instructions of the same tape. Record the
/// results before their slots are reused. Indices must be distinct, in
/// increasing order and in range.
pub fn eval_grads_at(tape: &Tape, ctx: &Ctx, pts: &[V3], instructions: &[usize]) -> Vec<Vec<[f64; 4]>> {
    assert!(instructions.windows(2).all(|w| w[0] < w[1]) && instructions.last().is_none_or(|&i| i < tape.len()));
    if instructions.is_empty() {
        return Vec::new();
    }
    if instructions == [tape.len() - 1] {
        return vec![eval_grads(tape, ctx, pts)];
    }
    let mut picked: Vec<_> = instructions.iter().map(|_| Vec::with_capacity(pts.len())).collect();
    match pts.len() {
        1 => grad_batches::<1, true>(tape, ctx, pts, instructions, &mut picked),
        ..=8 => grad_batches::<8, true>(tape, ctx, pts, instructions, &mut picked),
        _ => grad_batches::<L, true>(tape, ctx, pts, instructions, &mut picked),
    };
    picked
}

/// `eval_grads` in chunks of N points.
fn grads_in<const N: usize>(tape: &Tape, ctx: &Ctx, pts: &[V3]) -> Vec<[f64; 4]> {
    grad_batches::<N, false>(tape, ctx, pts, &[], &mut [])
}

#[cfg(test)]
thread_local! { pub(super) static GRAD_INSTRUCTIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) }; }

/// Capture the requested gradients, or return just the final gradient.
fn grad_batches<const N: usize, const PICK: bool>(tape: &Tape, ctx: &Ctx, pts: &[V3], instructions: &[usize], picked: &mut [Vec<[f64; 4]>]) -> Vec<[f64; 4]> {
    let mut out = if PICK { Vec::new() } else { Vec::with_capacity(pts.len()) };
    let mut s = vec![[[0.0f64; N]; 4]; tape.nslots];
    let end = if PICK { instructions[instructions.len() - 1] + 1 } else { tape.len() };
    for chunk in pts.chunks(N) {
        #[cfg(test)]
        GRAD_INSTRUCTIONS.with(|work| work.set(work.get() + end));
        let mut next = 0;
        for (i, ins) in tape.code[..end].iter().enumerate() {
            let o = ins.out as usize;
            match ins.op {
                Op::Const(v) => s[o] = [[v; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Var(v @ (Var::X | Var::Y | Var::Z)) => {
                    let axis = match v {
                        Var::X => 0,
                        Var::Y => 1,
                        _ => 2,
                    };
                    let r = &mut s[o];
                    *r = [[0.0; N]; 4];
                    for (k, p) in chunk.iter().enumerate() {
                        r[0][k] = p[axis];
                    }
                    r[axis + 1] = [1.0; N];
                }
                Op::Var(Var::T) => s[o] = [[ctx.t; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Var(Var::Param(p)) => s[o] = [[ctx.params[p as usize]; N], [0.0; N], [0.0; N], [0.0; N]],
                Op::Un(op) => grad_lanes(&mut s, o, [ins.a as usize; 3], |x| {
                    let (v, fa) = dun(op, x[0][0]);
                    [v, mul0(fa, x[0][1]), mul0(fa, x[0][2]), mul0(fa, x[0][3])]
                }),
                Op::Bin(op) => {
                    let args = [ins.a as usize, ins.b as usize, ins.b as usize];
                    // Dispatch before the lane loop for the common ops,
                    // so its arithmetic can be vectorized.
                    match op {
                        BinOp::Add => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Add, a, b)),
                        BinOp::Sub => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Sub, a, b)),
                        BinOp::Mul => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Mul, a, b)),
                        BinOp::Div => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Div, a, b)),
                        BinOp::Min => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Min, a, b)),
                        BinOp::Max => bin_grad(&mut s, o, args, |a, b| dbin(BinOp::Max, a, b)),
                        _ => bin_grad(&mut s, o, args, |a, b| dbin(op, a, b)),
                    }
                }
                Op::Select => grad_lanes(&mut s, o, [ins.a as usize, ins.b as usize, ins.c as usize], |x| if x[0][0] < 0.0 { x[1] } else { x[2] }),
                Op::Read => {
                    // An op of its coordinates: its partials (its form's
                    // gradient there) chained with theirs.
                    let f = &tape.code[tape.ssa[i].b as usize];
                    let (cx, cy, cz) = (s[f.a as usize], s[f.b as usize], s[f.c as usize]);
                    let at: Vec<V3> = (0..chunk.len()).map(|k| [cx[0][k], cy[0][k], cz[0][k]]).collect();
                    let gr = eval_grads(&tape.forms[ins.c as usize], ctx, &at);
                    let v = s[ins.a as usize][0];
                    let r = &mut s[o];
                    r[0] = v;
                    for (k, gr) in gr.iter().enumerate() {
                        for j in 1..4 {
                            r[j][k] = chain([gr[1], gr[2], gr[3]], [cx[j][k], cy[j][k], cz[j][k]]);
                        }
                    }
                }
                Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Real => s[o] = s[ins.a as usize],
                Op::Within(_) => s[o] = [s[ins.c as usize][0], [0.0; N], [0.0; N], [0.0; N]],
                Op::Frame => s[o] = [[0.0; N]; 4],
                Op::Grid(gi, m) => {
                    // Chain rule through the coordinate expressions.
                    let (cx, cy, cz) = (s[ins.a as usize], s[ins.b as usize], s[ins.c as usize]);
                    let g = &tape.grids[gi as usize];
                    let r = &mut s[o];
                    for k in 0..N {
                        let (v, gr) = g.component_grad([cx[0][k], cy[0][k], cz[0][k]], m);
                        r[0][k] = v;
                        for j in 1..4 {
                            r[j][k] = chain(gr, [cx[j][k], cy[j][k], cz[j][k]]);
                        }
                    }
                }
                Op::CornerReach(fa, fb) => {
                    let rv = s[ins.a as usize][0];
                    let gv = s[ins.c as usize][0];
                    let f = &tape.code[tape.ssa[i].b as usize];
                    let (cx, cy, cz) = (s[f.a as usize][0], s[f.b as usize][0], s[f.c as usize][0]);
                    let (ta, tb) = (&tape.forms[fa as usize], &tape.forms[fb as usize]);
                    let mut out = [[0.0; N]; 4];
                    for k in 0..chunk.len() {
                        if gv[k] < 0.0 {
                            out[0][k] = corner_reach_point(ta, tb, ctx, rv[k], [cx[k], cy[k], cz[k]]);
                        }
                    }
                    s[o] = out;
                }
            }
            if PICK && instructions.get(next) == Some(&i) {
                let r = &s[o];
                picked[next].extend((0..chunk.len()).map(|k| [r[0][k], r[1][k], r[2][k], r[3][k]]));
                next += 1;
            }
        }
        if !PICK {
            let r = &s[tape.out_slot()];
            for k in 0..chunk.len() {
                out.push([r[0][k], r[1][k], r[2][k], r[3][k]]);
            }
        }
    }
    out
}

#[inline(always)]
fn bin_grad<const N: usize>(s: &mut [[[f64; N]; 4]], o: usize, args: [usize; 3], partials: impl Fn(f64, f64) -> (f64, f64, f64)) {
    grad_lanes(s, o, args, |x| {
        let (v, fa, fb) = partials(x[0][0], x[1][0]);
        let d = |j: usize| defined(mul0(fa, x[0][j]) + mul0(fb, x[1][j]));
        [v, d(1), d(2), d(3)]
    });
}

/// `lanes` for values and their gradients (`grads_in`'s slots): each
/// lane's value and partials as `f` of its operands'. An output that reuses
/// an operand's slot is computed whole, then stored, as `lanes` does:
/// written lane by lane through the shared slot, each lane's operands were
/// copied out one at a time and read back wider than written, which stalled
/// (half of `grads_in`'s time on the Falcon swerve check).
#[inline(always)]
fn grad_lanes<const N: usize>(s: &mut [[[f64; N]; 4]], o: usize, args: [usize; 3], f: impl Fn([[f64; 4]; 3]) -> [f64; 4]) {
    let lane = |x: &[[f64; N]; 4], k: usize| [x[0][k], x[1][k], x[2][k], x[3][k]];
    if args.contains(&o) {
        let x = args.map(|a| &s[a]);
        let mut out = [[0.0; N]; 4];
        for k in 0..N {
            let r = f([lane(x[0], k), lane(x[1], k), lane(x[2], k)]);
            (0..4).for_each(|j| out[j][k] = r[j]);
        }
        s[o] = out;
        return;
    }
    let (below, rest) = s.split_at_mut(o);
    let (out, above) = rest.split_first_mut().unwrap();
    let x = args.map(|a| if a < o { &below[a] } else { &above[a - o - 1] });
    for k in 0..N {
        let r = f([lane(x[0], k), lane(x[1], k), lane(x[2], k)]);
        (0..4).for_each(|j| out[j][k] = r[j]);
    }
}

/// The chain rule's product: `0 * inf = 0` (a zero partial means "no
/// dependence"), as the value rules have it.
#[inline]
fn mul0(f: f64, d: f64) -> f64 {
    if f == 0.0 || d == 0.0 { 0.0 } else { f * d }
}

/// The chain rule through three coordinates: partials `f` times their
/// derivatives `d`, summed in order under the value rules (each sum
/// defined), as graph::Graph::diff writes it.
#[inline]
fn chain(f: V3, d: V3) -> f64 {
    defined(defined(mul0(f[0], d[0]) + mul0(f[1], d[1])) + mul0(f[2], d[2]))
}

/// Value and every partial derivative at one point:
/// `[f, d/dx, d/dy, d/dz, d/dt, d/dparam0, ...]`.
pub fn eval_dual(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<f64> {
    let dirs: Vec<usize> = (0..4 + ctx.params.len()).collect();
    dual_along(tape, ctx, p, &dirs)
}

/// `eval_dual` with the partials along ctx's ranged params alone (those
/// its derivations own and read among them), every other one 0: what a
/// witness search pushing those params reads (ink-kernel's
/// scope::pushed). Each partial is computed from its own direction's
/// alone, so these are eval_dual's to the bit, at a cost in the ranged
/// params, not in all of a model's (thousands in an assembly).
pub fn eval_dual_ranged(tape: &Tape, ctx: &Ctx, p: V3) -> Vec<f64> {
    let mut ks = ctx.p_spans.along();
    ks.retain(|&k| (k as usize) < ctx.params.len());
    let dirs: Vec<usize> = ks.iter().map(|&k| 4 + k as usize).collect();
    let d = dual_along(tape, ctx, p, &dirs);
    let mut out = vec![0.0; 5 + ctx.params.len()];
    out[0] = d[0];
    for (j, &dir) in dirs.iter().enumerate() {
        out[1 + dir] = d[1 + j];
    }
    out
}

/// Value and the partials along `dirs` (eval_dual's numbering less one:
/// 0..3 x y z, 3 t, 4 + i param i; distinct), in that order after it.
fn dual_along(tape: &Tape, ctx: &Ctx, p: V3, dirs: &[usize]) -> Vec<f64> {
    let w = dirs.len() + 1;
    // Each direction's column (0: not carried).
    let mut col = vec![0usize; 4 + ctx.params.len()];
    for (j, &dir) in dirs.iter().enumerate() {
        col[dir] = 1 + j;
    }
    let (mut s, mut r) = (vec![0.0f64; tape.nslots * w], vec![0.0; w]);
    for (i, ins) in tape.code.iter().enumerate() {
        let o = ins.out as usize * w;
        r.fill(0.0);
        match ins.op {
            Op::Const(v) => r[0] = v,
            Op::Var(v) => {
                let (val, dir) = match v {
                    Var::X => (p[0], 0),
                    Var::Y => (p[1], 1),
                    Var::Z => (p[2], 2),
                    Var::T => (ctx.t, 3),
                    Var::Param(i) => (ctx.params[i as usize], 4 + i as usize),
                };
                r[0] = val;
                if col[dir] != 0 {
                    r[col[dir]] = 1.0;
                }
            }
            Op::Un(op) => {
                let a = &s[ins.a as usize * w..][..w];
                let (v, fa) = dun(op, a[0]);
                r[0] = v;
                for j in 1..w {
                    r[j] = mul0(fa, a[j]);
                }
            }
            Op::Bin(op) => {
                let a = &s[ins.a as usize * w..][..w];
                let b = &s[ins.b as usize * w..][..w];
                let (v, fa, fb) = dbin(op, a[0], b[0]);
                r[0] = v;
                for j in 1..w {
                    r[j] = defined(mul0(fa, a[j]) + mul0(fb, b[j]));
                }
            }
            Op::Select => {
                let c = s[ins.a as usize * w];
                let src = if c < 0.0 { ins.b } else { ins.c } as usize * w;
                r.copy_from_slice(&s[src..src + w]);
            }
            Op::Read => {
                // As eval_grads.
                let f = &tape.code[tape.ssa[i].b as usize];
                let (a, b, c) = (&s[f.a as usize * w..][..w], &s[f.b as usize * w..][..w], &s[f.c as usize * w..][..w]);
                let gr = eval_grads(&tape.forms[ins.c as usize], ctx, &[[a[0], b[0], c[0]]])[0];
                r[0] = s[ins.a as usize * w];
                for j in 1..w {
                    r[j] = chain([gr[1], gr[2], gr[3]], [a[j], b[j], c[j]]);
                }
            }
            Op::Tag(_) | Op::Lipschitz | Op::Deriv(_) | Op::Real => r.copy_from_slice(&s[ins.a as usize * w..][..w]),
            Op::Within(_) => r[0] = s[ins.c as usize * w],
            Op::Frame => {}
            Op::Grid(gi, m) => {
                let (a, b, c) = (&s[ins.a as usize * w..][..w], &s[ins.b as usize * w..][..w], &s[ins.c as usize * w..][..w]);
                let (v, gr) = tape.grids[gi as usize].component_grad([a[0], b[0], c[0]], m);
                r[0] = v;
                for j in 1..w {
                    r[j] = chain(gr, [a[j], b[j], c[j]]);
                }
            }
            Op::CornerReach(fa, fb) => {
                let gv = s[ins.c as usize * w];
                if gv < 0.0 {
                    let f = &tape.code[tape.ssa[i].b as usize];
                    let q = [s[f.a as usize * w], s[f.b as usize * w], s[f.c as usize * w]];
                    r[0] = corner_reach_point(&tape.forms[fa as usize], &tape.forms[fb as usize], ctx, s[ins.a as usize * w], q);
                }
            }
        }
        s[o..o + w].copy_from_slice(&r);
    }
    let o = tape.out_slot() * w;
    s[o..o + w].to_vec()
}
