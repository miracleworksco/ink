use super::joins::rejoins;
use super::{BOTH, FIRST, Op, SECOND, Tape};
use crate::eval::{Ctx, eval_interval_clamped, eval_point};
use crate::graph::{BinOp, Graph, UnOp, Var};
use crate::interval::I;

/// A subset's search field, max(a + p0, -(b + p1)), rises with p0 and
/// falls with p1 wherever it is; read through |.|, or both ways, p2 moves
/// it either way; p3 it does not read.
#[test]
fn a_param_s_way_is_the_one_every_read_moves_the_output() {
    let mut g = Graph::new();
    let (x, y) = (g.var(Var::X), g.var(Var::Y));
    let [p0, p1, p2] = [0, 1, 2].map(|k| g.var(Var::Param(k)));
    let a = g.add(x, p0);
    let b = g.add(y, p1);
    let nb = g.un(UnOp::Neg, b);
    let f = g.max(a, nb);
    let c = g.add(x, p2);
    let both = g.un(UnOp::Abs, c);
    let out = g.min(f, both);
    let t = Tape::new(&g, out);
    assert_eq!([0, 1, 2, 3].map(|k| t.param_way(k)), [Some(1), Some(-1), Some(0), None]);
}

/// A min of a and a + c reads c's slot to order them (Tape::apart), so a
/// jump in c joins no earlier than that min, however `monotone` is
/// decided: here with every value but the jumps taken as monotone (a
/// wrong rule), c = atan2(x, z - 2.64) and a = |y| would join at a + c,
/// and the min after it would order a and a + c by the slot of c as the
/// second side left it.
#[test]
fn a_jump_joins_after_a_min_that_reads_it_to_order_two_values() {
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let a = g.un(UnOp::Abs, y);
    let k = g.constant(2.64);
    let zk = g.sub(z, k);
    let c = g.bin(BinOp::Atan2, x, zk);
    let v = g.add(a, c);
    let m = g.min(a, v);
    let fz = g.un(UnOp::Floor, z);
    let inner = g.max(v, m);
    let out = g.max(inner, fz);
    let t = Tape::new(&g, out);
    let ordered = (0..t.ssa.len()).find(|&i| Tape::apart_of(&t.ssa, i).is_some()).expect("a min of a and a + c");
    let wrong: Vec<bool> = (0..t.ssa.len()).map(|i| !t.jumps.contains(&(i as u32))).collect();
    let r = rejoins(&t.ssa, &t.jumps, &wrong);
    let at = t.jumps.iter().find(|&&j| matches!(t.ssa[j as usize].op, Op::Bin(BinOp::Atan2))).copied().unwrap();
    assert!(r[at as usize] as usize >= ordered, "atan2 at {at} joins at {}, before the min at {ordered} that reads it", r[at as usize]);
}

/// A simplified tape keeps its order but drops edges, so a jump's window
/// can hold another jump that feeds nothing the window joins: here
/// floor y joins at max(sin floor y, min(sin floor z + 10, floor y)),
/// the min decided for floor y over the box (and, after floor x is
/// taken as the case 0, for floor y plus 20 floor x), with floor z
/// inside the window, read past it by cos. Its window ends where its
/// value's last reader does; a window reaching floor z's join there
/// read the max as the second side left it, -9 to -9 where f is -10.
#[test]
fn a_simplified_tape_joins_each_jump_where_its_value_is_last_read() {
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let j = g.un(UnOp::Floor, y);
    let sj = g.un(UnOp::Sin, j);
    let k = g.un(UnOp::Floor, z);
    let sk = g.un(UnOp::Sin, k);
    let c10 = g.constant(10.0);
    let q = g.add(sk, c10);
    let fx = g.un(UnOp::Floor, x);
    let c20 = g.constant(20.0);
    let wx = g.mul(fx, c20);
    let shifted = g.add(j, wx);
    let ck = g.un(UnOp::Cos, k);
    let c11 = g.constant(11.0);
    let h = g.sub(ck, c11);
    let b = [I::new(-1.0, 1.0), I::new(0.5, 1.5), I::new(0.2, 0.4)];
    let ctx = Ctx::new(0.0, vec![]);
    // f(0.5, 0.7, 0.3) = sin 0 max (sin 0 + 10) min 0, plus cos 0 - 11: -10.
    let p = [0.5, 0.7, 0.3];
    let bound_of = |t: &Tape| {
        let mut ch = vec![BOTH; t.len()];
        crate::eval::eval_interval(t, &ctx, b, Some(&mut ch));
        let s = t.simplify(&ch);
        (crate::eval::eval_interval(&s, &ctx, b, None), eval_point(&s, &ctx, p))
    };
    let plain = {
        let m = g.min(q, j);
        let r = g.max(sj, m);
        let out = g.add(r, h);
        bound_of(&Tape::new(&g, out))
    };
    let cased = {
        let m = g.min(q, shifted);
        let r = g.max(sj, m);
        let out = g.add(r, h);
        let t = Tape::new(&g, out);
        let at = (0..t.len()).find(|&i| matches!(t.ssa[i].op, Op::Un(UnOp::Floor)) && matches!(t.ssa[t.ssa[i].a as usize].op, Op::Var(Var::X))).unwrap();
        bound_of(&t.with_consts(&[(at, 0.0)]))
    };
    for (r, v) in [plain, cased] {
        assert!(v == -10.0 && r.lo <= v && v <= r.hi, "{r:?} misses {v}");
    }
}

/// A term shared by many is one conjunct, whatever the paths to it: a
/// sphere's field under twenty emaxes of itself with itself is one term
/// (walked as a tree, 2^20 of them; thirty levels, 16 GiB).
#[test]
fn conjuncts_are_as_many_as_the_instructions_at_most() {
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let (xx, yy, zz) = (g.un(UnOp::Square, x), g.un(UnOp::Square, y), g.un(UnOp::Square, z));
    let xy = g.add(xx, yy);
    let r2 = g.add(xy, zz);
    let r = g.un(UnOp::Sqrt, r2);
    let one = g.constant(1.0);
    let mut f = g.sub(r, one);
    for _ in 0..20 {
        f = g.bin(BinOp::Emax, f, f);
    }
    let t = Tape::new(&g, f);
    let terms = t.conjuncts();
    assert!(terms.len() <= t.len() && terms.len() == 1, "{} terms of a tape {} long", terms.len(), t.len());
}

/// A min of a and a - c, c > 0 (a blend below what it rounds, named),
/// is the blend wherever c's bounds are above 0, though a's and the
/// blend's bounds overlap (each reads a, here x^2 over a box about
/// 0).
#[test]
fn a_min_of_a_and_a_less_c_is_ordered_by_c() {
    let mut g = Graph::new();
    let (x, y) = (g.var(Var::X), g.var(Var::Y));
    let a = g.un(UnOp::Square, x);
    let yy = g.un(UnOp::Square, y);
    let one = g.constant(1.0);
    let c = g.add(yy, one);
    let less = g.sub(a, c);
    let blend = g.tag("blend", less);
    for f in [g.min(a, blend), g.min(blend, a)] {
        let t = Tape::new(&g, f);
        let mut ch = vec![BOTH; t.len()];
        let b = [I::new(-1.0, 1.0), I::new(-0.5, 0.5), I::point(0.0)];
        let r = crate::eval::eval_interval(&t, &Ctx::new(0.0, vec![]), b, Some(&mut ch));
        // The blend's own bound: x^2 - (y^2 + 1) over the box.
        assert!(r.lo >= -1.25 - 1e-12 && r.hi <= 0.0, "{r:?}");
        let i = t.len() - 1;
        let picked = match ch[i] {
            FIRST => t.ssa[i].a,
            SECOND => t.ssa[i].b,
            _ => panic!("undecided: {r:?}"),
        };
        assert!(matches!(t.ssa[picked as usize].op, Op::Tag(_)), "{:?}", t.ssa[i]);
    }
}

/// A tape simplified by its sign choices over a box has its sign (below,
/// at or above 0) at every point of the box: through mins, maxes, an
/// emax, a negation and a product with a positive factor, and where a
/// branch reads 0 throughout a slab (max(0, floor(z + 0.5)), which no
/// choice may drop as below 0).
#[test]
fn sign_choices_keep_the_sign_everywhere() {
    let mut g = Graph::new();
    let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
    let (quarter, half, zero, two) = (g.constant(0.25), g.constant(0.5), g.constant(0.0), g.constant(2.0));
    let zh = g.add(z, half);
    let step = g.un(UnOp::Floor, zh);
    let plateau = g.max(zero, step);
    let xq = g.sub(x, quarter);
    let left = g.min(xq, plateau);
    let yy = g.mul(y, y);
    let ring = g.sub(yy, quarter);
    let slab = g.sub(z, half);
    let cap = g.bin(BinOp::Emax, ring, slab);
    let out = g.neg(cap);
    let scaled = g.mul(two, out);
    let f = g.max(left, scaled);
    let t = Tape::new(&g, f);
    let ctx = Ctx::new(0.0, vec![]);
    let class = |v: f64| {
        if v < 0.0 {
            -1
        } else if v > 0.0 {
            1
        } else {
            0
        }
    };
    let mut made = 0;
    for i in 0..64 {
        let at = |k: usize| -1.0 + 0.5 * ((i >> (2 * k)) & 3) as f64;
        let b: [I; 3] = std::array::from_fn(|k| I::new(at(k), at(k) + 0.5));
        let mut vals = Vec::new();
        eval_interval_clamped(&t, &ctx, b, &[], &mut vals).unwrap();
        let mut ch = vec![BOTH; t.len()];
        if !t.sign_choices(&vals, &mut ch) {
            continue;
        }
        made += 1;
        let s = t.simplify(&ch);
        assert!(s.len() < t.len());
        for q in 0..216 {
            let p: crate::V3 = std::array::from_fn(|k| b[k].lo + (b[k].hi - b[k].lo) * ((q / [1, 6, 36][k]) % 6) as f64 / 5.0);
            assert_eq!(class(eval_point(&t, &ctx, p)), class(eval_point(&s, &ctx, p)), "at {p:?} in {b:?}");
        }
    }
    assert!(made > 8, "choices made in {made} boxes");
}
