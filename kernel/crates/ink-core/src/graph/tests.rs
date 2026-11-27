use super::{BinOp, Fingerprint, Graph, Node, NodeId, UnOp, Var, bin_f64, dun, un_f64};

#[test]
fn a_number_known_within_bounds_is_read_over_them() {
    // Ranges read all of its bounds, through arithmetic on constants
    // too (never folded to its value); a point reads its value; it is
    // constant in space.
    use crate::eval::{Ctx, eval_interval, eval_point};
    use crate::interval::I;
    let mut g = Graph::new();
    let x = g.var(Var::X);
    let k = g.within(1.0, 2.0, 1.5, 1);
    let two = g.constant(2.0);
    let k2 = g.mul(k, two);
    let f = g.sub(x, k2);
    let t = crate::tape::Tape::new(&g, f);
    let ctx = Ctx::new(0.0, vec![]);
    let p = I::point(0.0);
    assert_eq!(eval_interval(&t, &ctx, [p, p, p], None), I::new(-4.0, -2.0));
    assert_eq!(eval_point(&t, &ctx, [0.0; 3]), -3.0);
    let d = g.diff(f, Var::X);
    assert_eq!(g.as_const(d), Some(1.0));
    assert_eq!(g.within(1.0, 1.0, 1.0, 2), g.constant(1.0));
}

#[test]
fn numbers_known_within_the_same_bounds_are_two_numbers() {
    // Two exact numbers can round to the same bounds: two identities,
    // their difference bounded, never folded to 0; one identity is one
    // number, less itself 0.
    use crate::eval::{Ctx, eval_interval};
    use crate::interval::I;
    let mut g = Graph::new();
    let (a, b) = (g.within(1.0, 2.0, 1.5, 1), g.within(1.0, 2.0, 1.5, 2));
    assert_eq!(g.within(1.0, 2.0, 1.5, 1), a);
    assert_ne!(a, b);
    let d = g.sub(a, b);
    assert_eq!(g.as_const(d), None);
    let t = crate::tape::Tape::new(&g, d);
    let p = I::point(0.0);
    assert_eq!(eval_interval(&t, &Ctx::new(0.0, vec![]), [p, p, p], None), I::new(-1.0, 1.0));
    let z = g.sub(a, a);
    assert_eq!(g.as_const(z), Some(0.0));
}

#[test]
fn sealed_names_change_the_fingerprint() {
    // A sealed name (provenance stops there: a fillet's face reports the
    // fillet) means different face names, so results keyed by content
    // must not be shared with the same expression under a plain name.
    let print = |sealed: bool| {
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let r = if sealed { g.tag_opaque("f", x) } else { g.tag("f", x) };
        let mut fp = Fingerprint::default();
        g.fingerprint(&[r], &[], true, &mut fp);
        fp.value()
    };
    assert_ne!(print(true), print(false));
}

/// Two operands of a sum that differ only in which param they read
/// are ordered by where the question first reads each, not by where
/// the graph made them: a model that declares or imports its params in
/// another order asks the same question.
#[test]
fn a_sum_of_two_params_has_one_key_whatever_their_order() {
    let key = |q_first: bool| {
        let mut g = Graph::new();
        let (p, q) = if q_first {
            let q = g.var(Var::Param(1));
            (g.var(Var::Param(0)), q)
        } else {
            let p = g.var(Var::Param(0));
            (p, g.var(Var::Param(1)))
        };
        let (sp, sq) = if q_first {
            let sq = g.un(UnOp::Sin, q);
            (g.un(UnOp::Sin, p), sq)
        } else {
            let sp = g.un(UnOp::Sin, p);
            (sp, g.un(UnOp::Sin, q))
        };
        let sum = g.add(sp, sq);
        let root = g.sub(p, sum);
        let tape = crate::tape::Tape::new(&g, root);
        let ctx = crate::eval::Ctx::new(0.0, vec![0.25, 0.5]);
        let names = crate::tape::Names::of(&[&tape], &ctx);
        let mut fp = Fingerprint::default();
        tape.fingerprint(&mut fp, &names);
        fp.value()
    };
    assert_eq!(key(false), key(true));
}

/// Neither read before, two params' terms are ordered by the params'
/// names, not by which the model declared or the graph made first: the
/// question, its values included, has one key.
#[test]
fn a_sum_of_two_unread_params_has_one_key_whatever_their_order() {
    let key = |q_first: bool| {
        let mut g = Graph::new();
        // q_first: q is param 0, declared and built first.
        let (kp, kq) = if q_first { (1, 0) } else { (0, 1) };
        g.name_param(kp, "p");
        g.name_param(kq, "q");
        let (sp, sq) = if q_first {
            let q = g.var(Var::Param(kq));
            let sq = g.un(UnOp::Sin, q);
            let p = g.var(Var::Param(kp));
            (g.un(UnOp::Sin, p), sq)
        } else {
            let p = g.var(Var::Param(kp));
            let sp = g.un(UnOp::Sin, p);
            let q = g.var(Var::Param(kq));
            (sp, g.un(UnOp::Sin, q))
        };
        let root = g.add(sp, sq);
        let tape = crate::tape::Tape::new(&g, root);
        let mut values = vec![0.0; 2];
        (values[kp as usize], values[kq as usize]) = (0.25, 0.5);
        let ctx = crate::eval::Ctx::new(0.0, values);
        let names = crate::tape::Names::of(&[&tape], &ctx);
        let mut fp = Fingerprint::default();
        tape.fingerprint(&mut fp, &names);
        assert!(ctx.fingerprint(&mut fp, &names));
        fp.value()
    };
    assert_eq!(key(false), key(true));
}

/// Terms alike but for their params are ordered by their own terms in
/// order, inner ties first: x - (sin (p + zz) + sin (q + r)) has one
/// key whichever of p and zz the graph made first (an unrelated
/// `(def unused zz)` makes zz first).
#[test]
fn nested_ties_are_ordered_from_the_inside_out() {
    let key = |zz_first: bool| {
        let mut g = Graph::new();
        for (k, n) in ["p", "zz", "q", "r"].iter().enumerate() {
            g.name_param(k as u32, n);
        }
        let (p, zz) = if zz_first {
            let zz = g.var(Var::Param(1));
            (g.var(Var::Param(0)), zz)
        } else {
            let p = g.var(Var::Param(0));
            (p, g.var(Var::Param(1)))
        };
        let (q, r) = (g.var(Var::Param(2)), g.var(Var::Param(3)));
        let pz = g.add(p, zz);
        let qr = g.add(q, r);
        let (a, b) = (g.un(UnOp::Sin, pz), g.un(UnOp::Sin, qr));
        let sum = g.add(a, b);
        let x = g.var(Var::X);
        let root = g.sub(x, sum);
        let tape = crate::tape::Tape::new(&g, root);
        let ctx = crate::eval::Ctx::new(0.0, vec![0.1, 0.2, 0.3, 0.4]);
        let names = crate::tape::Names::of(&[&tape], &ctx);
        let mut fp = Fingerprint::default();
        tape.fingerprint(&mut fp, &names);
        assert!(ctx.fingerprint(&mut fp, &names));
        fp.value()
    };
    assert_eq!(key(false), key(true));
}

/// Library functions give musl's bits, which no CPU dispatch changes.
/// glibc's differ in the last bit at these inputs (by CPU: sin at p,
/// where 1e16 (sin p + 0.3775604968083285) is 0 or -0.56, so a model
/// would be proven empty on one machine and disproven on another, and
/// the result cache would hand one's answer to the other).
#[test]
fn library_functions_are_the_same_on_every_machine() {
    let p = -21.60398819379739;
    assert_eq!(un_f64(UnOp::Sin, p), -0.3775604968083285);
    assert_eq!(un_f64(UnOp::Sin, -1499.931286), 0.9839855162426463);
    assert_eq!(un_f64(UnOp::Cos, -1499.877923), -0.23047820418116294);
    assert_eq!(un_f64(UnOp::Tan, -1499.88304), -4.320450143326324);
    assert_eq!(bin_f64(BinOp::Atan2, -1499.548242, 1.0), -1.5701294593852633);
    assert_eq!(un_f64(UnOp::Exp, -2.142849832857143), 0.11732002370048947);
    assert_eq!(un_f64(UnOp::Ln, 1499.959064), 7.313193096051238);
    // The derivatives' values, the constants a graph folds and the
    // range rules' ends are the same functions.
    assert_eq!(dun(UnOp::Sin, -1499.931286).0, 0.9839855162426463);
    let mut g = Graph::default();
    let c = g.constant(-1499.931286);
    let s = g.un(UnOp::Sin, c);
    assert_eq!(g.as_const(s), Some(0.9839855162426463));
    let r = crate::interval::un(UnOp::Exp, crate::interval::I::point(-2.142849832857143));
    assert_eq!(r, crate::interval::widen(crate::interval::I::point(0.11732002370048947), 4));
}

#[test]
fn a_rewrite_keeps_a_promise_only_where_it_still_holds() {
    // (lipschitz 1 (- x p)) holds for every p; with p warped to 1000 x
    // it is -999 x: kept, the promise would prove a clearance of 2
    // across a 0.9 gap.
    use crate::eval::{Ctx, eval_interval_grad};
    use crate::interval::{I, lipschitz};
    let slope = |g: &Graph, f: NodeId| {
        let t = crate::tape::Tape::new(g, f);
        lipschitz(&eval_interval_grad(&t, &Ctx::new(0.0, vec![0.5]), [I::new(-1.0, 1.0); 3]))
    };
    let mut g = Graph::new();
    let (x, p, t) = (g.var(Var::X), g.var(Var::Param(0)), g.var(Var::T));
    let (one, k) = (g.constant(1.0), g.constant(1000.0));
    let d = g.sub(x, p);
    let f = g.lipschitz(d, one);
    let kx = g.mul(k, x);
    let warped = g.subst(f, &[(Var::Param(0), kx)]);
    let replaced = g.replace(f, &[(p, kx)].into());
    assert!(slope(&g, warped) >= 999.0 && slope(&g, replaced) >= 999.0);
    // x (cos^2 p + sin^2 p) is x; with cos^2 p replaced by 1000, a
    // non-spatial node, it is not 1-Lipschitz.
    let (c, s) = (g.un(UnOp::Cos, p), g.un(UnOp::Sin, p));
    let (c2, s2) = (g.un(UnOp::Square, c), g.un(UnOp::Square, s));
    let sum = g.add(c2, s2);
    let e = g.mul(x, sum);
    let f2 = g.lipschitz(e, one);
    let swapped = g.replace(f2, &[(c2, k)].into());
    assert!(slope(&g, swapped) >= 1000.0);
    // Kept: x warped (the frame stretches with it), p to an expression
    // of t.
    let two = g.constant(2.0);
    let (x2, t2) = (g.mul(two, x), g.mul(two, t));
    for w in [g.subst(f, &[(Var::X, x2)]), g.subst(f, &[(Var::Param(0), t2)]), g.replace(f, &[(x, x2)].into())] {
        assert!(matches!(g.nodes[w as usize], Node::Lipschitz(..)));
    }
}
