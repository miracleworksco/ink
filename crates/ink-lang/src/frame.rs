//! Integrals over a placed body in its own frame: its shape there stays
//! put where its placement moves over a claim's scope (the motion enters
//! only the integrands), and is the same for every copy of the body
//! (claims::estimate_parts).

use crate::Placed;
use crate::view::ModelView;
use ink_core::graph::Var;
use ink_core::interval::{self, I};
use ink_core::tape::Tape;

/// A part integrated in its body's own frame (`in_frame`): the body's
/// shape, its integrands at R y + t, and 1 / |det R| over the scope.
pub(crate) struct InFrame {
    pub shape: u32,
    pub gs: Vec<Option<u32>>,
    pub inv_det: I,
}

/// A part integrated in its frame `f` (its integrands `gs` in the world),
/// placed by p = R y + t. A moving body's world field is wherever the body
/// may be; in its frame the shape stays put, and the motion enters only
/// the integrands. Copies of one body are one shape there, and one
/// integral where their integrands there are one field.
///
/// The world shape is S(A (x - t)), A = R^T (`Doc::items`), so where det R
/// is not 0, x = A^-1 y + t takes S's material onto the world shape's,
/// and the integral of g over that is the integral over S's material of
/// g(A^-1 y + t) / |det R|. Floats make a turn only nearly a rotation
/// (A^-1 = R^-T is not quite R), so this is used where the placement is a
/// rotation to within rounding (rigid::turn_bounds: its slack at most
/// 1e-9), or where no integrand reads x, y or z (the determinant alone
/// carries a stretch); elsewhere the world is used. Every entry of R and t
/// must be finite over the scope, and 1 / |det R| finite and above 0. An
/// entry of R or t that is one number over the scope, exactly
/// (rigid::constant_over: a turn's axis), is written as that number.
pub(crate) fn in_frame(v: &mut ModelView, f: &Placed, gs: &[Option<u32>]) -> Option<InFrame> {
    use ink_core::graph::BinOp;
    // The change of variables is one of numbers: an entry past the floats
    // over the scope (a move by e^1000, inf) places nothing in the world.
    let finite = |n: u32| {
        let x = ink_core::eval::eval_interval(&Tape::new(&v.doc.g, n), &v.ctx, [I::ENTIRE; 3], None);
        x.lo.is_finite() && x.hi.is_finite()
    };
    if !f.at.r.iter().flatten().chain(&f.at.t).all(|&n| finite(n)) {
        return None;
    }
    let turn = ink_kernel::rigid::turn_bounds(&v.doc.g, &f.at.r, &v.ctx)?;
    let det = match turn.det {
        d if d.lo > 0.0 => d,
        d if d.hi < 0.0 => I::new(-d.hi, -d.lo),
        _ => return None,
    };
    let slack = turn.slack()?;
    let inv_det = interval::bin(BinOp::Div, I::point(1.0), det).0;
    if !(inv_det.lo > 0.0 && inv_det.hi.is_finite()) {
        return None;
    }
    if slack > 1e-9 && gs.iter().flatten().any(|&g| v.doc.g.depends_on(g, |x| matches!(x, Var::X | Var::Y | Var::Z))) {
        return None;
    }
    let mut entry = |n: u32| match ink_kernel::rigid::constant_over(&v.doc.g, n, &v.ctx) {
        Some(c) => v.doc.g.constant(c),
        None => n,
    };
    let (r, t) = (f.at.r.map(|row| row.map(&mut entry)), f.at.t.map(&mut entry));
    let g = &mut v.doc.g;
    let p = [Var::X, Var::Y, Var::Z].map(|x| g.var(x));
    let map: Vec<(Var, u32)> = [Var::X, Var::Y, Var::Z]
        .into_iter()
        .enumerate()
        .map(|(i, x)| {
            let mut e = t[i];
            for (k, &pk) in p.iter().enumerate() {
                let m = g.mul(r[i][k], pk);
                e = g.add(e, m);
            }
            (x, e)
        })
        .collect();
    let gs = gs.iter().map(|h| h.map(|h| g.subst(h, &map))).collect();
    Some(InFrame { shape: f.shape, gs, inv_det })
}
