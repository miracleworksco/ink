use super::*;

/// Field `f` read at m q (m a signed permutation): x, y and z replaced by
/// signed x, y and z.
pub(super) fn relabel(g: &mut ink_core::graph::Graph, f: u32, m: &[[f64; 3]; 3]) -> u32 {
    let xyz = [Var::X, Var::Y, Var::Z];
    let map: Vec<(Var, u32)> = (0..3)
        .map(|i| {
            let (j, s) = (0..3).find_map(|j| (m[i][j] != 0.0).then_some((j, m[i][j]))).expect("a signed permutation");
            let q = g.var(xyz[j]);
            (xyz[i], if s > 0.0 { q } else { g.neg(q) })
        })
        .collect();
    g.subst(f, &map)
}

/// Field `f` with its fixed params read as the floats they are, as
/// `carried` reads centres.
fn fixed_params(v: &mut ModelView, f: u32) -> u32 {
    let ks: Vec<u32> = v
        .doc
        .g
        .reachable(&[f])
        .into_iter()
        .filter_map(|n| match v.doc.g.nodes[n as usize] {
            Node::Var(Var::Param(k)) if v.ctx.p_spans.get(k).is_none() => Some(k),
            _ => None,
        })
        .collect();
    let fixed: Vec<(Var, u32)> = ks.into_iter().map(|k| (Var::Param(k), v.doc.g.constant(v.ctx.params[k as usize]))).collect();
    if fixed.is_empty() { f } else { v.doc.g.subst(f, &fixed) }
}

/// Field `f` placed by `turn` (an isometry, a proper turn about the
/// vertical axis through `centre` of its own frame, then a move), its parts
/// every such turn leaves unchanged (symmetry::carried: a pin, a bore about
/// the axis) placed by the move alone, by `moved`: for such a part F,
/// F(R^T (p - T)) = F(p - (T + R l - l)) exactly, l the centre, and
/// `moved` is T + R l - l (T itself about the origin). A placement is a
/// substitution of x, y and z, which passes through min, max and negation,
/// so each part may be placed on its own. Its fixed params read as the
/// floats they are (`fixed_params`). None: no part is unchanged.
pub(super) fn split_about(v: &mut ModelView, f: u32, turn: &Xf, centre: [f64; 2], moved: [u32; 3]) -> Option<u32> {
    let f = fixed_params(v, f);
    let invariant = ink_core::symmetry::carried(&v.doc.g, f, centre, None).invariant;
    if invariant.is_empty() {
        return None;
    }
    let (z, o) = (v.doc.g.constant(0.0), v.doc.g.constant(1.0));
    let moved = Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: moved };
    let (tm, mm) = (v.doc.inverse_map(turn), v.doc.inverse_map(&moved));
    Some(split(&mut v.doc.g, f, &invariant, &tm, &mm))
}

/// `f` read through `tm`, its parts in `inv` (and what is made of them by min, max and negation)
/// through `mm`: each a substitution of x, y and z.
fn split(g: &mut ink_core::graph::Graph, f: u32, inv: &std::collections::HashSet<u32>, tm: &[(Var, u32)], mm: &[(Var, u32)]) -> u32 {
    fn go(g: &mut ink_core::graph::Graph, n: u32, inv: &std::collections::HashSet<u32>, tm: &[(Var, u32)], mm: &[(Var, u32)]) -> u32 {
        if inv.contains(&n) {
            return g.subst(n, mm);
        }
        match g.nodes[n as usize] {
            Node::Bin(op @ (BinOp::Min | BinOp::Max), a, b) => {
                let (a, b) = (go(g, a, inv, tm, mm), go(g, b, inv, tm, mm));
                g.bin(op, a, b)
            }
            Node::Un(UnOp::Neg, a) => {
                let a = go(g, a, inv, tm, mm);
                g.neg(a)
            }
            // The same value without its name or promise.
            Node::Tag(_, a) | Node::Lipschitz(a, ..) => go(g, a, inv, tm, mm),
            _ => g.subst(n, tm),
        }
    }
    go(g, f, inv, tm, mm)
}

/// Field `f` placed relative to a frame by `rel`, where the frame's origin is one point `k` of the
/// shape's own frame over the scope (`relative`: the shape does not turn in the world, the frame
/// turns about its own origin): read at R^T p + k, the frame's turn backwards and then one fixed
/// move, which is rel's own inverse R^T (p - T) in real numbers (T is R of -k) without the far
/// point's turn written out and taken back: a belt seen from the wheel at its far end is turned about
/// that wheel's centre, 260 mm from where the belt is written, and R^T (p - T) over a band of the
/// angle is loose by that length times the band, more than a mesh's room. Where the turn is about z
/// (`about_z`) and k is two floats, its parts every turn about k leaves unchanged are placed by the
/// move alone (`split_about`).
pub(super) fn about_own(v: &mut ModelView, f: u32, rel: &Xf, k: [u32; 3], about_z: bool) -> u32 {
    let g = &mut v.doc.g;
    let p = [Var::X, Var::Y, Var::Z].map(|x| g.var(x));
    let tm: Vec<(Var, u32)> = [Var::X, Var::Y, Var::Z]
        .into_iter()
        .enumerate()
        .map(|(i, var)| {
            let mut acc = g.constant(0.0);
            for (j, &pj) in p.iter().enumerate() {
                let m = g.mul(rel.r[j][i], pj);
                acc = g.add(acc, m);
            }
            (var, g.add(acc, k[i]))
        })
        .collect();
    let centre = (g.as_const(k[0]), g.as_const(k[1]));
    if let (true, (Some(x), Some(y))) = (about_z, centre) {
        let f = fixed_params(v, f);
        let invariant = ink_core::symmetry::carried(&v.doc.g, f, [x, y], None).invariant;
        if !invariant.is_empty() {
            let g = &mut v.doc.g;
            let mm = vec![(Var::X, g.add(p[0], k[0])), (Var::Y, g.add(p[1], k[1])), tm[2]];
            return split(g, f, &invariant, &tm, &mm);
        }
    }
    v.doc.g.subst(f, &tm)
}

/// A shape's placement relative to a frame, a turn about z, read as a
/// turn about a pivot: the point `l` of the shape's own frame (z 0) that
/// the motion holds within `delta` of the frame's point `c` at every
/// moment and param value of the scope. The shape so placed, p -> R (q -
/// l) + c, is the shape as placed moved by R l + T - c, one translation no
/// longer than `delta` (`charge`).
pub(super) struct Pivot {
    pub(super) l: [f64; 2],
    pub(super) c: [f64; 3],
    pub(super) delta: f64,
}

/// The pivot of `rel` (an isometry turning about z, its move not one
/// point over `span`), shown to drift less than `tol`: a linkage's
/// coupler seen from its rocker turns about the pin they share, though
/// the coupler is placed through its own angle and the pin through the
/// rocker's, so they meet only by the loop's closure, which no bound of
/// the two placements over a piece of the motion sees. The point is
/// fitted at the moments an estimate reads (claims::SAMPLES: the one
/// point two of their turns carry alike), taken at the centres the
/// shape's parts are written about (symmetry::centres, nearest on each
/// axis), so its parts round it are unchanged by the turn
/// (`split_about`); `c` is its image at the span's middle. Its drift
/// R l + T - c is bounded over pieces of the span, each t's band in the
/// view's scope (mean value bounds in it and the ranged params, what
/// solves on t set derived from it over the piece), the worst piece
/// halved until MAX_STEPS pieces: the largest bound, the sum of its axes'
/// magnitudes, rounded up. None: no turn over the span, no point fits,
/// or its drift is not shown under `tol`.
pub(super) fn pivot(v: &mut ModelView, f: u32, rel: &Xf, span: (f64, f64), tol: f64) -> Option<Pivot> {
    // Read with t as its param k, each piece of the span k's band in the
    // view's scope, so what solves on t set derived from it follows the
    // piece (Ctx::with_band).
    let k = v.t_param();
    let tk = v.doc.g.var(Var::Param(k));
    let rel = Xf { r: rel.r.map(|row| row.map(|n| v.doc.g.subst(n, &[(Var::T, tk)]))), t: rel.t.map(|n| v.doc.g.subst(n, &[(Var::T, tk)])) };
    let scope = v.ctx.clone();
    let over = |s: (f64, f64)| Ctx { t: ink_core::eval::span_mid(s), t_span: None, ..scope.with_band(k, s) };
    let tape = |v: &ModelView, n: u32| Tape::new(&v.doc.g, n);
    let (r, t) = (rel.r.map(|row| row.map(|n| tape(v, n))), rel.t.map(|n| tape(v, n)));
    let at = |s: f64| {
        let c = over((s, s));
        let e = |x: &Tape| ink_core::eval::eval_interval_scope(x, &c, [interval::I::point(0.0); 3], None).mid();
        (r.each_ref().map(|row| row.each_ref().map(e)), t.each_ref().map(e))
    };
    let n = claims::SAMPLES - 1;
    let moments: Vec<_> = (0..=n).map(|i| at(span.0 + (span.1 - span.0) * i as f64 / n as f64)).collect();
    // The two moments turned furthest apart: (R_0 - R_1) l = T_1 - T_0.
    let (r0, t0) = moments[0];
    let apart = |m: &([[f64; 3]; 3], [f64; 3])| (m.0[0][0] - r0[0][0]).hypot(m.0[1][0] - r0[1][0]);
    let (r1, t1) = *moments.iter().max_by(|x, y| apart(x).total_cmp(&apart(y)))?;
    let (a, b, c, d) = (r0[0][0] - r1[0][0], r0[0][1] - r1[0][1], r0[1][0] - r1[1][0], r0[1][1] - r1[1][1]);
    let det = a * d - b * c;
    if !(det.abs() > 0.0 && det.is_finite()) {
        return None;
    }
    let (u, w) = (t1[0] - t0[0], t1[1] - t0[1]);
    let fit = [(d * u - b * w) / det, (a * w - c * u) / det];
    let f = fixed_params(v, f);
    let centres = ink_core::symmetry::centres(&v.doc.g, f);
    let l: [f64; 2] = std::array::from_fn(|i| *centres[i].iter().min_by(|x, y| (*x - fit[i]).abs().total_cmp(&(*y - fit[i]).abs())).expect("0 is a centre"));
    let image = |m: &([[f64; 3]; 3], [f64; 3])| -> [f64; 3] { std::array::from_fn(|k| m.0[k][0] * l[0] + m.0[k][1] * l[1] + m.1[k]) };
    let c = image(&at(ink_core::eval::span_mid(span)));
    // Not near one point at the moments: no pivot, and no bound sought.
    if moments.iter().any(|m| (0..3).map(|k| (image(m)[k] - c[k]).abs()).sum::<f64>() >= tol) {
        return None;
    }
    let g = &mut v.doc.g;
    let drift: Vec<Tape> = (0..3)
        .map(|k| {
            let (lx, ly, ck) = (g.constant(l[0]), g.constant(l[1]), g.constant(c[k]));
            let (x, y) = (g.bin_real(BinOp::Mul, rel.r[k][0], lx), g.bin_real(BinOp::Mul, rel.r[k][1], ly));
            let s = g.bin_real(BinOp::Add, x, y);
            let s = g.bin_real(BinOp::Add, s, rel.t[k]);
            let e = g.bin_real(BinOp::Sub, s, ck);
            Tape::new(g, e)
        })
        .collect();
    let bound = |s: (f64, f64)| -> f64 {
        let at = over(s);
        let b = drift.iter().fold(0.0, |b, e| interval::add_up(b, interval::mag(ink_core::eval::eval_interval_scope(e, &at, [interval::I::point(0.0); 3], None))));
        if b.is_nan() { f64::INFINITY } else { b }
    };
    // Bounds are at least 0: their bits order as they do.
    let mut pieces = std::collections::BinaryHeap::from([(bound(span).to_bits(), span.0.to_bits(), span.1.to_bits())]);
    while pieces.len() < MAX_STEPS {
        let &(b, s0, s1) = pieces.peek().expect("a piece");
        let (s0, s1) = (f64::from_bits(s0), f64::from_bits(s1));
        let m = ink_core::eval::span_mid((s0, s1));
        if f64::from_bits(b) == 0.0 || !(s0 < m && m < s1) {
            break;
        }
        pieces.pop();
        for h in [(s0, m), (m, s1)] {
            pieces.push((bound(h).to_bits(), h.0.to_bits(), h.1.to_bits()));
        }
    }
    let delta = f64::from_bits(pieces.peek().expect("a piece").0);
    (delta < tol).then_some(Pivot { l, c, delta })
}

/// What a pair claim asked of shapes moved (each node of `moved` its
/// world shape moved by one translation no longer than its bound, the
/// rest exact) is charged to say the world's (`ranged_in` asks it so):
/// None where the moved shapes do not keep it within that. A clearance's
/// shapes each moved whole (an operand that is one moved shape or reads
/// none): the gap is off by at most the sum. Disjoint, subset and
/// overlaps ask whether a region of their shapes' points (a ∩ b; a less
/// b) holds a cube of half-side the accuracy: where it is an intersection
/// of parts each moved whole (`literal`), a cube of half-side s in the
/// world's region holds one of s - d in the region asked, d the largest
/// move, and a cube asked of s + d one of s in the world's at each moment.
/// Not where a moved part is one of a union: its members moving apart can
/// open a gap that closes deep inside the cube.
pub(super) fn charge(g: &ink_core::graph::Graph, claim: &Claim, moved: &std::collections::HashMap<u32, f64>) -> Option<f64> {
    let reads = |n: u32| g.reachable(&[n]).into_iter().any(|m| moved.contains_key(&m));
    // A name or a promise wraps the same shape.
    let whole = |mut n: u32| loop {
        match (moved.get(&n), g.nodes[n as usize]) {
            (Some(&d), _) => return Some(d),
            (None, Node::Tag(_, a) | Node::Lipschitz(a, ..)) => n = a,
            (None, _) => return (!reads(n)).then_some(0.0),
        }
    };
    // The largest move of the parts `n` (a region where `inside`, else its
    // complement) is the intersection of, each moved whole.
    fn literal(g: &ink_core::graph::Graph, n: u32, inside: bool, whole: &dyn Fn(u32) -> Option<f64>) -> Option<f64> {
        if let Some(d) = whole(n) {
            return Some(d);
        }
        match g.nodes[n as usize] {
            Node::Bin(BinOp::Max, a, b) if inside => Some(literal(g, a, inside, whole)?.max(literal(g, b, inside, whole)?)),
            Node::Bin(BinOp::Min, a, b) if !inside => Some(literal(g, a, inside, whole)?.max(literal(g, b, inside, whole)?)),
            Node::Un(UnOp::Neg, a) => literal(g, a, !inside, whole),
            Node::Tag(_, a) | Node::Lipschitz(a, ..) => literal(g, a, inside, whole),
            _ => None,
        }
    }
    match *claim {
        Claim::Clearance(a, b, _) => Some(interval::add_up(whole(a)?, whole(b)?)),
        Claim::Disjoint(a, b) | Claim::Overlaps(a, b) => Some(literal(g, a, true, &whole)?.max(literal(g, b, true, &whole)?)),
        Claim::Subset(a, b) => Some(literal(g, a, true, &whole)?.max(literal(g, b, false, &whole)?)),
        _ => None,
    }
}
