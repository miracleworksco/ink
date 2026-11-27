use super::{Group, M3, Model, Num, ROUNDING, V3, mul, mv, tr};
use ink_core::eval::Ctx;
use ink_core::graph::{Graph, NodeId, Var};
use ink_core::interval::{I, add_o, mul_o, sub_o};
use ink_core::v3::{cross, dot, norm};
use ink_kernel::rigid::Moments;
use ink_lang::Value;
use ink_lang::view::ModelView;
use std::collections::HashSet;

impl Group {
    /// Whether it moves through a solve: nothing of it is bounded over a
    /// step, so it is never read with its spin at rest (Model::round_enough).
    pub(super) fn solved(&self) -> bool {
        let numbers =
            self.frame.iter().chain(&self.dt).chain(&self.dt2).chain(self.dq.iter().flat_map(|(_, d)| d)).chain([&self.rate, &self.acc]).chain(self.along.iter().map(|(_, n)| n));
        numbers.into_iter().any(|n| matches!(n, Num::S(..)))
    }

    /// The most its departure from round, `(h, d2, d1)` (Model::departure),
    /// can move the effort along each coordinate over ctx's scope (a step):
    /// N mm or N per unit of the coordinate. Its points q in the frame at
    /// rest (spun by the angle) are at R (S q) + T; with the frame's
    /// angular velocity w, acceleration w', partial angular velocity a_k
    /// and partial velocity v_k, its origin's acceleration a_0 and gravity
    /// g, all in the frame, and the spin's rate, acceleration and partial
    /// along the axis e: nu = a_k + theta_k e, W = w + theta' e, L = w' +
    /// theta'' e + theta' w x e. The effort along k is linear in the mass's
    /// moments; the departure's part of it is dh . c - nu^T dJ L - W^T dJ
    /// (nu x W), c = (a_0 - g) x nu + v_k x L + W x (W x v_k), dh across e
    /// and dJ the second moment's part a turn about e moves. So it is at
    /// most h |c across| + d2 |x across| |y across| +
    /// d1 (|x across| |y_e| + |y across| |x_e|) for each of the two pairs (x, y), whatever the
    /// angle: from bounds over the scope (interval arithmetic). A frame that
    /// does not move gives h |g across| |theta_k|: a weight off a level
    /// shaft's axis.
    pub(super) fn moves(&self, (h, d2, d1): (f64, f64, f64), ctx: &Ctx, gravity: [f64; 3], nq: usize) -> Vec<f64> {
        use ink_core::interval::mag;
        // (Through a solve it has no bound over a step; at one moment it has its value.)
        if self.solved() && ctx.t_span.is_some() {
            return vec![f64::INFINITY; nq];
        }
        let m3 = |x: &[Num; 12]| -> [[I; 3]; 3] { std::array::from_fn(|p| std::array::from_fn(|q| x[3 * p + q].over(ctx))) };
        let v3 = |x: &[Num; 12]| -> [I; 3] { std::array::from_fn(|p| x[9 + p].over(ctx)) };
        let r = m3(&self.frame);
        let w = ivee(itm(&r, &m3(&self.dt)));
        let wd = ivee(itm(&r, &m3(&self.dt2)));
        let lin = isub(itv(&r, &v3(&self.dt2)), itv(&r, &gravity.map(|x| I::point(1000.0 * x))));
        let e: [I; 3] = std::array::from_fn(|k| I::point(f64::from(u8::from(k == self.axis))));
        let (rate, acc) = (self.rate.over(ctx), self.acc.over(ctx));
        let big = iadd(w, iscale(rate, e));
        let l = iadd(iadd(wd, iscale(acc, e)), iscale(rate, icross(w, e)));
        let (i, j, a) = ((self.axis + 1) % 3, (self.axis + 2) % 3, self.axis);
        let across = |x: [I; 3]| mag(x[i]).hypot(mag(x[j]));
        // A departure of 0 moves nothing, however the motion is bounded
        // (0 times an unbounded rate is 0, not NaN).
        let times = |d: f64, x: f64| if d == 0.0 { 0.0 } else { d * x };
        let pair = |x: [I; 3], y: [I; 3]| times(d2, across(x) * across(y)) + times(d1, across(x) * mag(y[a]) + across(y) * mag(x[a]));
        let zero = [I::point(0.0); 3];
        (0..nq)
            .map(|k| {
                let (ak, vk) = match self.dq.iter().find(|(l, _)| *l == k) {
                    Some((_, d)) => (ivee(itm(&r, &m3(d))), itv(&r, &v3(d))),
                    None => (zero, zero),
                };
                let tk = self.along.iter().find(|(l, _)| *l == k).map_or(I::point(0.0), |(_, n)| n.over(ctx));
                let nu = iadd(ak, iscale(tk, e));
                let c = iadd(iadd(icross(lin, nu), icross(vk, l)), icross(big, icross(big, vk)));
                // What no bound was found for (NaN) has none: infinite.
                let b = 1e-6 * (times(h, across(c)) + pair(nu, l) + pair(big, icross(nu, big)));
                if b.is_nan() { f64::INFINITY } else { b }
            })
            .collect()
    }
}

impl Group {
    /// The most its departure from round, `(h, d2, d1)` (Model::departure),
    /// can move the wrench on what carries it over ctx's scope (a step), at
    /// any angle of its spin: in the frame, with W, L and a_0 - g as in
    /// `moves`, the departure's first moment h_r (h, turned by the spin)
    /// moves the force by L x h_r + W x (W x h_r), at most h (|L| + |W|^2),
    /// and the moment about the frame's origin by h_r x (a_0 - g) - dS L - W
    /// x dS W, dS the part of the second moment a turn moves: its part in
    /// the plane across the axis (eigenvalues +-d2) acts on the parts of L
    /// and W across the axis alone, the rest (+-d1) on all of them, so at
    /// most h |a_0 - g| + d2 (|L across| + |W| |W across|) + d1 (|L| +
    /// |W|^2) (a steady spin on a still frame moves no moment by its d2):
    /// N and N mm. The weight and the origin's acceleration move no force
    /// (the departure has no mass).
    pub(super) fn shakes(&self, (h, d2, d1): (f64, f64, f64), ctx: &Ctx, gravity: [f64; 3]) -> (f64, f64) {
        use ink_core::interval::mag;
        if self.solved() && ctx.t_span.is_some() {
            return (f64::INFINITY, f64::INFINITY);
        }
        let m3 = |x: &[Num; 12]| -> [[I; 3]; 3] { std::array::from_fn(|p| std::array::from_fn(|q| x[3 * p + q].over(ctx))) };
        let v3 = |x: &[Num; 12]| -> [I; 3] { std::array::from_fn(|p| x[9 + p].over(ctx)) };
        let r = m3(&self.frame);
        let w = ivee(itm(&r, &m3(&self.dt)));
        let wd = ivee(itm(&r, &m3(&self.dt2)));
        let lin = isub(itv(&r, &v3(&self.dt2)), itv(&r, &gravity.map(|x| I::point(1000.0 * x))));
        let e: [I; 3] = std::array::from_fn(|k| I::point(f64::from(u8::from(k == self.axis))));
        let (rate, acc) = (self.rate.over(ctx), self.acc.over(ctx));
        let big = iadd(w, iscale(rate, e));
        let l = iadd(iadd(wd, iscale(acc, e)), iscale(rate, icross(w, e)));
        let size = |x: [I; 3]| x.iter().map(|c| mag(*c) * mag(*c)).sum::<f64>().sqrt();
        let (i, j) = ((self.axis + 1) % 3, (self.axis + 2) % 3);
        let across = |x: [I; 3]| mag(x[i]).hypot(mag(x[j]));
        let spin = size(l) + size(big) * size(big);
        // No departure moves nothing, however wide the bounds of the motion.
        let times = |a: f64, b: f64| if a == 0.0 { 0.0 } else { a * b };
        let moment = times(h, size(lin)) + times(d2, across(l) + size(big) * across(big)) + times(d1, spin);
        // What no bound was found for (NaN) has none: infinite.
        let none = |b: f64| if b.is_nan() { f64::INFINITY } else { b };
        (none(1e-6 * times(h, spin)), none(1e-6 * moment))
    }
}

/// Interval vectors and matrices, rounded outward: a + b, a - b, s a, a x b.
fn iadd(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| add_o(a[k], b[k]))
}

fn isub(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| sub_o(a[k], b[k]))
}

fn iscale(s: I, a: [I; 3]) -> [I; 3] {
    a.map(|x| mul_o(s, x))
}

fn icross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|k| {
        let (p, q) = ((k + 1) % 3, (k + 2) % 3);
        sub_o(mul_o(a[p], b[q]), mul_o(a[q], b[p]))
    })
}

/// a^T b, a^T v, and the axial vector of a matrix's skew part.
fn itm(a: &[[I; 3]; 3], b: &[[I; 3]; 3]) -> [[I; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).fold(I::point(0.0), |s, k| add_o(s, mul_o(a[k][i], b[k][j])))))
}

fn itv(a: &[[I; 3]; 3], v: &[I; 3]) -> [I; 3] {
    std::array::from_fn(|i| (0..3).fold(I::point(0.0), |s, k| add_o(s, mul_o(a[k][i], v[k]))))
}

fn ivee(w: [[I; 3]; 3]) -> [I; 3] {
    let half = I::point(0.5);
    [mul_o(half, sub_o(w[2][1], w[1][2])), mul_o(half, sub_o(w[0][2], w[2][0])), mul_o(half, sub_o(w[1][0], w[0][1]))]
}

/// A body's own spin: the lowest of its links that moves turns it, by a
/// joint, or by a WHERE that rides one (placed on a frame of what the
/// joint turns: `ride`) or is exactly a rotation about the child's own
/// axis through its origin (interfere's own turn, its translation free to
/// move); and the links below are fixed and a rotation to rounding.
/// Placed by X = Y S(angle), S the turn about the axis, its effort reads Y
/// and the angle's derivatives, and the angle only through its departure
/// from round about the axis. It is shown round where its mass seen in the
/// turn's frame is: its shape's and density's turn symmetry
/// (symmetry::turn_symmetry) every turn or of order 3 or more, which puts
/// its centre on the axis and gives it one inertia about every direction
/// across it, exactly; or, `balanced` (it has a spin-inertia), the axis is
/// its own z axis. Else it joins its link's group (Group).
pub(super) struct Spin {
    /// The body's placement with the spin at rest (Y): R row major, then T.
    pub(super) rest: [NodeId; 12],
    /// The spin's angle: a joint's drive (degrees), or the radians a
    /// WHERE's turn winds (rigid::winding); and radians per unit of it.
    pub(super) angle: NodeId,
    pub(super) radians: f64,
    /// The axis in the body's frame: a point on it, and its direction.
    pub(super) axis: (V3, V3),
    pub(super) round: bool,
    /// The link it spins with (its instance path); the frame that link
    /// turns in, at rest, in the world, and the axis of it the link turns
    /// about; the body's placement there (a body point p is R p + T).
    pub(super) link: String,
    pub(super) frame: ink_lang::Xf,
    pub(super) turn_axis: usize,
    pub(super) within: (M3, V3),
}

/// A turn joint whose drive moves: its instance path and joint, and, from
/// an item it places, the world placement of its child's frame (the links
/// down to it) and of the frame it turns in at rest (those above it, then
/// its WHERE).
pub(super) struct Turning {
    path: String,
    joint: ink_lang::Joint,
    child: ink_lang::Xf,
    frame: ink_lang::Xf,
}

pub(super) fn turning(v: &mut ModelView, items: &[ink_lang::Item], changes: &dyn Fn(&Graph, NodeId) -> bool) -> Vec<Turning> {
    let mut out: Vec<Turning> = Vec::new();
    for it in items {
        for (i, l) in it.links.iter().enumerate() {
            let Some(j) = l.joint.as_ref().filter(|j| j.turn && changes(&v.doc.g, j.drive) && !out.iter().any(|t| t.path == l.path)) else { continue };
            let mut above = v.doc.identity();
            for l in &it.links[..i] {
                above = v.doc.compose(&l.at, &above);
            }
            let child = v.doc.compose(&l.at, &above);
            let frame = v.doc.compose(&j.anchor, &above);
            out.push(Turning { path: l.path.clone(), joint: j.clone(), child, frame });
        }
    }
    out
}

/// What a spinning body spins with: the link's path, the axis it turns
/// about, the frame it turns in at rest, the body's placement there, the
/// angle and radians per unit of it.
type Turn = (String, usize, ink_lang::Xf, ink_lang::Xf, NodeId, f64);

pub(super) fn spin(
    v: &mut ModelView,
    it: &ink_lang::Item,
    density: Option<NodeId>,
    balanced: bool,
    changes: &dyn Fn(&Graph, NodeId) -> bool,
    turning: &[Turning],
    eps: &[u32],
    solved: &[u32],
) -> Option<Spin> {
    use ink_core::symmetry::{Order, turn_symmetry};
    let entries = |x: &ink_lang::Xf| -> [NodeId; 12] { std::array::from_fn(|k| if k < 9 { x.r[k / 3][k % 3] } else { x.t[k - 9] }) };
    let moves = |g: &Graph, x: &ink_lang::Xf| entries(x).iter().any(|&e| changes(g, e));
    let i = (0..it.links.len()).rev().find(|&k| moves(&v.doc.g, &it.links[k].at))?;
    let link = &it.links[i];
    let mut above = v.doc.identity();
    for l in &it.links[..i] {
        above = v.doc.compose(&l.at, &above);
    }
    // The body in its link's frame: the links below.
    let mut c = v.doc.identity();
    for l in it.links[i + 1..].iter().rev() {
        c = v.doc.compose(&c, &l.at);
    }
    // The turn: its axis, the frame it turns in at rest, the body there
    // (the links below, then `below`), its angle.
    let turn: Option<Turn> = match &link.joint {
        Some(j) if j.turn => Some((link.path.clone(), j.axis, v.doc.compose(&j.anchor, &above), v.doc.compose(&c, &j.below), j.drive, std::f64::consts::PI / 180.0)),
        Some(_) => return None,
        None => 'own: {
            if let Some(t) = ride(v, &it.instance.as_ref()?.1, turning, eps, solved) {
                break 'own Some(t);
            }
            let (g, r) = (&v.doc.g, link.at.r);
            // A turn through a solve: read as polynomials, a solved value
            // is a fixed number, so its angle is not shown.
            if r.iter().flatten().any(|&e| g.depends_on(e, |x| matches!(x, Var::Param(k) if solved.contains(&k)))) {
                break 'own None;
            }
            if !ink_kernel::rigid::turn_bounds(g, &r, &v.ctx).is_some_and(|b| b.skew == 0.0 && b.det.lo == 1.0 && b.det.hi == 1.0) {
                break 'own None;
            }
            let Some(k) = (0..3).find(|&k| {
                let (a, b) = ((k + 1) % 3, (k + 2) % 3);
                g.constant_value(r[k][k]) == Some(1.0) && [r[k][a], r[k][b], r[a][k], r[b][k]].iter().all(|&e| g.constant_value(e) == Some(0.0)) && changes(g, r[a][a])
            }) else {
                break 'own None;
            };
            let (a, b) = ((k + 1) % 3, (k + 2) % 3);
            let Some(wound) = ink_kernel::rigid::winding(g, r[a][a], r[b][a], &v.ctx) else { break 'own None };
            let id = v.doc.identity();
            let rest = v.doc.compose(&ink_lang::Xf { r: id.r, t: link.at.t }, &above);
            let g = &mut v.doc.g;
            let mut angle = g.constant(0.0);
            for (a, m) in wound {
                let m = g.constant(f64::from(m));
                let x = g.mul(m, a);
                angle = g.add(angle, x);
            }
            Some((link.path.clone(), k, rest, c, angle, 1.0))
        }
    };
    let (path, axis, frame, c, angle, radians) = turn?;
    if moves(&v.doc.g, &c) || !ink_kernel::rigid::turn_bounds(&v.doc.g, &c.r, &v.ctx).is_some_and(|b| b.skew <= ROUNDING) {
        return None;
    }
    let at = |n: NodeId| v.doc.scalar(n, &v.ctx);
    let (r, t): (M3, V3) = (std::array::from_fn(|i| std::array::from_fn(|j| at(c.r[i][j]))), std::array::from_fn(|i| at(c.t[i])));
    let point = mv(&tr(&r), &t).map(|x| -x);
    // A body with a published spin-inertia is a balanced rotor about its
    // own z axis: turned about that line, it is round.
    let near = ROUNDING * (1.0 + norm(t));
    let balanced_about = balanced && [r[axis][0], r[axis][1], point[0], point[1]].iter().all(|x| x.abs() <= near);
    // Its shape and density there, the axis turned onto z.
    let [x, y, z] = [Var::X, Var::Y, Var::Z].map(|k| v.doc.g.var(k));
    let onto_z = match axis {
        0 => [(Var::X, z), (Var::Y, x), (Var::Z, y)],
        1 => [(Var::X, y), (Var::Y, z), (Var::Z, x)],
        _ => [(Var::X, x), (Var::Y, y), (Var::Z, z)],
    };
    let mut round = |f: NodeId| {
        let Value::Field(h, _) = v.doc.to_world(&Value::Field(f, None), &c) else { return false };
        let h = v.doc.g.subst(h, &onto_z);
        matches!(turn_symmetry(&v.doc.g, h).order, Order::Every | Order::Of(3..))
    };
    let round = balanced_about || (round(it.shape) && density.is_none_or(round));
    let rest = entries(&v.doc.compose(&c, &frame));
    Some(Spin { rest, angle, radians, axis: (point, r[axis]), round, link: path, frame, turn_axis: axis, within: (r, t) })
}

/// The turn joint a body rides, placed on a frame of what it turns (a
/// pulley on a motor's shaft frame, its grub screw on the pulley's): one
/// whose drive its placement `x` reads, and whose child's frame sees it at
/// one placement to within rounding whatever t and the coordinates
/// (rigid::nearly_constant, the coordinates' params read as variables).
fn ride(v: &mut ModelView, x: &ink_lang::Xf, turning: &[Turning], eps: &[u32], solved: &[u32]) -> Option<Turn> {
    let entries = |x: &ink_lang::Xf| -> Vec<NodeId> { x.r.iter().flatten().chain(&x.t).copied().collect() };
    // Read as polynomials a solved value is a fixed number: a placement
    // through a solve is not shown to ride.
    if entries(x).iter().any(|&e| v.doc.g.depends_on(e, |k| matches!(k, Var::Param(k) if solved.contains(&k)))) {
        return None;
    }
    let reach: HashSet<NodeId> = v.doc.g.reachable(&entries(x)).into_iter().collect();
    let over = Ctx { p_spans: ink_core::eval::ParamSpans::new(eps.iter().map(|&k| (k, (v.ctx.params[k as usize], v.ctx.params[k as usize])))), ..v.ctx.clone() };
    for t in turning.iter().filter(|t| reach.contains(&t.joint.drive)) {
        let back = v.doc.inverse(&t.child);
        let seen = v.doc.compose(x, &back);
        let Some(fixed) = entries(&seen).iter().map(|&e| ink_kernel::rigid::nearly_constant(&v.doc.g, e, &over)).collect::<Option<Vec<(f64, f64)>>>() else { continue };
        let size = 1.0 + fixed[9..].iter().map(|f| f.0 * f.0).sum::<f64>().sqrt();
        if fixed.iter().enumerate().any(|(k, f)| f.1 > ROUNDING * if k < 9 { 1.0 } else { size }) {
            continue;
        }
        let g = &mut v.doc.g;
        let k: Vec<NodeId> = fixed.iter().map(|f| g.constant(f.0)).collect();
        let seen = ink_lang::Xf { r: std::array::from_fn(|i| std::array::from_fn(|j| k[3 * i + j])), t: [k[9], k[10], k[11]] };
        return Some((t.path.clone(), t.joint.axis, t.frame, v.doc.compose(&seen, &t.joint.below), t.joint.drive, std::f64::consts::PI / 180.0));
    }
    None
}

/// A spinning body's centre and inertia about it (in its frame) where its shape repeats about a
/// line of its own along its spin axis (`axis`: a point of the spin's, and its direction, which
/// must be one of the body's coordinates), but for small parts (symmetry::nearly_about): the part
/// that repeats has its centre on its own line exactly and, about it, no product of inertia with
/// the axis (and the same inertia every way across it, where it repeats three times or more a
/// turn), whatever its measure reads; the small parts are as measured. With them, what its first
/// and second moments are known to: what the small parts' own measures are known to (0 with none).
/// None: its shape or density is not shown so.
pub(super) fn own_round(v: &mut ModelView, shape: NodeId, rho: Option<NodeId>, m: &ink_kernel::query::Measured, (_, u): (V3, V3)) -> Option<(V3, M3, f64, f64, bool)> {
    use ink_core::symmetry::{Order, nearly_about};
    let k = (0..3).find(|&k| (u[k].abs() - 1.0).abs() <= ROUNDING)?;
    if rho.is_some_and(|r| v.doc.g.depends_on(r, |x| matches!(x, Var::X | Var::Y | Var::Z))) {
        return None;
    }
    let near = nearly_about(&mut v.doc.g, shape, k)?;
    let every = match near.symmetry.order {
        Order::Every | Order::Of(3..) => true,
        Order::Of(2) => false,
        Order::Of(_) => return None,
    };
    // The core's moments: the body's, less what is added and with what is cut.
    let times = |m: Moments, k: f64| Moments { about: m.about, m: m.m.map(|x| k * x) };
    let mut core = m.moments();
    let mut parts: Vec<(Moments, f64)> = Vec::new();
    let (mut first, mut second) = (0.0, 0.0);
    for (part, sign) in [(near.added, 1.0), (near.cut, -1.0)] {
        let Some(part) = part else { continue };
        // A part with no box is no part.
        let Some(pm) = v.inertia(part, rho, None) else { continue };
        let of = pm.moments();
        core = core.plus(&times(of, -sign));
        // What the part's own measure is known to, as a whole body's is
        // (its proven bounds, Measured::error), carried to the line's point.
        let size = |r: std::ops::Range<usize>| r.map(|j| pm.error(j).powi(2)).sum::<f64>().sqrt();
        let (d, mass) = (norm(std::array::from_fn(|k| of.about[k] - near.at[k])), pm.error(0).abs());
        first += size(1..4) + d * mass;
        second += size(4..10) + 2.0 * d * size(1..4) + d * d * mass;
        parts.push((of, sign));
    }
    let p = ink_kernel::rigid::mass_props(&core)?;
    let mut line = [0.0; 3];
    line[k] = 1.0;
    // Its centre onto its own line, its inertia as a shape so repeating has it.
    let along = (0..3).map(|j| (p.center[j] - near.at[j]) * line[j]).sum::<f64>();
    let centre: V3 = std::array::from_fn(|j| near.at[j] + along * line[j]);
    let inertia = match every {
        true => about_axis(p.mass, centre, &p.inertia, (near.at, line)).1,
        false => std::array::from_fn(|a| std::array::from_fn(|b| if (a == k) != (b == k) { 0.0 } else { p.inertia[a][b] })),
    };
    // Back to moments about its centre, with the small parts as measured.
    let half = 0.5 * (inertia[0][0] + inertia[1][1] + inertia[2][2]);
    let about = |a: usize, b: usize| (if a == b { half } else { 0.0 }) - inertia[a][b];
    let mut whole = Moments { about: centre, m: [p.mass, 0.0, 0.0, 0.0, about(0, 0), about(1, 1), about(2, 2), about(0, 1), about(1, 2), about(2, 0)] };
    for (of, sign) in parts {
        whole = whole.plus(&times(of, sign));
    }
    let p = ink_kernel::rigid::mass_props(&whole)?;
    Some((p.center, p.inertia, first, second + 2.0 * norm(std::array::from_fn(|k| p.center[k] - near.at[k])) * first, every))
}

/// An inertia's part round about unit u (the mean over a turn about it), about the same point.
pub(super) fn round_part(inertia: &M3, u: V3) -> M3 {
    let (iu, trace) = (dot(u, mv(inertia, &u)), inertia[0][0] + inertia[1][1] + inertia[2][2]);
    std::array::from_fn(|a| std::array::from_fn(|b| iu * u[a] * u[b] + 0.5 * (trace - iu) * (if a == b { 1.0 } else { 0.0 } - u[a] * u[b])))
}

/// Mass properties averaged over every turn about the line through o along
/// unit u (all in the body's frame): a body's centre and inertia about it,
/// as a body round about that line has them. The centre comes onto the
/// line; the inertia keeps its moment about u and takes the mean of the
/// others across it; the ring the centre's offset sweeps adds m r^2 (1 +
/// u u^T) / 2. The moment about the line is unchanged.
pub(super) fn about_axis(mass: f64, c: V3, inertia: &M3, (o, u): (V3, V3)) -> (V3, M3) {
    let d: V3 = std::array::from_fn(|k| c[k] - o[k]);
    let along = dot(d, u);
    let r: V3 = std::array::from_fn(|k| d[k] - along * u[k]);
    let (iu, trace) = (dot(u, mv(inertia, &u)), inertia[0][0] + inertia[1][1] + inertia[2][2]);
    let ring = mass * dot(r, r);
    let i = std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            let (uu, id) = (u[a] * u[b], if a == b { 1.0 } else { 0.0 });
            iu * uu + 0.5 * (trace - iu) * (id - uu) + 0.5 * ring * (id + uu)
        })
    });
    (std::array::from_fn(|k| o[k] + along * u[k]), i)
}

/// A spinning body's mass properties (its centre, its inertia about it, in
/// its frame) with each part of its departure from round about the line
/// through o along unit u that its measure does not know made round, as a
/// round body's are (about_axis): its first moment across the line, within
/// what its first moments are known to (`known[0]`, g mm: twice what the
/// last refinement of its measure moved them), its centre onto the line;
/// its second moment's parts about its centre a turn moves twice ((xx -
/// yy)/2 and xy across the line) and once (xz and yz), within what they
/// are known to there (`known[1]`, carried to the centre, g mm^2), the mean of a turn. Each within the rounding of it too
/// (ROUNDING of its reach). A measure's error is no imbalance: a wheel
/// measured off its axis by a millionth of its size (its centre known to a
/// ten-thousandth) spins round, and its reaction reads no m |c| w^2 from it.
/// Whether a part that was not 0 was made round: of its first moment, of its second.
pub(super) fn known_round(mass: f64, c: V3, inertia: &M3, (o, u): (V3, V3), known: [f64; 2]) -> (V3, M3, [bool; 2]) {
    let d: V3 = std::array::from_fn(|k| c[k] - o[k]);
    let along = dot(d, u);
    let r: V3 = std::array::from_fn(|k| d[k] - along * u[k]);
    // Its second moment about its centre, and two directions across u.
    let half = 0.5 * (inertia[0][0] + inertia[1][1] + inertia[2][2]);
    let mut s: M3 = std::array::from_fn(|p| std::array::from_fn(|q| if p == q { half } else { 0.0 } - inertia[p][q]));
    let spread = s[0][0] + s[1][1] + s[2][2];
    let k = (0..3).min_by(|&a, &b| u[a].abs().total_cmp(&u[b].abs())).unwrap_or(0);
    let e1 = {
        let x: V3 = std::array::from_fn(|i| f64::from(u8::from(i == k)) - u[k] * u[i]);
        let n = norm(x);
        x.map(|v| v / n)
    };
    let e2 = cross(u, e1);
    let sa = |a: V3, b: V3| dot(a, mv(&s, &b));
    let (twice, once) = ([0.5 * (sa(e1, e1) - sa(e2, e2)), sa(e1, e2)], [sa(e1, u), sa(e2, u)]);
    let (reach1, reach2) = (mass * norm(d) + (mass * spread.max(0.0)).sqrt(), spread.abs() + mass * dot(d, d));
    let (known1, known2) = (known[0], known[1]);
    let within = |v: f64, reach: f64, known: f64| v <= (ROUNDING * reach).max(known);
    let mut made = [false; 2];
    let mut center = c;
    let h = mass * norm(r);
    if within(h, reach1, known1) {
        made[0] |= h > 0.0;
        center = std::array::from_fn(|i| o[i] + along * u[i]);
    }
    // Take a part out of the second moment: x (a b^T + b a^T) for each of
    // its components x along the pair (a, b).
    let mut take = |x: f64, a: V3, b: V3| {
        (0..3).for_each(|p| (0..3).for_each(|q| s[p][q] -= x * (a[p] * b[q] + b[p] * a[q])));
    };
    if within(twice[0].hypot(twice[1]), reach2, known2) {
        made[1] |= twice != [0.0, 0.0];
        // (xx - yy)/2 off e1 e1^T and onto e2 e2^T, xy off both crosses.
        take(0.5 * twice[0], e1, e1);
        take(-0.5 * twice[0], e2, e2);
        take(twice[1], e1, e2);
    }
    if within(once[0].hypot(once[1]), reach2, known2) {
        made[1] |= once != [0.0, 0.0];
        take(once[0], e1, u);
        take(once[1], e2, u);
    }
    // The inertia from the second moment: trace(S) - S (S's trace kept).
    let trace = s[0][0] + s[1][1] + s[2][2];
    (center, std::array::from_fn(|p| std::array::from_fn(|q| if p == q { trace } else { 0.0 } - s[p][q])), made)
}

/// A published moment j about a body's own z axis (a rotor's, whose
/// catalog body is its shaft) in place of the measured one about it, the
/// centre's offset from the axis counted: the mass and the rest stay.
pub(super) fn spun(mass: f64, center: &V3, inertia: &mut M3, j: f64) {
    inertia[2][2] = j - mass * (center[0] * center[0] + center[1] * center[1]);
}

impl Model {
    /// A group's departure from round that is not known: its members' parts read as none within what
    /// their measures are known to (`Body::unknown`), summed as the largest they may be, in `departure`'s
    /// terms. It may lie at any angle, so it is carried as a departure read at rest is, by its bound
    /// (`Group::moves`, `Group::shakes`), at every moment, whether the group's spin is resolved or not.
    fn unknown(&self, g: &Group) -> (f64, f64, f64) {
        g.members.iter().fold((0.0, 0.0, 0.0), |(h, d2, d1), &(_, b, _, _)| {
            let u = self.bodies[b].unknown;
            (h + u[0], d2 + u[1], d1 + u[1])
        })
    }

    /// What the unknown departures (`unknown`) of the groups reaction `r` carries may move it by over
    /// `ctx`'s scope (a moment, or a step): each one's force, its moment about its frame's origin, and
    /// that origin, as `shaken` gives a group read at rest.
    pub(super) fn unknown_shakes(&self, r: usize, ctx: &Ctx, at: &Ctx) -> Vec<(f64, f64, V3)> {
        (0..self.groups.len())
            .filter(|&g| self.carries(r, g))
            .filter_map(|g| {
                let gr = &self.groups[g];
                let dep = self.unknown(gr);
                (dep != (0.0, 0.0, 0.0)).then(|| {
                    // Over the scope where the motion has a bound over it; else at the moment itself, as
                    // the reaction is read (a simulated or solved motion has no bound over a step).
                    let (f, mo) = Some(gr.shakes(dep, ctx, self.gravity)).filter(|(f, mo)| f.is_finite() && mo.is_finite()).unwrap_or_else(|| gr.shakes(dep, at, self.gravity));
                    (f, mo, std::array::from_fn(|k| gr.frame[9 + k].at(at)))
                })
            })
            .collect()
    }

    /// `unknown_shakes` at moment t alone.
    pub fn unknown_at(&self, r: usize, t: f64) -> Vec<(f64, f64, V3)> {
        let ctx = Ctx::new(t, self.params.clone());
        self.unknown_shakes(r, &ctx, &ctx)
    }

    /// What the unknown departures of every group may move the effort along each coordinate by over
    /// `ctx`'s scope (`Group::moves`), summed.
    pub(super) fn unknown_moves(&self, ctx: &Ctx) -> Vec<f64> {
        let nq = self.coords.len();
        // (At the scope's own moment where the motion has no bound over it, as `unknown_shakes`.)
        let at = Ctx { t_span: None, ..ctx.clone() };
        self.groups.iter().fold(vec![0.0; nq], |sum, gr| {
            let dep = self.unknown(gr);
            if dep == (0.0, 0.0, 0.0) {
                return sum;
            }
            let over = Some(gr.moves(dep, ctx, self.gravity, nq)).filter(|b| b.iter().all(|x| x.is_finite())).unwrap_or_else(|| gr.moves(dep, &at, self.gravity, nq));
            sum.iter().zip(over).map(|(x, y)| x + y).collect()
        })
    }

    /// Whether reaction r carries group g.
    pub(super) fn carries(&self, r: usize, g: usize) -> bool {
        self.groups[g].members.iter().any(|m| self.bodies[m.1].of.contains(&r))
    }

    /// Over each step, what each group in `which` can move the effort
    /// along each coordinate (Group::moves): [step][group][coordinate].
    pub(super) fn bounds(&self, steps: &[(f64, f64)], which: &[bool]) -> Vec<Vec<Vec<f64>>> {
        use rayon::prelude::*;
        let ctx = Ctx::new(0.0, self.params.clone());
        let nq = self.coords.len();
        let dep: Vec<(f64, f64, f64)> = self.groups.iter().map(|g| self.departure(g)).collect();
        steps
            .par_iter()
            .map(|&s| {
                let over = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..ctx.clone() };
                self.groups.iter().zip(&dep).zip(which).filter(|(_, w)| **w).map(|((g, d), _)| g.moves(*d, &over, self.gravity, nq)).collect()
            })
            .collect()
    }

    /// A group's departure from round about its axis, its members' masses
    /// summed in its frame: the first moment across the axis (g mm), and
    /// the parts of the second moment about the frame's origin a turn
    /// about the axis turns twice and once (g mm^2): of ((xx - yy)/2, xy)
    /// and of (xz, yz), the axis z. Averaged over every turn they are 0;
    /// what is left is the round mass the effort would read at every
    /// angle. A member's part within what its measure is known to is none
    /// already (known_round); one within the rounding of it (ROUNDING of
    /// its mass at its reach, of its second moment) is none, and so is
    /// their sum
    /// within ROUNDING of the parts summed: 0, which moves nothing however
    /// its spin is bounded (Group::moves, Group::shakes). A light member's
    /// real imbalance is kept beside a heavy one.
    pub(super) fn departure(&self, g: &Group) -> (f64, f64, f64) {
        let (i, j, a) = ((g.axis + 1) % 3, (g.axis + 2) % 3, g.axis);
        let (mut sum, mut summed) = ([[0.0f64; 2]; 3], [0.0f64; 3]);
        for &(_, b, r, t) in &g.members {
            let b = &self.bodies[b];
            let c: V3 = std::array::from_fn(|k| mv(&r, &b.center)[k] + t[k]);
            // Its second moment about its centre, from its inertia there,
            // and about the frame's origin.
            let half = 0.5 * (b.inertia[0][0] + b.inertia[1][1] + b.inertia[2][2]);
            let sc: M3 = std::array::from_fn(|p| std::array::from_fn(|q| if p == q { half } else { 0.0 } - b.inertia[p][q]));
            let sw = mul(&mul(&r, &sc), &tr(&r));
            let so = |p: usize, q: usize| sw[p][q] + b.mass * c[p] * c[q];
            let spread = sc[0][0] + sc[1][1] + sc[2][2];
            let (first, second) = (b.mass * norm(c) + (b.mass * spread.max(0.0)).sqrt(), spread.abs() + b.mass * dot(c, c));
            let parts = [([b.mass * c[i], b.mass * c[j]], first), ([0.5 * (so(i, i) - so(j, j)), so(i, j)], second), ([so(i, a), so(j, a)], second)];
            for (k, (v, reach)) in parts.into_iter().enumerate() {
                let m = v[0].hypot(v[1]);
                if m > ROUNDING * reach {
                    sum[k] = [sum[k][0] + v[0], sum[k][1] + v[1]];
                    summed[k] += m;
                }
            }
        }
        let past = |k: usize| {
            let x = sum[k][0].hypot(sum[k][1]);
            if x > ROUNDING * summed[k] { x } else { 0.0 }
        };
        (past(0), past(1), past(2))
    }
}
