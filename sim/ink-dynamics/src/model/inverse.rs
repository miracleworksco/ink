use super::couplings::unfed_need;
use super::drives::effort_split;
use super::spin::about_axis;
use super::{Coupled, M3, Model, Moment, Num, ROUNDING, Rides, Strip, V3, Wrench, dot_n, mul, mv, rot, split, tr, vee};
use ink_core::eval::Ctx;
use ink_core::v3::{cross, dot, norm};
use std::sync::Arc;

/// What the bodies sum into at a moment: the effort along each coordinate
/// (N mm per unit of it), its rounding's scale, and the generalized
/// inertia (g mm^2 per unit^2) and its diagonal's heft (Moment::heft;
/// None: not kept); each reaction's wrench (Wrench: force N, moment N mm
/// about the world's origin) and the scales of its rounding; gravity in
/// mm/s^2.
struct Sums<'a> {
    f: &'a mut [f64],
    scale: &'a mut [f64],
    gm: &'a mut [Vec<f64>],
    heft: Option<&'a mut [f64]>,
    wr: &'a mut [Wrench],
    gv: V3,
}

/// A share `x` at unit mass, at mass `err`: none where it is none, whatever
/// `err` (infinite: a mass not known).
fn by(err: f64, x: f64) -> f64 {
    if x == 0.0 { 0.0 } else { err * x }
}

/// Magnitudes of the terms actually projected onto this coordinate.
fn projected_size(direction: V3, magnitude: V3) -> f64 {
    dot(direction.map(f64::abs), magnitude.map(f64::abs))
}

impl Sums<'_> {
    /// What the rest of the model exerts on a body of the reactions `of`:
    /// its force and moment about the world's origin, at x, from its
    /// m (a - g) and I alpha + omega x I omega (sizes: their terms').
    fn react(&mut self, of: &[usize], x: V3, force: V3, torque: V3, size: (f64, f64)) {
        let moment = cross(x, force);
        for &r in of {
            let w = &mut self.wr[r];
            for k in 0..3 {
                w.force[k] += force[k];
                w.moment[k] += torque[k] + moment[k];
            }
            w.noise.0 += size.0;
            w.noise.1 += size.1 + norm(x) * size.0;
        }
    }

    /// A point mass's share: its mass (g), position, acceleration, partial
    /// velocities along the coordinates that move it, and the reactions
    /// it is of.
    fn point(&mut self, mass: f64, x: V3, a: V3, dq: &[(usize, V3)], of: &[usize]) {
        let force: V3 = std::array::from_fn(|k| 1e-6 * mass * (a[k] - self.gv[k]));
        let size = 1e-6 * mass * (norm(a) + norm(self.gv));
        for &(k, jv) in dq {
            self.f[k] += dot(jv, force);
            self.scale[k] += 1e-6 * mass * projected_size(jv, std::array::from_fn(|i| a[i].abs() + self.gv[i].abs()));
            for &(l, kv) in dq {
                self.gm[k][l] += mass * dot(jv, kv);
            }
            if let Some(h) = self.heft.as_deref_mut() {
                h[k] += mass * dot(jv, jv);
            }
        }
        self.react(of, x, force, [0.0; 3], (size, 0.0));
    }

    /// A point mass's inertia (its mass, g, at x with acceleration a) on
    /// reactions that carry what holds it, each by its weight: a belt's
    /// material on the parts its span or wrap is on (Model::belts).
    fn held(&mut self, mass: f64, x: V3, a: V3, on: &[(usize, f64)]) {
        let force: V3 = std::array::from_fn(|k| 1e-6 * mass * (a[k] - self.gv[k]));
        let size = 1e-6 * mass * (norm(a) + norm(self.gv));
        for &(r, w) in on {
            self.react(&[r], x, force.map(|f| w * f), [0.0; 3], (w.abs() * size, 0.0));
        }
    }

    /// A load's share: its force (N) and torque (N mm) in the world, the
    /// point it acts at, and that point's partial velocities and angular
    /// ones along the coordinates (those of what it acts on): an effort
    /// it supplies is one the actuators need not.
    fn load(&mut self, force: V3, torque: V3, x: V3, dq: &[(usize, V3, V3)], of: &[usize]) {
        for &(k, jv, jw) in dq {
            self.f[k] -= dot(jv, force) + dot(jw, torque);
            self.scale[k] += projected_size(jv, force) + projected_size(jw, torque);
        }
        self.react(of, x, force.map(|v| -v), torque.map(|v| -v), (norm(force), norm(torque)));
    }

    /// A rigid body's share: its mass (g), centre and inertia about it in
    /// its frame, at the moment its placement (Posed), and the reactions it
    /// is of.
    fn rigid(&mut self, posed: &Posed, mass: f64, c: V3, inertia: &M3, of: &[usize]) {
        let Posed { r, rd, rdd, pdd, w, .. } = posed;
        let Share { x, a, iw, force, torque, partial } = posed.share(mass, c, inertia, self.gv);
        let (w, al) = (*w, posed.al);
        let diagonal = (0..3).map(|i| iw[i][i].abs()).sum::<f64>();
        let big = iw.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
        let size = (1e-6 * mass * (norm(a) + norm(self.gv)), 1e-6 * big * (norm(al) + norm(w) * norm(w)));
        // A steady spin's angular acceleration is a cancellation inside
        // R'' Rᵀ + R' R'ᵀ. Keep the magnitudes of those products along
        // each axis before cancellation, just as for the final projection.
        let al_terms = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            0.5 * (projected_size(rdd[j], r[k]) + projected_size(rdd[k], r[j])) + projected_size(rd[j], rd[k])
        });
        let force_terms = std::array::from_fn(|i| 1e-6 * mass * (projected_size(rdd[i], c) + pdd[i].abs() + self.gv[i].abs()));
        let angular = iw.map(|row| projected_size(row, al_terms));
        let momentum = iw.map(|row| projected_size(row, w));
        let torque_terms = std::array::from_fn(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            1e-6 * (angular[i] + w[j].abs() * momentum[k] + w[k].abs() * momentum[j])
        });
        let turned: Vec<V3> = partial.iter().map(|(_, _, kw)| mv(&iw, kw)).collect();
        for &(k, jv, jw) in &partial {
            self.f[k] += dot(jv, force) + dot(jw, torque);
            self.scale[k] += projected_size(jv, force_terms) + projected_size(jw, torque_terms);
            for (&(l, kv, _), iwkw) in partial.iter().zip(&turned) {
                self.gm[k][l] += mass * dot(jv, kv) + dot(jw, *iwkw);
            }
            if let Some(h) = self.heft.as_deref_mut() {
                h[k] += mass * dot(jv, jv) + dot(jw, jw) * diagonal;
            }
        }
        self.react(of, x, force, torque, size);
    }
}

/// A rigid body's placement at a moment, as Sums::rigid reads it, kept
/// for its share at other masses, centres and inertias (what each of its
/// measured moments changes, Body::measure): the terms that rest on the
/// placement alone, worked out once.
struct Posed<'a> {
    r: M3,
    rd: M3,
    rt: M3,
    rdd: M3,
    pdd: V3,
    p: &'a [f64],
    w: V3,
    al: V3,
    /// Each partial: its coordinate, its rotation's and translation's
    /// partials, and its angular partial velocity.
    partial: Vec<(usize, M3, V3, V3)>,
}

impl<'a> Posed<'a> {
    fn new(r: &[f64], dr: &[f64], p: &'a [f64], dt2: &[f64], dq: &[(usize, &[f64])]) -> Self {
        let (r, rd) = (rot(r), rot(dr));
        let (rdd, pdd) = split(dt2);
        let rt = tr(&r);
        let w = vee(&mul(&rd, &rt));
        let al = {
            let (x, y) = (mul(&rdd, &rt), mul(&rd, &tr(&rd)));
            vee(&std::array::from_fn(|i| std::array::from_fn(|j| x[i][j] + y[i][j])))
        };
        let partial = dq
            .iter()
            .map(|(k, d)| {
                let (dr, dp) = split(d);
                (*k, dr, dp, vee(&mul(&dr, &rt)))
            })
            .collect();
        Posed { r, rd, rt, rdd, pdd, p, w, al, partial }
    }

    /// The body's share at this mass (g), centre and inertia in its frame,
    /// under gravity `gv` (mm/s^2): its centre's place and acceleration,
    /// its inertia in the world's axes, the force and torque it needs (N,
    /// N mm: g mm/s^2 and g mm^2/s^2 are 1e-6 of them), and its partial
    /// velocities and angular ones along each coordinate.
    fn share(&self, mass: f64, c: V3, inertia: &M3, gv: V3) -> Share {
        let Posed { r, rt, rdd, pdd, p, w, al, .. } = self;
        let (rc, rddc) = (mv(r, &c), mv(rdd, &c));
        let x: V3 = std::array::from_fn(|k| rc[k] + p[k]);
        let a: V3 = std::array::from_fn(|k| rddc[k] + pdd[k]);
        let iw = mul(&mul(r, inertia), rt);
        let force: V3 = std::array::from_fn(|k| 1e-6 * mass * (a[k] - gv[k]));
        let spin = cross(*w, mv(&iw, w));
        let ial = mv(&iw, al);
        let torque: V3 = std::array::from_fn(|k| 1e-6 * (ial[k] + spin[k]));
        let partial = self
            .partial
            .iter()
            .map(|&(k, dr, dp, jw)| {
                let drc = mv(&dr, &c);
                (k, std::array::from_fn(|i| drc[i] + dp[i]), jw)
            })
            .collect();
        Share { x, a, iw, force, torque, partial }
    }

    /// Sums::rigid's effort along each coordinate (added to `f`) and force
    /// and moment about the world's origin of the body alone at this mass
    /// (g), centre and inertia, worked as it works them: the same numbers,
    /// without the inertia, the heft and the roundings' scales, which no
    /// change of a measure reads.
    fn alone(&self, mass: f64, c: V3, inertia: &M3, gv: V3, f: &mut [f64]) -> (V3, V3) {
        let Share { x, force, torque, partial, .. } = self.share(mass, c, inertia, gv);
        for (k, jv, jw) in partial {
            f[k] += dot(jv, force) + dot(jw, torque);
        }
        // (As Sums::react sums them into a wrench at none.)
        let moment = cross(x, force);
        (std::array::from_fn(|k| 0.0 + force[k]), std::array::from_fn(|k| 0.0 + (torque[k] + moment[k])))
    }
}

/// A rigid body's share at a moment (Posed::share).
struct Share {
    x: V3,
    a: V3,
    iw: M3,
    force: V3,
    torque: V3,
    partial: Vec<(usize, V3, V3)>,
}

/// Gauss-Legendre nodes on -1..1 and weights: 2 integrate a line's
/// quadratic exactly (a point's acceleration and partials are affine along
/// it); 4 an arc's trig terms of degree 2 over a quarter turn to 5e-6 of the
/// arc's term.
const GAUSS2: [(f64, f64); 2] = [(-0.5773502691896257, 1.0), (0.5773502691896257, 1.0)];
const GAUSS4: [(f64, f64); 4] =
    [(-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461), (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538)];

/// A body's placement at a moment (Placement's entries there): its rotation
/// and translation, their partials along each coordinate and each joint
/// released for springs (by its row, after the coordinates: its unit
/// motion), and their first and second rates along t.
struct Placed {
    r: Vec<f64>,
    p: Vec<f64>,
    dq: Vec<(usize, Vec<f64>)>,
    dr: Vec<f64>,
    dt2: Vec<f64>,
}

impl Model {
    fn placed(&self, b: &super::Body, ctx: &Ctx, frames: &[(V3, V3, bool)]) -> Placed {
        let vals = |n: &[Num]| -> Vec<f64> { n.iter().map(|e| e.at(ctx)).collect() };
        let at = &b.at;
        let (r, p) = (vals(&at.r), vals(&at.p));
        let mut dq: Vec<(usize, Vec<f64>)> = at.dq.iter().map(|(k, d)| (*k, vals(d))).collect();
        dq.extend(b.released.iter().map(|&i| (self.coords.len() + i, self.unit_motion(frames, i, &r, &p).to_vec())));
        Placed { r, p, dq, dr: vals(&at.dr), dt2: vals(&at.dt2) }
    }

    /// At a state of a simulation (Model::state_near's), what each measured
    /// moment of each body (Body::measure, the bodies in order), alone at
    /// the far side of its proven bounds, changes the equations of motion
    /// by along each coordinate and each joint released for springs (n
    /// rows, Moment::f's then fr's): the effort the body needs there at the
    /// state's accelerations, and its inertia's along the released joints
    /// accelerating at `released` (no state of the model reads them: their
    /// inertia enters by M). Each is linear in the moments, so each moment
    /// moves them by its own change, whatever the others' (as Moment::mass).
    /// Infinite along what a moment not measured moves. After a spinning
    /// body's moments, its departure from round within what its measure is
    /// known to (Body::unknown), which may lie at any angle across its
    /// spin's axis: its centre that far off the axis along each of two ways
    /// across it, and its second moment not round about it along each of
    /// four (any angle's change a combination of theirs, each at most 1,
    /// so within the sum of their sizes). After the bodies, each carried
    /// body's mass per mm at the far side of what it is known to
    /// (Strip::err). Err: the state is not found (its loop does not close).
    pub(crate) fn mass_sources(&self, t: f64, q: &[[f64; 3]], released: &[f64]) -> Result<Vec<Vec<f64>>, String> {
        let ctx = self.near_ctx(t, q)?;
        let nq = self.coords.len();
        let n = nq + released.len();
        let gv = self.gravity.map(|x| 1000.0 * x);
        let frames = self.spring_frames(&ctx);
        let mut out = Vec::new();
        for b in self.bodies.iter().filter(|b| !b.measure.is_empty() || b.unknown != [0.0; 2]) {
            let Placed { r, p, dq, dr, dt2 } = self.placed(b, &ctx, &frames);
            let dq: Vec<(usize, &[f64])> = dq.iter().map(|(k, d)| (*k, &d[..])).collect();
            let posed = Posed::new(&r, &dr, &p, &dt2, &dq);
            // What it needs along each row, and that sum's rounding (ROUNDING of
            // the magnitudes summed): a change within 16 of it, as the
            // accelerations' own residual is held to (simulate::forces), is none
            // (a body falling free, m (a - g), moves by none of its mass).
            let alone = |mass: f64, c: V3, inertia: &M3| -> Vec<(f64, f64)> {
                let (mut f, mut scale, mut gm) = (vec![0.0; n], vec![0.0; n], vec![vec![0.0; n]; n]);
                let mut wr = [Wrench { force: [0.0; 3], moment: [0.0; 3], frame: [0.0; 12], noise: (0.0, 0.0) }];
                Sums { f: &mut f, scale: &mut scale, gm: &mut gm, heft: None, wr: &mut wr, gv }.rigid(&posed, mass, c, inertia, &[0]);
                (0..n)
                    .map(|k| {
                        let terms = released.iter().enumerate().map(|(j, a)| 1e-6 * gm[k][nq + j] * a);
                        (f[k] + terms.clone().sum::<f64>(), ROUNDING * (scale[k] + terms.map(f64::abs).sum::<f64>()))
                    })
                    .collect()
            };
            let base = alone(b.mass, b.center, &b.inertia);
            let change =
                |moved: Vec<(f64, f64)>| -> Vec<f64> { moved.iter().zip(&base).map(|(x, y)| if (x.0 - y.0).abs() <= 16.0 * x.1.max(y.1) { 0.0 } else { x.0 - y.0 }).collect() };
            for (dm, dc, di) in &b.measure {
                if !dm.is_finite() {
                    out.push((0..n).map(|k| if dq.iter().any(|d| d.0 == k) { f64::INFINITY } else { 0.0 }).collect());
                    continue;
                }
                out.push(change(alone(b.mass + dm, std::array::from_fn(|k| b.center[k] + dc[k]), &std::array::from_fn(|i| std::array::from_fn(|j| b.inertia[i][j] + di[i][j])))));
            }
            if let Some((_, (_, u))) = b.spin.filter(|_| b.unknown != [0.0; 2] && b.mass > 0.0) {
                let u = u.map(|x| x / norm(u));
                let least = (0..3).min_by(|&i, &j| u[i].abs().total_cmp(&u[j].abs())).expect("an axis");
                let a = cross(u, std::array::from_fn(|k| f64::from(u8::from(k == least))));
                let a = a.map(|x| x / norm(a));
                let w = cross(u, a);
                let [first, second] = b.unknown;
                let pair = |x: V3, y: V3| -> M3 { std::array::from_fn(|i| std::array::from_fn(|j| second * (x[i] * y[j] + y[i] * x[j]))) };
                let off = |e: V3| -> V3 { std::array::from_fn(|k| b.center[k] + first / b.mass * e[k]) };
                let plus = |d: M3| -> M3 { std::array::from_fn(|i| std::array::from_fn(|j| b.inertia[i][j] + d[i][j])) };
                let flat: M3 = std::array::from_fn(|i| std::array::from_fn(|j| second * (a[i] * a[j] - w[i] * w[j])));
                let departed = [alone(b.mass, off(a), &b.inertia), alone(b.mass, off(w), &b.inertia)]
                    .into_iter()
                    .chain([pair(a, u), pair(w, u), flat, pair(a, w)].map(|d| alone(b.mass, b.center, &plus(d))));
                out.extend(departed.map(change));
            }
        }
        // A carried body's mass per mm at the far side of what it is known to (Strip::err): its share at
        // that much, linear in it, within 16 of the share's own rounding none. Non-finite/negative error
        // fails closed to infinity rather than reading as exact.
        for s in self.strips.iter().filter(|s| s.err != 0.0) {
            let err = if s.err >= 0.0 { s.err } else { f64::INFINITY };
            let (f, scale, gm, _) = self.carried_alone(s, &ctx, n, gv);
            out.push(
                (0..n)
                    .map(|k| {
                        let terms = released.iter().enumerate().map(|(j, a)| 1e-6 * gm[k][nq + j] * a);
                        let (x, rounding) = (f[k] + terms.clone().sum::<f64>(), ROUNDING * (scale[k] + terms.map(f64::abs).sum::<f64>()));
                        if x.abs() <= 16.0 * rounding { 0.0 } else { by(err, x) }
                    })
                    .collect(),
            );
        }
        Ok(out)
    }

    /// A carried body's share at `mass` g per mm of its material: its
    /// material, path stations [shift, shift + length] (a loop's round its
    /// length, an open path's first and last pieces running on), on each
    /// piece the material covers, at x mm along the piece.
    fn carry(&self, s: &Strip, ctx: &Ctx, mass: f64, sums: &mut Sums) {
        let v3 = |n: &[Num; 3]| -> V3 { n.each_ref().map(|e| e.at(ctx)) };
        let (shift, length, total) = (s.shift.at(ctx), s.length.at(ctx), s.total.at(ctx));
        let n = s.pieces.len();
        let turns = if s.closed && total > 0.0 { (shift / total).floor() as i64..=((shift + length) / total).floor() as i64 } else { 0..=0 };
        for i in 0..n {
            let len = s.lengths[i].at(ctx);
            if len <= 0.0 {
                continue;
            }
            // Its nodes (x along it, and the mm of material each stands for).
            let (mut pts, mut ws) = (Vec::new(), Vec::new());
            for turn in turns.clone() {
                let start = s.starts[i].at(ctx) + turn as f64 * total;
                let lo = if !s.closed && i == 0 { f64::NEG_INFINITY } else { start };
                let hi = if !s.closed && i + 1 == n { f64::INFINITY } else { start + len };
                let (a, b) = (lo.max(shift), hi.min(shift + length));
                if !(b > a) {
                    continue;
                }
                // A line's two nodes; an arc's four per quarter turn.
                let (rule, parts): (&[(f64, f64)], usize) = match &s.sweeps[i] {
                    None => (&GAUSS2, 1),
                    Some(sw) => (&GAUSS4, ((b - a) * sw.at(ctx).abs() / len / std::f64::consts::FRAC_PI_2).ceil().max(1.0) as usize),
                };
                let h = (b - a) / parts as f64;
                for part in 0..parts {
                    let mid = a + h * (part as f64 + 0.5);
                    for &(node, weight) in rule {
                        pts.push([mid + 0.5 * h * node - start, 0.0, 0.0]);
                        ws.push(0.5 * h * weight);
                    }
                }
            }
            if pts.is_empty() {
                continue;
            }
            let p = &s.pieces[i];
            let mut col = Vec::new();
            // Each number at every node: by node, then number.
            let mut at = |n: &[Num; 3]| -> Vec<V3> {
                let mut out = vec![[0.0; 3]; pts.len()];
                for (k, e) in n.iter().enumerate() {
                    e.at_points(ctx, &pts, &mut col);
                    for (o, v) in out.iter_mut().zip(&col) {
                        o[k] = *v;
                    }
                }
                out
            };
            let (pos, acc) = (at(&p.pos), at(&p.acc));
            let dq: Vec<(usize, Vec<V3>)> = p.dq.iter().map(|(k, d)| (*k, at(d))).collect();
            for (j, w) in ws.iter().enumerate() {
                let dq: Vec<(usize, V3)> = dq.iter().map(|(k, d)| (*k, d[j])).collect();
                sums.point(w * mass, pos[j], acc[j], &dq, &s.of);
            }
            // Its material on the parts its span or wrap is on: a wrap's
            // where it is; a span's on its two ends, at them, by where it
            // is along it (a string's load shared between its ends), and
            // past an end (an open belt's tail at its clamp) on that end.
            if let Some(&rides) = s.on.get(i) {
                let nr = sums.wr.len();
                let on = |c: &Coupled, e: usize, x: f64| -> Vec<(usize, f64)> { (0..nr).filter(|&r| c.carries[r][e] && !s.of.contains(&r)).map(|r| (r, x)).collect() };
                for (j, w) in ws.iter().enumerate() {
                    match rides {
                        Rides::Span(c, l) => {
                            let c = &self.couplings[c];
                            let (q, u, _, len) = &c.lines[l];
                            let (q, u, len) = (v3(q), v3(u), len.at(ctx));
                            let n = norm(u);
                            let along = (0..3).map(|k| (pos[j][k] - q[k]) * u[k] / n).sum::<f64>() / len;
                            if (0.0..=1.0).contains(&along) {
                                let a_end: V3 = std::array::from_fn(|k| q[k] + u[k] / n * len);
                                sums.held(w * mass, a_end, acc[j], &on(c, 0, along));
                                sums.held(w * mass, q, acc[j], &on(c, 1, 1.0 - along));
                            } else {
                                sums.held(w * mass, pos[j], acc[j], &on(c, usize::from(along < 0.0), 1.0));
                            }
                        }
                        Rides::Wrap(c, _, e) => sums.held(w * mass, pos[j], acc[j], &on(&self.couplings[c], e, 1.0)),
                    }
                }
            }
        }
    }

    /// A carried body's share alone at 1 g per mm of its material (Model::carry), along `n` rows
    /// (the coordinates, then the joints released for springs): its effort, that sum's rounding's
    /// scale, its generalized inertia and each reaction's wrench. Linear in the mass.
    fn carried_alone(&self, s: &Strip, ctx: &Ctx, n: usize, gv: V3) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<Wrench>) {
        let (mut f, mut scale, mut gm) = (vec![0.0; n], vec![0.0; n], vec![vec![0.0; n]; n]);
        let mut wr: Vec<Wrench> = self.reactions.iter().map(|_| Wrench { force: [0.0; 3], moment: [0.0; 3], frame: [0.0; 12], noise: (0.0, 0.0) }).collect();
        self.carry(s, ctx, 1.0, &mut Sums { f: &mut f, scale: &mut scale, gm: &mut gm, heft: None, wr: &mut wr, gv });
        (f, scale, gm, wr)
    }

    /// The model at moment t as `at` gives it, not kept: a client reading
    /// each of many moments once (a viewer's frames) keeps none of them.
    pub fn at_once(&self, t: f64) -> Moment {
        self.compute(t, &vec![false; self.groups.len()])
    }

    /// The model at moment t, every body at its own mass properties.
    pub fn at(&self, t: f64) -> Arc<Moment> {
        self.moment(t, &vec![false; self.groups.len()])
    }

    /// The model at moment t, the groups `rest` holds read round (each
    /// moment computed once).
    pub(super) fn moment(&self, t: f64, rest: &[bool]) -> Arc<Moment> {
        let key = (t.to_bits(), rest.to_vec());
        if let Some(m) = self.moments.lock().unwrap().get(&key) {
            return m.clone();
        }
        let mut m = self.compute(t, rest);
        // What is not known of the groups' departures bounds each effort at this moment too.
        let unknown = self.unknown_moves(&Ctx::new(t, self.params.clone()));
        if unknown.iter().any(|b| *b != 0.0) && m.pm.iter().all(|p| p.len() == unknown.len()) {
            m.bound_f = m.bound_f.iter().zip(&unknown).map(|(x, y)| x + y).collect();
            m.bound = m.pm.iter().map(|p| p.iter().zip(&m.bound_f).map(|(p, d)| p.abs() * d).sum()).collect();
            m.bound_u = m.uncovered.iter().map(|q| q.iter().zip(&m.bound_f).map(|(q, d)| q * d).sum()).collect();
        }
        let m = Arc::new(m);
        self.moments.lock().unwrap().insert(key, m.clone());
        m
    }

    fn compute(&self, t: f64, rest: &[bool]) -> Moment {
        self.compute_in(t, self.ctx(t), rest)
    }

    /// The model at moment t in context `ctx` (its solves followed there).
    pub(super) fn compute_in(&self, t: f64, ctx: Result<Ctx, String>, rest: &[bool]) -> Moment {
        let nq = self.coords.len();
        let nr = self.springs.as_ref().map_or(0, |s| s.joints.len());
        let ctx = match ctx {
            Ok(ctx) => ctx,
            Err(why) => {
                let (na, nan) = (self.actuators.len(), f64::NAN);
                return Moment {
                    t,
                    q: vec![[nan; 3]; nq],
                    f: vec![nan; nq],
                    unpowered: vec![None; nq],
                    effort: vec![nan; na],
                    rate: vec![nan; na],
                    lo: vec![nan; na],
                    hi: vec![nan; na],
                    reflected: vec![nan; na],
                    noise_f: vec![0.0; nq],
                    noise: vec![0.0; na],
                    pm: vec![vec![0.0; nq]; na],
                    bound_f: vec![0.0; nq],
                    bound: vec![0.0; na],
                    mass: vec![0.0; na],
                    mass_w: vec![Vec::new(); self.reactions.len()],
                    mass_t: Vec::new(),
                    uncovered: vec![vec![0.0; nq]; nq],
                    bound_u: vec![0.0; nq],
                    reactions: self.reactions.iter().map(|_| Wrench { force: [nan; 3], moment: [nan; 3], frame: [nan; 12], noise: (0.0, 0.0) }).collect(),
                    tensions: Vec::new(),
                    tension_noise: Vec::new(),
                    frictions: Vec::new(),
                    coupling: None,
                    lossless: Vec::new(),
                    loss_bound: vec![0.0; na],
                    f_ideal: vec![nan; nq],
                    apart: vec![nan; nq],
                    worst: vec![0.0; na],
                    bm: vec![vec![nan; nq]; na],
                    sites: Vec::new(),
                    gm: vec![vec![nan; nq + nr]; nq + nr],
                    heft: vec![nan; nq + nr],
                    command: vec![None; na],
                    fr: vec![nan; nr],
                    noise_r: vec![0.0; nr],
                    springs: None,
                    unsolved: Some(why),
                };
            }
        };
        let gv = self.gravity.map(|x| 1000.0 * x);
        // Along the coordinates, then the joints released for springs (each
        // body they carry moved by its unit motion: springs.rs).
        let n = nq + nr;
        let frames = self.spring_frames(&ctx);
        let mut carried = vec![0.0; nr];
        let mut f = vec![0.0; n];
        // How large the products summed into each f are (|J| |F|: a dot
        // product of vectors at right angles rounds to that's 1e-16, not
        // to 0; and F's and N's own terms, m a and m g, I alpha and
        // omega x I omega, which cancel to rounding in a round body
        // spinning steadily): its rounding's scale.
        let mut scale = vec![0.0; n];
        // Generalized inertia, g mm^2 per unit^2 of the coordinates.
        let mut gm = vec![vec![0.0; n]; n];
        let mut heft = vec![0.0; n];
        let vals = |n: &[Num]| -> Vec<f64> { n.iter().map(|e| e.at(&ctx)).collect() };
        let v3 = |n: &[Num; 3]| -> V3 { n.each_ref().map(|e| e.at(&ctx)) };
        let mut reactions: Vec<Wrench> =
            self.reactions.iter().map(|r| Wrench { force: [0.0; 3], moment: [0.0; 3], frame: std::array::from_fn(|k| r.frame[k].at(&ctx)), noise: (0.0, 0.0) }).collect();
        let mut sums = Sums { f: &mut f, scale: &mut scale, gm: &mut gm, heft: Some(&mut heft), wr: &mut reactions, gv };
        // What the unseparated bodies need along each coordinate.
        let mut apart = vec![0.0; nq];
        let mut add_apart = |before: &[f64], after: &[f64]| {
            for k in 0..nq {
                apart[k] += after[k] - before[k];
            }
        };
        // Each coordinate's effort's reach for what the bodies' measures are known to (Body::measure).
        let mut mass_f = vec![0.0; n];
        // Each body's measured moment's change, as the reactions carrying it take it.
        let mut changes: Vec<Vec<(usize, V3, V3)>> = Vec::new();
        let (mut base_f, mut moved_f) = (vec![0.0; n], vec![0.0; n]);
        for b in &self.bodies {
            let before = if b.apart { sums.f.to_vec() } else { Vec::new() };
            // A member of a group read round: its mass averaged about the
            // axis, as a round body's is.
            let round = |mass: f64, c: V3, inertia: &M3| match b.spin {
                Some((g, axis)) if rest[g] => about_axis(mass, c, inertia, axis),
                _ => (c, *inertia),
            };
            let (c, inertia) = round(b.mass, b.center, &b.inertia);
            let Placed { r, p: pp, dq, dr, dt2 } = self.placed(b, &ctx, &frames);
            for &i in &b.released {
                carried[i] += b.mass;
            }
            let dq: Vec<(usize, &[f64])> = dq.iter().map(|(k, d)| (*k, &d[..])).collect();
            let posed = Posed::new(&r, &dr, &pp, &dt2, &dq);
            sums.rigid(&posed, b.mass, c, &inertia, &b.of);
            if !b.measure.is_empty() {
                base_f.fill(0.0);
                let base = posed.alone(b.mass, c, &inertia, gv, &mut base_f);
                for (dm, dc, di) in &b.measure {
                    // (A moment not measured: nothing it moves is known.)
                    if !dm.is_finite() {
                        dq.iter().for_each(|(k, _)| mass_f[*k] = f64::INFINITY);
                        changes.push(b.of.iter().map(|&r| (r, [f64::INFINITY; 3], [f64::INFINITY; 3])).collect());
                        continue;
                    }
                    let (c, inertia) =
                        round(b.mass + dm, std::array::from_fn(|k| b.center[k] + dc[k]), &std::array::from_fn(|i| std::array::from_fn(|j| b.inertia[i][j] + di[i][j])));
                    moved_f.fill(0.0);
                    let (force, moment) = posed.alone(b.mass + dm, c, &inertia, gv, &mut moved_f);
                    for (x, (y, z)) in mass_f.iter_mut().zip(moved_f.iter().zip(&base_f)) {
                        *x += (y - z).abs();
                    }
                    let gap = |a: V3, b: V3| -> V3 { std::array::from_fn(|k| a[k] - b[k]) };
                    changes.push(b.of.iter().map(|&r| (r, gap(force, base.0), gap(moment, base.1))).collect());
                }
            }
            if b.apart {
                add_apart(&before, sums.f);
            }
        }
        // A load: its force and torque, in the world's axes, at its frame's
        // origin x, which moves with what it acts on (placed R p + T: its
        // partial velocity dR R^T (x - T) + dT).
        for l in &self.loads {
            let x: V3 = std::array::from_fn(|k| l.at[k].at(&ctx));
            let (force, torque) = (v3(&l.force), v3(&l.torque));
            let on = vals(&l.on);
            let (ro, o) = split(&on);
            let arm: V3 = std::array::from_fn(|k| x[k] - o[k]);
            let mut partials: Vec<(usize, Vec<f64>)> = l.dq.iter().map(|(k, d)| (*k, vals(d))).collect();
            partials.extend(l.released.iter().map(|&i| (nq + i, self.unit_motion(&frames, i, &on[..9], &on[9..]).to_vec())));
            let dq: Vec<(usize, V3, V3)> = partials
                .iter()
                .map(|(k, d)| {
                    let (dr, dp) = split(d);
                    let w = mul(&dr, &tr(&ro));
                    (*k, std::array::from_fn(|i| mv(&w, &arm)[i] + dp[i]), vee(&w))
                })
                .collect();
            let before = if l.apart { sums.f.to_vec() } else { Vec::new() };
            sums.load(force, torque, x, &dq, &l.of);
            if l.apart {
                add_apart(&before, sums.f);
            }
        }
        // A carried body: its material at its mass per mm (Model::carry), and what that is known to
        // (Strip::err) moving it by its share at that much, linear in it.
        for s in &self.strips {
            let before = if s.apart { sums.f.to_vec() } else { Vec::new() };
            self.carry(s, &ctx, s.mass, &mut sums);
            if s.err != 0.0 {
                let err = if s.err >= 0.0 { s.err } else { f64::INFINITY };
                let (f, _, _, wr) = self.carried_alone(s, &ctx, n, gv);
                for (x, y) in mass_f.iter_mut().zip(&f) {
                    *x += by(err, *y).abs();
                }
                changes.push(
                    wr.iter()
                        .enumerate()
                        .filter(|(_, w)| w.force != [0.0; 3] || w.moment != [0.0; 3])
                        .map(|(r, w)| (r, w.force.map(|x| by(err, x)), w.moment.map(|x| by(err, x))))
                        .collect(),
                );
            }
            if s.apart {
                add_apart(&before, sums.f);
            }
        }
        let fr = f.split_off(nq);
        let noise_r: Vec<f64> = scale.split_off(nq).iter().map(|s| ROUNDING * s).collect();
        let springs = self.springs.as_ref().map(|_| Arc::new(self.springs_at(&ctx, &frames, carried)));
        let f_ideal = f.clone();
        let mut frictions = Vec::new();
        let mut tension_noise = Vec::new();
        let mut mass_t = Vec::new();
        let (tensions, coupling, lossless, sites) = self.couple(&ctx, &mut reactions, &mut f, &mut scale, &mut frictions, &mut tension_noise, &mut changes, &mut mass_t);
        let coupling = coupling.or_else(|| self.belts_off(&ctx));
        // B: each actuator's drive along each coordinate, in radians for a
        // turn. The actuators' efforts tau solve B^T tau = f by least
        // squares, of least sum of (tau / span)^2, each actuator's span its
        // limits' width there: unit-free, and actuators alike share alike.
        // With C = W B (W the spans), tau = W C (C^T C)^+ f = P f.
        let na = self.actuators.len();
        let at = |n: &Num| n.at(&ctx);
        let (lo, hi): (Vec<f64>, Vec<f64>) = self.actuators.iter().map(|a| (at(&a.lo), at(&a.hi))).unzip();
        let span: Vec<f64> = (0..na).map(|a| 0.5 * hi[a] - 0.5 * lo[a]).collect();
        let bm: Vec<Vec<f64>> = self
            .actuators
            .iter()
            .map(|a| {
                let s = if a.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
                a.along.iter().map(|x| s * x.at(&ctx)).collect()
            })
            .collect();
        let fixed: Vec<_> = (0..na).map(|a| if lo[a] == hi[a] { lo[a] } else { 0.0 }).collect();
        let (pm, offset) = effort_split(&bm, &span, &fixed, nq);
        // An actuator moving no coordinate (held at a pose a shape follows)
        // has no effort found. One moving simulated coordinates alone
        // supplies its command, clamped to its limits, or none (simulate.rs):
        // the motion is what that effort makes, so what the motion needs
        // differs from it by the simulation's error alone (unpowered, below).
        let simulated = |a: usize| {
            let moved: Vec<usize> = (0..nq).filter(|&k| bm[a][k] != 0.0).collect();
            !moved.is_empty() && moved.iter().all(|k| self.simulated.iter().any(|s| s.0 == *k))
        };
        let effort: Vec<f64> = (0..na)
            .map(|a| match &self.actuators[a] {
                act if act.coupled.is_empty() => f64::NAN,
                act if simulated(a) => act.command.as_ref().map_or(0.0, |c| at(c).max(lo[a]).min(hi[a])),
                _ => offset[a] + dot_n(&pm[a], &f),
            })
            .collect();
        let noise_f: Vec<f64> = scale.iter().map(|s| ROUNDING * s).collect();
        let noise: Vec<f64> = (0..na).map(|a| (0..nq).map(|k| pm[a][k].abs() * noise_f[k]).sum()).collect();
        // Along a simulated coordinate the motion is the simulation's, whose
        // every state its forces balance (simulate.rs): what no actuator
        // supplies there is how far its interpolation between knots is from
        // them, no unpowered effort. The residual is read to its own
        // rounding: the effort's along k, and that of the efforts the
        // actuators supply there (each one's from every coordinate it
        // balances, and their sum's): drives that cancel along k, as a
        // steer motor cancels a drive motor geared to its heading, leave
        // that alone where k carries no load.
        let unpowered = (0..nq)
            .map(|k| {
                let (supplied, rounding) = (0..na)
                    .filter(|&a| bm[a][k] != 0.0)
                    .map(|a| (bm[a][k] * effort[a], bm[a][k].abs() * (noise[a] + ROUNDING * effort[a].abs())))
                    .fold((0.0, noise_f[k]), |(s, n), (x, e)| (s + x, n + e));
                let r = f[k] - supplied;
                (self.needs_drive(k) && !(r.abs() <= rounding)).then_some(r)
            })
            .collect();
        let uncovered = (0..nq)
            .map(|k| (0..nq).map(|l| (f64::from(u8::from(k == l)) - (0..na).filter(|&a| bm[a][k] != 0.0).map(|a| bm[a][k] * pm[a][l]).sum::<f64>()).abs()).collect())
            .collect();
        // The inertia reflected to an actuator: the effort it needs per unit
        // of its own acceleration, the coordinates accelerating by the least
        // that gives it (B^+ = (B^T B)^+ B^T).
        let unit_shares = effort_split(&bm, &vec![1.0; na], &vec![0.0; na], nq).0;
        let reflected = (0..na).map(|a| (0..nq).map(|k| (0..nq).map(|l| pm[a][k] * gm[k][l] * unit_shares[a][l]).sum::<f64>()).sum()).collect();
        // The most the loss charged to each actuator can miss: the share of
        // the power its lossy couplings would pass for what the unseparated
        // bodies need along its coordinates (their pull through a coupling
        // the balances do not see), its worst efficiency's.
        // Stages in series compound: the worst of each stage moving its
        // coordinates, multiplied.
        let worst: Vec<f64> = (0..na)
            .map(|a| {
                let act = &self.actuators[a];
                let mut stages: Vec<(usize, f64)> = Vec::new();
                for s in sites.iter().filter(|s| act.coupled.iter().any(|&k| s.weights[k] != 0.0)) {
                    let w = s.worst;
                    match stages.iter_mut().find(|x| x.0 == s.stage) {
                        Some(x) => x.1 = x.1.max(w),
                        None => stages.push((s.stage, w)),
                    }
                }
                stages.iter().fold(1.0, |p, x| p * (1.0 + x.1)) - 1.0
            })
            .collect();
        let loss_bound = (0..na).map(|a| worst[a] * unfed_need(&pm[a], &apart, &sites, None)).collect();
        // Each reaction's changes, and where couplings lose power, the reactions and tensions widened
        // by their worst efficiency's share, as an effort's is.
        let lossy = 1.0 + sites.iter().map(|s| s.worst).fold(0.0, f64::max);
        let mut mass_w: Vec<Vec<(V3, V3)>> = vec![Vec::new(); self.reactions.len()];
        for ch in &changes {
            let mut each: Vec<(usize, V3, V3)> = Vec::new();
            for &(r, df, dm) in ch {
                match each.iter_mut().find(|e| e.0 == r) {
                    Some(e) => (0..3).for_each(|k| (e.1[k], e.2[k]) = (e.1[k] + df[k], e.2[k] + dm[k])),
                    None => each.push((r, df, dm)),
                }
            }
            for (r, df, dm) in each {
                mass_w[r].push((df.map(|x| lossy * x), dm.map(|x| lossy * x)));
            }
        }
        let mass_t: Vec<f64> = mass_t.iter().map(|x| lossy * x).collect();
        // (An effort that is its command rests on no mass.)
        let mass = (0..na)
            .map(|a| match simulated(a) {
                true => 0.0,
                false => (1.0 + worst[a]) * (0..nq).filter(|&k| pm[a][k] != 0.0).map(|k| pm[a][k].abs() * mass_f[k]).sum::<f64>(),
            })
            .collect();
        Moment {
            t,
            q: self.coords.iter().map(|c| [at(&c.q), at(&c.dq), at(&c.ddq)]).collect(),
            f,
            unpowered,
            rate: self.actuators.iter().map(|a| at(&a.rate)).collect(),
            lo,
            hi,
            effort,
            reflected,
            noise_f,
            noise,
            pm,
            bound_f: vec![0.0; nq],
            bound: vec![0.0; na],
            mass,
            mass_w,
            mass_t,
            loss_bound,
            uncovered,
            bound_u: vec![0.0; nq],
            reactions: reactions.into_iter().map(|w| Wrench { noise: (ROUNDING * w.noise.0, ROUNDING * w.noise.1), ..w }).collect(),
            tensions,
            tension_noise,
            frictions,
            coupling,
            lossless,
            f_ideal,
            apart,
            worst,
            bm,
            sites,
            command: self.actuators.iter().map(|a| a.command.as_ref().map(at)).collect(),
            gm,
            heft,
            fr,
            noise_r,
            springs,
            unsolved: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{M3, Posed, Sums, V3, Wrench};

    /// What a measure changes is read off the body's share alone
    /// (Posed::alone): the same numbers, to the bit, as Sums::rigid sums
    /// for that body by itself, over placements of every turn, spin,
    /// offset and size.
    #[test]
    fn a_body_alone_is_its_rigid_share_to_the_bit() {
        let mut seed = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = move |scale: f64| {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            scale * ((seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5)
        };
        for case in 0..500 {
            let size = [1e-3, 1.0, 1e3, 1e6][case % 4];
            let (u, v) = ([next(2.0), next(2.0), next(2.0)], next(6.3));
            let n = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
            let u = u.map(|x| x / n);
            let (s, c) = v.sin_cos();
            let r: Vec<f64> = (0..9)
                .map(|e| (e / 3, e % 3))
                .map(|(i, j)| u[i] * u[j] * (1.0 - c) + if i == j { c } else { [0.0, -u[2], u[1], u[2], 0.0, -u[0], -u[1], u[0], 0.0][3 * i + j] * s })
                .collect();
            let dr: Vec<f64> = (0..9).map(|_| next(10.0)).collect();
            let p: Vec<f64> = (0..3).map(|_| next(size)).collect();
            let dt2: Vec<f64> = (0..12).map(|_| next(1e4)).collect();
            let parts: Vec<(usize, Vec<f64>)> = (0..1 + case % 4).map(|k| (k, (0..12).map(|_| next(2.0)).collect())).collect();
            let dq: Vec<(usize, &[f64])> = parts.iter().map(|(k, d)| (*k, &d[..])).collect();
            let posed = Posed::new(&r, &dr, &p, &dt2, &dq);
            let mass = next(1e3).abs();
            let centre: V3 = std::array::from_fn(|_| next(size));
            let inertia: M3 = std::array::from_fn(|_| std::array::from_fn(|_| next(size * size)));
            let gv: V3 = [0.0, 0.0, -9810.0];
            let n = 4;
            let (mut f, mut scale, mut gm) = (vec![0.0; n], vec![0.0; n], vec![vec![0.0; n]; n]);
            let mut wr = [Wrench { force: [0.0; 3], moment: [0.0; 3], frame: [0.0; 12], noise: (0.0, 0.0) }];
            Sums { f: &mut f, scale: &mut scale, gm: &mut gm, heft: None, wr: &mut wr, gv }.rigid(&posed, mass, centre, &inertia, &[0]);
            let mut g = vec![0.0; n];
            let (force, moment) = posed.alone(mass, centre, &inertia, gv, &mut g);
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<u64>>();
            assert_eq!(bits(&f), bits(&g), "case {case}: effort");
            assert_eq!(bits(&wr[0].force), bits(&force), "case {case}: force");
            assert_eq!(bits(&wr[0].moment), bits(&moment), "case {case}: moment");
        }
    }
}
