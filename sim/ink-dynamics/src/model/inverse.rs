use super::couplings::unfed_need;
use super::drives::effort_split;
use super::spin::about_axis;
use super::{Coupled, M3, Model, Moment, Num, ROUNDING, Rides, V3, Wrench, dot_n, mul, mv, rot, split, tr, vee};
use ink_core::eval::Ctx;
use ink_core::v3::{cross, dot, norm};
use std::sync::Arc;

/// What the bodies sum into at a moment: the effort along each coordinate
/// (N mm per unit of it), its rounding's scale, and the generalized
/// inertia (g mm^2 per unit^2); each reaction's wrench (Wrench: force N,
/// moment N mm about the world's origin) and the scales of its rounding;
/// gravity in mm/s^2.
struct Sums<'a> {
    f: &'a mut [f64],
    scale: &'a mut [f64],
    gm: &'a mut [Vec<f64>],
    wr: &'a mut [Wrench],
    gv: V3,
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
    /// its frame, at the moment its placement's rotation, the rotation's
    /// rate and its translation, its numbers' second t-derivatives, their
    /// partials along the coordinates that move it, and the reactions it is
    /// of.
    fn rigid(&mut self, mass: f64, c: V3, inertia: &M3, r: &[f64], dr: &[f64], p: &[f64], dt2: &[f64], dq: &[(usize, &[f64])], of: &[usize]) {
        let (r, rd) = (rot(r), rot(dr));
        let (rdd, pdd) = split(dt2);
        let rt = tr(&r);
        let x: V3 = std::array::from_fn(|k| mv(&r, &c)[k] + p[k]);
        let a: V3 = std::array::from_fn(|k| mv(&rdd, &c)[k] + pdd[k]);
        let w = vee(&mul(&rd, &rt));
        let al = {
            let (x, y) = (mul(&rdd, &rt), mul(&rd, &tr(&rd)));
            vee(&std::array::from_fn(|i| std::array::from_fn(|j| x[i][j] + y[i][j])))
        };
        let iw = mul(&mul(&r, inertia), &rt);
        // N, and N mm: g mm/s^2 and g mm^2/s^2 are 1e-6 of them.
        let force: V3 = std::array::from_fn(|k| 1e-6 * mass * (a[k] - self.gv[k]));
        let spin = cross(w, mv(&iw, &w));
        let torque: V3 = std::array::from_fn(|k| 1e-6 * (mv(&iw, &al)[k] + spin[k]));
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
        let partial: Vec<(usize, V3, V3)> = dq
            .iter()
            .map(|(k, d)| {
                let (dr, dp) = split(d);
                let jv: V3 = std::array::from_fn(|i| mv(&dr, &c)[i] + dp[i]);
                (*k, jv, vee(&mul(&dr, &rt)))
            })
            .collect();
        for &(k, jv, jw) in &partial {
            self.f[k] += dot(jv, force) + dot(jw, torque);
            self.scale[k] += projected_size(jv, force_terms) + projected_size(jw, torque_terms);
            for &(l, kv, kw) in &partial {
                self.gm[k][l] += mass * dot(jv, kv) + dot(jw, mv(&iw, &kw));
            }
        }
        self.react(of, x, force, torque, size);
    }
}

/// Gauss-Legendre nodes on -1..1 and weights: 2 integrate a line's
/// quadratic exactly (a point's acceleration and partials are affine along
/// it); 4 an arc's trig terms of degree 2 over a quarter turn to 5e-6 of the
/// arc's term.
const GAUSS2: [(f64, f64); 2] = [(-0.5773502691896257, 1.0), (0.5773502691896257, 1.0)];
const GAUSS4: [(f64, f64); 4] =
    [(-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461), (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538)];

impl Model {
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
        let (mut heft, mut carried) = (vec![0.0; nr], vec![0.0; nr]);
        let mut f = vec![0.0; n];
        // How large the products summed into each f are (|J| |F|: a dot
        // product of vectors at right angles rounds to that's 1e-16, not
        // to 0; and F's and N's own terms, m a and m g, I alpha and
        // omega x I omega, which cancel to rounding in a round body
        // spinning steadily): its rounding's scale.
        let mut scale = vec![0.0; n];
        // Generalized inertia, g mm^2 per unit^2 of the coordinates.
        let mut gm = vec![vec![0.0; n]; n];
        let vals = |n: &[Num]| -> Vec<f64> { n.iter().map(|e| e.at(&ctx)).collect() };
        let v3 = |n: &[Num; 3]| -> V3 { n.each_ref().map(|e| e.at(&ctx)) };
        let mut reactions: Vec<Wrench> =
            self.reactions.iter().map(|r| Wrench { force: [0.0; 3], moment: [0.0; 3], frame: std::array::from_fn(|k| r.frame[k].at(&ctx)), noise: (0.0, 0.0) }).collect();
        let mut sums = Sums { f: &mut f, scale: &mut scale, gm: &mut gm, wr: &mut reactions, gv };
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
        let mut alone = (vec![0.0; n], vec![0.0; n], vec![vec![0.0; n]; n], [Wrench { force: [0.0; 3], moment: [0.0; 3], frame: [0.0; 12], noise: (0.0, 0.0) }]);
        for b in &self.bodies {
            let before = if b.apart { sums.f.to_vec() } else { Vec::new() };
            // A member of a group read round: its mass averaged about the
            // axis, as a round body's is.
            let round = |mass: f64, c: V3, inertia: &M3| match b.spin {
                Some((g, axis)) if rest[g] => about_axis(mass, c, inertia, axis),
                _ => (c, *inertia),
            };
            let (c, inertia) = round(b.mass, b.center, &b.inertia);
            let p = &b.at;
            let (r, pp) = (vals(&p.r), vals(&p.p));
            let mut dq: Vec<(usize, Vec<f64>)> = p.dq.iter().map(|(k, d)| (*k, vals(d))).collect();
            for &i in &b.released {
                dq.push((nq + i, self.unit_motion(&frames, i, &r, &pp).to_vec()));
                carried[i] += b.mass;
                let (_, o, turn) = frames[i];
                heft[i] += if turn {
                    let x: V3 = std::array::from_fn(|k| mv(&rot(&r), &c)[k] + pp[k] - o[k]);
                    b.mass * dot(x, x) + inertia[0][0] + inertia[1][1] + inertia[2][2]
                } else {
                    b.mass
                };
            }
            let dq: Vec<(usize, &[f64])> = dq.iter().map(|(k, d)| (*k, &d[..])).collect();
            let (dr, dt2) = (vals(&p.dr), vals(&p.dt2));
            sums.rigid(b.mass, c, &inertia, &r, &dr, &pp, &dt2, &dq, &b.of);
            if !b.measure.is_empty() {
                let mut along = |mass: f64, c: V3, inertia: &M3| -> (Vec<f64>, V3, V3) {
                    alone.0.fill(0.0);
                    alone.3[0].force = [0.0; 3];
                    alone.3[0].moment = [0.0; 3];
                    let mut s = Sums { f: &mut alone.0, scale: &mut alone.1, gm: &mut alone.2, wr: &mut alone.3, gv };
                    s.rigid(mass, c, inertia, &r, &dr, &pp, &dt2, &dq, &[0]);
                    (alone.0.clone(), alone.3[0].force, alone.3[0].moment)
                };
                let base = along(b.mass, c, &inertia);
                for (dm, dc, di) in &b.measure {
                    // (A moment not measured: nothing it moves is known.)
                    if !dm.is_finite() {
                        dq.iter().for_each(|(k, _)| mass_f[*k] = f64::INFINITY);
                        changes.push(b.of.iter().map(|&r| (r, [f64::INFINITY; 3], [f64::INFINITY; 3])).collect());
                        continue;
                    }
                    let (c, inertia) =
                        round(b.mass + dm, std::array::from_fn(|k| b.center[k] + dc[k]), &std::array::from_fn(|i| std::array::from_fn(|j| b.inertia[i][j] + di[i][j])));
                    let (f, force, moment) = along(b.mass + dm, c, &inertia);
                    for (x, (y, z)) in mass_f.iter_mut().zip(f.iter().zip(&base.0)) {
                        *x += (y - z).abs();
                    }
                    let gap = |a: V3, b: V3| -> V3 { std::array::from_fn(|k| a[k] - b[k]) };
                    changes.push(b.of.iter().map(|&r| (r, gap(force, base.1), gap(moment, base.2))).collect());
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
        // A carried body: its material, path stations
        // [shift, shift + length] (a loop's round its length, an open path's
        // first and last pieces running on), on each piece the material
        // covers, at x mm along the piece.
        for s in &self.strips {
            let before = if s.apart { sums.f.to_vec() } else { Vec::new() };
            let (shift, length, total) = (s.shift.at(&ctx), s.length.at(&ctx), s.total.at(&ctx));
            let n = s.pieces.len();
            let turns = if s.closed && total > 0.0 { (shift / total).floor() as i64..=((shift + length) / total).floor() as i64 } else { 0..=0 };
            for i in 0..n {
                let len = s.lengths[i].at(&ctx);
                if len <= 0.0 {
                    continue;
                }
                // Its nodes (x along it, and the mm of material each stands for).
                let (mut pts, mut ws) = (Vec::new(), Vec::new());
                for turn in turns.clone() {
                    let start = s.starts[i].at(&ctx) + turn as f64 * total;
                    let lo = if !s.closed && i == 0 { f64::NEG_INFINITY } else { start };
                    let hi = if !s.closed && i + 1 == n { f64::INFINITY } else { start + len };
                    let (a, b) = (lo.max(shift), hi.min(shift + length));
                    if !(b > a) {
                        continue;
                    }
                    // A line's two nodes; an arc's four per quarter turn.
                    let (rule, parts): (&[(f64, f64)], usize) = match &s.sweeps[i] {
                        None => (&GAUSS2, 1),
                        Some(sw) => (&GAUSS4, ((b - a) * sw.at(&ctx).abs() / len / std::f64::consts::FRAC_PI_2).ceil().max(1.0) as usize),
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
                        e.at_points(&ctx, &pts, &mut col);
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
                    sums.point(w * s.mass, pos[j], acc[j], &dq, &s.of);
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
                                let (q, u, len) = (v3(q), v3(u), len.at(&ctx));
                                let n = norm(u);
                                let along = (0..3).map(|k| (pos[j][k] - q[k]) * u[k] / n).sum::<f64>() / len;
                                if (0.0..=1.0).contains(&along) {
                                    let a_end: V3 = std::array::from_fn(|k| q[k] + u[k] / n * len);
                                    sums.held(w * s.mass, a_end, acc[j], &on(c, 0, along));
                                    sums.held(w * s.mass, q, acc[j], &on(c, 1, 1.0 - along));
                                } else {
                                    sums.held(w * s.mass, pos[j], acc[j], &on(c, usize::from(along < 0.0), 1.0));
                                }
                            }
                            Rides::Wrap(c, _, e) => sums.held(w * s.mass, pos[j], acc[j], &on(&self.couplings[c], e, 1.0)),
                        }
                    }
                }
            }
            if s.apart {
                add_apart(&before, sums.f);
            }
        }
        let fr = f.split_off(nq);
        let noise_r: Vec<f64> = scale.split_off(nq).iter().map(|s| ROUNDING * s).collect();
        let springs = self.springs.as_ref().map(|_| Arc::new(self.springs_at(&ctx, &frames, heft, carried)));
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
            fr,
            noise_r,
            springs,
            unsolved: None,
        }
    }
}
