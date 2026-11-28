use super::sampling::Tracker;
use super::{M3, Model, Moment, Num, Point, Strip, V3, mv, tr};
use ink_carried::{Carried, PathPiece};
use ink_core::eval::Ctx;
use ink_core::graph::{Graph, NodeId, Var};
use ink_core::v3::{cross, norm};
use ink_lang::Doc;
use ink_lang::view::ModelView;

/// The numbers a carried body's material moves by: its path's pieces,
/// stations and length, and its shift.
pub(super) fn path_numbers(c: &Carried) -> Vec<NodeId> {
    let mut ns: Vec<NodeId> = c.pieces.iter().flat_map(PathPiece::numbers).collect();
    ns.extend(c.starts.iter().chain(&c.lengths).chain([&c.total, &c.shift]));
    ns
}

/// A carried body placed by xf, its material of density rho: its profile
/// cut to its length measured (None: not proven bounded), spread evenly
/// along it on its path line, and on each piece of its path the point a
/// (x) mm along it. Material station m is a = m + shift - start along the
/// piece, so its derivatives at fixed m are the chain rule's: d/dt =
/// (d/da) a' + d/dt at fixed a, d/dq = (d/da) da/dq + d/dq at fixed a.
/// Tracked: its circulation, the shift's speed and acceleration (the
/// path's shape moves with the bodies it wraps).
pub(super) fn strip(
    v: &mut ModelView,
    _label: &str,
    c: &Carried,
    xf: &ink_lang::Xf,
    rho: NodeId,
    along: &[usize],
    eps: &[u32],
    of: Vec<usize>,
    tracker: &mut Tracker,
) -> Option<(Strip, Vec<[[NodeId; 3]; 2]>, ink_kernel::query::Measured)> {
    let length = v.doc.scalar(c.length, &v.ctx);
    let g = &mut v.doc.g;
    let x = g.var(Var::X);
    let l = g.constant(length);
    let (before, after) = (g.neg(x), g.sub(x, l));
    let ends = g.max(before, after);
    let cut = g.max(c.profile, ends);
    let m = v.inertia(cut, Some(rho), None)?;
    let d = &mut v.doc;
    let speed = d.derivative(c.shift, Var::T);
    let acc = d.derivative(speed, Var::T);
    tracker.add(d, &[speed], None, None, &[]);
    tracker.add(d, &[acc], None, None, &[]);
    let mut pieces = Vec::new();
    let mut ends = Vec::new();
    for i in 0..c.pieces.len() {
        let g = &mut d.g;
        let placed = |g: &mut Graph, p: [NodeId; 3]| -> [NodeId; 3] {
            std::array::from_fn(|k| {
                let terms: [NodeId; 3] = std::array::from_fn(|j| g.mul(xf.r[k][j], p[j]));
                let s = g.add(terms[0], terms[1]);
                let s = g.add(s, terms[2]);
                g.add(s, xf.t[k])
            })
        };
        let p = c.point(g, i, x);
        let w = placed(g, p);
        let zero = g.constant(0.0);
        let (p0, p1) = (c.point(g, i, zero), c.point(g, i, c.lengths[i]));
        ends.push([placed(g, p0), placed(g, p1)]);
        let rel = g.sub(c.shift, c.starts[i]);
        let rate = d.derivative(rel, Var::T);
        let material = |d: &mut Doc, e: NodeId| {
            let (ea, et) = (d.derivative(e, Var::X), d.derivative(e, Var::T));
            let p = d.g.mul(ea, rate);
            d.g.add(p, et)
        };
        let acc: [NodeId; 3] = w.map(|e| {
            let de = material(d, e);
            material(d, de)
        });
        let mut dq = Vec::new();
        for &k in along {
            let da = d.derivative(rel, Var::Param(eps[k]));
            let dk: [NodeId; 3] = w.map(|e| {
                let (ea, eq) = (d.derivative(e, Var::X), d.derivative(e, Var::Param(eps[k])));
                let p = d.g.mul(ea, da);
                d.g.add(p, eq)
            });
            if dk.iter().any(|&x| d.g.constant_value(x) != Some(0.0)) {
                dq.push((k, dk.map(|n| Num::of(&d.g, n))));
            }
        }
        pieces.push(Point { pos: w.map(|n| Num::of(&d.g, n)), acc: acc.map(|n| Num::of(&d.g, n)), dq });
    }
    let g = &d.g;
    let num = |n: NodeId| Num::of(g, n);
    let strip = Strip {
        label: _label.to_string(),
        mass: if length > 0.0 { m.sums[0] / length } else { 0.0 },
        err: 0.0,
        shift: num(c.shift),
        length: num(c.length),
        closed: c.closed,
        total: num(c.total),
        starts: c.starts.iter().map(|&n| num(n)).collect(),
        lengths: c.lengths.iter().map(|&n| num(n)).collect(),
        sweeps: c.pieces.iter().map(|p| if let PathPiece::Arc { sweep, .. } = p { Some(num(*sweep)) } else { None }).collect(),
        pieces,
        ends: ends.iter().map(|e| e.map(|p| p.map(num))).collect(),
        on: Vec::new(),
        of,
        apart: false,
    };
    Some((strip, ends, m))
}

impl Model {
    /// What reaction r's path carries, for one who adds up its balance apart from the solve: each rigid
    /// body's mass (g) and its centre of mass in the world at moment t (mm), by its instance path; and the
    /// carried instances that are no weighed rigid body here (a belt counted along its path, a body whose
    /// shape moves, one left out for having no mass), which such a sum cannot count.
    #[allow(clippy::type_complexity)]
    pub fn carried_masses(&self, r: usize, t: f64) -> Option<(Vec<(f64, V3)>, Vec<&str>)> {
        let ctx = self.ctx(t).ok()?;
        let carries = &self.reactions[r].carries;
        let mut bodies = Vec::new();
        for b in self.bodies.iter().filter(|b| carries.contains(&b.label)) {
            let (rot, p) = (b.at.r.each_ref().map(|n| n.at(&ctx)), b.at.p.each_ref().map(|n| n.at(&ctx)));
            let c: V3 = std::array::from_fn(|i| rot[3 * i] * b.center[0] + rot[3 * i + 1] * b.center[1] + rot[3 * i + 2] * b.center[2] + p[i]);
            bodies.push((b.mass, c));
        }
        let other = carries.iter().filter(|l| !self.bodies.iter().any(|b| b.label == **l)).chain(&self.reactions[r].belts).map(String::as_str).collect();
        Some((bodies, other))
    }

    /// The angular momentum about the world's origin of the rigid bodies reaction r's path carries, at moment
    /// t (N mm s): each one's mass at its centre's velocity, and its own inertia at its turning rate, both by
    /// differences of where it is and how it is turned at t - h and t + h. For one who adds up the joint's
    /// balance of moments apart from the solve.
    pub fn carried_momentum(&self, r: usize, t: f64, h: f64) -> Option<V3> {
        let carries = &self.reactions[r].carries;
        let (before, after) = (self.ctx(t - h).ok()?, self.ctx(t + h).ok()?);
        let mut sum = [0.0; 3];
        for b in self.bodies.iter().filter(|b| carries.contains(&b.label)) {
            let placed = |ctx: &Ctx| -> (M3, V3) {
                let (rot, p) = (b.at.r.each_ref().map(|n| n.at(ctx)), b.at.p.each_ref().map(|n| n.at(ctx)));
                let rot: M3 = std::array::from_fn(|i| std::array::from_fn(|j| rot[3 * i + j]));
                let c = mv(&rot, &b.center);
                (rot, std::array::from_fn(|i| c[i] + p[i]))
            };
            let ((r0, c0), (r1, c1)) = (placed(&before), placed(&after));
            let (c, v): (V3, V3) = (std::array::from_fn(|i| 0.5 * (c0[i] + c1[i])), std::array::from_fn(|i| (c1[i] - c0[i]) / (2.0 * h)));
            // Its turning rate: R1 R0' is the turn made over 2 h, whose skew part is the rate's, to its order.
            let turn: M3 = std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| r1[i][k] * r0[j][k]).sum()));
            let w: V3 = [(turn[2][1] - turn[1][2]) / (4.0 * h), (turn[0][2] - turn[2][0]) / (4.0 * h), (turn[1][0] - turn[0][1]) / (4.0 * h)];
            // Its inertia in the world at t: the mean of the two turnings', to the same order.
            let spin = |rot: &M3| mv(rot, &mv(&b.inertia, &mv(&tr(rot), &w)));
            let (s0, s1) = (spin(&r0), spin(&r1));
            let orbit = cross(c, v);
            for i in 0..3 {
                sum[i] += 1e-6 * (b.mass * orbit[i] + 0.5 * (s0[i] + s1[i]));
            }
        }
        Some(sum)
    }

    /// The moment about the world's origin of what acts on what reaction r's path carries at moment m besides
    /// its joint (N mm): each body's weight at its centre under `gravity` (m/s2), each load on it (its torque,
    /// and its force at its point), and the couplings' lines from outside (as `coupling_pull`).
    pub fn carried_moment(&self, r: usize, m: &Moment, gravity: V3) -> Option<V3> {
        let ctx = self.ctx(m.t).ok()?;
        let (bodies, _) = self.carried_masses(r, m.t)?;
        let mut sum = [0.0; 3];
        for (mass, c) in bodies {
            let w = cross(c, gravity.map(|g| mass * 1e-3 * g));
            (0..3).for_each(|i| sum[i] += w[i]);
        }
        for load in self.loads.iter().filter(|l| l.of.contains(&r)) {
            let (f, tq, at) = (load.force.each_ref().map(|n| n.at(&ctx)), load.torque.each_ref().map(|n| n.at(&ctx)), load.at.each_ref().map(|n| n.at(&ctx)));
            let lever = cross(at, f);
            (0..3).for_each(|i| sum[i] += tq[i] + lever[i]);
        }
        if self.couplings.iter().any(|c| c.side[r] != 0) {
            if m.coupling.is_some() || m.tensions.is_empty() {
                return None;
            }
            for c in self.couplings.iter().filter(|c| c.side[r] != 0) {
                for (i, (p, u, pitch, _)) in c.lines.iter().enumerate() {
                    let (p, u) = (p.each_ref().map(|e| e.at(&ctx)), u.each_ref().map(|e| e.at(&ctx)));
                    let t = m.tensions[c.first + i] / norm(u);
                    let f = u.map(|x| f64::from(c.side[r]) * t * x);
                    let lever = cross(p, f);
                    (0..3).for_each(|k| sum[k] += lever[k] + pitch.at(&ctx) * f[k]);
                }
            }
        }
        for fr in m.frictions.iter().filter(|fr| self.bears(fr, r)) {
            (0..3).for_each(|k| sum[k] += fr.moment[k]);
        }
        Some(sum)
    }

    /// Body `label`'s mass properties set to these (mass g, centre mm and
    /// inertia about it g mm^2, in its frame): references stated with
    /// exact masses are checked at them. Another crate's tests turn it on
    /// with the feature `exact-masses` (Cargo.toml).
    #[cfg(any(test, feature = "exact-masses"))]
    pub fn set_mass(&mut self, label: &str, mass: f64, center: V3, inertia: M3) {
        let b = self.bodies.iter_mut().find(|b| b.label == label).unwrap_or_else(|| panic!("no moving body {label}"));
        (b.mass, b.center, b.inertia) = (mass, center, inertia);
        self.moments.lock().unwrap().clear();
        self.chosen.lock().unwrap().clear();
        self.swept.lock().unwrap().clear();
    }

    /// How many measured moments its bodies have (Model::mass_sources'
    /// sources).
    #[cfg(test)]
    pub(crate) fn mass_source_count(&self) -> usize {
        self.bodies.iter().map(|b| b.measure.len()).sum()
    }

    /// Source j (Model::mass_sources' order) moved `by` times the far side
    /// of its bounds: the body's mass, centre and inertia changed by that
    /// much of its change, the motion a rerun reads against the band.
    #[cfg(test)]
    pub(crate) fn move_source(&mut self, j: usize, by: f64) {
        let (b, k) = self.bodies.iter().flat_map(|b| (0..b.measure.len()).map(move |k| (b.label.clone(), k))).nth(j).expect("a source");
        let b = self.bodies.iter_mut().find(|x| x.label == b).expect("its body");
        let (dm, dc, di) = b.measure[k];
        b.mass += by * dm;
        (0..3).for_each(|i| b.center[i] += by * dc[i]);
        (0..3).for_each(|i| (0..3).for_each(|l| b.inertia[i][l] += by * di[i][l]));
        self.moments.lock().unwrap().clear();
        self.swept.lock().unwrap().clear();
    }

    /// Carried body `label`'s mass per mm (g), what that is known to (g),
    /// and its length at the model's moment (mm).
    #[cfg(test)]
    pub(crate) fn strip(&self, label: &str) -> (f64, f64, f64) {
        let s = self.strips.iter().find(|s| s.label == label).unwrap_or_else(|| panic!("no carried body {label}"));
        (s.mass, s.err, s.length.at(&Ctx::new(0.0, self.params.clone())))
    }

    /// Carried body `label`'s mass per mm set to this (g).
    #[cfg(test)]
    pub(crate) fn set_strip_mass(&mut self, label: &str, mass: f64) {
        self.strips.iter_mut().find(|s| s.label == label).unwrap_or_else(|| panic!("no carried body {label}")).mass = mass;
        self.moments.lock().unwrap().clear();
    }
}
