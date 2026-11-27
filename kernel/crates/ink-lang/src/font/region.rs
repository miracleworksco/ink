//! A glyph's filled region as simple loops. A font fills where its
//! contours wind round a point (nonzero), and a variable font's contours
//! overlap and cross themselves (strokes drawn over each other, a notch
//! closed by a loop), which `polygon` refuses. So: split every edge where
//! another crosses or touches it, keep the pieces with filled on one side
//! and empty on the other (judged by the winding number just off each),
//! turned so that filled is on their left, and join them into loops,
//! keeping touching regions apart. Each loop is simple: counter-clockwise
//! an outside, clockwise a hole.

type P = [f64; 2];

fn sub(a: P, b: P) -> P {
    [a[0] - b[0], a[1] - b[1]]
}
fn add(a: P, b: P) -> P {
    [a[0] + b[0], a[1] + b[1]]
}
fn mul(a: P, s: f64) -> P {
    [a[0] * s, a[1] * s]
}
fn dot(a: P, b: P) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: P, b: P) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn norm(a: P) -> f64 {
    a[0].hypot(a[1])
}

/// A point within this of an edge (in ems) lies on it: what rounding
/// leaves of where two arcs meet, nearly tangent (the square root of a
/// float's precision).
const EPS: f64 = 1e-8;
/// Points closer than this are one (40 nm on letters 40 mm to the em).
const MERGE: f64 = 1e-6;
/// How far off an edge its sides are judged (ems).
const SIDE: f64 = 4e-6;

/// An edge from a to b: a line where `bu` is 0, else an arc of sweep
/// 4 atan(bu), positive counter-clockwise (lying right of the chord).
#[derive(Clone, Copy, Debug)]
pub(super) struct Edge {
    pub a: P,
    pub b: P,
    pub bu: f64,
}

/// An arc's circle: centre, radius, the angle of a from it, and the sweep.
struct Circle {
    c: P,
    r: f64,
    start: f64,
    sweep: f64,
}

impl Edge {
    fn line(&self) -> bool {
        self.bu.abs() < 1e-12
    }

    fn circle(&self) -> Circle {
        let ch = sub(self.b, self.a);
        let l = norm(ch);
        let sweep = 4.0 * self.bu.atan();
        let n = [-ch[1] / l, ch[0] / l];
        let c = add(mul(add(self.a, self.b), 0.5), mul(n, l / (2.0 * (sweep / 2.0).tan())));
        let u = sub(self.a, c);
        Circle { c, r: norm(u), start: u[1].atan2(u[0]), sweep }
    }

    /// The point a fraction u along it.
    fn at(&self, u: f64) -> P {
        if self.line() {
            return add(self.a, mul(sub(self.b, self.a), u));
        }
        let k = self.circle();
        let t = k.start + k.sweep * u;
        add(k.c, [k.r * t.cos(), k.r * t.sin()])
    }

    /// Its unit tangent a fraction u along it.
    fn tangent(&self, u: f64) -> P {
        if self.line() {
            let d = sub(self.b, self.a);
            return mul(d, 1.0 / norm(d));
        }
        let k = self.circle();
        let t = k.start + k.sweep * u;
        let s = k.sweep.signum();
        [-t.sin() * s, t.cos() * s]
    }

    /// Where p lies along it (a fraction), if it lies on it within eps.
    fn param_of(&self, p: P, eps: f64) -> Option<f64> {
        if self.line() {
            let d = sub(self.b, self.a);
            let l2 = dot(d, d);
            let u = dot(sub(p, self.a), d) / l2;
            let q = add(self.a, mul(d, u));
            let l = l2.sqrt();
            return (norm(sub(p, q)) <= eps && u * l >= -eps && (u - 1.0) * l <= eps).then_some(u.clamp(0.0, 1.0));
        }
        let k = self.circle();
        let v = sub(p, k.c);
        if (norm(v) - k.r).abs() > eps {
            return None;
        }
        let mut phi = v[1].atan2(v[0]) - k.start;
        let tau = std::f64::consts::TAU;
        if k.sweep > 0.0 {
            phi = phi.rem_euclid(tau);
            if phi > k.sweep + (tau - k.sweep) / 2.0 {
                phi -= tau;
            }
        } else {
            phi = -(-phi).rem_euclid(tau);
            if phi < k.sweep - (tau + k.sweep) / 2.0 {
                phi += tau;
            }
        }
        let u = phi / k.sweep;
        (u * k.sweep.abs() * k.r >= -eps && (u - 1.0) * k.sweep.abs() * k.r <= eps).then_some(u.clamp(0.0, 1.0))
    }

    /// The piece from fraction u0 to u1.
    fn piece(&self, u0: f64, u1: f64, a: P, b: P) -> Edge {
        Edge { a, b, bu: if self.line() { 0.0 } else { (4.0 * self.bu.atan() * (u1 - u0) / 4.0).tan() } }
    }

    fn bbox(&self) -> (P, P) {
        let (mut lo, mut hi) = ([self.a[0].min(self.b[0]), self.a[1].min(self.b[1])], [self.a[0].max(self.b[0]), self.a[1].max(self.b[1])]);
        if !self.line() {
            for k in 1..16 {
                let p = self.at(k as f64 / 16.0);
                for i in 0..2 {
                    lo[i] = lo[i].min(p[i]);
                    hi[i] = hi[i].max(p[i]);
                }
            }
            // The arc between samples bulges at most r (1 - cos(sweep / 32)).
            let k = self.circle();
            let m = k.r * (1.0 - (k.sweep / 32.0).cos()) + EPS;
            lo = [lo[0] - m, lo[1] - m];
            hi = [hi[0] + m, hi[1] + m];
        }
        (lo, hi)
    }

    /// The angle it turns round p (radians, counter-clockwise positive).
    fn winding_angle(&self, p: P) -> f64 {
        let (u, v) = (sub(self.a, p), sub(self.b, p));
        let chord = cross(u, v).atan2(dot(u, v));
        if self.line() {
            return chord;
        }
        // p between the chord and the arc: the arc goes round it the long
        // way.
        let k = self.circle();
        let side = cross(sub(self.b, self.a), sub(p, self.a));
        let inside = norm(sub(p, k.c)) < k.r && side * k.sweep < 0.0;
        if inside { chord + std::f64::consts::TAU * k.sweep.signum() } else { chord }
    }
}

/// Where edges e and f meet: points on both.
fn meets(e: &Edge, f: &Edge) -> Vec<P> {
    let mut out = Vec::new();
    // Each one's ends lying on the other (touching, and overlapping runs).
    for (p, g) in [(e.a, f), (e.b, f), (f.a, e), (f.b, e)] {
        if g.param_of(p, EPS * 10.0).is_some() {
            out.push(p);
        }
    }
    match (e.line(), f.line()) {
        (true, true) => {
            let (d, g) = (sub(e.b, e.a), sub(f.b, f.a));
            let den = cross(d, g);
            if den.abs() > 1e-15 * norm(d) * norm(g) {
                let w = sub(f.a, e.a);
                let (t, u) = (cross(w, g) / den, cross(w, d) / den);
                if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
                    out.push(add(e.a, mul(d, t)));
                }
            }
        }
        (true, false) | (false, true) => {
            let (l, a) = if e.line() { (e, f) } else { (f, e) };
            let k = a.circle();
            let d = sub(l.b, l.a);
            let w = sub(l.a, k.c);
            let (qa, qb, qc) = (dot(d, d), 2.0 * dot(w, d), dot(w, w) - k.r * k.r);
            let disc = qb * qb - 4.0 * qa * qc;
            if disc >= 0.0 {
                let s = disc.sqrt();
                for t in [(-qb - s) / (2.0 * qa), (-qb + s) / (2.0 * qa)] {
                    if (0.0..=1.0).contains(&t) {
                        let p = add(l.a, mul(d, t));
                        if a.param_of(p, EPS * 10.0).is_some() {
                            out.push(p);
                        }
                    }
                }
            }
        }
        (false, false) => {
            let (k1, k2) = (e.circle(), f.circle());
            let d = norm(sub(k2.c, k1.c));
            if d > EPS && d <= k1.r + k2.r && d >= (k1.r - k2.r).abs() {
                let a = (k1.r * k1.r - k2.r * k2.r + d * d) / (2.0 * d);
                let h = (k1.r * k1.r - a * a).max(0.0).sqrt();
                let u = mul(sub(k2.c, k1.c), 1.0 / d);
                let m = add(k1.c, mul(u, a));
                for s in [-1.0, 1.0] {
                    let p = add(m, mul([-u[1], u[0]], s * h));
                    if e.param_of(p, EPS * 10.0).is_some() && f.param_of(p, EPS * 10.0).is_some() {
                        out.push(p);
                    }
                }
            }
        }
    }
    out
}

/// The winding number of `edges` round p.
fn winding(edges: &[Edge], p: P) -> i32 {
    (edges.iter().map(|e| e.winding_angle(p)).sum::<f64>() / std::f64::consts::TAU).round() as i32
}

/// The region `loops` fill (nonzero), as simple loops of edges, filled on
/// their left.
pub(super) fn filled(loops: &[Vec<Edge>]) -> Vec<Vec<Edge>> {
    let edges: Vec<Edge> = loops.iter().flatten().copied().filter(|e| norm(sub(e.b, e.a)) > MERGE).collect();
    // Points, merged within EPS.
    let mut pts: Vec<P> = Vec::new();
    let id = |p: P, pts: &mut Vec<P>| -> usize {
        if let Some(i) = pts.iter().position(|q| norm(sub(*q, p)) <= MERGE) {
            return i;
        }
        pts.push(p);
        pts.len() - 1
    };
    // Every edge's cut points: its ends and where others meet it.
    let boxes: Vec<(P, P)> = edges.iter().map(|e| e.bbox()).collect();
    let mut cuts: Vec<Vec<(f64, P)>> = edges.iter().map(|e| vec![(0.0, e.a), (1.0, e.b)]).collect();
    for i in 0..edges.len() {
        for j in i + 1..edges.len() {
            let ((la, ha), (lb, hb)) = (boxes[i], boxes[j]);
            if la[0] > hb[0] + EPS || lb[0] > ha[0] + EPS || la[1] > hb[1] + EPS || lb[1] > ha[1] + EPS {
                continue;
            }
            for p in meets(&edges[i], &edges[j]) {
                for k in [i, j] {
                    if let Some(u) = edges[k].param_of(p, EPS * 100.0) {
                        cuts[k].push((u, p));
                    }
                }
            }
        }
    }
    // The pieces, each between consecutive cuts.
    let mut pieces: Vec<(usize, usize, Edge)> = Vec::new();
    for (e, cs) in edges.iter().zip(&mut cuts) {
        cs.sort_by(|x, y| x.0.total_cmp(&y.0));
        let mut run: Vec<(f64, usize)> = Vec::new();
        for &(u, p) in cs.iter() {
            let k = id(p, &mut pts);
            if run.last().is_none_or(|l| l.1 != k) {
                run.push((u, k));
            }
        }
        for w in run.windows(2) {
            let ((u0, a), (u1, b)) = (w[0], w[1]);
            pieces.push((a, b, e.piece(u0, u1, pts[a], pts[b])));
        }
    }
    // Kept: filled on exactly one side; turned to have it on the left; one
    // of each pair that coincide.
    let all: Vec<Edge> = pieces.iter().map(|p| p.2).collect();
    let mut kept: Vec<(usize, usize, Edge)> = Vec::new();
    for &(a, b, e) in &pieces {
        let (m, t) = (e.at(0.5), e.tangent(0.5));
        let n = [-t[1], t[0]];
        let (left, right) = (winding(&all, add(m, mul(n, SIDE))) != 0, winding(&all, sub(m, mul(n, SIDE))) != 0);
        if left == right {
            continue;
        }
        let k = if left { (a, b, e) } else { (b, a, Edge { a: e.b, b: e.a, bu: -e.bu }) };
        if !kept.iter().any(|q| q.0 == k.0 && q.1 == k.1 && (q.2.bu - k.2.bu).abs() < 1e-9) {
            kept.push(k);
        }
    }
    // Joined into loops: at a point where several leave, the one turning
    // most to the left, which keeps regions that touch there apart.
    let mut used = vec![false; kept.len()];
    let mut out = Vec::new();
    for s in 0..kept.len() {
        if used[s] {
            continue;
        }
        used[s] = true;
        let mut lp = vec![kept[s].2];
        let mut at = kept[s].1;
        while at != kept[s].0 {
            let din = lp.last().expect("an edge").tangent(1.0);
            let next = (0..kept.len()).filter(|&i| !used[i] && kept[i].0 == at).max_by(|&i, &j| {
                let turn = |k: usize| {
                    let d = kept[k].2.tangent(0.0);
                    cross(din, d).atan2(dot(din, d))
                };
                turn(i).total_cmp(&turn(j))
            });
            let Some(i) = next else { break };
            used[i] = true;
            lp.push(kept[i].2);
            at = kept[i].1;
        }
        if at == kept[s].0 && lp.len() >= 2 {
            out.push(lp);
        }
    }
    out
}
