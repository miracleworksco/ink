//! polygon against the exact signed distance to its lines and arcs, at
//! random points and densely off its edges: the sign is exact and the value
//! never exceeds the distance.

use super::*;

fn ev(t: &T, x: f64, y: f64) -> f64 {
    match t {
        T::Lin(a, b, c, _) => a * x + b * y + c,
        T::Disc(cx, cy, r, _) => (x - cx).hypot(y - cy) - r,
        T::Max(v) => v.iter().map(|t| ev(t, x, y)).fold(f64::MIN, f64::max),
        T::Min(v) => v.iter().map(|t| ev(t, x, y)).fold(f64::MAX, f64::min),
        T::Neg(a) => -ev(a, x, y),
    }
}

/// The exact signed distance to a profile of lines and arcs: to its nearest
/// edge (a segment, or an arc: off its circle within its sweep, else to its
/// nearer end), negative inside, by the crossings of a ray along +x with its
/// edges (an arc in pieces monotone in y, split at its circle's top and
/// bottom; each piece's ends half-open in y, as a segment's).
fn exact(pts: &[[f64; 3]], x: f64, y: f64) -> f64 {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let (n, p) = (pts.len(), [x, y]);
    let mut d = f64::MAX;
    let mut inside = false;
    let mut cross_at = |a: P, q: P, xq: &dyn Fn() -> f64| {
        if (a[1] > y) != (q[1] > y) && x < xq() {
            inside = !inside;
        }
    };
    for i in 0..n {
        let (a, q, b) = ([pts[i][0], pts[i][1]], [pts[(i + 1) % n][0], pts[(i + 1) % n][1]], pts[i][2]);
        if b == 0.0 {
            cross_at(a, q, &|| a[0] + (y - a[1]) * (q[0] - a[0]) / (q[1] - a[1]));
            let e = sub(q, a);
            let t = (dot(sub(p, a), e) / dot(e, e)).clamp(0.0, 1.0);
            d = d.min(norm(sub(p, add(a, mul(e, t)))));
            continue;
        }
        // From a, turning by the sweep (b > 0: counter-clockwise).
        let arc = arc_of(a, q, b);
        let (u, ua) = (sub(p, arc.c), sub(a, arc.c));
        let (t0, sweep) = (ua[1].atan2(ua[0]), 4.0 * b.atan());
        let turn = (b.signum() * cross(ua, u)).atan2(dot(ua, u)).rem_euclid(TAU);
        d = d.min(if turn <= sweep.abs() { (norm(u) - arc.r).abs() } else { norm(sub(p, a)).min(norm(sub(p, q))) });
        // Its pieces, split where it passes its circle's top or bottom.
        let mut cuts: Vec<(f64, P)> = vec![(0.0, a)];
        let first = ((t0 - FRAC_PI_2) / PI).floor() as i64 - 1;
        for k in first..first + 6 {
            let th = FRAC_PI_2 + k as f64 * PI;
            let along = (th - t0) * sweep.signum();
            if along > 0.0 && along < sweep.abs() {
                cuts.push((along, [arc.c[0], arc.c[1] + if k.rem_euclid(2) == 0 { arc.r } else { -arc.r }]));
            }
        }
        cuts.sort_by(|x, y| x.0.total_cmp(&y.0));
        cuts.push((sweep.abs(), q));
        for w in cuts.windows(2) {
            let mid = t0 + sweep.signum() * 0.5 * (w[0].0 + w[1].0);
            let side = mid.cos().signum();
            cross_at(w[0].1, w[1].1, &|| arc.c[0] + side * (arc.r * arc.r - (y - arc.c[1]).powi(2)).max(0.0).sqrt());
        }
    }
    if inside { -d } else { d }
}

/// The field's sign against the exact section's, densely: at points off each
/// edge (along its normal, at distances from 1e-6 to 1 of the profile's
/// size) and over its box. Its sign is the section's but where both are
/// within polygon's allowance (100 eps) of 0; it is never more than the
/// distance.
fn signs(f: impl Fn(f64, f64) -> f64, sec: &[[f64; 3]], what: &str) {
    let n = sec.len();
    let (lo, hi) = (0..2).fold(([f64::MAX; 2], [f64::MIN; 2]), |(mut lo, mut hi), k| {
        for p in sec {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
        (lo, hi)
    });
    let size = norm(sub(hi, lo));
    let allow = 100.0 * 1e-9 * size;
    let probe = |x: f64, y: f64| {
        let (v, d) = (f(x, y), exact(sec, x, y));
        assert!((v < 0.0) == (d < 0.0) && v != 0.0 || (d.abs() <= allow && v.abs() <= allow), "{what}: sign at ({x}, {y}): {v} against {d}");
        assert!(v.abs() <= d.abs() * (1.0 + 1e-9) + allow, "{what}: over the distance at ({x}, {y}): {v} against {d}");
    };
    for i in 0..n {
        let (a, q, b) = ([sec[i][0], sec[i][1]], [sec[(i + 1) % n][0], sec[(i + 1) % n][1]], sec[i][2]);
        for k in 1..24 {
            let u = k as f64 / 24.0;
            let (at, nrm) = if b == 0.0 {
                (add(a, mul(sub(q, a), u)), right(a, q))
            } else {
                let arc = arc_of(a, q, b);
                let ua = sub(a, arc.c);
                let t = ua[1].atan2(ua[0]) + 4.0 * b.atan() * u;
                let dir = [t.cos(), t.sin()];
                (add(arc.c, mul(dir, arc.r)), dir)
            };
            for e in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0] {
                for side in [-1.0, 1.0] {
                    let p = add(at, mul(nrm, side * e * size));
                    probe(p[0], p[1]);
                }
            }
        }
    }
    for i in 0..=40 {
        for j in 0..=40 {
            probe(lo[0] - 0.1 * size + (hi[0] - lo[0] + 0.2 * size) * i as f64 / 40.0, lo[1] - 0.1 * size + (hi[1] - lo[1] + 0.2 * size) * j as f64 / 40.0);
        }
    }
}

struct Rng(u64);
impl Rng {
    fn f(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.f()
    }
}

/// A star of n points, radii in [r0, r1], with bulges on some edges.
fn star(rng: &mut Rng, n: usize, r0: f64, r1: f64, bulge: f64, cw: bool) -> Vec<[f64; 3]> {
    let mut v: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let a = std::f64::consts::TAU * (i as f64 + rng.range(-0.2, 0.2)) / n as f64;
            let r = rng.range(r0, r1);
            let b = if rng.f() < 0.4 { rng.range(-bulge, bulge) } else { 0.0 };
            [r * a.cos(), r * a.sin(), b]
        })
        .collect();
    if cw {
        let m = v.len();
        v = (0..m).map(|k| [v[m - 1 - k][0], v[m - 1 - k][1], -v[(2 * m - 2 - k) % m][2]]).collect();
    }
    v
}

fn simple(pts: &[[f64; 3]]) -> bool {
    // Densely sampled, as `exact` does, then no two segments cross.
    let n = pts.len();
    let mut poly: Vec<P> = Vec::new();
    for i in 0..n {
        let (a, q, b) = ([pts[i][0], pts[i][1]], [pts[(i + 1) % n][0], pts[(i + 1) % n][1]], pts[i][2]);
        for s in 0..40 {
            let u = s as f64 / 40.0;
            poly.push(if b == 0.0 {
                add(a, mul(sub(q, a), u))
            } else {
                let arc = arc_of(a, q, b);
                let a0 = (a[1] - arc.c[1]).atan2(a[0] - arc.c[0]);
                add(arc.c, mul([(a0 + 4.0 * b.atan() * u).cos(), (a0 + 4.0 * b.atan() * u).sin()], arc.r))
            });
        }
    }
    let m = poly.len();
    let lp = Loop { pts: poly, bs: vec![0.0; m], eps: 1e-9, defs: vec![Def::Half; m], ids: (0..m as u32).collect(), next: m as u32 };
    lp.crossing().is_none()
}

fn check(pts: &[[f64; 3]], rng: &mut Rng, what: &str) {
    let t = polygon(pts).unwrap_or_else(|e| panic!("{what}: {e}: {pts:?}"));
    let (lo, hi) = (0..2).fold(([f64::MAX; 2], [f64::MIN; 2]), |(mut lo, mut hi), k| {
        for p in pts {
            lo[k] = lo[k].min(p[k] - 8.0);
            hi[k] = hi[k].max(p[k] + 8.0);
        }
        (lo, hi)
    });
    for _ in 0..200 {
        let (x, y) = (rng.range(lo[0], hi[0]), rng.range(lo[1], hi[1]));
        let (f, d) = (ev(&t, x, y), exact(pts, x, y));
        // The dense sampling is within ~1e-5 of the arcs.
        if d.abs() > 1e-4 {
            assert_eq!(f < 0.0, d < 0.0, "{what}: sign at ({x}, {y}): {f} against {d}: {pts:?}");
        }
        assert!(f.abs() <= d.abs() + 1e-4, "{what}: over the distance at ({x}, {y}): {f} against {d}: {pts:?}");
    }
    // On each half-plane's line and each disc's circle, where two terms
    // could cancel to a zero off the boundary (a crease: an interval proof
    // cannot cross it, a ray can stop on it), the sign is strict and right.
    fn leaves(t: &T, out: &mut Vec<T>) {
        match t {
            T::Lin(..) | T::Disc(..) => out.push(t.clone()),
            T::Max(v) | T::Min(v) => v.iter().for_each(|t| leaves(t, out)),
            T::Neg(a) => leaves(a, out),
        }
    }
    let mut ls = Vec::new();
    leaves(&t, &mut ls);
    let span = norm(sub(hi, lo));
    for l in &ls {
        for k in 0..8 {
            let s = rng.f();
            let [x, y] = match *l {
                T::Lin(a, b, c, _) => add(mul([a, b], -c), mul([-b, a], (s - 0.5) * 2.0 * span)),
                T::Disc(cx, cy, r, _) => add([cx, cy], mul([(s * std::f64::consts::TAU).cos(), (s * std::f64::consts::TAU).sin()], r)),
                _ => unreachable!(),
            };
            if !(lo[0]..hi[0]).contains(&x) || !(lo[1]..hi[1]).contains(&y) {
                continue;
            }
            let (f, d) = (ev(&t, x, y), exact(pts, x, y));
            if d.abs() > 1e-3 {
                assert!(f.abs() > 1e-9 && (f < 0.0) == (d < 0.0), "{what}: on a term's zero ({k}) at ({x}, {y}): {f} against {d}: {pts:?}");
            }
        }
    }
}

#[test]
fn collinear_arc_chords_keep_their_bumps() {
    let b = QUARTER;
    let pts = [
        [-22.75, -13.75, b],
        [-20.25, -16.25, -b],
        [-16.25, -20.25, b],
        [-13.75, -22.75, 0.0],
        [13.75, -22.75, b],
        [16.25, -20.25, -b],
        [20.25, -16.25, b],
        [22.75, -13.75, 0.0],
        [22.75, 13.75, b],
        [20.25, 16.25, -b],
        [16.25, 20.25, b],
        [13.75, 22.75, 0.0],
        [-13.75, 22.75, b],
        [-16.25, 20.25, -b],
        [-20.25, 16.25, b],
        [-22.75, 13.75, 0.0],
    ];
    let t = polygon(&pts).unwrap();
    let v = ev(&t, -22.017766952966, -15.517766952966);
    assert!(v.abs() <= 100.0 * allowance(&pts), "on the convex r2.5 arc: {v}");
    let v = ev(&t, -21.95, -15.45);
    assert!(v < 0.0, "inside the convex r2.5 arc: {v}");
    signs(|x, y| ev(&t, x, y), &pts, "collinear tangent arcs");
    // Straight, convex and concave edges may share the same hull line.
    // Their curves still bound separate circular segments.
    for bs in [[-0.2, -0.2], [0.2, 0.2], [0.2, 0.0], [0.0, 0.2], [0.2, -0.2]] {
        let pts = [[0.0, 0.0, bs[0]], [4.0, 0.0, bs[1]], [8.0, 0.0, 0.0], [8.0, 6.0, 0.0], [0.0, 6.0, 0.0]];
        let t = polygon(&pts).unwrap();
        signs(|x, y| ev(&t, x, y), &pts, &format!("collinear chords {bs:?}"));
    }
}

#[test]
fn a_circular_wave_strip_keeps_its_trough() {
    // An arc next to a line it crosses: 1561T4's native tangent circular strip. The lower and
    // upper boundaries at x=-4.7625 are 4.639204166667 and 4.969404166667.
    let pts = [
        [-9.525, 6.2055375, -0.160222927551],
        [-7.14375, 5.422370833, 0.328888888889],
        [-2.38125, 5.422370833, -0.328888888889],
        [2.38125, 5.422370833, 0.328888888889],
        [7.14375, 5.422370833, -0.160222927551],
        [9.525, 6.2055375, 0.0],
        [9.525, 6.5357375, 0.160222927551],
        [6.947752424, 5.688109408, -0.328888888889],
        [2.577247576, 5.688109408, 0.328888888889],
        [-2.577247576, 5.688109408, -0.328888888889],
        [-6.947752424, 5.688109408, 0.160222927551],
        [-9.525, 6.5357375, 0.0],
    ];
    let t = polygon(&pts).unwrap();
    let v = ev(&t, -4.7625, 4.8);
    assert!(v < 0.0, "inside the wave strip's trough: {v}");
    signs(|x, y| ev(&t, x, y), &pts, "wave strip");
}

#[test]
fn random_tangent_arc_and_line_strips_have_the_winding_sign() {
    let mut rng = Rng(0xa703_6824_771e_903d);
    for case in 0..32 {
        let angle = rng.range(0.15, 0.7);
        let thickness = rng.range(0.1, 0.6);
        let count = 2 + (rng.f() * 5.0) as usize;
        let (mut p, mut heading) = ([0.0, 0.0], 0.0f64);
        let (mut lower, mut upper) = (Vec::new(), Vec::new());
        // Integrate arcs and tangent lines; the other boundary is their
        // normal offset, with the same sweep and a changed radius.
        let mut edge = |turn: f64, size: f64| {
            let b = (turn / 4.0).tan();
            lower.push([p[0], p[1], b]);
            upper.push([p[0] - thickness * heading.sin(), p[1] + thickness * heading.cos(), b]);
            if turn == 0.0 {
                p = add(p, mul([heading.cos(), heading.sin()], size));
            } else {
                let c = add(p, mul([-heading.sin(), heading.cos()], size * turn.signum()));
                heading += turn;
                p = sub(c, mul([-heading.sin(), heading.cos()], size * turn.signum()));
            }
        };
        edge(0.0, rng.range(0.5, 3.0));
        edge(-angle, rng.range(2.0, 8.0));
        for k in 0..count {
            // With no intervening line these chords are collinear, even
            // with unequal radii: precisely the case star profiles miss.
            if case % 3 == 0 {
                edge(0.0, rng.range(0.2, 2.0));
            }
            edge(if k % 2 == 0 { 2.0 * angle } else { -2.0 * angle }, rng.range(2.0, 8.0));
        }
        edge(if count.is_multiple_of(2) { angle } else { -angle }, rng.range(2.0, 8.0));
        edge(0.0, rng.range(0.5, 3.0));
        edge(0.0, 0.0); // The cap, then the other boundary back to the start.
        let mut pts = lower;
        for k in (0..upper.len()).rev() {
            pts.push([upper[k][0], upper[k][1], if k == 0 { 0.0 } else { -upper[k - 1][2] }]);
        }
        // Either winding, in an arbitrary frame.
        if case % 2 == 0 {
            let n = pts.len();
            pts = (0..n).map(|k| [pts[n - 1 - k][0], pts[n - 1 - k][1], -pts[(2 * n - 2 - k) % n][2]]).collect();
        }
        let a = rng.range(0.0, std::f64::consts::TAU);
        let shift = [rng.range(-10.0, 10.0), rng.range(-10.0, 10.0)];
        for p in &mut pts {
            let q = add([a.cos() * p[0] - a.sin() * p[1], a.sin() * p[0] + a.cos() * p[1]], shift);
            (p[0], p[1]) = (q[0], q[1]);
        }
        let t = polygon(&pts).unwrap_or_else(|e| panic!("tangent strip {case}: {e}: {pts:?}"));
        signs(|x, y| ev(&t, x, y), &pts, &format!("tangent strip {case}"));
    }
}

#[test]
fn polygons_and_arcs() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let q = (std::f64::consts::PI / 8.0).tan();
    // A rounded rectangle and a slot.
    let rr = [[4.0, -6.0, q], [6.0, -4.0, 0.0], [6.0, 4.0, q], [4.0, 6.0, 0.0], [-4.0, 6.0, q], [-6.0, 4.0, 0.0], [-6.0, -4.0, q], [-4.0, -6.0, 0.0]];
    let slot = [[-8.0, -2.0, 0.0], [8.0, -2.0, 1.0], [8.0, 2.0, 0.0], [-8.0, 2.0, 1.0]];
    for (p, name) in [(&rr[..], "rounded rectangle"), (&slot[..], "slot")] {
        check(p, &mut rng, name);
    }
    let mut tried = 0;
    while tried < 150 {
        let n = 3 + (rng.f() * 24.0) as usize;
        let cw = rng.f() < 0.5;
        let (r0, bulge) = if tried % 3 == 0 { (10.0, 0.0) } else { (rng.range(3.0, 9.0), rng.range(0.0, 0.4)) };
        let pts = star(&mut rng, n, r0, 12.0, bulge, cw);
        if !simple(&pts) {
            continue;
        }
        tried += 1;
        check(&pts, &mut rng, "star");
    }
}

/// 1678-26a-0409's outer loop: a circle with six lobes and concave fillets
/// between, 24 tangent arcs, each fillet meeting a lobe or the circle at an
/// inflection. Its hull is the lobes' common tangents; its pockets' fillets
/// are not the segments added past their chords, which would make the field
/// 0 on a chord inside a fillet's disc (2e-5 at (154.7818, 306.9589),
/// outside; -4e-5 at (153.7563, 304.1839), inside).
const LOBES: [[f64; 3]; 24] = [
    [160.076496, 320.509146, -0.071943],
    [163.615709, 322.552512, 0.219551],
    [164.786824, 323.788159, -0.717273],
    [173.033176, 323.788159, 0.219551],
    [174.204291, 322.552512, -0.071943],
    [177.743504, 320.509146, 0.219551],
    [179.399163, 320.112755, -0.717273],
    [183.522339, 312.971204, 0.219551],
    [183.037794, 311.339165, -0.071943],
    [183.037794, 307.252435, 0.219551],
    [183.522339, 305.620396, -0.717273],
    [179.399163, 298.478845, 0.219551],
    [177.743504, 298.082454, -0.071943],
    [174.204291, 296.039088, 0.219551],
    [173.033176, 294.803441, -0.717273],
    [164.786824, 294.803441, 0.219551],
    [163.615709, 296.039088, -0.071943],
    [160.076496, 298.082454, 0.219551],
    [158.420837, 298.478845, -0.717273],
    [154.297661, 305.620396, 0.219551],
    [154.782206, 307.252435, -0.071943],
    [154.782206, 311.339165, 0.219551],
    [154.297661, 312.971204, -0.717273],
    [158.420837, 320.112755, 0.219551],
];

/// A plate's corners with ears (a lobe between two fillets; its first
/// contact is discarded when renumbered), and a run of S bends.
const EARS: [[f64; 3]; 16] = [
    [91.2812, 35.3019, -0.12284188688065735],
    [91.8392, 37.5387, 0.7401667088289183],
    [87.545, 41.833, -0.12284275700161461],
    [85.3082, 41.275, 0.0],
    [-85.3082, 41.275, -0.12284188688065707],
    [-87.545, 41.833, 0.7400962681121066],
    [-91.8392, 37.5387, -0.12284275700161412],
    [-91.2812, 35.3019, 0.0],
    [-91.2812, -35.3019, -0.12284188688065735],
    [-91.8392, -37.5387, 0.7401667088289183],
    [-87.545, -41.833, -0.12284275700161461],
    [-85.3082, -41.275, 0.0],
    [85.3082, -41.275, -0.12284188688065707],
    [87.545, -41.833, 0.7400962681121066],
    [91.8392, -37.5387, -0.12284275700161412],
    [91.2812, -35.3019, 0.0],
];
const BENDS: [[f64; 3]; 24] = [
    [-14.4432, 15.49135, -0.131652],
    [-15.08125, 17.8726, -0.339454],
    [-11.55138, 22.47282, -0.339454],
    [-6.1943, 20.25385, 0.0],
    [-5.49926, 19.05, 0.267949],
    [0.0, 15.875, -0.267949],
    [13.74815, 7.9375, 0.267949],
    [19.24741, 4.7625, 0.0],
    [20.6375, 4.7625, -0.414214],
    [25.4, 0.0, -0.414214],
    [20.6375, -4.7625, 0.0],
    [19.24741, -4.7625, 0.267949],
    [13.74815, -7.9375, -0.267949],
    [0.0, -15.875, 0.267949],
    [-5.49926, -19.05, 0.0],
    [-6.1943, -20.25385, -0.339454],
    [-11.55138, -22.47282, -0.339454],
    [-15.08125, -17.8726, -0.131652],
    [-14.4432, -15.49135, 0.0],
    [-13.74815, -14.2875, 0.267949],
    [-13.74815, -7.9375, -0.131652],
    [-15.875, 0.0, -0.131652],
    [-13.74815, 7.9375, 0.267949],
    [-13.74815, 14.2875, 0.0],
];

/// Stars of arcs with every corner rounded, as their solved sections.
fn rounded_stars(rng: &mut Rng, count: usize) -> Vec<Vec<[f64; 3]>> {
    let mut out = Vec::new();
    while out.len() < count {
        let n = 3 + (rng.f() * 12.0) as usize;
        let cw = rng.f() < 0.5;
        let pts: Vec<[f64; 4]> = star(rng, n, 5.0, 12.0, 0.3, cw).iter().map(|p| [p[0], p[1], p[2], rng.range(0.3, 1.5)]).collect();
        let Ok(sec) = walls::Walls::new(&pts, &vec![0.0; n]).and_then(|w| w.solve(0.0)) else { continue };
        if simple(&sec.0) {
            out.push(sec.0);
        }
    }
    out
}

#[test]
fn tangent_arcs_leave_no_crease() {
    // Arcs meeting tangent with either curvature: 0409's lobes, a plate's
    // corners with ears (a lobe between two fillets), a run of S bends, and
    // stars of arcs with every corner rounded. The sign is exact off the
    // boundary at random points and on every term's line or circle.
    let mut rng = Rng(7);
    let t = polygon(&LOBES).unwrap();
    for (x, y) in [(154.7818, 306.9589), (153.7563, 304.1839)] {
        let (f, d) = (ev(&t, x, y), exact(&LOBES, x, y));
        assert!(f.abs() > 1e-4 && (f < 0.0) == (d < 0.0), "0409 at ({x}, {y}): {f} against {d}");
    }
    check(&LOBES, &mut rng, "0409's outer loop");
    // Its lobes are split once each side, where the hull's common tangents
    // touch them: the hull and six pockets, 114 terms (with the hull split
    // at its farthest points again and again, 19 pieces and 198 terms).
    let (t, _, tr) = build(&LOBES, allowance(&LOBES)).unwrap();
    fn terms(t: &T) -> usize {
        match t {
            T::Lin(..) | T::Disc(..) => 1,
            T::Max(v) | T::Min(v) => v.iter().map(terms).sum(),
            T::Neg(a) => terms(a),
        }
    }
    let pieces = tr.pieces;
    assert!(pieces <= 7 && terms(&t) <= 114, "0409's outer loop: {pieces} pieces, {} terms", terms(&t));
    check(&EARS, &mut rng, "ears");
    check(&BENDS, &mut rng, "S bends");
    for sec in rounded_stars(&mut rng, 40) {
        check(&sec, &mut rng, "rounded star");
    }
}

static CORE: crate::Registry = crate::Registry::new(&[]);

/// The field `src` evaluates to, at points.
fn field(src: &str) -> impl Fn([f64; 3]) -> f64 + use<> {
    let d = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(f, _) = d.def("f").unwrap().clone() else { panic!() };
    let t = crate::tape::Tape::new(&d.g, f);
    move |p| crate::eval::eval_point(&t, &crate::eval::Ctx::new(0.0, vec![]), p)
}

fn ink_points(pts: &[[f64; 4]]) -> String {
    let v: Vec<String> = pts.iter().map(|p| format!("(vec {} {} {} {})", p[0], p[1], p[2], p[3])).collect();
    format!("(vec {})", v.join(" "))
}

/// An annular window between circles of radii 11 and 20 and two spokes
/// 2.4 from its axis, each corner rounded to 1.5 (a moulded wheel's
/// windows), counter-clockwise and not; an L (a pocket under a lid),
/// rounded inside and at the corners away from its pocket's mouth, and
/// rounded at its mouth too (the lid tangent to the round); and a hexagon
/// with a bump between fillets on one side (the hull's lines from two
/// corners touch the bump; its joins are smooth, so it is drafted alike).
/// Each with its drafts.
fn profiles() -> Vec<(&'static str, Vec<[f64; 4]>, Vec<f64>)> {
    let (ri, ro, hs, ang, r) = (11.0f64, 20.0f64, 2.4f64, std::f64::consts::PI / 5.0, 1.5);
    let (n, v) = ([-ang.cos(), ang.sin()], [ang.sin(), ang.cos()]);
    let on = |rad: f64| add(mul(n, hs), mul(v, (rad * rad - hs * hs).sqrt()));
    let (a, b) = (on(ri), on(ro));
    let half = |p: P| (p[0].atan2(p[1]) / 2.0).tan();
    let w = vec![[a[0], a[1], 0.0, r], [b[0], b[1], half(b), r], [-b[0], b[1], 0.0, r], [-a[0], a[1], -half(a), r]];
    let k = w.len();
    let cw: Vec<[f64; 4]> = (0..k).map(|i| [w[k - 1 - i][0], w[k - 1 - i][1], -w[(2 * k - 2 - i) % k][2], w[k - 1 - i][3]]).collect();
    let l = vec![[0.0, 0.0, 0.0, 1.0], [12.0, 0.0, 0.0, 1.0], [12.0, 4.0, 0.0, 0.0], [4.0, 4.0, 0.0, 1.0], [4.0, 10.0, 0.0, 0.0], [0.0, 10.0, 0.0, 1.0]];
    let mouth = vec![[0.0, 0.0, 0.0, 1.0], [12.0, 0.0, 0.0, 1.0], [12.0, 4.0, 0.0, 1.5], [4.0, 4.0, 0.0, 1.0], [4.0, 10.0, 0.0, 1.0], [0.0, 10.0, 0.0, 1.0]];
    let bump = vec![
        [7.3323, 0.0, 0.0, 0.0],
        [3.6662, 6.35, 0.0, 0.0],
        [2.4539, 6.35, -0.19272798, 0.0],
        [1.2269, 6.8411, 0.40032568, 0.0],
        [-1.2269, 6.8411, -0.19272798, 0.0],
        [-2.4539, 6.35, 0.0, 0.0],
        [-3.6662, 6.35, 0.0, 0.0],
        [-7.3323, 0.0, 0.0, 0.0],
        [-3.6662, -6.35, 0.0, 0.0],
        [3.6662, -6.35, 0.0, 0.0],
    ];
    let varied = |p: &Vec<[f64; 4]>| (0..p.len()).map(|i| [3.0, -2.0, 1.0, 4.0, 2.0, -1.0][i % 6]).collect();
    let alike = vec![2.0; bump.len()];
    vec![("window", varied(&w), w), ("window, clockwise", varied(&cw), cw), ("L", varied(&l), l), ("L rounded at its mouth", varied(&mouth), mouth), ("bump", alike, bump)]
        .into_iter()
        .map(|(n, d, p)| (n, p, d))
        .collect()
}

#[test]
fn rounded_corners_are_tangent_arcs() {
    // A rectangle rounded at its corners is the rectangle of lines and
    // quarter arcs.
    let q = (std::f64::consts::PI / 8.0).tan();
    let w = walls::Walls::new(&[[-6.0, -6.0, 0.0, 2.0], [6.0, -6.0, 0.0, 2.0], [6.0, 6.0, 0.0, 2.0], [-6.0, 6.0, 0.0, 2.0]], &[0.0; 4]).unwrap();
    let (got, _) = w.solve(0.0).unwrap();
    let want = [[-4.0, -6.0, 0.0], [4.0, -6.0, q], [6.0, -4.0, 0.0], [6.0, 4.0, q], [4.0, 6.0, 0.0], [-4.0, 6.0, q], [-6.0, 4.0, 0.0], [-6.0, -4.0, q]];
    for k in 0..8 {
        let g = got[(k + 1) % 8];
        assert!((0..3).all(|i| (g[i] - want[k][i]).abs() < 1e-12), "{got:?}");
    }
    let mut rng = Rng(5);
    for (name, p, _) in profiles() {
        let (pts, _) = walls::Walls::new(&p, &vec![0.0; p.len()]).unwrap().solve(0.0).unwrap();
        check(&pts, &mut rng, name);
    }
}

#[test]
fn drafted_profiles_lean_their_walls_and_keep_their_rounds() {
    // Each wall leans by its own draft and each round keeps its radius:
    // against the profile solved at each height, the field's sign is the
    // same and its value never over the distance, at the span's ends and
    // between (corners, rounds and the hull's contacts on their curves;
    // within the dense sampling's 1e-4); beyond the span the profile is
    // straight. The field is 1-Lipschitz in x, y and z, out to far beyond
    // its box.
    let mut rng = Rng(11);
    for (name, p, drafts) in profiles() {
        let (z0, z1) = (-6.0, 8.0);
        let ds: Vec<String> = drafts.iter().map(|d| d.to_string()).collect();
        let f = field(&format!("(polygon {} (vec {}) {z0} {z1})", ink_points(&p), ds.join(" ")));
        let w = walls::Walls::new(&p, &drafts.iter().map(|d| d.to_radians().tan()).collect::<Vec<_>>()).unwrap();
        for k in 0..40 {
            let z = match k {
                0 => z0,
                1 => z1,
                2 => z1 + 3.0,
                _ => rng.range(z0, z1),
            };
            let (sec, _) = w.solve(z.clamp(z0, z1)).unwrap();
            let tol = 1e-7;
            for _ in 0..60 {
                let (x, y) = (rng.range(-24.0, 24.0), rng.range(-4.0, 26.0));
                let (v, d) = (f([x, y, z]), exact(&sec, x, y));
                if d.abs() > tol + 1e-4 {
                    assert_eq!(v < 0.0, d < 0.0, "{name}: sign at ({x}, {y}, {z}): {v} against {d}");
                }
                assert!(v.abs() <= d.abs() + tol + 1e-4, "{name}: over the distance at ({x}, {y}, {z}): {v} against {d}");
            }
        }
        for _ in 0..4000 {
            let p: [f64; 3] = [rng.range(-80.0, 80.0), rng.range(-80.0, 80.0), rng.range(-40.0, 40.0)];
            let s = rng.range(1e-3, 30.0);
            let q = [p[0] + rng.range(-s, s), p[1] + rng.range(-s, s), p[2] + rng.range(-s, s)];
            let gap = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            assert!((f(p) - f(q)).abs() <= gap * (1.0 + 1e-9) + 1e-12, "{name}: not 1-Lipschitz from {p:?} to {q:?}");
        }
    }
}

#[test]
fn a_drafted_circle_leans_out_either_way_round() {
    // A circle of two half arcs has no chord area: its turn is its arcs'
    // (clockwise, it leaned in: 1678-26c-1311's drafted bore).
    for b in [1.0, -1.0] {
        let f = field(&format!("(polygon (vec (vec -4 0 {b}) (vec 4 0 {b})) 10 -5 5)"));
        for z in [-5.0, 0.0, 5.0] {
            let r = 4.0 + z * 10f64.to_radians().tan();
            assert!(f([r - 1e-3, 0.0, z]) < 0.0 && f([r + 1e-3, 0.0, z]) > 0.0, "bulge {b} at z = {z}");
        }
    }
}

#[test]
fn drafts_that_change_a_profile_are_refused() {
    let err = |src: &str| Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).err().unwrap().message;
    // A thin slot's walls meet when they lean in far enough, at z = 0.5 /
    // tan 10 = 2.836: it closes there, empty beyond (a convex profile is
    // the meet of its walls' half-planes).
    let f = field("(polygon (vec (vec 0 0) (vec 10 0) (vec 10 1) (vec 0 1)) (vec -10 0 -10 0) 0 10)");
    assert!(f([5.0, 0.5, 2.8]) < 0.0 && f([5.0, 0.5, 2.87]) > 0.0 && f([5.0, 0.5, 9.0]) > 0.0);
    // Tangent walls that lean differently part.
    let e = err("(polygon (vec (vec -5 -2) (vec 5 -2 1) (vec 5 2) (vec -5 2 1)) (vec 0 -2 0 -2) 0 5)");
    assert!(e.contains("no longer meet"), "{e}");
    // A round on a smooth corner.
    let e = err("(polygon (vec (vec -5 -2) (vec 5 -2 1 1) (vec 5 2) (vec -5 2 1)))");
    assert!(e.contains("smooth"), "{e}");
    // (A round at a pocket's mouth, the lid tangent to it, builds: its
    // contact moves on its curve, not the chord of its path:
    // drafted_profiles_lean_their_walls_and_keep_their_rounds.)
    // Walls crossing at an end are refused saying where it holds:
    // WCP-1708-002's wedge as written (6 decimals): its triangle closes at
    // z = -1.130932, 1e-6 inside the span's end, six times the allowance for
    // its size, and its walls cross from there on.
    let e = err("(polygon (vec (vec 13.752523 -6.746875 0.004319) (vec 14.516837 -6.502706) (vec 14.516837 -6.746875)) (vec 0 34.051877 0) -1.130933 1.130933)");
    assert!(e.contains("walls cross, at z = -1.130933") && e.contains("it holds from z = -1.13093"), "{e}");
}

#[test]
fn a_curves_derivative_holds_its_value() {
    // The field's bounds on a curve's rate are its derivative's intervals:
    // they hold the real value. d(t / 3) is not folded to the float nearest
    // 1/3, a point below it.
    use curve::{Ex, kx, xdiv, xlin};
    use ink_core::interval::I;
    let d = Ex::d(&xdiv(&xlin(0.0, 1.0), &kx(3.0))).at(I::point(0.5));
    let third = 1.0 / 3.0; // below 1/3
    assert!(d.lo <= third && d.hi > third, "{d:?}");
}

/// A REV omni hub's moulded slab: a 44-edge loop drafted 1 and 2 degrees.
/// A pocket's arc ends tangent to its lid at the lid's corner, so it
/// reaches the lid exactly at every height: its reach is a length times a
/// moving radius, and taken at their extremes over the span it would be
/// 0.63 past the lid (refused at z0 itself). Decided exactly, it holds, and the
/// field is the profile solved at each height (the_sign_is_the_sections_densely).
fn slab() -> (Vec<[f64; 4]>, Vec<f64>, f64, f64) {
    let p: Vec<[f64; 4]> = vec![
        [-15.5315, 17.0928, 0.0, 0.0],
        [-18.0788, 21.5048, 0.475177, 0.0],
        [-19.2932, 21.7014, 0.161624, 0.0],
        [-28.4406, 5.8577, 0.475177, 0.0],
        [-27.6631, 4.90429, 0.0, 0.0],
        [-22.5685, 4.90429, 0.0, 0.0],
        [-19.737, 0.0, 0.0, 0.0],
        [-22.5685, -4.90429, 0.0, 0.0],
        [-27.6631, -4.90429, 0.475177, 0.0],
        [-28.4406, -5.8577, 0.161624, 0.0],
        [-19.2932, -21.7014, 0.475177, 0.0],
        [-18.0788, -21.5048, 0.0, 0.0],
        [-15.5315, -17.0928, 0.0, 0.0],
        [-9.86851, -17.0928, 0.0, 0.0],
        [-7.03702, -21.997, 0.0, 0.0],
        [-9.58429, -26.409, 0.475177, 0.0],
        [-9.14736, -27.5591, 0.0802909, 0.0],
        [-0.0, -29.0375, 0.0802909, 0.0],
        [9.14736, -27.5591, 0.475177, 0.0],
        [9.58429, -26.409, 0.0, 0.0],
        [7.03702, -21.997, 0.0, 0.0],
        [9.86851, -17.0928, 0.0, 0.0],
        [15.5315, -17.0928, 0.0, 0.0],
        [18.0788, -21.5048, 0.475177, 0.0],
        [19.2932, -21.7014, 0.161624, 0.0],
        [28.4406, -5.8577, 0.475177, 0.0],
        [27.6631, -4.90429, 0.0, 0.0],
        [22.5685, -4.90429, 0.0, 0.0],
        [19.737, -0.0, 0.0, 0.0],
        [22.5685, 4.90429, 0.0, 0.0],
        [27.6631, 4.90429, 0.475177, 0.0],
        [28.4406, 5.8577, 0.161624, 0.0],
        [19.2932, 21.7014, 0.475177, 0.0],
        [18.0788, 21.5048, 0.0, 0.0],
        [15.5315, 17.0928, 0.0, 0.0],
        [9.86851, 17.0928, 0.0, 0.0],
        [7.03702, 21.997, 0.0, 0.0],
        [9.58429, 26.409, 0.475177, 0.0],
        [9.14736, 27.5591, 0.0802909, 0.0],
        [0.0, 29.0375, 0.0802909, 0.0],
        [-9.14736, 27.5591, 0.475177, 0.0],
        [-9.58429, 26.409, 0.0, 0.0],
        [-7.03702, 21.997, 0.0, 0.0],
        [-9.86851, 17.0928, 0.0, 0.0],
    ];
    let drafts = [2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1].map(f64::from);
    (p, drafts.to_vec(), -0.746125, 0.841375)
}

#[test]
fn the_sign_is_the_sections_densely() {
    // The sign everywhere off the boundary is the solved section's, densely
    // (`signs`): the profiles of tangent arcs, and drafted, at many heights,
    // its curves lowered as the proof reads them (a round's centre on its
    // curve under the tangent from it: the drafted window, its field +4e-6
    // at a point 1e-6 inside were the tangent to read the centre's chord;
    // and a corner where a line meets an arc wall nearly tangentially, whose
    // field would be 0 were its rate bounded over fixed parts of the span).
    // The points those would get wrong: 1e-6 inside the window (-1.0e-6 is
    // the distance), and in the near-tangent profile (-1.088).
    let (_, w, _) = profiles().remove(0);
    let f = field(&format!("(polygon {} (vec 7.5 -5 2.5 10) -6 8)", ink_points(&w)));
    let v = f([3.383960840828903, 10.281094040663824, 1.0]);
    assert!(v < 0.0, "the window 1e-6 inside: {v}");
    let f = field("(polygon (vec (vec -1 -1 1) (vec 1 -1) (vec 1.001 -0.98) (vec 1.001 3) (vec -1 3)) 10 0 1)");
    let v = f([0.0, 0.0, 0.5]);
    assert!(v < -0.1, "the near-tangent profile: {v}");
    let mut rng = Rng(13);
    let mut flat: Vec<(String, Vec<[f64; 3]>)> = vec![("0409's outer loop".into(), LOBES.to_vec()), ("ears".into(), EARS.to_vec()), ("S bends".into(), BENDS.to_vec())];
    flat.extend(rounded_stars(&mut rng, 40).into_iter().map(|s| ("rounded star".into(), s)));
    for (name, p, _) in profiles() {
        flat.push((name.into(), walls::Walls::new(&p, &vec![0.0; p.len()]).unwrap().solve(0.0).unwrap().0));
    }
    use rayon::prelude::*;
    flat.par_iter().for_each(|(name, sec)| {
        let t = polygon(sec).unwrap_or_else(|e| panic!("{name}: {e}"));
        signs(|x, y| ev(&t, x, y), sec, name);
    });
}

#[test]
fn the_sign_of_a_drafted_profile_is_its_sections_at_every_height() {
    // The profiles above, drafted, at 24 heights.
    use rayon::prelude::*;
    let window = profiles().remove(0).1;
    let mut drafted: Vec<(String, Vec<[f64; 4]>, Vec<f64>, f64, f64)> = profiles().into_iter().map(|(n, p, d)| (n.to_string(), p, d, -6.0, 8.0)).collect();
    drafted.push(("window drafted 7.5, -5, 2.5, 10".into(), window, vec![7.5, -5.0, 2.5, 10.0], -6.0, 8.0));
    let (p, d, z0, z1) = slab();
    drafted.push(("omni slab".into(), p, d, z0, z1));
    let near = vec![[-1.0, -1.0, 1.0, 0.0], [1.0, -1.0, 0.0, 0.0], [1.001, -0.98, 0.0, 0.0], [1.001, 3.0, 0.0, 0.0], [-1.0, 3.0, 0.0, 0.0]];
    drafted.push(("a line meeting an arc wall nearly tangentially".into(), near, vec![10.0; 5], 0.0, 1.0));
    drafted.par_iter().for_each(|(name, p, drafts, z0, z1)| drafted_signs(p, drafts, *z0, *z1, 24, name).unwrap());
}

/// A drafted profile's field against its section solved at each of k + 1
/// heights over its span (`signs`; at an end where the profile degenerates,
/// 1e-6 of the span inside it); where polygon refuses it, why.
fn drafted_signs(p: &[[f64; 4]], drafts: &[f64], z0: f64, z1: f64, k: usize, name: &str) -> Result<(), String> {
    let ds: Vec<String> = drafts.iter().map(|d| d.to_string()).collect();
    let src = format!("(polygon {} (vec {}) {z0} {z1})", ink_points(p), ds.join(" "));
    if let Err(e) = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]) {
        return Err(format!("{name}: {}: {src}", e.message));
    }
    let f = field(&src);
    let w = walls::Walls::new(p, &drafts.iter().map(|d| d.to_radians().tan()).collect::<Vec<_>>()).unwrap();
    use rayon::prelude::*;
    (0..=k).into_par_iter().for_each(|i| {
        let mut z = z0 + (z1 - z0) * i as f64 / k as f64;
        if w.solve(z).is_err() && (i == 0 || i == k) {
            z += (z1 - z0) * if i == 0 { 1e-6 } else { -1e-6 };
        }
        let sec = w.solve(z).unwrap().0;
        signs(|x, y| f([x, y, z]), &sec, &format!("{name} at z = {z}"));
    });
    Ok(())
}

/// Sections of moulded hubs as step2ink reads them (6 decimals, bulges 5):
/// REV-21-2476's rim between two roller pockets, and AndyMark am-3340's
/// hub half at two heights. Each has rounds meeting arcs nearly tangent,
/// some at a reflex corner (a r 48.9 arc into a r 1.54 round, 0.02 rad;
/// a r 73 arc into a r 0.13 round, 0.0044 rad: a bridge 13 eps deep).
const HUBS: [&[[f64; 3]]; 3] = [
    &[
        [0.0, 0.0, 0.0],
        [50.057767, -36.288534, 0.0],
        [61.056055, -9.736318, 0.0],
        [63.617511, -9.401372, 0.0],
        [71.567189, -12.694236, -0.020852],
        [71.864856, -12.832024, 0.507218],
        [72.047687, -12.740877, 0.087739],
        [72.0477, 12.740801, 0.506242],
        [71.864851, 12.832055, -0.020942],
        [71.567189, 12.694236, 0.0],
        [63.617511, 9.401372, 0.0],
        [61.056055, 9.736318, 0.0],
        [50.057767, 36.288534, 0.0],
        [51.632147, 38.3366, 0.0],
        [59.581825, 41.629464, -0.020852],
        [59.889738, 41.742516, 0.507218],
        [59.954568, 41.936247, 0.087739],
        [41.93631, 59.954524, 0.0],
    ],
    &[
        [-10.1487, 14.3474, 0.0],
        [-19.9321, 31.3669, 0.26694],
        [-21.2992, 32.163, -0.67498],
        [-24.5013, 40.2358, 0.0],
        [-24.5755, 40.3649, 0.43156],
        [-26.822, 40.9011, 0.07137],
        [-37.2391, 31.7108, 0.29185],
        [-37.4985, 30.0769, 0.29134],
        [-36.1675, 29.1004, -0.43797],
        [-29.0969, 19.9684, 0.31266],
        [-28.328, 18.3318, 0.0],
        [-14.3091, 10.2027, 0.36167],
        [-12.3078, 10.5431, -0.0395],
        [-10.4966, 12.3474, 0.36167],
    ],
    &[
        [-36.1872, 29.0904, -0.43913],
        [-29.1103, 19.9586, 0.31308],
        [-28.3419, 18.3197, 0.0],
        [-14.3053, 10.1803, 0.3619],
        [-12.3031, 10.5217, -0.03991],
        [-10.4753, 12.3427, 0.3619],
        [-10.1263, 14.3435, 0.0],
        [-19.9187, 31.3784, 0.26751],
        [-21.2891, 32.1746, -0.67647],
        [-24.4867, 40.2337, 0.0],
        [-24.5716, 40.3813, 0.43148],
        [-26.8177, 40.9178, 0.07155],
        [-37.2631, 31.7006, 0.29214],
        [-37.5211, 30.0652, 0.29163],
    ],
];

#[test]
fn rounded_data_is_read_exactly() {
    // Tangent arcs and rounds whose numbers are rounded, as a STEP reader
    // writes them, leave joins a little off tangent either way: hull lines
    // that miss an arc by less than its sag, reflex corners a hair deep.
    // Each builds (the hull, followed round the arcs, keeps one turned past
    // its start; a reflex corner needs no cut the contacts leave within its
    // gap), its sign exact off the boundary and its value never over the
    // distance (`check`, `signs`): the hubs, the stars rounded to 2, 4 and
    // 6 decimals, and a repeated point (the edge on from it is its last
    // copy's, not the first's bulge).
    let mut rng = Rng(21);
    for (k, sec) in HUBS.iter().enumerate() {
        let t = polygon(sec).unwrap_or_else(|e| panic!("hub {k}: {e}"));
        check(sec, &mut rng, "hub");
        signs(|x, y| ev(&t, x, y), sec, "hub");
    }
    let stars = rounded_stars(&mut rng, 150);
    for (k, sec) in stars.iter().enumerate() {
        let q = [1e2, 1e4, 1e6][k % 3];
        let r: Vec<[f64; 3]> = sec.iter().map(|p| [(p[0] * q).round() / q, (p[1] * q).round() / q, (p[2] * q * 10.0).round() / q / 10.0]).collect();
        if !simple(&r) {
            continue;
        }
        let t = polygon(&r).unwrap_or_else(|e| panic!("rounded star: {e}: {r:?}"));
        check(&r, &mut rng, "rounded star");
        signs(|x, y| ev(&t, x, y), &r, "rounded star");
    }
    let repeated = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, 0.0, 0.4], [10.0, 10.0, 0.0], [0.0, 10.0, 0.0]];
    check(&repeated, &mut rng, "a repeated point");
}

#[test]
fn edges_in_an_arcs_segment_are_kept_out() {
    // The profile is its chord polygon with each arc's segment added or
    // taken away, so nothing else may reach into a segment. A quarter turn
    // of radius 20 closed by two lines meeting 3 inside it (a crescent: its
    // chord triangle is simple and turns the other way, which would make
    // the field +2.08 at a point 0.91 inside); and the same arc over a
    // polyline of 8 lines g inside it (a curved band: the polyline crosses
    // the arc's chord, so the chord polygon's area would turn the loop over;
    // at g = 1e-4 the arc is halved into 256 pieces before its chords clear
    // the polyline, past a cap of 8 n + 64 points).
    let mut rng = Rng(9);
    let q = (std::f64::consts::PI / 8.0).tan();
    let mut profiles = vec![vec![[20.0, 0.0, q], [0.0, 20.0, 0.0], [12.0, 12.0, 0.0]]];
    for g in [1.0, 1e-2, 1e-4] {
        let mut band = vec![[20.0, 0.0, q], [0.0, 20.0, 0.0]];
        band.extend((0..=8).map(|k| {
            let a = std::f64::consts::FRAC_PI_2 * (8 - k) as f64 / 8.0;
            [(20.0 - g) * a.cos(), (20.0 - g) * a.sin(), 0.0]
        }));
        profiles.push(band);
    }
    for p in &profiles {
        let t = polygon(p).unwrap_or_else(|e| panic!("{e}: {p:?}"));
        check(p, &mut rng, "arc closed in on");
        signs(|x, y| ev(&t, x, y), p, "arc closed in on");
    }
}

#[test]
fn drafts_that_change_a_construction_are_read_in_pieces() {
    // As the walls lean, a profile's construction can change partway: a
    // point crosses a hull line, a contact runs off its arc's end. Each such
    // span is read in pieces, not refused ("built differently"), each
    // proven over its part, and the field is the section's at every height.
    // A notch whose inner corners are rounded, drafted in 45 deg: the hull
    // line from the far corner touches the round, then runs to its end.
    let notch = vec![
        [10.0, -10.0, 0.0, 0.0],
        [10.0, -5.0, 0.0, 0.0],
        [6.45, -3.0, 0.0, 1.75],
        [6.45, 3.0, 0.0, 1.75],
        [10.0, 5.0, 0.0, 0.0],
        [10.0, 10.0, 0.0, 0.0],
        [-10.0, 10.0, 0.0, 0.0],
        [-10.0, -10.0, 0.0, 0.0],
    ];
    for (z0, z1) in [(0.0, 2.0), (1.0, 2.0)] {
        drafted_signs(&notch, &[-45.0; 8], z0, z1, 32, "notch").unwrap();
    }
    // AndyMark am-3340's hub (8 in MK mecanum) as moulded.py reads its
    // section at x 16.98 (4 decimals), drafted 2 deg: joins tangent only to
    // that rounding: a reflex join a hair deep is not refused as its corner
    // changing kind, nor two arcs meeting nearly tangent (r 49 and 1.5,
    // 6e-4 rad) as an arc's corner turning out past its neighbour by the
    // larger radius's sag.
    for (k, hub) in HUB_LOOPS.iter().enumerate() {
        drafted_signs(hub, &vec![2.0; hub.len()], -2.075, 2.075, 8, &format!("hub loop {k}")).unwrap();
    }
}

/// Loops of am-3340's hub section (drafts_that_change_a_construction_are_read_in_pieces).
const HUB_LOOPS: [&[[f64; 4]]; 3] = [
    &[
        [-5.5098, 14.8959, -0.04722, 0.0],
        [-2.6173, 15.665, 0.36603, 0.0],
        [-1.2919, 17.2273, 0.0, 0.0],
        [-1.2912, 17.5968, 0.58805, 0.0],
        [-3.6305, 18.9039, -0.1393, 0.0],
        [-6.1631, 18.1843, 0.0, 0.0],
        [-6.249, 18.185, 0.58268, 0.0],
        [-7.6055, 15.888, 0.0, 0.0],
        [-7.4361, 15.5934, 0.36603, 0.0],
    ],
    &[
        [38.6014, 2.0965, -0.27734, 0.0],
        [42.545, 4.4662, -0.23715, 0.0],
        [46.1289, 2.6651, 0.35135, 0.0],
        [47.9561, 2.1228, 0.35211, 0.0],
        [48.9853, 3.7328, 0.0688, 0.0],
        [46.1347, 16.8844, 0.29739, 0.0],
        [44.811, 17.9178, 0.29688, 0.0],
        [43.3027, 17.1875, -0.25714, 0.0],
        [36.8451, 13.6313, -0.18169, 0.0],
        [31.8128, 15.5224, 0.32052, 0.0],
        [29.9765, 15.7041, 0.0, 0.0],
        [15.5934, 7.4361, 0.36603, 0.0],
        [14.8959, 5.5098, -0.04722, 0.0],
        [15.665, 2.6173, 0.36603, 0.0],
        [17.2273, 1.2919, 0.0, 0.0],
        [37.197, 1.2542, 0.27786, 0.0],
    ],
    &[
        [42.545, -4.4662, -0.27023, 0.0],
        [38.6583, -2.2, 0.26974, 0.0],
        [37.28, -1.3946, 0.0, 0.0],
        [28.4581, -1.378, 0.58268, 0.0],
        [27.1383, -3.6963, -0.42831, 0.0],
        [25.2671, -10.4137, 0.58805, 0.0],
        [25.2193, -13.0929, 0.0, 0.0],
        [29.8709, -15.7903, 0.32558, 0.0],
        [31.7297, -15.5962, -0.18545, 0.0],
        [36.8451, -13.6313, -0.25714, 0.0],
        [43.3027, -17.1874, 0.29688, 0.0],
        [44.811, -17.9178, 0.29739, 0.0],
        [46.1347, -16.8844, 0.0739, 0.0],
        [49.0509, -2.7389, 0.38244, 0.0],
        [47.7295, -1.0852, 0.24345, 0.0],
        [46.56, -1.9562, -0.28686, 0.0],
    ],
];

#[test]
fn rounded_data_is_drafted_exactly() {
    // Stars of tangent arcs and rounds, their numbers rounded as a STEP
    // reader writes them (2, 4 and 6 decimals: joins tangent only to that),
    // drafted at angles each way over spans of several heights, some edges
    // leaning more than others (so hulls and pockets change as they lean).
    // At heights over its span each one's field is its section's. Each
    // builds, including where the profile built again where its
    // construction stops holding is built the same way. The stars side by
    // side (each alone takes the test 15 s).
    use rayon::prelude::*;
    let mut rng = Rng(23);
    let stars = rounded_stars(&mut rng, 60);
    let results: Vec<Result<(), String>> = stars
        .par_iter()
        .enumerate()
        .filter_map(|(k, sec)| {
            let q = [1e2, 1e4, 1e6][k % 3];
            let p: Vec<[f64; 4]> = sec.iter().map(|p| [(p[0] * q).round() / q, (p[1] * q).round() / q, (p[2] * q * 10.0).round() / q / 10.0, 0.0]).collect();
            let r: Vec<[f64; 3]> = p.iter().map(|p| [p[0], p[1], p[2]]).collect();
            if !simple(&r) {
                return None;
            }
            let a = [-7.0, -2.0, 2.0, 7.0][k % 4];
            let drafts: Vec<f64> = (0..p.len()).map(|i| if k % 5 == 0 { a * [1.0, 0.5, 1.5][i % 3] } else { a }).collect();
            let (z0, z1) = [(-3.0, 3.0), (0.0, 5.0), (-5.0, 1.0)][k % 3];
            // Drafted in, a thin part can close: only spans solved simple.
            let w = walls::Walls::new(&p, &drafts.iter().map(|d| d.to_radians().tan()).collect::<Vec<_>>()).unwrap();
            if (0..=8).any(|i| w.solve(z0 + (z1 - z0) * i as f64 / 8.0).map_or(true, |s| !simple(&s.0))) {
                return None;
            }
            Some(drafted_signs(&p, &drafts, z0, z1, 4, &format!("star {k}")))
        })
        .collect();
    let tried = results.len();
    let refused: Vec<String> = results.into_iter().filter_map(Result::err).collect();
    assert!(tried >= 30 && refused.is_empty(), "{} of {tried} refused: {refused:#?}", refused.len());
}

#[test]
// Its points are rounded data: a -6.28 there is a coordinate, not tau.
#[allow(clippy::approx_constant)]
fn a_construction_built_alike_where_it_stops_is_built_again_farther_on() {
    // Where a drafted profile's construction stops holding, the profile
    // built again there can be built the same way (its build's tolerance
    // takes the change as not yet made), and each is not refused ("built
    // again there it does not carry on"): rounded data whose arc runs out past a
    // hull line; a nut strip, its 45 degree chamfers drafted from
    // a slot's rounds; REV's MPM saddle mount's plate, chamfered 45 degrees
    // into an undrafted fillet: at its bottom face the chamfer is tangent to
    // the fillet, so above it the corner turns reflex at a rate without
    // bound. The next construction is built where the profile builds
    // differently, and the field is the section's at every height.
    let rounded: Vec<[f64; 4]> = [
        [10.02, -1.53, 0.456],
        [10.72, 0.51, 0.0],
        [7.85, 5.47, 0.124],
        [7.48, 5.86, 0.0],
        [3.3, 8.44, 0.021],
        [3.27, 8.46, 0.0],
        [-2.05, 11.18, 0.328],
        [-2.94, 11.03, 0.107],
        [-6.97, 4.62, 0.13],
        [-7.01, 4.14, 0.0],
        [-6.26, -0.28, 0.003],
        [-6.26, -0.28, 0.0],
        [-5.74, -3.11, 0.22],
        [-5.02, -4.13, 0.0],
        [-2.38, -5.65, 0.125],
        [-1.7, -5.85, -0.113],
        [2.16, -6.84, 0.419],
        [3.88, -6.28, 0.0],
        [5.35, -3.32, -0.215],
        [5.8, -2.95, 0.03],
    ]
    .map(|p| [p[0], p[1], p[2], 0.0])
    .to_vec();
    drafted_signs(&rounded, &vec![7.0; rounded.len()], -3.0, 3.0, 12, "rounded data").unwrap();
    let strip: Vec<[f64; 4]> = [
        [-22.225, 6.1, 0.0],
        [22.225, 6.1, 0.0],
        [22.225, 2.0193, 1.0],
        [22.225, -2.0193, 0.0],
        [22.225, -6.1, 0.0],
        [5.735063, -6.1, -0.190935],
        [5.519564, -5.55625, 0.0],
        [5.519564, -4.92125, 0.414214],
        [4.725814, -4.1275, 0.0],
        [-15.832306, -4.1275, -0.145114],
        [-16.264929, -3.999239, 0.082904],
        [-21.62422, -1.663528, -0.343985],
        [-22.225, -0.893591, 0.0],
    ]
    .map(|p| [p[0], p[1], p[2], 0.0])
    .to_vec();
    let drafts = [45.0, 0.0, 0.0, 0.0, 45.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    drafted_signs(&strip, &drafts, -0.25, 0.25, 12, "nut strip").unwrap();
    // Beyond the span the profile is straight, and its field is its end
    // section's: the strip's last piece, built just past where a hull line's
    // contact is born, has that line turning so fast its scale is 4e-6, and
    // the bar's body beyond read -1.6e-6 at 4.1 deep.
    let ds: Vec<String> = drafts.iter().map(|d| d.to_string()).collect();
    let f = field(&format!("(polygon {} (vec {}) -0.25 0.25)", ink_points(&strip), ds.join(" ")));
    for z in [-3.0, 0.5, 5.5] {
        let v = f([-10.0, 0.0, z]);
        assert!(v <= -(z.abs() - 0.25).min(4.1275) + 1e-9, "the strip's body at z = {z}: {v}");
    }
    let plate: Vec<[f64; 4]> = [
        [-11.46338, -51.554063, -0.198912],
        [-11.807997, -51.411317, 0.0],
        [-23.471317, -39.747997, -0.198912],
        [-23.614062, -39.40338, 0.0],
        [-23.614062, 39.403432, -0.198912],
        [-23.471317, 39.74805, 0.0],
        [-20.69805, 42.521317, -0.198912],
        [-20.353432, 42.664062, 0.0],
        [20.35338, 42.664062, -0.198912],
        [20.697997, 42.521317, 0.0],
        [23.662613, 39.556701, -0.009867],
        [23.675942, 39.542836, 0.0],
        [32.053351, 30.476652, -0.403232],
        [32.7025, 29.802599, 0.0],
        [32.7025, -11.088539, -0.279351],
        [32.340152, -11.686481, 0.0],
        [23.673373, -39.467099, -0.121339],
        [23.552742, -39.666572, 0.0],
        [11.807997, -51.411317, -0.198912],
        [11.46338, -51.554063, 0.0],
        [6.661305, -51.554062, -0.199387],
        [6.183469, -51.7525, 0.0],
        [-6.183469, -51.7525, -0.199387],
        [-6.661305, -51.554062, 0.0],
    ]
    .map(|p| [p[0], p[1], p[2], 0.0])
    .to_vec();
    let chamfers = [-45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, 0.0, 0.0, 0.0, -45.0, -45.0, -45.0, -45.0, -45.0, 0.0, 0.0, 0.0, -45.0];
    drafted_signs(&plate, &chamfers, -0.198437, 0.198437, 12, "MPM saddle mount").unwrap();
}

#[test]
fn a_drafted_profile_may_degenerate_at_its_ends() {
    // A chamfer drafted to where it runs out at a corner: at the span's end
    // an edge, or the whole section, reaches zero (in rounded data, a hair
    // past it), and none is refused there ("its walls cross"). Each holds
    // over the span to within the allowance of its ends, and the field is
    // the section's at every height (at an end, just inside it): 1678's
    // wedge (part 1115), a triangle whose 45 degree wall meets its far
    // corner at z0; part-44's strip, chamfered 45 degrees both sides to a
    // ridge; REV-21-3005-p05's needle; and AndyMark's winch spool's tooth,
    // whose 45 degree chamfered arc shrinks to the point its sides meet at,
    // 5e-8 past z0 as written.
    let wedge = [[36.83, 319.0875, 0.0, 0.0], [36.83, 317.5, 0.0, 0.0], [38.4175, 317.5, 0.0, 0.0]];
    drafted_signs(&wedge, &[45.0, 0.0, 0.0], -1.5875, 1.5875, 8, "wedge").unwrap();
    let ridge = [[-393.133162, -92.10675, 0.0, 0.0], [-392.68415, -92.42425, 0.0, 0.0], [-392.68415, -47.27575, 0.0, 0.0], [-393.133162, -47.59325, 0.0, 0.0]];
    drafted_signs(&ridge, &[30.0, 45.0, 30.0, 45.0], -0.224506, 0.224506, 8, "ridge").unwrap();
    let needle =
        [[-27.609593, -16.002, 0.0, 0.0], [-27.609593, -21.139, 0.0, 0.0], [-23.947809, -24.800784, 0.0, 0.0], [-28.206038, -21.139, 0.0, 0.0], [-28.206038, -16.002, 0.0, 0.0]];
    drafted_signs(&needle, &[0.0, 0.0, 21.156354, 30.691133, 0.0], -1.004882, 1.004882, 8, "needle").unwrap();
    let tooth = [[-0.805667, 6.94052, -0.057587, 0.0], [0.805813, 6.940436, 0.0, 0.0], [1.979946, 6.26255, 0.154314, 0.0], [-1.979946, 6.26255, 0.0, 0.0]];
    drafted_signs(&tooth, &[-45.0, 0.0, 0.0, 0.0], -0.418794, 0.418794, 8, "spool tooth").unwrap();
}

/// REV's MPM saddle mount's plate (`a_construction_built_alike_where_it_stops_is_built_again_farther_on`):
/// chamfered 45 degrees over z0..z1 into undrafted fillets at a tab.
fn mpm_plate() -> (Vec<[f64; 4]>, Vec<f64>, f64, f64) {
    let plate: Vec<[f64; 4]> = [
        [-11.46338, -51.554063, -0.198912],
        [-11.807997, -51.411317, 0.0],
        [-23.471317, -39.747997, -0.198912],
        [-23.614062, -39.40338, 0.0],
        [-23.614062, 39.403432, -0.198912],
        [-23.471317, 39.74805, 0.0],
        [-20.69805, 42.521317, -0.198912],
        [-20.353432, 42.664062, 0.0],
        [20.35338, 42.664062, -0.198912],
        [20.697997, 42.521317, 0.0],
        [23.662613, 39.556701, -0.009867],
        [23.675942, 39.542836, 0.0],
        [32.053351, 30.476652, -0.403232],
        [32.7025, 29.802599, 0.0],
        [32.7025, -11.088539, -0.279351],
        [32.340152, -11.686481, 0.0],
        [23.673373, -39.467099, -0.121339],
        [23.552742, -39.666572, 0.0],
        [11.807997, -51.411317, -0.198912],
        [11.46338, -51.554063, 0.0],
        [6.661305, -51.554062, -0.199387],
        [6.183469, -51.7525, 0.0],
        [-6.183469, -51.7525, -0.199387],
        [-6.661305, -51.554062, 0.0],
    ]
    .map(|p| [p[0], p[1], p[2], 0.0])
    .to_vec();
    let chamfers =
        [-45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, -45.0, 0.0, 0.0, 0.0, -45.0, -45.0, -45.0, -45.0, -45.0, 0.0, 0.0, 0.0, -45.0].to_vec();
    (plate, chamfers, -0.198437, 0.198437)
}

/// The distance from (x, y, z) to a drafted profile's walls (straight beyond
/// its span), from its sections every 1/200 of the span.
fn wall_distance(w: &walls::Walls, z0: f64, z1: f64) -> impl Fn(f64, f64, f64) -> f64 {
    let zs: Vec<f64> = (0..=200).map(|i| z0 + (z1 - z0) * i as f64 / 200.0).collect();
    let secs: Vec<Vec<[f64; 3]>> = zs.iter().map(|&z| w.solve(z).unwrap().0).collect();
    move |x, y, z| {
        let at = |k: usize, dz: f64| exact(&secs[k], x, y).abs().hypot(dz);
        let along = (0..zs.len()).map(|k| at(k, z - zs[k])).fold(f64::MAX, f64::min);
        let beyond = if z < z0 {
            at(0, 0.0)
        } else if z > z1 {
            at(200, 0.0)
        } else {
            f64::MAX
        };
        along.min(beyond)
    }
}

#[test]
fn a_drafted_wall_reads_its_plane() {
    // A drafted line wall is its plane, signed and linear: inside a square
    // drafted 45 degrees the field is the 3D distance to the nearer leaning
    // wall, and past its corner the nearer plane's (the walls mitre, as an
    // undrafted polygon's half-planes do). A wall scaled by how fast its ends
    // slide along it (sqrt 2 per mm of z here) reads 0.816 of the distance.
    let f = field("(polygon (vec (vec -5 -5) (vec 5 -5) (vec 5 5) (vec -5 5)) -45 0 2)");
    for (p, want) in [([3.8, 0.0, 1.0], -0.2 * 0.5f64.sqrt()), ([3.5, 0.0, 1.0], -0.5 * 0.5f64.sqrt()), ([4.5, 4.5, 1.0], 0.5 * 0.5f64.sqrt())] {
        let v = f(p);
        assert!((v - want).abs() < 1e-9, "at {p:?}: {v}, the distance {want}");
    }
    // So past a 1 degree square's corner (the distance to the box its wall
    // sweeps reads 0.885 of that plane's there).
    let f = field("(polygon (vec (vec -5 -5) (vec 5 -5) (vec 5 5) (vec -5 5)) 1 0 10)");
    let a = 1f64.to_radians();
    let want = (5.8 - (5.0 + 5.0 * a.tan())) * a.cos();
    let v = f([5.8, 5.8, 5.0]);
    assert!((v - want).abs() < 1e-9, "past the corner: {v}, the nearer wall's plane {want}");
}

#[test]
fn a_drafted_field_reads_the_distance_to_its_walls() {
    // MPM's plate, chamfered 45 degrees into undrafted fillets: the lines its
    // construction adds (hull lines over its pockets, its arcs' chords) turn
    // and grow as the walls lean. Each reads the distance across it, scaled
    // by its own rate near it (scaled once over the profile's box, the field
    // would read 0.001 of the distance inside: -0.0063 at 6.4 deep, a mesh
    // of it 0.3 off). Within 0.3 of the walls, at
    // heights its construction holds through, the field is mostly the 3D
    // distance to the walls (a median of 0.9 at each height), at least a
    // tenth of it where a corner slides fast along a fillet it meets nearly
    // tangentially (0.14 there), never over it, and its slope at most 1.
    // One construction tree stays within 4000 instructions; a separate
    // distance tree signed by the construction would duplicate that work.
    let (p, d, z0, z1) = mpm_plate();
    let ds: Vec<String> = d.iter().map(|d| d.to_string()).collect();
    let src = format!("(polygon {} (vec {}) {z0} {z1})", ink_points(&p), ds.join(" "));
    let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    assert!(tape.len() <= 4000, "{} instructions", tape.len());
    let w = walls::Walls::new(&p, &d.iter().map(|d| d.to_radians().tan()).collect::<Vec<_>>()).unwrap();
    let dist = wall_distance(&w, z0, z1);
    let mut rng = Rng(5);
    for z in [0.0, 0.1, 0.15] {
        let (sec, _) = w.solve(z).unwrap();
        let mut pts = Vec::new();
        while pts.len() < 60 {
            let (x, y) = (rng.range(-24.0, 33.0), rng.range(-52.0, 43.0));
            let e = exact(&sec, x, y);
            if e.abs() < 0.3 && e.abs() > 0.005 {
                pts.push(([x, y, z], e.signum()));
            }
        }
        let gs = crate::eval::eval_grads(&tape, &crate::eval::Ctx::new(0.0, vec![]), &pts.iter().map(|p| p.0).collect::<Vec<_>>());
        let mut ratios: Vec<f64> = pts
            .iter()
            .zip(&gs)
            .map(|(&([x, y, z], side), g)| {
                let (v, t) = (g[0], dist(x, y, z));
                assert!(v * side >= 0.1 * t && v.abs() <= t * (1.0 + 1e-9) + 1e-9, "({x}, {y}, {z}): {v} against {t} from the walls");
                assert!(g[1].hypot(g[2]).hypot(g[3]) <= 1.0 + 1e-9, "({x}, {y}, {z}): slope {g:?}");
                v * side / t
            })
            .collect();
        ratios.sort_by(f64::total_cmp);
        assert!(ratios[30] >= 0.9, "at z = {z}: the median {} of the distance", ratios[30]);
    }
    let v = crate::eval::eval_point(&tape, &crate::eval::Ctx::new(0.0, vec![]), [0.0, -45.0, 0.0]);
    assert!(v < -0.15, "6.4 deep: {v}");
}

#[test]
fn a_drafted_field_s_slope_is_at_most_1_where_its_construction_changes() {
    // MPM's plate where two drafted pieces hand over: weighted across their
    // overlap, (f_next - f_prev) times the weight's slope (1 / its 0.01 mm)
    // read -2.66 and -0.79 4e-12 apart in z, a slope over 4e11.
    let (p, d, z0, z1) = mpm_plate();
    let ds: Vec<String> = d.iter().map(|d| d.to_string()).collect();
    let src = format!("(polygon {} (vec {}) {z0} {z1})", ink_points(&p), ds.join(" "));
    let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let at = |z: f64| crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [29.0, 23.0, z].map(ink_core::interval::I::point), None);
    let (a, b) = (at(-0.168249779095), at(-0.168249779091));
    let apart = (b.lo - a.hi).max(a.lo - b.hi).max(0.0);
    assert!(apart <= (-0.168249779091f64 - -0.168249779095f64), "{a:?} {b:?}");
}

#[test]
fn a_polygon_across_any_magnitudes_is_built() {
    // From 1e-100 to 1e70: its small differences' products, at one scale
    // for the whole loop, fall below the floats; it was refused as
    // doubling back.
    let src = "(def f (polygon (vec (vec 0 0) (vec 1e-100 0) (vec 2e-100 1e-100) (vec 1e70 1e70) (vec 0 1e70))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let at = |p: [f64; 3]| crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), p.map(ink_core::interval::I::point), None);
    assert!(at([1e60, 1e69, 0.0]).hi < 0.0 && at([1e69, 1e60, 0.0]).lo > 0.0);
}

#[test]
fn a_straight_line_is_at_most_1_lipschitz_in_real_numbers() {
    // The triangle (0, 0) (3, 4) (0, 4): its slanted edge scaled by 1 / 5
    // rounded to nearest (above 1 / 5) read 1 + 2^-54 between (3/2, 9/4) and
    // (7/4, 33/16), in exact arithmetic.
    let doc = Doc::load_str_in(&CORE, "(def f (polygon (vec (vec 0 0) (vec 3 4) (vec 0 4))))", "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    // |f(q) - f(p)| <= 5/16.
    assert!(exact::differs_by_at_most(&doc.g, fid, [1.5, 2.25], [1.75, 2.0625], 0.3125), "slope above 1");
}

#[test]
fn a_straight_polygon_does_not_drift_over_many_corners() {
    // A bottom of 8192 segments, corner i at (i / 8192, 1.2e-9 i ln i): each
    // corner within polygon's allowance of the line its neighbours' merges
    // left, they drifted to a chord 3.4e-6 off the drawn run, and a point
    // 1e-6 inside read outside. Built exactly, it reads inside.
    let n = 8192;
    let mut pts: Vec<String> = (0..=n)
        .map(|i| {
            let x = i as f64 / n as f64;
            let y = if i > 1 { 1.2e-9 * i as f64 * (i as f64).ln() } else { 0.0 };
            format!("(vec {x:?} {y:?})")
        })
        .collect();
    pts.extend(["(vec 1.0 1.0)".into(), "(vec 0.0 1.0)".into()]);
    let doc = Doc::load_str_in(&CORE, &format!("(def f (polygon (vec {})))", pts.join(" ")), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let i = (n / 2) as f64;
    let (x, y) = (0.5, 1.2e-9 * i * i.ln() + 1e-6);
    let v = crate::eval::eval_interval(
        &tape,
        &crate::eval::Ctx::new(0.0, vec![]),
        [ink_core::interval::I::point(x), ink_core::interval::I::point(y), ink_core::interval::I::point(0.0)],
        None,
    );
    assert!(v.hi < 0.0 && (v.lo + 1e-6).abs() < 1e-12, "{v:?}");
}

#[test]
fn a_drafted_line_is_read_about_its_middle() {
    // A 10 square, corners rounded 0.02, its bottom edge an arc bowed by
    // a bulge of 2.9e-9 (at the size where a bow is kept, not dropped),
    // drafted 45 degrees or mixed: the lines its construction adds there
    // slope alike read about their middle or an end, and read about an
    // end, the bound at a point 0.3 microns outside spanned 0 (lane E's
    // arcbound probes: 90 of 160 went unknown, 30 more excluded the side
    // they are on). Each point's bound there is on its side.
    let outline = |drafts: &str| format!("(polygon (vec (vec 0 0 -2.9e-09 0.02) (vec 10 0 0 0.02) (vec 10 10 0 0.02) (vec 0 10 0 0.02)) (vec {drafts}) 0 0.5)");
    let ctx = crate::eval::Ctx::new(0.0, vec![]);
    for (src, p) in [
        (outline("45 45 45 45"), [10.494142347710815, -0.4941423507183473, 0.5]),
        (outline("-17 7 23 -11"), [0.005857652285083452, 0.005857652319153936, 0.0]),
        (outline("-17 7 23 -11"), [0.032487754636732445, 0.047742755832634313, 0.137]),
        (outline("-17 7 23 -11"), [0.10304780685314548, 0.1587229936052508, 0.5]),
    ] {
        let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let Value::Field(f, _) = doc.def("f").unwrap().clone() else { panic!() };
        let t = crate::tape::Tape::new(&doc.g, f);
        let v = crate::eval::eval_interval(&t, &ctx, p.map(ink_core::interval::I::point), None);
        assert!(v.lo > 0.0, "{src} at {p:?}: [{:e}, {:e}] spans 0, 3e-7 outside", v.lo, v.hi);
    }
}

#[test]
fn a_drafted_profile_prunes_to_a_short_formula_in_small_cells() {
    // A rectangle stepped 0.2 down its bottom by a small concave arc,
    // drafted 45 degrees: the lines its construction adds there (the arc's
    // chord, lines to it) run between points under 0.05 apart, each moving
    // about 1 per unit z, so a normal read from two of them has no bound
    // over a cell taller than they are apart, and neither had the field
    // there however small the cell (MPM's saddle plate: 233 of its 3,636
    // instructions left at 0.006, three quarters of the cells its surface
    // was not near). Each line's promised slope bounds it instead. Followed
    // down an octree, the cells the surface may cross, each's bound at most
    // its diagonal wide and its pruned formula short.
    use ink_core::interval::I;
    let src = "(polygon (vec (vec -10 -0.2) (vec -10 10) (vec 10 10) (vec 10 -0.4) (vec 0 -0.4 -0.2) (vec -0.5 -0.2)) -45 -0.2 0.2)";
    let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let ctx = crate::eval::Ctx::new(0.0, vec![]);
    let root = std::sync::Arc::new(crate::tape::Tape::new(&doc.g, fid));
    // Cubes from one 24 across about it, 2^-12 of it at the finest.
    let mut cells = vec![([-12.0, -7.0, -12.0], 24.0, root)];
    let mut rng = Rng(0x9e3779b97f4a7c15);
    let (mut widest, mut lens) = (0.0f64, Vec::new());
    for level in 1..=12 {
        let mut next = Vec::new();
        for (o, s, t) in &cells {
            let h = 0.5 * s;
            for k in 0..8 {
                let a: [f64; 3] = std::array::from_fn(|j| o[j] + if k >> j & 1 == 1 { h } else { 0.0 });
                let mut ch = vec![crate::tape::BOTH; t.len()];
                let v = crate::eval::eval_interval(t, &ctx, a.map(|x| I::new(x, x + h)), Some(&mut ch));
                if v.lo <= 0.0 && v.hi >= 0.0 {
                    if level == 12 {
                        widest = widest.max((v.hi - v.lo) / (h * 3f64.sqrt()));
                        lens.push(t.simplify(&ch).len());
                    }
                    next.push((a, h, std::sync::Arc::new(t.simplify(&ch))));
                }
            }
        }
        while next.len() > 64 {
            let i = (rng.f() * next.len() as f64) as usize;
            next.swap_remove(i);
        }
        cells = next;
    }
    lens.sort();
    let median = lens[lens.len() / 2];
    assert!(widest <= 1.0 + 1e-9 && median <= 32, "widest {widest} of a cell's diagonal, median {median} instructions");
}

/// A double-double number (hi + lo, about 32 digits): enough to see a
/// slope past 1 by an ulp.
#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct Dd(f64, f64);

impl Dd {
    fn k(v: f64) -> Dd {
        Dd(v, 0.0)
    }
    fn norm(hi: f64, lo: f64) -> Dd {
        let s = hi + lo;
        Dd(s, lo - (s - hi))
    }
    fn add(self, o: Dd) -> Dd {
        let s = self.0 + o.0;
        let bb = s - self.0;
        let e = (self.0 - (s - bb)) + (o.0 - bb);
        Dd::norm(s, e + self.1 + o.1)
    }
    fn neg(self) -> Dd {
        Dd(-self.0, -self.1)
    }
    fn mul(self, o: Dd) -> Dd {
        let p = self.0 * o.0;
        Dd::norm(p, self.0.mul_add(o.0, -p) + self.0 * o.1 + self.1 * o.0)
    }
    fn div(self, o: Dd) -> Dd {
        let q = self.0 / o.0;
        let r = self.add(o.mul(Dd::k(q)).neg());
        Dd::norm(q, r.0 / o.0)
    }
    fn sqrt(self) -> Dd {
        if self.0 <= 0.0 {
            return Dd::k(0.0);
        }
        let x = self.0.sqrt();
        let r = self.add(Dd::k(x).mul(Dd::k(x)).neg());
        Dd::norm(x, r.0 / (2.0 * x))
    }
}

/// The value and gradient of a tape at p, in double-double arithmetic.
fn dd_grad(t: &crate::tape::Tape, p: [f64; 3]) -> (Dd, [Dd; 3]) {
    use crate::graph::{BinOp, UnOp, Var};
    use crate::tape::Op;
    let z = Dd::k(0.0);
    let mut v: Vec<(Dd, [Dd; 3])> = Vec::with_capacity(t.ssa.len());
    for ins in &t.ssa {
        let (a, b) = (ins.a as usize, ins.b as usize);
        let lin = |f: &dyn Fn(Dd, Dd) -> Dd, (x, dx): (Dd, [Dd; 3]), (y, dy): (Dd, [Dd; 3])| (f(x, y), [0, 1, 2].map(|k| f(dx[k], dy[k])));
        let r = match ins.op {
            Op::Const(c) => (Dd::k(c), [z; 3]),
            // A number known within bounds: its point value, as a point reads.
            Op::Within(_) => (v[ins.c as usize].0, [z; 3]),
            Op::Var(Var::X | Var::Y | Var::Z) => {
                let j = match ins.op {
                    Op::Var(Var::X) => 0,
                    Op::Var(Var::Y) => 1,
                    _ => 2,
                };
                (Dd::k(p[j]), [0, 1, 2].map(|k| Dd::k(if k == j { 1.0 } else { 0.0 })))
            }
            Op::Tag(_) | Op::Lipschitz | Op::Real => v[a],
            Op::Frame => (z, [z; 3]),
            Op::Bin(BinOp::Add) => lin(&|x, y| x.add(y), v[a], v[b]),
            Op::Bin(BinOp::Sub) => lin(&|x, y| x.add(y.neg()), v[a], v[b]),
            Op::Bin(BinOp::Mul) => {
                let ((x, dx), (y, dy)) = (v[a], v[b]);
                (x.mul(y), [0, 1, 2].map(|k| dx[k].mul(y).add(dy[k].mul(x))))
            }
            Op::Bin(BinOp::Div) => {
                let ((x, dx), (y, dy)) = (v[a], v[b]);
                let q = x.div(y);
                (q, [0, 1, 2].map(|k| dx[k].add(dy[k].mul(q).neg()).div(y)))
            }
            Op::Bin(BinOp::Min) => {
                if v[a].0 <= v[b].0 {
                    v[a]
                } else {
                    v[b]
                }
            }
            Op::Bin(BinOp::Max) => {
                if v[a].0 >= v[b].0 {
                    v[a]
                } else {
                    v[b]
                }
            }
            Op::Un(UnOp::Neg) => (v[a].0.neg(), v[a].1.map(Dd::neg)),
            Op::Un(UnOp::Square) => {
                let (x, dx) = v[a];
                (x.mul(x), dx.map(|d| d.mul(x).mul(Dd::k(2.0))))
            }
            Op::Un(UnOp::Sqrt) => {
                let (x, dx) = v[a];
                let r = x.sqrt();
                (r, dx.map(|d| if r.0 > 0.0 { d.div(r.mul(Dd::k(2.0))) } else { z }))
            }
            op => panic!("{op:?} in a drafted field"),
        };
        v.push(r);
    }
    v[v.len() - 1]
}

/// Each promised leaf of the drafted field `src`, read in real numbers from
/// its constants as computed: its worst slope past its promise, and where.
fn promise_excess(src: &str) -> Option<(Dd, [f64; 3])> {
    let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).ok()?;
    let mut worst: Option<(Dd, [f64; 3])> = None;
    for node in &doc.g.nodes {
        let crate::graph::Node::Lipschitz(a, k, _) = *node else { continue };
        let k = doc.g.as_const(k).expect("a drafted promise's constant");
        let t = crate::tape::Tape::new(&doc.g, a);
        for x in [-10.0, -1.0, -0.3, 0.0, 0.004, 1.0] {
            for y in [-1.0, -0.004, 0.0, 1.0] {
                for z in [-0.02, -0.01, -0.005, 0.0, 0.005, 0.01, 0.02] {
                    let (_, g) = dd_grad(&t, [x, y, z]);
                    let over = g[0].mul(g[0]).add(g[1].mul(g[1])).add(g[2].mul(g[2])).sqrt().add(Dd::k(-k));
                    if worst.is_none_or(|w| over > w.0) {
                        worst = Some((over, [x, y, z]));
                    }
                }
            }
        }
    }
    worst
}

#[test]
fn a_drafted_profiles_promises_hold_in_real_numbers() {
    // Each promised leaf of a drafted profile, read in real numbers from its
    // constants as computed, has a slope of at most its promise: its scale
    // from its rate's bound rounded up, and down itself (rounded to nearest,
    // each of these went over by 1e-17 to 1.4e-16). Steps by a small
    // concave arc, drafted alike or not.
    for src in [
        "(polygon (vec (vec -96.70866030515135 -0.41653411565695303) (vec -96.70866030515135 96.70866030515135) (vec 96.70866030515135 96.70866030515135) (vec 96.70866030515135 -0.8330682313139061) (vec 0.0 -0.8330682313139061 -0.2939028637685043) (vec -2.5553740565951535 -0.41653411565695303)) (vec -2.6115328216160507 0 -2.6115328216160507 -2.6115328216160507 -2.6115328216160507 -2.6115328216160507) -0.1277687028297577 0.1277687028297577)",
        "(polygon (vec (vec -112.21413049019851 -0.12157340354488907) (vec -112.21413049019851 112.21413049019851) (vec 112.21413049019851 112.21413049019851) (vec 112.21413049019851 -0.24314680708977815) (vec 0.0 -0.24314680708977815 -0.23441377092754592) (vec -0.17593704088099904 -0.12157340354488907)) -0.9288834944677701 -0.008796852044049952 0.008796852044049952)",
        "(polygon (vec (vec -0.016565269827376064 -3.174596157293424e-06) (vec -0.016565269827376064 0.016565269827376064) (vec 0.016565269827376064 0.016565269827376064) (vec 0.016565269827376064 -6.349192314586848e-06) (vec 0.0 -6.349192314586848e-06 -0.5341988189136128) (vec -2.799769430685777e-05 -3.174596157293424e-06)) -18.48080005810317 -1.3998847153428887e-06 1.3998847153428887e-06)",
    ] {
        let (over, at) = promise_excess(src).expect("promised leaves");
        // Double-double's own rounding is under 1e-28.
        assert!(over.0 <= 1e-28, "a promise exceeded by {} at {at:?}", over.0);
    }
}

#[test]
fn a_drafted_profiles_planes_and_cones_slope_at_most_1() {
    // A drafted line wall's plane, scaled coefficient by coefficient by 1 /
    // hypot(1, its slope), read 11258999068426241 / 11258999068426240 at
    // the slanted edge of (0, 0) (3, 4) (0, 4) drafted 0 over 0..1; each
    // plane is scaled once from all its coefficients, rounded down. A cone
    // drafted 1e-8 degrees, its rate's hypot rounded to 1, was promised
    // unscaled: its scale is rounded down too.
    let doc = Doc::load_str_in(&CORE, "(def f (polygon (vec (vec 0 0) (vec 3 4) (vec 0 4)) 0 0 1))", "test", std::path::Path::new("."), &[]).unwrap();
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let t = crate::tape::Tape::new(&doc.g, fid);
    for x in [-1.0, 0.5, 1.5, 2.0, 3.5] {
        for y in [-1.0, 1.0, 2.25, 3.0, 4.5] {
            for z in [0.25, 0.5, 0.75] {
                let (_, g) = dd_grad(&t, [x, y, z]);
                let slope = g[0].mul(g[0]).add(g[1].mul(g[1])).add(g[2].mul(g[2])).sqrt();
                // Double-double's own rounding is under 1e-28.
                assert!(slope <= Dd::k(1.0 + 1e-28), "sloped {} + {} at {:?}", slope.0, slope.1, [x, y, z]);
            }
        }
    }
    let (over, at) = promise_excess("(polygon (vec (vec -4 0 1) (vec 4 0 1)) 0.00000001 -1 1)").expect("promised leaves");
    assert!(over.0 <= 1e-28, "a promise exceeded by {} at {at:?}", over.0);
}

#[test]
fn a_flat_arc_is_its_circle() {
    // A square 3.1e-60 across, its bottom an arc of
    // bulge 2.8e-9 (sag 4.4e-69). A point 3.2e-70 inside, below the chord
    // and above the arc, read outside (polygon's float disc); an arc is its
    // circle exactly, so it reads inside.
    let s = 3.111507638930571e-60;
    let src = format!("(def f (polygon (vec (vec 0 0 2.8e-09) (vec {s:?} 0) (vec {s:?} {s:?}) (vec 0 {s:?}))))");
    let doc = Doc::load_str_in(&CORE, &src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let p = |v: f64| ink_core::interval::I::point(v);
    let v = crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [p(3.111508012991262e-61), p(-1.2446030555722282e-69), p(0.0)], None);
    assert!(v.hi < 0.0, "{v:?}");
}

#[test]
fn an_arc_flatter_than_a_float_step_is_its_circle() {
    // A 10 mm square turned 0.1 rad, 1e12 from the origin,
    // its bottom an arc of bulge 2e-8 (sag 1e-7, far below a float's step
    // in its coordinates' scale). A point 3e-6 outside its bottom read
    // 7.4e-6 inside (polygon's allowance); its kite's corner is found
    // along the chord's normal, so it reads outside.
    let src = "(def f (polygon (vec (vec 1000000000000 0 2e-08) (vec 1000000000009.9501 0.99833416646828155 0) (vec 1000000000008.9518 10.94837581924854 0) (vec 999999999999.00171 9.9500416527802589 0))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let p = |v: f64| ink_core::interval::I::point(v);
    let v = crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [p(1000000000002.4875), p(0.24958353166699715), p(0.0)], None);
    assert!(v.lo > 0.0, "{v:?}");
}

#[test]
fn an_arc_far_from_the_origin_is_its_circle() {
    // A 1 mm square 2e6 mm out, its bottom an arc of bulge 0.3. Its circle
    // is read about the arc's end, its centre and radius each known within
    // bounds on their exact values: a point 1.7e-10 inside the arc reads
    // inside (polygon's construction read the square empty: 1.07e-3).
    let src = "(def f (polygon (vec (vec 1048576.3 2097152.7 0.3) (vec 1048576.301 2097152.7) (vec 1048576.301 2097152.701) (vec 1048576.3 2097152.701))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let p = |v: f64| ink_core::interval::I::point(v);
    let v = crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [p(1048576.3008818447), p(2097152.699934159), p(0.0)], None);
    assert!(v.hi < 0.0 && v.lo > -1.7e-10, "{v:?}");
}

#[test]
fn a_round_is_its_circle() {
    // A triangle 1e5 across, each corner rounded 0.02 (the one at
    // (1e5, 0) turns 168 degrees). A point 9.5e-12 outside that round read
    // inside (its tangent points solved in floats, polygon's allowance);
    // the corners' polygon is exact and the round's centre known within
    // bounds, so it reads outside.
    let src = "(def f (polygon (vec (vec 0 0 0 0.02) (vec 100000 0 0 0.02) (vec 20000 18000 0 0.02))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let p = |v: f64| ink_core::interval::I::point(v);
    let v = crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [p(99999.83985368833), p(0.017585241216589442), p(0.0)], None);
    assert!(v.lo > 9.5e-12 && v.hi < 9.53e-12, "{v:?}");
}

#[test]
fn a_round_beside_an_arc_is_its_circle() {
    // A 10 mm square 2048 mm out, its bottom an arc of bulge 1e-9, each
    // corner rounded 0.2: the bottom corners' rounds meet the arc. A point
    // 1e-8 outside the bottom right round read inside (its tangent points
    // solved in floats); its centre is bounded where its offset line and
    // circle meet, so it reads outside.
    let src = "(def f (polygon (vec (vec 2048 0 1e-09 0.2) (vec 2058 0 0 0.2) (vec 2058 10 0 0.2) (vec 2048 10 0 0.2))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let p = |v: f64| ink_core::interval::I::point(v);
    let v = crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [p(2057.997537678002), p(0.16871310507353643), p(0.0)], None);
    assert!(v.lo > 9.9e-9 && v.hi < 1.01e-8, "{v:?}");
}

#[test]
fn a_filled_round_keeps_its_profile_in_its_box() {
    // An L whose inner corner is rounded (a fill, which adds material):
    // the field is positive over a large box beside it, so the body is
    // proven bounded.
    let src = "(def f (polygon (vec (vec 0 0) (vec 2 0) (vec 2 1) (vec 1 1 0 0.3) (vec 1 2) (vec 0 2))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let ctx = crate::eval::Ctx::new(0.0, vec![]);
    let i = ink_core::interval::I::new;
    for b in [[i(3.0, 1e6), i(-1e6, 1e6)], [i(-1e6, 1e6), i(3.0, 1e6)], [i(-1e6, -1.0), i(-1e6, 1e6)]] {
        let v = crate::eval::eval_interval(&tape, &ctx, [b[0], b[1], i(0.0, 0.0)], None);
        assert!(v.lo > 0.0, "{b:?}: {v:?}");
    }
    let v = crate::eval::eval_point(&tape, &ctx, [1.05, 1.05, 0.0]);
    assert!(v < 0.0, "the fill: {v}");
}

#[test]
fn a_fill_is_held_to_its_corner() {
    // A leadscrew thread's section, its two inner corners rounded 0.12
    // (fills). Each fill's triangle must be proven apart from the other
    // edges by its corners' hull in order: read as the corners in the
    // order found, one fill's triangle took in the other flank, and
    // (0.98, -0.12), inside the other fill's circle (outside the profile),
    // read inside.
    let src = "(def f (polygon (vec (vec 0.8666265877365439 -0.02900635094611581) (vec 1.0 -0.06474364905389207 0 0.12) (vec 1.0 -0.1875) (vec 2.0 -0.1875) (vec 2.0 0.1875) (vec 1.0 0.1875) (vec 1.0 0.06474364905389207 0 0.12) (vec 0.8666265877365439 0.02900635094611581))))";
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let ctx = crate::eval::Ctx::new(0.0, vec![]);
    for (p, outside) in [([0.98, -0.12], true), ([0.98, 0.12], true), ([0.99, -0.08], false), ([0.96, 0.0], false)] {
        let v = crate::eval::eval_point(&tape, &ctx, [p[0], p[1], 0.0]);
        assert_eq!(v > 0.0, outside, "{p:?}: {v}");
    }
}

fn interval_at(src: &str, p: [f64; 2]) -> ink_core::interval::I {
    let doc = Doc::load_str_in(&CORE, src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let i = ink_core::interval::I::point;
    crate::eval::eval_interval(&tape, &crate::eval::Ctx::new(0.0, vec![]), [i(p[0]), i(p[1]), i(0.0)], None)
}

#[test]
fn a_round_s_triangle_is_bounded() {
    // A round 0.006 beside an arc of bulge 0.37, at a sharp corner: the
    // line cutting its wedge beyond the corner crosses one of the wedge's
    // rays behind its centre, so the wedge runs on past it, across the
    // whole profile. Its triangle's bound (the hull of where the normals
    // meet that line) held none of that, and (1, 0.68), deep inside, read
    // outside.
    let src = "(def f (polygon (vec (vec 0.023 -0.054 0.37) (vec 0 0 0 0.006) (vec 0.55 -1.56) (vec 2 -2) (vec 2 2) (vec 0.1 2) (vec 0.1 -0.1))))";
    let v = interval_at(src, [1.0, 0.68]);
    assert!(v.hi < -0.8, "{v:?}");
}

#[test]
fn a_fill_s_box_holds_every_arc() {
    // An L 1e8 out, its bottom an arc of bulge 2e-9 (2e-9 below its chord,
    // under an ulp of y there), its inner corner a fill: the profile's box
    // must hold the arc, and the chord's middle, 2e-9 inside, reads inside.
    let src = "(def f (polygon (vec (vec 0 100000000 2e-9) (vec 2 100000000) (vec 2 100000001) (vec 1 100000001 0 0.3) (vec 1 100000002) (vec 0 100000002))))";
    let v = interval_at(src, [1.0, 1e8]);
    assert!(v.hi < 0.0, "{v:?}");
}

#[test]
fn a_round_s_hold_is_rounded_out_at_its_exact_corners() {
    // A round 0.25 at (0, 1e8) between square edges: its centre less the
    // corner, (0.25, 0.25), is exact, and its triangle reaches
    // y = 1e8 + 0.25. Its bound took both ends of each exact coordinate
    // down, so it stopped an ulp short, and a tab whose lower edge lies
    // there passed as apart: (-0.5, 1e8 + 0.25), inside the tab, read 0.
    let src = "(def f (polygon (vec (vec 0 100000002.0) (vec 0 100000000.0 0 0.25) (vec 2 100000000.0) (vec 2 100000003.0) (vec -1 100000003.0) (vec -1 100000000.24999999) (vec -0.4 100000000.24999999) (vec -0.4 100000000.25000001) (vec -0.9 100000000.25000001) (vec -0.9 100000002.0))))";
    let v = interval_at(src, [-0.5, 100000000.25]);
    assert!(v.hi < 0.0, "{v:?}");
}

#[test]
fn two_arcs_whose_circles_round_alike_are_two_circles() {
    // Two unit squares whose bottom arcs' bulges differ by 2^-55: their
    // circles' bounds round to the same floats, but they are two exact
    // circles, so their difference at (0.5, 0), -2^-56 there, is not 0.
    let src = "(def f (- (polygon (vec (vec 0 0 0.24876511400000006) (vec 1 0) (vec 1 1) (vec 0 1))) (polygon (vec (vec 0 0 0.24876511400000004) (vec 1 0) (vec 1 1) (vec 0 1)))))";
    let v = interval_at(src, [0.5, 0.0]);
    let truth = -(0.24876511400000006f64 - 0.24876511400000004) / 2.0;
    assert!(v.lo <= truth && truth <= v.hi, "{v:?} {truth}");
}

/// Two numbers +1/2 and -1/2 exactly, from two unit squares whose arcs'
/// bulges differ by 2^-55 (their circles' bounds the same floats): each
/// square reads -b/2 at (0.5, 0), scaled up by 2^56 about the mean.
const TWO_ARC_NUMBERS: &str = "(def a (+ 0.5 (* 72057594037927936 (+ (warp (polygon (vec (vec 0 0 0.24876511400000004) (vec 1 0) (vec 1 1) (vec 0 1))) 0.5 0 0) 0.12438255700000002)))) (def b (+ 0.5 (* 72057594037927936 (+ (warp (polygon (vec (vec 0 0 0.24876511400000006) (vec 1 0) (vec 1 1) (vec 0 1))) 0.5 0 0) 0.12438255700000002)))) (def c (select a -1 1)) (def d (select b -1 1))";

#[test]
fn one_arc_number_twice_and_two_arc_numbers_key_apart() {
    // c c is a square (at least 0), c d is not (it is -1): a key that read
    // their bounds alone gave both one key, so a grid baked of one was
    // read with the other's error bound.
    let key = |f: &str| {
        let doc = Doc::load_str_in(&CORE, &format!("{TWO_ARC_NUMBERS} (def f {f})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
        let mut fp = crate::graph::Fingerprint::default();
        doc.g.fingerprint(&[fid], &[], false, &mut fp);
        fp.value()
    };
    assert_ne!(key("(* c c)"), key("(* c d)"));
}

#[test]
fn a_select_between_two_arc_numbers_may_switch() {
    // select(x - y, x + a, x + b), a = 1/2 and b = -1/2: along y = z = 0
    // it changes sign at -1/2, at its jump at 0, and at 1/2. Its branches
    // read as one expression (a and b alike by their bounds), so it was
    // taken as continuous there, one sheet.
    let doc =
        Doc::load_str_in(&CORE, &format!("{TWO_ARC_NUMBERS} (def f (select (- x y) (+ x a) (+ x b)))"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let t = crate::tape::Tape::new(&doc.g, fid);
    let region = [ink_core::interval::I::new(-2.0, 2.0), ink_core::interval::I::point(0.0), ink_core::interval::I::point(0.0)];
    let mut choices = vec![crate::tape::BOTH; t.len()];
    crate::eval::eval_interval(&t, &crate::eval::Ctx::new(0.0, vec![]), region, Some(&mut choices));
    assert!(t.switches(&choices));
}

#[test]
fn a_select_of_an_arc_number_and_its_double_is_continuous() {
    // select(h - j, h, j) with j = 2 h is min(h, j), continuous whatever
    // the arc's number a in h = x + a: a number known within bounds is
    // one constant, so the select is one field, no seam.
    let src = format!(
        "{TWO_ARC_NUMBERS} (def h (+ x (warp (polygon (vec (vec 0 0 0.24876511400000004) (vec 1 0) (vec 1 1) (vec 0 1))) 0.5 0 0))) (def j (* 2 h)) (def f (select (- h j) h j))"
    );
    let doc = Doc::load_str_in(&CORE, &src, "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let t = crate::tape::Tape::new(&doc.g, fid);
    let region = [ink_core::interval::I::new(-1.0, 1.0), ink_core::interval::I::point(0.0), ink_core::interval::I::point(0.0)];
    let mut choices = vec![crate::tape::BOTH; t.len()];
    crate::eval::eval_interval(&t, &crate::eval::Ctx::new(0.0, vec![]), region, Some(&mut choices));
    assert!(!t.switches(&choices));
}

#[test]
fn an_arc_past_a_hull_line_over_a_pocket_keeps_its_segment() {
    // A hole of SHERPA-MINI's core: its 143 degree arc's segment pokes past
    // a hull line with a pocket under it. That line giving way to the
    // segment let the pocket (which reads past its own lid) take part of
    // the segment: (-0.34, 5.39), 1.1 inside, read 0.34 outside.
    let pts = "(vec (vec -6.5 -0.0 -0.71552) (vec 5.1456 3.9715 0.37636) (vec 6.4396 3.7179 -0.08405) (vec 7.7838 4.2218 0.0) (vec 7.8162 4.23 0.0) (vec 23.2851 9.0561 -0.50583) (vec 24.5829 8.1015 -0.03748) (vec 24.5717 7.9522 0.0) (vec 22.5925 -5.1567 -0.37098) (vec 21.4962 -6.1 0.0) (vec 9.3 -6.1 -0.41421) (vec 7.8 -4.6 0.0) (vec 7.8 -4.4665 0.34513) (vec 7.5717 -4.1752 -0.07135) (vec 6.4396 -3.7179 0.37636) (vec 5.1456 -3.9715 -0.71552))";
    let doc = Doc::load_str_in(&CORE, &format!("(def f (polygon {pts}))"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = crate::tape::Tape::new(&doc.g, fid);
    let v = crate::eval::eval_point(&tape, &crate::eval::Ctx::new(0.0, vec![]), [-0.3396011352538846, 5.392036247253436, 0.0]);
    assert!((v + 1.0973).abs() < 1e-3, "{v}");
}

#[test]
fn a_mesh_of_a_drafted_chamfer_lies_on_its_walls() {
    // A uniform mesh of MPM's plate where its chamfer meets the undrafted
    // fillet of its tab: each vertex, placed by the field's values at its
    // cell's corners, is within a sixth of a cell of the walls (a field read
    // at a thousandth of the distance inside and in full outside would set
    // vertices up to a cell inside them).
    let (p, d, z0, z1) = mpm_plate();
    let ds: Vec<String> = d.iter().map(|d| d.to_string()).collect();
    let src = format!("(max (polygon {} (vec {}) {z0} {z1}) (- (abs z) {z1}))", ink_points(&p), ds.join(" "));
    let doc = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
    let Value::Field(fid, _) = doc.def("f").unwrap().clone() else { panic!() };
    let tape = std::sync::Arc::new(crate::tape::Tape::new(&doc.g, fid));
    let h = 0.05;
    let m = ink_kernel::mesh::mesh(&tape, &crate::eval::Ctx::new(0.0, vec![]), [-8.0, -52.5, -0.3], [-5.0, -50.5, 0.3], h);
    let w = walls::Walls::new(&p, &d.iter().map(|d| d.to_radians().tan()).collect::<Vec<_>>()).unwrap();
    let dist = wall_distance(&w, z0, z1);
    assert!(m.verts.len() > 1000, "{} vertices", m.verts.len());
    let worst = m.verts.iter().map(|&[x, y, z]| dist(x, y, z).min((z.abs() - z1).abs())).fold(0.0, f64::max);
    assert!(worst <= h / 6.0, "a vertex {worst} from the plate's surface, at cells of {h}");
}

#[test]
fn a_drafted_wall_has_its_normal_on_it() {
    // WCP's ball-bearing block chamfers its edges with a drafted polygon
    // whose undrafted walls are axis-aligned: a mesh's crossings land exactly
    // on them. There the field is 0 and its slope is across the wall (a
    // distance along the wall, unsigned, is 0 there with no slope, and the
    // vertices lose the wall's normal: 0.42 off at cells of 0.45).
    let src = "(polygon (vec (vec 23.7109 66.675) (vec 3.9751 66.675) (vec 3.9751 69.85) (vec 23.7109 69.85 -0.414214) (vec 28.1051 65.4558) (vec 28.1051 34.925) (vec 24.9301 34.925) (vec 24.9301 65.4558 0.414214)) (vec 0 -45 0 0 0 -45 0 0) -3.175 3.175)";
    let d = Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).unwrap();
    let Value::Field(fid, _) = d.def("f").unwrap().clone() else { panic!() };
    let t = crate::tape::Tape::new(&d.g, fid);
    let pts: Vec<[f64; 3]> = [10.0, 15.115, 20.0].iter().flat_map(|&x| [-2.5, 0.0, 2.667534].map(|z| [x, 66.675, z])).collect();
    for (p, g) in pts.iter().zip(crate::eval::eval_grads(&t, &crate::eval::Ctx::new(0.0, vec![]), &pts)) {
        assert!(g[0] == 0.0 && g[2] < -0.1 && g[1].abs() < 1e-9 && g[3].abs() < 1e-9, "on the wall at {p:?}: {g:?}");
    }
}

/// Why polygon refuses `src`, if it does.
fn refused(src: &str) -> Option<String> {
    Doc::load_str_in(&CORE, &format!("(def f {src})"), "test", std::path::Path::new("."), &[]).err().map(|e| e.message)
}

#[test]
fn a_drafted_profile_keeps_motion_small_beside_its_coordinates() {
    // A rounded square 1e7 from the origin, drafted 45 degrees over 8e-6:
    // its walls move 8e-6, 570 times its allowance (100 eps, 1.4e-6) and
    // under the 1e-5 its paths read as rounding beside a coordinate of 1e7.
    // Frozen, its rounds would stand still and its box clip its right wall
    // at z = 0. At z = 6e-6 the right wall is at x0 + 10 + 6e-6, the lower
    // right round's centre at (x0 + 9 + z, 1 - z).
    let x0 = 1e7;
    let src = format!("(polygon (vec (vec {x0} 0 0 1) (vec {} 0 0 1) (vec {} 10 0 1) (vec {x0} 10 0 1)) 45 0 0.000008)", x0 + 10.0, x0 + 10.0);
    if let Some(e) = refused(&src) {
        panic!("refused: {e}");
    }
    let f = field(&src);
    let z = 6e-6;
    assert!(f([x0 + 10.0 + 3e-6, 5.0, z]) < 0.0, "3e-6 inside the moving right wall");
    assert!(f([x0 + 10.0 + 9e-6, 5.0, z]) > 0.0, "3e-6 outside the moving right wall");
    // Along the round's diagonal, 4e-6 inside it (outside the frozen one).
    let (c, u) = ([x0 + 9.0 + z, 1.0 - z], std::f64::consts::FRAC_1_SQRT_2);
    assert!(f([c[0] + u * (1.0 - 4e-6), c[1] - u * (1.0 - 4e-6), z]) < 0.0, "4e-6 inside the moving round");
    assert!(f([c[0] + u * (1.0 + 4e-6), c[1] - u * (1.0 + 4e-6), z]) > 0.0, "4e-6 outside the moving round");
}

#[test]
fn a_drafted_shallow_arc_keeps_its_radius_growing() {
    // An arc of radius 1e8 (sagitta 1.25e-7) drafted 45 degrees over 5e-5:
    // its radius grows by 5e-5, 35 times the allowance (1.4e-6) and under
    // the 1e-4 its path reads as rounding beside the radius. Frozen, the
    // arc would stand still while its corners moved on. At z = 2.5e-5 its
    // bottom at x = 0 is at y = -1.25e-7 - z (the probe 1.26e-5 inside it is
    // inside, as a claim of accuracy 0).
    let src = "(polygon (vec (vec -5 0 0.000000025) (vec 5 0) (vec 5 10) (vec -5 10)) (vec 45 0 0 0) 0 0.00005)";
    if let Some(e) = refused(src) {
        panic!("refused: {e}");
    }
    let f = field(src);
    let (z, bottom) = (2.5e-5, -1.25e-7 - 2.5e-5);
    assert!(f([0.0, -1.25e-5, z]) < 0.0, "1.26e-5 inside the drafted arc");
    assert!(f([0.0, bottom + 5e-6, z]) < 0.0, "5e-6 inside the drafted arc");
    assert!(f([0.0, bottom - 5e-6, z]) > 0.0, "5e-6 outside the drafted arc");
}

#[test]
fn a_drafted_profile_far_from_the_origin_is_read() {
    // 1e12 from the origin. A rounded square drafted 45
    // degrees over 0.5: with its motion read as rounding beside 1e12 its
    // rounds would stand still while its walls moved, and the point 0.125
    // inside its lower right round at z = 0.25 (centre (x0 + 9.25, 0.75),
    // radius 1) read outside, as a claim of accuracy 0. Its walls are axis-aligned, so
    // its paths are exact: it is read right.
    let x0 = 1e12;
    let f = field("(polygon (vec (vec 1000000000000 0 0 1) (vec 1000000000010 0 0 1) (vec 1000000000010 10 0 1) (vec 1000000000000 10 0 1)) 45 0 0.5)");
    assert!(f([x0 + 10.125, 0.75, 0.25]) < 0.0, "0.125 inside the moving round");
    assert!(f([x0 + 10.375, 0.75, 0.25]) > 0.0, "0.125 outside the moving round");
    // A strip whose sides lean in to meet at z = 0.25: frozen, it would be
    // read as a stationary rectangle. It closes there, empty beyond, as it
    // does moved to the origin (the same geometry: 0.5 is exact).
    for x in [1e12, 0.0] {
        let src = format!("(polygon (vec (vec {x} 0) (vec {} 0) (vec {} 10) (vec {x} 10)) (vec 0 -45 0 -45) 0 0.5)", x + 0.5, x + 0.5);
        let f = field(&src);
        assert!(f([x + 0.25, 5.0, 0.24]) < 0.0 && f([x + 0.25, 5.0, 0.26]) > 0.0 && f([x + 0.4, 5.0, 0.2]) > 0.0, "{src}");
    }
    // The square turned (its sides along 3-4-5 triangles), so its walls'
    // data round: 1e10 from the origin a coordinate's rounding (1.9e-6) is
    // the allowance's 2e-6, and its paths could not be shown on its walls
    // within it. Built exactly (each wall the plane through its drawn edge,
    // each corner and round moving with its walls' bounds), it is read
    // right: 0.01 either side of its first wall, offset 0.25 at z = 0.25.
    for x0 in [1e9, 1e10, 1e12] {
        let p = |dx: f64, dy: f64| format!("(vec {} {dy} 0 1)", x0 + dx);
        let f = field(&format!("(polygon (vec {} {} {} {}) 45 0 0.5)", p(0.0, 0.0), p(8.0, 6.0), p(2.0, 14.0), p(-6.0, 8.0)));
        let at = |d: f64| f([x0 + 4.0 + 0.6 * d, 3.0 - 0.8 * d, 0.25]);
        assert!(at(0.24) < 0.0 && at(0.26) > 0.0, "{x0}: {} {}", at(0.24), at(0.26));
        assert!(f([x0 + 1.0, 7.0, 0.25]) < -3.0, "{x0}");
    }
}

#[test]
fn a_drafted_polygon_s_box_is_outward() {
    // The box a drafted polygon is cut by holds its emitted walls over the
    // span, grown by the allowance, rounded outward: its right side past
    // the real 10 + 0.1 tan(45) (above the float 10.1).
    let w = walls::Walls::new(&[[0.0, 0.0, 0.0, 0.0], [10.0, 0.0, 0.0, 0.0], [10.0, 10.0, 0.0, 0.0], [0.0, 10.0, 0.0, 0.0]], &[0.0, 45f64.to_radians().tan(), 0.0, 0.0]).unwrap();
    let d = w.drafted(0.0, 0.1).unwrap();
    assert!(d.hi[0] > 10.1 && d.hi[0] < 10.1 + 1e-5 && d.lo[0] < 0.0 && d.lo[1] < 0.0 && d.hi[1] > 10.0, "{:?} {:?}", d.lo, d.hi);
}

#[test]
fn a_contact_s_identities_are_taken_only_where_lowering_keeps_them() {
    // The graph folds an operation on two constants in floats as it is
    // lowered; `touch` folds so too, and says so: its tangency is an
    // identity of what the field reads only where nothing folded.
    use curve::{kx, pk, xlin};
    let fixed = curve::touch(&pk([0.0, 0.0]), &kx(1.0), &pk([3.0, 1.0]), &kx(0.0), 1.0);
    assert!(fixed.1, "fixed circles fold");
    assert!(fixed.0.iter().all(|x| matches!(**x, curve::Ex::K(_))), "folded to constants: {:?}", fixed.0);
    let moving = curve::touch(&pk([0.0, 0.0]), &kx(1.0), &[xlin(3.0, 4.0), kx(1.0)], &kx(0.0), 1.0);
    assert!(!moving.1, "a moving point does not");
}

#[test]
fn a_polygon_that_crosses_itself_is_refused_where_it_does() {
    // An arc crossing the line beside it, near (0.08, 0): its chords do not
    // cross, and no corner is inside its segment.
    let pts = [[5.04004, -1.33967668, 0.132669958832], [0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [10.0, -5.0, 0.0], [5.04004, -5.0, 0.0]];
    let e = polygon(&pts).expect_err("refused");
    assert!(e.contains("crosses itself at (0.08"), "{e}");
    // A fillet's outside, as STEP data writes it: a line, then a quarter
    // arc leaving tangent to it, back along it. Its bulge rounded to 1e-6
    // lifts the arc over the line, meeting it again 3.6e-6 from the corner
    // across a sliver 1e-12 thick, far within the allowance. Accepted.
    let cusp = [[-8.045452, 2.4257, 0.0], [-5.619752, 2.4257, 0.414214], [-8.045452, 0.0, 0.0]];
    if let Err(e) = polygon(&cusp) {
        panic!("a rounded cusp: {e}");
    }
    // A quarter arc crossing the line beside it far out on a profile 1e6
    // across, whose allowance (100 eps) is 0.1 and then 0.21: its sliver
    // (the arc's sag over the chord to the crossing) is 0.1056 and 0.2401,
    // over the allowance though d^2 / 8 r reads 0.1 and 0.2113. Refused.
    let q = -(std::f64::consts::PI / 8.0).tan();
    for (pts, what) in [
        (vec![[1.0, 1.0, q], [0.0, 0.0, 0.0], [2.0, 1.0, 0.0], [1e6, 1.0, 0.0], [1e6, 100.0, 0.0], [1.0, 100.0, 0.0]], "sag 0.1056"),
        (
            vec![[1.0, 1.0, q], [0.0, 0.0, 0.0], [1.9758289399641862, 1.6900000000000002, 0.0], [2112500.0, 1.6900000000000002, 0.0], [2112500.0, 100.0, 0.0], [1.0, 100.0, 0.0]],
            "sag 0.2401",
        ),
    ] {
        let e = polygon(&pts).expect_err(what);
        assert!(e.contains("crosses itself"), "{what}: {e}");
    }
    // Stars with arcs whose edges, densely sampled, cross: every one is
    // refused.
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut tried = 0;
    while tried < 300 {
        let n = 3 + (rng.f() * 12.0) as usize;
        let (r0, bulge, cw) = (rng.range(1.0, 9.0), rng.range(0.2, 1.5), rng.f() < 0.5);
        let pts = star(&mut rng, n, r0, 12.0, bulge, cw);
        if simple(&pts) {
            continue;
        }
        tried += 1;
        if let Ok(t) = polygon(&pts) {
            panic!("accepted a crossing polygon: {pts:?}: {:?}", ev(&t, 0.0, 0.0));
        }
    }
}

#[test]
fn a_drafted_profile_keeps_its_depth_where_its_construction_changes() {
    // A tooth in a pocket under a square's top edge, its flanks drafted 45
    // degrees: its tip rises through the edge's line at z = 0.0452, where
    // the hull changes (the tip joins it). Read as pieces cut by slabs, a
    // point 5 deep read -6e-5 there; one construction below the change and
    // the other above (the change's height known within bounds), it reads
    // its depth to the walls.
    let tooth = [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [7.0, 10.0], [7.0, 9.0], [5.0, 9.95], [3.0, 9.0], [3.0, 10.0], [0.0, 10.0]].map(|p| [p[0], p[1], 0.0, 0.0]);
    let drafts = [0.0, 0.0, 0.0, 0.0, 45.0, 45.0, 0.0, 0.0, 0.0];
    let ds: Vec<String> = drafts.iter().map(|d| d.to_string()).collect();
    let f = field(&format!("(polygon {} (vec {}) 0 0.2)", ink_points(&tooth), ds.join(" ")));
    for z in [0.04, 0.0451, 0.045163875, 0.0451639, 0.04517, 0.05] {
        assert!(f([5.0, 5.0, z]) < -0.4, "z = {z}: {}", f([5.0, 5.0, z]));
    }
    // Each piece's slope is at most 1 (pairs 1e-4 apart, not across the
    // change's height, where the pieces' values differ).
    let mut rng = Rng(7);
    for _ in 0..20000 {
        let p = [rng.range(-1.0, 11.0), rng.range(-1.0, 11.0), rng.range(0.0, 0.2)];
        let q = [p[0] + rng.range(-1e-4, 1e-4), p[1] + rng.range(-1e-4, 1e-4), p[2] + rng.range(-1e-4, 1e-4)];
        if (p[2] - 0.045164) * (q[2] - 0.045164) <= 0.0 {
            continue;
        }
        let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
        assert!((f(p) - f(q)).abs() <= d * (1.0 + 1e-9), "{p:?} {q:?}: {} {}", f(p), f(q));
    }
    drafted_signs(&tooth, &drafts, 0.0, 0.2, 40, "tooth").unwrap();
    drafted_signs(&tooth, &drafts, 0.0451, 0.0453, 20, "tooth at its change").unwrap();
}

#[test]
fn drafted_rounds_between_straight_walls_move_with_their_corners() {
    // A round between two straight walls keeps its radius as they lean, so
    // it moves with its corner: built exactly, its sign the section's at
    // every height. Rounded convex and reflex corners: a square, an L, a
    // notch, drafted out, in, and each wall its own way.
    let sq = [[0.0, 0.0, 0.0, 1.0], [10.0, 0.0, 0.0, 2.0], [10.0, 10.0, 0.0, 0.5], [0.0, 10.0, 0.0, 1.5]];
    let l = [[0.0, 0.0, 0.0, 0.5], [12.0, 0.0, 0.0, 1.0], [12.0, 4.0, 0.0, 0.5], [4.0, 4.0, 0.0, 1.0], [4.0, 10.0, 0.0, 0.5], [0.0, 10.0, 0.0, 0.5]];
    let notch = [
        [0.0, 0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0, 0.3],
        [4.0, 3.0, 0.0, 0.4],
        [6.0, 3.0, 0.0, 0.4],
        [6.0, 0.0, 0.0, 0.3],
        [10.0, 0.0, 0.0, 0.0],
        [10.0, 10.0, 0.0, 1.0],
        [0.0, 10.0, 0.0, 1.0],
    ];
    let cases: Vec<(&str, Vec<[f64; 4]>, Vec<f64>)> = vec![
        ("square out", sq.to_vec(), vec![10.0; 4]),
        ("square in", sq.to_vec(), vec![-10.0; 4]),
        ("L", l.to_vec(), vec![5.0, -3.0, 8.0, 0.0, 2.0, 12.0]),
        ("notch", notch.to_vec(), vec![0.0, 10.0, -20.0, 10.0, 0.0, 0.0, 3.0, 3.0]),
    ];
    for (name, p, d) in &cases {
        assert!(super::drafted::exactly(p, d, -1.0, 2.0), "{name}");
        drafted_signs(p, d, -1.0, 2.0, 24, name).unwrap();
    }
}

#[test]
fn a_drafted_edge_that_vanishes_is_trimmed() {
    // A square's corner chamfered, the chamfer drafted 45 degrees: it runs
    // out at the corner at z = 1 / sqrt 2. Past there its wall leaves and
    // its neighbours meet (a CAD trim), so the sharp corner is inside; it
    // was refused. Below z = 0 it leans in.
    let ch = [[0.0, 0.0], [10.0, 0.0], [10.0, 9.0], [9.0, 10.0], [0.0, 10.0]].map(|p| [p[0], p[1], 0.0, 0.0]);
    let d = [0.0, 0.0, 45.0, 0.0, 0.0];
    assert!(super::drafted::exactly(&ch, &d, -1.0, 1.0));
    let f = field(&format!("(polygon {} (vec 0 0 45 0 0) -1 1)", ink_points(&ch)));
    let at = |z: f64, s: f64| f([9.5 + s, 9.5 + s, z]);
    // The chamfer x + y = 19 + z sqrt 2 until it reaches the corner.
    for z in [-0.5, 0.0, 0.3, 0.7] {
        let s = 0.5 * z * std::f64::consts::SQRT_2;
        assert!(at(z, s - 0.01) < 0.0 && at(z, s + 0.01) > 0.0, "z = {z}");
    }
    for z in [0.71, 0.9, 1.0, 1.5] {
        assert!(f([9.99, 9.99, z]) < 0.0 && f([10.01, 9.99, z]) > 0.0, "z = {z}");
    }
    drafted_signs(&ch, &d, -1.0, 0.7, 16, "chamfer").unwrap();
}

#[test]
fn a_drafted_corner_turning_straight_and_an_edge_vanishing_are_read() {
    // Eight walls 1e6 from the origin, drafted -5 to 45 degrees: corner 3
    // turns from convex to reflex at z = 0.032 and edge 2 vanishes by
    // z = 0.064. Read by the allowance, a point 4.8e-5 inside the section
    // at z = 0.064 (by its walls solved in 50-digit decimals) read 0.041
    // outside.
    let p = [
        [1000000.888829729, -999999.2914697953],
        [1000000.7675199523, -999998.0240828756],
        [1000000.5257067526, -999998.4308712197],
        [1000000.1292991528, -999999.0041193124],
        [999999.8914277705, -999998.4010801787],
        [999998.111227277, -999999.8957132017],
        [999998.1988513363, -1000001.5895861739],
        [1000001.3794245238, -1000000.6128118325],
    ]
    .map(|q| [q[0], q[1], 0.0, 0.0]);
    let d = [-5.0, 0.0, 45.0, -20.0, 5.0, 0.0, 5.0, 0.0];
    assert!(super::drafted::exactly(&p, &d, 0.0, 0.14168026157146713));
    let ds: Vec<String> = d.iter().map(|v| v.to_string()).collect();
    let f = field(&format!("(polygon {} (vec {}) 0 0.14168026157146713)", ink_points(&p), ds.join(" ")));
    assert!(f([1000000.7578525974, -999997.982737937, 0.0639919310276438]) < 0.0);
    assert!(f([1000000.6606148512, -999998.0498147813, 0.10581535163270141]) < 0.0);
}

#[test]
fn a_drafted_arc_between_straight_walls_is_built_exactly() {
    // A rectangle whose bottom edge is an arc, bulging out or in, drafted
    // alike or each wall its own way: built exactly (its corners beside the
    // arc on curves, its segment its cone within a kite), its sign the
    // section's at every height.
    for (b, dr) in [(0.2, vec![10.0; 4]), (-0.2, vec![10.0; 4]), (0.2, vec![5.0, 0.0, 20.0, 0.0]), (0.414, vec![3.0; 4])] {
        let p = [[0.0, 0.0, b, 0.0], [10.0, 0.0, 0.0, 0.0], [10.0, 6.0, 0.0, 0.0], [0.0, 6.0, 0.0, 0.0]];
        assert!(super::drafted::exactly(&p, &dr, -1.0, 2.0), "{b} {dr:?}");
        drafted_signs(&p, &dr, -1.0, 2.0, 24, &format!("arc {b} {dr:?}")).unwrap();
    }
}

#[test]
fn drafted_arcs_meeting_each_other_are_built_exactly() {
    // Two arcs meeting at a corner (bulging alike or not), and a triangle of
    // arcs, drafted alike or each wall its own way: built exactly (each
    // corner of two arcs where their cones meet), its sign the section's at
    // every height.
    let two = |b1: f64, b2: f64| vec![[0.0, 0.0, b1, 0.0], [5.0, -0.5, b2, 0.0], [10.0, 0.0, 0.0, 0.0], [10.0, 6.0, 0.0, 0.0], [0.0, 6.0, 0.0, 0.0]];
    let cases = [
        (two(0.2, 0.15), vec![10.0; 5]),
        (two(0.2, -0.15), vec![10.0; 5]),
        (two(-0.2, -0.3), vec![4.0, 12.0, 0.0, 8.0, 2.0]),
        // Two upright arcs between leaning walls: their corner stays put.
        (two(0.2, 0.15), vec![0.0, 0.0, 10.0, 10.0, 10.0]),
        (vec![[0.0, 0.0, 0.2, 0.0], [10.0, 0.0, 0.2, 0.0], [5.0, 8.0, 0.2, 0.0]], vec![6.0; 3]),
        (vec![[0.0, 0.0, 0.2, 0.0], [10.0, 0.0, -0.1, 0.0], [5.0, 8.0, 0.3, 0.0]], vec![6.0, 2.0, 9.0]),
    ];
    for (p, dr) in cases {
        assert!(super::drafted::exactly(&p, &dr, -1.0, 2.0), "{p:?} {dr:?}");
        drafted_signs(&p, &dr, -1.0, 2.0, 24, &format!("arcs {p:?} {dr:?}")).unwrap();
    }
}

#[test]
fn a_drafted_arc_moving_past_its_own_sag_is_built_exactly() {
    // A slot's end (after a sprocket's hub, step2ink's numbers): one shallow
    // arc wall drafted 45 degrees over +-0.254, so it moves seven times its
    // sag. Its kite moves with it (scaled about the arc's centre).
    let p = [[9.236959, -0.79375, 0.0, 0.0], [7.897713, -0.79375, 0.050126, 0.0], [7.897713, 0.79375, 0.0, 0.0], [9.236959, 0.79375, -0.042887, 0.0]];
    for dr in [[0.0, 0.0, 0.0, 45.0], [0.0, 30.0, 0.0, 45.0], [10.0, 10.0, 10.0, 10.0]] {
        assert!(super::drafted::exactly(&p, &dr, -0.254, 0.254), "{dr:?}");
        drafted_signs(&p, &dr, -0.254, 0.254, 24, &format!("slot end {dr:?}")).unwrap();
    }
}

#[test]
fn a_drafted_arc_s_kite_holds_its_segment_at_every_height() {
    // Five arcs drafted each its own way (random, found by checking against
    // sections): kites of neighbouring arcs overlapped, and a point across
    // a chord brought in toward it at z = 0 was not across it at the span's
    // end (kite points move), so a point 0.9 inside read outside there.
    let p = [
        [3.817476960510768, 1.4132037021445196, 0.3988299480498855, 0.0],
        [0.2132427884221318, 4.46646300085881, 0.3183941415556502, 0.0],
        [-2.4973677977119135, 1.6286268469258767, 0.4700696716958046, 0.0],
        [-2.3981339470084455, -1.8865535235948039, -0.03400373223432407, 0.0],
        [1.8946033478533517, -3.3376178328806096, -0.13959255589737474, 0.0],
    ];
    let dr = [8.665801747981483, -8.165957131037231, 10.677595788435706, 10.519288390416126, -6.688330681847795];
    drafted_signs(&p, &dr, -1.142905790518473, 0.685698611180729, 24, "five arcs").unwrap();
}
