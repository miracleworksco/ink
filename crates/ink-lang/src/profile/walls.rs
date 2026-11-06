//! A profile as walls, for rounded corners and drafts: `(polygon pts)`
//! whose points carry a corner round, `(vec u v bulge round)`, and
//! `(polygon pts drafts z0 z1)`, whose walls lean (a mould's draft).
//!
//! A wall is an edge's line or circle. Drafted by angle a it is offset
//! outward by z tan a at height z: a plane, or a cone, tilted by a (a
//! sketch's arc is a wall, so a drafted extrude's rounds grow, as CAD's
//! do). A corner rounded to r is the circle of radius r tangent to both
//! its walls, its centre r inside both, solved again as they move: a
//! mould's round, which keeps its radius as the walls lean.
//!
//! A drafted profile is solved and built (`build`) at a height, and each
//! point of that construction's loop moves as its definition does, in
//! closed form (curve.rs): a corner where its walls meet, a round's ends
//! where it touches its walls, a halving at its direction from its arc's
//! centre, where a hull line touches an arc. `Moving` proves over how much
//! of the span the construction holds for that moving loop (each assumption
//! it recorded, `Ob`, the loop simple, its curves defined), to within
//! polygon's allowance for its data; where it stops holding, the profile is
//! built again there, and the pieces overlap (`Walls::drafted`). Each leaf
//! moves with the loop (`Doc::lower_profile`): a wall's half-plane is a
//! plane tilted by its draft, an arc wall's disc a cone, a round's disc its
//! circle about its moving centre; a line the construction adds (an arc's
//! chord, a pocket's lid) passes through its two moving points. So at every
//! height the field is the profile solved there.

use super::*;
use curve::{Ex, Sl, Tape, X, XP};
use ink_core::interval::I;

#[derive(Clone, Copy, Debug)]
enum Wall {
    /// n . p = k, n the outward unit normal.
    Line(P, f64),
    /// |p - c| = r, the profile inside the circle along this wall (s = 1)
    /// or outside it (s = -1).
    Circle(P, f64, f64),
}

impl Wall {
    /// Offset outward by d.
    fn moved(self, d: f64) -> Wall {
        match self {
            Wall::Line(n, k) => Wall::Line(n, k + d),
            Wall::Circle(c, r, s) => Wall::Circle(c, r + s * d, s),
        }
    }

    /// The point of the wall nearest p.
    fn foot(self, p: P) -> P {
        match self {
            Wall::Line(n, k) => sub(p, mul(n, dot(n, p) - k)),
            Wall::Circle(c, r, _) => add(c, mul(unit(sub(p, c)), r)),
        }
    }
}

/// Where two walls meet, nearest `near` (walls that coincide: near's foot).
fn meet(a: Wall, b: Wall, near: P, eps: f64) -> Option<P> {
    let pick = |c: [P; 2]| if norm(sub(c[0], near)) <= norm(sub(c[1], near)) { c[0] } else { c[1] };
    match (a, b) {
        (Wall::Line(n1, k1), Wall::Line(n2, k2)) => {
            let det = cross(n1, n2);
            if det.abs() <= 1e-12 {
                return (dot(n1, n2) > 0.0 && (k1 - k2).abs() <= eps).then(|| a.foot(near));
            }
            Some([(k1 * n2[1] - k2 * n1[1]) / det, (n1[0] * k2 - n2[0] * k1) / det])
        }
        (Wall::Line(n, k), Wall::Circle(c, r, _)) | (Wall::Circle(c, r, _), Wall::Line(n, k)) => {
            let h = k - dot(n, c);
            let foot = add(c, mul(n, h));
            let q = r * r - h * h;
            if q < -2.0 * r * eps {
                return None;
            }
            let s = q.max(0.0).sqrt();
            let t = [-n[1], n[0]];
            Some(pick([add(foot, mul(t, s)), sub(foot, mul(t, s))]))
        }
        (Wall::Circle(c1, r1, _), Wall::Circle(c2, r2, _)) => {
            let e = sub(c2, c1);
            let d = norm(e);
            if d <= eps {
                return ((r1 - r2).abs() <= eps).then(|| a.foot(near));
            }
            let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
            let q = r1 * r1 - along * along;
            if q < -2.0 * r1 * eps {
                return None;
            }
            let (u, s) = (mul(e, 1.0 / d), q.max(0.0).sqrt());
            let m = add(c1, mul(u, along));
            let t = [-u[1], u[0]];
            Some(pick([add(m, mul(t, s)), sub(m, mul(t, s))]))
        }
    }
}

/// A profile as walls, counter-clockwise: its corners, each edge's bulge
/// and draft slope (tan of its angle), each corner's round and turn (1
/// convex, -1 reflex, 0 smooth), and each corner's index as given.
#[derive(Clone)]
pub(crate) struct Walls {
    pts: Vec<P>,
    bs: Vec<f64>,
    slopes: Vec<f64>,
    rounds: Vec<f64>,
    walls: Vec<Wall>,
    turn: Vec<f64>,
    given: Vec<usize>,
    eps: f64,
}

/// Twice the area between the arc from a to q of bulge b and its chord,
/// signed as the bulge (to the chord's right): the arc's term in a loop's
/// shoelace sum beside its chord's.
fn segment(a: P, q: P, b: f64) -> f64 {
    if b == 0.0 {
        return 0.0;
    }
    let (th, r) = (4.0 * b.abs().atan(), arc_of(a, q, b).r);
    b.signum() * r * r * (th - th.sin())
}

/// The tangents of the edge from a to q of bulge b at its two ends.
fn tangents(a: P, q: P, b: f64) -> (P, P) {
    let t = unit(sub(q, a));
    let h = 2.0 * b.atan(); // half the sweep
    let rot = |v: P, th: f64| [v[0] * th.cos() - v[1] * th.sin(), v[0] * th.sin() + v[1] * th.cos()];
    (rot(t, -h), rot(t, h))
}

impl Walls {
    /// From (u v bulge round) points and each edge's draft slope.
    pub(crate) fn new(input: &[[f64; 4]], slopes: &[f64]) -> Result<Walls, String> {
        let scale = {
            let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
            for p in input {
                for k in 0..2 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            norm(sub(hi, lo)).max(1e-300)
        };
        let eps = 1e-9 * scale;
        let (mut pts, mut bs, mut rounds, mut ss, mut given) = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for (i, p) in input.iter().enumerate() {
            if let Some(&last) = pts.last()
                && norm(sub([p[0], p[1]], last)) <= eps
            {
                continue;
            }
            pts.push([p[0], p[1]]);
            bs.push(p[2]);
            rounds.push(p[3]);
            ss.push(slopes[i]);
            given.push(i);
        }
        while pts.len() > 1 && norm(sub(pts[0], pts[pts.len() - 1])) <= eps {
            pts.pop();
            bs.pop();
            rounds.pop();
            ss.pop();
            given.pop();
        }
        let n = pts.len();
        if n < 2 || (n < 3 && bs.iter().all(|&b| b == 0.0)) {
            return Err(format!("a polygon needs 3 corners (or arcs to close it); it has {n}"));
        }
        // Its turn: its chord polygon's area and each arc's segment (a circle
        // of two half arcs has no chord area).
        let area: f64 = (0..n).map(|i| cross(pts[i], pts[(i + 1) % n]) + segment(pts[i], pts[(i + 1) % n], bs[i])).sum();
        if area < 0.0 {
            // Corner k is old corner n-1-k; edge k (from it) is old edge n-2-k, reversed.
            let edge = |k: usize| (2 * n - 2 - k) % n;
            pts.reverse();
            rounds.reverse();
            given.reverse();
            bs = (0..n).map(|k| -bs[edge(k)]).collect();
            ss = (0..n).map(|k| ss[edge(k)]).collect();
        }
        let walls: Vec<Wall> = (0..n)
            .map(|i| {
                let (a, q) = (pts[i], pts[(i + 1) % n]);
                if bs[i] == 0.0 {
                    let nn = right(a, q);
                    Wall::Line(nn, dot(nn, a))
                } else {
                    let arc = arc_of(a, q, bs[i]);
                    Wall::Circle(arc.c, arc.r, bs[i].signum())
                }
            })
            .collect();
        let turn: Vec<f64> = (0..n)
            .map(|i| {
                let j = (i + n - 1) % n;
                let (_, t_in) = tangents(pts[j], pts[i], bs[j]);
                let (t_out, _) = tangents(pts[i], pts[(i + 1) % n], bs[i]);
                let s = cross(t_in, t_out).atan2(dot(t_in, t_out));
                if s.abs() <= 1e-6 { 0.0 } else { s.signum() }
            })
            .collect();
        let w = Walls { pts, bs, slopes: ss, rounds, walls, turn, given, eps };
        for i in 0..n {
            if w.rounds[i] < 0.0 {
                return Err(format!("the polygon's corner at {} has a negative round", w.at(i)));
            }
            if w.rounds[i] > 0.0 && w.turn[i] == 0.0 {
                return Err(format!("the polygon's corner at {} is smooth (its edges are tangent): it has no corner to round", w.at(i)));
            }
        }
        Ok(w)
    }

    fn at(&self, i: usize) -> String {
        let p = self.pts[i];
        format!("({}, {}) (point {})", p[0], p[1], self.given[i] + 1)
    }

    /// The profile at height z: its points, (u v bulge), each with where it
    /// comes from.
    pub(crate) fn solve(&self, z: f64) -> Result<(Vec<[f64; 3]>, Vec<Src>), String> {
        let n = self.pts.len();
        let w: Vec<Wall> = (0..n).map(|i| self.walls[i].moved(z * self.slopes[i])).collect();
        let lost = |i: usize, what: &str| format!("the polygon's corner at {}: {what} at z = {z}", self.at(i));
        // Each corner's entry and exit points, and its round's bulge.
        let mut ends: Vec<(P, P, f64)> = Vec::with_capacity(n);
        for i in 0..n {
            let (a, b) = (w[(i + n - 1) % n], w[i]);
            let corner = meet(a, b, self.pts[i], self.eps)
                .ok_or_else(|| lost(i, "its two walls no longer meet (tangent walls that lean differently part: round the corner, or give them one draft)"))?;
            let r = self.rounds[i];
            if r == 0.0 {
                ends.push((corner, corner, 0.0));
                continue;
            }
            let s = self.turn[i];
            let c = meet(a.moved(-s * r), b.moved(-s * r), corner, self.eps).ok_or_else(|| lost(i, "its round does not fit between its walls"))?;
            let (t0, t1) = (a.foot(c), b.foot(c));
            let (u0, u1) = (sub(t0, c), sub(t1, c));
            let sweep = cross(u0, u1).atan2(dot(u0, u1));
            if sweep * s <= 0.0 {
                return Err(lost(i, "its round does not fit between its walls"));
            }
            ends.push((t0, t1, (sweep / 4.0).tan()));
        }
        let (mut out, mut of) = (Vec::with_capacity(2 * n), Vec::with_capacity(2 * n));
        for i in 0..n {
            let (t0, t1, rb) = ends[i];
            let next = ends[(i + 1) % n].0;
            let short = || format!("the polygon's edge from {} to {}: its corners' rounds overlap, or its walls cross, at z = {z}", self.at(i), self.at((i + 1) % n));
            // The edge from this corner's exit to the next one's entry.
            let b = match w[i] {
                Wall::Line(..) => {
                    if dot(sub(next, t1), sub(self.pts[(i + 1) % n], self.pts[i])) <= 0.0 {
                        return Err(short());
                    }
                    0.0
                }
                Wall::Circle(c, _, _) => {
                    let sweep0 = 4.0 * self.bs[i].atan();
                    let (u, v) = (sub(t1, c), sub(next, c));
                    let mut sweep = cross(u, v).atan2(dot(u, v));
                    if sweep0 > 0.0 && sweep <= 0.0 {
                        sweep += std::f64::consts::TAU;
                    } else if sweep0 < 0.0 && sweep >= 0.0 {
                        sweep -= std::f64::consts::TAU;
                    }
                    if (sweep - sweep0).abs() >= std::f64::consts::FRAC_PI_2 {
                        return Err(short());
                    }
                    (sweep / 4.0).tan()
                }
            };
            if rb != 0.0 {
                out.push([t0[0], t0[1], rb]);
                of.push(Src::RoundIn(i));
            }
            out.push([t1[0], t1[1], b]);
            of.push(if rb != 0.0 { Src::RoundOut(i) } else { Src::Corner(i) });
        }
        Ok((out, of))
    }

    /// The profile drafted over [z0, z1]: polygon's construction, built at
    /// a height and moved with its walls (`Moving`: each point on its path,
    /// each leaf with them) over the part of the span it is proven to hold
    /// over. Where it stops holding (as the walls lean, a point crosses a
    /// hull line, a contact runs off its arc's end: the hull or a pocket
    /// changes), the next construction is built where the profile builds
    /// differently (`next`) and must hold back past where the last stopped:
    /// the two pieces overlap, each cut to its part of the span (`Piece`),
    /// so the field (the least of them, each with its part's slab) has no
    /// crease. A piece proven to within the allowance of the span's end
    /// (its points move less than it between) covers the rest: an end where
    /// the profile degenerates (a chamfer running out at a corner, its edge
    /// or its whole section reaching zero there) is read as the profile
    /// just inside it.
    pub(crate) fn drafted(&self, z0: f64, z1: f64) -> Result<Drafted, String> {
        // One allowance over the span: its larger end's (an end whose
        // profile degenerates gives none: the other's, else the profile's as
        // given).
        let eps = [z0, z1].into_iter().filter_map(|z| self.solve(z).ok()).map(|(q, _)| allowance(&q)).reduce(f64::max).unwrap_or(self.eps);
        let slack = 100.0 * eps;
        let mut pieces: Vec<Piece> = Vec::new();
        let mut why = String::new();
        let (mut from, mut at) = (z0, z0);
        loop {
            let c = match self.next(pieces.last(), [from, z1], at, eps) {
                Ok(c) => c,
                Err(e) => {
                    // Stopped within the allowance of the span's end: it covers it.
                    if let Some(p) = pieces.last_mut()
                        && p.mv.still(p.held[1], z1, slack)
                    {
                        p.range[1] = z1;
                        break;
                    }
                    return Err(e.unwrap_or(why));
                }
            };
            let Try { t, mv, lo, up, fail, which, h } = c;
            let bottom = match pieces.last_mut() {
                None => z0,
                Some(p) => {
                    let third = (p.range[1] - lo) / 3.0;
                    p.range[1] -= third;
                    p.held[1] = p.range[1];
                    lo + third
                }
            };
            // The heights its leaves read, their t (as the field reads it:
            // `Moving::t_at`) within where it is proven.
            let mut hi = mv.height(up);
            while mv.t_at(hi).hi > up {
                hi = hi.next_down();
            }
            let done = up >= mv.t_at(z1).hi;
            if let Some(w) = which {
                why = w;
            }
            let top = if done { z1 } else { hi };
            let held = [if pieces.is_empty() { lo } else { bottom }, top];
            pieces.push(Piece { t, range: [bottom, top], held, mv });
            if done {
                break;
            }
            (from, at) = (h, fail);
        }
        let mut box_ = ([f64::MAX; 2], [f64::MIN; 2]);
        for p in &pieces {
            let (lo, hi) = p.mv.extent(p.held);
            for k in 0..2 {
                box_.0[k] = box_.0[k].min(lo[k]);
                box_.1[k] = box_.1[k].max(hi[k]);
            }
        }
        Ok(Drafted { pieces, lo: box_.0, hi: box_.1 })
    }

    /// The next piece: the construction built where the last stopped (at),
    /// else the nearest built farther on (halving back from the span's end
    /// toward there), the first that `attempt` takes; else why the profile
    /// does not build where the last stopped or nearest past there (the
    /// first piece: else why none built at the span's start holds), if it
    /// does not.
    fn next(&self, last: Option<&Piece>, window: [f64; 2], at: f64, eps: f64) -> Result<Try, Option<String>> {
        let farther = std::iter::successors(Some(window[1]), |&h| {
            let m = 0.5 * (at + h);
            (m != at && m != h).then_some(m)
        });
        let (mut here, mut near, mut held, mut inside) = (None, None, None, f64::MAX);
        for (k, h) in std::iter::once(at).chain(farther).enumerate() {
            match self.attempt(last, window, h, eps) {
                Ok(c) => return Ok(c),
                Err(Not::Built(e)) if k == 0 => here = Some(e),
                Err(Not::Built(e)) => near = Some(e),
                Err(Not::Held(e)) if k == 0 && last.is_none() => held = Some(e),
                Err(Not::Held(_)) => {}
                Err(Not::Short(lo)) => inside = inside.min(lo),
                Err(Not::Cert(e)) => return Err(Some(e)),
            }
        }
        let why = here.or(near).or(held);
        Err(if last.is_none() && inside < window[1] { why.map(|e| format!("{e}. Built inside the span, it holds from z = {inside}. Start the span there")) } else { why })
    }

    /// The construction built at h, moving over the window, where it builds
    /// there plainly (no two of its points, nor an arc and its chord, within
    /// the allowance: merged there, a point would move as only one of its
    /// definitions), differs from the last piece's, holds back past where
    /// that one was proven to (the first: to within the allowance of the
    /// window's start) and on past there.
    fn attempt(&self, last: Option<&Piece>, window: [f64; 2], h: f64, eps: f64) -> Result<Try, Not> {
        let (q, src) = self.solve(h).map_err(Not::Built)?;
        let n = q.len();
        for i in 0..n {
            let (a, b) = ([q[i][0], q[i][1]], [q[(i + 1) % n][0], q[(i + 1) % n][1]]);
            let l = norm(sub(b, a));
            if l <= eps || (q[i][2] != 0.0 && (q[i][2] * l / 2.0).abs() <= eps) {
                let what = match src[i] {
                    Src::Corner(c) | Src::RoundOut(c) => format!("edge from {} to {}", self.at(c), self.at((c + 1) % self.pts.len())),
                    Src::RoundIn(c) => format!("round at {}", self.at(c)),
                };
                return Err(Not::Built(format!("the polygon's {what} vanishes at z = {h} (under its allowance, {eps}): shorten the span or the drafts")));
            }
        }
        let (t, l, tr) = build(&q, eps).map_err(|e| Not::Built(format!("{e} (drafted, at z = {h})")))?;
        if last.is_some_and(|p| p.mv.l.defs == l.defs && same(&p.t, &t)) {
            return Err(Not::Held(String::new()));
        }
        let mv = Moving::new(self, l, src, window, h).map_err(Not::Built)?;
        let conds = mv.conditions(&tr.obs);
        let tb = (h - window[0]) / (window[1] - window[0]);
        // An identity of the paths not proven where the construction holds
        // otherwise: refused (the paths built again elsewhere are the same).
        let uncertified = |which: Option<usize>, upto: f64| {
            which.filter(|&w| conds[w].cert).map(|w| {
                Not::Cert(format!(
                    "the drafted polygon near its corner at {}: ink cannot certify that {} past z = {}: the paths the field reads are not proven within polygon's allowance, {:e}, of its walls there (rounding at coordinates far larger than the profile, or walls meeting nearly tangentially, can pass it): refused, not built unproven",
                    mv.corner(conds[w].at),
                    conds[w].what,
                    mv.height(upto),
                    100.0 * eps
                ))
            })
        };
        // Built at the window's start, it holds back to there.
        let below = last.map_or(window[0], |p| p.range[1]);
        let lo = if tb > 0.0 {
            let (back, _, which) = mv.reach(&conds, tb, 0.0);
            let mut lo = mv.height(back);
            while mv.t_at(lo).lo < back {
                lo = lo.next_up();
            }
            if lo >= below
                && let Some(e) = uncertified(which, back)
            {
                return Err(e);
            }
            lo
        } else {
            window[0]
        };
        if !(lo < below || (last.is_none() && (lo <= window[0] || mv.still(window[0], lo, 100.0 * eps)))) {
            return Err(if last.is_none() { Not::Short(lo) } else { Not::Held(String::new()) });
        }
        // On to the t the field reads at the window's end.
        let end = mv.t_at(window[1]).hi;
        let (up, fail, which) = mv.reach(&conds, tb, end);
        if let Some(e) = uncertified(which, up) {
            return Err(e);
        }
        let which = which.map(|w| {
            format!(
                "the drafted polygon near its corner at {}: {} only up to z = {} (as its walls lean its hull, a pocket or an arc's part changes, and built again past there it does not carry on): shorten the span or the drafts",
                mv.corner(conds[w].at),
                conds[w].what,
                mv.height(up)
            )
        });
        if up < end && mv.height(up) <= below {
            return Err(Not::Held(which.unwrap_or_default()));
        }
        Ok(Try { t, h, lo, up, fail: mv.height(fail), which, mv })
    }
}

/// Why a construction built at a height is not taken: the profile does not
/// build plainly there, what it builds does not hold far enough, (the
/// first) holds back only to a height short of the span's start, or its
/// paths are not proven on its walls (the profile is refused).
enum Not {
    Built(String),
    Held(String),
    Short(f64),
    Cert(String),
}

/// A construction built at height h and proven over its window: the
/// heights it holds back to (lo) and on to (up, in t), where it failed on
/// (fail), and why it stopped.
struct Try {
    t: T,
    mv: Moving,
    h: f64,
    lo: f64,
    up: f64,
    fail: f64,
    which: Option<String>,
}

/// Whether two trees are one construction: the same shape, their leaves
/// from the same loop edges.
fn same(a: &T, b: &T) -> bool {
    match (a, b) {
        (T::Lin(.., x), T::Lin(.., y)) => x == y,
        (T::Disc(.., x), T::Disc(.., y)) => x == y,
        (T::Neg(a), T::Neg(b)) => same(a, b),
        (T::Max(a), T::Max(b)) | (T::Min(a), T::Min(b)) => a.len() == b.len() && a.iter().zip(b).all(|(a, b)| same(a, b)),
        _ => false,
    }
}

/// Height z of t in [0, 1] over a window from z0 to z1.
fn mv_height(z0: f64, z1: f64, t: f64) -> f64 {
    z0 + t * (z1 - z0)
}

/// A drafted profile: its pieces, and the box that holds it over its span.
pub(crate) struct Drafted {
    pub(crate) pieces: Vec<Piece>,
    pub(crate) lo: P,
    pub(crate) hi: P,
}

/// One construction of a drafted profile: its tree as built, the part of
/// the span it is cut to, the part its leaves read (where it is proven: at
/// a span's end, it covers the allowance beyond as the profile there), and
/// its loop moving with its walls.
pub(crate) struct Piece {
    pub(crate) t: T,
    pub(crate) range: [f64; 2],
    pub(crate) held: [f64; 2],
    pub(crate) mv: Moving,
}

/// Where a point of the solved profile comes from: corner i (where its two
/// walls meet; wall i's edge follows), or where corner i's round meets wall
/// i - 1 (the round's arc follows) or wall i (wall i's edge follows).
#[derive(Clone, Copy, Debug)]
pub(crate) enum Src {
    Corner(usize),
    RoundIn(usize),
    RoundOut(usize),
}

/// The circle an edge lies on: wall i's, or corner i's round's.
#[derive(Clone, Copy, Debug)]
enum On {
    Wall,
    Round(usize),
}

/// A condition on the moving loop: the quantities it reads, and whether it
/// holds given their bounds; what it says, and the loop point it is near.
struct Cond {
    qs: Tape,
    holds: Box<dyn Fn(&[I]) -> bool>,
    what: &'static str,
    at: usize,
    /// An identity of the paths (`identities`): where it is not proven the
    /// profile is refused, not built again.
    cert: bool,
}

/// a op b, rounded outward.
fn iv(op: BinOp, a: I, b: I) -> I {
    ink_core::interval::bin(op, a, b).0
}

impl Cond {
    /// Whether it holds over all of t, by the quantities' bounds in slope
    /// form.
    fn over(&self, t: I) -> bool {
        let b: Vec<I> = self.qs.at(Sl::of(t)).into_iter().map(|q| q.e).collect();
        (self.holds)(&b)
    }
}

/// A drafted profile's construction moving with its walls over a window of
/// the span (t from 0 at z[0] to 1 at z[1]; built at height h). Each point
/// moves as its definition does, in closed form (curve.rs): a corner of two
/// walls where they meet (two planes: linearly; a plane or cone and a cone:
/// on a curve), a round's ends at its feet on its walls from its centre
/// (where its walls, offset by its radius, meet), a halving at a fixed
/// direction from its arc's centre, a contact where the hull line from its
/// point, or common to its two arcs, touches its arc (linearly where its
/// normal stays the same). Each edge's arc reads its one circle (a wall's,
/// its centre fixed and its radius moving linearly, or a round's), as do
/// the points defined on it.
///
/// What polygon's construction assumed of its loop (`Ob`), and what it
/// needs of any loop (simple, no corner folding back, each arc's part
/// under a half turn), and that each curve's square roots stay positive
/// (its walls meet, its tangent exists), are conditions on the signs of
/// expressions in t (`Cond`), proven over parts of the window by their
/// bounds in second-order slope form (`Sl`), to within polygon's allowance
/// for its data, 100 eps (where one is broken by less, the field differs
/// from the construction's by a sliver as thin). What the construction
/// makes true by its definitions is proven too, on the paths the field
/// reads (`identities`): in real arithmetic each is an identity, but the
/// paths are folded in floats and some are replaced by a path between
/// their ends, so each is checked, and where it is not proven the profile
/// is refused. Each condition reads the paths through exact operations
/// (`curve::exact`), so its bounds hold of the real values of the
/// expressions the field lowers. The identities pair t with the height at
/// which the field reads it (`z_of_t`), exactly.
pub(crate) struct Moving {
    w: Walls,
    /// The build height, in t.
    tb: f64,
    l: Loop,
    pub(crate) z: [f64; 2],
    src: Vec<Src>,
    /// The circle each input edge lies on.
    edges: Vec<On>,
    /// Each loop point, moving.
    pub(crate) pts: Vec<XP>,
    /// Each loop edge's circle (a line's unused).
    pub(crate) circles: Vec<(XP, X)>,
    /// The line walls each loop point lies on (a line through two points on
    /// one is that wall: its plane), and each line wall's n, k and slope.
    pub(crate) on_walls: Vec<Vec<usize>>,
    pub(crate) planes: Vec<Option<(P, f64, f64)>>,
    /// Each contact moving as `curve::touch` builds it (not replaced): the
    /// point its line runs to, whether that is a contact common to two arcs,
    /// and its side.
    touched: Vec<Option<(usize, bool, f64)>>,
}

impl Moving {
    fn new(w: &Walls, l: Loop, src: Vec<Src>, z: [f64; 2], h: f64) -> Result<Moving, String> {
        let n = l.n();
        let edges = src
            .iter()
            .map(|s| match *s {
                Src::Corner(_) | Src::RoundOut(_) => On::Wall,
                Src::RoundIn(i) => On::Round(i),
            })
            .collect();
        let planes = w.walls.iter().zip(&w.slopes).map(|(wl, &s)| if let Wall::Line(n, k) = *wl { Some((n, k, s)) } else { None }).collect();
        let tb = (h - z[0]) / (z[1] - z[0]);
        let mut m = Moving { w: w.clone(), tb, l, z, src, edges, pts: Vec::new(), circles: Vec::new(), on_walls: Vec::new(), planes, touched: vec![None; n] };
        m.circles = (0..n).map(|e| m.edge_circle(e, tb)).collect();
        let mut pts: Vec<Option<XP>> = vec![None; n];
        let mut touched = vec![None; n];
        for i in 0..n {
            m.point(i, tb, &mut pts, &mut touched)?;
        }
        m.touched = touched;
        m.pts = pts.into_iter().map(|p| p.expect("every point moves")).collect();
        // The line walls each point lies on.
        m.on_walls = (0..n)
            .map(|i| match m.l.defs[i] {
                Def::Given(k) => {
                    let nw = m.w.walls.len();
                    let ws = match m.src[k] {
                        Src::Corner(c) => vec![(c + nw - 1) % nw, c],
                        Src::RoundIn(c) => vec![(c + nw - 1) % nw],
                        Src::RoundOut(c) => vec![c],
                    };
                    ws.into_iter().filter(|&w| matches!(m.w.walls[w], Wall::Line(..))).collect()
                }
                _ => vec![],
            })
            .collect();
        Ok(m)
    }

    /// polygon's eps for the construction (its allowance is 100 eps).
    pub(crate) fn eps(&self) -> f64 {
        self.l.eps
    }

    pub(crate) fn height(&self, t: f64) -> f64 {
        mv_height(self.z[0], self.z[1], t)
    }

    /// Wall i, moved by the draft over the window: at z0 and at z1.
    fn wall(&self, i: usize, by: f64) -> (Wall, Wall) {
        let at = |z: f64| self.w.walls[i].moved(z * self.w.slopes[i]).moved(by);
        (at(self.z[0]), at(self.z[1]))
    }

    /// Where two walls, each at both ends, meet, on the side nearest p at
    /// t; None where both are lines (they meet linearly).
    fn meet_curve(a: (Wall, Wall), b: (Wall, Wall), p: P, t: f64) -> Option<XP> {
        let side = |s: f64| match (a, b) {
            ((Wall::Line(n, k0), Wall::Line(_, k1)), (Wall::Circle(cc, r0, _), Wall::Circle(_, r1, _)))
            | ((Wall::Circle(cc, r0, _), Wall::Circle(_, r1, _)), (Wall::Line(n, k0), Wall::Line(_, k1))) => Some(curve::line_circle(n, [k0, k1], cc, [r0, r1], s)),
            ((Wall::Circle(c1, r10, _), Wall::Circle(_, r11, _)), (Wall::Circle(c2, r20, _), Wall::Circle(_, r21, _))) => {
                Some(curve::circle_circle(c1, [r10, r11], c2, [r20, r21], s))
            }
            _ => None,
        };
        side(1.0)?;
        Some(curve::nearer(|s| side(s).unwrap(), p, t).0)
    }

    /// Whether two walls touch tangentially at both ends (a smooth join):
    /// they stay tangent between (the gap under their meeting's square
    /// root is linear in t and 0 at both ends), meeting at their foot,
    /// which moves linearly. Or two arcs of one circle.
    fn tangent_walls(a: (Wall, Wall), b: (Wall, Wall), eps: f64) -> bool {
        let gap = |(a, b): (Wall, Wall)| match (a, b) {
            (Wall::Line(n, k), Wall::Circle(c, r, _)) | (Wall::Circle(c, r, _), Wall::Line(n, k)) => (r * r - (k - dot(n, c)).powi(2), r),
            (Wall::Circle(c1, r1, _), Wall::Circle(c2, r2, _)) => {
                let d = norm(sub(c2, c1));
                if d <= eps {
                    return (0.0, 1.0);
                }
                let along = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
                (r1 * r1 - along * along, r1.max(r2))
            }
            _ => (1.0, 0.0),
        };
        [(a.0, b.0), (a.1, b.1)].iter().all(|&w| {
            let (q, r) = gap(w);
            q.abs() <= 2.0 * r * 100.0 * eps
        })
    }

    /// Where corner i's walls meet at height z (near the given corner).
    fn corner_at(&self, i: usize, z: f64, by: f64) -> Option<P> {
        let nw = self.w.walls.len();
        let w = |k: usize| self.w.walls[k].moved(z * self.w.slopes[k]);
        let corner = meet(w((i + nw - 1) % nw), w(i), self.w.pts[i], self.w.eps)?;
        if by == 0.0 {
            return Some(corner);
        }
        meet(w((i + nw - 1) % nw).moved(by), w(i).moved(by), corner, self.w.eps)
    }

    /// Corner i's round's centre, moving: its walls offset by its radius
    /// meet there.
    fn round_centre(&self, i: usize, t: f64) -> XP {
        let nw = self.w.walls.len();
        let (r, s) = (self.w.rounds[i], self.w.turn[i]);
        let (a, b) = (self.wall((i + nw - 1) % nw, -s * r), self.wall(i, -s * r));
        let (c0, c1) = (self.corner_at(i, self.z[0], -s * r).unwrap_or(self.w.pts[i]), self.corner_at(i, self.z[1], -s * r).unwrap_or(self.w.pts[i]));
        let at_t = add(c0, mul(sub(c1, c0), t));
        let line = curve::plin(c0, c1);
        match Moving::meet_curve(a, b, at_t, t) {
            Some(x) if !(Moving::tangent_walls(a, b, self.w.eps) && self.stays_on(&line, &[((i + nw - 1) % nw, -s * r), (i, -s * r)])) => x,
            _ => line,
        }
    }

    /// How far p is off wall i offset out by `by` at the height of t: past
    /// its plane n . p = k + s z + by, or off its cone |p - c| = r + sg (s z
    /// + by) (polygon's own data, read exactly: `curve::exact`).
    fn off(&self, i: usize, by: f64, p: &XP) -> X {
        use curve::exact::{kx, pdot, pk, pnorm, psub, xadd, xmul, xsub};
        let z = self.z_of_t();
        let s = self.w.slopes[i];
        match self.w.walls[i] {
            Wall::Line(nv, k) => xsub(&pdot(&pk(nv), p), &xadd(&xadd(&kx(k), &xmul(&kx(s), &z)), &kx(by))),
            Wall::Circle(c, r, sg) => xsub(&pnorm(&psub(p, &pk(c))), &xadd(&xadd(&kx(r), &xmul(&kx(sg * s), &z)), &kx(sg * by))),
        }
    }

    /// The height of t: z0 + t (z1 - z0), its span's length rounded as
    /// `xlin` rounds it. The field reads t at height z as (z - z0) / that
    /// length, in the graph's operations (`Doc::t_of`), so this is the
    /// height at which it reads t, exactly.
    /// The t at which the field reads height z, in outward bounds.
    pub(crate) fn t_at(&self, z: f64) -> I {
        use ink_core::interval::bin;
        bin(BinOp::Div, bin(BinOp::Sub, I::point(z), I::point(self.z[0])).0, I::point(self.z[1] - self.z[0])).0
    }

    fn z_of_t(&self) -> X {
        curve::xlin(self.z[0], self.z[1])
    }

    /// Whether a path stays within the allowance (100 eps) of walls (each
    /// offset out by its `by`) over the whole window, by its bounds over
    /// parts of it, halved until they settle (`curve::reach`): a curve is
    /// replaced by the path between its ends (where walls touch
    /// tangentially at both ends) only where that is proven.
    fn stays_on(&self, p: &XP, walls: &[(usize, f64)]) -> bool {
        let slack = 100.0 * self.l.eps;
        let tape = Tape::new(&walls.iter().map(|&(i, by)| self.off(i, by, p)).collect::<Vec<_>>());
        curve::reach(|t| tape.at(Sl::of(t)).iter().all(|q| q.e.lo >= -slack && q.e.hi <= slack), 0.0, 1.0).0 >= 1.0
    }

    /// Where corner i's round meets wall k (its foot from the round's
    /// centre): off a line by the radius along its normal, on a circle at
    /// its radius along the way to the centre.
    fn round_foot(&self, i: usize, k: usize, t: f64) -> XP {
        let c = self.round_centre(i, t);
        let (w0, w1) = self.wall(k, 0.0);
        match (w0, w1) {
            (Wall::Line(n, k0), Wall::Line(_, k1)) => {
                // c - (n . c - k) n
                let off = curve::xsub(&curve::pdot(&curve::pk(n), &c), &curve::xlin(k0, k1));
                curve::psub(&c, &curve::pscale(&curve::pk(n), &off))
            }
            (Wall::Circle(cc, r0, sg), Wall::Circle(_, r1, _)) => {
                // cc + R / (R - sg s r) (c - cc), |c - cc| being R - sg s r.
                let (r, s) = (self.w.rounds[i], self.w.turn[i]);
                let big = curve::xlin(r0, r1);
                let far = curve::xlin(r0 - sg * s * r, r1 - sg * s * r);
                curve::padd(&curve::pk(cc), &curve::pscale(&curve::psub(&c, &curve::pk(cc)), &curve::xdiv(&big, &far)))
            }
            _ => unreachable!("a wall is a line or a circle at both ends"),
        }
    }

    /// The input point loop edge e runs from (its arc's): from the given
    /// point at or before it.
    fn given(&self, e: usize) -> usize {
        let n = self.l.n();
        let mut k = e % n;
        loop {
            if let Def::Given(g) = self.l.defs[k] {
                return g;
            }
            k = (k + n - 1) % n;
        }
    }

    /// Edge e's arc's circle: a round's (its centre where its walls, offset
    /// by its radius, meet: on a curve where one is an arc) or a wall's (its
    /// centre fixed, its radius moving linearly). Every part of an arc
    /// (halved, split at a contact) reads that one circle, as do the points
    /// defined on it (a contact, a round's foot, a part's end), so that what
    /// the construction makes true of them is an identity of the
    /// expressions the field reads.
    fn edge_circle(&self, e: usize, t: f64) -> (XP, X) {
        let g = self.given(e);
        if let On::Round(i) = self.edges[g] {
            return (self.round_centre(i, t), curve::kx(self.w.rounds[i]));
        }
        let (Src::Corner(w) | Src::RoundOut(w)) = self.src[g] else { unreachable!("a round's edge reads its round") };
        if let (Wall::Circle(c, r0, _), Wall::Circle(_, r1, _)) = self.wall(w, 0.0) {
            return (curve::pk(c), curve::xlin(r0, r1));
        }
        // A line's: unused.
        (curve::pk([0.0, 0.0]), curve::kx(0.0))
    }

    /// Loop point i's path (`pts`: each found once; a contact reads its
    /// partner's), and that it is where the construction put it at t.
    fn point(&self, i: usize, t: f64, pts: &mut Vec<Option<XP>>, touched: &mut Vec<Option<(usize, bool, f64)>>) -> Result<(), String> {
        if pts[i].is_some() {
            return Ok(());
        }
        let nw = self.w.walls.len();
        let at = self.l.pts[i];
        let x = match self.l.defs[i] {
            Def::Given(k) => match self.src[k] {
                Src::Corner(c) => {
                    let (a, b) = (self.wall((c + nw - 1) % nw, 0.0), self.wall(c, 0.0));
                    let ends = (self.corner_at(c, self.z[0], 0.0), self.corner_at(c, self.z[1], 0.0));
                    let line = if let (Some(p0), Some(p1)) = ends { Some(curve::plin(p0, p1)) } else { None };
                    // Two lines meet linearly; where a curve's walls touch
                    // tangentially at both ends, the path between its ends
                    // where that is proven on them (`stays_on`).
                    match (Moving::meet_curve(a, b, at, t), line) {
                        (None, Some(l)) => l,
                        (Some(_), Some(l)) if Moving::tangent_walls(a, b, self.w.eps) && self.stays_on(&l, &[((c + nw - 1) % nw, 0.0), (c, 0.0)]) => l,
                        (Some(x), _) => x,
                        _ => {
                            return Err(format!(
                                "the drafted polygon near its corner at {}: its walls do not meet at both ends of z = {}..{}",
                                self.corner(i),
                                self.z[0],
                                self.z[1]
                            ));
                        }
                    }
                }
                Src::RoundIn(c) => self.round_foot(c, (c + nw - 1) % nw, t),
                Src::RoundOut(c) => self.round_foot(c, c, t),
            },
            Def::Half => {
                // At its direction from its arc's centre where it was built.
                let (c, r) = &self.circles[i];
                let u = unit(sub(at, curve::at(c, t)));
                curve::padd(c, &curve::pscale(&curve::pk(u), r))
            }
            Def::Touch(to) => {
                let j = self.l.ids.iter().position(|&d| d == to).ok_or("a drafted polygon's contact has lost its line")?;
                let (c, r) = self.circles[i].clone();
                let mutual = self.l.defs[j] == Def::Touch(self.l.ids[i]);
                let (c2, r2) = if mutual {
                    self.circles[j].clone()
                } else {
                    self.point(j, t, pts, touched)?;
                    (pts[j].clone().expect("found"), curve::kx(0.0))
                };
                let (x, side) = curve::nearer(|s| curve::touch(&c, &r, &c2, &r2, s).0, at, t);
                let rounded = curve::touch(&c, &r, &c2, &r2, side).1;
                // Its normal the same at both ends, from a point or circle
                // moving linearly: the tangent moves linearly (to rounding,
                // which the conditions bound: it is not `touched`).
                let linear = |p: &XP, r: &X| p.iter().all(curve::linear) && curve::linear(r);
                let normal = |s: f64| unit(sub(curve::at(&x, s), curve::at(&c, s)));
                if linear(&c, &r) && linear(&c2, &r2) && norm(sub(normal(0.0), normal(1.0))) <= 1e-12 {
                    curve::plin(curve::at(&x, 0.0), curve::at(&x, 1.0))
                } else {
                    // Its identities hold of what the field reads only where
                    // no fold rounded it (`curve::touch`).
                    if !rounded {
                        touched[i] = Some((j, mutual, side));
                    }
                    x
                }
            }
        };
        // Every condition is proven on the paths (each on its walls, over
        // the window: `identities`), so this only guards against a wrong
        // one (a curve's other side): a point defined by walls or circles
        // meeting within the allowance of tangent shifts along them by up to
        // the half-chord whose sag is the allowance (at the profile's size,
        // 1e9 eps).
        if !(norm(sub(curve::at(&x, t), at)) <= (200.0 * self.l.eps * 1e9 * self.l.eps).sqrt()) {
            return Err(format!("the drafted polygon near its corner at {}: a point moving on its path is not where the profile has it at z = {}", self.corner(i), self.height(t)));
        }
        pts[i] = Some(x);
        Ok(())
    }

    /// Whether loop point p is a contact moving as `curve::touch` builds it
    /// whose line runs to loop point q (a contact common to two arcs: each
    /// built so, on opposite sides, so with one normal).
    fn exact_contact(&self, p: usize, q: usize) -> bool {
        let n = self.l.n();
        let (p, q) = (p % n, q % n);
        match self.touched[p] {
            Some((j, false, _)) => j == q,
            Some((j, true, s)) => j == q && self.touched[q] == Some((p, true, -s)),
            None => false,
        }
    }

    fn p(&self, i: usize) -> XP {
        self.pts[i % self.l.n()].clone()
    }

    fn normal(&self, n: Normal) -> XP {
        use curve::exact::{kx, psub, xsub};
        match n {
            Normal::Line(a, b) => {
                let d = psub(&self.p(b), &self.p(a));
                [d[1].clone(), xsub(&kx(0.0), &d[0])]
            }
            Normal::Arc(e, p) => psub(&self.p(p), &self.circles[e % self.l.n()].0),
            Normal::Toward(e, p) => psub(&self.circles[e % self.l.n()].0, &self.p(p)),
        }
    }

    /// The conditions the construction holds by over the window.
    fn conditions(&self, obs: &[Ob]) -> Vec<Cond> {
        use curve::exact::{kx, pcross, pdot, pnorm, psub, xadd, xdiv, xmul};
        let n = self.l.n();
        let eps = self.l.eps;
        let slack = 100.0 * eps;
        let mut out: Vec<Cond> = Vec::new();
        let cond = |qs: Vec<X>, holds: Box<dyn Fn(&[I]) -> bool>, what: &'static str, at: usize| Cond { qs: Tape::new(&qs), holds, what, at, cert: false };
        // A curve's square roots stay positive: its walls meet, its tangent
        // exists (where one reaches 0 its rate has no bound).
        let mut seen = std::collections::HashSet::new();
        for i in 0..n {
            let mut rads = Vec::new();
            self.pts[i].iter().chain(self.circles[i].0.iter()).for_each(|x| Ex::radicands(x, &mut rads));
            for r in rads {
                if seen.insert(std::rc::Rc::as_ptr(&r)) {
                    out.push(cond(vec![r], Box::new(|b| b[0].lo > 0.0), "its walls meet, and its hull lines touch their arcs, nearly tangentially", i));
                }
            }
        }
        // No corner folds back (its edges overlapping), and each arc's part
        // stays under a half turn (its centre on its side of its chord).
        for k in 0..n {
            let (u, v) = (psub(&self.p(k), &self.p(k + n - 1)), psub(&self.p(k + 1), &self.p(k)));
            let (uf, vf) = (curve::at(&u, self.tb), curve::at(&v, self.tb));
            let q0 = cross(uf, vf);
            let class = if q0.abs() > eps * (norm(uf) + norm(vf)) { q0.signum() } else { 0.0 };
            out.push(cond(
                vec![pdot(&u, &v), xmul(&kx(class), &pcross(&u, &v))],
                Box::new(move |b| b[0].lo > 0.0 || (class != 0.0 && b[1].lo > 0.0)),
                "its corner does not fold back",
                k,
            ));
            if self.l.bs[k] != 0.0 {
                let s = self.l.bs[k].signum();
                let chord = psub(&self.p(k + 1), &self.p(k));
                let side = xmul(&kx(s), &pcross(&chord, &psub(&self.circles[k].0, &self.p(k))));
                out.push(cond(vec![side], Box::new(|b| b[0].lo > 0.0), "an arc's part stays under a half turn", k));
            }
        }
        // No two edges meet but at their shared corner: apart by the boxes
        // their ends sweep, or by the signs of the four crosses.
        let boxes: Vec<[I; 2]> = (0..n).map(|i| curve::bounds(&self.pts[i], 0.0, 1.0)).collect();
        for i in 0..n {
            for j in i + 2..n {
                if i == 0 && j == n - 1 {
                    continue;
                }
                let (a, b, c, d) = (boxes[i], boxes[(i + 1) % n], boxes[j], boxes[(j + 1) % n]);
                if (0..2).any(|k| a[k].hi.max(b[k].hi) + eps < c[k].lo.min(d[k].lo) || c[k].hi.max(d[k].hi) + eps < a[k].lo.min(b[k].lo)) {
                    continue;
                }
                let (a, b, c, d) = (self.p(i), self.p(i + 1), self.p(j), self.p(j + 1));
                let (s, o) = (psub(&b, &a), psub(&d, &c));
                let mut qs = vec![pcross(&s, &psub(&c, &a)), pcross(&s, &psub(&d, &a)), pcross(&o, &psub(&a, &c)), pcross(&o, &psub(&b, &c))];
                qs.extend([a, b, c, d].into_iter().flatten());
                out.push(cond(
                    qs,
                    Box::new(move |q| {
                        let same = |x: I, y: I| (x.lo > 0.0 && y.lo > 0.0) || (x.hi < 0.0 && y.hi < 0.0);
                        let (a, b, c, d) = ([q[4], q[5]], [q[6], q[7]], [q[8], q[9]], [q[10], q[11]]);
                        let boxed = (0..2).any(|k| a[k].hi.max(b[k].hi) + eps < c[k].lo.min(d[k].lo) || c[k].hi.max(d[k].hi) + eps < a[k].lo.min(b[k].lo));
                        boxed || same(q[0], q[1]) || same(q[2], q[3])
                    }),
                    "the profile does not cross itself",
                    i,
                ));
            }
        }
        // What the construction assumed.
        for ob in obs {
            match *ob {
                Ob::Side { a, b, p, on } => {
                    let e = psub(&self.p(b), &self.p(a));
                    let d = xdiv(&pcross(&e, &psub(&self.p(p), &self.p(a))), &pnorm(&e));
                    let what = if on { "a point on a hull edge stays on it" } else { "a point stays inside the convex hull of its piece" };
                    out.push(cond(vec![d], Box::new(move |b| b[0].lo >= -slack && (!on || b[0].hi <= slack)), what, p));
                }
                Ob::Along { p, q, a, b } => {
                    let d = pdot(&psub(&self.p(q), &self.p(p)), &psub(&self.p(b), &self.p(a)));
                    out.push(cond(vec![d], Box::new(|b| b[0].lo > 0.0), "a run of edges along a hull edge stays in order", p));
                }
                Ob::Reach { e, a, b } => {
                    // Its ends, and its disc's extreme along the line's
                    // normal where that lies within its sweep, at most 100
                    // eps past the line.
                    let nv = self.normal(Normal::Line(a, b));
                    let nn = pnorm(&nv);
                    let past = |x: XP| xdiv(&pdot(&nv, &psub(&x, &self.p(a))), &nn);
                    let ends: Vec<usize> = [e % n, (e + 1) % n].into_iter().filter(|&x| x != a % n && x != b % n).collect();
                    let mut qs: Vec<X> = ends.iter().map(|&x| past(self.p(x))).collect();
                    let k = qs.len();
                    let (c, r) = self.circles[e % n].clone();
                    // A line drawn to touch this arc where `touch` builds its
                    // contact x (on its circle: `touch`; else, a contact moving
                    // on the path between its ends, it is checked as any
                    // other line, in bounds): it is at right angles
                    // to n, and the circle's centre c = x - r n is r n . n' = r
                    // sigma behind it (n' its unit normal, sigma = n . n' = +-1),
                    // so its disc reaches r (1 - sigma) past it. sigma is fixed
                    // over the window (w > 0, its radicands' condition), so it
                    // is 1 (the disc reaches 0 past) where at the build height
                    // r > 0 and n' . (x - c) = r sigma > 0, in bounds that
                    // hold of the real values.
                    let ab = [a % n, b % n];
                    let tangent = [e % n, (e + 1) % n].into_iter().any(|x| {
                        let other = if x == ab[0] { ab[1] } else { ab[0] };
                        ab.contains(&x) && self.exact_contact(x, other) && {
                            let at = |v: &X| v.at(I::point(self.tb));
                            at(&r).lo > 0.0 && at(&pdot(&nv, &psub(&self.p(x), &c))).lo > 0.0
                        }
                    });
                    if !tangent {
                        let turn = kx(self.l.bs[e % n].signum());
                        qs.push(xmul(&turn, &pcross(&psub(&self.p(e), &c), &nv)));
                        qs.push(xmul(&turn, &pcross(&nv, &psub(&self.p(e + 1), &c))));
                        qs.push(xadd(&past(c), &r));
                    }
                    out.push(cond(
                        qs,
                        Box::new(move |q| q[..k].iter().all(|x| x.hi <= slack) && (tangent || q[k].hi < 0.0 || q[k + 1].hi < 0.0 || q[k + 2].hi <= slack)),
                        "an arc stays within the hull of the piece above it",
                        e,
                    ));
                }
                Ob::Turn { n: nu, m: nm } => {
                    // A contact `touch` builds and its line's: its arc's normal
                    // there is r n, the line at right angles to n
                    // (`curve::touch`), so their cross is 0 and the condition
                    // below holds.
                    let contact = |x: Normal, y: Normal| match (x, y) {
                        (Normal::Arc(e, p), Normal::Line(a, b)) => {
                            (p == a || p == b) && [e % n, (e + 1) % n].contains(&(p % n)) && self.exact_contact(p, if p == a { b } else { a })
                        }
                        _ => false,
                    };
                    if contact(nu, nm) || contact(nm, nu) {
                        continue;
                    }
                    // Turned out by at most turn_out's angle, each arc
                    // bulging out at its largest over the part of the span.
                    // A round's normal and its wall's at its foot are checked
                    // so too: parallel in real arithmetic, their cross is the
                    // paths' rounding, not taken as parallel by
                    // construction.
                    let r = [nu, nm].map(|k| match k {
                        Normal::Arc(e, _) => self.circles[e % n].1.clone(),
                        Normal::Toward(..) | Normal::Line(..) => kx(0.0),
                    });
                    let (u, v) = (self.normal(nu), self.normal(nm));
                    let (q, uv) = (pcross(&u, &v), xmul(&pnorm(&u), &pnorm(&v)));
                    let at = match nu {
                        Normal::Arc(_, p) | Normal::Toward(_, p) | Normal::Line(p, _) => p,
                    };
                    out.push(cond(
                        vec![q, uv, pdot(&u, &v), r[0].clone(), r[1].clone()],
                        Box::new(move |b| {
                            // Convex, or turned in by under the allowance's
                            // turn for each arc's largest radius over the part
                            // (a lower bound: `turn_out`), over the whole part.
                            let s = turn_out([b[3].hi, b[4].hi], slack);
                            b[0].lo >= 0.0 || (iv(BinOp::Add, b[0], iv(BinOp::Mul, I::point(s), b[1])).lo >= 0.0 && b[2].lo > 0.0)
                        }),
                        "an arc's corners stay convex",
                        at,
                    ));
                }
                Ob::Sliver => out.push(cond(vec![], Box::new(|_| false), "a piece is not a sliver", 0)),
            }
        }
        self.identities(&mut out);
        out
    }

    /// What the construction makes true by its definitions, as conditions
    /// on the paths the field reads: each path within the allowance (100
    /// eps) of the wall or circle it is defined on, at every height of the
    /// window. Walls are polygon's own data (n . p = k + s z, a plane, or
    /// |p - c| = r + sg s z, a cone; a round's walls offset by its radius):
    /// the field reads a line wall as that plane (`lower_profile`), and so
    /// each leaf it reads through these points or about these circles lies
    /// within the allowance of the specified profile's. In real arithmetic
    /// each is an identity (a corner is where its walls meet, a round's
    /// centre where they meet offset by its radius, its feet its centre's
    /// projections onto them, a halving or a contact c + r u on its arc's
    /// circle); in floats it holds to rounding where the paths are as
    /// defined, and fails where they are not: motion dropped as rounding,
    /// a replacement by the path between its ends (tangent walls, a
    /// contact whose normal stays put) that strays, coordinates so large
    /// beside the profile that their rounding alone passes the allowance.
    fn identities(&self, out: &mut Vec<Cond>) {
        use curve::exact::{kx, pnorm, psub, xadd, xmul, xsub};
        let (n, nw, slack) = (self.l.n(), self.w.walls.len(), 100.0 * self.l.eps);
        let z = self.z_of_t();
        let off = |i: usize, by: f64, p: &XP| self.off(i, by, p);
        let cert = |qs: Vec<X>, what: &'static str, at: usize| Cond {
            qs: Tape::new(&qs),
            holds: Box::new(move |b| b.iter().all(|q| q.lo >= -slack && q.hi <= slack)),
            what,
            at,
            cert: true,
        };
        // A corner on both its walls; a round's foot on its wall (and on its
        // round: its arc's end, below).
        for i in 0..n {
            let Def::Given(g) = self.l.defs[i] else { continue };
            let (ws, what) = match self.src[g] {
                Src::Corner(c) => (vec![(c + nw - 1) % nw, c], "a corner stays on its walls"),
                Src::RoundIn(c) => (vec![(c + nw - 1) % nw], "a round's end stays on its wall"),
                Src::RoundOut(c) => (vec![c], "a round's end stays on its wall"),
            };
            out.push(cert(ws.iter().map(|&w| off(w, 0.0, &self.pts[i])).collect(), what, i));
        }
        for e in 0..n {
            let g = self.given(e);
            let (c, r) = &self.circles[e];
            match (self.edges[g], self.src[g]) {
                // A round's circle tangent to both its walls: its centre on
                // each offset in by its radius (r inside a line, r off a
                // circle on its side), so its disc meets each wall at its
                // foot, to within the allowance.
                (On::Round(k), _) => {
                    let by = -self.w.turn[k] * self.w.rounds[k];
                    out.push(cert(vec![off((k + nw - 1) % nw, by, c), off(k, by, c)], "a round stays tangent to its walls", e));
                }
                // An arc wall's circle (its centre the wall's own): its
                // radius the cone's.
                (On::Wall, Src::Corner(w) | Src::RoundOut(w)) if matches!(self.w.walls[w], Wall::Circle(..)) => {
                    let Wall::Circle(_, r0, sg) = self.w.walls[w] else { unreachable!() };
                    let cone = xadd(&kx(r0), &xmul(&kx(sg * self.w.slopes[w]), &z));
                    out.push(cert(vec![xsub(r, &cone)], "an arc wall stays its cone", e));
                }
                _ => {}
            }
            // An arc's ends on the circle the field reads it about. A
            // contact `touch` builds on it is: c + r n with n . n = 1, so
            // |p - c| = |r| (`curve::touch`; its arc's parts read one circle,
            // built alike: `edge_circle`), which is r where r >= 0.
            if self.l.bs[e] != 0.0 {
                let on = |p: &XP| xsub(&pnorm(&psub(p, c)), r);
                let ends: Vec<X> = [e, e + 1].into_iter().filter(|&x| self.touched[x % n].is_none()).map(|x| on(&self.p(x))).collect();
                if ends.len() < 2 {
                    out.push(Cond {
                        qs: Tape::new(std::slice::from_ref(r)),
                        holds: Box::new(move |b| b[0].lo >= -0.5 * slack),
                        what: "an arc's radius stays positive",
                        at: e,
                        cert: true,
                    });
                }
                if !ends.is_empty() {
                    out.push(cert(ends, "an arc's ends stay on its circle", e));
                }
            }
        }
    }

    /// How far from t0 toward t1 every condition is proven: the furthest t,
    /// where the first to fail beyond it fails, and which.
    fn reach(&self, conds: &[Cond], t0: f64, t1: f64) -> (f64, f64, Option<usize>) {
        let (mut upto, mut fail, mut which) = (t1, t1, None);
        for (k, c) in conds.iter().enumerate() {
            let (u, f) = curve::reach(|t| c.over(t), t0, fail);
            if u != fail {
                (upto, fail, which) = (u, f, Some(k));
            }
        }
        (upto, fail, which)
    }

    /// Whether every point of the loop and every circle (its centre and
    /// radius) moves less than `by` between heights a and b in its window,
    /// by their bounds there.
    fn still(&self, a: f64, b: f64, by: f64) -> bool {
        let [ta, tb] = [a, b].map(|z| (z - self.z[0]) / (self.z[1] - self.z[0]));
        let (lo, hi) = (ta.min(tb), ta.max(tb));
        let narrow = |i: Option<I>| i.is_some_and(|i| i.hi - i.lo <= by);
        self.pts.iter().chain(self.circles.iter().map(|c| &c.0)).all(|p| p.iter().all(|x| narrow(curve::hull(|t| x.at(Sl::of(t)).e, lo, hi))))
            && self.circles.iter().all(|(_, r)| narrow(curve::hull(|t| r.at(Sl::of(t)).e, lo, hi)))
    }

    /// A box that holds the loop the field reads over the heights [z0, z1]
    /// in its window, and the profile it is proven to within the allowance
    /// of there: each point's path, and each arc's circle's extreme along
    /// each axis its sweep may hold over those heights (an arc's part under
    /// a half turn: its conditions), in outward bounds over the heights'
    /// t (`Doc::t_of`, in outward bounds); grown by the allowance, rounded
    /// outward. (A box of float values at the ends clipped real material by
    /// its rounding.)
    fn extent(&self, z: [f64; 2]) -> (P, P) {
        use curve::exact::{kx, pcross, psub, xadd, xmul};
        let n = self.l.n();
        let dz = self.z[1] - self.z[0];
        let t = |z: f64| ink_core::interval::bin(BinOp::Div, ink_core::interval::bin(BinOp::Sub, I::point(z), I::point(self.z[0])).0, I::point(dz)).0;
        let (ta, tb) = (t(z[0].min(z[1])), t(z[0].max(z[1])));
        let (a, b) = if dz > 0.0 { (ta.lo, tb.hi) } else { (tb.lo, ta.hi) };
        let hull = |x: &X| curve::hull(|t| x.at(Sl::of(t)).e, a, b).unwrap_or(I::ENTIRE);
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        let mut take = |k: usize, v: I| {
            lo[k] = lo[k].min(v.lo);
            hi[k] = hi[k].max(v.hi);
        };
        for p in &self.pts {
            (0..2).for_each(|k| take(k, hull(&p[k])));
        }
        for e in (0..n).filter(|&e| self.l.bs[e] != 0.0) {
            let (c, r) = &self.circles[e];
            let turn = kx(self.l.bs[e].signum());
            for (k, sg) in [(0, 1.0), (1, 1.0), (0, -1.0), (1, -1.0)] {
                let mut u = [kx(0.0), kx(0.0)];
                u[k] = kx(sg);
                // Along u from its centre, outside its sweep throughout?
                let from = xmul(&turn, &pcross(&psub(&self.p(e), c), &u));
                let to = xmul(&turn, &pcross(&u, &psub(&self.p(e + 1), c)));
                if hull(&from).hi < 0.0 || hull(&to).hi < 0.0 {
                    continue;
                }
                take(k, hull(&xadd(&c[k], &xmul(&kx(sg), r))));
            }
        }
        let slack = I::point(100.0 * self.l.eps);
        let out = |v: f64, op: BinOp| ink_core::interval::bin(op, I::point(v), slack).0;
        ([out(lo[0], BinOp::Sub).lo, out(lo[1], BinOp::Sub).lo], [out(hi[0], BinOp::Add).hi, out(hi[1], BinOp::Add).hi])
    }

    /// The reference corner nearest loop point i, for messages.
    fn corner(&self, i: usize) -> String {
        let p = self.l.pts[i % self.l.n()];
        let k = (0..self.w.pts.len()).min_by(|&a, &b| norm(sub(self.w.pts[a], p)).total_cmp(&norm(sub(self.w.pts[b], p)))).unwrap_or(0);
        let q = self.w.pts[k];
        format!("({}, {}) (point {})", q[0], q[1], self.w.given[k] + 1)
    }
}
