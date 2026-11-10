//! A drafted `polygon` of straight walls, rounded or not, built exactly:
//! the flat construction's tree (its convex hull less the pockets between
//! the two, each pocket read the same way) over corners that move with
//! height, and each round as the flat construction reads one (its disc in
//! its triangle), moving with its corner: both its walls are planes moving
//! parallel to themselves, so its centre and tangent points are the
//! corner plus fixed offsets.
//!
//! Each wall is a plane, m . (p - a) = s |m| z (m the edge's outward normal
//! unnormalised, s the draft's slope), exactly the drawn edge at z = 0. A
//! corner is where its two walls meet: c(z) = p + z V, p the drawn point
//! (exact), V known within bounds (|m| is a square root). So every decision
//! of the construction (which side of a line through two corners a third
//! lies on, which of two runs farther along a wall) is the sign of a
//! polynomial in z of degree two at most, its coefficients known within
//! bounds (the one at z = 0 exactly), and the tree is built over a range of
//! heights only where each decision it takes has one sign over the whole
//! range, so it is the profile's construction at every height there.
//!
//! The span is cut where a decision cannot be: halved until it can, each
//! part's tree built over it, parts with one tree joined. Where the tree
//! changes, the height it changes at lies between two heights (floats a
//! float step or so apart) and the one decision that changes there is
//! proven to change once (its slope has one sign between): the field is the
//! one tree below that height and the other above (`select`, on a height
//! known within those two, `Graph::within`), so a range over the band reads
//! both. Lines of the tree are walls (planes) or, through two corners of
//! different walls, lines that turn with height, exact through both
//! corners in real numbers, their constants known within bounds.
//!
//! The span is read outward from z = 0, where the profile is drawn (below
//! 0 mirrored: its slopes negated). Where an edge's length (linear in z)
//! reaches zero, its wall leaves and its neighbours meet (a CAD trim): the
//! profile beyond reads a corner drawn from the vanished edge's point plus
//! an offset within bounds, and the band the edge vanishes in is a change
//! like any other, the construction below holding until it (each decision
//! it takes that degenerates there is the edge's length times a factor
//! proven one sign) and the trimmed one built over it.

use super::exact::{E, diff, mul, orient, safe, sign, sum};
use super::*;
use ink_core::interval::{self as iv, I, add_o, mul_o, sub_o};
use std::rc::Rc;

fn tot(e: &[f64]) -> I {
    e.iter().fold(I::point(0.0), |s, &t| add_o(s, I::point(t)))
}
fn ineg(a: I) -> I {
    I::new(-a.hi, -a.lo)
}
fn icross(u: [I; 2], v: [I; 2]) -> I {
    sub_o(mul_o(u[0], v[1]), mul_o(u[1], v[0]))
}
fn idot(u: [I; 2], v: [I; 2]) -> I {
    add_o(mul_o(u[0], v[0]), mul_o(u[1], v[1]))
}
fn isub(u: [I; 2], v: [I; 2]) -> [I; 2] {
    [sub_o(u[0], v[0]), sub_o(u[1], v[1])]
}
fn mag(a: I) -> f64 {
    a.lo.abs().max(a.hi.abs())
}
fn decided(v: I) -> Option<i8> {
    if v.lo > 0.0 {
        Some(1)
    } else if v.hi < 0.0 {
        Some(-1)
    } else {
        None
    }
}

/// c0 + c1 z + c2 z^2, its coefficients within bounds, and the sign of c0
/// where it is known exactly (2 where not); or, where a corner it reads
/// moves on a curve (`Path`), that decision's enclosures over a range of
/// heights (`Curve`), its coefficients unused.
#[derive(Clone)]
struct Q {
    c: [I; 3],
    s0: i8,
    f: Option<Curve>,
}

/// A decision of corners on curves: its value and its slope over a range of
/// heights, each bounded outward.
#[derive(Clone)]
struct Curve {
    value: Rc<dyn Fn(I) -> I>,
    slope: Rc<dyn Fn(I) -> I>,
}

impl Q {
    fn poly(c: [I; 3], s0: i8) -> Q {
        Q { c, s0, f: None }
    }

    /// Its slope over heights z.
    fn slope(&self, z: I) -> I {
        match &self.f {
            Some(f) => (f.slope)(z),
            None => add_o(self.c[1], mul_o(mul_o(I::point(2.0), self.c[2]), z)),
        }
    }

    /// Its sign over the heights z (an interval not holding 0 inside), or
    /// none where it is not one sign there. At an end z = 0 where c0 is
    /// exactly 0 it is z (c1 + c2 z): the sign over the rest (0 at the end,
    /// where the construction either side is the limit of one that holds).
    fn over(&self, z: I) -> Option<i8> {
        if let Some(f) = &self.f {
            // 0 at z = 0 exactly: its slope's sign over the rest. Else about
            // its middle: q(zm) + (z - zm) q'(z).
            if self.s0 == 0 && (z.lo == 0.0 || z.hi == 0.0) {
                let s = decided((f.slope)(z))?;
                return Some(if z.hi > 0.0 { s } else { -s });
            }
            let m = I::point(z.mid());
            return decided(add_o((f.value)(m), mul_o(sub_o(z, m), (f.slope)(z))));
        }
        let [c0, c1, c2] = self.c;
        if self.s0 == 0 && (z.lo == 0.0 || z.hi == 0.0) {
            let s = decided(add_o(c1, mul_o(c2, z)))?;
            return Some(if z.hi > 0.0 { s } else { -s });
        }
        // About its middle: q(zm) + (z - zm) (c1 + c2 (z + zm)).
        let m = I::point(z.mid());
        let at = add_o(c0, mul_o(m, add_o(c1, mul_o(c2, m))));
        let slope = add_o(c1, mul_o(c2, add_o(z, m)));
        decided(add_o(at, mul_o(sub_o(z, m), slope)))
    }

    /// Whether its slope has one sign over z (so it is 0 once at most).
    fn monotone(&self, z: I) -> bool {
        decided(self.slope(z)).is_some()
    }
}

/// A decision: what it reads (for matching across a change), its
/// polynomial and the sign it took.
#[derive(Clone)]
struct Pred {
    key: (u8, [usize; 3]),
    /// What it reads: an orientation of three corners (0), a run along
    /// (1).
    raw: (u8, [usize; 3]),
    q: Q,
    sign: i8,
}

/// The tree: a line through two corners (positive right of the first to
/// the second), a max, a negation.
#[derive(Clone, Debug, PartialEq)]
enum Tr {
    Line(usize, usize),
    Max(Vec<Tr>),
    /// A pocket, and the hull line it lies under (its lid's line).
    Neg(Box<Tr>, (usize, usize)),
}

#[derive(Clone)]
struct Wall {
    a: P,
    d: [E; 2],
    /// s |m|: how far its plane moves along m per unit z, within bounds.
    sl: I,
    /// Walls on one plane share a class.
    class: usize,
    /// Which wall of the drawn profile it is.
    id: u64,
    /// An arc's wall (a cone): its circle about its drawn start.
    arc: Option<ArcWall>,
}

/// A drafted arc's wall: G = a - C (a its drawn start, C its centre) and
/// R0, within bounds of their exact values (`exact::circle`), sigma (1 where
/// the profile is inside its circle) and its draft's slope: its radius at
/// height z is R0 + sigma s z.
#[derive(Clone, Copy)]
struct ArcWall {
    g: [I; 2],
    r0: I,
    sigma: f64,
    s: f64,
}

/// The profile: its walls (edge i from corner i to i + 1, counter-
/// clockwise) and corners (corner i where walls i - 1 and i meet).
pub(super) struct Drafted {
    walls: Vec<Wall>,
    /// Corner i is p + o + z v: p a drawn point, o an offset within bounds
    /// (0 but where an edge was trimmed), v its speed within bounds.
    p: Vec<P>,
    o: Vec<[I; 2]>,
    v: Vec<[I; 2]>,
    /// Each corner's round (0: none).
    r: Vec<f64>,
    /// Edges `simple` leaves out (`simple_but`).
    skip: Vec<usize>,
    /// What defines it (its drawn points and slopes): the identity its
    /// numbers known within bounds take theirs from.
    key: u64,
    /// The walls trimmed from the drawn profile to make it, in order.
    removed: Vec<u64>,
    /// Each corner's path where it is not p + o + z v (a corner on an
    /// arc's wall moves on a curve).
    curve: Vec<Option<Path>>,
}

/// A corner's path on a curve: its position and its velocity over a range
/// of heights, each bounded outward. At z = 0 it is its drawn point exactly.
#[derive(Clone)]
struct Path {
    at: Rc<dyn Fn(I) -> ([I; 2], [I; 2])>,
    /// What it is, to write it.
    meet: Meet,
}

/// Which walls a corner on a curve is the meeting of.
#[derive(Clone)]
enum Meet {
    /// A line wall and an arc's cone.
    Line(Rc<ConeCorner>),
    /// Two arcs' cones.
    Cone(Rc<ConeCone>),
}

/// One part of the span and its tree, with the decisions it took.
struct Part {
    z: [f64; 2],
    tree: Tr,
    preds: Vec<Pred>,
}

impl Drafted {
    /// From (u v) points (no arcs, no rounds) and each edge's slope; none
    /// where this construction does not read it (repeated points, walls
    /// that are parallel at a corner but not one plane, a profile that
    /// crosses itself at z = 0).
    pub(super) fn new(input: &[[f64; 4]], slopes: &[f64]) -> Option<Drafted> {
        let n0 = input.len();
        // Arcs under half a turn and not flat (flat ones are read about their
        // end throughout: not yet), with no round beside one.
        let arc = |b: f64| b != 0.0;
        if n0 < 3
            || input.iter().any(|p| !(p[3] >= 0.0 && p[3].is_finite()) || (arc(p[2]) && !(p[2].abs() >= 1.0 / 1024.0 && p[2].abs() < 1.0)))
            || (0..n0).any(|i| input[i][..2] == input[(i + 1) % n0][..2])
            || (0..n0).any(|i| input[i][3] != 0.0 && (arc(input[i][2]) || arc(input[(i + n0 - 1) % n0][2])))
        {
            return None;
        }
        let mut pts: Vec<P> = input.iter().map(|p| [p[0], p[1]]).collect();
        let mut ss: Vec<f64> = slopes.to_vec();
        let mut rr: Vec<f64> = input.iter().map(|p| p[3]).collect();
        let mut bs: Vec<f64> = input.iter().map(|p| p[2]).collect();
        // Counter-clockwise (by its lowest point's turn), each slope kept
        // with its edge.
        let lo = (0..n0).fold(0, |m, i| if pts[i][0] < pts[m][0] || (pts[i][0] == pts[m][0] && pts[i][1] < pts[m][1]) { i } else { m });
        match orient(pts[(lo + n0 - 1) % n0], pts[lo], pts[(lo + 1) % n0]) {
            0 => return None,
            o if o < 0 => {
                pts.reverse();
                rr.reverse();
                ss = (0..n0).map(|k| ss[if k + 1 < n0 { n0 - 2 - k } else { n0 - 1 }]).collect();
                bs = (0..n0).map(|k| -bs[if k + 1 < n0 { n0 - 2 - k } else { n0 - 1 }]).collect();
            }
            _ => {}
        }
        // A corner on the line between its neighbours whose walls are one
        // plane (one slope) is no corner; with two slopes, or doubling
        // back, its walls do not meet.
        loop {
            let n = pts.len();
            if n < 3 {
                return None;
            }
            let at = (0..n).find(|&i| bs[(i + n - 1) % n] == 0.0 && bs[i] == 0.0 && orient(pts[(i + n - 1) % n], pts[i], pts[(i + 1) % n]) == 0);
            let Some(i) = at else { break };
            let h = (i + n - 1) % n;
            if super::exact::onward(pts[h], pts[i], pts[(i + 1) % n]) <= 0 || ss[h] != ss[i] || rr[i] != 0.0 {
                return None;
            }
            pts.remove(i);
            ss.remove(i);
            rr.remove(i);
            bs.remove(i);
        }
        let n = pts.len();
        let mut walls: Vec<Wall> = (0..n)
            .map(|i| {
                let (a, b) = (pts[i], pts[(i + 1) % n]);
                let d = [diff(b[0], a[0]), diff(b[1], a[1])];
                let l = iv::un(UnOp::Sqrt, add_o(mul_o(tot(&d[0]), tot(&d[0])), mul_o(tot(&d[1]), tot(&d[1]))));
                // An arc's wall: its circle about its start, never one plane
                // with another (its class its own).
                let arc = (bs[i] != 0.0).then(|| {
                    let c = super::exact::circle(a, b, bs[i]);
                    ArcWall { g: [ineg(I::new(c[0][0], c[0][1])), ineg(I::new(c[1][0], c[1][1]))], r0: I::new(c[2][0], c[2][1]), sigma: bs[i].signum(), s: ss[i] }
                });
                let class = if arc.is_some() { usize::MAX - i } else { i };
                Wall { a, d, sl: mul_o(I::point(ss[i]), l), class, id: i as u64, arc }
            })
            .collect();
        // Planes: one line, one way round, one slope.
        for j in 0..n {
            for k in 0..j {
                let (a, b) = (pts[k], pts[(k + 1) % n]);
                let (dj, dk) = (&walls[j].d, &walls[k].d);
                let exact = safe(&[dj[0].clone(), dj[1].clone(), dk[0].clone(), dk[1].clone()]);
                if walls[j].arc.is_none()
                    && walls[k].arc.is_none()
                    && ss[j] == ss[k]
                    && exact
                    && orient(a, b, pts[j]) == 0
                    && orient(a, b, pts[(j + 1) % n]) == 0
                    && dsign(dj, dk) > 0
                {
                    walls[j].class = walls[k].class;
                    break;
                }
            }
        }
        let mut v = Vec::with_capacity(n);
        let mut curve: Vec<Option<Path>> = vec![None; n];
        for i in 0..n {
            let (wh, wi) = (&walls[(i + n - 1) % n], &walls[i]);
            match (wh.arc, wi.arc) {
                (None, None) => v.push(meet(wh, wi, [w_sl(wh), w_sl(wi)])?),
                // Both walls upright: the corner stays where it is drawn.
                _ if ss[(i + n - 1) % n] == 0.0 && ss[i] == 0.0 => v.push([I::point(0.0); 2]),
                (Some(a1), Some(a2)) => {
                    // Two cones: G at u is the first arc's end's (its start's
                    // plus its chord) and the second's start's.
                    let g1 = [add_o(a1.g[0], tot(&wh.d[0])), add_o(a1.g[1], tot(&wh.d[1]))];
                    let cc = Rc::new(ConeCone::new(pts[i], g1, a2.g, [a1.r0, a2.r0], [a1.sigma * a1.s, a2.sigma * a2.s])?);
                    let c2 = cc.clone();
                    curve[i] = Some(Path { at: Rc::new(move |z: I| c2.at(z)), meet: Meet::Cone(cc) });
                    v.push([I::ENTIRE; 2]);
                }
                (ah, ai) => {
                    // A line and an arc: the corner on the cone (`ConeCorner`),
                    // G at it: the arc's start's, or its end's (plus its chord).
                    let (line, arc, g) = match (ah, ai) {
                        (Some(a), None) => (wi, a, [add_o(a.g[0], tot(&wh.d[0])), add_o(a.g[1], tot(&wh.d[1]))]),
                        (None, Some(a)) => (wh, a, a.g),
                        _ => unreachable!(),
                    };
                    let d = [tot(&line.d[0]), tot(&line.d[1])];
                    let m = m_of(line);
                    let l = iv::un(UnOp::Sqrt, idot(m, m));
                    let mhat = [iv::bin(BinOp::Div, m[0], l).0, iv::bin(BinOp::Div, m[1], l).0];
                    let s_l = if ah.is_some() { ss[i] } else { ss[(i + n - 1) % n] };
                    let cc = Rc::new(ConeCorner::new(pts[i], d, mhat, s_l, g, arc.r0, arc.sigma, arc.s)?);
                    let c2 = cc.clone();
                    curve[i] = Some(Path { at: Rc::new(move |z: I| c2.at(z)), meet: Meet::Line(cc) });
                    v.push([I::ENTIRE; 2]);
                }
            }
        }
        let mut f = ink_core::graph::Fingerprint::default();
        f.str("drafted polygon");
        for (p, s) in input.iter().zip(slopes) {
            p.iter().chain([s]).for_each(|v| _ = f.word(v.to_bits()));
        }
        Some(Drafted { walls, p: pts, o: vec![[I::point(0.0); 2]; n], v, r: rr, skip: Vec::new(), key: f.value() as u64, removed: Vec::new(), curve })
    }

    /// The profile beyond where edge e vanished: its wall left out, walls
    /// e - 1 and e + 1 meeting at a corner drawn from corner e's point (on
    /// wall e - 1 at z = 0) plus an offset. None where they are parallel,
    /// fewer than three walls are left, or a vanishing corner is rounded.
    fn trim(&self, e: usize) -> Option<Drafted> {
        let n = self.n();
        let (h, e1, f) = ((e + n - 1) % n, (e + 1) % n, (e + 2) % n);
        // (Trims beside a corner on a curve: not yet.)
        if n < 4 || self.r[e] != 0.0 || self.r[e1] != 0.0 || self.curved(&[h, e, e1, f]) {
            return None;
        }
        let (wh, wf) = (&self.walls[h], &self.walls[e1]);
        // m_h . d = 0, m_f . d = m_f . (a_f - (p_e + o_e)).
        let pe = [add_o(I::point(self.p[e][0]), self.o[e][0]), add_o(I::point(self.p[e][1]), self.o[e][1])];
        let rhs = idot(m_of(wf), isub([I::point(wf.a[0]), I::point(wf.a[1])], pe));
        let dlt = meet(wh, wf, [I::point(0.0), rhs])?;
        let v = meet(wh, wf, [w_sl(wh), w_sl(wf)])?;
        // In loop order from the new corner: walls e + 1 .. e - 1.
        let order: Vec<usize> = (0..n - 1).map(|k| (f + k) % n).collect();
        let mut removed = self.removed.clone();
        removed.push(self.walls[e].id);
        let mut d = Drafted { walls: Vec::new(), p: Vec::new(), o: Vec::new(), v: Vec::new(), r: Vec::new(), skip: Vec::new(), key: self.key, removed, curve: Vec::new() };
        d.walls.push(wf.clone());
        d.p.push(self.p[e]);
        d.o.push([add_o(self.o[e][0], dlt[0]), add_o(self.o[e][1], dlt[1])]);
        d.v.push(v);
        d.r.push(0.0);
        d.curve.push(None);
        for &k in &order[..n - 2] {
            d.walls.push(self.walls[k].clone());
            d.p.push(self.p[k]);
            d.o.push(self.o[k]);
            d.v.push(self.v[k]);
            d.r.push(self.r[k]);
            d.curve.push(self.curve[k].clone());
        }
        // Walls: e + 1, then e + 2 .. e - 1; corners: the new one (between
        // e - 1 and e + 1), then e + 2 .. e - 1.
        Some(d)
    }

    fn n(&self) -> usize {
        self.p.len()
    }

    fn clone_of(d: &Drafted) -> Drafted {
        Drafted {
            walls: d.walls.clone(),
            p: d.p.clone(),
            o: d.o.clone(),
            v: d.v.clone(),
            r: d.r.clone(),
            skip: d.skip.clone(),
            key: d.key,
            removed: d.removed.clone(),
            curve: d.curve.clone(),
        }
    }

    /// A number's identity (Graph::within) from what defines it: this
    /// profile's data, what it is, and what it is of (walls, a corner's
    /// drawn point, an axis).
    fn ident(&self, what: &str, of: &[u64]) -> u64 {
        let mut f = ink_core::graph::Fingerprint::default();
        f.word(self.key).str(what);
        of.iter().for_each(|&w| _ = f.word(w));
        f.value() as u64
    }

    /// The identity of the height where decision p changes sign in band g
    /// (one root of it there: `change`): what defines it, the profile as
    /// trimmed so far, what the decision reads (its corners by their
    /// walls, a plane by its class), and which of its roots it is (its
    /// slope's sign there; a decision is of degree two at most).
    fn change_id(&self, p: &Pred, g: I) -> Option<(u64, I)> {
        let (kind, ix) = p.key;
        let up = decided(p.q.slope(g))? > 0;
        let mut of: Vec<u64> = self.removed.clone();
        of.push(u64::MAX);
        of.push(kind as u64);
        match kind {
            3 => {
                of.push(ix[0] as u64);
                of.extend(self.cid(ix[1]));
            }
            _ => ix.iter().for_each(|&c| of.extend(self.cid(c))),
        }
        of.push(u64::from(up));
        // A decision of corners on curves: no closed root; its bracket is the
        // band it was found in, and the band is in its identity.
        if p.q.f.is_some() {
            of.extend([g.lo.to_bits(), g.hi.to_bits()]);
            return Some((self.ident("change height", &of), g));
        }
        // Its bracket from what the identity names alone (the decision's own
        // polynomial), never from the band it was found in (the span's).
        let root = if kind == 3 {
            // A corner against a plane's line: its signed distance, linear
            // in z, m . (c_o - a) - s|m| z, of the instance's first wall on
            // that plane.
            let w = self.walls.iter().filter(|w| w.class == ix[0]).min_by_key(|w| w.id)?;
            let o = ix[1];
            let m = m_of(w);
            let at0 = isub([add_o(I::point(self.p[o][0]), self.o[o][0]), add_o(I::point(self.p[o][1]), self.o[o][1])], [I::point(w.a[0]), I::point(w.a[1])]);
            linear_root(idot(m, at0), sub_o(idot(m, self.v[o]), w.sl))?
        } else {
            quadratic_root(&p.q, up)?
        };
        Some((self.ident("change height", &of), root))
    }

    /// The identity of the height where edge e vanishes (its length is
    /// linear in z: one root): the profile as trimmed so far and e's wall.
    fn vanish_id(&self, e: usize) -> Option<(u64, I)> {
        let mut of: Vec<u64> = self.removed.clone();
        of.push(u64::MAX);
        of.push(self.walls[e].id);
        let q = self.wall_q(e, (e + 1) % self.n(), e);
        if q.f.is_some() {
            return None;
        }
        Some((self.ident("vanishing height", &of), linear_root(q.c[0], q.c[1])?))
    }

    /// Corner i's walls (a corner is where its two walls meet).
    fn cid(&self, i: usize) -> [u64; 2] {
        let n = self.n();
        [self.walls[(i + n - 1) % n].id, self.walls[i].id]
    }

    /// Whether any of corners cs moves on a curve.
    fn curved(&self, cs: &[usize]) -> bool {
        cs.iter().any(|&c| self.curve[c].is_some())
    }

    /// Corner i's position and velocity over a range of heights.
    fn motion(&self, i: usize) -> Rc<dyn Fn(I) -> ([I; 2], [I; 2])> {
        if let Some(c) = &self.curve[i] {
            return c.at.clone();
        }
        let (p, v) = ([0, 1].map(|k| add_o(I::point(self.p[i][k]), self.o[i][k])), self.v[i]);
        Rc::new(move |z: I| ([0, 1].map(|k| add_o(p[k], mul_o(z, v[k]))), v))
    }

    /// A decision of corners cs some of which move on curves: its value
    /// from their positions, its slope from their positions and
    /// velocities, over a range of heights (`Curve`).
    fn path_q(&self, cs: &[usize], s0: i8, value: impl Fn(&[[I; 2]]) -> I + 'static, slope: impl Fn(&[[I; 2]], &[[I; 2]]) -> I + 'static) -> Q {
        let ms: Rc<Vec<Rc<dyn Fn(I) -> ([I; 2], [I; 2])>>> = Rc::new(cs.iter().map(|&c| self.motion(c)).collect());
        let m2 = ms.clone();
        let value = Rc::new(move |z: I| value(&ms.iter().map(|m| m(z).0).collect::<Vec<_>>()));
        let slope = Rc::new(move |z: I| {
            let (ps, vs): (Vec<[I; 2]>, Vec<[I; 2]>) = m2.iter().map(|m| m(z)).unzip();
            slope(&ps, &vs)
        });
        Q { c: [I::ENTIRE; 3], s0, f: Some(Curve { value, slope }) }
    }

    /// Whether corners' offsets are all exactly 0 (their points are drawn
    /// points, so decisions at z = 0 are exact predicates on them).
    fn drawn(&self, cs: &[usize]) -> bool {
        cs.iter().all(|&c| self.o[c].iter().all(|v| v.lo == 0.0 && v.hi == 0.0))
    }

    /// c_a - c_b at z = 0 (exact two-term differences, plus offsets).
    fn d0(&self, a: usize, b: usize) -> [I; 2] {
        [0, 1].map(|k| add_o(tot(&diff(self.p[a][k], self.p[b][k])), sub_o(self.o[a][k], self.o[b][k])))
    }

    /// The planes (classes) corner i lies on.
    fn classes(&self, i: usize) -> [usize; 2] {
        let n = self.n();
        [self.walls[(i + n - 1) % n].class, self.walls[i].class]
    }

    /// Whether corners i and j lie on one plane, and which wall's.
    fn shared(&self, i: usize, j: usize) -> Option<usize> {
        let n = self.n();
        let ws = |k: usize| [(k + n - 1) % n, k];
        // (An arc's wall is no plane: corners on one have no line in common.)
        ws(i).into_iter().find(|&a| self.walls[a].arc.is_none() && ws(j).iter().any(|&b| self.walls[a].class == self.walls[b].class))
    }

    /// Corner i over heights z.
    fn at(&self, i: usize, z: I) -> [I; 2] {
        if self.curve[i].is_some() {
            return self.motion(i)(z).0;
        }
        [0, 1].map(|k| add_o(add_o(I::point(self.p[i][k]), self.o[i][k]), mul_o(z, self.v[i][k])))
    }

    /// (c_b - c_a) x (c_p - c_a).
    fn orient_q(&self, a: usize, b: usize, p: usize) -> Q {
        if self.curved(&[a, b, p]) {
            let s0 = if self.drawn(&[a, b, p]) { orient(self.p[a], self.p[b], self.p[p]) } else { 2 };
            return self.path_q(
                &[a, b, p],
                s0,
                |c| icross(isub(c[1], c[0]), isub(c[2], c[0])),
                |c, v| add_o(icross(isub(v[1], v[0]), isub(c[2], c[0])), icross(isub(c[1], c[0]), isub(v[2], v[0]))),
            );
        }
        let (u0, w0) = (self.d0(b, a), self.d0(p, a));
        let (u1, w1) = (isub(self.v[b], self.v[a]), isub(self.v[p], self.v[a]));
        let s0 = if self.drawn(&[a, b, p]) { orient(self.p[a], self.p[b], self.p[p]) } else { 2 };
        Q::poly([icross(u0, w0), add_o(icross(u0, w1), icross(u1, w0)), icross(u1, w1)], s0)
    }

    /// (c_a - c_b) . (c_b - c_o): positive where o -> b -> a runs on.
    fn along_q(&self, a: usize, b: usize, o: usize) -> Q {
        if self.curved(&[a, b, o]) {
            return self.path_q(
                &[a, b, o],
                2,
                |c| idot(isub(c[0], c[1]), isub(c[1], c[2])),
                |c, v| add_o(idot(isub(v[0], v[1]), isub(c[1], c[2])), idot(isub(c[0], c[1]), isub(v[1], v[2]))),
            );
        }
        let e = |u: P, w: P| [diff(u[0], w[0]), diff(u[1], w[1])];
        let (u0, w0) = (e(self.p[a], self.p[b]), e(self.p[b], self.p[o]));
        let (u1, w1) = (isub(self.v[a], self.v[b]), isub(self.v[b], self.v[o]));
        let s0 = if self.drawn(&[a, b, o]) && safe(&[u0[0].clone(), u0[1].clone(), w0[0].clone(), w0[1].clone()]) { dsign(&u0, &w0) } else { 2 };
        let (u0, w0) = (self.d0(a, b), self.d0(b, o));
        Q::poly([idot(u0, w0), add_o(idot(u0, w1), idot(u1, w0)), idot(u1, w1)], s0)
    }

    /// (c_a - c_b) . g, g a direction no profile is symmetric about.
    fn dir_q(&self, a: usize, b: usize) -> Q {
        const G: [f64; 2] = [1.0, 0.618_033_988_749_894_9];
        if self.curved(&[a, b]) {
            let g = [I::point(G[0]), I::point(G[1])];
            return self.path_q(&[a, b], 2, move |c| idot(isub(c[0], c[1]), g), move |_, v| idot(isub(v[0], v[1]), g));
        }
        let e = [diff(self.p[a][0], self.p[b][0]), diff(self.p[a][1], self.p[b][1])];
        let s0 = if self.drawn(&[a, b]) && safe(&e) { sign(&sum(&mul(&e[0], &[G[0]]), &mul(&e[1], &[G[1]]))) } else { 2 };
        let g = [I::point(G[0]), I::point(G[1])];
        Q::poly([idot(self.d0(a, b), g), idot(isub(self.v[a], self.v[b]), g), I::point(0.0)], s0)
    }

    /// (c_b - c_a) . d_w: corner b on from corner a along wall w's way.
    fn wall_q(&self, a: usize, b: usize, w: usize) -> Q {
        if self.curved(&[a, b]) {
            let d = &self.walls[w].d;
            let dd = [tot(&d[0]), tot(&d[1])];
            return self.path_q(&[a, b], 2, move |c| idot(isub(c[1], c[0]), dd), move |_, v| idot(isub(v[1], v[0]), dd));
        }
        let e = [diff(self.p[b][0], self.p[a][0]), diff(self.p[b][1], self.p[a][1])];
        let d = &self.walls[w].d;
        let s0 = if self.drawn(&[a, b]) && safe(&[e[0].clone(), e[1].clone(), d[0].clone(), d[1].clone()]) { dsign(&e, d) } else { 2 };
        let dd = [tot(&d[0]), tot(&d[1])];
        Q::poly([idot(self.d0(b, a), dd), idot(isub(self.v[b], self.v[a]), dd), I::point(0.0)], s0)
    }
}

/// A wall's m = (dy, -dx), its outward normal unnormalised.
fn m_of(w: &Wall) -> [I; 2] {
    [tot(&w.d[1]), ineg(tot(&w.d[0]))]
}

fn w_sl(w: &Wall) -> I {
    w.sl
}

/// x with m_a . x = r[0] and m_b . x = r[1] (two walls' normals); none
/// where they are parallel or it is not finite.
fn meet(wa: &Wall, wb: &Wall, r: [I; 2]) -> Option<[I; 2]> {
    let (m0, m1) = (m_of(wa), m_of(wb));
    let det = icross(m0, m1);
    decided(det)?;
    let div = |x: I| iv::bin(BinOp::Div, x, det).0;
    let x = [div(sub_o(mul_o(r[0], m1[1]), mul_o(r[1], m0[1]))), div(sub_o(mul_o(m0[0], r[1]), mul_o(m1[0], r[0])))];
    x.iter().all(|v| v.lo.is_finite() && v.hi.is_finite()).then_some(x)
}

/// The root of c0 + c1 z, bounded outward; none where c1 may be 0.
fn linear_root(c0: I, c1: I) -> Option<I> {
    decided(c1)?;
    let r = iv::bin(BinOp::Div, ineg(c0), c1).0;
    (r.lo.is_finite() && r.hi.is_finite()).then_some(r)
}

/// The root of q (degree two at most) where its slope is positive (`up`)
/// or negative, bounded outward: the root of a line, or of a quadratic by
/// the form that does not cancel (t = -(c1 + sign(c1) sqrt(disc)) / 2, the
/// roots t / c2 and c0 / t), the one whose slope 2 c2 z + c1 has that
/// sign. None where its bounds cannot tell.
fn quadratic_root(q: &Q, up: bool) -> Option<I> {
    let [c0, c1, c2] = q.c;
    if c2.lo == 0.0 && c2.hi == 0.0 {
        return linear_root(c0, c1).filter(|_| decided(c1) == Some(if up { 1 } else { -1 }));
    }
    let disc = sub_o(mul_o(c1, c1), mul_o(mul_o(I::point(4.0), c0), c2));
    if !(disc.lo > 0.0) {
        return None;
    }
    let sq = iv::un(UnOp::Sqrt, disc);
    let half = |x: I| mul_o(I::point(-0.5), x);
    let ts: Vec<I> = match decided(c1) {
        Some(1) => vec![half(add_o(c1, sq))],
        Some(_) => vec![half(sub_o(c1, sq))],
        None => vec![half(add_o(c1, sq)), half(sub_o(c1, sq))],
    };
    let div = |a: I, b: I| decided(b).map(|_| iv::bin(BinOp::Div, a, b).0);
    let mut cands: Vec<I> = Vec::new();
    for t in ts {
        cands.extend(div(t, c2));
        cands.extend(div(c0, t));
    }
    let want = if up { 1 } else { -1 };
    let mut hit = cands.into_iter().filter(|r| r.lo.is_finite() && r.hi.is_finite() && decided(add_o(mul_o(mul_o(I::point(2.0), c2), *r), c1)) == Some(want));
    let r = hit.next()?;
    // Both candidates of that slope must be the one root (a quadratic has one).
    hit.all(|o| o.lo <= r.hi && r.lo <= o.hi).then_some(r)
}

/// The sign of u . v, exactly.
fn dsign(u: &[E; 2], v: &[E; 2]) -> i8 {
    sign(&sum(&mul(&u[0], &v[0]), &mul(&u[1], &v[1])))
}

/// Decisions over a range of heights, recorded.
struct Oracle<'a> {
    d: &'a Drafted,
    z: I,
    preds: Vec<Pred>,
}

impl Oracle<'_> {
    fn take(&mut self, key: (u8, [usize; 3]), raw: (u8, [usize; 3]), q: Q) -> Option<i8> {
        let s = q.over(self.z)?;
        self.preds.push(Pred { key, raw, q, sign: s });
        Some(s)
    }

    /// The side of a -> b corner p lies on: 0 only on one plane with both.
    fn orient(&mut self, a: usize, b: usize, p: usize) -> Option<i8> {
        if a == p || b == p || a == b {
            return Some(0);
        }
        let (ca, cb, cp) = (self.d.classes(a), self.d.classes(b), self.d.classes(p));
        if ca.iter().any(|c| cb.contains(c) && cp.contains(c)) {
            return Some(0);
        }
        // One key for the three corners in any order (the sign flips with
        // an odd reordering).
        let mut k = [a, b, p];
        let mut flip = 1;
        for i in 0..3 {
            for j in 0..2 - i {
                if k[j] > k[j + 1] {
                    k.swap(j, j + 1);
                    flip = -flip;
                }
            }
        }
        // Two of them on one plane: the third's side of that plane's line,
        // which is 0 at the same heights whichever two they are.
        let pair = [(0, 1, 2), (0, 2, 1), (1, 2, 0)].into_iter().find_map(|(i, j, o)| self.d.shared(k[i], k[j]).map(|w| (self.d.walls[w].class, k[o])));
        let key = match pair {
            Some((c, o)) => (3, [c, o, 0]),
            None => (0, k),
        };
        let s = self.take(key, (0, k), self.d.orient_q(k[0], k[1], k[2]))?;
        Some(s * flip)
    }

    fn along(&mut self, a: usize, b: usize, o: usize) -> Option<i8> {
        self.take((1, [a, b, o]), (1, [a, b, o]), self.d.along_q(a, b, o))
    }

    /// Which of corners a and b lies farther along a fixed direction, to
    /// start the hull's wrap at: any hull corner serves (the wrap proves
    /// each edge of its chain, which is then the hull whichever corner it
    /// started from), so not a decision the tree holds by, and where it is
    /// not one sign over the range, the middle's.
    fn dir(&mut self, a: usize, b: usize) -> Option<i8> {
        let q = self.d.dir_q(a, b);
        Some(q.over(self.z).unwrap_or_else(|| {
            let m = I::point(self.z.mid());
            let v = match &q.f {
                Some(f) => (f.value)(m),
                None => add_o(q.c[0], mul_o(q.c[1], m)),
            };
            if v.mid() > 0.0 { 1 } else { -1 }
        }))
    }
}

impl Drafted {
    /// The tree over heights z, and the decisions it took; none where one
    /// is not one sign over them.
    fn tree(&self, z: I) -> Option<(Tr, Vec<Pred>)> {
        let mut o = Oracle { d: self, z, preds: Vec::new() };
        let root: Vec<usize> = (0..self.n()).collect();
        let t = piece(&mut o, &root, false)?;
        Some((t, o.preds))
    }

    /// Whether the loop is simple over heights z, its edges running their
    /// walls' way: no edge vanishes or turns back, no corner folds, no two
    /// edges meet, it winds counter-clockwise.
    fn simple(&self, z: I) -> bool {
        self.clauses(z).iter().all(|c| c.iter().any(|all| all.iter().all(|(q, sg)| q.over(z) == Some(*sg))))
    }

    /// `simple` as clauses over heights z, each an OR of ANDs of decisions
    /// and the sign each must take (pairs of edges whose boxes over z are
    /// apart need none).
    fn clauses(&self, z: I) -> Vec<Vec<Vec<(Q, i8)>>> {
        let n = self.n();
        let mut out: Vec<Vec<Vec<(Q, i8)>>> = Vec::new();
        for i in 0..n {
            let near = self.skip.iter().any(|&e| [e, (e + 1) % n].contains(&i));
            if !self.skip.contains(&i) {
                out.push(vec![vec![(self.wall_q(i, (i + 1) % n, i), 1)]]);
            }
            // No corner folds: it turns, or its edges run on (a corner
            // passing through straight as it changes from convex to
            // reflex).
            if !near {
                let o = self.orient_q((i + n - 1) % n, i, (i + 1) % n);
                out.push(vec![vec![(o.clone(), 1)], vec![(o, -1)], vec![(self.along_q((i + 1) % n, i, (i + n - 1) % n), 1)]]);
            }
        }
        let bx: Vec<[I; 2]> = (0..n).map(|i| self.at(i, z)).collect();
        let seg = |i: usize| {
            let (a, b) = (bx[i], bx[(i + 1) % n]);
            [I::new(a[0].lo.min(b[0].lo), a[0].hi.max(b[0].hi)), I::new(a[1].lo.min(b[1].lo), a[1].hi.max(b[1].hi))]
        };
        let mut ix: Vec<usize> = (0..n).collect();
        ix.sort_by(|&i, &j| seg(i)[0].lo.total_cmp(&seg(j)[0].lo));
        for (k, &i) in ix.iter().enumerate() {
            let si = seg(i);
            for &j in &ix[k + 1..] {
                let sj = seg(j);
                if sj[0].lo > si[0].hi {
                    break;
                }
                if sj[1].lo > si[1].hi || si[1].lo > sj[1].hi || (i + 1) % n == j || (j + 1) % n == i {
                    continue;
                }
                // The edges beside a vanishing one (`simple_but`).
                if self.skip.iter().any(|&e| [(e + n - 1) % n, (e + 1) % n].contains(&i) && [(e + n - 1) % n, (e + 1) % n].contains(&j)) {
                    continue;
                }
                let (a, b, c, d) = (i, (i + 1) % n, j, (j + 1) % n);
                // Apart: one edge wholly on one side of the other's line,
                // or (on one plane) one after the other along it.
                let (o1, o2, o3, o4) = (self.orient_q(a, b, c), self.orient_q(a, b, d), self.orient_q(c, d, a), self.orient_q(c, d, b));
                let mut cl = vec![vec![(o1.clone(), 1), (o2.clone(), 1)], vec![(o1, -1), (o2, -1)], vec![(o3.clone(), 1), (o4.clone(), 1)], vec![(o3, -1), (o4, -1)]];
                if self.walls[i].class == self.walls[j].class {
                    cl.push(vec![(self.wall_q(b, j, i), 1)]);
                    cl.push(vec![(self.wall_q(d, a, i), 1)]);
                }
                out.push(cl);
            }
        }
        // Counter-clockwise: twice its area, about corner 0, positive.
        let all: Vec<usize> = (0..n).collect();
        if self.curved(&all) {
            let area = |c: &[[I; 2]]| (1..c.len() - 1).fold(I::point(0.0), |a, i| add_o(a, icross(isub(c[i], c[0]), isub(c[i + 1], c[0]))));
            let slope = |c: &[[I; 2]], v: &[[I; 2]]| {
                (1..c.len() - 1).fold(I::point(0.0), |a, i| add_o(a, add_o(icross(isub(v[i], v[0]), isub(c[i + 1], c[0])), icross(isub(c[i], c[0]), isub(v[i + 1], v[0])))))
            };
            out.push(vec![vec![(self.path_q(&all, 2, area, slope), 1)]]);
            return out;
        }
        let mut area = [I::point(0.0); 3];
        for i in 1..n - 1 {
            let (u0, w0) = (self.d0(i, 0), self.d0(i + 1, 0));
            let (u1, w1) = (isub(self.v[i], self.v[0]), isub(self.v[i + 1], self.v[0]));
            area[0] = add_o(area[0], icross(u0, w0));
            area[1] = add_o(area[1], add_o(icross(u0, w1), icross(u1, w0)));
            area[2] = add_o(area[2], icross(u1, w1));
        }
        out.push(vec![vec![(Q::poly(area, 2), 1)]]);
        out
    }

    /// Whether the loop is simple over heights z but for edge e, which
    /// shrinks there (its length falling linearly), its corners turning
    /// as their other edges and d_e do.
    fn simple_but(&self, z: I, es: &[usize]) -> bool {
        if !es.iter().all(|&e| self.shrinks(z, e)) {
            return false;
        }
        let mut t = Drafted::clone_of(self);
        t.skip = es.to_vec();
        t.simple(z)
    }

    /// Edge e's part of `simple_but`: its length falls, its corners turn
    /// as their other edges and d_e do, and the edges beside it do not
    /// cross before it vanishes.
    fn shrinks(&self, z: I, e: usize) -> bool {
        let n = self.n();
        let (h, e1, f) = ((e + n - 1) % n, (e + 1) % n, (e + 2) % n);
        let lam = self.wall_q(e, e1, e);
        if self.curved(&[h, e, e1, f]) || !(lam.c[1].hi < 0.0) {
            return false;
        }
        let d = [tot(&self.walls[e].d[0]), tot(&self.walls[e].d[1])];
        let lin = |a: usize, b: usize, sg: bool| {
            // (c_b - c_a) x d_e, or d_e x (c_b - c_a).
            let c = [icross(self.d0(b, a), d), icross(isub(self.v[b], self.v[a]), d)];
            let q = Q::poly([c[0], c[1], I::point(0.0)], 2);
            let q = if sg { q } else { Q::poly([ineg(q.c[0]), ineg(q.c[1]), I::point(0.0)], 2) };
            q.over(z).is_some()
        };
        if !(lin(h, e, true) && lin(e1, f, false)) {
            return false;
        }
        // Edges e - 1 and e + 1 lie on lines meeting at x (where all three
        // walls meet as e vanishes, so c_e - x and c_e1 - x are 0 there and
        // linear in z): they cross only both holding x, so they do not
        // where, below the root, x lies past c_e along wall e - 1 (its
        // slope (v_e - v_x) . d_(e-1) positive) or short of c_e1 along wall
        // e + 1 ((v_x - v_e1) . d_(e+1) positive).
        let Some(vx) = meet(&self.walls[h], &self.walls[e1], [self.walls[h].sl, self.walls[e1].sl]) else { return false };
        let dh = [tot(&self.walls[h].d[0]), tot(&self.walls[h].d[1])];
        let df = [tot(&self.walls[e1].d[0]), tot(&self.walls[e1].d[1])];
        idot(isub(self.v[e], vx), dh).lo > 0.0 || idot(isub(vx, self.v[e1]), df).lo > 0.0
    }

    /// The edge whose vanishing explains band g after part a (the tree
    /// below), if one does: its length alone among the edges' is not
    /// positive over the band and falls there, the loop is otherwise
    /// simple, and each of a's decisions not one sign over the band reads
    /// both its corners, so is its length times a factor proven one sign
    /// (it holds until the edge vanishes).
    fn vanish(&self, a: &Part, g: I) -> Option<usize> {
        let n = self.n();
        let mut es = (0..n).filter(|&e| self.wall_q(e, (e + 1) % n, e).over(g) != Some(1));
        let e = es.next()?;
        if es.next().is_some() || !self.simple_but(g, &[e]) {
            return None;
        }
        a.preds.iter().all(|p| p.q.over(g) == Some(p.sign) || self.degenerates(p, e, g)).then_some(e)
    }

    /// Whether decision p reads both corners of edge e, so is its length
    /// times a factor, and the factor is one sign over heights g: d_e x
    /// (c_x - c_e), d_e . (c_b - c_o).
    fn degenerates(&self, p: &Pred, e: usize, g: I) -> bool {
        let e1 = (e + 1) % self.n();
        if p.q.f.is_some() || self.curved(&p.raw.1) {
            return false;
        }
        let d = [tot(&self.walls[e].d[0]), tot(&self.walls[e].d[1])];
        let both = |c: &[usize]| c.contains(&e) && c.contains(&e1);
        let lin = |a: usize, b: usize, cross: bool| {
            let (u0, u1) = (self.d0(a, b), isub(self.v[a], self.v[b]));
            let c = if cross { [icross(d, u0), icross(d, u1)] } else { [idot(d, u0), idot(d, u1)] };
            Q::poly([c[0], c[1], I::point(0.0)], 2).over(g).is_some()
        };
        let (kind, ix) = p.raw;
        match kind {
            0 if both(&ix) => ix.iter().find(|&&c| c != e && c != e1).is_some_and(|&x| lin(x, e, true)),
            1 if both(&ix[..2]) => lin(ix[1], ix[2], false),
            1 if both(&ix[1..]) => lin(ix[0], ix[1], false),
            _ => false,
        }
    }

    /// Whether the tree of part a holds over heights g (the last of the
    /// span): each decision it took, and each of `simple`'s, one sign over
    /// g, or (where the profile degenerates just past the span's end)
    /// monotone over g and of its sign at both ends, so 0 nowhere in it.
    fn holds_to(&self, a: &Part, g: I) -> bool {
        let ends = [I::point(g.lo), I::point(g.hi)];
        let ok = |q: &Q, sg: i8| q.over(g) == Some(sg) || (q.monotone(g) && ends.iter().all(|&e| q.over(e) == Some(sg)));
        self.clauses(g).iter().all(|c| c.iter().any(|all| all.iter().all(|(q, sg)| ok(q, *sg)))) && a.preds.iter().all(|p| ok(&p.q, p.sign))
    }

    /// Cover heights r with parts, halving where the tree is not built, up
    /// to the first band (at a float step or so) where it is not: whether
    /// it stopped there.
    fn cover(&self, r: [f64; 2], out: &mut Vec<Result<Part, [f64; 2]>>, budget: &mut usize) -> Option<bool> {
        if *budget == 0 {
            return None;
        }
        *budget -= 1;
        let z = I::new(r[0], r[1]);
        if self.simple(z)
            && let Some((tree, preds)) = self.tree(z)
        {
            out.push(Ok(Part { z: r, tree, preds }));
            return Some(false);
        }
        let m = 0.5 * (r[0] + r[1]);
        if !(m > r[0] && m < r[1]) || (r[1] - r[0]) <= 4.0 * f64::EPSILON * r[0].abs().max(r[1].abs()) {
            out.push(Err(r));
            return Some(true);
        }
        if self.cover([r[0], m], out, budget)? {
            return Some(true);
        }
        self.cover([m, r[1]], out, budget)
    }
}

/// The profile over heights 0..top (its own: a span below 0 is read
/// mirrored, its slopes negated, so it too starts at the drawn profile):
/// its parts, each with the profile it reads (an edge trimmed where it
/// vanished) and its tree, and between two the band the change lies in.
/// None where a part is not built or a change not proven.
fn half(d0: Drafted, top: f64, end: &mut bool) -> Option<(Vec<(Rc<Drafted>, Part)>, Vec<Band>)> {
    // A convex profile is the meet of its walls' half-planes at every
    // height: an edge vanishing is its half-plane no longer cutting (its
    // neighbours meet inside it, the trim), the section closing is them
    // meeting nowhere. So it is the max of its walls' planes throughout.
    let n0 = d0.n();
    // (Not with arcs: a small arc's disc does not hold the profile.)
    if d0.walls.iter().all(|w| w.arc.is_none())
        && d0.r.iter().all(|&r| r == 0.0)
        && d0.drawn(&(0..n0).collect::<Vec<_>>())
        && (0..n0).all(|i| orient(d0.p[(i + n0 - 1) % n0], d0.p[i], d0.p[(i + 1) % n0]) > 0)
    {
        let tree = Tr::Max((0..n0).map(|i| Tr::Line(i, (i + 1) % n0)).collect());
        return Some((vec![(Rc::new(d0), Part { z: [0.0, top], tree, preds: Vec::new() })], Vec::new()));
    }
    let mut d = Rc::new(d0);
    let mut parts: Vec<(Rc<Drafted>, Part)> = Vec::new();
    let mut gaps: Vec<Band> = Vec::new();
    // Ranges tried (each a tree built and checked), at most: a change
    // is found to a float step in 60 or so. With corners on curves a tree
    // costs ten times as much, and a change its decisions cannot show (a
    // corner a hair off a line, its decision near a double root) is given
    // up on only when the budget runs out: there it is in corners times
    // ranges, as much work as 400 ranges of 40 corners, so a large profile
    // is left to the old builder soon and a small one searches as long.
    let mut budget = if d.curve.iter().any(Option::is_some) { (16000 / d.n()).clamp(400, 4000) } else { 4000 };
    let mut from = 0.0;
    // A band the construction changes in, before the next part; a trim's
    // band (and its height's identity), joining the last part to the
    // trimmed profile's first.
    let (mut band, mut join): (Option<[f64; 2]>, Option<Band>) = (None, None);
    loop {
        let mut segs = Vec::new();
        d.cover([from, top], &mut segs, &mut budget)?;
        let mut stop = None;
        for s in segs {
            let p = match s {
                Ok(p) => p,
                Err(g) => {
                    stop = Some(g);
                    break;
                }
            };
            if let Some(j) = join.take() {
                // The trimmed profile, built over the band (from its start).
                if p.z[0] != j.band[0] || p.z[1] < j.band[1].min(top) {
                    return None;
                }
                gaps.push(j);
                parts.push((d.clone(), p));
                continue;
            }
            if let Some(g) = band.take() {
                let (_, last) = parts.last_mut()?;
                let gi = I::new(g[0], g[1]);
                if !d.simple(gi) {
                    return None;
                }
                let at = change(last, &p, gi)?;
                if last.tree == p.tree {
                    last.z[1] = p.z[1];
                    last.preds.extend(p.preds);
                } else {
                    let (id, root) = d.change_id(&at?, gi)?;
                    gaps.push(Band { band: g, root, id });
                    parts.push((d.clone(), p));
                }
                continue;
            }
            match parts.last_mut() {
                Some((_, l)) if l.tree == p.tree => {
                    l.z[1] = p.z[1];
                    l.preds.extend(p.preds);
                }
                // Two trees certified at one height take the same decisions
                // there, so are one tree: not reached.
                Some(_) => return None,
                None => parts.push((d.clone(), p)),
            }
        }
        let Some(g) = stop else { break };
        if join.is_some() {
            return None;
        }
        let g = band.take().map_or(g, |b| [b[0], g[1]]);
        // An edge vanishing at or just past the band (its length falls
        // linearly; its decisions degenerate a few steps before it does):
        // the band reaches its root, bounded; at the top, past it, so
        // above the top is the part below.
        let (_, last) = parts.last()?;
        // Its edges vanishing only past the top: the tree holds to it.
        if d.holds_to(last, I::new(g[0], top)) {
            let (_, last) = parts.last_mut().expect("a part");
            last.z[1] = top;
            break;
        }
        let n = d.n();
        let near = (g[1] - g[0]).max(64.0 * f64::EPSILON * g[1].abs());
        let cand = (0..n)
            .filter_map(|e| {
                let q = d.wall_q(e, (e + 1) % n, e);
                let r = iv::bin(BinOp::Div, ineg(q.c[0]), q.c[1]).0;
                (q.f.is_none() && q.c[1].hi < 0.0 && r.lo.is_finite() && r.hi.is_finite() && r.hi >= g[0] && r.lo <= g[1] + near).then_some((e, r))
            })
            .min_by(|a, b| a.1.lo.total_cmp(&b.1.lo));
        if let Some((e, r)) = cand {
            let hi = if r.hi >= top { top.next_up() } else { g[1].max(r.hi) };
            let b = [g[0], hi];
            let gone = d.vanish(last, I::new(b[0], b[1].min(top))) == Some(e);
            if gone && let Some(t) = d.trim(e) {
                let (id, root) = d.vanish_id(e)?;
                d = Rc::new(t);
                join = Some(Band { band: b, root, id });
                from = g[0];
                continue;
            }
        }
        // Else the tree changes: on past the band, twice as wide each time
        // the next part does not start at its end (a decision changing
        // sign stays within its bounds of 0 for a while past its root).
        if !d.simple(I::new(g[0], g[1])) || g[1] >= top {
            // Within polygon's allowance (1e-7 of its diagonal) of its end:
            // it degenerates there, as rounded data a hair from where an
            // edge or the section closes does.
            *end = top - g[0] <= 100.0 * allowance(&d.p.iter().map(|p| [p[0], p[1], 0.0]).collect::<Vec<_>>());
            return None;
        }
        let next = (g[1] + (g[1] - g[0])).min(top);
        band = Some([g[0], next]);
        from = next;
    }
    if join.is_some() || band.is_some() || parts.is_empty() {
        return None;
    }
    Some((parts, gaps))
}

/// Where a span's construction changes: the band it was found in (the
/// heights the parts either side are read over), its height's bracket from
/// what defines it (an enclosure of the true root, the same in any span:
/// `change_id`, `vanish_id`) and its identity.
#[derive(Clone, Copy)]
struct Band {
    band: [f64; 2],
    root: I,
    id: u64,
}

/// Whether the tree below band g holds up to a height in it and the tree
/// above from there, and the decision that changes there: every decision of each that is not one sign over the
/// band is one (the same for both, matched by what it reads) that changes
/// sign once at most there.
/// That decision, if one is not one sign there.
fn change(a: &Part, b: &Part, g: I) -> Option<Option<Pred>> {
    let mut at: Option<Pred> = None;
    for p in a.preds.iter().chain(&b.preds) {
        if p.q.over(g) == Some(p.sign) {
            continue;
        }
        if !p.q.monotone(g) || at.is_some_and(|k| k.key != p.key) {
            return None;
        }
        at = Some(p.clone());
    }
    Some(at)
}

/// The tree of the piece over corners vs (counter-clockwise; in a pocket
/// the edge from its last to its first is its lid, left out).
fn piece(o: &mut Oracle, vs: &[usize], pocket: bool) -> Option<Tr> {
    let m = vs.len();
    let h = hull(o, vs, pocket)?;
    let nh = h.len();
    let mut terms = Vec::new();
    let mut pockets = Vec::new();
    for k in 0..nh {
        let (i, j) = (h[k], h[(k + 1) % nh]);
        let (u, v) = (vs[i], vs[j]);
        if pocket && i == m - 1 && j == 0 {
            continue;
        }
        terms.push(Tr::Line(u, v));
        // The run under it, cut where it touches the line.
        let len = (j + m - i) % m;
        let mut seg = vec![u];
        for s in 1..=len {
            let p = vs[(i + s) % m];
            seg.push(p);
            let side = if s == len { 0 } else { o.orient(u, v, p)? };
            if side < 0 {
                return None;
            }
            if side == 0 {
                if seg.len() > 2 {
                    seg.reverse();
                    pockets.push((seg.clone(), (u, v)));
                }
                seg = vec![p];
            }
        }
    }
    for (p, lid) in pockets {
        terms.push(Tr::Neg(Box::new(piece(o, &p, true)?), lid));
    }
    Some(Tr::Max(terms))
}

/// The strict convex hull of corners vs, as places in vs in order: by
/// wrapping (each next corner the one every other lies left of, or on and
/// short of), from the lid's start in a pocket, else the corner farthest
/// along a fixed direction.
fn hull(o: &mut Oracle, vs: &[usize], pocket: bool) -> Option<Vec<usize>> {
    let m = vs.len();
    if m < 3 {
        return None;
    }
    let start = if pocket {
        m - 1
    } else {
        let mut best = 0;
        for k in 1..m {
            if o.dir(vs[k], vs[best])? > 0 {
                best = k;
            }
        }
        best
    };
    let mut h = vec![start];
    let mut cur = start;
    loop {
        let mut q = (cur + 1) % m;
        for r in 0..m {
            if r == cur || r == q {
                continue;
            }
            let s = o.orient(vs[cur], vs[q], vs[r])?;
            if s < 0 || (s == 0 && o.along(vs[r], vs[q], vs[cur])? > 0) {
                q = r;
            }
        }
        if q == start {
            break;
        }
        if h.len() > m {
            return None;
        }
        h.push(q);
        cur = q;
    }
    h.sort_unstable();
    Some(h)
}

impl Doc {
    /// `(polygon pts drafts z0 z1)` of straight walls, built exactly (the
    /// module's comment); none where it is not read so.
    /// `(polygon pts drafts z0 z1)` of straight walls built exactly, or
    /// whether it is not only where it degenerates within polygon's
    /// allowance of an end (`half`).
    pub(super) fn exact_drafted(&mut self, input: &[[f64; 4]], slopes: &[f64], z0: f64, z1: f64) -> Result<NodeId, bool> {
        let mut end = false;
        self.exact_drafted_in(input, slopes, z0, z1, &mut end).ok_or(end)
    }

    fn exact_drafted_in(&mut self, input: &[[f64; 4]], slopes: &[f64], z0: f64, z1: f64, end: &mut bool) -> Option<NodeId> {
        // Far from 1 in size, built at a power-of-two scale (exact, slopes
        // kept) putting its size near 1, read at that scale.
        let size = (0..2).map(|c| input.iter().fold(f64::MIN, |m, p| m.max(p[c])) - input.iter().fold(f64::MAX, |m, p| m.min(p[c]))).fold(0.0f64, f64::max);
        let big = input.iter().fold(z0.abs().max(z1.abs()), |m, p| m.max(p[0].abs()).max(p[1].abs()));
        if !(size > 0.0 && size.is_finite() && big.is_finite()) {
            return None;
        }
        let k = if (2f64.powi(-200)..=2f64.powi(200)).contains(&size) && big <= 2f64.powi(400) {
            0
        } else {
            (-size.log2().ceil()).min(400.0 - big.log2().ceil()).clamp(-1000.0, 1000.0) as i32
        };
        let sc = 2f64.powi(k);
        // A span the floats lose at the profile's size (a span 1e-301 high
        // on a profile 1e60 across): its walls move far below its corners'
        // rounding, read by the allowance.
        if !(z0 * sc < z1 * sc && (z1 * sc - z0 * sc).is_finite()) {
            *end = true;
            return None;
        }
        // Every length: coordinates, rounds (bulges and slopes are ratios).
        let scaled: Vec<[f64; 4]> = input.iter().map(|p| [p[0] * sc, p[1] * sc, p[2], p[3] * sc]).collect();
        let (x0, y0, z00) = (self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z));
        let xyz = if k == 0 {
            [x0, y0, z00]
        } else {
            let c = self.g.constant(sc);
            [x0, y0, z00].map(|v| self.g.mul(c, v))
        };
        let f = self.drafted_at(&scaled, slopes, z0 * sc, z1 * sc, xyz, end)?;
        if k == 0 {
            return Some(f);
        }
        let back = self.g.constant(1.0 / sc);
        Some(self.g.mul(back, f))
    }

    fn drafted_at(&mut self, input: &[[f64; 4]], slopes: &[f64], z0: f64, z1: f64, [x, y, z]: [NodeId; 3], end: &mut bool) -> Option<NodeId> {
        if !(z0.is_finite() && z1.is_finite() && z0 < z1) {
            return None;
        }
        // Each side of z = 0 read outward from the drawn profile: above as
        // drawn, below mirrored (its slopes negated, its height -z). Each
        // over its own heights lo..top.
        let neg: Vec<f64> = slopes.iter().map(|s| -s).collect();
        let mut halves = Vec::new();
        if z1 > 0.0 {
            halves.push((false, Drafted::new(input, slopes)?, z0.max(0.0), z1));
        }
        if z0 < 0.0 {
            halves.push((true, Drafted::new(input, &neg)?, (-z1).max(0.0), -z0));
        }
        let mut built = Vec::new();
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for (down, d0, from, top) in halves {
            let rounds = d0.rounds(from, top)?;
            // Its arcs' kites (none without arcs), over its heights.
            let kites = if d0.walls.iter().any(|w| w.arc.is_some()) { d0.kites(top)? } else { Vec::new() };
            for k in &kites {
                for p in [k.t, k.t2] {
                    // (Where it moves to over the span.)
                    let (q, _) = d0.kp(p)(I::new(0.0, top));
                    for c in 0..2 {
                        lo[c] = lo[c].min(q[c].lo);
                        hi[c] = hi[c].max(q[c].hi);
                    }
                }
            }
            let d0r = Rc::new(d0);
            let (parts, gaps) = half(Drafted::clone_of(&d0r), top, end)?;
            // Rounds only where no edge was trimmed.
            if !rounds.is_empty() && parts.iter().any(|(d, _)| d.n() != d0r.n()) {
                return None;
            }
            // The parts reaching lo..top.
            let first = parts.iter().position(|(_, p)| p.z[1] >= from)?;
            let parts: Vec<(Rc<Drafted>, Part)> = parts.into_iter().skip(first).collect();
            let gaps: Vec<Band> = gaps.into_iter().skip(first).collect();
            // Each height's bracket lies within the parts either side of it
            // (past the span's top, which it never reads, it may reach on):
            // there one reads the profile below the root, the other above.
            for (k, b) in gaps.iter().enumerate() {
                let (lo, hi) = (parts[k].1.z[0], parts[k + 1].1.z[1]);
                if !(b.root.lo >= lo && (b.root.hi <= hi || hi >= top)) {
                    return None;
                }
            }
            // The box holding every corner over its heights (outward bounds).
            for (k, (d, p)) in parts.iter().enumerate() {
                let zr = I::new(
                    if k > 0 { gaps[k - 1].band[0].min(gaps[k - 1].root.lo) } else { p.z[0].max(from) },
                    if k < gaps.len() { gaps[k].band[1].max(gaps[k].root.hi).min(top) } else { p.z[1] },
                );
                for i in 0..d.n() {
                    let c = d.at(i, zr);
                    for k in 0..2 {
                        lo[k] = lo[k].min(c[k].lo);
                        hi[k] = hi[k].max(c[k].hi);
                    }
                }
            }
            built.push((down, d0r, rounds, parts, gaps, kites, top));
        }
        if !(lo.iter().chain(&hi).all(|v| v.is_finite())) {
            return None;
        }
        // Turning lines read x and y clamped to the box grown by its size.
        let size = add_o(I::point(hi[0] - lo[0]), I::point(hi[1] - lo[1])).hi;
        let blo = [sub_o(I::point(lo[0]), I::point(size)).lo, sub_o(I::point(lo[1]), I::point(size)).lo];
        let bhi = [add_o(I::point(hi[0]), I::point(size)).hi, add_o(I::point(hi[1]), I::point(size)).hi];
        let clamp = |d: &mut Doc, v: NodeId, lo: f64, hi: f64| {
            let (lo, hi) = (d.g.constant(lo), d.g.constant(hi));
            let m = d.g.min(v, hi);
            d.g.max(m, lo)
        };
        let zc = clamp(self, z, z0, z1);
        let (xc, yc) = (clamp(self, x, blo[0], bhi[0]), clamp(self, y, blo[1], bhi[1]));
        let mut sides = Vec::new();
        for (down, d0, rounds, parts, gaps, kites, top) in &built {
            let w = if *down { self.g.neg(zc) } else { zc };
            // With corners on curves, its own heights from 0 (a range across
            // 0 reads the drawn profile there, not a path past its proof).
            let w = if d0.curve.iter().any(Option::is_some) {
                let z0 = self.g.constant(0.0);
                self.g.max(w, z0)
            } else {
                w
            };
            let mut planes: HashMap<(usize, usize, bool), NodeId> = HashMap::new();
            let mut fields = Vec::with_capacity(parts.len());
            // Each part's heights, with the bands either side.
            let zrs: Vec<I> = (0..parts.len())
                .map(|k| {
                    I::new(
                        if k > 0 { gaps[k - 1].band[0].min(gaps[k - 1].root.lo) } else { parts[k].1.z[0] },
                        if k < gaps.len() { gaps[k].band[1].max(gaps[k].root.hi).min(*top) } else { parts[k].1.z[1] },
                    )
                })
                .collect();
            // The arcs read at their chords' leaves in every part (the same
            // profile throughout: no edge trimmed), the rest within kites.
            let leafed: Vec<Vec<Leafed>> = kites
                .iter()
                .filter_map(|k| parts.iter().zip(&zrs).map(|((d, p), &zr)| if d.n() == d0.n() { d.at_leaf(&p.tree, k.i, k.t, zr) } else { None }).collect::<Option<Vec<Leafed>>>())
                .collect();
            let placed: Vec<usize> = leafed.iter().map(|l| l[0].i).collect();
            for (k, (d, p)) in parts.iter().enumerate() {
                let here: Vec<Leafed> = leafed.iter().map(|l| l[k].clone()).collect();
                let f = self.lower_drafted(d, &p.tree, zrs[k], [blo, bhi], [x, y, w, xc, yc], &mut planes, &here)?;
                fields.push(f);
            }
            // Below each change its lower tree, above it its upper.
            let mut f = *fields.last().expect("a part");
            for k in (0..gaps.len()).rev() {
                let Band { root, id, .. } = gaps[k];
                let c = self.g.within(root.lo, root.hi, root.mid().clamp(root.lo, root.hi), id);
                let cond = self.g.sub(w, c);
                f = self.g.select(cond, fields[k], f);
            }
            let f = self.drafted_arcs(d0, kites, &placed, f, [x, y, w, xc, yc], [blo, bhi], I::new(0.0, *top));
            sides.push((*down, self.drafted_rounds(d0, rounds, f, [x, y, w])));
        }
        // Below z = 0 the mirrored side, above it the drawn.
        let f = match sides.as_slice() {
            [(_, f)] => *f,
            [(false, up), (true, down)] => self.g.select(zc, *down, *up),
            _ => return None,
        };
        let mut terms = vec![f];
        for (v, k) in [(x, 0), (y, 1)] {
            terms.push(self.axis_face(v, hi[k], true));
            terms.push(self.axis_face(v, lo[k], false));
        }
        Some(self.reduce(BinOp::Max, &terms))
    }

    fn lower_drafted(&mut self, d: &Drafted, t: &Tr, zr: I, b: [P; 2], v: [NodeId; 5], planes: &mut HashMap<(usize, usize, bool), NodeId>, placed: &[Leafed]) -> Option<NodeId> {
        match t {
            Tr::Max(ts) => {
                let mut ids = Vec::with_capacity(ts.len());
                for t in ts {
                    let id = self.lower_drafted(d, t, zr, b, v, planes, placed)?;
                    if !ids.contains(&id) {
                        ids.push(id);
                    }
                }
                Some(self.reduce(BinOp::Max, &ids))
            }
            Tr::Neg(a, lid) => {
                let mut a = self.lower_drafted(d, a, zr, b, v, planes, placed)?;
                // Under a line giving way to an arc's segment: it leaves the
                // segment out.
                for l in placed.iter().filter(|l| l.relax.contains(lid)) {
                    let seg = self.arc_segment(d, l, zr, b, v, planes)?;
                    let ns = self.g.neg(seg);
                    a = self.g.max(a, ns);
                }
                Some(self.g.neg(a))
            }
            &Tr::Line(u, w) => {
                let n = d.n();
                // Arcs read at leaves (`Drafted::at_leaf`): at its chord's,
                // the line or its segment (out of the piece), or the line
                // less it (into the piece); at a line giving way to one, the
                // line or the segment (its cone and beyond its chord).
                if placed.iter().any(|l| l.dir == (u, w) || l.relax.contains(&(u, w))) {
                    let mut h = self.lower_drafted(d, t, zr, b, v, planes, &[])?;
                    // Its own arc first, then those it gives way to (in the
                    // other order its arc's disc would take theirs back).
                    for l in placed.iter().filter(|l| l.dir == (u, w)) {
                        let cone = self.arc_cone(d, l.i, [v[0], v[1], v[2]]);
                        h = if l.out {
                            self.g.min(h, cone)
                        } else {
                            let nc = self.g.neg(cone);
                            self.g.max(h, nc)
                        };
                    }
                    for l in placed.iter().filter(|l| l.relax.contains(&(u, w))) {
                        let seg = self.arc_segment(d, l, zr, b, v, planes)?;
                        h = self.g.min(h, seg);
                    }
                    return Some(h);
                }
                // A wall's own edge, or two corners on one plane: the plane,
                // its sign by which way the line runs along it.
                // (An arc's chord is no wall's plane.)
                let plane = if w == (u + 1) % n && d.walls[u].arc.is_none() {
                    Some((u, true))
                } else if u == (w + 1) % n && d.walls[w].arc.is_none() {
                    Some((w, false))
                } else if let Some(k) = d.shared(u, w) {
                    Some((k, d.wall_q(u, w, k).over(zr)? > 0))
                } else {
                    None
                };
                if let Some((k, out)) = plane {
                    let k0 = d.walls[k].class;
                    let rep = (0..n).find(|&j| d.walls[j].class == k0).expect("a wall of its class");
                    // (Its class names the plane in any profile trimmed from one.)
                    if let Some(&id) = planes.get(&(k0, 0, out)) {
                        return Some(id);
                    }
                    let id = self.wall_plane(d, &d.walls[rep], v[0], v[1], v[2]);
                    let id = if out { id } else { self.g.neg(id) };
                    planes.insert((k0, 0, out), id);
                    return Some(id);
                }
                if d.curved(&[u, w]) {
                    let (cu, cw) = (self.corner_xy(d, u, v[2]), self.corner_xy(d, w, v[2]));
                    let (mu, mw) = (d.motion(u), d.motion(w));
                    return Some(self.moving_line([cu, cw], [mu, mw], zr, b, [v[3], v[4]]));
                }
                Some(self.turning_line(d, u, w, zr, b, [v[2], v[3], v[4]]))
            }
        }
    }

    /// The profile f with its arcs' segments: each bulging out min(f, K),
    /// each bulging in max(f, -K), K its cone within its kite (the cone's
    /// field and the kite's four sides', each positive outside): the kites
    /// are apart, so their order does not matter.
    fn drafted_arcs(&mut self, d: &Drafted, kites: &[ArcKite], placed: &[usize], f: NodeId, [x, y, z, xc, yc]: [NodeId; 5], b: [P; 2], zr: I) -> NodeId {
        let n = d.n();
        let mut f = f;
        for k in kites.iter().filter(|k| !placed.contains(&k.i)) {
            let w = d.walls[k.i].arc.expect("an arc");
            let cone = self.arc_cone(d, k.i, [x, y, z]);
            // The kite, counter-clockwise: apex right of the chord where the
            // arc bulges out.
            let (ca, cq) = (Kp::C(k.i), Kp::C((k.i + 1) % n));
            let poly = if w.sigma > 0.0 { [ca, k.t, cq, k.t2] } else { [ca, k.t2, cq, k.t] };
            let mut terms = vec![cone];
            for e in 0..4 {
                let (p0, p1) = (poly[e], poly[(e + 1) % 4]);
                let pos = |dc: &mut Doc, p: Kp| match p {
                    Kp::C(i) => dc.corner_xy(d, i, z),
                    Kp::M(q, v) => [0, 1].map(|c| {
                        let (x, vc) = (dc.g.constant(q[c]), dc.g.constant(v[c]));
                        let vz = dc.g.mul(vc, z);
                        dc.g.add(x, vz)
                    }),
                };
                let (c0, c1) = (pos(self, p0), pos(self, p1));
                terms.push(self.moving_line([c0, c1], [d.kp(p0), d.kp(p1)], zr, b, [xc, yc]));
            }
            let kd = self.reduce(BinOp::Max, &terms);
            f = if w.sigma > 0.0 {
                self.g.min(f, kd)
            } else {
                let nk = self.g.neg(kd);
                self.g.max(f, nk)
            };
        }
        f
    }

    /// A leaf-read arc's segment: its cone and its chord's side.
    fn arc_segment(&mut self, d: &Drafted, l: &Leafed, zr: I, b: [P; 2], v: [NodeId; 5], planes: &mut HashMap<(usize, usize, bool), NodeId>) -> Option<NodeId> {
        let cone = self.arc_cone(d, l.i, [v[0], v[1], v[2]]);
        let c = self.lower_drafted(d, &Tr::Line(l.dir.0, l.dir.1), zr, b, v, planes, &[])?;
        let side = if l.out { self.g.neg(c) } else { c };
        Some(self.g.max(cone, side))
    }

    /// Arc i's cone, |p - C| - (R0 + sigma s z), C = a - G, scaled to slope
    /// 1: negative inside its circle at height z.
    fn arc_cone(&mut self, d: &Drafted, i: usize, [x, y, z]: [NodeId; 3]) -> NodeId {
        {
            let w = d.walls[i].arc.expect("an arc");
            let a0 = d.p[i];
            let id = d.walls[i].id;
            let gx = self.g.within(w.g[0].lo, w.g[0].hi, w.g[0].mid().clamp(w.g[0].lo, w.g[0].hi), d.ident("arc centre", &[id, 0]));
            let gy = self.g.within(w.g[1].lo, w.g[1].hi, w.g[1].mid().clamp(w.g[1].lo, w.g[1].hi), d.ident("arc centre", &[id, 1]));
            let r0 = self.g.within(w.r0.lo, w.r0.hi, w.r0.mid().clamp(w.r0.lo, w.r0.hi), d.ident("arc radius", &[id]));
            let (ax, ay) = (self.g.constant(a0[0]), self.g.constant(a0[1]));
            let (px, py) = (self.g.sub(x, ax), self.g.sub(y, ay));
            let (dx, dy) = (self.g.add(px, gx), self.g.add(py, gy));
            let (sx, sy) = (self.g.un(UnOp::Square, dx), self.g.un(UnOp::Square, dy));
            let q = self.g.add(sx, sy);
            let dist = self.g.un(UnOp::Sqrt, q);
            let rate = self.g.constant(w.sigma * w.s);
            let rz = self.g.mul(rate, z);
            let r = self.g.add(r0, rz);
            let cone = self.g.sub(dist, r);
            let kk = self.g.constant(iv::inv_norm_down(&[1.0, w.s.abs()]));
            let cone = self.g.mul(kk, cone);
            let one = self.g.constant(1.0);
            self.g.lipschitz(cone, one)
        }
    }

    /// Corner i's position as the graph reads it at height z: p + o + z v,
    /// or on its cone (`cone_xy`).
    fn corner_xy(&mut self, d: &Drafted, i: usize, z: NodeId) -> [NodeId; 2] {
        match d.curve[i].as_ref().map(|c| c.meet.clone()) {
            Some(Meet::Line(cc)) => return self.cone_xy(d, i, &cc, z),
            Some(Meet::Cone(cc)) => return self.cones_xy(d, i, &cc, z),
            None => {}
        }
        let c = d.cid(i);
        let pb = [d.p[i][0].to_bits(), d.p[i][1].to_bits()];
        [0, 1].map(|k| {
            let mut x = self.g.constant(d.p[i][k]);
            let o = d.o[i][k];
            if !(o.lo == 0.0 && o.hi == 0.0) {
                let o = self.g.within(o.lo, o.hi, o.mid().clamp(o.lo, o.hi), d.ident("offset", &[c[0], c[1], pb[0], pb[1], k as u64]));
                x = self.g.add(x, o);
            }
            let v = d.v[i][k];
            let v = self.g.within(v.lo, v.hi, v.mid().clamp(v.lo, v.hi), d.ident("speed", &[c[0], c[1], k as u64]));
            let vz = self.g.mul(v, z);
            self.g.add(x, vz)
        })
    }

    /// A corner on a cone (`ConeCorner`) as the graph reads it at height z:
    /// u + z w + tau(z) d, tau in the form that does not cancel, its
    /// numbers known within bounds (named by the corner's walls and role).
    fn cone_xy(&mut self, d: &Drafted, i: usize, cc: &ConeCorner, z: NodeId) -> [NodeId; 2] {
        let c = d.cid(i);
        let wi = |g: &mut Graph, x: I, role: u64| g.within(x.lo, x.hi, x.mid().clamp(x.lo, x.hi), d.ident("cone corner", &[c[0], c[1], role]));
        let (two, sg) = (I::point(2.0), I::point(cc.sigma));
        let div = |a: I, b: I| iv::bin(BinOp::Div, a, b).0;
        let one = self.g.constant(1.0);
        let tau = if cc.s[0] == cc.s[1] {
            let s = I::point(cc.s[0]);
            let a = wi(&mut self.g, mul_o(mul_o(two, sg), mul_o(s, cc.b)), 0);
            let dn = wi(&mut self.g, mul_o(cc.dd, cc.n), 1);
            let k = wi(&mut self.g, div(mul_o(mul_o(two, sg), s), cc.n), 2);
            let kz = self.g.mul(k, z);
            let f = self.g.add(one, kz);
            let sf = self.g.un(UnOp::Sqrt, f);
            let den0 = self.g.add(one, sf);
            let den = self.g.mul(dn, den0);
            let az = self.g.mul(a, z);
            self.g.bin(BinOp::Div, az, den)
        } else {
            let (k2, k1) = (wi(&mut self.g, cc.k[0], 3), wi(&mut self.g, cc.k[1], 4));
            let (b, b2, dd) = (wi(&mut self.g, cc.b, 5), wi(&mut self.g, mul_o(cc.b, cc.b), 6), wi(&mut self.g, cc.dd, 7));
            let zz = self.g.un(UnOp::Square, z);
            let (t2, t1) = (self.g.mul(k2, zz), self.g.mul(k1, z));
            let kk = self.g.add(t2, t1);
            let ddk = self.g.mul(dd, kk);
            let disc = self.g.sub(b2, ddk);
            let sq = self.g.un(UnOp::Sqrt, disc);
            let sq = if decided(cc.b) == Some(-1) { self.g.neg(sq) } else { sq };
            let den = self.g.add(b, sq);
            let num = self.g.neg(kk);
            self.g.bin(BinOp::Div, num, den)
        };
        [0, 1].map(|k| {
            let u = self.g.constant(cc.u[k]);
            let w = wi(&mut self.g, cc.w[k], 8 + k as u64);
            let dk = wi(&mut self.g, cc.d[k], 10 + k as u64);
            let wz = self.g.mul(w, z);
            let td = self.g.mul(tau, dk);
            let a = self.g.add(u, wz);
            self.g.add(a, td)
        })
    }

    /// A corner of two cones (`ConeCone`) as the graph reads it at height
    /// z: u + a e^ + tau f, tau in the factored form where it has one.
    fn cones_xy(&mut self, d: &Drafted, i: usize, cc: &ConeCone, z: NodeId) -> [NodeId; 2] {
        let c = d.cid(i);
        let wi = |g: &mut Graph, x: I, role: u64| g.within(x.lo, x.hi, x.mid().clamp(x.lo, x.hi), d.ident("cone cone corner", &[c[0], c[1], role]));
        let one = self.g.constant(1.0);
        let zz = self.g.un(UnOp::Square, z);
        let (a1, a2) = (wi(&mut self.g, cc.a[0], 0), wi(&mut self.g, cc.a[1], 1));
        let (t1, t2) = (self.g.mul(a1, z), self.g.mul(a2, zz));
        let a = self.g.add(t1, t2);
        let tau = match cc.fac {
            Some([r1, x0, q]) => {
                let four = I::point(4.0);
                let div = |a: I, b: I| iv::bin(BinOp::Div, a, b).0;
                let cq = wi(&mut self.g, div(mul_o(mul_o(four, r1), cc.b), q), 2);
                let kq = wi(&mut self.g, div(mul_o(four, r1), q), 3);
                let x0 = wi(&mut self.g, x0, 4);
                let r1 = self.g.constant(r1.lo);
                let rz = self.g.mul(r1, z);
                let xr = self.g.add(x0, rz);
                let t = self.g.mul(z, xr);
                let kt = self.g.mul(kq, t);
                let ff = self.g.add(one, kt);
                let sf = self.g.un(UnOp::Sqrt, ff);
                let den = self.g.add(one, sf);
                let num = self.g.mul(cq, t);
                self.g.bin(BinOp::Div, num, den)
            }
            None => {
                let (b, b2) = (wi(&mut self.g, cc.b, 6), wi(&mut self.g, mul_o(cc.b, cc.b), 7));
                // K = z (k0 + z (k1 + z (k2 + z k3))).
                let mut kk = wi(&mut self.g, cc.k[3], 8);
                for (j, kj) in cc.k[..3].iter().enumerate().rev() {
                    let kj = wi(&mut self.g, *kj, 14 + j as u64);
                    let t = self.g.mul(kk, z);
                    kk = self.g.add(kj, t);
                }
                let kk = self.g.mul(kk, z);
                let disc = self.g.sub(b2, kk);
                let sq = self.g.un(UnOp::Sqrt, disc);
                let sq = if decided(cc.b) == Some(-1) { self.g.neg(sq) } else { sq };
                let den = self.g.add(b, sq);
                let num = self.g.neg(kk);
                self.g.bin(BinOp::Div, num, den)
            }
        };
        [0, 1].map(|k| {
            let u = self.g.constant(cc.u[k]);
            let e = wi(&mut self.g, cc.eh[k], 10 + k as u64);
            let fk = wi(&mut self.g, cc.f[k], 12 + k as u64);
            let ae = self.g.mul(a, e);
            let tf = self.g.mul(tau, fk);
            let x = self.g.add(u, ae);
            self.g.add(x, tf)
        })
    }

    /// The line through two moving points (graph positions c, motions m),
    /// positive right of the first to the second: (q - c0) x (c1 - c0) k, q
    /// (x, y) clamped to the box, k at most 1 / its gradient's length over
    /// heights zr (as `turning_line`, from the points' positions and
    /// velocities over zr: read in full over the box, or about its middle
    /// and faded beyond rho = |D| / |D'|, whichever stronger).
    fn moving_line(&mut self, c: [[NodeId; 2]; 2], m: [Rc<dyn Fn(I) -> ([I; 2], [I; 2])>; 2], zr: I, b: [P; 2], [xc, yc]: [NodeId; 2]) -> NodeId {
        use ink_core::interval::{add_up, mul_up};
        let dd = [self.g.sub(c[1][0], c[0][0]), self.g.sub(c[1][1], c[0][1])];
        let qq = [self.g.sub(xc, c[0][0]), self.g.sub(yc, c[0][1])];
        let (a, cr) = (self.g.mul(qq[0], dd[1]), self.g.mul(qq[1], dd[0]));
        let f = self.g.sub(a, cr);
        let len = |v: [I; 2]| iv::norm([mag(v[0]), mag(v[1]), 0.0]);
        let ((p0, v0), (p1, v1)) = (m[0](zr), m[1](zr));
        let (dl, d1l) = (len(isub(p1, p0)), len(isub(v1, v0)));
        let vm = len([mul_o(I::point(0.5), add_o(v0[0], v1[0])), mul_o(I::point(0.5), add_o(v0[1], v1[1]))]);
        let mid = |k: usize| mul_o(I::point(0.5), add_o(p0[k], p1[k]));
        let reach = len([sub_o(I::new(b[0][0], b[1][0]), mid(0)), sub_o(I::new(b[0][1], b[1][1]), mid(1))]);
        let one = self.g.constant(1.0);
        let boxed = iv::inv_norm_down(&[dl, add_up(mul_up(d1l, reach), mul_up(dl, vm))]);
        let rho = if d1l > 0.0 { dl / d1l } else { f64::INFINITY };
        let sat = if rho < reach { iv::inv_norm_down(&[dl, add_up(mul_up(2.0 * rho, d1l), mul_up(dl, vm))]) } else { 0.0 };
        if !(boxed.is_finite() || sat.is_finite()) {
            let k = self.g.constant(0.0);
            return self.g.mul(k, f);
        }
        if sat > boxed {
            let half = self.g.constant(0.5);
            let (hx, hy) = (self.g.mul(half, dd[0]), self.g.mul(half, dd[1]));
            let (mx, my) = (self.g.sub(qq[0], hx), self.g.sub(qq[1], hy));
            let (sx, sy) = (self.g.un(UnOp::Square, mx), self.g.un(UnOp::Square, my));
            let q2 = self.g.add(sx, sy);
            let r = self.g.un(UnOp::Sqrt, q2);
            let kr = self.g.constant(rho);
            let far = self.g.max(r, kr);
            let wgt = self.g.bin(BinOp::Div, kr, far);
            let f = self.g.mul(f, wgt);
            let k = self.g.constant(sat);
            let f = self.g.mul(k, f);
            return self.g.lipschitz(f, one);
        }
        let k = self.g.constant(boxed);
        let f = self.g.mul(k, f);
        self.g.lipschitz(f, one)
    }

    /// Wall w's plane, (m . (p - a) - s |m| z) k, k at most 1 / its
    /// gradient's length: positive outside, exact on the drawn edge at
    /// z = 0, its offset per unit z known within bounds.
    fn wall_plane(&mut self, d: &Drafted, w: &Wall, x: NodeId, y: NodeId, z: NodeId) -> NodeId {
        let (ux, uy) = (self.g.constant(w.a[0]), self.g.constant(w.a[1]));
        let (xu, yu) = (self.g.sub(x, ux), self.g.sub(y, uy));
        let mut acc: Option<NodeId> = None;
        for (var, e, sg) in [(xu, &w.d[1], 1.0), (yu, &w.d[0], -1.0)] {
            for &c in e.iter().rev() {
                let c = self.g.constant(sg * c);
                let term = self.g.mul(var, c);
                acc = Some(match acc {
                    None => term,
                    Some(a) => self.g.add(a, term),
                });
            }
        }
        let mut f = acc.expect("a wall has length");
        if !(w.sl.lo == 0.0 && w.sl.hi == 0.0) {
            let s = self.g.within(w.sl.lo, w.sl.hi, w.sl.mid().clamp(w.sl.lo, w.sl.hi), d.ident("wall's offset per height", &[w.id]));
            let sz = self.g.mul(s, z);
            f = self.g.sub(f, sz);
        }
        let size = |e: &E| ink_core::interval::sum_up(e.iter().map(|v| v.abs()));
        let k = ink_core::interval::inv_norm_down(&[size(&w.d[0]), size(&w.d[1]), mag(w.sl)]);
        let k = self.g.constant(k);
        let f = self.g.mul(k, f);
        let one = self.g.constant(1.0);
        self.g.lipschitz(f, one)
    }

    /// The line through corners u and w, positive right of u -> w:
    /// (c_w - c_u) x (q - c_u) k, q (x, y) clamped to the box, its constants
    /// within bounds, k at most 1 / its gradient's length over the box and
    /// heights zr.
    fn turning_line(&mut self, d: &Drafted, u: usize, w: usize, zr: I, b: [P; 2], [z, xc, yc]: [NodeId; 3]) -> NodeId {
        let (pu, pw) = (d.p[u], d.p[w]);
        let wi = |g: &mut Graph, i: I, id: u64| g.within(i.lo, i.hi, i.mid().clamp(i.lo, i.hi), id);
        let d1 = isub(d.v[w], d.v[u]);
        let (cu, cw) = (d.cid(u), d.cid(w));
        // D = (pw - pu) + z (Vw - Vu), each axis; q - c_u = (q - pu) - z Vu.
        let mut dd = [0; 2];
        let mut qq = [0; 2];
        for k in 0..2 {
            let e = diff(pw[k], pu[k]);
            let mut s: Option<NodeId> = None;
            for &c in e.iter().rev() {
                let c = self.g.constant(c);
                s = Some(match s {
                    None => c,
                    Some(a) => self.g.add(a, c),
                });
            }
            let dv = wi(&mut self.g, d1[k], d.ident("speed difference", &[cu[0], cu[1], cw[0], cw[1], k as u64]));
            let dz = self.g.mul(dv, z);
            dd[k] = match s {
                Some(s) => self.g.add(s, dz),
                None => dz,
            };
            // Offsets (where an edge was trimmed).
            let od = sub_o(d.o[w][k], d.o[u][k]);
            if !(od.lo == 0.0 && od.hi == 0.0) {
                let pb = |i: usize| [d.p[i][0].to_bits(), d.p[i][1].to_bits()];
                let od = wi(&mut self.g, od, d.ident("offset difference", &[cu[0], cu[1], pb(u)[0], pb(u)[1], cw[0], cw[1], pb(w)[0], pb(w)[1], k as u64]));
                dd[k] = self.g.add(dd[k], od);
            }
            let pk = self.g.constant(pu[k]);
            let mut q = self.g.sub(if k == 0 { xc } else { yc }, pk);
            if !(d.o[u][k].lo == 0.0 && d.o[u][k].hi == 0.0) {
                let ou = wi(&mut self.g, d.o[u][k], d.ident("offset", &[cu[0], cu[1], pu[0].to_bits(), pu[1].to_bits(), k as u64]));
                q = self.g.sub(q, ou);
            }
            let vu = wi(&mut self.g, d.v[u][k], d.ident("speed", &[cu[0], cu[1], k as u64]));
            let vz = self.g.mul(vu, z);
            qq[k] = self.g.sub(q, vz);
        }
        // Right of u -> w: (q - c_u) x D = qx Dy - qy Dx, which is
        // (q - m) x D about its middle m = c_u + D / 2.
        let (a, c) = (self.g.mul(qq[0], dd[1]), self.g.mul(qq[1], dd[0]));
        let f = self.g.sub(a, c);
        // Its gradient: across, |D|; along z, (q - m) x D' - Vm x D, at
        // most |D'| |q - m| + |D| |Vm| (Vm the middle's speed).
        use ink_core::interval::{add_up, mul_up};
        let len = |v: [I; 2]| iv::norm([mag(v[0]), mag(v[1]), 0.0]);
        let d0 = d.d0(w, u);
        let dz_i = [add_o(d0[0], mul_o(zr, d1[0])), add_o(d0[1], mul_o(zr, d1[1]))];
        let (dl, d1l) = (len(dz_i), len(d1));
        let vm = len([mul_o(I::point(0.5), add_o(d.v[u][0], d.v[w][0])), mul_o(I::point(0.5), add_o(d.v[u][1], d.v[w][1]))]);
        let (cu, cw) = (d.at(u, zr), d.at(w, zr));
        let mid = |k: usize| mul_o(I::point(0.5), add_o(cu[k], cw[k]));
        let reach = len([sub_o(I::new(b[0][0], b[1][0]), mid(0)), sub_o(I::new(b[0][1], b[1][1]), mid(1))]);
        // Read in full over the box; or within rho of m and as rho / r of
        // itself beyond (r = |q - m|: there it is rho |D| sin of its angle
        // to D, whose slope across is at most |D| rho / r and along z at
        // most 2 rho |D'| + |D| |Vm|). The larger scale.
        let boxed = iv::inv_norm_down(&[dl, add_up(mul_up(d1l, reach), mul_up(dl, vm))]);
        let rho = if d1l > 0.0 { dl / d1l } else { f64::INFINITY };
        let sat = if rho < reach { iv::inv_norm_down(&[dl, add_up(mul_up(2.0 * rho, d1l), mul_up(dl, vm))]) } else { 0.0 };
        let one = self.g.constant(1.0);
        if sat > boxed {
            let half = self.g.constant(0.5);
            let (hx, hy) = (self.g.mul(half, dd[0]), self.g.mul(half, dd[1]));
            let (mx, my) = (self.g.sub(qq[0], hx), self.g.sub(qq[1], hy));
            let (sx, sy) = (self.g.un(UnOp::Square, mx), self.g.un(UnOp::Square, my));
            let q2 = self.g.add(sx, sy);
            let r = self.g.un(UnOp::Sqrt, q2);
            let kr = self.g.constant(rho);
            let far = self.g.max(r, kr);
            let wgt = self.g.bin(BinOp::Div, kr, far);
            let f = self.g.mul(f, wgt);
            let k = self.g.constant(sat);
            let f = self.g.mul(k, f);
            return self.g.lipschitz(f, one);
        }
        let k = self.g.constant(boxed);
        let f = self.g.mul(k, f);
        self.g.lipschitz(f, one)
    }
}

/// A point moving linearly with height, p + z v, p and v within bounds.
#[derive(Clone, Copy)]
struct Mp {
    p: [I; 2],
    v: [I; 2],
}

impl Mp {
    fn plus(self, o: [I; 2]) -> Mp {
        Mp { p: [add_o(self.p[0], o[0]), add_o(self.p[1], o[1])], v: self.v }
    }
}

/// (b - a) x (c - a), of moving points.
fn orient_mp(a: Mp, b: Mp, c: Mp) -> Q {
    let (u0, w0, u1, w1) = (isub(b.p, a.p), isub(c.p, a.p), isub(b.v, a.v), isub(c.v, a.v));
    Q::poly([icross(u0, w0), add_o(icross(u0, w1), icross(u1, w0)), icross(u1, w1)], 2)
}

/// (b - a) . d.
fn along_mp(a: Mp, b: Mp, d: [I; 2]) -> Q {
    Q::poly([idot(isub(b.p, a.p), d), idot(isub(b.v, a.v), d), I::point(0.0)], 2)
}

/// Whether convex polygons A and B (counter-clockwise; a segment is two
/// points) are apart over heights z: an edge's line has the other strictly
/// outside it throughout.
fn apart_mp(a: &[Mp], b: &[Mp], z: I) -> bool {
    let sep = |x: &[Mp], y: &[Mp]| {
        let n = x.len();
        (0..n).any(|i| {
            let (p, q) = (x[i], x[(i + 1) % n]);
            y.iter().all(|&r| orient_mp(p, q, r).over(z) == Some(-1)) || (n == 2 && y.iter().all(|&r| orient_mp(q, p, r).over(z) == Some(-1)))
        })
    };
    sep(a, b) || sep(b, a)
}

/// Whether `ok` holds over heights z, halving them up to 64 parts.
fn throughout(z: I, depth: u32, ok: &dyn Fn(I) -> bool) -> bool {
    ok(z)
        || (depth > 0 && {
            let m = z.mid();
            m > z.lo && m < z.hi && throughout(I::new(z.lo, m), depth - 1, ok) && throughout(I::new(m, z.hi), depth - 1, ok)
        })
}

/// A round at a corner, moving with it: its centre corner + w, convex (in
/// the profile) or not, its triangle (the centre, and where the normals to
/// its walls through it meet the line m through corner - w / 8 square to w),
/// counter-clockwise.
struct Round {
    i: usize,
    r: f64,
    convex: bool,
    w: [I; 2],
    tri: [Mp; 3],
}

impl Drafted {
    fn corner(&self, i: usize) -> Mp {
        Mp { p: [0, 1].map(|k| add_o(I::point(self.p[i][k]), self.o[i][k])), v: self.v[i] }
    }

    /// Wall j's outward unit normal and its direction (unnormalised).
    fn normal(&self, j: usize) -> ([I; 2], [I; 2]) {
        let d = [tot(&self.walls[j].d[0]), tot(&self.walls[j].d[1])];
        let l = iv::un(UnOp::Sqrt, idot(d, d));
        ([iv::bin(BinOp::Div, d[1], l).0, ineg(iv::bin(BinOp::Div, d[0], l).0)], d)
    }

    /// Each round, moving with its corner, certified over the span: its
    /// corner keeps its turn, its tangent points stay inside their edges
    /// (past the next round's), its triangle apart from every edge but its
    /// own two and from every other round's. None where one is not.
    fn rounds(&self, z0: f64, z1: f64) -> Option<Vec<Round>> {
        let n = self.n();
        let span = I::new(z0, z1);
        let mut out: Vec<Round> = Vec::new();
        // Tangent points: on wall i - 1 and on wall i.
        let mut tangents: Vec<Option<[Mp; 2]>> = vec![None; n];
        for i in (0..n).filter(|&i| self.r[i] > 0.0) {
            let h = (i + n - 1) % n;
            let convex = self.orient_q(h, i, (i + 1) % n).s0 > 0;
            let k = if convex { self.r[i] } else { -self.r[i] };
            let ((n0, _), (n1, _)) = (self.normal(h), self.normal(i));
            // n0 . w = -k = n1 . w.
            let det = icross(n0, n1);
            decided(det)?;
            let kk = I::point(-k);
            let div = |x: I| iv::bin(BinOp::Div, x, det).0;
            let w = [div(sub_o(mul_o(kk, n1[1]), mul_o(kk, n0[1]))), div(sub_o(mul_o(n0[0], kk), mul_o(n1[0], kk)))];
            let c = self.corner(i).plus(w);
            let t = |nn: [I; 2]| c.plus([mul_o(I::point(k), nn[0]), mul_o(I::point(k), nn[1])]);
            tangents[i] = Some([t(n0), t(n1)]);
            // Where the normal through c along nn meets m (through corner -
            // w / 8, square to w): c + s nn, s = -(9 / 8) |w|^2 / (nn . w).
            let ww = idot(w, w);
            let meet = |nn: [I; 2]| {
                let s = iv::bin(BinOp::Div, mul_o(I::point(-1.125), ww), idot(nn, w)).0;
                c.plus([mul_o(s, nn[0]), mul_o(s, nn[1])])
            };
            let mut tri = [c, meet(n0), meet(n1)];
            match orient_mp(tri[0], tri[1], tri[2]).over(I::point(0.5 * (z0 + z1)))? {
                1 => {}
                _ => tri.swap(1, 2),
            }
            if !tri.iter().flat_map(|m| m.p.iter().chain(&m.v)).all(|v| v.lo.is_finite() && v.hi.is_finite()) {
                return None;
            }
            let turn = if convex { 1 } else { -1 };
            if !throughout(span, 6, &|z| self.orient_q(h, i, (i + 1) % n).over(z) == Some(turn)) {
                return None;
            }
            out.push(Round { i, r: self.r[i], convex, w, tri });
        }
        for rd in &out {
            let i = rd.i;
            let (h, j) = ((i + n - 1) % n, (i + 1) % n);
            let [t0, t1] = tangents[i].expect("a round's tangents");
            let (d0, d1) = (self.normal(h).1, self.normal(i).1);
            let from = tangents[h].map_or(self.corner(h), |t| t[1]);
            let to = tangents[j].map_or(self.corner(j), |t| t[0]);
            let inside = |z: I| {
                along_mp(from, t0, d0).over(z) == Some(1)
                    && along_mp(t0, self.corner(i), d0).over(z) == Some(1)
                    && along_mp(self.corner(i), t1, d1).over(z) == Some(1)
                    && along_mp(t1, to, d1).over(z) == Some(1)
            };
            if !throughout(span, 6, &inside) {
                return None;
            }
            for e in (0..n).filter(|&e| e != h && e != i) {
                let seg = [self.corner(e), self.corner((e + 1) % n)];
                if !throughout(span, 6, &|z| apart_mp(&rd.tri, &seg, z)) {
                    return None;
                }
            }
        }
        for (a, ra) in out.iter().enumerate() {
            for rb in &out[a + 1..] {
                if !throughout(span, 6, &|z| apart_mp(&ra.tri, &rb.tri, z)) {
                    return None;
                }
            }
        }
        Some(out)
    }
}

impl Doc {
    /// (q - (o + z v)) . dir k, k at most 1 / its gradient's length: o the
    /// drawn corner p plus an offset within bounds, dir within bounds.
    fn moving_plane(&mut self, p: P, off: [I; 2], v: [I; 2], dir: [I; 2], [x, y, z]: [NodeId; 3], ids: [[u64; 2]; 3]) -> NodeId {
        let mut acc: Option<NodeId> = None;
        for (k, var) in [(0, x), (1, y)] {
            let q = self.moving_axis(p[k], off[k], v[k], var, z, [ids[0][k], ids[1][k]]);
            let dk = self.g.within(dir[k].lo, dir[k].hi, dir[k].mid().clamp(dir[k].lo, dir[k].hi), ids[2][k]);
            let t = self.g.mul(q, dk);
            acc = Some(match acc {
                None => t,
                Some(a) => self.g.add(a, t),
            });
        }
        let along = mag(add_o(mul_o(dir[0], v[0]), mul_o(dir[1], v[1])));
        let k = iv::inv_norm_down(&[mag(dir[0]), mag(dir[1]), along]);
        let k = self.g.constant(k);
        let f = self.g.mul(k, acc.expect("two axes"));
        let one = self.g.constant(1.0);
        self.g.lipschitz(f, one)
    }

    /// var - (p + off + z v) along one axis (off's and v's identities).
    fn moving_axis(&mut self, p: f64, off: I, v: I, var: NodeId, z: NodeId, ids: [u64; 2]) -> NodeId {
        let wi = |g: &mut Graph, i: I, id: u64| g.within(i.lo, i.hi, i.mid().clamp(i.lo, i.hi), id);
        let pc = self.g.constant(p);
        let mut q = self.g.sub(var, pc);
        if !(off.lo == 0.0 && off.hi == 0.0) {
            let o = wi(&mut self.g, off, ids[0]);
            q = self.g.sub(q, o);
        }
        if !(v.lo == 0.0 && v.hi == 0.0) {
            let vv = wi(&mut self.g, v, ids[1]);
            let vz = self.g.mul(vv, z);
            q = self.g.sub(q, vz);
        }
        q
    }

    /// The sharp profile f with its rounds: each convex one max(f, min(D,
    /// -T)), each other min(f, max(-D, T)), D its disc (a cylinder leaning
    /// with its corner), T its triangle's field (positive outside: the
    /// normals through its centre, positive toward the corner, and m beyond
    /// the corner).
    fn drafted_rounds(&mut self, d: &Drafted, rounds: &[Round], f: NodeId, xyz: [NodeId; 3]) -> NodeId {
        let n = d.n();
        let (mut cut, mut fill) = (Vec::new(), Vec::new());
        for rd in rounds {
            let (i, h) = (rd.i, (rd.i + n - 1) % n);
            let (p, v) = (d.p[i], d.v[i]);
            let w = [add_o(d.o[i][0], rd.w[0]), add_o(d.o[i][1], rd.w[1])];
            let c = d.cid(i);
            let ax2 = |what: &str| [0, 1].map(|k| d.ident(what, &[c[0], c[1], p[0].to_bits(), p[1].to_bits(), rd.r.to_bits(), k]));
            let (centre, speed) = (ax2("round centre"), [0, 1].map(|k| d.ident("speed", &[c[0], c[1], k])));
            // The disc: |q - c(z)| - r, slope at most sqrt(1 + |v|^2).
            let ax = [0, 1].map(|k| self.moving_axis(p[k], w[k], v[k], xyz[k], xyz[2], [centre[k], speed[k]]));
            let (sx, sy) = (self.g.un(UnOp::Square, ax[0]), self.g.un(UnOp::Square, ax[1]));
            let q = self.g.add(sx, sy);
            let dist = self.g.un(UnOp::Sqrt, q);
            let r = self.g.constant(rd.r);
            let disc = self.g.sub(dist, r);
            let k = self.g.constant(iv::inv_norm_down(&[1.0, iv::norm([mag(v[0]), mag(v[1]), 0.0])]));
            let disc = self.g.mul(k, disc);
            let one = self.g.constant(1.0);
            let disc = self.g.lipschitz(disc, one);
            // Its normals, positive toward the corner, negated; m.
            let (d0, d1) = (d.normal(h).1, d.normal(i).1);
            let mut sides = Vec::with_capacity(3);
            let dirs = [(d0, d.ident("direction", &[d.walls[h].id, 0])), ([ineg(d1[0]), ineg(d1[1])], d.ident("direction", &[d.walls[i].id, 1]))];
            for (dir, did) in dirs {
                let s = self.moving_plane(p, w, v, dir, xyz, [centre, speed, [0, 1].map(|k| d.ident("axis", &[did, k]))]);
                sides.push(self.g.neg(s));
            }
            let mw = [ineg(rd.w[0]), ineg(rd.w[1])];
            let eighth = [add_o(d.o[i][0], mul_o(I::point(0.125), mw[0])), add_o(d.o[i][1], mul_o(I::point(0.125), mw[1]))];
            sides.push(self.moving_plane(p, eighth, v, mw, xyz, [ax2("round cut"), speed, ax2("round cut's normal")]));
            let tri = self.reduce(BinOp::Max, &sides);
            if rd.convex {
                let out = self.g.neg(tri);
                cut.push(self.g.min(disc, out));
            } else {
                let nd = self.g.neg(disc);
                fill.push(self.g.max(nd, tri));
            }
        }
        let mut f = f;
        if !cut.is_empty() {
            cut.push(f);
            f = self.reduce(BinOp::Max, &cut);
        }
        if !fill.is_empty() {
            fill.push(f);
            f = self.reduce(BinOp::Min, &fill);
        }
        f
    }
}

/// Whether a drafted profile is built exactly (here), not by the allowance.
#[cfg(test)]
pub(super) fn exactly(input: &[[f64; 4]], drafts: &[f64], z0: f64, z1: f64) -> bool {
    let slopes: Vec<f64> = drafts.iter().map(|a| a.to_radians().tan()).collect();
    let neg: Vec<f64> = slopes.iter().map(|s| -s).collect();
    let side = |s: &[f64], from: f64, top: f64| {
        Drafted::new(input, s)
            .is_some_and(|d| d.rounds(from, top).is_some() && (d.walls.iter().all(|w| w.arc.is_none()) || d.kites(top).is_some()) && half(d, top, &mut false).is_some())
    };
    (z1 <= 0.0 || side(&slopes, z0.max(0.0), z1)) && (z0 >= 0.0 || side(&neg, (-z1).max(0.0), -z0))
}

/// Where a line wall meets an arc's wall (a cone) at their drawn corner u,
/// moving: u + z w + tau(z) d, w = s_l m^ on the line moved out by its
/// draft (m^ its outward unit normal, d its direction), tau the root through
/// 0 of |d|^2 tau^2 + 2 B tau + K(z) = 0, B = d . G, G = u - C, and
/// K = z^2 (s_l^2 - s_a^2) + 2 sigma z (R0 (s_l - s_a) - s_l E), the cone's
/// radius R0 + sigma s_a z (sigma 1 where the profile is inside its circle).
/// No term of the circle's size cancels: E = R0 - sigma m^ . G, small at a
/// joint near tangent, is B^2 / (|d|^2 N), N = R0 + sigma m^ . G (|G| = R0
/// in real numbers is never formed). With one draft (s_l = s_a = s) the
/// discriminant factors: tau = 2 sigma s B z / (|d|^2 N (1 + sqrt F)),
/// F = 1 + 2 sigma s z / N, which holds while F > 0, however tangent the
/// joint (no division by B).
#[derive(Clone, Debug)]
struct ConeCorner {
    u: P,
    d: [I; 2],
    dd: I,
    w: [I; 2],
    b: I,
    n: I,
    s: [f64; 2],
    sigma: f64,
    k: [I; 2],
}

impl ConeCorner {
    /// From the corner u, the line's direction d (exact two-term parts) and
    /// outward unit normal m^, its slope, and the arc's G = u - C, R0
    /// (bounds), sigma and slope. None where N is not shown positive.
    fn new(u: P, d: [I; 2], mhat: [I; 2], s_l: f64, g: [I; 2], r0: I, sigma: f64, s_a: f64) -> Option<ConeCorner> {
        let dd = idot(d, d);
        let b = idot(d, g);
        let sg = I::point(sigma);
        let n = add_o(r0, mul_o(sg, idot(mhat, g)));
        if !(n.lo > 0.0 && dd.lo > 0.0) {
            return None;
        }
        let e = iv::bin(BinOp::Div, mul_o(b, b), mul_o(dd, n)).0;
        let (sl, sa) = (I::point(s_l), I::point(s_a));
        let k2 = sub_o(mul_o(sl, sl), mul_o(sa, sa));
        let k1 = mul_o(mul_o(I::point(2.0), sg), sub_o(mul_o(r0, sub_o(sl, sa)), mul_o(sl, e)));
        let w = [mul_o(sl, mhat[0]), mul_o(sl, mhat[1])];
        Some(ConeCorner { u, d, dd, w, b, n, s: [s_l, s_a], sigma, k: [k2, k1] })
    }

    /// Its position and velocity over heights z; unbounded where its root
    /// is not shown to exist there (the walls may part).
    fn at(&self, z: I) -> ([I; 2], [I; 2]) {
        let none = ([I::ENTIRE; 2], [I::ENTIRE; 2]);
        let (sg, two) = (I::point(self.sigma), I::point(2.0));
        let div = |a: I, b: I| iv::bin(BinOp::Div, a, b).0;
        let (tau, dtau) = if self.s[0] == self.s[1] {
            let s = I::point(self.s[0]);
            let f = add_o(I::point(1.0), div(mul_o(mul_o(two, sg), mul_o(s, z)), self.n));
            if !(f.lo > 0.0) {
                return none;
            }
            let sf = iv::un(UnOp::Sqrt, f);
            let num = mul_o(mul_o(two, sg), mul_o(s, self.b));
            let tau = div(mul_o(num, z), mul_o(mul_o(self.dd, self.n), add_o(I::point(1.0), sf)));
            let dtau = div(mul_o(mul_o(sg, s), self.b), mul_o(mul_o(self.dd, self.n), sf));
            (tau, dtau)
        } else {
            let [k2, k1] = self.k;
            let kz = add_o(mul_o(k2, mul_o(z, z)), mul_o(k1, z));
            let dk = add_o(mul_o(mul_o(two, k2), z), k1);
            let disc = sub_o(mul_o(self.b, self.b), mul_o(self.dd, kz));
            let Some(sb) = decided(self.b) else { return none };
            if !(disc.lo > 0.0) {
                return none;
            }
            let sq = mul_o(I::point(sb as f64), iv::un(UnOp::Sqrt, disc));
            (div(ineg(kz), add_o(self.b, sq)), div(ineg(dk), mul_o(two, sq)))
        };
        let pos = [0, 1].map(|k| add_o(add_o(I::point(self.u[k]), mul_o(z, self.w[k])), mul_o(tau, self.d[k])));
        let vel = [0, 1].map(|k| add_o(self.w[k], mul_o(dtau, self.d[k])));
        if pos.iter().chain(&vel).any(|v| !(v.lo.is_finite() && v.hi.is_finite())) {
            return none;
        }
        (pos, vel)
    }
}

#[cfg(test)]
mod cone_tests {
    use super::*;

    fn arc_g(u: P, q: P, b: f64) -> ([I; 2], I) {
        let c = super::super::exact::circle(u, q, b);
        ([ineg(I::new(c[0][0], c[0][1])), ineg(I::new(c[1][0], c[1][1]))], I::new(c[2][0], c[2][1]))
    }

    #[test]
    fn a_line_tangent_to_an_arc_stays_tangent_when_both_lean_alike() {
        // A line along y = 0 into a quarter circle about (0, 3), tangent at
        // the origin, both drafted 10 degrees: the line moves down by s z, the
        // circle grows by s z, so they touch at (0, -s z) at every height.
        let s = 10f64.to_radians().tan();
        let (g, r0) = arc_g([0.0, 0.0], [3.0, 3.0], (22.5f64).to_radians().tan());
        let cc = ConeCorner::new([0.0, 0.0], [I::point(5.0), I::point(0.0)], [I::point(0.0), I::point(-1.0)], s, g, r0, 1.0, s).unwrap();
        for z in [-1.0, -0.1, 0.0, 0.3, 2.0] {
            let (p, v) = cc.at(I::point(z));
            // (The bulge tan 22.5 deg is a float: tangent to its rounding.)
            assert!(p[0].lo.abs() < 1e-12 && p[0].hi.abs() < 1e-12, "{z}: {p:?}");
            assert!(p[1].lo <= -s * z + 1e-15 && -s * z - 1e-15 <= p[1].hi, "{z}: {p:?}");
            assert!((v[1].mid() + s).abs() < 1e-12, "{v:?}");
        }
    }

    #[test]
    fn a_corner_on_a_cone_is_where_its_walls_meet() {
        // A line along y = 0 into an arc of bulge 0.2 to (3, 3), the line
        // drafted 0.1 and the arc 0.3 (slopes): at each height the corner is
        // where the line, moved down 0.1 z, meets the circle grown by 0.3 z.
        let (u, q, bulge) = ([0.0, 0.0], [3.0, 3.0], 0.2);
        let (g, r0) = arc_g(u, q, bulge);
        let cc = ConeCorner::new(u, [I::point(5.0), I::point(0.0)], [I::point(0.0), I::point(-1.0)], 0.1, g, r0, 1.0, 0.3).unwrap();
        let c = super::super::arc_of(u, q, bulge).c;
        let r = super::super::arc_of(u, q, bulge).r;
        for z in [-0.5, 0.0, 0.25, 1.0] {
            let (p, v) = cc.at(I::point(z));
            let y = -0.1 * z;
            let rr = r + 0.3 * z;
            let h = (rr * rr - (y - c[1]).powi(2)).sqrt();
            // The meeting nearer the drawn corner.
            let x = if (c[0] - h).abs() < (c[0] + h).abs() { c[0] - h } else { c[0] + h };
            assert!((p[0].mid() - x).abs() < 1e-12 && (p[1].mid() - y).abs() < 1e-15, "{z}: {p:?} vs ({x}, {y})");
            // Its velocity is its path's slope.
            let h = 1e-6;
            let (p2, _) = cc.at(I::point(z + h));
            assert!(((p2[0].mid() - p[0].mid()) / h - v[0].mid()).abs() < 1e-5, "{z}: {v:?}");
        }
    }
}

#[cfg(test)]
mod arc_tree_tests {
    use super::*;

    #[test]
    fn a_drafted_profile_with_an_arc_builds_its_tree_on_its_moving_corners() {
        // A rectangle whose bottom edge is an arc bulging out or in, drafted
        // alike or each wall its own way: its corners beside the arc move on
        // curves, and its tree is built over the span on decisions read from
        // their paths (the convex shortcut, the meet of its walls, is not
        // taken: a small arc's disc does not hold the profile).
        for (b, dr) in [(0.2, [10.0; 4]), (-0.2, [10.0; 4]), (0.2, [5.0, 0.0, 20.0, 0.0]), (0.414, [3.0; 4])] {
            let p = [[0.0, 0.0, b, 0.0], [10.0, 0.0, 0.0, 0.0], [10.0, 6.0, 0.0, 0.0], [0.0, 6.0, 0.0, 0.0]];
            let s: Vec<f64> = dr.iter().map(|a: &f64| a.to_radians().tan()).collect();
            let d = Drafted::new(&p, &s).expect("built");
            assert!(d.curve.iter().filter(|c| c.is_some()).count() == 2, "{b}");
            let (parts, _) = half(d, 2.0, &mut false).expect("its tree");
            assert!(parts.iter().all(|(_, p)| !p.preds.is_empty()), "{b} {dr:?}");
        }
    }
}

/// An arc read at its chord's leaf (`Drafted::at_leaf`): its chord as its
/// piece runs it, whether its segment is out of the piece, and the lines
/// that give way to it.
#[derive(Clone, Debug)]
struct Leafed {
    i: usize,
    dir: (usize, usize),
    out: bool,
    relax: Vec<(usize, usize)>,
}

/// A point a kite reads: a corner, or a float point moving on a line, q +
/// z v (its place at z = 0 and its velocity, floats), so a kite follows
/// its arc however far it moves.
#[derive(Clone, Copy, PartialEq)]
enum Kp {
    C(usize),
    M(P, P),
}

/// An arc's kite, moving with its ends: corners i and i + 1 (the arc's
/// chord), the apex t beyond where its end tangents meet and the point t2
/// across its chord. Its disc within it is its segment.
#[derive(Clone, Copy)]
struct ArcKite {
    i: usize,
    t: Kp,
    t2: Kp,
}

impl Drafted {
    fn kp(&self, p: Kp) -> Rc<dyn Fn(I) -> ([I; 2], [I; 2])> {
        match p {
            Kp::C(i) => self.motion(i),
            Kp::M(q, v) => {
                let (q, v) = ([I::point(q[0]), I::point(q[1])], [I::point(v[0]), I::point(v[1])]);
                Rc::new(move |z: I| ([add_o(q[0], mul_o(z, v[0])), add_o(q[1], mul_o(z, v[1]))], v))
            }
        }
    }

    /// Arc i's chord middle and where its end tangents meet, at heights 0
    /// and top (floats: a guess the proofs then check).
    fn kite_guide(&self, i: usize, top: f64) -> Option<[[P; 2]; 2]> {
        let n = self.n();
        let w = self.walls[i].arc.expect("an arc");
        let c = [self.p[i][0] - w.g[0].mid(), self.p[i][1] - w.g[1].mid()];
        let at = |z: f64| -> Option<[P; 2]> {
            let (a, q) = (self.motion(i)(I::point(z)).0, self.motion((i + 1) % n)(I::point(z)).0);
            let (a, q) = ([a[0].mid(), a[1].mid()], [q[0].mid(), q[1].mid()]);
            let (na, nq) = ([a[0] - c[0], a[1] - c[1]], [q[0] - c[0], q[1] - c[1]]);
            // x . na = a . na, x . nq = q . nq.
            let det = na[0] * nq[1] - na[1] * nq[0];
            let (ra, rq) = (a[0] * na[0] + a[1] * na[1], q[0] * nq[0] + q[1] * nq[1]);
            let x = [(ra * nq[1] - rq * na[1]) / det, (na[0] * rq - nq[0] * ra) / det];
            let m = [0.5 * (a[0] + q[0]), 0.5 * (a[1] + q[1])];
            (m.iter().chain(&x).all(|v| v.is_finite())).then_some([m, x])
        };
        Some([at(0.0)?, at(top)?])
    }

    /// The side of a -> b point c lies on over heights z (centred form), or
    /// 0 where c is a or b.
    fn side(&self, a: Kp, b: Kp, c: Kp, z: I) -> Option<i8> {
        if c == a || c == b {
            return Some(0);
        }
        let (fa, fb, fc) = (self.kp(a), self.kp(b), self.kp(c));
        let m = I::point(z.mid());
        let ((pa, _), (pb, _), (pc, _)) = (fa(m), fb(m), fc(m));
        let ((qa, va), (qb, vb), (qc, vc)) = (fa(z), fb(z), fc(z));
        let at = icross(isub(pb, pa), isub(pc, pa));
        let slope = add_o(icross(isub(vb, va), isub(qc, qa)), icross(isub(qb, qa), isub(vc, va)));
        decided(add_o(at, mul_o(sub_o(z, m), slope)))
    }

    /// Whether convex polygons x and y (counter-clockwise; a segment is two
    /// points) have no interior point in common over heights z: an edge's
    /// line has the other on its closed outer side throughout.
    fn apart_k(&self, x: &[Kp], y: &[Kp], z: I) -> bool {
        let sep = |x: &[Kp], y: &[Kp]| {
            let n = x.len();
            (0..n).any(|i| {
                let (p, q) = (x[i], x[(i + 1) % n]);
                y.iter().all(|&r| matches!(self.side(p, q, r, z), Some(0) | Some(-1))) || (n == 2 && y.iter().all(|&r| matches!(self.side(q, p, r, z), Some(0) | Some(-1))))
            })
        };
        sep(x, y) || sep(y, x)
    }

    /// Each arc's kite (`ArcKite`), certified over heights 0..top: its apex
    /// on the arc's side beyond both end tangents (the triangle holds the
    /// segment), its centre across its chord (under half a turn), its
    /// triangles apart from every other edge and kite (closed: a shared
    /// corner touches). Each point moves on a line, from its place at z =
    /// 0 to its place at top, guided by where the end tangents meet there:
    /// the apex pushed out past it until it holds, the point across halved
    /// toward the chord until apart. None where one is not.
    fn kites(&self, top: f64) -> Option<Vec<ArcKite>> {
        let n = self.n();
        let span = I::new(0.0, top);
        let arcs: Vec<usize> = (0..n).filter(|&i| self.walls[i].arc.is_some()).collect();
        let mut out: Vec<ArcKite> = Vec::new();
        // Each kite's guide (its apex's) and how far across its point is.
        let mut across: Vec<([[P; 2]; 2], f64)> = Vec::new();
        // A point f(m, x) at heights 0 and top, moving on the line between.
        let moving = |g: [[P; 2]; 2], f: &dyn Fn(P, P) -> P| {
            let (p0, p1) = (f(g[0][0], g[0][1]), f(g[1][0], g[1][1]));
            Kp::M(p0, if top > 0.0 { [(p1[0] - p0[0]) / top, (p1[1] - p0[1]) / top] } else { [0.0; 2] })
        };
        let out_by = |l: f64| move |m: P, x: P| [m[0] + l * (x[0] - m[0]), m[1] + l * (x[1] - m[1])];
        for &i in &arcs {
            let w = self.walls[i].arc.expect("an arc");
            let a0 = self.p[i];
            let g = self.kite_guide(i, top)?;
            let arc_side: i8 = if w.sigma > 0.0 { -1 } else { 1 }; // right of a -> q where it bulges out
            let (ca, cq) = (Kp::C(i), Kp::C((i + 1) % n));
            let ctr = [add_o(I::point(a0[0]), ineg(w.g[0])), add_o(I::point(a0[1]), ineg(w.g[1]))];
            // (T - c) . (c - C) over z, for an end corner c.
            let beyond = |t: Kp, c: usize, z: I| {
                let (fc, ft) = (self.motion(c), self.kp(t));
                let m = I::point(z.mid());
                let ((pz, vz), (tz, vt)) = (fc(z), ft(z));
                let ((pm, _), (tm, _)) = (fc(m), ft(m));
                let at = idot(isub(tm, pm), isub(pm, ctr));
                let slope = add_o(idot(isub(vt, vz), isub(pz, ctr)), idot(isub(tz, pz), vz));
                decided(add_o(at, mul_o(sub_o(z, m), slope))) == Some(1)
            };
            let holds = |t: Kp| {
                throughout(span, 8, &|z| {
                    self.side(ca, cq, t, z) == Some(arc_side) && beyond(t, i, z) && beyond(t, (i + 1) % n, z) && {
                        // Under half a turn: the centre across the chord.
                        let (fa, fq) = (self.motion(i), self.motion((i + 1) % n));
                        let m = I::point(z.mid());
                        let ((pa, _), (pq, _)) = (fa(m), fq(m));
                        let ((qa, va), (qq, vq)) = (fa(z), fq(z));
                        let at = icross(isub(pq, pa), isub(ctr, pa));
                        let slope = add_o(icross(isub(vq, va), isub(ctr, qa)), icross(isub(qq, qa), ineg2(va)));
                        decided(add_o(at, mul_o(sub_o(z, m), slope))) == Some(-arc_side)
                    }
                })
            };
            // The apex past where the end tangents meet, pushed out until it
            // holds.
            let l = (0..8).map(|j| 1.25 * 2f64.powi(j)).find(|&l| holds(moving(g, &out_by(l))))?;
            let t = moving(g, &out_by(l));
            // Its triangle (counter-clockwise) apart from every other edge.
            let tri = if arc_side < 0 { vec![ca, t, cq] } else { vec![ca, cq, t] };
            for e in (0..n).filter(|&e| e != i) {
                let seg = [Kp::C(e), Kp::C((e + 1) % n)];
                if !throughout(span, 8, &|z| self.apart_k(&tri, &seg, z)) {
                    return None;
                }
            }
            // The point across, toward the chord until its triangle is apart.
            let g = g.map(|[m, x]| [m, out_by(l)(m, x)]);
            let mut k = 1.0;
            let mut t2 = None;
            for _ in 0..30 {
                let p = moving(g, &out_by(-k));
                if self.across_holds(i, p, span) {
                    t2 = Some(p);
                    break;
                }
                k *= 0.5;
            }
            out.push(ArcKite { i, t, t2: t2? });
            across.push((g, k));
        }
        // Kites apart from each other: where two are not, both points across
        // are brought halfway to their chords, shown again to be across and
        // apart from the other edges, and tried again.
        let poly = |k: &ArcKite| {
            let (ca, cq) = (Kp::C(k.i), Kp::C((k.i + 1) % n));
            // Counter-clockwise: apex right of the chord where it bulges out.
            if self.walls[k.i].arc.expect("an arc").sigma > 0.0 { vec![ca, k.t, cq, k.t2] } else { vec![ca, k.t2, cq, k.t] }
        };
        let mut done = vec![false; out.len() * out.len()];
        for _ in 0..30 {
            let mut closer = vec![false; out.len()];
            for x in 0..out.len() {
                for y in x + 1..out.len() {
                    if done[x * out.len() + y] {
                        continue;
                    }
                    // (A kite is convex where its two triangles are.)
                    let (px, py) = (poly(&out[x]), poly(&out[y]));
                    if throughout(span, 8, &|z| {
                        self.apart_k(&px[..3], &py[..3], z)
                            && self.apart_k(&px[..3], &[py[2], py[3], py[0]], z)
                            && self.apart_k(&[px[2], px[3], px[0]], &py[..3], z)
                            && self.apart_k(&[px[2], px[3], px[0]], &[py[2], py[3], py[0]], z)
                    }) {
                        done[x * out.len() + y] = true;
                    } else {
                        closer[x] = true;
                        closer[y] = true;
                    }
                }
            }
            if !closer.contains(&true) {
                return Some(out);
            }
            for j in (0..out.len()).filter(|&j| closer[j]) {
                let (g, k) = &mut across[j];
                *k *= 0.5;
                let p = moving(*g, &out_by(-*k));
                if !self.across_holds(out[j].i, p, span) {
                    return None;
                }
                out[j].t2 = p;
                for x in 0..out.len() {
                    done[x.min(j) * out.len() + x.max(j)] = false;
                }
            }
        }
        None
    }
}

impl Drafted {
    /// How arc i is read at its chord's leaf in tree t over heights zr, if
    /// it can be: lane E's flat construction (`exact::piece`), its reading
    /// proven over the heights. Its segment lies out of its piece (the
    /// piece gains it) or into it (loses it); its triangle (chord and apex,
    /// which holds the segment) is proven on the closed inner side of each
    /// line of its piece and the pieces above, or else: past a line of a
    /// piece of the other parity than the segment's (out of a piece of
    /// depth k it is of parity k, into it k + 1) it lies where its side is
    /// read already; past one of the same parity that line gives way to it
    /// (the line or the segment) and each pocket under it leaves it out; a
    /// pocket's lid (its parent's line, of the pocket's parity) cannot give
    /// way, so it is not read so.
    fn at_leaf(&self, t: &Tr, i: usize, apex: Kp, zr: I) -> Option<Leafed> {
        let n = self.n();
        let q = (i + 1) % n;
        // The lines of the piece holding the chord and of those above it,
        // each with its piece's depth and whether it is a lid; and the
        // depth of the chord's piece.
        type L = ((usize, usize), usize, bool);
        fn find(t: &Tr, c: [(usize, usize); 2], depth: usize, anc: &mut Vec<L>) -> Option<usize> {
            let Tr::Max(ts) = t else { return None };
            let base = anc.len();
            anc.extend(ts.iter().filter_map(|t| if let Tr::Line(u, w) = t { Some(((*u, *w), depth, false)) } else { None }));
            if anc[base..].iter().any(|l| c.contains(&l.0)) {
                return Some(depth);
            }
            for t in ts {
                if let Tr::Neg(a, lid) = t {
                    anc.push((*lid, depth + 1, true));
                    if let Some(k) = find(a, c, depth + 1, anc) {
                        return Some(k);
                    }
                    anc.pop();
                }
            }
            anc.truncate(base);
            None
        }
        let mut lines: Vec<L> = Vec::new();
        let depth = find(t, [(i, q), (q, i)], 0, &mut lines)?;
        let dir = lines.iter().find(|l| !l.2 && l.1 == depth && (l.0 == (i, q) || l.0 == (q, i))).expect("its chord").0;
        let out = (dir.0 == i) == (self.walls[i].arc.expect("an arc").sigma > 0.0);
        let parity = (depth + usize::from(!out)) % 2;
        let tri = [Kp::C(i), apex, Kp::C(q)];
        let inner = |(u, w): (usize, usize)| throughout(zr, 8, &|z| tri.iter().all(|&p| matches!(self.side(Kp::C(u), Kp::C(w), p, z), Some(0) | Some(1))));
        let mut relax = Vec::new();
        for &(l, dl, lid) in lines.iter().filter(|l| l.2 || l.0 != dir) {
            if dl % 2 != parity || inner(l) {
                continue;
            }
            if lid {
                return None;
            }
            relax.push(l);
        }
        Some(Leafed { i, dir, out, relax })
    }

    /// Whether p lies across arc i's chord from the arc over heights z and
    /// its triangle with the chord is apart from every other edge.
    fn across_holds(&self, i: usize, p: Kp, span: I) -> bool {
        let n = self.n();
        let arc_side: i8 = if self.walls[i].arc.expect("an arc").sigma > 0.0 { -1 } else { 1 };
        let (ca, cq) = (Kp::C(i), Kp::C((i + 1) % n));
        let below = if arc_side < 0 { vec![ca, cq, p] } else { vec![ca, p, cq] };
        throughout(span, 8, &|z| self.side(ca, cq, p, z) == Some(-arc_side))
            && (0..n).filter(|&e| e != i).all(|e| throughout(span, 8, &|z| self.apart_k(&below, &[Kp::C(e), Kp::C((e + 1) % n)], z)))
    }
}

fn ineg2(v: [I; 2]) -> [I; 2] {
    [ineg(v[0]), ineg(v[1])]
}

/// Where two arcs' walls (cones) meet at their drawn corner u, moving:
/// u + a(z) e^ + tau(z) f, e^ the unit vector along e = G1 - G2 (Gk = u -
/// Ck, so from the second centre to the first) and f = e^ turned a quarter
/// left. Circle k grows at rate rk = sigma_k s_k: |delta + Gk|^2 = (Rk +
/// rk z)^2, that is |delta|^2 + 2 delta . Gk = hk = 2 rk Rk z + rk^2 z^2
/// (|Gk| = Rk is never formed). Their difference puts the corner on the
/// moving line a = (h1 - h2) / 2L (L = |e|); tau is the root through 0 of
/// tau^2 + 2 B tau + K = 0, B = f . G1 = (G1 x G2) / L, K = a^2 + 2 a e^.G1
/// - h1, real while D = B^2 - K > 0.
///
/// D is the half chord of the grown circles squared: (P^2 - L^2)(L^2 -
/// M^2) / 4L^2, P and M their radii's sum and difference. Where the rates
/// are equal (r2 = r1) M stays put, where opposite (r2 = -r1) P does: with
/// X the one that moves (X = X0 + 2 r1 z, X0 = R1 + eps R2, eps = r2 / r1),
/// D = B^2 (X^2 - L^2) / Q, Q = X0^2 - L^2 = 2 (eps R1 R2 + G1 . G2). So a
/// joint tangent in intent (B tiny, never 0 in the data) is shown real by
/// B's sign and F = 1 + 4 r1 z (X0 + r1 z) / Q > 0, no term cancelling:
/// tau = 4 r1 B z (X0 + r1 z) / (Q (1 + sqrt F)). Otherwise the plain form.
#[derive(Clone, Debug)]
struct ConeCone {
    u: P,
    eh: [I; 2],
    f: [I; 2],
    b: I,
    /// a = a1 z + a2 z^2.
    a: [I; 2],
    /// The plain form: K = z (k0 + k1 z + k2 z^2 + k3 z^3) (z taken out,
    /// so it is tight over a range of heights).
    k: [I; 4],
    /// Equal or opposite rates: r1, X0, Q.
    fac: Option<[I; 3]>,
}

impl ConeCone {
    fn new(u: P, g1: [I; 2], g2: [I; 2], r0: [I; 2], r: [f64; 2]) -> Option<ConeCone> {
        let e = isub(g1, g2);
        let el = iv::un(UnOp::Sqrt, idot(e, e));
        if !(el.lo > 0.0) {
            return None;
        }
        let div = |a: I, b: I| iv::bin(BinOp::Div, a, b).0;
        let eh = [div(e[0], el), div(e[1], el)];
        let f = [ineg(eh[1]), eh[0]];
        let b = div(sub_o(mul_o(g1[0], g2[1]), mul_o(g1[1], g2[0])), el);
        decided(b)?;
        let (r1, r2, two) = (I::point(r[0]), I::point(r[1]), I::point(2.0));
        let a1 = div(sub_o(mul_o(r1, r0[0]), mul_o(r2, r0[1])), el);
        let a2 = div(sub_o(mul_o(r1, r1), mul_o(r2, r2)), mul_o(two, el));
        let h1 = [mul_o(mul_o(two, r1), r0[0]), mul_o(r1, r1)];
        let fac = (r[0] != 0.0 && (r[1] == r[0] || r[1] == -r[0])).then(|| {
            let eps = I::point(r[1] / r[0]);
            let q = mul_o(two, add_o(mul_o(eps, mul_o(r0[0], r0[1])), idot(g1, g2)));
            [r1, add_o(r0[0], mul_o(eps, r0[1])), q]
        });
        if fac.is_some_and(|[_, _, q]| decided(q).is_none()) {
            return None;
        }
        // K = a^2 + 2 a p - h1, p = e^ . G1, a = z (a1 + a2 z).
        let p = idot(eh, g1);
        let k = [sub_o(mul_o(mul_o(two, a1), p), h1[0]), sub_o(add_o(mul_o(a1, a1), mul_o(mul_o(two, a2), p)), h1[1]), mul_o(mul_o(two, a1), a2), mul_o(a2, a2)];
        Some(ConeCone { u, eh, f, b, a: [a1, a2], k, fac })
    }

    /// Its position and velocity over heights z; unbounded where its root
    /// is not shown to exist there (the walls may part).
    fn at(&self, z: I) -> ([I; 2], [I; 2]) {
        let none = ([I::ENTIRE; 2], [I::ENTIRE; 2]);
        let (two, four) = (I::point(2.0), I::point(4.0));
        let div = |a: I, b: I| iv::bin(BinOp::Div, a, b).0;
        let Some(sb) = decided(self.b) else { return none };
        let [a1, a2] = self.a;
        let a = add_o(mul_o(a1, z), mul_o(a2, mul_o(z, z)));
        let da = add_o(a1, mul_o(mul_o(two, a2), z));
        let (tau, dtau) = match self.fac {
            Some([r1, x0, q]) => {
                // t = z (X0 + r1 z), F = 1 + 4 r1 t / Q.
                let t = mul_o(z, add_o(x0, mul_o(r1, z)));
                let dt = add_o(x0, mul_o(mul_o(two, r1), z));
                let ff = add_o(I::point(1.0), div(mul_o(mul_o(four, r1), t), q));
                if !(ff.lo > 0.0) {
                    return none;
                }
                let sf = iv::un(UnOp::Sqrt, ff);
                let c = div(mul_o(mul_o(four, r1), self.b), q);
                (div(mul_o(c, t), add_o(I::point(1.0), sf)), div(mul_o(mul_o(I::point(0.5), c), dt), sf))
            }
            None => {
                let [k0, k1, k2, k3] = self.k;
                let k = mul_o(z, add_o(k0, mul_o(z, add_o(k1, mul_o(z, add_o(k2, mul_o(z, k3)))))));
                let dk = add_o(k0, mul_o(z, add_o(mul_o(two, k1), mul_o(z, add_o(mul_o(I::point(3.0), k2), mul_o(mul_o(four, k3), z))))));
                let disc = sub_o(mul_o(self.b, self.b), k);
                if !(disc.lo > 0.0) {
                    return none;
                }
                let sq = mul_o(I::point(sb as f64), iv::un(UnOp::Sqrt, disc));
                (div(ineg(k), add_o(self.b, sq)), div(ineg(dk), mul_o(two, sq)))
            }
        };
        let pos = [0, 1].map(|c| add_o(add_o(I::point(self.u[c]), mul_o(a, self.eh[c])), mul_o(tau, self.f[c])));
        let vel = [0, 1].map(|c| add_o(mul_o(da, self.eh[c]), mul_o(dtau, self.f[c])));
        if pos.iter().chain(&vel).any(|v| !(v.lo.is_finite() && v.hi.is_finite())) {
            return none;
        }
        (pos, vel)
    }
}

#[cfg(test)]
mod cone_cone_tests {
    use super::*;

    /// Arc k's G at u and R0, from lane E's circle (drawn a -> q, bulge b).
    fn arc_g(a: P, q: P, b: f64, u: P) -> ([I; 2], I) {
        let c = super::super::exact::circle(a, q, b);
        let g = [ineg(I::new(c[0][0], c[0][1])), ineg(I::new(c[1][0], c[1][1]))];
        let d = [u[0] - a[0], u[1] - a[1]];
        ([add_o(g[0], I::point(d[0])), add_o(g[1], I::point(d[1]))], I::new(c[2][0], c[2][1]))
    }

    /// Where the circles about c1, c2 of radii r1, r2 meet, nearest u.
    fn meet(c1: P, r1: f64, c2: P, r2: f64, u: P) -> P {
        let (dx, dy) = (c2[0] - c1[0], c2[1] - c1[1]);
        let d = dx.hypot(dy);
        let al = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
        let h = (r1 * r1 - al * al).max(0.0).sqrt();
        let m = [c1[0] + dx / d * al, c1[1] + dy / d * al];
        let cs = [[m[0] - dy / d * h, m[1] + dx / d * h], [m[0] + dy / d * h, m[1] - dx / d * h]];
        let dist = |p: P| (p[0] - u[0]).hypot(p[1] - u[1]);
        if dist(cs[0]) < dist(cs[1]) { cs[0] } else { cs[1] }
    }

    fn check(q1: P, b1: f64, u: P, q2: P, b2: f64, r: [f64; 2], tol: f64) {
        let (g1, r1) = arc_g(q1, u, b1, u);
        let (g2, r2) = arc_g(u, q2, b2, u);
        let cc = ConeCone::new(u, g1, g2, [r1, r2], r).unwrap();
        let c1 = [u[0] - g1[0].mid(), u[1] - g1[1].mid()];
        let c2 = [u[0] - g2[0].mid(), u[1] - g2[1].mid()];
        for z in [-0.4, -0.01, 0.0, 0.3] {
            let (p, v) = cc.at(I::point(z));
            let x = meet(c1, r1.mid() + r[0] * z, c2, r2.mid() + r[1] * z, u);
            assert!((p[0].mid() - x[0]).abs() < tol && (p[1].mid() - x[1]).abs() < tol, "{z}: {p:?} vs {x:?}");
            assert!(p.iter().all(|c| c.hi - c.lo < 1e-9), "{z}: {p:?}");
            let hh = 1e-6;
            let (p2, _) = cc.at(I::point(z + hh));
            for k in 0..2 {
                assert!(((p2[k].mid() - p[k].mid()) / hh - v[k].mid()).abs() < 1e-4, "{z}: {v:?}");
            }
        }
    }

    #[test]
    fn a_corner_of_two_cones_is_where_their_circles_meet() {
        // Two arcs meeting at the origin at an angle: drafts unequal (the
        // plain form), equal, and opposite (the factored one).
        let (q1, u, q2) = ([-4.0, 1.0], [0.0, 0.0], [3.0, 2.0]);
        check(q1, -0.3, u, q2, 0.2, [0.1, 0.25], 1e-11);
        check(q1, -0.3, u, q2, 0.2, [0.1, 0.1], 1e-11);
        check(q1, -0.3, u, q2, 0.2, [0.1, -0.1], 1e-11);
    }

    #[test]
    fn two_arcs_tangent_in_intent_meet_at_every_height_with_equal_drafts() {
        // Two arcs tangent at the origin (centres on the y axis, one inside
        // the other's circle, and back to back), their tangency rounded off
        // in the data: the grown circles stay tangent, so the corner moves
        // radially, and the factored form shows it at every height where
        // the plain one cannot.
        let u = [0.0, 0.0];
        // Circle 1 about (0, 5), radius 5, through u and tangent to y = 0
        // there; circle 2 of radius 2 through u, its centre turned 1e-9 off
        // the y axis (above it, then below it, the profile outside it).
        let at = |c: P, r: f64, t: f64| [c[0] + r * t.sin(), c[1] - r * t.cos()];
        let q1 = at([0.0, 5.0], 5.0, -0.8);
        let (b1, b2) = ((0.8f64 / 4.0).tan(), (1.1f64 / 4.0).tan());
        let (g1, r1) = arc_g(q1, u, b1, u);
        let s = 0.05;
        for (side, r) in [(1.0, [s, s]), (-1.0, [s, -s])] {
            let phi = 1e-9f64;
            let c2 = [2.0 * phi.sin(), side * 2.0 * phi.cos()];
            // u's angle about c2, then on 1.1 rad, anticlockwise where the
            // bulge is positive.
            let tu = (-c2[1]).atan2(-c2[0]);
            let q2 = [c2[0] + 2.0 * (tu + side * 1.1).cos(), c2[1] + 2.0 * (tu + side * 1.1).sin()];
            let (g2, r2) = arc_g(u, q2, side * b2, u);
            let cc = ConeCone::new(u, g1, g2, [r1, r2], r).unwrap();
            let plain = ConeCone { fac: None, ..cc.clone() };
            let mut shown = 0;
            for z in [-1.0, -0.2, 0.5, 3.0] {
                let (p, _) = cc.at(I::point(z));
                assert!(p[0].hi - p[0].lo < 1e-9 && p[0].mid().abs() < 1e-7, "{side} {z}: {p:?}");
                assert!((p[1].mid() + s * z).abs() < 1e-7, "{side} {z}: {p:?}");
                shown += plain.at(I::point(z)).0[0].hi.is_finite() as usize;
            }
            assert!(shown < 4, "{side}: the plain form shows it too");
        }
    }
}
