//! `(polygon pts)`: a closed profile of lines and arcs in the xy plane, as
//! half-planes and discs, so interval arithmetic can prune it per region
//! like any min or max (a nearest-edge choice among edges cannot be
//! pruned where their distances' intervals overlap), at about 4
//! instructions an edge.
//!
//! The profile is its convex hull less the pockets between the two, each
//! pocket read the same way (a tree of convex pieces). The hull is the true
//! one, of the corners and of the arcs bulging out: each arc is split where
//! the hull touches it (a line tangent to it, in closed form: from a corner,
//! or common to two arcs), so the chord polygon's hull is the hull, and a
//! pocket's own hull the same way. A piece is the max of its edges'
//! half-planes; a pocket leaves out the hull edge over its mouth (its lid:
//! beyond it the parent is empty already), so an offset does not close the
//! mouth. An arc belongs to the piece its chord is an edge of: one bulging
//! out of the piece is its chord's half-plane or its disc; one bulging in
//! adds the disc's outside. Where an arc bulging out turns out past the
//! next hull edge and the boundary goes on tangent to it (an inflection: a
//! concave fillet's end, or a pocket starting along the tangent; no hull
//! follows there but through a pocket at every scale), that edge gives way
//! on the arc's side of the tangent: its half-plane or the tangent's. So
//! the sign is exact (a term is zero off the boundary only where another is
//! not), and the value is 1-Lipschitz and never more than the distance: the
//! distance along each edge, a bound past corners (an offset mitres those
//! between lines, but past a line's join with an arc an outward offset
//! keeps a sliver: an arc's chord reaches beyond the grown arc; std's
//! `polygon-distance` is the exact distance outside, whose offsets round
//! them). Rounded corners and drafts: walls.rs.
//!
//! Each arc is its chord with the segment between them added or taken
//! away, so arcs over a quarter turn are halved first, and arcs are halved
//! while another edge reaches into a segment or across a chord.

/// A field of x and y, to be lowered into the graph.
#[derive(Clone, Debug)]
pub(crate) enum T {
    /// a x + b y + c, and what line it is in the loop (for a drafted
    /// profile's moving leaves): through loop points a -> b (on its left),
    /// or tangent to edge e's arc at loop point p.
    Lin(f64, f64, f64, Option<Normal>),
    /// |(x, y) - (cx, cy)| - r, promised 1-Lipschitz, and the loop edge
    /// whose arc's circle it is.
    Disc(f64, f64, f64, Option<usize>),
    Max(Vec<T>),
    Neg(Box<T>),
    Min(Vec<T>),
}

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
fn unit(a: P) -> P {
    mul(a, 1.0 / norm(a))
}
/// The unit normal on the right of a -> b (outward, for a counter-clockwise loop).
fn right(a: P, b: P) -> P {
    let t = unit(sub(b, a));
    [t[1], -t[0]]
}

/// n . (x - p): the signed distance past the line through p with unit normal n.
fn half(n: P, p: P, line: [usize; 2]) -> T {
    T::Lin(n[0], n[1], -dot(n, p), Some(Normal::Line(line[0], line[1])))
}

/// What the construction assumed of its loop's points, to be proven over
/// a drafted profile's span (walls.rs). Indices are the loop's points; a
/// line a -> b has its piece on its left.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Ob {
    /// p is left of the line a -> b by more than eps, or within eps of it.
    Side { a: usize, b: usize, p: usize, on: bool },
    /// q - p runs forward along a -> b.
    Along { p: usize, q: usize, a: usize, b: usize },
    /// Edge e's arc reaches at most 100 eps past the line a -> b.
    Reach { e: usize, a: usize, b: usize },
    /// From normal n to normal m turns counter-clockwise (convex), or
    /// clockwise by at most `turn_out`'s angle.
    Turn { n: Normal, m: Normal },
    /// A piece found empty (a sliver).
    Sliver,
}

/// An outward normal, unnormalised: right of the line a -> b, from edge
/// e's arc's centre to point p, or from point p to edge e's arc's centre.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Normal {
    Line(usize, usize),
    Arc(usize, usize),
    Toward(usize, usize),
}

/// What the construction (`build`) assumed, and how many pieces it read.
#[derive(Default)]
pub(crate) struct Trace {
    pub(crate) obs: Vec<Ob>,
    pieces: usize,
}
/// An edge of the loop: from its point to the next, a line where b = 0,
/// else an arc of bulge b = tan(sweep / 4), positive turning
/// counter-clockwise (lying right of the chord).
#[derive(Clone, Copy, Debug)]
struct Arc {
    c: P,
    r: f64,
}

/// The arc from a to q of bulge b: its centre and radius.
fn arc_of(a: P, q: P, b: f64) -> Arc {
    let e = sub(q, a);
    let l = norm(e);
    let s = b * l / 2.0; // signed sagitta, positive to the right
    let r = l * (1.0 + b * b) / (4.0 * b.abs());
    let rn = [e[1] / l, -e[0] / l];
    let apex = add(mul(add(a, q), 0.5), mul(rn, s));
    Arc { c: sub(apex, mul(rn, s.signum() * r)), r }
}

/// The arc's midpoint and the bulge of each half.
fn halve(a: P, q: P, b: f64) -> (P, f64) {
    let e = sub(q, a);
    let m = add(mul(add(a, q), 0.5), mul([e[1], -e[0]], b / 2.0));
    (m, b / (1.0 + (1.0 + b * b).sqrt()))
}

const QUARTER: f64 = 0.414_213_562_373_095_1; // tan(pi / 8): a quarter turn's bulge

/// The turn on a circle of radius r whose sag is d: an arc turned out past a
/// line by less pokes past it by under d (polygon's allowance for its data,
/// 100 eps, as a sliver that thin).
fn sag_turn(r: f64, d: f64) -> f64 {
    if d >= 2.0 * r { std::f64::consts::PI } else { (1.0 - d / r).acos() }
}

/// The sag of a circle of radius r over a chord d (r - sqrt(r^2 - d^2 / 4),
/// as (d / 2)^2 / (r + sqrt(r^2 - d^2 / 4)), which does not cancel), in an
/// outward bound: at least the true sag.
fn sag_over(r: f64, d: f64) -> f64 {
    use ink_core::interval::{self, I};
    let b = |op: BinOp, x: I, y: I| interval::bin(op, x, y).0;
    let q = b(BinOp::Mul, I::point(0.5 * d), I::point(0.5 * d));
    let rest = b(BinOp::Sub, b(BinOp::Mul, I::point(r), I::point(r)), q);
    let root = if rest.lo > 0.0 { interval::un(UnOp::Sqrt, rest).lo } else { 0.0 };
    b(BinOp::Div, q, b(BinOp::Add, I::point(r), I::point(root))).hi.min(r)
}

/// A lower bound on the sine of how far an arc bulging out of a piece may
/// turn out past its neighbour where they meet, read as it is (no pocket
/// bridging the corner): until the corner lies d (polygon's allowance)
/// inside the line that bridges it, tangent to both sides (1 where the turn
/// reaches a quarter). Each side that curves away from that line, an arc
/// bulging out of radius r (else 0: a line, or an arc bulging in), turns by
/// sag_turn(r, d) to its contact, whose cosine is 1 - d / r and sine
/// sqrt((d / r) (2 - d / r)): the sum's sine and cosine, in outward bounds
/// (no rounded acos or sin decides what is within the allowance).
fn turn_out(r: [f64; 2], d: f64) -> f64 {
    use ink_core::interval::{self, I};
    let b = |op: BinOp, x: I, y: I| interval::bin(op, x, y).0;
    let mut sum: Option<(I, I)> = None;
    for &r in r.iter().filter(|&&r| r > 0.0) {
        if d >= 2.0 * r {
            return 1.0;
        }
        let q = b(BinOp::Div, I::point(d), I::point(r));
        let (c, s) = (b(BinOp::Sub, I::point(1.0), q), interval::un(UnOp::Sqrt, b(BinOp::Mul, q, b(BinOp::Sub, I::point(2.0), q))));
        sum = Some(match sum {
            None => (c, s),
            Some((c0, s0)) => (b(BinOp::Sub, b(BinOp::Mul, c0, c), b(BinOp::Mul, s0, s)), b(BinOp::Add, b(BinOp::Mul, s0, c), b(BinOp::Mul, c0, s))),
        });
    }
    match sum {
        None => 0.0,
        Some((c, _)) if c.hi <= 0.0 => 1.0,
        Some((_, s)) => s.lo.max(0.0),
    }
}

pub(crate) fn polygon(input: &[[f64; 3]]) -> Result<T, String> {
    let (f, lp, _) = build(input, allowance(input))?;
    Ok(T::Max(vec![f, lp.bounds()]))
}

/// polygon's eps: 1e-9 of the size of its points' box (its allowance for its
/// data is 100 eps).
fn allowance(input: &[[f64; 3]]) -> f64 {
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in input {
        for k in 0..2 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    1e-9 * norm(sub(hi, lo)).max(1e-300)
}

/// The construction, reading points within eps as one (`allowance`): the
/// profile's tree (without its bounding box), its loop as split, and what it
/// assumed.
fn build(input: &[[f64; 3]], eps: f64) -> Result<(T, Loop, Trace), String> {
    // Points, and each edge's bulge; repeated points dropped (the edge on
    // from a repeated point, and so the point, is its last copy's).
    let mut pts: Vec<P> = Vec::new();
    let mut bs: Vec<f64> = Vec::new();
    let mut defs: Vec<Def> = Vec::new();
    for (k, p) in input.iter().enumerate() {
        if let Some(&last) = pts.last()
            && norm(sub([p[0], p[1]], last)) <= eps
        {
            *bs.last_mut().expect("a point before") = p[2];
            *defs.last_mut().expect("a point before") = Def::Given(k);
            continue;
        }
        pts.push([p[0], p[1]]);
        bs.push(p[2]);
        defs.push(Def::Given(k));
    }
    while pts.len() > 1 && norm(sub(pts[0], pts[pts.len() - 1])) <= eps {
        pts.pop();
        bs.pop();
        defs.pop();
    }
    // An arc too flat to tell from its chord (or no arc) is a line.
    for i in 0..pts.len() {
        let l = norm(sub(pts[(i + 1) % pts.len()], pts[i]));
        if (bs[i] * l / 2.0).abs() <= eps {
            bs[i] = 0.0;
        }
    }
    if pts.len() < 2 || (pts.len() < 3 && bs.iter().all(|&b| b == 0.0)) {
        return Err(format!("a polygon needs 3 corners (or arcs to close it); it has {}", pts.len()));
    }
    let ids = (0..pts.len() as u32).collect();
    let next = pts.len() as u32;
    let mut lp = Loop { pts, bs, eps, defs, ids, next };
    lp.split(|_, b| b.abs() > QUARTER * (1.0 + 1e-12));
    lp.merge_straight();
    if let Some((i, j, p)) = lp.meets() {
        let (a, b) = (lp.pts[i], lp.pts[j]);
        return Err(format!("the polygon crosses itself at ({}, {}): its edges from ({}, {}) and from ({}, {}) meet there", p[0], p[1], a[0], a[1], b[0], b[1]));
    }
    // The profile is its chord polygon with each arc's segment (between it
    // and its chord) added or taken away, so long as no other edge reaches
    // into a segment or across a chord. Each round halves the arcs that
    // break that (a simple loop's arcs, halved, close in on themselves until
    // nothing else is near enough; an arc within eps of its chord is that
    // chord, so the loop crosses itself), or splits arcs where the hulls of
    // one more level of pieces touch them (a piece is read once its hull is
    // its chord polygon's; each cut is a new point, so it ends). The loop
    // turns as its chord polygon does once that is simple.
    loop {
        if let Some((i, j)) = lp.crossing() {
            let sag = |e: usize| (lp.bs[e] * norm(sub(lp.pts[(e + 1) % lp.n()], lp.pts[e])) / 2.0).abs();
            let halve = [i, j].map(|e| (sag(e) > eps).then_some(e));
            if halve.iter().any(Option::is_some) {
                lp.split(|e, _| halve.contains(&Some(e)));
                continue;
            }
            let (a, b) = (lp.pts[i], lp.pts[j]);
            return Err(format!("the polygon crosses itself: its edges from ({}, {}) and from ({}, {}) meet", a[0], a[1], b[0], b[1]));
        }
        if lp.area() < 0.0 {
            lp.reverse();
        }
        let mut cuts = Vec::new();
        let mut tr = Trace::default();
        let root: Vec<usize> = (0..lp.pts.len()).collect();
        let f = lp.piece(&root, false, &[], &mut cuts, &mut tr)?;
        if cuts.is_empty() {
            return Ok((f, lp, tr));
        }
        if !lp.split_at(&cuts) {
            return Err("a polygon's arcs reach past the hull of their corners where the hull touches them (is the profile self-touching?)".into());
        }
    }
}

pub(crate) struct Loop {
    pts: Vec<P>,
    bs: Vec<f64>,
    eps: f64,
    /// Where each point comes from, its id (kept as points are added),
    /// and the next id.
    pub(crate) defs: Vec<Def>,
    ids: Vec<u32>,
    next: u32,
}

/// Where a loop point comes from (a drafted profile's points move as their
/// definitions do).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Def {
    /// Point k of the input.
    Given(usize),
    /// Where an arc is halved.
    Half,
    /// Where a hull's line from it to the point with this id touches its
    /// arc (from a corner, or common to two arcs: each point the other's).
    Touch(u32),
}

/// Where a hull touches an arc: its loop edge, the arc's outward normal
/// there, and where the hull's line through it runs: to a loop point, or to
/// the cut on another arc it is common to.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Cut {
    e: usize,
    n: P,
    to: To,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum To {
    Pt(usize),
    Cut(usize),
}

/// A piece of a piece's boundary, from a, with its outward normal at each
/// end.
struct Seg {
    a: P,
    n0: P,
    n1: P,
    kind: Kind,
    /// The hull edge it lies under (its chord's line, for an arc), and
    /// that edge's place in the hull.
    line: [usize; 2],
    hull: usize,
    /// The loop edge it is, if it is one.
    edge: Option<usize>,
    /// Where n0 and n1 come from.
    src: [Normal; 2],
}

#[derive(Clone, Copy)]
enum Kind {
    Line,
    /// A pocket's own lid: no term.
    Lid,
    /// An arc bulging out of the piece, and its chord's normal.
    Out(Arc, P),
    /// An arc bulging into the piece, and its chord's normal.
    Cut(Arc, P),
}

impl Loop {
    fn n(&self) -> usize {
        self.pts.len()
    }

    /// Halve every arc `which` picks, repeatedly while it still picks it.
    fn split(&mut self, which: impl Fn(usize, f64) -> bool) {
        let (mut pts, mut bs, mut defs, mut ids) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for i in 0..self.n() {
            let (a, q) = (self.pts[i], self.pts[(i + 1) % self.n()]);
            let mut stack = vec![(a, q, self.bs[i], which(i, self.bs[i]))];
            let mut out = Vec::new();
            while let Some((a, q, b, go)) = stack.pop() {
                if b != 0.0 && go {
                    let (m, h) = halve(a, q, b);
                    let more = h.abs() > QUARTER * (1.0 + 1e-12);
                    stack.push((m, q, h, more));
                    stack.push((a, m, h, more));
                } else {
                    out.push((a, b));
                }
            }
            for (k, (a, b)) in out.into_iter().enumerate() {
                pts.push(a);
                bs.push(b);
                if k == 0 {
                    defs.push(self.defs[i]);
                    ids.push(self.ids[i]);
                } else {
                    defs.push(Def::Half);
                    ids.push(self.next);
                    self.next += 1;
                }
            }
        }
        (self.pts, self.bs, self.defs, self.ids) = (pts, bs, defs, ids);
    }

    /// Split each arc `cuts` names where its outward normal (from its
    /// centre) is the one given: where a hull touches it. A cut within 10
    /// eps of sag of its arc's end, or of another cut, falls on that point.
    /// Whether any was made.
    fn split_at(&mut self, cuts: &[Cut]) -> bool {
        let n = self.n();
        let ang = |u: P, v: P| cross(u, v).atan2(dot(u, v));
        // Each cut's turn from its edge's start, and the id it gets.
        let mut place: Vec<Option<f64>> = vec![None; cuts.len()];
        let mut id: Vec<u32> = vec![0; cuts.len()];
        for i in 0..n {
            let (a, q, b) = (self.pts[i], self.pts[(i + 1) % n], self.bs[i]);
            let mut here: Vec<(f64, usize)> = cuts.iter().enumerate().filter(|c| c.1.e == i && b != 0.0).map(|(k, c)| (ang(sub(a, arc_of(a, q, b).c), c.n), k)).collect();
            if here.is_empty() {
                continue;
            }
            let (arc, sweep) = (arc_of(a, q, b), 4.0 * b.atan());
            let gap = sag_turn(arc.r, 10.0 * self.eps);
            here.sort_by(|x, y| (x.0 * sweep.signum()).total_cmp(&(y.0 * sweep.signum())));
            let mut last: Option<(f64, u32)> = None;
            for (t, k) in here {
                id[k] = if t * sweep.signum() <= gap {
                    self.ids[i]
                } else if (sweep - t) * sweep.signum() <= gap {
                    self.ids[(i + 1) % n]
                } else if let Some((lt, li)) = last
                    && (t - lt).abs() <= gap
                {
                    li
                } else {
                    place[k] = Some(t);
                    self.next += 1;
                    last = Some((t, self.next - 1));
                    self.next - 1
                };
            }
        }
        if place.iter().all(Option::is_none) {
            return false;
        }
        let (mut pts, mut bs, mut defs, mut ids) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            pts.push(self.pts[i]);
            defs.push(self.defs[i]);
            ids.push(self.ids[i]);
            let (a, q, b) = (self.pts[i], self.pts[(i + 1) % n], self.bs[i]);
            let mut at: Vec<(f64, usize)> = (0..cuts.len()).filter(|&k| cuts[k].e == i).filter_map(|k| place[k].map(|t| (t, k))).collect();
            if at.is_empty() {
                bs.push(b);
                continue;
            }
            let (arc, sweep) = (arc_of(a, q, b), 4.0 * b.atan());
            at.sort_by(|x, y| (x.0 * sweep.signum()).total_cmp(&(y.0 * sweep.signum())));
            let u = unit(sub(a, arc.c));
            let mut from = 0.0;
            for (t, k) in at {
                bs.push(((t - from) / 4.0).tan());
                pts.push(add(arc.c, mul([u[0] * t.cos() - u[1] * t.sin(), u[0] * t.sin() + u[1] * t.cos()], arc.r)));
                defs.push(Def::Touch(match cuts[k].to {
                    To::Pt(j) => self.ids[j],
                    To::Cut(c) => id[c],
                }));
                ids.push(id[k]);
                from = t;
            }
            bs.push(((sweep - from) / 4.0).tan());
        }
        (self.pts, self.bs, self.defs, self.ids) = (pts, bs, defs, ids);
        true
    }

    /// Corners on a straight run of lines add nothing.
    fn merge_straight(&mut self) {
        let mut k = 0;
        while k < self.n() && self.n() > 3 {
            let n = self.n();
            let (ip, inx) = ((k + n - 1) % n, (k + 1) % n);
            let (a, b, c) = (self.pts[ip], self.pts[k], self.pts[inx]);
            let d = sub(c, a);
            if self.bs[ip] == 0.0
                && self.bs[k] == 0.0
                && norm(d) > self.eps
                && cross(d, sub(b, a)).abs() / norm(d) <= self.eps
                && dot(sub(b, a), d) > 0.0
                && dot(sub(c, b), d) > 0.0
            {
                self.pts.remove(k);
                self.bs.remove(k);
                self.defs.remove(k);
                self.ids.remove(k);
                k = k.saturating_sub(1);
            } else {
                k += 1;
            }
        }
    }

    fn area(&self) -> f64 {
        (0..self.n()).map(|i| cross(self.pts[i], self.pts[(i + 1) % self.n()])).sum()
    }

    fn reverse(&mut self) {
        let n = self.n();
        let pts: Vec<P> = (0..n).map(|k| self.pts[n - 1 - k]).collect();
        // Edge k now runs from old point n-1-k to n-2-k: old edge n-2-k, reversed.
        let bs: Vec<f64> = (0..n).map(|k| -self.bs[(2 * n - 2 - k) % n]).collect();
        self.pts = pts;
        self.bs = bs;
        self.defs.reverse();
        self.ids.reverse();
    }

    /// Two edges as drawn (lines and arcs, each at most a quarter turn) that
    /// meet, or come within eps, other than beside a corner they share,
    /// within the allowance for the data (100 eps): the edges and a point where
    /// they meet. Each pair's candidates are where their lines and circles
    /// cross (beside a shared corner, the other crossing, which stays exact
    /// where they are tangent there), each one's ends and midpoint (a
    /// touching corner, an overlap).
    fn meets(&self) -> Option<(usize, usize, P)> {
        let (n, eps) = (self.n(), self.eps);
        let drawn = |e: usize| {
            let (a, q, b) = (self.pts[e], self.pts[(e + 1) % n], self.bs[e]);
            (a, q, (b != 0.0).then(|| arc_of(a, q, b)), if b == 0.0 { mul(add(a, q), 0.5) } else { halve(a, q, b).0 })
        };
        // Within eps of edge e as drawn (by unit directions: a profile at
        // any scale, its products unable to underflow).
        let on = |e: usize, p: P| {
            let (a, q, arc, _) = drawn(e);
            match arc {
                None => {
                    let (d, l) = (unit(sub(q, a)), norm(sub(q, a)));
                    let t = dot(sub(p, a), d);
                    cross(d, sub(p, a)).abs() <= eps && t >= -eps && t <= l + eps
                }
                Some(o) => {
                    let (ua, uq, u) = (unit(sub(a, o.c)), unit(sub(q, o.c)), sub(p, o.c));
                    let (sg, w) = (cross(ua, uq).signum(), unit(u));
                    (norm(u) - o.r).abs() <= eps && sg * cross(ua, w) >= -eps / o.r && sg * cross(w, uq) >= -eps / o.r && dot(w, add(ua, uq)) > 0.0
                }
            }
        };
        let reach = |e: usize| {
            let (a, q, _, _) = drawn(e);
            let pad = (self.bs[e] * norm(sub(q, a)) / 2.0).abs() + eps;
            ([a[0].min(q[0]) - pad, a[1].min(q[1]) - pad], [a[0].max(q[0]) + pad, a[1].max(q[1]) + pad])
        };
        for i in 0..n {
            let (lo, hi) = reach(i);
            for j in i + 1..n {
                let (lo2, hi2) = reach(j);
                if lo[0] > hi2[0] || lo2[0] > hi[0] || lo[1] > hi2[1] || lo2[1] > hi[1] {
                    continue;
                }
                let shared: Vec<P> = [(j == i + 1).then(|| self.pts[j]), (i == (j + 1) % n).then(|| self.pts[i])].into_iter().flatten().collect();
                let (ei, ej) = (drawn(i), drawn(j));
                let mut cands = vec![ei.0, ei.1, ei.3, ej.0, ej.1, ej.3];
                match (ei.2, ej.2) {
                    (None, None) => {
                        let (d1, d2) = (sub(ei.1, ei.0), sub(ej.1, ej.0));
                        let den = cross(d1, d2);
                        if shared.is_empty() && den != 0.0 {
                            cands.push(add(ei.0, mul(d1, cross(sub(ej.0, ei.0), d2) / den)));
                        }
                    }
                    (None, Some(o)) | (Some(o), None) => {
                        let (a, q) = if ei.2.is_none() { (ei.0, ei.1) } else { (ej.0, ej.1) };
                        let u = unit(sub(q, a));
                        if let Some(&s) = shared.first() {
                            cands.push(sub(s, mul(u, 2.0 * dot(sub(s, o.c), u))));
                        } else {
                            let foot = add(a, mul(u, dot(sub(o.c, a), u)));
                            let h = norm(sub(foot, o.c));
                            if h <= o.r + eps {
                                let w = (o.r * o.r - h * h).max(0.0).sqrt();
                                cands.extend([add(foot, mul(u, w)), sub(foot, mul(u, w))]);
                            }
                        }
                    }
                    (Some(o1), Some(o2)) => {
                        let dc = sub(o2.c, o1.c);
                        let d = norm(dc);
                        if d <= eps {
                            // One circle: an overlap shows at an end or a midpoint.
                        } else if let Some(&s) = shared.first() {
                            let u = mul(dc, 1.0 / d);
                            let v = sub(s, o1.c);
                            cands.push(add(o1.c, sub(mul(u, 2.0 * dot(v, u)), v)));
                        } else if d <= o1.r + o2.r + eps && d >= (o1.r - o2.r).abs() - eps {
                            let u = mul(dc, 1.0 / d);
                            let a = (d * d + o1.r * o1.r - o2.r * o2.r) / (2.0 * d);
                            let h = (o1.r * o1.r - a * a).max(0.0).sqrt();
                            let m = add(o1.c, mul(u, a));
                            cands.extend([add(m, mul([-u[1], u[0]], h)), sub(m, mul([-u[1], u[0]], h))]);
                        }
                    }
                }
                // Beside a shared corner, a meeting whose sliver between the
                // two is within the allowance for the data (100 eps): a
                // tangent joint whose points were rounded. The sliver is at
                // most the sum of each arc's sag over the chord from the
                // corner (sag_over, rounded up).
                let arcs: Vec<f64> = [ei.2, ej.2].iter().flatten().map(|o| o.r).collect();
                let apart = |p: P, s: P| {
                    let d = norm(sub(p, s));
                    d > 100.0 * eps && (arcs.is_empty() || arcs.iter().map(|&r| sag_over(r, d)).sum::<f64>() > 100.0 * eps)
                };
                if let Some(p) = cands.into_iter().find(|&p| shared.iter().all(|&s| apart(p, s)) && on(i, p) && on(j, p)) {
                    return Some((i, j, p));
                }
            }
        }
        None
    }

    /// Two chord edges that meet other than at a shared corner, or an arc
    /// and an edge reaching into its segment by more than eps: a corner
    /// inside it (the arc twice), or another arc dipping across its chord.
    fn crossing(&self) -> Option<(usize, usize)> {
        let n = self.n();
        if n < 3 {
            return Some((0, 1 % n));
        }
        let eps = self.eps;
        let side = |a: P, b: P, c: P| cross(sub(b, a), sub(c, a));
        for i in 0..n {
            let (a, b) = (self.pts[i], self.pts[(i + 1) % n]);
            for j in i + 1..n {
                let adjacent = j == i + 1 || (i == 0 && j == n - 1);
                let (c, d) = (self.pts[j], self.pts[(j + 1) % n]);
                if a[0].min(b[0]) > c[0].max(d[0]) + eps || c[0].min(d[0]) > a[0].max(b[0]) + eps || a[1].min(b[1]) > c[1].max(d[1]) + eps || c[1].min(d[1]) > a[1].max(b[1]) + eps
                {
                    continue;
                }
                if adjacent {
                    // Adjacent edges meet only at their corner: not folded back on each other.
                    let (p, q, r) = if j == i + 1 { (a, b, d) } else { (c, a, b) };
                    let (u, v) = (sub(p, q), sub(r, q));
                    if cross(u, v).abs() <= eps * norm(u).max(norm(v)) && dot(u, v) > 0.0 {
                        return Some((i, j));
                    }
                    continue;
                }
                let (d1, d2, d3, d4) = (side(a, b, c), side(a, b, d), side(c, d, a), side(c, d, b));
                if d1 * d2 <= 0.0 && d3 * d4 <= 0.0 {
                    return Some((i, j));
                }
            }
        }
        for i in (0..n).filter(|&i| self.bs[i] != 0.0) {
            let (a, q) = (self.pts[i], self.pts[(i + 1) % n]);
            let (arc, out) = (arc_of(a, q, self.bs[i]), mul(right(a, q), self.bs[i].signum()));
            let depth = |p: P| dot(sub(p, a), out).min(arc.r - norm(sub(p, arc.c)));
            for j in 0..n {
                if j != i && j != (i + 1) % n && depth(self.pts[j]) > eps {
                    return Some((i, i));
                }
                if j == i || self.bs[j] == 0.0 {
                    continue;
                }
                // Arc j's circle across the chord's line, where both
                // crossings are on the chord and on arc j: the arc between.
                let o = arc_of(self.pts[j], self.pts[(j + 1) % n], self.bs[j]);
                let (d, f) = (sub(q, a), sub(a, o.c));
                let (aa, bb, cc) = (dot(d, d), 2.0 * dot(f, d), dot(f, f) - o.r * o.r);
                let disc = bb * bb - 4.0 * aa * cc;
                if disc <= 0.0 {
                    continue;
                }
                let (ua, uq) = (sub(self.pts[j], o.c), sub(self.pts[(j + 1) % n], o.c));
                let sg = cross(ua, uq).signum();
                let on = |t: f64| {
                    let u = sub(add(a, mul(d, t)), o.c);
                    t * norm(d) > eps && (1.0 - t) * norm(d) > eps && sg * cross(ua, u) > 0.0 && sg * cross(u, uq) > 0.0
                };
                let (t1, t2) = ((-bb - disc.sqrt()) / (2.0 * aa), (-bb + disc.sqrt()) / (2.0 * aa));
                if on(t1) && on(t2) {
                    let m = unit(add(sub(add(a, mul(d, t1)), o.c), sub(add(a, mul(d, t2)), o.c)));
                    if depth(add(o.c, mul(m, o.r))) > eps {
                        return Some((i, j));
                    }
                }
            }
        }
        None
    }

    /// The most of n . x over loop edge e (its arc, or its chord).
    fn reach(&self, e: usize, n: P) -> f64 {
        let (a, q, b) = (self.pts[e], self.pts[(e + 1) % self.n()], self.bs[e]);
        let mut best = dot(n, a).max(dot(n, q));
        if b != 0.0 {
            let arc = arc_of(a, q, b);
            let (ua, uq) = (sub(a, arc.c), sub(q, arc.c));
            let sg = cross(ua, uq).signum();
            if sg * cross(ua, n) >= 0.0 && sg * cross(n, uq) >= 0.0 {
                best = best.max(dot(n, arc.c) + arc.r);
            }
        }
        best
    }

    /// Edge index between consecutive loop points i, j (either way round).
    fn edge(&self, i: usize, j: usize) -> (usize, bool) {
        if (i + 1) % self.n() == j { (i, true) } else { (j, false) }
    }

    /// The convex piece over vs (loop point indices, counter-clockwise for
    /// the piece; consecutive ones joined by loop edges, the last to the
    /// first by a loop edge too, or, in a pocket, by its lid, left out),
    /// within the hull of the piece above (`within`: its lines' outward
    /// normals, a point and the loop points they pass through). Where an
    /// arc bulging out of the piece reaches past the hull of its chord
    /// polygon, that is not the piece's hull: the arcs are cut where the
    /// hull touches them (`cuts`) and the piece is read again next round.
    fn piece(&self, vs: &[usize], pocket: bool, within: &[(P, P, [usize; 2])], cuts: &mut Vec<Cut>, tr: &mut Trace) -> Result<T, String> {
        let p = |k: usize| self.pts[vs[k]];
        let m = vs.len();
        tr.pieces += 1;
        // Its loop edges, as places in vs, and each one's bulge as the
        // piece runs (positive: bulging out of it).
        let mut chain: Vec<(usize, usize)> = (1..m).map(|k| (k - 1, k)).collect();
        if !pocket {
            chain.push((m - 1, 0));
        }
        let bulge = |a: usize, b: usize| {
            let (e, fwd) = self.edge(vs[a], vs[b]);
            (e, if fwd { self.bs[e] } else { -self.bs[e] })
        };
        for &(a, b) in &chain {
            let (e, bl) = bulge(a, b);
            if bl != 0.0 {
                tr.obs.extend(within.iter().map(|&(_, _, [a, b])| Ob::Reach { e, a, b }));
            }
        }
        let hull = {
            let pts: Vec<P> = (0..m).map(p).collect();
            let mut h = hull(&pts, self.eps);
            h.sort_unstable();
            h
        };
        if hull.len() < 3 {
            // A sliver within 1e-9 of the profile's size of a line: empty.
            tr.obs.push(Ob::Sliver);
            return Ok(T::Lin(0.0, 0.0, 1.0, None));
        }
        // Each point stays on its side of each hull line (or on it), and
        // each hull corner stays a corner.
        let h = hull.len();
        for k in 0..h {
            let (a, b, c) = (vs[hull[k]], vs[hull[(k + 1) % h]], vs[hull[(k + h - 1) % h]]);
            let side = |q: usize| cross(sub(self.pts[b], self.pts[a]), sub(self.pts[q], self.pts[a])).abs() <= self.eps * norm(sub(self.pts[b], self.pts[a]));
            tr.obs.extend(vs.iter().filter(|&&q| q != a && q != b).map(|&q| Ob::Side { a, b, p: q, on: side(q) }));
            if h > 3 {
                tr.obs.push(Ob::Side { a: b, b: c, p: a, on: false });
            }
        }
        let lines: Vec<(P, P, [usize; 2])> = (0..h).map(|k| (right(p(hull[k]), p(hull[(k + 1) % h])), p(hull[k]), [vs[hull[k]], vs[hull[(k + 1) % h]]])).collect();
        // Its boundary in order: loop edges, lids over pockets, and a
        // pocket's own lid (no term).
        let mut segs: Vec<Seg> = Vec::new();
        let mut pockets: Vec<Vec<usize>> = Vec::new();
        for k in 0..h {
            let (s, f) = (hull[k], hull[(k + 1) % h]);
            let wrap = k + 1 == h;
            let run: Vec<usize> = if wrap { (s..m).chain(0..=f).collect() } else { (s..=f).collect() };
            let (a, b) = (p(s), p(f));
            let n = right(a, b);
            let line = [vs[s], vs[f]];
            let flat = [Normal::Line(vs[s], vs[f]); 2];
            if pocket && wrap && s == m - 1 && f == 0 {
                segs.push(Seg { a, n0: n, n1: n, kind: Kind::Lid, line, hull: k, edge: None, src: flat });
                continue;
            }
            // A convex arc owns its chord's half-plane: another edge on
            // that line would clip its bump back to the chord. Such a run
            // goes through the true hull's contacts and pockets below.
            // Inward arcs only subtract discs, so may share a hull line.
            let on_line = (run.len() == 2 || run.windows(2).all(|w| bulge(w[0], w[1]).1 <= 0.0))
                && run[1..run.len() - 1].iter().all(|&i| dot(sub(p(i), a), n).abs() <= self.eps)
                && run.windows(2).all(|w| dot(sub(p(w[1]), p(w[0])), sub(b, a)) > 0.0);
            if !on_line {
                // A pocket under a lid: the lid here, the pocket read the
                // same way (counter-clockwise: reversed).
                pockets.push(run.iter().rev().map(|&i| vs[i]).collect());
                segs.push(Seg { a, n0: n, n1: n, kind: Kind::Line, line, hull: k, edge: None, src: flat });
                continue;
            }
            for w in run.windows(2) {
                let (i, j) = (vs[w[0]], vs[w[1]]);
                tr.obs.push(Ob::Along { p: i, q: j, a: vs[s], b: vs[f] });
                let (pi, pj) = (self.pts[i], self.pts[j]);
                let (e, bl) = bulge(w[0], w[1]);
                let seg = if bl == 0.0 {
                    Seg { a: pi, n0: n, n1: n, kind: Kind::Line, line, hull: k, edge: Some(e), src: flat }
                } else {
                    let arc = arc_of(pi, pj, bl);
                    if bl > 0.0 {
                        let (n0, n1) = (unit(sub(pi, arc.c)), unit(sub(pj, arc.c)));
                        Seg { a: pi, n0, n1, kind: Kind::Out(arc, n), line, hull: k, edge: Some(e), src: [Normal::Arc(e, i), Normal::Arc(e, j)] }
                    } else {
                        // Bulging in, its outward normals point into its disc.
                        let (n0, n1) = (unit(sub(arc.c, pi)), unit(sub(arc.c, pj)));
                        Seg { a: pi, n0, n1, kind: Kind::Cut(arc, n), line, hull: k, edge: Some(e), src: [Normal::Toward(e, i), Normal::Toward(e, j)] }
                    }
                };
                segs.push(seg);
            }
        }
        let ns = segs.len();
        // Where an arc bulging out meets one bulging in and turns out past
        // its chord, tangent or turning in (an inflection: no hull follows
        // the two but through a pocket at every scale), the neighbour's
        // chord gives way on the arc's side of the arc's tangent there.
        // So does a hull edge over a pocket whose first edge is a line
        // tangent to the arc there (a lobe, a straight run, a concave
        // fillet): the line's half-plane, which the pocket lies beyond.
        // Each arc's end that meets its neighbour so: the normal it meets.
        let mut relax: Vec<Vec<(P, P, Normal)>> = vec![Vec::new(); ns];
        let mut free: Vec<[Option<(P, Normal)>; 2]> = vec![[None; 2]; ns];
        let slack = 100.0 * self.eps;
        let place = |v: usize| vs.iter().position(|&x| x == v).expect("a seg's point is the piece's");
        for k in 0..ns {
            let Kind::Out(arc, _) = segs[k].kind else { continue };
            let (ip, inx) = ((k + ns - 1) % ns, (k + 1) % ns);
            // Tangent or turning in but for a sliver, and turned out past
            // the chord by more.
            let meets = |u: P, v: P, r: f64| cross(u, v) >= -turn_out([r, 0.0], slack) && dot(u, v) > 0.0;
            let past = -sag_turn(arc.r, slack).sin();
            match segs[inx].kind {
                Kind::Cut(b, c) if meets(segs[k].n1, segs[inx].n0, arc.r.max(b.r)) && cross(segs[k].n1, c) < past => {
                    relax[inx].push((segs[k].n1, segs[inx].a, segs[k].src[1]));
                    free[k][1] = Some((segs[inx].n0, segs[inx].src[0]));
                }
                Kind::Line if segs[inx].edge.is_none() && cross(segs[k].n1, segs[inx].n0) < past => {
                    // The pocket's first edge, from this arc's end.
                    let at = place(segs[inx].line[0]);
                    let (_, bl) = bulge(at, (at + 1) % m);
                    let n = right(p(at), p((at + 1) % m));
                    if bl == 0.0 && meets(segs[k].n1, n, arc.r) {
                        let src = Normal::Line(vs[at], vs[(at + 1) % m]);
                        relax[inx].push((n, p(at), src));
                        free[k][1] = Some((n, src));
                    }
                }
                _ => {}
            }
            match segs[ip].kind {
                Kind::Cut(b, c) if meets(segs[ip].n1, segs[k].n0, arc.r.max(b.r)) && cross(c, segs[k].n0) < past => {
                    relax[ip].push((segs[k].n0, segs[k].a, segs[k].src[0]));
                    free[k][0] = Some((segs[ip].n1, segs[ip].src[1]));
                }
                Kind::Line if segs[ip].edge.is_none() && cross(segs[ip].n0, segs[k].n0) < past => {
                    // The pocket's last edge, to this arc's start.
                    let at = place(segs[ip].line[1]);
                    let (_, bl) = bulge((at + m - 1) % m, at);
                    let n = right(p((at + m - 1) % m), p(at));
                    if bl == 0.0 && meets(n, segs[k].n0, arc.r) {
                        let src = Normal::Line(vs[(at + m - 1) % m], vs[at]);
                        relax[ip].push((n, p(at), src));
                        free[k][0] = Some((n, src));
                    }
                }
                _ => {}
            }
        }
        // An arc bulging out past a hull line (but its own chord's, and a
        // neighbour's that gives way to it) by more than polygon's allowance
        // says this is not the hull: such arcs are cut where the hull of the
        // piece's corners and arcs touches them.
        let pokes: Vec<usize> = chain
            .iter()
            .filter_map(|&(a, b)| {
                let (e, bl) = bulge(a, b);
                let own = segs.iter().position(|s| s.edge == Some(e));
                (bl > 0.0
                    && lines.iter().enumerate().any(|(k, &(n, q, _))| {
                        let skip = own.is_some_and(|s| {
                            segs[s].hull == k || (free[s][1].is_some() && segs[(s + 1) % ns].hull == k) || (free[s][0].is_some() && segs[(s + ns - 1) % ns].hull == k)
                        });
                        !skip && self.reach(e, n) - dot(n, q) > 100.0 * self.eps
                    }))
                .then_some(e)
            })
            .collect();
        // An arc whose end turns out past its neighbour's normal there (its
        // chord's, or its own where it gives way, or an arc's) by more than
        // `turn_out` (as `Ob::Turn` reads it) is cut where the hull bridging
        // the corner touches it; where that is within its ends' gap the
        // corner's pocket is too (the piece is right to the allowance: a
        // reflex join of arcs each inside the other's chord is read right by
        // their terms).
        let mut turns: Vec<usize> = Vec::new();
        for k in 0..ns {
            let Kind::Out(arc, _) = segs[k].kind else { continue };
            for end in [0, 1] {
                let o = if end == 1 { (k + 1) % ns } else { (k + ns - 1) % ns };
                let nb = &segs[o];
                let (theirs, r) = match (free[k][end], nb.kind) {
                    (Some((n, _)), _) => (n, 0.0),
                    (None, Kind::Cut(_, c)) => (c, 0.0),
                    (None, Kind::Out(b, _)) => (if end == 1 { nb.n0 } else { nb.n1 }, b.r),
                    (None, _) => (nb.n0, 0.0),
                };
                let mine = if end == 1 { segs[k].n1 } else { segs[k].n0 };
                let sin = if end == 1 { cross(mine, theirs) } else { cross(theirs, mine) };
                // Turned in by no more than the allowance's turn (none
                // where no arc bulges out), in floats: each normal is a unit
                // vector from a difference, a hypot and a division, each
                // within a few ulps of the emitted geometry's, so their cross
                // is within TOL of its (and a turn just under 0 is read as
                // turning in).
                const TOL: f64 = 64.0 * f64::EPSILON;
                let s = turn_out([arc.r, r], slack);
                if !(sin >= TOL || (sin >= TOL - s && dot(mine, theirs) > 0.0)) {
                    turns.extend(segs[k].edge);
                    if let Kind::Out(..) = nb.kind {
                        turns.extend(nb.edge);
                    }
                }
            }
        }
        if !pokes.is_empty() || !turns.is_empty() {
            let pts: Vec<P> = (0..m).map(p).collect();
            let arcs: Vec<(usize, Arc, usize, usize)> = chain
                .iter()
                .filter_map(|&(a, b)| {
                    let (e, bl) = bulge(a, b);
                    (bl > 0.0).then(|| (e, arc_of(p(a), p(b), bl), a, b))
                })
                .collect();
            let touch = hull_contacts(&pts, &arcs, self.eps);
            // The arcs that poke or turn, cut where the hull touches them,
            // and the arcs their hull lines are common to; where the hull
            // touches none that poke (they lie inside it, under a line that
            // misses another arc's contact), every arc it touches.
            let mut keep: Vec<bool> = touch.iter().map(|c| pokes.contains(&c.e) || turns.contains(&c.e)).collect();
            if !pokes.is_empty() && !keep.contains(&true) {
                keep.iter_mut().for_each(|k| *k = true);
            }
            for k in 0..touch.len() {
                if keep[k]
                    && let To::Cut(o) = touch[k].to
                {
                    keep[o] = true;
                }
            }
            if keep.contains(&true) {
                let base = cuts.len();
                // Each kept contact's place among those kept.
                let renum: Vec<usize> = keep
                    .iter()
                    .scan(0, |n, &k| {
                        let at = *n;
                        *n += k as usize;
                        Some(at)
                    })
                    .collect();
                cuts.extend(touch.iter().zip(&keep).filter(|c| *c.1).map(|(c, _)| Cut {
                    e: c.e,
                    n: c.n,
                    to: match c.to {
                        To::Pt(i) => To::Pt(vs[i]),
                        To::Cut(o) => To::Cut(base + renum[o]),
                    },
                }));
                return Ok(T::Lin(0.0, 0.0, 1.0, None));
            }
            if !pokes.is_empty() {
                return Err("a polygon's arcs reach past the hull of their corners, which touches none of them (is the profile self-touching?)".into());
            }
        }
        let mut children = Vec::with_capacity(pockets.len());
        for c in &pockets {
            children.push(self.piece(c, true, &lines, cuts, tr)?);
        }
        // Each arc left bulging out keeps its convex corners (where its
        // neighbour gives way, its neighbour's own tangent).
        let from = |q: &Seg, end: usize| if let Kind::Cut(..) = q.kind { Normal::Line(q.line[0], q.line[1]) } else { q.src[end] };
        for k in 0..ns {
            if let Kind::Out(..) = segs[k].kind {
                let (pv, nx) = (&segs[(k + ns - 1) % ns], &segs[(k + 1) % ns]);
                tr.obs.push(Ob::Turn { n: free[k][0].map_or(from(pv, 1), |f| f.1), m: segs[k].src[0] });
                tr.obs.push(Ob::Turn { n: segs[k].src[1], m: free[k][1].map_or(from(nx, 0), |f| f.1) });
            }
        }
        let mut terms: Vec<T> = Vec::with_capacity(2 * ns + children.len());
        for (s, r) in segs.iter().zip(&relax) {
            match s.kind {
                Kind::Lid => {}
                Kind::Line if r.is_empty() => terms.push(half(s.n0, s.a, s.line)),
                Kind::Line => {
                    let mut give = vec![half(s.n0, s.a, s.line)];
                    give.extend(r.iter().map(|&(t, q, src)| T::Lin(t[0], t[1], -dot(t, q), Some(src))));
                    terms.push(T::Min(give));
                }
                Kind::Cut(arc, n) => {
                    // Bulging in: the chord's half-plane less the disc.
                    terms.push(T::Neg(Box::new(T::Disc(arc.c[0], arc.c[1], arc.r, s.edge))));
                    let chord = half(n, s.a, s.line);
                    if r.is_empty() {
                        terms.push(chord);
                    } else {
                        let mut give = vec![chord];
                        give.extend(r.iter().map(|&(t, q, src)| T::Lin(t[0], t[1], -dot(t, q), Some(src))));
                        terms.push(T::Min(give));
                    }
                }
                Kind::Out(arc, n) => {
                    // Bulging out: the chord's half-plane or the disc.
                    terms.push(T::Min(vec![half(n, s.a, s.line), T::Disc(arc.c[0], arc.c[1], arc.r, s.edge)]));
                }
            }
        }
        terms.extend(children.into_iter().map(|c| T::Neg(Box::new(c))));
        Ok(T::Max(terms))
    }

    /// The bounding box of the profile (its arcs' extremes too), as four
    /// half-planes: interval arithmetic proves it empty far away.
    fn bounds(&self) -> T {
        let (lo, hi) = self.extent();
        T::Max(vec![T::Lin(1.0, 0.0, -hi[0], None), T::Lin(-1.0, 0.0, lo[0], None), T::Lin(0.0, 1.0, -hi[1], None), T::Lin(0.0, -1.0, lo[1], None)])
    }

    /// The profile's bounding box, its arcs' extremes too.
    fn extent(&self) -> (P, P) {
        extent(&self.pts, &self.bs, 100.0 * self.eps)
    }
}

/// A box that holds a loop of points and bulges, its arcs' extremes too,
/// grown by d (polygon's allowance) and rounded outward: each arc's extreme
/// (its circle's, as the field reads it) in outward bounds, and where its
/// sweep holds an axis within rounding, the box grown by the allowance
/// holds it anyway (an extreme just outside the sweep is under its ends').
fn extent(pts: &[P], bs: &[f64], d: f64) -> (P, P) {
    use ink_core::interval::{self, I};
    let n = pts.len();
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    let mut take = |p: [I; 2]| {
        for k in 0..2 {
            lo[k] = lo[k].min(p[k].lo);
            hi[k] = hi[k].max(p[k].hi);
        }
    };
    let iv = |op: BinOp, a: f64, b: f64| interval::bin(op, I::point(a), I::point(b)).0;
    for i in 0..n {
        let (a, q, b) = (pts[i], pts[(i + 1) % n], bs[i]);
        take([I::point(a[0]), I::point(a[1])]);
        if b != 0.0 {
            let arc = arc_of(a, q, b);
            let (ua, uq) = (sub(a, arc.c), sub(q, arc.c));
            for d in [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]] {
                // d within the arc's (under a half turn) sweep from ua to uq.
                let s = cross(ua, uq).signum();
                if s * cross(ua, d) >= 0.0 && s * cross(d, uq) >= 0.0 {
                    take([iv(BinOp::Add, arc.c[0], d[0] * arc.r), iv(BinOp::Add, arc.c[1], d[1] * arc.r)]);
                }
            }
        }
    }
    let grow = |v: f64, op: BinOp| iv(op, v, d);
    ([grow(lo[0], BinOp::Sub).lo, grow(lo[1], BinOp::Sub).lo], [grow(hi[0], BinOp::Add).hi, grow(hi[1], BinOp::Add).hi])
}

/// Indices of the strict convex hull of pts, counter-clockwise.
fn hull(pts: &[P], eps: f64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..pts.len()).collect();
    order.sort_by(|&a, &b| pts[a][0].total_cmp(&pts[b][0]).then(pts[a][1].total_cmp(&pts[b][1])));
    // a leaves the hull where o, a, b turn the wrong way, or where it lies
    // within eps of the edge o-b between its ends (on it; not where it lies
    // beyond o, as a point a hair across from o does where the points run
    // almost along the axis they are sorted by).
    let drop = |o: usize, a: usize, b: usize| {
        let (u, v) = (sub(pts[a], pts[o]), sub(pts[b], pts[o]));
        let c = cross(u, v) / norm(v).max(1e-300);
        c <= 0.0 || (c <= eps && dot(u, v) > 0.0 && dot(u, v) < dot(v, v))
    };
    let mut lower: Vec<usize> = Vec::new();
    for &i in &order {
        while lower.len() >= 2 && drop(lower[lower.len() - 2], lower[lower.len() - 1], i) {
            lower.pop();
        }
        lower.push(i);
    }
    let mut upper: Vec<usize> = Vec::new();
    for &i in order.iter().rev() {
        while upper.len() >= 2 && drop(upper[upper.len() - 2], upper[upper.len() - 1], i) {
            upper.pop();
        }
        upper.push(i);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);
    lower
}

/// Where the convex hull of pts and arcs touches each arc inside it (not
/// within 10 eps of sag of its ends): its loop edge, its outward normal
/// there, and where the hull's line through it runs. arcs are (loop edge,
/// circle, places in pts of its start and end), each bulging out, turning
/// counter-clockwise. The hull's line of normal t touches an arc where the
/// arc's support along t is the greatest of every site's. Each other site
/// beats it over one interval of t: a point, where the arc faces it, out to
/// the tangents from it; an arc, within its sweep, out to the circles'
/// common tangents (beyond its sweep its ends stand for it). So an arc is on
/// the hull over its sweep less those intervals, and touched where they
/// start and end, the line running to the site that beats it there. A site
/// that beats it by eps at most does not. Each arc is read on its own, so
/// contacts cannot be missed by the order the hull is followed in, and two
/// arcs' common tangent is found from both (the same t).
fn hull_contacts(pts: &[P], arcs: &[(usize, Arc, usize, usize)], eps: f64) -> Vec<Cut> {
    use std::f64::consts::{PI, TAU};
    let angle = |v: P| v[1].atan2(v[0]);
    // Each arc's outward normals: from t0, counter-clockwise by sw.
    let span: Vec<(f64, f64)> = arcs
        .iter()
        .map(|(_, c, a, b)| {
            let (u, v) = (sub(pts[*a], c.c), sub(pts[*b], c.c));
            (angle(u), cross(u, v).atan2(dot(u, v)).rem_euclid(TAU))
        })
        .collect();
    // The normals from `from`, counter-clockwise by len, as turns from t0
    // within 0..sw.
    let part = |from: f64, len: f64, t0: f64, sw: f64| {
        let s = (from - t0).rem_euclid(TAU);
        [s, s - TAU].map(|s| (s.max(0.0), (s + len).min(sw))).into_iter().filter(|(lo, hi)| lo < hi)
    };
    // Each contact: its arc, its turn from the arc's start, and the site.
    let mut out: Vec<(usize, f64, To)> = Vec::new();
    for (k, &(_, c, a, b)) in arcs.iter().enumerate() {
        let (t0, sw) = span[k];
        let mut beat: Vec<(f64, f64, To)> = Vec::new();
        for (j, &q) in pts.iter().enumerate() {
            let v = sub(q, c.c);
            if j != a && j != b && norm(v) - c.r > eps {
                let w = (c.r / norm(v)).acos();
                beat.extend(part(angle(v) - w, 2.0 * w, t0, sw).map(|(lo, hi)| (lo, hi, To::Pt(j))));
            }
        }
        for (k2, &(_, d, ..)) in arcs.iter().enumerate() {
            let e = sub(d.c, c.c);
            if k2 == k || norm(e) + d.r - c.r <= eps {
                continue;
            }
            let w = if d.r - c.r >= norm(e) { PI } else { ((c.r - d.r) / norm(e)).clamp(-1.0, 1.0).acos() };
            for (lo, hi) in part(angle(e) - w, 2.0 * w, t0, sw) {
                // Within k2's sweep.
                for (l2, h2) in part(span[k2].0, span[k2].1, t0, sw) {
                    if lo.max(l2) < hi.min(h2) {
                        beat.push((lo.max(l2), hi.min(h2), To::Cut(k2)));
                    }
                }
            }
        }
        // The ends of the intervals' union inside the sweep.
        beat.sort_by(|x, y| x.0.total_cmp(&y.0));
        let gap = sag_turn(c.r, 10.0 * eps);
        let mut at = |u: f64, to: To| {
            if u > gap && u < sw - gap {
                out.push((k, t0 + u, to));
            }
        };
        let mut run: Option<(f64, To)> = None;
        for &(lo, hi, to) in &beat {
            match run {
                Some((h, _)) if lo <= h => {
                    if hi > h {
                        run = Some((hi, to));
                    }
                }
                _ => {
                    if let Some((h, t)) = run {
                        at(h, t);
                    }
                    at(lo, to);
                    run = Some((hi, to));
                }
            }
        }
        if let Some((h, t)) = run {
            at(h, t);
        }
    }
    // A common tangent's other contact (on arc j, at the same normal); if it
    // is within its arc's gap, that arc's nearer end.
    out.iter()
        .map(|&(k, t, to)| Cut {
            e: arcs[k].0,
            n: [t.cos(), t.sin()],
            to: match to {
                To::Pt(i) => To::Pt(i),
                To::Cut(j) => out
                    .iter()
                    .position(|&(o, ot, oto)| o == j && matches!(oto, To::Cut(x) if x == k) && ((ot - t + PI).rem_euclid(TAU) - PI).abs() <= 1e-9)
                    .map_or_else(|| To::Pt(if (t - span[j].0).rem_euclid(TAU) < span[j].1 / 2.0 { arcs[j].2 } else { arcs[j].3 }), To::Cut),
            },
        })
        .collect()
}

use super::*;
use ink_core::interval::I;
use walls::Walls;

/// How to write a profile whose corners move (a mould draft: walls that
/// tilt along z).
const MOVING: &str = "constants and params work here; x, y, z and t do not. A profile whose walls lean with z (a mould's draft) \
                      is its profile at z = 0 drafted: (polygon pts drafts z0 z1), each wall leaning out by its angle (degrees; one \
                      for all, or one an edge) over z0..z1; a corner rounded (vec u v bulge round) keeps its radius as its walls lean";

/// How a drafted piece's leaves move: t over its window (0 at z0, 1 at
/// z1), z clamped to its part of the span (t over [ta, tb]), and x and y to
/// its box grown by its size (for the lines that turn).
struct Span<'a> {
    z0: f64,
    z1: f64,
    ta: f64,
    tb: f64,
    zc: NodeId,
    xc: NodeId,
    yc: NodeId,
    lo: P,
    hi: P,
    mv: &'a walls::Moving,
}

/// A drafted profile whose field has no bound on its rate somewhere.
fn unbounded(sp: &Span) -> String {
    format!(
        "the drafted polygon over z = {}..{}: a curve's rate has no bound over part of the span (its walls meet, or a hull line touches its arc, too nearly tangentially there): give them one draft, or shorten the span",
        sp.z0, sp.z1
    )
}

/// Which side of 0 a field of (x, y) is on wherever x >= 0, proven: `Below`
/// is at most 0 there and 0 at most on the axis x = 0 (an edge on the axis,
/// -k x + c with k > 0 and c <= 0), `Above` the same of its negation. A
/// `Below` term of a max, or an `Above` one of a min, decides nothing off
/// the axis, and on it reads 0 where the max (min) is inside (outside).
#[derive(Clone, Copy, PartialEq, Debug)]
enum Side {
    Below,
    Above,
}

impl Side {
    fn flip(self) -> Side {
        if self == Side::Below { Side::Above } else { Side::Below }
    }
}

/// `revolve-z`'s rewrite of one profile: its terms' coefficients, sides and
/// rewrites, each found once per node (a profile's expressions are shared).
#[derive(Default)]
struct OffAxis {
    affine: HashMap<NodeId, Option<[I; 4]>>,
    side: HashMap<NodeId, Option<Side>>,
    done: HashMap<NodeId, NodeId>,
}

/// A constant above 0 (a scale), if `n` is one.
fn positive(g: &ink_core::graph::Graph, n: NodeId) -> Option<f64> {
    g.as_const(n).filter(|&k| k > 0.0 && k.is_finite())
}

impl OffAxis {
    /// `node` as a + b x + c y + d z of constants alone, if it is one: each
    /// coefficient an interval holding its value in real numbers (outward
    /// where the arithmetic rounds), so a sign or a zero read from it is
    /// proven.
    fn affine(&mut self, g: &ink_core::graph::Graph, node: NodeId) -> Option<[I; 4]> {
        use ink_core::graph::{BinOp as B, Node, UnOp as U};
        use ink_core::interval::{add_o, bin, mul_o, sub_o};
        if let Some(&a) = self.affine.get(&node) {
            return a;
        }
        let (o, l) = (I::point(0.0), I::point(1.0));
        let a = (|| {
            Some(match g.nodes[node as usize] {
                Node::Const(v) => [I::point(v), o, o, o],
                Node::Var(Var::X) => [o, l, o, o],
                Node::Var(Var::Y) => [o, o, l, o],
                Node::Var(Var::Z) => [o, o, o, l],
                Node::Un(U::Neg, a) => self.affine(g, a)?.map(|v| I::new(-v.hi, -v.lo)),
                Node::Tag(_, a) | Node::Real(a) => self.affine(g, a)?,
                Node::Bin(op @ (B::Add | B::Sub), a, b) => {
                    let (a, b) = (self.affine(g, a)?, self.affine(g, b)?);
                    [0, 1, 2, 3].map(|i| if op == B::Add { add_o(a[i], b[i]) } else { sub_o(a[i], b[i]) })
                }
                Node::Bin(B::Mul, a, b) => match (g.as_const(a), g.as_const(b)) {
                    (Some(k), _) => self.affine(g, b)?.map(|v| mul_o(I::point(k), v)),
                    (_, Some(k)) => self.affine(g, a)?.map(|v| mul_o(I::point(k), v)),
                    _ => return None,
                },
                Node::Bin(B::Div, a, b) => {
                    let k = g.as_const(b).filter(|&k| k != 0.0)?;
                    self.affine(g, a)?.map(|v| bin(B::Div, v, I::point(k)).0)
                }
                _ => return None,
            })
        })();
        self.affine.insert(node, a);
        a
    }

    /// The side `node` is on wherever x >= 0 (`Side`), proven: of an edge
    /// on the axis, and through names, promises, neg, min, max, a scale
    /// and a select whose branches are on one side.
    fn side(&mut self, g: &ink_core::graph::Graph, node: NodeId) -> Option<Side> {
        use ink_core::graph::{BinOp as B, Node, UnOp as U};
        if let Some(&s) = self.side.get(&node) {
            return s;
        }
        let zero = |v: I| v.lo == 0.0 && v.hi == 0.0;
        let s = match self.affine(g, node) {
            Some([c, kx, ky, kz]) if zero(ky) && zero(kz) => {
                if kx.hi < 0.0 && c.hi <= 0.0 {
                    Some(Side::Below)
                } else if kx.lo > 0.0 && c.lo >= 0.0 {
                    Some(Side::Above)
                } else {
                    None
                }
            }
            Some(_) => None,
            None => match g.nodes[node as usize] {
                Node::Un(U::Neg, a) => self.side(g, a).map(Side::flip),
                Node::Tag(_, a) | Node::Lipschitz(a, ..) => self.side(g, a),
                // The max of two terms is below if both are, above if either
                // is (0 only where that one is: on the axis); the min the
                // other way round.
                Node::Bin(op @ (B::Max | B::Min), a, b) => {
                    let (sa, sb) = (self.side(g, a), self.side(g, b));
                    let (all, any) = if op == B::Max { (Side::Below, Side::Above) } else { (Side::Above, Side::Below) };
                    if sa == Some(any) || sb == Some(any) {
                        Some(any)
                    } else if sa == Some(all) && sb == Some(all) {
                        Some(all)
                    } else {
                        None
                    }
                }
                Node::Bin(B::Mul, a, b) => match (positive(g, a), positive(g, b)) {
                    (Some(_), _) => self.side(g, b),
                    (_, Some(_)) => self.side(g, a),
                    _ => None,
                },
                Node::Bin(B::Div, a, b) if positive(g, b).is_some() => self.side(g, a),
                Node::Select(_, a, b) => Some(self.side(g, a)?).filter(|&s| self.side(g, b) == Some(s)),
                _ => None,
            },
        };
        self.side.insert(node, s);
        s
    }

    /// `a` without the terms on the axis (`Side`) of its maxes and mins,
    /// once per node (a shared expression is rewritten once), through what
    /// keeps a sign: names, promises (dropped where the field changes),
    /// neg, a scale, and a select (its test too, so it picks the same
    /// branch off the axis). Off the axis every sign is as it was.
    fn rewrite(&mut self, g: &mut ink_core::graph::Graph, a: NodeId) -> NodeId {
        use ink_core::graph::{BinOp as B, Node, UnOp as U};
        if let Some(&n) = self.done.get(&a) {
            return n;
        }
        let n = match g.nodes[a as usize] {
            // A term on the axis (a group of them: `side` classes it) is
            // left out, unless both are: then the max (min) is one, for
            // what holds it to leave out.
            Node::Bin(op @ (B::Max | B::Min), p, q) => {
                let off = Some(if op == B::Max { Side::Below } else { Side::Above });
                match (self.side(g, p) == off, self.side(g, q) == off) {
                    (true, true) => a,
                    (true, false) => self.rewrite(g, q),
                    (false, true) => self.rewrite(g, p),
                    (false, false) => {
                        let (p, q) = (self.rewrite(g, p), self.rewrite(g, q));
                        g.bin(op, p, q)
                    }
                }
            }
            Node::Un(U::Neg, p) => {
                let p = self.rewrite(g, p);
                g.neg(p)
            }
            Node::Tag(t, p) => {
                let p = self.rewrite(g, p);
                let name = g.tag_names[t as usize].clone();
                g.tag(&name, p)
            }
            // A promise made of the whole holds of it, not of fewer terms.
            Node::Lipschitz(p, _, _) => {
                let q = self.rewrite(g, p);
                if q == p { a } else { q }
            }
            Node::Bin(op @ (B::Mul | B::Div), p, q) if positive(g, q).is_some() => {
                let p = self.rewrite(g, p);
                g.bin(op, p, q)
            }
            Node::Bin(B::Mul, p, q) if positive(g, p).is_some() => {
                let q = self.rewrite(g, q);
                g.mul(p, q)
            }
            Node::Select(c, p, q) => {
                let (c, p, q) = (self.rewrite(g, c), self.rewrite(g, p), self.rewrite(g, q));
                g.select(c, p, q)
            }
            _ => a,
        };
        self.done.insert(a, n);
        n
    }
}

impl Doc {
    /// `(revolve-z g)`: the 2D profile g(x, y) turned about the z axis, its
    /// x the radius and its y the height. A profile reaching the axis has
    /// an edge on it, a term -k x (k > 0) of a max: read at the radius r,
    /// never below 0, it is below 0 but on the axis, where it reads 0 (no
    /// surface: the axis is inside). So the profile is read without its
    /// terms on the axis (`OffAxis::rewrite`): off the axis every sign stays
    /// as it was; values change only inside, near an edge left out, where
    /// they read the walls that remain.
    pub(super) fn eval_revolve(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let [g] = args else {
            return Err(self.err_hint(
                e,
                "revolve-z takes one profile",
                Some("usage: (revolve-z g): the 2D profile g(x, y) turned about the z axis, x the radius and y the height".into()),
            ));
        };
        let g = self.field(g, env)?;
        let g = OffAxis::default().rewrite(&mut self.g, g);
        let (x, y, z) = (self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z));
        let (xx, yy) = (self.g.un(UnOp::Square, x), self.g.un(UnOp::Square, y));
        let s = self.g.add(xx, yy);
        let r = self.g.un(UnOp::Sqrt, s);
        let zero = self.g.constant(0.0);
        Ok(Value::Field(self.g.subst(g, &[(Var::X, r), (Var::Y, z), (Var::Z, zero)]), None))
    }
}

impl Doc {
    /// `(polygon pts [drafts z0 z1])`: pts are (vec u v) corners, (vec u v
    /// b) where the edge to the next point is an arc of bulge b, or (vec u v
    /// b r), the corner rounded to radius r. Drafted, each edge's wall leans
    /// out by its angle over z0..z1 (walls.rs). Its numbers are read at load
    /// time (their params are baked: derivatives through them are taken by
    /// differences).
    pub(super) fn eval_polygon(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let usage = "usage: (polygon (vec (vec u v) ...) [drafts z0 z1]): a closed profile in the xy plane through those corners; \
                     a point (vec u v b) makes the edge to the next point an arc of bulge b = tan(sweep / 4), \
                     positive turning counter-clockwise, and (vec u v b r) rounds the corner there to radius r. \
                     Drafted, the profile is pts at z = 0 with each edge's wall leaning out by its angle (degrees: one for all \
                     edges, or a vec, one an edge, from each point to the next) as z goes from 0, solved over z0..z1 and \
                     straight beyond; its rounds keep their radius";
        let (l, draft) = match args {
            [l] => (l, None),
            [l, d, z0, z1] => (l, Some((d, z0, z1))),
            _ => return Err(self.err_hint(e, "polygon takes a list of points, and drafts with their span", Some(usage.into()))),
        };
        let items = self.list(l, env)?;
        let mut pts = Vec::with_capacity(items.len());
        for it in &items {
            let Value::Vec(c) = it else {
                return Err(self.err_hint(l, format!("a polygon's point is a vec, not a {}", it.kind()), Some(usage.into())));
            };
            if !(2..=4).contains(&c.len()) {
                return Err(self.err_hint(l, format!("a polygon's point is (vec u v), (vec u v bulge) or (vec u v bulge round), not {} numbers", c.len()), Some(usage.into())));
            }
            let mut p = [0.0; 4];
            for (k, v) in c.iter().enumerate() {
                let Value::Field(f, _) = v else { return Err(self.err_hint(l, "a polygon's point holds numbers", Some(usage.into()))) };
                p[k] = self.load_number(l, *f)?.ok_or_else(|| self.err_hint(l, "a polygon's points must be known at load time", Some(MOVING.into())))?;
            }
            pts.push(p);
        }
        let bad = |d: &mut Doc, m: String| d.err_hint(e, m, Some(usage.into()));
        let Some((d, z0, z1)) = draft else {
            if pts.iter().all(|p| p[3] == 0.0) {
                let pts: Vec<[f64; 3]> = pts.iter().map(|p| [p[0], p[1], p[2]]).collect();
                let t = polygon(&pts).map_err(|m| bad(self, m))?;
                let f = self.lower_profile(&t, None).map_err(|m| bad(self, m))?;
                return Ok(Value::Field(f, None));
            }
            let t = Walls::new(&pts, &vec![0.0; pts.len()]).and_then(|w| polygon(&w.solve(0.0)?.0)).map_err(|m| bad(self, m))?;
            let f = self.lower_profile(&t, None).map_err(|m| bad(self, m))?;
            return Ok(Value::Field(f, None));
        };
        let angles: Vec<f64> = match self.eval(d, env)? {
            Value::Vec(v) => {
                if v.len() != pts.len() {
                    return Err(self.err_hint(
                        d,
                        format!("a polygon of {} points has {} edges, so {} drafts, not {}", pts.len(), pts.len(), pts.len(), v.len()),
                        Some(usage.into()),
                    ));
                }
                let mut a = Vec::with_capacity(v.len());
                for x in &v {
                    let Value::Field(f, _) = x else { return Err(self.err_hint(d, "a polygon's drafts are numbers (degrees)", Some(usage.into()))) };
                    a.push(self.load_number(d, *f)?.ok_or_else(|| self.err_hint(d, "a polygon's drafts must be known at load time", Some(MOVING.into())))?);
                }
                a
            }
            Value::Field(f, _) => vec![self.load_number(d, f)?.ok_or_else(|| self.err_hint(d, "a polygon's drafts must be known at load time", Some(MOVING.into())))?; pts.len()],
            other => return Err(self.err_hint(d, format!("a polygon's drafts are a number or a vec of numbers, not a {}", other.kind()), Some(usage.into()))),
        };
        if let Some(a) = angles.iter().find(|a| !(a.abs() < 90.0)) {
            return Err(self.err_hint(d, format!("a draft is an angle in degrees within (-90, 90), not {a}"), Some(usage.into())));
        }
        let (z0, z1) = (self.const_num(z0, env)?, self.const_num(z1, env)?);
        if !(z0 != z1 && z0.is_finite() && z1.is_finite()) {
            return Err(self.err_hint(e, format!("a drafted polygon's span z0..z1 must be two different heights, not {z0}..{z1}"), Some(usage.into())));
        }
        let (z0, z1) = (z0.min(z1), z0.max(z1));
        let slopes: Vec<f64> = angles.iter().map(|a| a.to_radians().tan()).collect();
        let w = Walls::new(&pts, &slopes).map_err(|m| bad(self, m))?;
        let d = w.drafted(z0, z1).map_err(|m| bad(self, m))?;
        let (x, y, z) = (self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z));
        let clamp = |d: &mut Doc, v: NodeId, lo: f64, hi: f64| {
            let (lo, hi) = (d.g.constant(lo), d.g.constant(hi));
            let m = d.g.min(v, hi);
            d.g.max(m, lo)
        };
        // A line that turns reads x and y clamped to the box grown by its
        // size, so its slope along z is bounded, and its value is the line's
        // wherever an offset or a fillet of the profile reads it.
        let grow = norm(sub(d.hi, d.lo));
        let (blo, bhi) = ([d.lo[0] - grow, d.lo[1] - grow], [d.hi[0] + grow, d.hi[1] + grow]);
        let (xc, yc) = (clamp(self, x, blo[0], bhi[0]), clamp(self, y, blo[1], bhi[1]));
        // Each piece cut to its part of the span by a slab (but beyond the
        // span's ends, where the profile is straight), the least of them.
        let last = d.pieces.len() - 1;
        let mut fs = Vec::with_capacity(d.pieces.len());
        for (k, pc) in d.pieces.iter().enumerate() {
            let [w0, w1] = pc.mv.z;
            let t = |z: f64| (z - w0) / (w1 - w0);
            let sp = Span { z0: w0, z1: w1, ta: t(pc.held[0]), tb: t(pc.held[1]), zc: clamp(self, z, pc.held[0], pc.held[1]), xc, yc, lo: blo, hi: bhi, mv: &pc.mv };
            let mut f = self.lower_profile(&pc.t, Some(&sp)).map_err(|m| bad(self, m))?;
            if k > 0 {
                let r0 = self.g.constant(pc.range[0]);
                let below = self.g.sub(r0, z);
                f = self.g.max(f, below);
            }
            if k < last {
                let r1 = self.g.constant(pc.range[1]);
                let above = self.g.sub(z, r1);
                f = self.g.max(f, above);
            }
            fs.push(f);
        }
        // Beyond its ends the profile is straight: there it is also its end
        // section's flat polygon (whose leaves the walls' rates do not scale
        // down, as a line turning fast near where a construction changes
        // does its piece's), cut to the heights past that end.
        for (h, past) in [(d.pieces[0].held[0], -1.0), (d.pieces[last].held[1], 1.0)] {
            if let Ok(t) = w.solve(h).and_then(|(q, _)| polygon(&q)) {
                let end = self.lower_profile(&t, None).map_err(|m| bad(self, m))?;
                let cut = self.ramp(past * h, -past, z);
                fs.push(self.g.max(end, cut));
            }
        }
        let f = self.reduce(BinOp::Min, &fs);
        let fixed = T::Max(vec![T::Lin(1.0, 0.0, -d.hi[0], None), T::Lin(-1.0, 0.0, d.lo[0], None), T::Lin(0.0, 1.0, -d.hi[1], None), T::Lin(0.0, -1.0, d.lo[1], None)]);
        let b = self.lower_profile(&fixed, None).map_err(|m| bad(self, m))?;
        Ok(Value::Field(self.g.max(f, b), None))
    }

    /// A line the construction adds through points moving with the walls,
    /// or from a centre moving with them (curve.rs): n(t) . (p - c(t)), its
    /// unit normal n (m / |m|) and c a point on it (its segment's middle or
    /// an end, whichever bounds it least; or the point a tangent is drawn
    /// at), the points' paths lowered as the proof reads them: the distance
    /// across it. Its rate along z at p is how
    /// fast it turns times how far p is from c, and how fast c moves across
    /// it; each is bounded in intervals over parts of the piece's span
    /// (`curve::most`; refused where one has none), and the line is scaled by
    /// the least of its slope's bounds. A line whose normal strays from its
    /// middle one by under polygon's allowance over the box (it does not
    /// turn: a round's chord, or one through points sliding along it) is read
    /// with that normal fixed, n0 . p - n0 . c(t): a plane where c moves
    /// linearly. One that turns is read over the box (x and y clamped to it,
    /// as `Span` says), or, where it turns so fast that it would read weaker
    /// so, in full only within rho of c, and beyond as n . u rho (u the
    /// direction from c): its sign everywhere, its slope across it at most 1,
    /// its rate along z at most rho |n'| and how fast c moves (across it
    /// within rho; beyond, at r from c, rho / r of how fast it moves across
    /// u, which is at most its speed across n and half along it). rho makes
    /// the turn's part of the rate 1.
    fn curved_line(&mut self, line: Normal, sp: &Span) -> Result<NodeId, String> {
        use curve::{Ex, padd, pcross, pdot, pnorm, pscale, psub, xdiv};
        use ink_core::interval::I;
        let pt = |i: usize| sp.mv.pts[i].clone();
        let centre = |e: usize| sp.mv.circles[e].0.clone();
        // Its normal, and the points on it it may be read about (c).
        let (m, cs) = match line {
            Normal::Line(i, j) => {
                let e = psub(&pt(j), &pt(i));
                let mid = pscale(&padd(&pt(i), &pt(j)), &curve::kx(0.5));
                ([e[1].clone(), psub(&[curve::kx(0.0), curve::kx(0.0)], &e)[0].clone()], vec![mid, pt(i), pt(j)])
            }
            Normal::Arc(e, i) => (psub(&pt(i), &centre(e)), vec![pt(i)]),
            Normal::Toward(e, i) => (psub(&centre(e), &pt(i)), vec![pt(i)]),
        };
        let c = cs[0].clone();
        let n = pscale(&m, &xdiv(&curve::kx(1.0), &pnorm(&m)));
        let dz = sp.z1 - sp.z0;
        let sl = |x: &curve::X, t: I| x.at(curve::Sl::of(t)).e;
        // The most of |x| over 4 equal parts of the piece's span.
        let parts = |f: &dyn Fn(I) -> I| {
            let cut = |k: usize| sp.ta + (sp.tb - sp.ta) * k as f64 / 4.0;
            (0..4).map(|k| curve::most(f, cut(k), cut(k + 1))).try_fold(0.0f64, |a, b| b.map(|b| a.max(b)))
        };
        let most = |x: &curve::X| parts(&|t| sl(x, t));
        let ok = |l: f64| l.is_finite() && l > 0.0;
        let t = self.t_of(sp);
        let mut memo = HashMap::new();
        let (x, y) = (self.g.var(Var::X), self.g.var(Var::Y));
        // m = m0 + m1 t whose direction strays from its middle one by under
        // the allowance over the box: it does not turn.
        let lin = |x: &curve::X| match **x {
            Ex::K(v) => Some((v, 0.0)),
            Ex::Lin(a, b) => Some((a, b)),
            _ => None,
        };
        if let (Some(mx), Some(my), Some(cx), Some(cy)) = (lin(&m[0]), lin(&m[1]), lin(&c[0]), lin(&c[1])) {
            let at = |t: f64| [mx.0 + mx.1 * t, my.0 + my.1 * t];
            let mid = at(0.5 * (sp.ta + sp.tb));
            let n0 = unit(mid);
            // |n0 x m(t)| is linear in t, |m(t)| least where t is nearest its foot.
            let side = cross(n0, at(sp.ta)).abs().max(cross(n0, at(sp.tb)).abs());
            let d1 = [mx.1, my.1];
            let foot = if dot(d1, d1) > 0.0 { (-dot([mx.0, my.0], d1) / dot(d1, d1)).clamp(sp.ta, sp.tb) } else { sp.ta };
            let least = norm(at(foot)).min(norm(at(sp.ta))).min(norm(at(sp.tb)));
            let reach = norm(sub(sp.hi, sp.lo));
            // Over the box n0 . (p - c) strays from n . (p - c) by under
            // |n - n0| reach, |n - n0| = sin / cos(half the angle) and sin =
            // |n0 x m| / |m| <= side / least: with the angle under 1e-7
            // (the allowance over the box), 1 + 1e-9 bounds the cosine's
            // part and the floats' rounding of this test.
            if n0.iter().all(|v| v.is_finite()) && least > 0.0 && side / least * reach * (1.0 + 1e-9) <= 100.0 * sp.mv.eps() {
                // n0 . p - o(t), o = n0 . c = o0 + o1 t: a plane.
                let (o0, o1) = (n0[0] * cx.0 + n0[1] * cy.0, n0[0] * cx.1 + n0[1] * cy.1);
                let k = 1.0 / (o1 / dz).hypot(1.0);
                return Ok(self.lin3([n0[0] * k, n0[1] * k, -o1 / dz * k, -(o0 - o1 * sp.z0 / dz) * k], [x, y, sp.zc]));
            }
        }
        // Its turn (m x m' / |m|^2, per unit t).
        let dm = [Ex::d(&m[0]), Ex::d(&m[1])];
        let lm = pnorm(&m);
        let turn = most(&xdiv(&pcross(&m, &dm), &pdot(&m, &m)));
        // Over the box (within its diagonal of c), at its full slope across
        // it (|n| is 1, to rounding).
        let reach = norm(sub(sp.hi, sp.lo));
        // Read about each c: in full over the box, or in full within rho of c
        // (1 / its turn per unit z) and saturated beyond: within rho its rate
        // is rho times its turn and how fast c moves across it, beyond at most
        // that and half how fast c moves at all. The least of the slope bounds.
        let read = |c: &curve::XP| {
            let dc = [Ex::d(&c[0]), Ex::d(&c[1])];
            let (across, speed) = most(&xdiv(&pdot(&m, &dc), &lm)).zip(most(&pnorm(&dc)))?;
            let w = turn?;
            let boxed = ((w * reach + across) / dz).hypot(1.0);
            let rho = if w > 0.0 { dz / w } else { f64::INFINITY };
            let sat = ((w * rho + across.max(0.5 * (across + speed))) / dz).hypot(1.0);
            Some((boxed, rho, sat))
        };
        let (c, (boxed, rho, sat)) =
            cs.iter().filter_map(|c| read(c).map(|v| (c.clone(), v))).min_by(|a, b| a.1.0.min(a.1.2).total_cmp(&b.1.0.min(b.1.2))).ok_or_else(|| unbounded(sp))?;
        let boxed = Some(boxed).filter(|&l| ok(l));
        let sat = Some((rho, sat)).filter(|&(r, l)| ok(r) && ok(l));
        // n . (p - c), p clamped to the box or not.
        let (nx, ny) = (self.lower_ex(&n[0], t, &mut memo), self.lower_ex(&n[1], t, &mut memo));
        let (cx, cy) = (self.lower_ex(&c[0], t, &mut memo), self.lower_ex(&c[1], t, &mut memo));
        let across = |d: &mut Doc, px: NodeId, py: NodeId| {
            let (dx, dy) = (d.g.sub(px, cx), d.g.sub(py, cy));
            let (ax, ay) = (d.g.mul(nx, dx), d.g.mul(ny, dy));
            (d.g.add(ax, ay), dx, dy)
        };
        match (boxed, sat) {
            (Some(l), s) if s.is_none_or(|s| l <= s.1) => {
                let (f, _, _) = across(self, sp.xc, sp.yc);
                let k = self.g.constant(1.0 / l);
                Ok(self.g.mul(f, k))
            }
            (_, Some((rho, l))) => {
                let (f, dx, dy) = across(self, x, y);
                let (sx, sy) = (self.g.un(UnOp::Square, dx), self.g.un(UnOp::Square, dy));
                let q = self.g.add(sx, sy);
                let r = self.g.un(UnOp::Sqrt, q);
                // rho / r beyond rho.
                let kr = self.g.constant(rho);
                let far = self.g.max(r, kr);
                let w = self.g.bin(BinOp::Div, kr, far);
                let f = self.g.mul(f, w);
                let k = self.g.constant(1.0 / l);
                let f = self.g.mul(f, k);
                let one = self.g.constant(1.0);
                Ok(self.g.lipschitz(f, one))
            }
            _ => Err(unbounded(sp)),
        }
    }

    /// A drafted profile's disc, its circle as the proof reads it (a wall's
    /// cone, a round's centre moving linearly, or on a curve where its
    /// walls, offset by its radius, meet and one is an arc: curve.rs):
    /// |p - c(t)| - r(t), scaled so its gradient is at most 1: by its rate
    /// along z, its centre's and radius's, bounded in intervals over parts
    /// of the span (`curve::most`; refused where it has no bound).
    fn moving_disc(&mut self, c: &curve::XP, r: &curve::X, sp: &Span) -> Result<NodeId, String> {
        use curve::{Ex, Val};
        use ink_core::interval::{I, mag};
        let (dc, dr) = ([Ex::d(&c[0]), Ex::d(&c[1])], Ex::d(r));
        let dz = sp.z1 - sp.z0;
        let iv = <I as Val>::bin;
        let sl = |x: &curve::X, t: I| x.at(curve::Sl::of(t)).e;
        let rate = curve::most(
            |t| {
                let d = [sl(&dc[0], t), sl(&dc[1], t)];
                let speed = <I as Val>::sqrt(iv(BinOp::Add, iv(BinOp::Mul, d[0], d[0]), iv(BinOp::Mul, d[1], d[1])));
                iv(BinOp::Add, speed, I::point(mag(sl(&dr, t))))
            },
            sp.ta,
            sp.tb,
        )
        .ok_or_else(|| unbounded(sp))?;
        let (x, y) = (self.g.var(Var::X), self.g.var(Var::Y));
        let t = self.t_of(sp);
        let mut memo = HashMap::new();
        let (cx, cy, rr) = (self.lower_ex(&c[0], t, &mut memo), self.lower_ex(&c[1], t, &mut memo), self.lower_ex(r, t, &mut memo));
        let (dx, dy) = (self.g.sub(x, cx), self.g.sub(y, cy));
        let (sx, sy) = (self.g.un(UnOp::Square, dx), self.g.un(UnOp::Square, dy));
        let q = self.g.add(sx, sy);
        let d = self.g.un(UnOp::Sqrt, q);
        let mut f = self.g.sub(d, rr);
        if rate != 0.0 {
            let s = self.g.constant(1.0 / (rate / dz).hypot(1.0));
            f = self.g.mul(f, s);
        }
        let one = self.g.constant(1.0);
        Ok(self.g.lipschitz(f, one))
    }

    /// a x + b y + c z + d (c z left out where c is 0).
    fn lin3(&mut self, [a, b, c, d]: [f64; 4], [x, y, z]: [NodeId; 3]) -> NodeId {
        let (a, b, d) = (self.g.constant(a), self.g.constant(b), self.g.constant(d));
        let (ax, by) = (self.g.mul(a, x), self.g.mul(b, y));
        let mut s = self.g.add(ax, by);
        if c != 0.0 {
            let c = self.g.constant(c);
            let cz = self.g.mul(c, z);
            s = self.g.add(s, cz);
        }
        self.g.add(s, d)
    }

    /// A piece's t at height z (clamped to its held part): (z - z0) / (z1 -
    /// z0), the length rounded as the proof's height of t rounds it
    /// (`Moving::z_of_t`), in the graph's operations (each real, its range
    /// rounded outward), so that the field reads each path at the height
    /// the proof pairs with it, exactly.
    fn t_of(&mut self, sp: &Span) -> NodeId {
        let (z0, dz) = (self.g.constant(sp.z0), self.g.constant(sp.z1 - sp.z0));
        let d = self.g.sub(sp.zc, z0);
        self.g.bin(BinOp::Div, d, dz)
    }

    /// v0 + k z.
    fn ramp(&mut self, v0: f64, k: f64, z: NodeId) -> NodeId {
        let v0 = self.g.constant(v0);
        if k == 0.0 {
            return v0;
        }
        let k = self.g.constant(k);
        let kz = self.g.mul(k, z);
        self.g.add(kz, v0)
    }

    /// The tree as fields: its leaves, or with a span, each moving with
    /// the piece's walls, scaled so its gradient is at most 1: a line
    /// through two points of a line wall is the wall's plane, tilted by its
    /// draft; a disc, its circle (a cone, or a round's moving centre); a
    /// line the construction adds, through points or from a centre moving
    /// differently, turns (`curved_line`).
    fn lower_profile(&mut self, a: &T, sp: Option<&Span>) -> Result<NodeId, String> {
        let (x, y) = (self.g.var(Var::X), self.g.var(Var::Y));
        match a {
            &T::Lin(a0, b0, c0, line) => {
                let Some(sp) = sp.filter(|_| line.is_some()) else {
                    return Ok(self.lin3([a0, b0, 0.0, c0], [x, y, x]));
                };
                let line = line.expect("a drafted line is the loop's");
                if let Normal::Line(i, j) = line
                    && let Some(&w) = sp.mv.on_walls[i].iter().find(|w| sp.mv.on_walls[j].contains(w))
                {
                    // n . p = k + s z, its normal the leaf's way round.
                    let (n, k, s) = sp.mv.planes[w].expect("a line wall");
                    let (sg, h) = ((a0 * n[0] + b0 * n[1]).signum(), 1.0 / s.hypot(1.0));
                    return Ok(self.lin3([sg * n[0] * h, sg * n[1] * h, -sg * s * h, -sg * k * h], [x, y, sp.zc]));
                }
                self.curved_line(line, sp)
            }
            &T::Disc(cx, cy, r, edge) => {
                if let Some(sp) = sp {
                    let (c, r) = &sp.mv.circles[edge.expect("a drafted disc is a loop edge's")];
                    return self.moving_disc(c, r, sp);
                }
                let (cx, cy, r) = (self.g.constant(cx), self.g.constant(cy), self.g.constant(r));
                let (dx, dy) = (self.g.sub(x, cx), self.g.sub(y, cy));
                let (sx, sy) = (self.g.un(UnOp::Square, dx), self.g.un(UnOp::Square, dy));
                let q = self.g.add(sx, sy);
                let d = self.g.un(UnOp::Sqrt, q);
                let f = self.g.sub(d, r);
                let one = self.g.constant(1.0);
                Ok(self.g.lipschitz(f, one))
            }
            T::Neg(a) => {
                let a = self.lower_profile(a, sp)?;
                Ok(self.g.neg(a))
            }
            T::Max(v) | T::Min(v) => {
                // Nested maxes (mins) flattened into one balanced tree.
                let max = matches!(a, T::Max(_));
                fn flat<'a>(a: &'a T, max: bool, out: &mut Vec<&'a T>) {
                    match a {
                        T::Max(v) if max => v.iter().for_each(|a| flat(a, max, out)),
                        T::Min(v) if !max => v.iter().for_each(|a| flat(a, max, out)),
                        _ => out.push(a),
                    }
                }
                let mut leaves = Vec::new();
                v.iter().for_each(|a| flat(a, max, &mut leaves));
                let mut ids: Vec<NodeId> = Vec::with_capacity(leaves.len());
                for a in leaves {
                    let id = self.lower_profile(a, sp)?;
                    if !ids.contains(&id) {
                        ids.push(id);
                    }
                }
                Ok(self.reduce(if max { BinOp::Max } else { BinOp::Min }, &ids))
            }
        }
    }
}

mod curve;
mod walls;

#[cfg(test)]
mod tests;
