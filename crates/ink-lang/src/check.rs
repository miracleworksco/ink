//! `(check TEST VALUE)`: a condition a builder states on its inputs.
//! VALUE is returned unchanged; TEST is kept, with where it was stated, and
//! `ink check` reports a failing one as a warning (Doc::check_failures), at
//! the params' values and moment it is run at: every check loading ran,
//! in the checked file and in what it imports. Two tests:
//!
//! - `(value E CMP V)`: a number of params and t, compared as a value claim;
//! - `(path-thickness PTS >= R)`: the thickness of the polyline PTS (its
//!   reach, `path_thickness`). A body made by moving a section along a path
//!   is free of self-overlap when the thickness is at least the section's
//!   reach from the path.

use super::*;
use ink_core::v3::{dot, norm, sub};
use ink_kernel::ranged::Cmp;

/// A path whose first and last points are closer than this (mm) is a loop:
/// std.ink's path-closes, which its sweeps select their ends by.
pub const LOOP_GAP: f64 = 1e-9;

/// A condition stated by `check`.
#[derive(Clone, Debug)]
pub struct Check {
    pub test: Test,
    /// The test as written.
    pub text: String,
    /// The top-level form it was stated under (file:line (head name ...)),
    /// and the calls it was stated in, outermost first: `name (file:line)`.
    pub site: String,
    pub calls: Vec<String>,
}

#[derive(Clone, Debug)]
pub enum Test {
    Value {
        e: NodeId,
        cmp: Cmp,
        v: NodeId,
    },
    /// The path's points (x, y, z; z is 0 for a path in the plane).
    PathThickness {
        pts: Vec<[NodeId; 3]>,
        v: NodeId,
    },
}

impl Test {
    fn nodes(&self) -> Vec<NodeId> {
        match self {
            Test::Value { e, v, .. } => vec![*e, *v],
            Test::PathThickness { pts, v, .. } => pts.iter().flatten().copied().chain([*v]).collect(),
        }
    }

    /// The kind of test and its comparison, which with its nodes say what
    /// it tests.
    fn kind(&self) -> u8 {
        match self {
            Test::Value { cmp, .. } => 1 + *cmp as u8,
            Test::PathThickness { .. } => 0,
        }
    }
}

/// Where a path is thinnest: `radius` is half the distance between two
/// stretches (at stations `s` along it, at `at`), or the reach at which the
/// mitred ends of a segment cross (one station, the segment's middle).
#[derive(Clone, Debug, PartialEq)]
pub struct Thin {
    pub radius: f64,
    pub s: Vec<f64>,
    pub at: Vec<V3>,
}

/// How thin a path is where it is thinner than `cap`, if anywhere: the
/// tightest bend and the closest pair of distant stretches, and how many
/// segment pairs come closer than 2 cap. A path whose first and last points
/// coincide is a loop: it has no ends, and it bends where it closes.
///
/// The thickness of a curve (Gonzalez and Maddocks' global radius of
/// curvature; for polygons, Rawdon's) is the least of its radius of
/// curvature and half the distance between stretches that face each other:
/// a pair of points whose chord is perpendicular to the path at both ends
/// (doubly critical), which excludes pairs that are merely near along a
/// bend. A tube of radius r about the curve overlaps itself iff the
/// thickness is below r; a section of reach r (its farthest point from the
/// path) is then clear when the thickness is at least r.
///
/// A polyline bends at its joints, which std's sweeps mitre (the plane
/// halving the turn). Segment i, of length L along u, holds its two mitres
/// apart within reach r iff L >= r |m|, where m = n1 - n0, the mitre at its
/// end n1 = perp(u') / (1 + u . u') for the next direction u' (perp: the
/// part square to u), and n0 the same of the previous direction: a point
/// w of the section (square to u) meets the start mitre at w . n0 along u,
/// the end one at L + w . n1. Turns the same way add (a U-turn needs L >=
/// 2 r), a jog's opposite turns cancel, and a path's ends are cut square
/// (no term); a circle of radius rho sampled by n segments gives rho
/// cos(pi / n). Pairs are segments two or more apart; the closest points
/// of two segments are doubly critical when each is inside its segment, or
/// a vertex the path does not approach the other point from beyond (a
/// path's ends always are).
pub fn path_thickness(pts: &[V3], cap: f64) -> (Option<Thin>, Option<Thin>, usize) {
    let at = |a: V3, d: V3, s: f64| [a[0] + d[0] * s, a[1] + d[1] * s, a[2] + d[2] * s];
    // Repeated points (an arc sampled from the tangent point a run ends at)
    // make no segment.
    let scale = pts.iter().flatten().fold(1.0f64, |m, c| m.max(c.abs()));
    let eps = 1e-12 * scale;
    let mut p: Vec<V3> = Vec::with_capacity(pts.len());
    for &q in pts {
        if p.last().is_none_or(|&l| norm(sub(q, l)) > eps) {
            p.push(q);
        }
    }
    // A loop (three segments at least): segment i runs from p[i] to p[i + 1],
    // the last back to p[0]. Its ends coincide by std's rule (path-closes),
    // which the sweeps checked here build by.
    let closed = p.len() >= 4 && norm(sub(p[0], p[p.len() - 1])) < LOOP_GAP;
    if closed {
        p.pop();
    }
    let n = if closed { p.len() } else { p.len().saturating_sub(1) };
    let end = |i: usize| p[(i + 1) % p.len()];
    let d: Vec<V3> = (0..n).map(|i| sub(end(i), p[i])).collect();
    let len: Vec<f64> = d.iter().map(|&v| norm(v)).collect();
    let u: Vec<V3> = (0..n).map(|i| at([0.0; 3], d[i], 1.0 / len[i])).collect();
    let station: Vec<f64> = std::iter::once(0.0)
        .chain(len.iter().scan(0.0, |s, l| {
            *s += l;
            Some(*s)
        }))
        .collect();
    let prev = |i: usize| {
        if i > 0 {
            Some(i - 1)
        } else if closed {
            Some(n - 1)
        } else {
            None
        }
    };
    let next = |i: usize| {
        if i + 1 < n {
            Some(i + 1)
        } else if closed {
            Some(0)
        } else {
            None
        }
    };
    // The mitre between segment i and its neighbour j, as the shift along
    // u_i per unit of the section's offset: perp(u_j) / (1 + u_i . u_j).
    let mitre = |i: usize, j: Option<usize>| -> V3 {
        let Some(j) = j else { return [0.0; 3] };
        let c = dot(u[i], u[j]);
        let w = sub(u[j], at([0.0; 3], u[i], c));
        let den = 1.0 + c;
        if den <= 1e-15 { [f64::INFINITY; 3] } else { at([0.0; 3], w, 1.0 / den) }
    };
    let mut bend: Option<Thin> = None;
    for i in 0..n {
        let m = sub(mitre(i, next(i)), mitre(i, prev(i)));
        let m = norm(m);
        let r = if m > 0.0 { len[i] / m } else { f64::INFINITY };
        if r < cap && bend.as_ref().is_none_or(|b| r < b.radius) {
            bend = Some(Thin { radius: r, s: vec![station[i] + 0.5 * len[i]], at: vec![at(p[i], d[i], 0.5)] });
        }
    }
    // Whether point q, at parameter s of segment i, is critical for the
    // distance from the point o: inside the segment it is (the closest
    // points of two segments are a minimum over both); at a joint, moving
    // into the neighbouring segment must not bring it closer.
    let critical = |i: usize, s: f64, q: V3, o: V3| -> bool {
        let w = sub(q, o);
        let slack = -1e-9 * norm(w) - eps;
        match (s <= 0.0, s >= 1.0) {
            (true, _) => prev(i).is_none_or(|j| dot(w, at([0.0; 3], u[j], -1.0)) >= slack),
            (_, true) => next(i).is_none_or(|j| dot(w, u[j]) >= slack),
            _ => true,
        }
    };
    let lo: Vec<V3> = (0..n).map(|i| [0, 1, 2].map(|k| p[i][k].min(end(i)[k]))).collect();
    let hi: Vec<V3> = (0..n).map(|i| [0, 1, 2].map(|k| p[i][k].max(end(i)[k]))).collect();
    let mut pair: Option<Thin> = None;
    let mut close = 0;
    for i in 0..n {
        for j in i + 2..n {
            // Neighbours across a loop's closing point are not a pair.
            if closed && i == 0 && j == n - 1 {
                continue;
            }
            // Boxes at least 2 cap apart hold no pair closer.
            let gap2: f64 = (0..3).map(|k| (lo[j][k] - hi[i][k]).max(lo[i][k] - hi[j][k]).max(0.0).powi(2)).sum();
            if gap2 >= 4.0 * cap * cap {
                continue;
            }
            let (s, t) = closest(p[i], d[i], p[j], d[j]);
            let (a, b) = (at(p[i], d[i], s), at(p[j], d[j], t));
            let r = 0.5 * norm(sub(a, b));
            if r < cap && critical(i, s, a, b) && critical(j, t, b, a) {
                close += 1;
                if pair.as_ref().is_none_or(|q| r < q.radius) {
                    pair = Some(Thin { radius: r, s: vec![station[i] + s * len[i], station[j] + t * len[j]], at: vec![a, b] });
                }
            }
        }
    }
    (bend, pair, close)
}

/// Parameters (s, t) in [0, 1] of the closest points of segments a + s da
/// and b + t db. Parallel segments that overlap along their direction take
/// the middle of the overlap, inside both, so a joint is not chosen where
/// the inside of a segment is as close.
fn closest(a: V3, da: V3, b: V3, db: V3) -> (f64, f64) {
    let r = sub(a, b);
    let (aa, ee, bb) = (dot(da, da), dot(db, db), dot(da, db));
    let (c, f) = (dot(da, r), dot(db, r));
    let den = aa * ee - bb * bb;
    let t_of = |s: f64| ((bb * s + f) / ee).clamp(0.0, 1.0);
    if den <= 1e-12 * aa * ee {
        // b's ends projected onto a, as parameters of a.
        let (p0, p1) = (-c / aa, (bb - c) / aa);
        let (lo, hi) = (p0.min(p1).max(0.0), p0.max(p1).min(1.0));
        let s = if lo <= hi {
            0.5 * (lo + hi)
        } else if p0.max(p1) < 0.0 {
            0.0
        } else {
            1.0
        };
        let t = t_of(s);
        return (((bb * t - c) / aa).clamp(0.0, 1.0), t);
    }
    let s = ((bb * f - c * ee) / den).clamp(0.0, 1.0);
    let t = (bb * s + f) / ee;
    if t < 0.0 {
        ((-c / aa).clamp(0.0, 1.0), 0.0)
    } else if t > 1.0 {
        (((bb - c) / aa).clamp(0.0, 1.0), 1.0)
    } else {
        (s, t)
    }
}

impl Doc {
    /// `(check TEST VALUE)` (see the module).
    pub(super) fn eval_check(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let usage = "usage: (check TEST VALUE): VALUE, with a condition ink check reports when it fails; TEST is (value E CMP V), E a number \
                     of params and t, or (path-thickness PTS >= R), the thickness of the path PTS, a list of (vec x y z) or (vec x y)";
        let [test, value] = args else { return Err(self.err_hint(e, "bad check", Some(usage.into()))) };
        let Some([head, arg, op, v]) = test.list() else { return Err(self.err_hint(test, "bad check test", Some(usage.into()))) };
        let cmp = self.cmp_of(op)?;
        let v = self.field(v, env)?;
        let text = test.to_string();
        let test = match head.sym() {
            Some("value") => Test::Value { e: self.field(arg, env)?, cmp, v },
            Some("path-thickness") => {
                if cmp != Cmp::Ge {
                    return Err(self.err_hint(op, "path-thickness is tested with >=: the path must be at least as thick as R", Some(usage.into())));
                }
                let Value::Vec(items) = self.eval(arg, env)? else { return Err(self.err_hint(arg, "path-thickness takes a list of points", Some(usage.into()))) };
                let mut pts = Vec::with_capacity(items.len());
                for it in &items {
                    let c: Option<Vec<NodeId>> = match it {
                        Value::Vec(c) if (2..=3).contains(&c.len()) => c.iter().map(|x| if let Value::Field(f, _) = x { Some(*f) } else { None }).collect(),
                        _ => None,
                    };
                    let Some(c) = c else { return Err(self.err_hint(arg, "a path's point is a (vec x y z) or (vec x y) of numbers", Some(usage.into()))) };
                    let z = if c.len() == 3 { c[2] } else { self.g.constant(0.0) };
                    pts.push([c[0], c[1], z]);
                }
                Test::PathThickness { pts, v }
            }
            _ => return Err(self.err_hint(test, "unknown check test", Some(usage.into()))),
        };
        let nodes = test.nodes();
        if nodes.iter().any(|&n| self.g.depends_on(n, |v| matches!(v, Var::X | Var::Y | Var::Z))) {
            return Err(self.err_hint(e, "a check's test takes numbers (constants, params, t), not fields of x y z", Some(usage.into())));
        }
        let value = self.eval(value, env)?;
        // The same test stated again (a part used twice) is one condition,
        // reported where it was first stated.
        if self.checks_seen.insert((test.kind(), nodes)) {
            let calls = self.calls.iter().filter(|c| c.0 != "lambda").map(|(f, file, line)| format!("{f} ({file}:{line})")).collect();
            self.checks.push(Check { test, text, site: self.form_at.clone(), calls });
        }
        Ok(value)
    }

    /// The checks that fail at `ctx`, each as a warning saying where it was
    /// stated, its values and, for a path, where it is thin. A path's
    /// stretches may touch: they fail when they overlap by more than `tol`
    /// (the model's accuracy), so the test is against R - tol / 2. A path
    /// tested against several reaches (a hose's cover and bore) is reported
    /// for the largest, which fails wherever a smaller one does. Checks
    /// resting on a solve that failed at ctx (`failed`) are skipped: what
    /// rests on it is refused already.
    pub fn check_failures(&self, ctx: &Ctx, tol: f64, failed: &[usize]) -> Vec<String> {
        let num = |n: NodeId| self.g.as_const(n).unwrap_or_else(|| self.scalar(n, ctx));
        let r3 = |v: V3| format!("({})", v.iter().map(|c| format!("{}", (c * 1e3).round() / 1e3)).collect::<Vec<_>>().join(", "));
        let r = |c: f64| (c * 1e3).round() / 1e3;
        let reported: Vec<&Check> = self.checks.iter().filter(|x| !self.solves_read(&x.test.nodes()).iter().any(|i| failed.contains(i))).collect();
        // The largest reach each path is tested against.
        let mut largest: HashMap<&[[NodeId; 3]], f64> = HashMap::new();
        for x in &reported {
            if let Test::PathThickness { pts, v } = &x.test {
                let m = largest.entry(pts).or_insert(f64::NEG_INFINITY);
                *m = m.max(num(*v));
            }
        }
        let mut tested = HashSet::new();
        let mut out = Vec::new();
        for x in reported {
            let nodes = x.test.nodes();
            let calls = if x.calls.is_empty() { String::new() } else { format!(", in {}", x.calls.join(" > ")) };
            let mut values: Vec<String> = self.params_read(&nodes).into_iter().map(|k| format!("{} = {}", self.params[k].name, ctx.params[k])).collect();
            if nodes.iter().any(|&n| self.g.depends_on(n, |v| matches!(v, Var::T))) {
                values.insert(0, format!("t = {}", ctx.t));
            }
            let values = if values.is_empty() { String::new() } else { format!(" (at {})", values.join(", ")) };
            let head = format!("{}{calls}{values}", x.site);
            match &x.test {
                Test::Value { e, cmp, v } => {
                    let (e, v) = (num(*e), num(*v));
                    if !cmp.holds(e, v) {
                        out.push(format!("{head}: check {} fails: {e} {} {v} is false", x.text, cmp_text(*cmp)));
                    }
                }
                Test::PathThickness { pts, v } => {
                    let reach = num(*v);
                    if largest.get(pts.as_slice()).is_some_and(|&m| m > reach) || !tested.insert(pts.as_slice()) {
                        continue;
                    }
                    let pts: Vec<V3> = pts.iter().map(|p| p.map(num)).collect();
                    let (bend, pair, close) = path_thickness(&pts, reach - 0.5 * tol);
                    let swept = format!("what is swept along it reaches {} mm from the path", r(reach));
                    if let Some(b) = bend {
                        out.push(format!(
                            "{head}: the path turns too sharply for what is swept along it: at {} ({} mm along it) the mitres at the ends of its segment \
                             cross {} mm from the path, and {swept}: a section that reaches that far across the turn folds over itself there. Make the \
                             turn gentler (a bend radius of at least {} mm, or longer segments between turns) or the section smaller. A flat strip \
                             bent across its thickness reaches less far across the turn and may be clear; in a plane, build it as a band extruded \
                             square to the plane, checked on its half-thickness (docs/IR.md, Paths)",
                            r3(b.at[0]),
                            r(b.s[0]),
                            r(b.radius),
                            r(reach)
                        ));
                    }
                    if let Some(p) = pair {
                        let others = if close > 1 { format!(" ({close} pairs of segments in all)") } else { String::new() };
                        out.push(format!(
                            "{head}: the path runs through itself: its stretches at {} and {} ({} and {} mm along it) are {} mm apart, and {swept}, \
                             so they must be {} mm apart: the body overlaps itself there{others}. Reroute the path to keep them apart, or make the section smaller",
                            r3(p.at[0]),
                            r3(p.at[1]),
                            r(p.s[0]),
                            r(p.s[1]),
                            r(2.0 * p.radius),
                            r(2.0 * reach)
                        ));
                    }
                }
            }
        }
        out
    }
}

fn cmp_text(c: Cmp) -> &'static str {
    match c {
        Cmp::Lt => "<",
        Cmp::Le => "<=",
        Cmp::Gt => ">",
        Cmp::Ge => ">=",
    }
}
