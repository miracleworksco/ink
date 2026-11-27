//! std's `(carried PROFILE P LENGTH SHIFT)` (std.ink's `carried` and
//! `path-of`), as the packages that bend or move a carried body read it:
//! sharing's GLB strips and dynamics' moving material. To the language it
//! is an ordinary attribute; its layout is read here and nowhere else.

use ink_core::graph::{BinOp, Graph, NodeId, UnOp};
use ink_lang::Value;

/// A piece of std's path (std.ink, paths of lines and arcs), its numbers as
/// fields (they may read t).
#[derive(Clone, Debug, PartialEq)]
pub enum PathPiece {
    Line {
        a: [NodeId; 3],
        b: [NodeId; 3],
    },
    /// About centre c and unit axis u, radius r, from c + r f, sweep radians,
    /// rho mm of station per radian (r unless given).
    Arc {
        c: [NodeId; 3],
        u: [NodeId; 3],
        r: NodeId,
        f: [NodeId; 3],
        sweep: NodeId,
        rho: NodeId,
    },
}

impl PathPiece {
    /// Its numbers, in order.
    pub fn numbers(&self) -> Vec<NodeId> {
        match self {
            PathPiece::Line { a, b } => [&a[..], &b[..]].concat(),
            PathPiece::Arc { c, u, r, f, sweep, rho } => [&c[..], &u[..], &[*r], &f[..], &[*sweep, *rho]].concat(),
        }
    }
}

/// std's `(carried PROFILE P LENGTH SHIFT)`: a body that is PROFILE (a
/// solid in material coordinates, x the station along the material from 0
/// to LENGTH) carried along path P, material station m at path station
/// m + SHIFT.
#[derive(Clone, Debug, PartialEq)]
pub struct Carried {
    pub profile: NodeId,
    pub pieces: Vec<PathPiece>,
    /// The path's `up`, unit.
    pub up: [NodeId; 3],
    pub closed: bool,
    pub length: NodeId,
    pub shift: NodeId,
    /// Each piece's station at its start, and its length; the path's.
    pub starts: Vec<NodeId>,
    pub lengths: Vec<NodeId>,
    pub total: NodeId,
}

impl Carried {
    /// The `carried` among a body's attributes: None if it has none, an
    /// error if it is not std's.
    pub fn in_attrs(g: &Graph, attrs: &[(String, Value)]) -> Result<Option<Carried>, String> {
        let Some((_, v)) = attrs.iter().find(|(k, _)| k == "carried") else { return Ok(None) };
        let c = Carried::of(g, v);
        c.map(Some).ok_or_else(|| "attr carried is not std's (carried PROFILE P LENGTH SHIFT): PROFILE a field, P a (path ...) or (path-loop ...), LENGTH and SHIFT numbers".into())
    }

    /// The description in a `carried` attribute's value; None if it is not
    /// std's.
    fn of(g: &Graph, v: &Value) -> Option<Carried> {
        let point = |v: &Value| match v {
            Value::Vec(c) => match c.as_slice() {
                [Value::Field(a, _), Value::Field(b, _), Value::Field(c, _)] => Some([*a, *b, *c]),
                _ => None,
            },
            _ => None,
        };
        let field = |v: &Value| if let Value::Field(f, _) = v { Some(*f) } else { None };
        let Value::Vec(v) = v else { return None };
        let [Value::Field(profile, _), Value::Vec(path), Value::Field(length, _), Value::Field(shift, _)] = v.as_slice() else { return None };
        // path-of's (vec pieces up closed starts total)
        let [Value::Vec(pieces), up, Value::Field(closed, _), Value::Vec(starts), Value::Field(total, _)] = path.as_slice() else { return None };
        let (pieces, lengths): (Vec<PathPiece>, Vec<NodeId>) = pieces
            .iter()
            .map(|p| {
                let Value::Vec(p) = p else { return None };
                // path-line's (vec 0 length a b), path-arc's (vec 1 length c u r f sweep rho)
                let piece = match (p.first().and_then(field).and_then(|k| g.as_const(k)), p.as_slice()) {
                    (Some(0.0), [_, _, a, b]) => PathPiece::Line { a: point(a)?, b: point(b)? },
                    (Some(1.0), [_, _, c, u, r, f, sweep, rho]) => {
                        PathPiece::Arc { c: point(c)?, u: point(u)?, r: field(r)?, f: point(f)?, sweep: field(sweep)?, rho: field(rho)? }
                    }
                    _ => return None,
                };
                Some((piece, field(&p[1])?))
            })
            .collect::<Option<Vec<_>>>()?
            .into_iter()
            .unzip();
        let starts = starts.iter().map(field).collect::<Option<Vec<_>>>()?;
        if pieces.is_empty() || starts.len() != pieces.len() {
            return None;
        }
        Some(Carried { profile: *profile, pieces, up: point(up)?, closed: g.as_const(*closed)? == 1.0, length: *length, shift: *shift, starts, lengths, total: *total })
    }

    /// The angle arc i has turned `along` mm from its start (signed as its
    /// sweep: std's path-arc-angle).
    fn angle(&self, g: &mut Graph, i: usize, along: NodeId) -> NodeId {
        let PathPiece::Arc { sweep, .. } = self.pieces[i] else { return g.constant(0.0) };
        let a = g.mul(along, sweep);
        g.bin(BinOp::Div, a, self.lengths[i])
    }

    /// The point `along` mm from piece i's start, its line or circle run
    /// on past its ends (std's path-piece-move).
    pub fn point(&self, g: &mut Graph, i: usize, along: NodeId) -> [NodeId; 3] {
        match self.pieces[i] {
            PathPiece::Line { a, b } => {
                let d = sub(g, b, a).map(|c| g.bin(BinOp::Div, c, self.lengths[i]));
                let run = scale(g, d, along);
                add(g, a, run)
            }
            PathPiece::Arc { c, u, r, f, .. } => {
                let th = self.angle(g, i, along);
                let (co, si) = (g.un(UnOp::Cos, th), g.un(UnOp::Sin, th));
                let w = cross(g, u, f);
                let (rc, rs) = (g.mul(r, co), g.mul(r, si));
                let (fc, ws) = (scale(g, f, rc), scale(g, w, rs));
                let off = add(g, fc, ws);
                add(g, c, off)
            }
        }
    }
}

fn cross(g: &mut Graph, a: [NodeId; 3], b: [NodeId; 3]) -> [NodeId; 3] {
    std::array::from_fn(|k| {
        let (i, j) = ((k + 1) % 3, (k + 2) % 3);
        let (p, q) = (g.mul(a[i], b[j]), g.mul(a[j], b[i]));
        g.sub(p, q)
    })
}

fn add(g: &mut Graph, a: [NodeId; 3], b: [NodeId; 3]) -> [NodeId; 3] {
    std::array::from_fn(|k| g.add(a[k], b[k]))
}

fn sub(g: &mut Graph, a: [NodeId; 3], b: [NodeId; 3]) -> [NodeId; 3] {
    std::array::from_fn(|k| g.sub(a[k], b[k]))
}

fn scale(g: &mut Graph, a: [NodeId; 3], k: NodeId) -> [NodeId; 3] {
    a.map(|c| g.mul(c, k))
}

/// std's carried path is read by its pieces' names, and a value that is
/// not std's is refused, not taken as no path without a word.
#[test]
fn reads_stds_carried_path() {
    use ink_lang::{Doc, Registry};
    static NONE: Registry = Registry::new(&[]);
    let load = |body: &str| {
        let src = format!(
            concat!(
                "(import \"std\")\n",
                "(def P (path-loop (vec (path-line (vec 0 -10 0) (vec 50 -10 0)) (path-arc (vec 50 0 0) (vec 0 0 1) 10 (vec 0 -1 0) pi)\n",
                "  (path-line (vec 50 10 0) (vec 0 10 0)) (path-arc (vec 0 0 0) (vec 0 0 1) 10 (vec 0 1 0) pi)) (vec 0 0 1)))\n{}"
            ),
            body
        );
        let d = Doc::load_str_in(&NONE, &src, "t", std::path::Path::new("."), &[]).unwrap_or_else(|e| panic!("{}", e.message));
        let c = Carried::in_attrs(&d.g, &d.body("b").unwrap().attrs);
        (d, c)
    };
    let (d, c) = load("(body b (shape (sphere 1)) (attr carried (carried (box 1 1 1) P (path-length P) (* 10 (/ t 1s)))))");
    let c = c.unwrap().unwrap();
    assert!(c.closed && c.pieces.len() == 4 && d.g.as_const(c.length).is_some_and(|l| (l - (100.0 + 20.0 * std::f64::consts::PI)).abs() < 1e-9));
    let PathPiece::Arc { r, rho, c: centre, .. } = &c.pieces[1] else { panic!("{:?}", c.pieces[1]) };
    assert!(r == rho && d.g.as_const(*r) == Some(10.0) && d.g.as_const(centre[0]) == Some(50.0));
    assert!(matches!(c.pieces[0], PathPiece::Line { .. }));
    assert!(load("(body b (shape (sphere 1)) (attr carried (vec (box 1 1 1) 3)))").1.unwrap_err().contains("is not std's (carried PROFILE P LENGTH SHIFT)"));
    assert_eq!(load("(body b (shape (sphere 1)))").1, Ok(None));
}
