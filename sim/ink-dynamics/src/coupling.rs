//! `(couple LINES T0 [E [EB]])` in a relation: the clause this package owns
//! there (Package::relation_clauses), kept with the relation's kind and
//! placed with each relation under a target (`related`).

use ink_core::graph::{BinOp, NodeId, Var};
use ink_lang::sexp::Sexp;
use ink_lang::{Doc, InkError, Kind, Quantity, Reader, Value, Xf};
use std::any::Any;
use std::rc::Rc;
use std::sync::Arc;

const USAGE: &str = "(couple LINES T0 [E [EB]]): what A and B exert on each other along LINES, each (vec P U [H [L [K]]]): B pulled along U, A \
                     the other way, through point P, with a moment H (mm, a screw's lead / 2 pi) times the tension about U, L its length (mm, a \
                     belt span's), K its stiffness (N per mm its ends draw apart along it; not given: rigid); T0 the pretension (N) a tension the parts' balance leaves free is installed at, locked (the belt keeps its \
                     length; a span that would push slack), or (held T0), a constant-force tensioner holding the least line at T0; E the share \
                     of the power A gives that B receives, EB of what B gives that A receives (default E), not given: unknown; among \
                     them std's (tensile-stiffness EA), EA a force (16903N) the tension that would stretch a line by its length, makes each line a \
                     spring of EA / L";

/// `(couple LINES T0 [E [EB]])` in a relation: what A and B exert on each
/// other. Each line pulls B along its direction and A the other way, its
/// tension T >= 0 (a belt's span, a gear's line of action, a screw's axis
/// either way), with a moment of its pitch times T about it (a screw's lead
/// / 2 pi). A tension the parts' balance leaves free (a belt's spans
/// tightened together) is set by T0, installed and locked: the belt keeps
/// its length, its spans' tensions less T0, weighted by their lengths,
/// summing to 0, and a span that would push slack at 0; `(held T0)`, a
/// constant-force tensioner, holds its least line at T0. E is the share of
/// the power A gives that B receives, EB of what B gives that A receives
/// (None: not given, so unknown).
#[derive(Clone, Debug)]
pub struct Couple {
    pub lines: Vec<CoupleLine>,
    pub pretension: NodeId,
    pub held: bool,
    pub efficiency: Option<(NodeId, NodeId)>,
    /// Where its relation is related (file:line): what states its lines'
    /// stiffness.
    pub site: String,
}

/// A line of a coupling: a point on it and its direction (the relating
/// assembly's coordinates as read, the world's where placed), its pitch
/// (mm: the moment about it per N of its tension), its length (mm: a belt
/// span's, weighing its stretch; None: alike), and its stiffness (N per mm
/// its ends draw apart along it: a belt span's EA / L; None: rigid).
#[derive(Clone, Copy, Debug)]
pub struct CoupleLine {
    pub at: [NodeId; 3],
    pub dir: [NodeId; 3],
    pub pitch: NodeId,
    pub length: Option<NodeId>,
    pub stiffness: Option<NodeId>,
}

impl Couple {
    pub fn fields(&self) -> Vec<NodeId> {
        let lines = self.lines.iter().flat_map(|l| l.at.into_iter().chain(l.dir).chain([l.pitch]).chain(l.length).chain(l.stiffness));
        lines.chain([self.pretension]).chain(self.efficiency.iter().flat_map(|e| [e.0, e.1])).collect()
    }

    /// A kind's coupling, as read (the relating assembly's coordinates).
    pub fn of(kind: &Kind) -> Option<&Couple> {
        kind.clauses.iter().find(|c| c.package == "dynamics" && c.head == "couple").and_then(|c| c.kept.downcast_ref::<Couple>())
    }
}

/// A relation under a target (Doc::relations_under): its ends as full paths,
/// its kind, and its coupling's lines in the world, where it has one.
pub struct Related {
    pub a: String,
    pub b: String,
    pub kind: Option<Rc<Kind>>,
    pub couple: Option<Couple>,
}

/// The relations under `target`, each coupling's lines placed in the world:
/// a point by the relating assembly's placement, a direction turned.
pub fn related(doc: &mut Doc, target: &str) -> Vec<Related> {
    let mut out = Vec::new();
    for r in doc.relations_under(target) {
        let couple = r.kind.as_deref().and_then(Couple::of).cloned().map(|mut c| {
            for l in &mut c.lines {
                let at = turned(doc, &r.at, l.at);
                l.at = std::array::from_fn(|i| doc.g.add(at[i], r.at.t[i]));
                l.dir = turned(doc, &r.at, l.dir);
            }
            c.site = r.kind.as_ref().map_or(String::new(), |k| k.site.clone());
            c
        });
        out.push(Related { a: r.a, b: r.b, kind: r.kind, couple });
    }
    out
}

fn turned(doc: &mut Doc, xf: &Xf, v: [NodeId; 3]) -> [NodeId; 3] {
    std::array::from_fn(|i| {
        let a = doc.g.mul(xf.r[i][0], v[0]);
        let b = doc.g.mul(xf.r[i][1], v[1]);
        let c = doc.g.mul(xf.r[i][2], v[2]);
        let ab = doc.g.add(a, b);
        doc.g.add(ab, c)
    })
}

/// Read `(couple LINES T0 [E [EB]])` in a relation (Package::clause): LINES
/// a vec of lines, each (vec P U [H [L [K]]]), P a point on it and U its
/// direction (vecs of 3, in the relating assembly's coordinates), H its
/// pitch (mm, default 0); T0 the pretension (N); E, EB the efficiencies
/// forward and back (EB defaulting to E). Numbers of params, t and frames.
pub fn read(r: &mut Reader, c: &Sexp, args: &[Sexp]) -> Result<Arc<dyn Any + Send + Sync>, InkError> {
    let [lines, t0, eff @ ..] = args else { return Err(r.err_hint(c, "bad couple", Some(USAGE.into()))) };
    if eff.len() > 2 {
        return Err(r.err_hint(c, "bad couple", Some(USAGE.into())));
    }
    let number = |r: &mut Reader, v: &Value, at: &Sexp| -> Result<NodeId, InkError> {
        let Value::Field(f, _) = v else { return Err(r.err_hint(at, format!("a coupling's numbers are numbers, not a {}", v.kind()), Some(USAGE.into()))) };
        if r.graph().depends_on(*f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
            return Err(r.err(at, "a coupling's numbers are of constants, params, t and frames, not fields of x y z"));
        }
        Ok(*f)
    };
    // EA is a force (N), said in its unit (16903N, 10kN, 3.8lbf): a bare number
    // is refused, as a joint's (stiffness K) refuses one (Reader::field_of).
    let force = |r: &mut Reader, v: &Value, at: &Sexp| -> Result<NodeId, InkError> {
        const WHAT: &str = "a tensile stiffness (EA)";
        let f = number(r, v, at)?;
        let n = match v {
            Value::Field(_, Some(said)) => match said.read_as(Quantity::Force, WHAT) {
                Ok(1.0) => f,
                Ok(k) => {
                    let k = r.graph().constant(k);
                    r.graph().mul(f, k)
                }
                Err(why) => return Err(r.err_hint(at, why, Some(USAGE.into()))),
            },
            _ if r.graph().as_const(f) == Some(0.0) => f,
            _ => {
                return Err(r.err_hint(
                    at,
                    format!("{WHAT} takes a force, and it is a bare number: say its unit ((tensile-stiffness 16903N), or of a value worked out (* 1N EA))"),
                    Some(format!("a bare number is a length in mm or a pure number; a force says its unit: {}", ink_lang::units::units_of(Quantity::Force))),
                ));
            }
        };
        // Known at load, it is at least 0 and finite (a mode and a
        // simulation check one of params where they read it).
        match r.graph().as_const(n) {
            Some(k) if !(k >= 0.0 && k.is_finite()) => Err(r.err_hint(at, format!("{WHAT} is {k} N: a stiffness is at least 0 and finite"), Some(USAGE.into()))),
            _ => Ok(n),
        }
    };
    let three = |r: &mut Reader, v: &Value, at: &Sexp| -> Result<[NodeId; 3], InkError> {
        match v {
            Value::Vec(xs) if xs.len() == 3 => Ok([number(r, &xs[0], at)?, number(r, &xs[1], at)?, number(r, &xs[2], at)?]),
            _ => Err(r.err_hint(at, "a coupling line's point and direction are vecs of 3", Some(USAGE.into()))),
        }
    };
    let each = "each coupling line is (vec P U [H [L [K]]])";
    let Value::Vec(written) = r.value(lines)? else { return Err(r.err_hint(lines, "LINES is a vec of lines", Some(USAGE.into()))) };
    let mut out = Vec::new();
    for line in written.iter() {
        let Value::Vec(parts) = line else { return Err(r.err_hint(lines, each, Some(USAGE.into()))) };
        let (p, u, rest) = match parts.as_slice() {
            [p, u, rest @ ..] if rest.len() <= 3 => (p, u, rest),
            _ => return Err(r.err_hint(lines, each, Some(USAGE.into()))),
        };
        let (at, dir) = (three(r, p, lines)?, three(r, u, lines)?);
        let pitch = match rest.first() {
            Some(h) => number(r, h, lines)?,
            None => r.graph().constant(0.0),
        };
        let length = rest.get(1).map(|l| number(r, l, lines)).transpose()?;
        let stiffness = rest.get(2).map(|k| number(r, k, lines)).transpose()?;
        out.push(CoupleLine { at, dir, pitch, length, stiffness });
    }
    if out.is_empty() {
        return Err(r.err_hint(lines, "a coupling has a line", Some(USAGE.into())));
    }
    // (held T0) is std's (vec "held" T0).
    let (pretension, held) = match r.value_of(t0, Quantity::Force, "a coupling's pretension")? {
        Value::Vec(v) if matches!(v.as_slice(), [Value::Str(h), _] if h == "held") => (number(r, &v[1], t0)?, true),
        v => (number(r, &v, t0)?, false),
    };
    // E and EB, or a vec of them (a kind passing on its optional
    // arguments: none, E, or E and EB); and among them std's
    // (tensile-stiffness EA), (vec "tensile-stiffness" EA), which no
    // number can be read as.
    let tagged = |v: &Value| matches!(v, Value::Vec(t) if matches!(t.as_slice(), [Value::Str(h), _] if h == "tensile-stiffness"));
    let mut given = Vec::new();
    let mut ea = None;
    for x in eff {
        let vs = match r.value(x)? {
            v if tagged(&v) => vec![v],
            Value::Vec(vs) => std::rc::Rc::unwrap_or_clone(vs),
            v => vec![v],
        };
        for v in vs {
            match v {
                Value::Vec(t) if tagged(&Value::Vec(t.clone())) => {
                    if ea.is_some() {
                        return Err(r.err_hint(x, "a coupling's tensile stiffness is given twice", Some(USAGE.into())));
                    }
                    ea = Some(force(r, &t[1], x)?);
                }
                v => given.push((v, x)),
            }
        }
    }
    if let Some(ea) = ea {
        for l in &mut out {
            let Some(len) = l.length else {
                return Err(r.err_hint(lines, "a tensile stiffness makes each line a spring of EA / L: give each line its length, (vec P U H L)", Some(USAGE.into())));
            };
            if l.stiffness.is_some() {
                return Err(r.err_hint(lines, "a line states its stiffness K and the coupling its tensile stiffness: state one", Some(USAGE.into())));
            }
            l.stiffness = Some(r.graph().bin(BinOp::Div, ea, len));
        }
    }
    let efficiency = match given.as_slice() {
        [] => None,
        [(e, at)] => {
            let e = number(r, e, at)?;
            Some((e, e))
        }
        [(e, at), (b, bt)] => Some((number(r, e, at)?, number(r, b, bt)?)),
        _ => return Err(r.err_hint(&eff[0], "a coupling's efficiency is E, or E and EB", Some(USAGE.into()))),
    };
    Ok(Arc::new(Couple { lines: out, pretension, held, efficiency, site: String::new() }))
}
