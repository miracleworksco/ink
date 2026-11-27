//! `(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)])` in
//! an assembly: the clause this package owns (Package::clauses), kept with
//! the assembly and placed with each instance of it.

use ink_core::graph::{Graph, NodeId, Var};
use ink_lang::InkError;
use ink_lang::sexp::Sexp;
use ink_lang::{Doc, PlacedClause, Quantity, Reader, Xf};
use std::any::Any;
use std::sync::Arc;

const USAGE: &str = "(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)]) in an assembly: a force in N and a torque in N mm, in \
                     this assembly's axes, on what its instance PATH places (an instance read so far, or a path into one), acting at FRAME's \
                     origin (a frame in its coordinates, (frame-of rod tip); default PATH's own origin); components are fields of its \
                     arguments, params and t, and may read its instances' frames (a spring: along the line between two of them). Placed with \
                     each instance of the assembly, it is a load of every mechanism placing it (docs/DYNAMICS.md)";

/// An assembly's load: a force (N) and torque (N mm) in the assembly's axes
/// on what its instance path `path` places, acting at `at`'s origin (a
/// frame in its coordinates; None: the path's own origin), fields of its
/// arguments, params and t, which may read its instances' frames.
pub struct AssemblyLoad {
    pub name: String,
    pub path: String,
    pub force: [NodeId; 3],
    pub torque: [NodeId; 3],
    pub at: Option<Xf>,
}

/// An assembly's load under a target (loads_under): its name and path as
/// full paths (matching `Item::label`), its force and torque in the world's
/// axes, and the frame it acts at in the world; where it is written.
pub struct PlacedLoad {
    pub name: String,
    pub path: String,
    pub force: [NodeId; 3],
    pub torque: [NodeId; 3],
    pub at: Option<Xf>,
    pub file: String,
    pub line: u32,
}

/// Read `(load ...)` in an assembly (Package::clause).
pub fn read(r: &mut Reader, c: &Sexp, args: &[Sexp]) -> Result<Arc<dyn Any + Send + Sync>, InkError> {
    let [name, path, rest @ ..] = args else { return Err(r.err_hint(c, "bad load", Some(USAGE.into()))) };
    let name = r.name(name)?;
    let p = r.instance_path(path).map_err(|e| InkError { hint: Some(USAGE.into()), ..e })?;
    if r.kept("load").iter().filter_map(|k| k.downcast_ref::<AssemblyLoad>()).any(|l| l.name == name) {
        return Err(r.err(c, format!("this assembly has two loads named '{name}'")));
    }
    let zero = r.graph().constant(0.0);
    let (mut force, mut torque, mut at) = (None, None, None);
    for cl in rest {
        let l = cl.list().unwrap_or(&[]);
        match (l.first().and_then(|h| h.sym()), l.get(1..).unwrap_or(&[])) {
            (Some(k @ ("force" | "torque")), [a, b, cc]) => {
                let q = if k == "force" { Quantity::Force } else { Quantity::Torque };
                let mut v = [zero; 3];
                for (i, e) in [a, b, cc].into_iter().enumerate() {
                    v[i] = r.field_of(e, q, k)?;
                    if r.graph().depends_on(v[i], |x| matches!(x, Var::X | Var::Y | Var::Z)) {
                        return Err(r.err_hint(e, "a load's components are fields of params and t, not of x y z", Some(USAGE.into())));
                    }
                }
                let slot = if k == "force" { &mut force } else { &mut torque };
                if slot.replace(v).is_some() {
                    return Err(r.err(cl, format!("one ({k} ...) per load")));
                }
            }
            (Some("at"), [f]) if at.is_none() => at = Some(r.frame(f, "a load's (at FRAME)")?),
            _ => return Err(r.err_hint(cl, "bad load clause", Some(USAGE.into()))),
        }
    }
    if force.is_none() && torque.is_none() {
        return Err(r.err_hint(c, "a load needs (force fx fy fz) or (torque tx ty tz)", Some(USAGE.into())));
    }
    Ok(Arc::new(AssemblyLoad { name, path: p, force: force.unwrap_or([zero; 3]), torque: torque.unwrap_or([zero; 3]), at }))
}

/// Every assembly load placed under `target` (a top-level assembly, or a
/// path into one), its path in full and its force, torque and frame in the
/// world: turned and moved as its assembly is placed.
pub fn loads_under(doc: &mut Doc, target: &str) -> Vec<PlacedLoad> {
    let placed = doc.clauses_under(target, "dynamics");
    place(placed, &mut |xf, l| (turned(&mut doc.g, xf, l.force), turned(&mut doc.g, xf, l.torque), l.at.map(|f| doc.compose(&f, xf))))
}

/// loads_under, read through a package's reader (what they read recorded).
pub fn loads_read(r: &mut Reader, target: &str) -> Result<Vec<PlacedLoad>, String> {
    let placed = r.clauses_under(target, "dynamics")?;
    Ok(place(placed, &mut |xf, l| (turned(r.graph(), xf, l.force), turned(r.graph(), xf, l.torque), l.at.map(|f| r.compose(&f, xf)))))
}

/// `v` in the world's axes, its assembly placed by `xf`.
fn turned(g: &mut Graph, xf: &Xf, v: [NodeId; 3]) -> [NodeId; 3] {
    std::array::from_fn(|i| {
        let a = g.mul(xf.r[i][0], v[0]);
        let b = g.mul(xf.r[i][1], v[1]);
        let c = g.mul(xf.r[i][2], v[2]);
        let ab = g.add(a, b);
        g.add(ab, c)
    })
}

type Placed = ([NodeId; 3], [NodeId; 3], Option<Xf>);

fn place(placed: Vec<PlacedClause>, at: &mut dyn FnMut(&Xf, &AssemblyLoad) -> Placed) -> Vec<PlacedLoad> {
    let mut out = Vec::new();
    for p in placed {
        let Some(l) = p.clause.kept.downcast_ref::<AssemblyLoad>() else { continue };
        let (force, torque, frame) = at(&p.at, l);
        let prefix = &p.path;
        out.push(PlacedLoad {
            name: format!("{prefix}/{}", l.name),
            path: format!("{prefix}/{}", l.path),
            force,
            torque,
            at: frame,
            file: p.clause.file.to_string(),
            line: p.clause.line,
        });
    }
    out
}
