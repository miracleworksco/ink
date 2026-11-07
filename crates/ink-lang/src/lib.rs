//! The language (docs/CORE.md, docs/IR.md): reads its text (sexp.rs; numbers
//! in their units, units.rs), expands S-expressions into the field graph
//! and holds the document's structure (params, bodies, assemblies, claims);
//! the registry and reader packages plug their forms
//! into (ext.rs, reader.rs); load-time programming (author.rs) and kernel
//! services at load time (services.rs); results as fields (results.rs); the
//! materials library (materials.rs); a loaded model and its analysis
//! (view.rs), and checking a document's claims on it (claims.rs). It names
//! no package: a document reads the packages of the registry it is loaded
//! with.

// The crates it builds on, at the paths the language's code uses.
use ink_core::{V3, eval, graph, tape};
mod error;
pub use error::InkError;
use ink_kernel::{cache, mesh, query};
pub use units::Quantity;

/// Sampled fields (ink-core) and meshes read into fields (ink-io).
mod grid {
    pub use ink_core::grid::*;
    pub use ink_io::stl::{MeshDistance, read_stl};
}

use crate::eval::{Ctx, eval_point};
use crate::graph::{BinOp, Graph, NodeId, UnOp, Var};
use crate::sexp::{SK, Sexp, parse};
use crate::tape::Tape;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::rc::Rc;

pub mod answer;
mod author;
/// Conditions builders state on what they build (`check`).
pub mod check;
/// Checking a document's claims on a model view.
pub mod claims;
pub mod command;
mod ext;
pub mod follow;
mod frame;
mod inexact;
pub mod interfere;
pub mod look;
pub mod materials;
pub mod motion;
pub mod pairs;
mod profile;
mod rates;
mod reader;
mod records;
mod results;
/// Certified roots of a solve's equations (Krawczyk's test).
pub mod roots;
mod services;
pub mod sexp;
pub mod units;
pub use units::Given;
pub mod validity;
/// A loaded model and what analysis has learned of it.
pub mod view;

pub use self::ext::*;
pub use self::follow::Follow;
pub use self::inexact::*;
pub use self::reader::*;
pub use self::results::sampled_result;

pub const STD_SRC: &str = include_str!("../../../std/std.ink");

const MAX_CALL_DEPTH: usize = 256;
/// The window, where a shape not proven bounded is examined (see
/// `Doc::domain`): the 1 m cube at the origin.
pub const DEFAULT_DOMAIN: (V3, V3) = ([-500.0; 3], [500.0; 3]);

/// The accuracy an automatic proof works to, as a fraction of the window
/// (see `Doc::finest`): 1/20 000, 0.05 mm, below what parts are made to,
/// for every model. Surfaces that touch by design (teeth in mesh, threads,
/// a screw seated on a plate) never resolve, so the search stops at cells
/// half of it and proves them to the accuracy (twice the cell: this); a gap
/// narrower than about this is proven to the accuracy too, not outright. A
/// claim that needs finer (MEMS, optics) says so, `(accuracy x)` or
/// `(res r)`, as does `--res`.
pub const ACCURACY: f64 = 1.0 / 20_000.0;

/// Where analysis of a field looks (see `Doc::reach`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reach {
    pub lo: V3,
    pub hi: V3,
    /// The box holds all of the field's material. If not (the field is not
    /// proven bounded) it is the window, and a result that needs all of the
    /// shape (emptiness, an upper bound on volume, a box, a distance's lower
    /// bound) is unknown; one from a witness inside still holds.
    pub whole: bool,
}

#[derive(Clone, Debug)]
pub enum Value {
    /// A number expression, and the number written with a unit its
    /// quantity was said by, if any (units::Said).
    Field(NodeId, Option<Rc<units::Said>>),
    /// A point (`probe`): three numbers, constants and params. Claims read
    /// it as a point (its sign in a shape, exactly); `place` puts it at a
    /// frame (R p + T, forward); it is no field, and no warp moves it.
    Point([NodeId; 3]),
    Str(String),
    Vec(Vec<Value>),
    Data(Data),
    /// A function as a value (for map, reduce, apply, and higher-order fns).
    Fn(Callable),
    /// A rigid placement: an origin and x, y, z axes (`(at ...)`).
    Frame(Xf),
    /// What a kind gives `(relate A B KIND)`: `(relation NAME ...)`.
    Relation(Rc<RelationValue>),
}

/// `(relation NAME [(zone WHAT ZA [ZB])] (claim NAME VALUE LO HI) ...
/// [(data KEY VALUE ...)] [CLAUSE ...])`: a relation of two parts as a kind
/// writes it, CLAUSE a package's (Package::relation_clauses).
#[derive(Debug)]
pub struct RelationValue {
    pub kind: String,
    pub zone: Option<(String, Zone)>,
    /// Each LO <= VALUE <= HI, a value claim of the assembly relating them,
    /// and the accuracy it states, `(accuracy x)` in VALUE's units (None: 0).
    pub claims: Vec<(String, NodeId, NodeId, NodeId, Option<f64>)>,
    pub data: Vec<(String, Value)>,
    /// Its clauses that packages own, in order.
    pub clauses: Vec<OwnedClause>,
}

impl RelationValue {
    /// The fields of its zone and claims, and those its packages' clauses
    /// read (its data's are values).
    pub fn fields(&self) -> Vec<NodeId> {
        let zone = self.zone.iter().flat_map(|(_, z)| std::iter::once(z.a).chain(z.b));
        let owned = self.clauses.iter().flat_map(|c| c.reads.fields.iter().copied());
        zone.chain(self.claims.iter().flat_map(|c| [c.1, c.2, c.3])).chain(owned).collect()
    }
}

#[derive(Clone, Debug)]
pub enum Callable {
    Builtin(String),
    User(String, Rc<FnDef>),
    Closure(Rc<Closure>),
}

/// A lambda with the local names it captured.
#[derive(Debug)]
pub struct Closure {
    pub params: Vec<String>,
    /// The last param takes the remaining arguments, as a vec (`rest ...`).
    pub rest: bool,
    pub body: Sexp,
    pub env: Vec<(String, Value)>,
    pub file: Rc<str>,
    pub scope: Scope,
}

/// Where a function reads the names it does not bind (lexical scope): the
/// file it is written in (an index into Doc::files; 0 is std). It sees the
/// names of that file, of the files it imports (SourceFile::bindings), and std's:
/// never one of the file calling it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scope {
    pub file: u32,
}

impl Scope {
    pub const STD: Scope = Scope { file: 0 };
}

/// A file run: its name, its canonical path and the namespace it is
/// imported under (one file under two is run twice), its directory (the
/// paths it reads are relative to it), and the names it reads, as written
/// in it: its own and those of the files it imports, transitively (an
/// import `as p` adds them as `p.name`). `bindings`: each value's and
/// function's key; `objects`: each body's and assembly's name in the
/// document (`p.plate` for `plate` in a file imported as p: Doc::resolve).
/// `ahead`: where it defines each value, and `declared` every name its
/// forms declare (top_names: bodies and assemblies too). `given`: the
/// values its import gave its params (ImportValue).
struct SourceFile {
    name: Rc<str>,
    path: String,
    ns: Option<String>,
    dir: std::path::PathBuf,
    bindings: HashMap<String, String>,
    objects: HashMap<String, String>,
    ahead: HashMap<String, u32>,
    declared: HashSet<String>,
    given: Vec<ImportValue>,
}

impl SourceFile {
    fn new(name: Rc<str>, path: String, ns: Option<String>, dir: std::path::PathBuf) -> SourceFile {
        SourceFile { name, path, ns, dir, bindings: HashMap::new(), objects: HashMap::new(), ahead: HashMap::new(), declared: HashSet::new(), given: Vec::new() }
    }
}

/// A value an import gives a param of its file, `(import "a.ink" (base-t
/// 3.2 [(tol lo hi)]))`: the param's default for the whole load, as if
/// written there, and its band where it gives one; where it is written,
/// and the units its value and band ends are said in, each where
/// (ImportValue::fits).
#[derive(Clone)]
struct ImportValue {
    name: String,
    value: f64,
    said: Option<Rc<units::Said>>,
    tol: Option<(f64, f64)>,
    file: Rc<str>,
    line: u32,
    units: Vec<(Rc<units::Said>, u32, u32)>,
}

impl ImportValue {
    /// Refused, at the import, where its value or band is said in a unit
    /// of another quantity than `q`, its param's (`what`).
    fn fits(&self, q: Quantity, what: &str) -> R<()> {
        match self.units.iter().find_map(|(s, line, col)| Some((s.misfit(q, what)?, *line, *col))) {
            Some((message, line, col)) => Err(InkError { message, file: self.file.to_string(), line, col, hint: None, called_from: Vec::new() }),
            None => Ok(()),
        }
    }
}

/// JSON data (`read-json`, `data`, a material): a place in a document
/// shared by every value read from it, not a copy. A file is read once
/// per load however many files read it (a family's table, read by each of
/// its parts), and `get` moves the place.
#[derive(Clone, Debug)]
pub struct Data {
    root: Rc<serde_json::Value>,
    /// Its JSON pointer in `root`.
    at: Rc<str>,
}

impl Data {
    pub fn new(v: serde_json::Value) -> Data {
        Data { root: Rc::new(v), at: Rc::from("") }
    }

    /// The place `key` names in this object.
    fn get(&self, key: &str) -> Option<Data> {
        self.value().get(key)?;
        Some(Data { root: self.root.clone(), at: Rc::from(format!("{}/{}", self.at, key.replace('~', "~0").replace('/', "~1"))) })
    }

    fn nth(&self, i: usize) -> Data {
        Data { root: self.root.clone(), at: Rc::from(format!("{}/{i}", self.at)) }
    }

    pub fn value(&self) -> &serde_json::Value {
        self.root.pointer(&self.at).expect("a place made from the document")
    }
}

/// What value `v`'s quantities were said by (units::Said): a field's, or
/// a vec's fields'.
fn said(v: &Value) -> Vec<&units::Said> {
    match v {
        Value::Field(_, Some(s)) => vec![s],
        Value::Vec(items) => items.iter().flat_map(said).collect(),
        _ => Vec::new(),
    }
}

impl Value {
    /// What it is, for messages ("field", "frame").
    pub fn kind(&self) -> &'static str {
        match self {
            Value::Field(..) => "field",
            Value::Point(..) => "point",
            Value::Str(_) => "string",
            Value::Vec(_) => "vec",
            Value::Data(_) => "data",
            Value::Fn(_) => "function",
            Value::Frame(_) => "frame",
            Value::Relation(_) => "relation",
        }
    }

    /// The fields it holds, in order, onto `out`: a frame's, a point's, a
    /// vec's and a relation's; a closure's captured values only with
    /// `closures` (a promise's inputs treat a function as code).
    pub(crate) fn fields_into(&self, out: &mut Vec<NodeId>, closures: bool) {
        match self {
            Value::Field(f, _) => out.push(*f),
            Value::Point(p) => out.extend(p),
            Value::Vec(items) => items.iter().for_each(|x| x.fields_into(out, closures)),
            Value::Fn(Callable::Closure(c)) if closures => c.env.iter().for_each(|(_, x)| x.fields_into(out, closures)),
            Value::Frame(xf) => out.extend(xf.r.iter().flatten().chain(&xf.t)),
            Value::Relation(r) => {
                out.extend(r.fields());
                r.data.iter().for_each(|(_, x)| x.fields_into(out, closures));
            }
            Value::Str(_) | Value::Data(_) | Value::Fn(_) => {}
        }
    }
}

#[derive(Debug)]
pub struct FnDef {
    pub params: Vec<String>,
    /// The last param takes the remaining arguments, as a vec (`rest ...`).
    pub rest: bool,
    pub body: Sexp,
    pub file: Rc<str>,
    pub scope: Scope,
}

/// Whether a param is a (ranged lo hi) no def named: named by its form and
/// a number ("t-slot-reaches 0"), a space no written name has.
pub fn unnamed_range(name: &str) -> bool {
    name.contains(' ')
}

pub struct Param {
    pub name: String,
    pub default: f64,
    /// A tolerance, `(tol lo hi)`: deviations from the value (lo <= hi),
    /// so proofs cover every value in [value + lo, value + hi] (`band`).
    pub tol: Option<(f64, f64)>,
}

impl Param {
    /// Its tolerance band around `value`: the real sums value + lo and
    /// value + hi, each the float itself where the sum is exact, else
    /// rounded outward (so the band holds the real one). None without a
    /// tolerance, or where an end overflows.
    pub fn band(&self, value: f64) -> Option<(f64, f64)> {
        let (lo, hi) = self.tol?;
        let (a, b) = (add_rounded(value, lo, false), add_rounded(value, hi, true));
        (a.is_finite() && b.is_finite()).then_some((a, b))
    }

    /// Its band around its declared value, where that excludes the value
    /// itself (a drawing's basic size with a fit's deviations: 8 f9 is
    /// 7.951..7.987): the value is a size no real part has, which the
    /// commands that work at one value would measure. Why, as a warning.
    pub fn outside_band(&self) -> Option<String> {
        let (v, (a, b)) = (self.default, self.band(self.default)?);
        if a <= v && v <= b {
            return None;
        }
        let num = |x: f64| (x * 1e9).round() / 1e9;
        let (a, b, mid) = (num(a), num(b), num(0.5 * a + 0.5 * b));
        let (name, half) = (&self.name, num(0.5 * (b - a)));
        Some(format!(
            "param {name} is {v} with its band {a}..{b}, which excludes {v} (a basic size with a fit's deviations): assert and interfere cover the band, but the commands that work at one value (distance, mass, eval, render, mesh) measure {v}, a size no real part has, and say so (outside_band). \
             Write the value as a size in the band, e.g. (param {name} {mid} (tol -{half} {half})), or measure one with --param {name}=V"
        ))
    }
}

/// a + b rounded down (or up): the float sum, stepped outward by one ulp
/// where it is not exact (Knuth's TwoSum gives the rounding error exactly).
fn add_rounded(a: f64, b: f64, up: bool) -> f64 {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    match (up, err) {
        (false, e) if e < 0.0 => s.next_down(),
        (true, e) if e > 0.0 => s.next_up(),
        _ => s,
    }
}

/// The sine and cosine of a turn of `deg` degrees, every placement's (`at`
/// steps, joints, a motion's keys; std's rotate-x/y/z are warps of their
/// own, read in radians: off a quarter turn by rounding, about 1e-16, far
/// below any claim's tolerance). A constant turn is reduced to within
/// 45 degrees of a quarter turn in degrees, exactly, before it is read in
/// radians: a multiple of 90 is exactly 0 and +-1, where sin(180 pi/180)
/// is 1.2e-16 (an axis turned half round would be off by that, and a
/// claim that two axes are parallel would fail). Others are as near as
/// before. A turn reading params or t stays the graph's sin and cos.
pub(crate) fn turn_sin_cos(g: &mut Graph, deg: NodeId) -> [NodeId; 2] {
    if let Some(d) = g.as_const(deg).filter(|d| d.is_finite()) {
        let r = d % 360.0;
        let q = (r / 90.0).round();
        let (s, c) = (r - 90.0 * q).to_radians().sin_cos();
        let (s, c) = match q.rem_euclid(4.0) as u8 {
            0 => (s, c),
            1 => (c, -s),
            2 => (-s, -c),
            _ => (-c, s),
        };
        // + 0: no -0 entries.
        return [g.constant(s + 0.0), g.constant(c + 0.0)];
    }
    let k = g.constant(std::f64::consts::PI / 180.0);
    let rad = g.mul(deg, k);
    [g.un(UnOp::Sin, rad), g.un(UnOp::Cos, rad)]
}

/// Rigid placement: `p_world = R * p_local + t`. Entries are graph nodes, so
/// placements can depend on parameters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Xf {
    pub r: [[NodeId; 3]; 3],
    pub t: [NodeId; 3],
}

/// A placed body's own frame: its shape there, and its placement. Its world
/// shape is the shape at R^T (p - t) (`Doc::items`), and an integral over
/// that is one over this shape of the integrand at R p + t, with the
/// placement's determinant (claims: a body turned through a ranged angle
/// is integrated where it does not move).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Placed {
    pub shape: NodeId,
    pub at: Xf,
}

pub struct Body {
    pub name: String,
    /// Its arguments, with the values it was built with as written (as
    /// an assembly's: `Assembly::params`).
    pub params: Vec<(String, String)>,
    /// For one built with arguments: the body it is an instance of.
    pub of: Option<String>,
    pub shape: NodeId,
    pub attrs: Attrs,
    /// Declared frames, in the body's own coordinates.
    pub frames: Vec<(String, Xf)>,
    /// Marks, curves drawn on the part (a sheet's bend lines): each its
    /// name and a field pair (g, r), the curve where g = 0 inside r < 0,
    /// in the body's own coordinates.
    pub marks: Vec<(String, NodeId, NodeId)>,
    /// The file that defines it (as named to the loader, or an import's
    /// path): where it is written, or, built with arguments, where they
    /// were passed.
    pub file: Rc<str>,
}

/// A body's attributes, shared by every item that places it (flattening an
/// assembly of thousands of parts copies no catalog records).
pub type Attrs = Rc<Vec<(String, Value)>>;

pub struct Instance {
    pub name: String,
    pub child: String,
    /// The final placement, in the assembly's coordinates.
    pub at: Xf,
    pub joint: Option<Joint>,
    /// The instances (paths in its assembly) whose frames WHERE reads: the
    /// child is placed on them, joined at an interface (a screw on its tap,
    /// a carriage on its rail, a hinge with clearance), so connected to
    /// them, and what the instance carries loads them
    /// (a reaction). A sub-assembly's frame of one of its own instances
    /// names that instance (`s/output` for a servo's output frame). One
    /// placed in the assembly's own coordinates is placed on nothing.
    pub on: Vec<String>,
    /// The child's frame met face to face with WHERE (CFRAME), if one was.
    pub met: Option<String>,
}

/// An instance's one degree of freedom: `(turn AXIS deg)` or `(slide AXIS
/// mm)` about or along an axis of its WHERE frame (`anchor`), driven by an
/// expression of params and t (`text` as written), within `limits`.
#[derive(Clone)]
pub struct Joint {
    pub turn: bool,
    pub axis: usize,
    pub anchor: Xf,
    /// The child's placement below the motion, in the frame it moves:
    /// C^-1 then FLIP where a CFRAME is given, else the identity. The
    /// instance's placement is this, then the motion, then WHERE.
    pub below: Xf,
    pub drive: NodeId,
    pub text: String,
    pub limits: Option<(NodeId, NodeId)>,
    /// `(speed-limits LO HI)`: bounds on the drive's rate, d drive / dt.
    pub speed_limits: Option<(NodeId, NodeId)>,
    /// Its clauses that packages own (Package::joint_clauses), in order,
    /// each at where the motion is set (where the arguments were passed,
    /// else the instance), as its limits' claims are.
    pub clauses: Vec<OwnedClause>,
}

impl Joint {
    /// What package clause `head` kept, if the joint has one.
    pub fn clause<T: std::any::Any>(&self, head: &str) -> Option<&T> {
        self.clauses.iter().find(|c| c.head == head).and_then(|c| c.kept.downcast_ref::<T>())
    }
}

pub struct Assembly {
    pub name: String,
    /// Its parameters, with the values it was built with as written (the
    /// defaults, or an instance's arguments).
    pub params: Vec<(String, String)>,
    /// For one built with arguments: the assembly it is an instance of
    /// (named after the instance path that first passed them).
    pub of: Option<String>,
    pub instances: Vec<Instance>,
    /// Declared frames, in the assembly's coordinates.
    pub frames: Vec<(String, Xf)>,
    /// The instances (paths in it) each frame is of (Instance::on).
    pub frame_on: Vec<(String, Vec<String>)>,
    /// The file that defines it, as a body's (`Body::file`).
    pub file: Rc<str>,
    /// What is declared of pairs of its instances: `(contact A B ZONE)` and
    /// `(relate A B KIND)`.
    pub relations: Vec<Relation>,
    /// Its clauses that packages own (Package::clauses), in order: placed
    /// with each instance of it (Doc::clauses_under).
    pub clauses: Vec<OwnedClause>,
    /// Its published `(attr mass M)`, on the assembly that defines it.
    pub mass: Option<AssemblyMass>,
}

/// A clause a package owns (an assembly's, a joint's or a relation's:
/// Package::clauses): the package and the clause's head, what the package
/// kept of it (Package::clause), what it read reading it (handed on with
/// it: Reader::clauses_under), and where it is written.
#[derive(Clone, Debug)]
pub struct OwnedClause {
    pub package: &'static str,
    pub head: &'static str,
    pub kept: std::sync::Arc<dyn std::any::Any + Send + Sync>,
    pub reads: Reads,
    pub file: Rc<str>,
    pub line: u32,
}

/// A package's clause under a target (Doc::clauses_under): the path of the
/// assembly instance it is written in (`root/arm`, as `Item::label`s
/// begin), where that assembly is placed in the world, and the clause, in
/// that assembly's coordinates.
pub struct PlacedClause {
    pub path: String,
    pub at: Xf,
    pub clause: OwnedClause,
}

/// An assembly's published mass (grams) and the density it gives the
/// bodies it places with none of their own, the derived param
/// `NAME.mass-density`: None where every one has its own, so it is not
/// used.
pub struct AssemblyMass {
    pub mass: NodeId,
    pub density: Option<NodeId>,
}

/// A relation of instance paths `a` and `b`: overlap of them inside `zone`
/// is intended (`interfere` leaves it out of their pair). A contact is one
/// with a zone alone; `(relate A B KIND)` gives its kind too.
pub struct Relation {
    pub a: String,
    pub b: String,
    pub zone: Option<Zone>,
    pub kind: Option<Rc<Kind>>,
}

/// Where two parts are meant to overlap: the intersection of `a`, a field
/// in A's own coordinates, and `b`, one in B's (none: everywhere).
#[derive(Clone, Copy, Debug)]
pub struct Zone {
    pub a: NodeId,
    pub b: Option<NodeId>,
}

/// A relation's kind as `(relation ...)` gave it: its name, what its zone
/// leaves unchecked, the names of the claims it added and its data.
pub struct Kind {
    pub name: String,
    pub skipped: Option<String>,
    pub claims: Vec<String>,
    pub data: Vec<(String, Value)>,
    /// Its clauses that packages own, in the relating assembly's
    /// coordinates.
    pub clauses: Vec<OwnedClause>,
    /// Where it is related (file:line).
    pub site: String,
}

/// A relation under a target, as `interfere` reads it: its ends as full
/// paths (matching `Item::label`), its zone in world coordinates, its zone
/// as written (`local`: each part in its end's own coordinates), and where
/// the relating assembly is placed in the world (its kind's package
/// clauses are in its coordinates).
pub struct PlacedRelation {
    pub a: String,
    pub b: String,
    pub zone: Option<NodeId>,
    pub local: Option<Zone>,
    pub kind: Option<Rc<Kind>>,
    pub at: Xf,
}

/// What a target path names (`root`, `root/outer/inner`, an assembly's
/// or a body's definition name), resolved once (Doc::target): every
/// question of it reads this, so a path and the definition it places
/// are asked alike.
pub struct Target {
    pub path: String,
    /// The body or assembly it places.
    pub def: String,
    /// The file that defines `def`: whose span it is asked over.
    pub file: Rc<str>,
    /// Where it is placed in its top assembly's coordinates.
    pub at: Xf,
    /// The assemblies enclosing it, from the top: each one's path, its
    /// definition and where it is placed. Their relations between parts
    /// of the target are the target's too.
    pub above: Vec<(String, String, Xf)>,
}

const INSTANCE: &str = "(instance NAME CHILD WHERE [CFRAME] [MOTION] [(limits LO HI [(accuracy x)])] [(speed-limits LO HI [(accuracy x)])] \
    [CLAUSE ...]): CHILD a body or assembly, or (CHILD ARG ...) to pass its arguments; WHERE a frame in this assembly's \
    coordinates, e.g. (at (move 0 0 5)) or (frame-of b top); CFRAME the name of a frame CHILD declares, met face to face with \
    WHERE; MOTION (turn AXIS deg) or (slide AXIS mm), AXIS x, y or z of WHERE; limits bound the motion's value, speed limits \
    its rate (deg/s or mm/s, t in seconds); (accuracy x), how far past them the claim allows, in their units; CLAUSE a package's joint clause";

/// The assembly whose clauses are being read: its instances so far (what
/// a relative `frame-of` path starts at), the names of those placed
/// later (naming one is an error), and, built with arguments, where they
/// were passed: the claims its joints add are that file's, at that line.
/// `read`: the instance paths `frame-of` read since it was last emptied.
struct Building {
    asm: Assembly,
    later: Vec<String>,
    passed: Option<(Rc<str>, u32)>,
    /// Its `(attr mass M)` and `(attr mass-params ...)`, and where.
    mass: Vec<(String, Value, Sexp)>,
    read: Vec<String>,
}

/// A body or an assembly with arguments, `(body (NAME (P DEFAULT) ...)
/// ...)` or `(assembly (NAME (P DEFAULT) ...) ...)`: its clauses as
/// written, built again for each instance passing other values, in the
/// names of the file it is written in (Scope).
pub struct ArgDef {
    pub params: Vec<(String, Sexp)>,
    pub clauses: Vec<Sexp>,
    pub file: Rc<str>,
    pub scope: Scope,
}

/// A claim's comparison against its bound (the kernel's, which proves it
/// over intervals).
pub use ink_kernel::ranged::Cmp;

#[derive(Clone)]
pub enum Claim {
    /// The field's solid (where it is < 0) holds no point deeper than the
    /// claim's accuracy (none at all, to 0).
    Empty(NodeId),
    /// The two shapes' overlap holds no point deeper than the accuracy:
    /// Empty of their intersection, kept as a pair so motion can be checked
    /// as their distance.
    Disjoint(NodeId, NodeId),
    /// The first shape sticks out of the second nowhere by more than the
    /// accuracy: Empty of max(a, -b) (`within`: b a box).
    Subset(NodeId, NodeId),
    /// Values in a region (`range`): Empty of the field given, which is
    /// below 0 where a value is out of its range inside the region,
    /// decided exactly (a value has no length to be within).
    Range(NodeId),
    /// The shapes overlap deeper than the accuracy somewhere.
    Overlaps(NodeId, NodeId),
    /// The shapes sit against each other: they overlap nowhere deeper than
    /// the accuracy, and come within it of each other somewhere.
    Touches(NodeId, NodeId),
    /// A point is inside the shape (`inside`), or not: the field's sign
    /// there, decided exactly (a point has no size to be within an
    /// accuracy of). Its coordinates are numbers (constants and params).
    Point { f: NodeId, at: [NodeId; 3], inside: bool },
    /// No point of the second field is within distance d of the first
    /// (sound for any fields, via gradient bounds).
    Clearance(NodeId, NodeId, f64),
    /// A package's metric of one of its forms or of a shape compared with
    /// a bound (`(max-stress books < 30)`, `(min-thickness part > 1)`):
    /// measured by the package, numerically. `excluding`: what the
    /// metric leaves out (`Excluded::key`) that the claim leaves out on
    /// purpose, `(excluding KEY ...)`: it answers on the rest, saying so.
    Metric { metric: FormMetric, cmp: Cmp, value: f64, excluding: Vec<&'static str> },
    /// A number (any expression of params and t) compared with a bound.
    Value(NodeId, Cmp, f64),
    /// The sum over parts of the integral of a field over each part's
    /// shape, for each of one or more fields: `parts` are (shape, the
    /// fields), one per body of the target (a target with no bodies is one
    /// part, of an empty shape), each with the same number of fields, and
    /// `frames` each part's body's own frame, where it is a placed body.
    /// Each sum is compared with the bound.
    Integral { parts: Vec<(NodeId, Vec<NodeId>)>, frames: Vec<Option<Placed>>, cmp: Cmp, bound: f64 },
    /// The weighted mean of a field over parts, the sum of the integrals of
    /// w h over the sum of those of w, for each of one or more fields h:
    /// `parts` are (shape, w, the h's), one per body of the target (as for
    /// an integral, with their `frames`). Each mean is compared with the
    /// bound. Weighted by density, the mean of x is the centre of mass's x.
    Mean { parts: Vec<(NodeId, NodeId, Vec<NodeId>)>, frames: Vec<Option<Placed>>, cmp: Cmp, bound: f64 },
}

/// What a claim's accuracy, `(accuracy x)`, is in (Claim::accuracy).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Accuracy {
    /// A length: how deep a claim about shapes may be violated (by default
    /// the model's accuracy).
    Length,
    /// The claim's own units: how far a comparison's value may pass its
    /// bound (by default 0, exactly).
    Value,
    /// None: decided exactly (a point's sign, values in a region).
    Exact,
}

impl Claim {
    /// What the claim's accuracy is in: a length for the claims about
    /// shapes; the value's units for comparisons with a bound (values,
    /// volumes, integrals, means, metrics); points and ranges
    /// are decided exactly.
    pub fn accuracy(&self) -> Accuracy {
        match self {
            Claim::Empty(_) | Claim::Disjoint(..) | Claim::Subset(..) | Claim::Clearance(..) | Claim::Overlaps(..) | Claim::Touches(..) => Accuracy::Length,
            Claim::Value(..) | Claim::Integral { .. } | Claim::Mean { .. } | Claim::Metric { .. } => Accuracy::Value,
            Claim::Point { .. } | Claim::Range(_) => Accuracy::Exact,
        }
    }

    /// Whether the claim is checked to a length: the claims about shapes.
    pub fn accuracy_is_length(&self) -> bool {
        self.accuracy() == Accuracy::Length
    }

    /// The fields checking the claim reads: its shapes, integrands and
    /// values (a package's metric reads its own apart),
    /// and of a part integrated in its body's frame, the body's shape and
    /// its placement (what it is integrated over there).
    pub fn fields(&self) -> Vec<NodeId> {
        let framed = |frames: &[Option<Placed>]| -> Vec<NodeId> {
            frames.iter().flatten().flat_map(|f| std::iter::once(f.shape).chain(f.at.r.into_iter().flatten()).chain(f.at.t)).collect()
        };
        match self {
            Claim::Empty(f) | Claim::Range(f) | Claim::Value(f, ..) => vec![*f],
            Claim::Clearance(p, q, _) | Claim::Disjoint(p, q) | Claim::Subset(p, q) | Claim::Overlaps(p, q) | Claim::Touches(p, q) => vec![*p, *q],
            Claim::Point { f, at, .. } => std::iter::once(*f).chain(at.iter().copied()).collect(),
            Claim::Integral { parts, frames, .. } => parts.iter().flat_map(|(s, gs)| std::iter::once(*s).chain(gs.iter().copied())).chain(framed(frames)).collect(),
            Claim::Mean { parts, frames, .. } => parts.iter().flat_map(|(s, w, hs)| [*s, *w].into_iter().chain(hs.iter().copied())).chain(framed(frames)).collect(),
            Claim::Metric { metric, .. } => metric.shape.into_iter().collect(),
        }
    }
}

/// The names a file's forms define at top level, as written: values,
/// functions, bodies, assemblies, claims and package forms.
fn top_names(forms: &[Sexp], registry: &Registry) -> HashSet<String> {
    // Core forms that define a name; every extension form does.
    const DEFINING: &[&str] = &["def", "fn", "param", "body", "assembly", "assert"];
    let mut all = HashSet::new();
    for f in forms {
        let Some([h, n, ..]) = f.list() else { continue };
        // A body or an assembly with arguments is named by its list's head.
        let n = match n.list() {
            Some([n, ..]) if matches!(h.sym(), Some("assembly" | "body")) => n,
            _ => n,
        };
        let (Some(h), Some(n)) = (h.sym(), n.sym()) else { continue };
        if DEFINING.contains(&h) || registry.package(h).is_some() {
            all.insert(n.to_string());
        }
    }
    all
}

/// The directories an import not found beside its file is read from, in
/// order: `$INK_PATH`'s, separated as `$PATH`'s are. A host or a
/// distribution gives its libraries so (a parts catalog: `parts/` in one).
pub fn import_dirs() -> Vec<std::path::PathBuf> {
    let dirs = std::env::var_os("INK_PATH").unwrap_or_default();
    std::env::split_paths(&dirs).filter(|d| !d.as_os_str().is_empty()).map(|d| d.canonicalize().unwrap_or(d)).collect()
}

/// The file `(import "p")` in a file in `base` reads: `p` relative to
/// `base` (or absolute); else the first of `import_dirs` holding it (none:
/// `p` relative to `base`). And the files looked for before it, absent: a
/// load rests on their absence as on what it reads.
pub fn import_path(base: &Path, p: &str) -> (std::path::PathBuf, Vec<std::path::PathBuf>) {
    let beside = base.join(p);
    if beside.exists() || Path::new(p).is_absolute() {
        return (beside, vec![]);
    }
    let mut tried = vec![beside.clone()];
    for f in import_dirs().into_iter().map(|d| d.join(p)) {
        if f.exists() {
            return (f, tried);
        }
        tried.push(f);
    }
    tried.remove(0);
    (beside, tried)
}

/// A file as ink names it to the user (in errors, claims' sources, `ink
/// check`'s imports): one in an import directory (`import_dirs`) by its
/// path there, the path that imports it from anywhere (an agent of a
/// hosted ink opens it so); any other as given.
pub fn shown(path: &Path) -> String {
    let real = path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    match import_dirs().iter().find_map(|d| real.strip_prefix(d).ok().filter(|r| !r.as_os_str().is_empty())) {
        Some(r) => r.display().to_string(),
        None => path.display().to_string(),
    }
}

pub struct AssertDef {
    pub name: String,
    pub claim: Claim,
    pub text: String,
    /// The claim's bound as a field of params (its value is in `claim`, at
    /// load time).
    pub bound: Option<NodeId>,
    /// The cell size the claim says it needs, `(res r)`: its search reaches
    /// the finer of this and a res given to the run.
    pub res: Option<f64>,
    /// The accuracy the claim states, `(accuracy x)`, in its units
    /// (Claim::accuracy; None: its default, claims::claim_accuracy).
    pub accuracy: Option<f64>,
    /// What the claim's value and accuracy are in, where known (None: the
    /// model's own units): Doc::accuracy_quantity, or a joint's limits'.
    pub quantity: Option<Quantity>,
    /// Where it is written: the file (as named to the loader, or an
    /// import's path) and line.
    pub file: Rc<str>,
    pub line: u32,
    /// Named parts of a value claim, reported where its value is (a joint's
    /// drive and limits at the witness).
    pub terms: Vec<(String, NodeId)>,
    /// Ranged params its value was read with at load, each one number there
    /// (a relation's kind reading its parts' shapes by bounds): a run ranging
    /// one refuses the claim, not the file.
    pub fixed: Vec<usize>,
}

/// A body placed in world space (the result of resolving a target): its
/// field and attributes in its own frame, and its placement. One body
/// placed many times is one field.
pub struct Item {
    pub label: String,
    /// The body's field in its own frame (a def's as written). Where it
    /// lies in the world is Doc::placed_field, built only where read.
    pub shape: NodeId,
    /// The body's attributes in its frame (placed: Doc::placed_attrs).
    pub attrs: Attrs,
    /// For a placed body: the body's name and its placement, so queries can
    /// work once per distinct body in its own frame and transform the results
    /// (instancing).
    pub instance: Option<(String, Xf)>,
    /// The instances from the world down to it, top first, the target's
    /// own path among them (none for a body or a def): its placement is
    /// their composition, so two items are placed relative to each other by
    /// the links they do not share, and moved and stretched by those they
    /// do, whatever target named them (motion::over_joint, pairs::framed).
    pub links: Vec<Link>,
}

/// An instance on the way from a target to an item: its path, its
/// placement in its parent's coordinates, and its joint.
#[derive(Clone)]
pub struct Link {
    pub path: String,
    pub at: Xf,
    pub joint: Option<Joint>,
    /// The instances its WHERE reads (Instance::on).
    pub on: Vec<String>,
}

pub struct Doc {
    /// Results packages solved, by content key (a package's own type).
    pub solved: HashMap<u128, std::sync::Arc<dyn std::any::Any + Send + Sync>>,
    /// The packages whose forms this document can read.
    pub registry: &'static Registry,
    pub g: Graph,
    pub params: Vec<Param>,
    pub defs: Vec<(String, Value)>,
    def_index: HashMap<String, usize>,
    pub fns: HashMap<String, Rc<FnDef>>,
    pub bodies: Vec<Body>,
    pub assemblies: Vec<Assembly>,
    /// Bodies and assemblies with arguments, by name: what an instance
    /// passing arguments builds again.
    pub arg_defs: HashMap<String, Rc<ArgDef>>,
    /// Those built so far, by the definition and the values passed: the
    /// name of what they built (instances passing equal values share it).
    built: HashMap<String, String>,
    /// The values of the arguments of the body or assembly whose clauses
    /// are being read: what they read the arguments' names as.
    args: Env,
    pub asserts: Vec<AssertDef>,
    /// Extension forms (`study`, `optimize`, ...), as their packages read them:
    /// see `forms` and `form`.
    pub ext: Vec<ExtForm>,
    /// The forms of `ext` that provide params (solves), each after every
    /// one that sets a param it reads: the order providers are applied in
    /// (`providers`), at load time, at every context, and over ranged
    /// params (ModelView::scope_spans).
    pub solves: Vec<Solve>,
    /// Densities derived from a published mass (`(attr mass M)`): each
    /// param, and what it is the mass of, as answers name it.
    pub published: Vec<(usize, String)>,
    /// The params of solves on demand that something reads (Doc::want):
    /// provided at every context.
    wanted: HashSet<usize>,
    /// The last successful provision, by form: reuse it while
    /// its inputs and warm-start values are unchanged.
    provided: std::cell::RefCell<HashMap<usize, ext::Provided>>,
    /// The solves compiled, those on t followed along t (`follow`), until
    /// the solves or their load-time values change; and the steps they
    /// were followed by, kept across that.
    follow: std::cell::RefCell<Option<std::sync::Arc<follow::Follow>>>,
    paths: follow::Paths,
    /// The run file's `(t-span A B)`: the motion its claims that read t
    /// are checked over by default.
    pub t_span: Option<(f64, f64)>,
    /// Each file's `(t-span A B)`, the run file's and its imports', by the
    /// file's name (as `AssertDef::file`): a claim is checked over the span
    /// of the file it is written in, an assembly's motion over its file's
    /// (`span_of`).
    pub file_spans: HashMap<Rc<str>, (f64, f64)>,
    /// World fields of placed bodies (an instance, a construction in one,
    /// `place`), with the body's own field and the placement: a body that
    /// moves rigidly, which motion checks bound by how far its points move.
    /// Keyed by field, so one of the decompositions of equal fields: an
    /// item's own body and placement are its instance's and links'.
    pub placed: HashMap<NodeId, Placed>,
    /// Each file imported, by path and namespace ("std", a path, "PATH as
    /// ns"): run once, where it was first imported ("b.ink:2").
    imported: HashMap<String, String>,
    /// The values the next file run takes (SourceFile::given).
    giving: Vec<ImportValue>,
    /// Each param an import set (ImportValue), by index: its written
    /// default and where the import is ("b.ink:4").
    pub set_at: HashMap<usize, (f64, String)>,
    /// Files whose run under a namespace defined values and functions
    /// alone (a part family's file: its tables and its functions), reading
    /// nothing that differs by namespace (`varies`), by path: the file and
    /// the namespace they were read under. Imported under another
    /// namespace, such a file is not run again: its names are bound there
    /// to what that run defined (`library_under`).
    libraries: HashMap<String, (u32, String)>,
    /// How many reads of what differs by namespace the load has made: a
    /// file imported that is not std or a library. A run during which it
    /// grew is no library.
    varies: usize,
    /// The JSON files read (`read-json`), by path: read once, shared.
    json: HashMap<std::path::PathBuf, Data>,
    /// The fingerprint of each file's content the load read (`Doc::read`):
    /// see `sources`.
    sources_read: Vec<u128>,
    /// The files on disk the load depended on: see `files_read`.
    files_read: Vec<(std::path::PathBuf, Option<u128>)>,
    /// The namespace of the file being run (`import ... as ns`), if any.
    ns: Option<String>,
    /// Params a solve sets (by index): numerical at their values (a
    /// numerical method's answer), certified where a view derives them over
    /// ranged params (ModelView::set_param_spans).
    pub solved_params: HashSet<usize>,
    /// The variables `(ranged lo hi)` made (params at lo, no tolerance),
    /// each with [lo, hi]: a range of the question, which every claim
    /// reading one ranges over and no run pins (ModelView::declared_spans).
    pub ranged: Vec<(usize, (f64, f64))>,
    /// Each solved param's rate along each input asked for (Doc::rate):
    /// the param its rates' solve sets, or None where it does not vary.
    rates_of: HashMap<(usize, Var), Option<NodeId>>,
    /// Params a derivative held fixed along an input (no solve set them
    /// when it was taken), with where: a solve read later may not make one
    /// vary along it (Doc::varies_held).
    held: HashMap<(usize, Var), String>,
    /// Joints' speed-limit claims (index in `asserts`), each with its
    /// drive and limits: their rates are taken once the file is read
    /// (Doc::rate_joints), through every solve their drives read, however
    /// the file orders them (a loop closed over the assembly's frames
    /// comes after it).
    rated: Vec<(usize, NodeId, NodeId, NodeId)>,
    /// The solves read so far do not all solve at load time: the first
    /// that fails, as an error at its form. An error once the document is
    /// read, unless a solve read later makes them solve (Doc::solve_at_load).
    unsolved: Option<InkError>,
    /// The solves depending on t that do not solve at load time (t = 0),
    /// by index in `solves`, with why: an error once the document is
    /// read if a value resting on one is used at load time.
    unsolved_t: Vec<(usize, InkError)>,
    /// Directions that must have a length wherever the model is used
    /// (Doc::unit): a rotation's axis, a frame's, that is not a constant.
    pub conditions: Vec<Condition>,
    /// The bound of the claim being parsed (see `claim`).
    bound: Option<NodeId>,
    /// Warnings found while loading (undefined cases in constants).
    pub load_warnings: Vec<String>,
    /// What each param's quantity was said by, by index, where its value
    /// was (units::Said): a slot taking a known quantity checks it.
    param_said: HashMap<usize, Rc<units::Said>>,
    /// Conditions builders stated (`check`), which `ink check` reports
    /// where they fail.
    pub checks: Vec<check::Check>,
    /// Each check's test (its kind, comparison and the nodes it reads), to
    /// state each once.
    checks_seen: HashSet<(u8, Vec<NodeId>)>,
    /// The calls being evaluated, outermost first: the function's name
    /// and the file and line it is called at, for where a check was
    /// stated.
    calls: Vec<(String, Rc<str>, u32)>,
    /// Each field placed (Doc::warped), by its field and placement.
    placed_fields: HashMap<(NodeId, Xf), NodeId>,
    /// Which of x, y and z each field placed reads (Doc::xyz_reads).
    xyz_read: HashMap<NodeId, [bool; 3]>,
    /// While reading an integral's or a mean's integrand: the body being
    /// integrated (its label) and its attributes in world space, which
    /// `(attr NAME)` reads; for a target that is a field or has no bodies,
    /// why it has none.
    integrand_of: Option<(String, Result<Attrs, &'static str>)>,
    /// The file being run: while loading, the one whose forms are read;
    /// after, the file loaded.
    file: Rc<str>,
    /// The scope names are read in: the file's at top level, a function's
    /// in its body (Scope).
    scope: Scope,
    files: Vec<SourceFile>,
    /// The files being run, innermost last (an import cycle is an error).
    running: Vec<u32>,
    /// The file and line of each value key, and of the last public name
    /// (ordinary names are unique; std names may be shadowed per file).
    origin: HashMap<String, (u32, u32)>,
    param_keys: HashMap<String, usize>,
    depth: usize,
    /// Per call being evaluated, innermost last: how many entries its env
    /// starts with that are its inputs (what a closure captured, then the
    /// arguments; the body's lets follow). A promise made in the body holds
    /// with these fixed (`lipschitz`).
    inputs: Vec<usize>,
    /// How many imports deep the file being run is (0: the file given).
    import_depth: usize,
    /// Parameter values used for load-time computations such as (bounds f).
    overrides: Vec<(String, Given)>,
    /// Params whose values were consumed while loading (load-time numbers,
    /// sampled or solved fields, measurements, claim bounds): field
    /// derivatives do not see them.
    baked: std::cell::RefCell<HashSet<usize>>,
    /// Of those, the ones whose value became a number the model keeps
    /// (everything but a claim's bound, which is kept as an expression), and
    /// where first: the top-level form's file and line, and its head and
    /// name. A proof cannot cover such a param over a range.
    materialized: std::cell::RefCell<HashMap<usize, String>>,
    /// Undefined cases the current form's load-time numbers passed through
    /// (`load_eval`), for `at_form` to report.
    undefined_at_load: std::cell::RefCell<Vec<&'static str>>,
    /// `(shape-of path)`'s fields by path, while no top-level form but
    /// claims has run since: what is defined does not change between them,
    /// and a whole assembly's is a walk of every item (an assembly's claims
    /// read it again and again: the RC car's 100 tool-access claims, 3 s of
    /// its 3.4 s load).
    shape_of: HashMap<String, NodeId>,
    /// The top-level form being run: file:line (head name).
    form_at: String,
    /// The first number written outside ink's range and where (file:line):
    /// what the range's unknown cites (validity::numbers).
    pub(crate) wild_written: Option<(f64, String)>,
    /// Its name (the def's, the claim's), which a `ranged` it makes is named by.
    form_name: String,
    /// The assembly being read, if any.
    building: Option<Building>,
    /// Reading a relation end's FRAME, which is in the part's own
    /// coordinates: `frame-of` names none of the assembly's instances there.
    own_frame: bool,
    /// What the file run was loaded from, so it can be loaded again
    /// (`with_coordinates`).
    source: Option<Origin>,
    /// In a load `with_coordinates`: the assembly whose instances'
    /// arguments are built as coordinates.
    coordinates_of: Option<String>,
    /// In a load `with_states`: the coordinates that move are built as a
    /// state (Coord::jet), and this is the param of the moment t0 it is at.
    states: bool,
    pub state_t: Option<usize>,
    /// Those coordinates, in the order they were passed.
    pub coords: Vec<Coord>,
    /// Values that carry a param as their own coordinate's, each with it
    /// (Reader::as_coordinate): passed itself by the target's instance, the
    /// value is that coordinate, so what else reads it moves along it.
    pub coordinate_params: Vec<(NodeId, usize)>,
}

/// The source a document was loaded from (Doc::load_named's inputs).
#[derive(Clone)]
struct Origin {
    src: String,
    name: String,
    path: std::path::PathBuf,
    base: std::path::PathBuf,
}

/// A document's file run as it was loaded (Doc::reloader).
#[derive(Clone)]
pub struct Reload {
    registry: &'static Registry,
    origin: Origin,
    /// The param overrides it was loaded with.
    pub overrides: Vec<(String, Given)>,
}

impl Reload {
    /// Doc::with_coordinates.
    pub fn with_coordinates(&self, target: &str) -> R<Doc> {
        Doc::load_origin(self.registry, self.origin.clone(), &self.overrides, Some((target.to_string(), false))).map_err(|(e, _)| e)
    }

    /// Doc::with_states.
    pub fn with_states(&self, target: &str) -> R<Doc> {
        Doc::load_origin(self.registry, self.origin.clone(), &self.overrides, Some((target.to_string(), true))).map_err(|(e, _)| e)
    }
}

/// An argument an instance of an assembly passes (`INSTANCE.ARG`), in a
/// load `with_coordinates`: built as its value plus `param`, a param at 0,
/// so what reads the argument reads the param too, and a derivative along
/// the param is one along the argument. In a load `with_states`, one whose
/// value reads t is built as a state instead, `param + rate (t - t0) + acc
/// (t - t0)^2 / 2` (`jet`: rate's and acc's params; t0 is Doc::state_t),
/// so at t = t0 the mechanism is at that position, speed and acceleration
/// whatever the motion as written, and what reads its speed (a
/// torque-speed line's `(diff angle t)`) reads that speed.
#[derive(Clone, Debug)]
pub struct Coord {
    pub name: String,
    pub param: usize,
    /// The argument's value as written (no param): a field of params and t.
    pub value: NodeId,
    pub jet: Option<(usize, usize)>,
    /// What the instance is passed: value + param, or the state.
    pub passed: NodeId,
}

type Env = Vec<(String, Value)>;
type R<T> = Result<T, InkError>;

/// A direction that is not a constant (Doc::unit): its squared length,
/// which must be positive, and finite, wherever the model is used.
#[derive(Clone, Debug)]
pub struct Condition {
    pub len2: NodeId,
    /// What it is ("the axis of (rot a 0 0 60)") and the top-level form
    /// it is in (file:line (head name ...)).
    pub what: String,
    pub site: String,
    /// Where it is written, and how to give it a length: errors about it
    /// point there.
    pub at: InkError,
}

#[derive(Clone, Copy)]
enum Builtin {
    Un(UnOp),
    Bin(BinOp),
    Variadic(BinOp),
    Minus,
    Select,
}

/// The arithmetic forms, by name.
const BUILTINS: &[(&str, Builtin)] = {
    use Builtin::*;
    &[
        ("+", Variadic(BinOp::Add)),
        ("*", Variadic(BinOp::Mul)),
        ("min", Variadic(BinOp::Min)),
        ("max", Variadic(BinOp::Max)),
        ("emax", Variadic(BinOp::Emax)),
        ("-", Minus),
        ("/", Bin(BinOp::Div)),
        ("atan2", Bin(BinOp::Atan2)),
        ("pow", Bin(BinOp::Pow)),
        ("neg", Un(UnOp::Neg)),
        ("abs", Un(UnOp::Abs)),
        ("sqrt", Un(UnOp::Sqrt)),
        ("square", Un(UnOp::Square)),
        ("exp", Un(UnOp::Exp)),
        ("ln", Un(UnOp::Ln)),
        ("sin", Un(UnOp::Sin)),
        ("cos", Un(UnOp::Cos)),
        ("tan", Un(UnOp::Tan)),
        ("floor", Un(UnOp::Floor)),
        ("select", Select),
    ]
};

fn builtin(name: &str) -> Option<Builtin> {
    BUILTINS.iter().find(|(n, _)| *n == name).map(|&(_, b)| b)
}

const SPECIAL: &[&str] = &[
    "let",
    "vec",
    "unit",
    "nth",
    "data",
    "warp",
    "shape-of",
    "lambda",
    "map",
    "reduce",
    "apply",
    "range",
    "len",
    "concat",
    "char-codes",
    "slice",
    "if",
    "read-json",
    "get",
    "at",
    "axes",
    "frame-of",
    "place",
    "fold",
    "flatten",
    "bounds",
    "tag",
    "lipschitz",
    "grid-bake",
    "import-stl",
    "redistance",
    "study-field",
    "diff",
    "material",
    "by-body",
    "attr",
    "attr-of",
    "met-by",
    "check",
    "error",
    "profile-reach",
    "polygon",
    "revolve-z",
    "probe",
    "relation",
    "ranged",
];
/// The forms that are functions: each argument an expression, evaluated
/// once (no names, clauses or literal text), so they are values like a
/// user `fn`: `(map polygon outlines)`. The other forms are syntax.
const FORM_FUNCTIONS: &[&str] = &[
    "vec",
    "unit",
    "nth",
    "map",
    "reduce",
    "fold",
    "flatten",
    "apply",
    "range",
    "len",
    "concat",
    "char-codes",
    "slice",
    "get",
    "warp",
    "lipschitz",
    "bounds",
    "profile-reach",
    "polygon",
    "revolve-z",
    "probe",
    "ranged",
];
/// Builtins beyond the kernel ops: comparisons yielding 1 or 0 (as fields).
const COMPARISONS: &[&str] = &["<", "<=", ">", ">=", "=="];
const CONSTANTS: &[&str] = &["x", "y", "z", "t", "pi", "inf"];

/// Whether a value is, or holds, a point.
fn has_point(v: &Value) -> bool {
    match v {
        Value::Point(..) => true,
        Value::Vec(items) => items.iter().any(has_point),
        _ => false,
    }
}

fn is_reserved(name: &str) -> bool {
    builtin(name).is_some() || SPECIAL.contains(&name) || CONSTANTS.contains(&name) || COMPARISONS.contains(&name)
}

/// An error's hint naming what is_reserved keeps, words first (the
/// operators' symbols aside): every name hit is in it.
fn reserved_hint() -> String {
    let mut words: Vec<&str> = BUILTINS.iter().map(|(n, _)| *n).chain(SPECIAL.iter().copied()).filter(|n| n.starts_with(char::is_alphabetic)).collect();
    words.sort();
    format!("the built-in names: {}; {}", CONSTANTS.join(" "), words.join(" "))
}

/// `(fn NAME ARGS ...)` for a function written as a def: `(def (f x)
/// body)`, `(def f (x) body)` or `(def f (fn (x) body))`.
fn fn_written_as_def(args: &[Sexp]) -> Option<String> {
    let names = |e: &Sexp| e.list().filter(|l| l.iter().all(|a| a.sym().is_some())).map(|l| l.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(" "));
    let (name, ps) = match args {
        [head, _] if head.list().is_some_and(|l| l.len() > 1) => {
            let all = names(head)?;
            let (n, ps) = all.split_once(' ')?;
            (n.to_string(), ps.to_string())
        }
        [n, ps, _] => (n.sym()?.to_string(), names(ps)?),
        [n, e] => match e.list()? {
            [h, ps, _] if h.sym() == Some("fn") => (n.sym()?.to_string(), names(ps)?),
            _ => return None,
        },
        _ => return None,
    };
    Some(format!("(fn {name} ({ps}) body)"))
}

/// Typos between two names: letters inserted, deleted, changed, or two
/// swapped ('stduy' is one from 'study').
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let (mut before, mut prev): (Vec<usize>, Vec<usize>) = (Vec::new(), (0..=b.len()).collect());
    for i in 1..=a.len() {
        let mut cur = vec![i; b.len() + 1];
        for j in 1..=b.len() {
            let c = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + c);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                cur[j] = cur[j].min(before[j - 2] + 1);
            }
        }
        before = std::mem::replace(&mut prev, cur);
    }
    prev[b.len()]
}

/// std's forms, parsed once.
fn std_forms() -> &'static [Sexp] {
    static STD: std::sync::OnceLock<Vec<Sexp>> = std::sync::OnceLock::new();
    STD.get_or_init(|| parse(STD_SRC).expect("std parses"))
}

/// Whether std defines a function `name`.
fn in_std(name: &str) -> bool {
    static FNS: std::sync::OnceLock<std::collections::HashSet<String>> = std::sync::OnceLock::new();
    FNS.get_or_init(|| {
        std_forms()
            .iter()
            .filter_map(|f| match f.list() {
                Some([h, n, ..]) if h.sym() == Some("fn") => n.sym().map(String::from),
                _ => None,
            })
            .collect()
    })
    .contains(name)
}

/// The candidate nearest `name` by edit distance, as a hint, if a typo
/// could make one of the other: at most a third of the name's letters
/// differ (one, in a short name: 'hole' is not 'fold').
pub fn closest<'a>(name: &str, candidates: impl Iterator<Item = &'a str>) -> Option<String> {
    near(name, candidates).first().map(|c| format!("did you mean '{c}'?"))
}

/// What to pick instead of `name`, not one of `names` (`what`, e.g.
/// "its instances"): the nearest by `closest`, else all of them.
pub fn choices<'a>(name: &str, what: &str, names: impl Iterator<Item = &'a str> + Clone) -> String {
    let all = names.clone().collect::<Vec<_>>().join(", ");
    closest(name, names).unwrap_or_else(|| format!("{what}: {}", if all.is_empty() { "none" } else { &all }))
}

/// For an import `written` of a file `path` that cannot be read: the file
/// in its directory a typo (or an old extension) could make of its name.
fn sibling(path: &Path, written: &str) -> Option<String> {
    let name = path.file_name()?.to_str()?;
    let names: Vec<String> = std::fs::read_dir(path.parent()?).ok()?.filter_map(|e| e.ok()?.file_name().into_string().ok()).collect();
    let found = near(name, names.iter().map(String::as_str)).into_iter().next()?;
    let dir = written.rsplit_once('/').map_or("", |(d, _)| d);
    Some(format!("did you mean \"{dir}{}{found}\"?", if dir.is_empty() { "" } else { "/" }))
}

/// The candidates a typo could make of `name` (closest), nearest first:
/// each keeps a letter of it in place ('w' is no typo of 't').
fn near<'a>(name: &str, candidates: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let n = name.chars().count();
    let most = (n / 3).max(1);
    // An edit changes the length by one at most: a candidate longer or
    // shorter by more than `most` is further (not measured).
    let mut c: Vec<(usize, &str)> = candidates
        .map(|c| (c, c.chars().count()))
        .filter(|(_, k)| k.abs_diff(n) <= most)
        .filter_map(|(c, k)| Some((edit_distance(name, c), c)).filter(|&(d, _)| d <= most && d < n.max(k)))
        .collect();
    c.sort();
    c.into_iter().map(|(_, c)| c).collect()
}

/// Words agents bring from other CAD tools, with what std calls the same
/// thing, and how it is written where that is not a call by another name.
const CAD_WORDS: &[(&str, &[&str], &str)] = &[
    ("hole", &["subtract", "cylinder-z"], "a hole is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("bore", &["subtract", "cylinder-z"], "a bore is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("drill", &["subtract", "cylinder-z"], "a drilled hole is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("cut", &["subtract"], ""),
    ("difference", &["subtract"], ""),
    ("minus", &["subtract"], ""),
    ("fillet", &["fillet-union", "corner-fillet", "smooth-union", "rounded-box"], ""),
    ("round", &["floor", "fillet-union", "corner-fillet", "smooth-union", "rounded-box"], "a number rounds as (floor (+ v 0.5)); an edge, with a fillet"),
    ("chamfer", &["chamfer-union", "chamfer-intersect"], ""),
    ("bevel", &["chamfer-union", "chamfer-intersect"], ""),
    ("extrude", &["extrude-z", "sweep", "loft-z"], ""),
    ("prism", &["extrude-z"], ""),
    ("pad", &["extrude-z"], ""),
    ("revolve", &["revolve-z"], ""),
    ("lathe", &["revolve-z"], ""),
    ("translate", &["move"], ""),
    ("rotate", &["rotate-x", "rotate-y", "rotate-z"], ""),
    ("mirror", &["mirror-x", "mirror-y", "mirror-z"], ""),
    ("cube", &["box"], ""),
    ("cuboid", &["box"], ""),
    ("block", &["box", "rounded-box"], ""),
    ("ball", &["sphere"], ""),
    ("union", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("join", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("combine", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("intersect", &["max"], "an intersection is the max of its shapes: (max a b ...)"),
    ("intersection", &["max"], "an intersection is the max of its shapes: (max a b ...)"),
    ("offset", &["-"], "f grown by d is (- f d), shrunk by d (+ f d)"),
    ("array", &["pattern-linear", "pattern-polar", "repeat-linear", "repeat-x"], ""),
    ("hollow", &["shell"], ""),
    ("thicken", &["shell", "-"], ""),
    ("plane", &["half-space"], ""),
    ("ring", &["ring-2d", "torus-z"], ""),
    ("annulus", &["ring-2d"], ""),
    ("hex", &["hexagon-2d"], ""),
    ("hexagon", &["hexagon-2d"], ""),
    ("cone", &["cone-z"], ""),
    ("frustum", &["cone-z"], ""),
    ("taper", &["cone-z"], ""),
    ("countersink", &["cone-z"], ""),
    ("rectangle", &["rect-2d"], ""),
    ("screw", &["thread-iso"], ""),
    ("spring", &["helix-points", "sweep"], ""),
];

/// A hint for a call of an unknown function `name`: what std calls the
/// thing a CAD word names (CAD_WORDS: the whole name or a part of it
/// between dashes), the names one part of which a typo could make of it
/// ('cylnder': 'cylinder-z'), and the names a typo could make of it; `names` the
/// functions in scope, std's included where it is imported.
fn unknown_function_hint(name: &str, names: &[String], std_imported: bool) -> Option<String> {
    let known = |c: &str| names.iter().any(|n| n == c) || in_std(c);
    let mut pick: Vec<String> = Vec::new();
    let mut add = |c: &str| {
        if !pick.iter().any(|p| p == c) {
            pick.push(c.to_string());
        }
    };
    let words: Vec<&(&str, &[&str], &str)> = CAD_WORDS.iter().filter(|(w, ..)| name == *w || name.split('-').any(|s| s == *w || s == format!("{w}s"))).collect();
    for (_, std_names, _) in &words {
        std_names.iter().filter(|c| known(c)).for_each(|c| add(c));
    }
    if name.len() >= 3 {
        names.iter().filter(|c| c.contains('-') && !near(name, c.split('-')).is_empty()).for_each(|c| add(c));
    }
    near(name, names.iter().map(String::as_str)).into_iter().take(3).for_each(&mut add);
    if pick.is_empty() {
        return None;
    }
    pick.truncate(6);
    let quoted: Vec<String> = pick.iter().map(|c| format!("'{c}'")).collect();
    let mut hint = match quoted.as_slice() {
        [one] => format!("did you mean {one}?"),
        all => format!("did you mean one of {}?", all.join(", ")),
    };
    if let Some((.., how)) = words.iter().find(|(.., how)| !how.is_empty()) {
        hint += &format!(" ({how})");
    }
    if !std_imported && pick.iter().any(|c| in_std(c) && !names.contains(c)) {
        hint += " (std's: add (import \"std\"))";
    }
    Some(hint + " `ink std --section NAME` shows one.")
}

impl Doc {
    fn empty(registry: &'static Registry) -> Doc {
        Doc {
            solved: HashMap::new(),
            registry,
            g: Graph::new(),
            params: Vec::new(),
            defs: Vec::new(),
            def_index: HashMap::new(),
            fns: HashMap::new(),
            bodies: Vec::new(),
            assemblies: Vec::new(),
            arg_defs: HashMap::new(),
            built: HashMap::new(),
            args: Env::new(),
            asserts: Vec::new(),
            ext: Vec::new(),
            solves: Vec::new(),
            published: Vec::new(),
            wanted: HashSet::new(),
            provided: Default::default(),
            follow: Default::default(),
            paths: Default::default(),
            t_span: None,
            file_spans: HashMap::new(),
            placed: HashMap::new(),
            imported: HashMap::new(),
            giving: Vec::new(),
            set_at: HashMap::new(),
            libraries: HashMap::new(),
            varies: 0,
            json: HashMap::new(),
            sources_read: Vec::new(),
            files_read: Vec::new(),
            ns: None,
            load_warnings: Vec::new(),
            param_said: HashMap::new(),
            checks: Vec::new(),
            checks_seen: HashSet::new(),
            calls: Vec::new(),
            solved_params: HashSet::new(),
            ranged: Vec::new(),
            rates_of: HashMap::new(),
            held: HashMap::new(),
            rated: Vec::new(),
            unsolved: None,
            unsolved_t: Vec::new(),
            conditions: Vec::new(),
            bound: None,
            placed_fields: HashMap::new(),
            xyz_read: HashMap::new(),
            integrand_of: None,
            file: Rc::from(""),
            scope: Scope::STD,
            files: vec![SourceFile::new(Rc::from("std"), "std".into(), None, ".".into())],
            running: Vec::new(),
            origin: HashMap::new(),
            param_keys: HashMap::new(),
            depth: 0,
            inputs: Vec::new(),
            import_depth: 0,
            overrides: Vec::new(),
            baked: Default::default(),
            materialized: Default::default(),
            undefined_at_load: Default::default(),
            shape_of: HashMap::new(),
            form_at: String::new(),
            wild_written: None,
            form_name: String::new(),
            building: None,
            own_frame: false,
            source: None,
            coordinates_of: None,
            states: false,
            state_t: None,
            coords: Vec::new(),
            coordinate_params: Vec::new(),
        }
    }

    /// Load a file with the packages of `registry`, and parameter
    /// overrides applied to load-time computations (parameters stay
    /// symbolic in fields, so derivatives still work). The standard
    /// distribution's loaders (`ink::load_file`) pass `ink::registry()`.
    pub fn load_file_in(registry: &'static Registry, path: &Path, overrides: &[(String, Given)]) -> R<Doc> {
        Self::load_file_reads(registry, path, overrides).map_err(|(e, _)| e)
    }

    /// `load_file_in`, and for a load that fails, what it read before
    /// failing (`files_read`, the file it could not read among them): a
    /// failed load rests on those as a load rests on its `files_read`, so
    /// it fails the same way while each reads as it did.
    pub fn load_file_reads(registry: &'static Registry, path: &Path, overrides: &[(String, Given)]) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let src = std::fs::read_to_string(path).map_err(|e| (InkError::msg(format!("cannot read {}: {e}", shown(path))), vec![(path.to_path_buf(), None)]))?;
        let base = path.parent().unwrap_or(Path::new("."));
        Self::load_named(registry, &src, &shown(path), path, base, overrides)
    }

    /// Note a file the load read (a source it ran, std included, or data
    /// it read): its name as shown (`shown`: in the catalog, `parts/...`),
    /// where it is (std: none, it is the language's) and its content.
    fn read(&mut self, name: &str, path: Option<&Path>, bytes: &[u8]) {
        let mut fp = graph::Fingerprint::default();
        fp.str(name).word(bytes.len() as u64).raw(bytes);
        self.sources_read.push(fp.value());
        self.files_read.extend(path.map(|p| (p.to_path_buf(), Some(Doc::content_key(bytes)))));
    }

    /// Note a path the load used for a file it read already, maybe by
    /// another path (a symlink, `..`): the load rests on what each path it
    /// used names, so a path retargeted later is seen.
    pub(crate) fn read_again(&mut self, path: &Path) {
        if !self.files_read.iter().any(|(f, _)| f == path) {
            self.files_read.push((path.to_path_buf(), std::fs::read(path).ok().map(|b| Doc::content_key(&b))));
        }
    }

    /// A file's content as `files_read` keys it.
    pub fn content_key(bytes: &[u8]) -> u128 {
        graph::Fingerprint::default().word(bytes.len() as u64).raw(bytes).value()
    }

    /// The files on disk the load depended on: each it read (a source it
    /// ran or data it read) with its content's key (`content_key`), and
    /// each it found absent (None: a `parts/...` import not beside its
    /// file, read from the catalog; a file it could not read, where the
    /// load failed: `load_file_reads`). With the code, they are all a load
    /// rests on: where each reads as it did, the load is the same.
    pub fn files_read(&self) -> &[(std::path::PathBuf, Option<u128>)] {
        &self.files_read
    }

    /// The fingerprint of everything the model was read from: the name and
    /// content of its own file, of every file it imports, transitively
    /// (std's text included), and of every file it reads (`read-json`,
    /// `import-stl`), as a set. A file in the catalog is named by its path
    /// there, so the same catalog gives the same fingerprint wherever it
    /// is. With the kernel's code, it is all a load of the model rests on.
    pub fn sources(&self) -> u128 {
        let mut read = self.sources_read.clone();
        read.sort_unstable();
        read.dedup();
        let mut fp = graph::Fingerprint::default();
        fp.word(read.len() as u64);
        read.iter().for_each(|k| _ = fp.key(*k));
        fp.value()
    }

    /// Load source text with the packages of `registry`.
    pub fn load_str_in(registry: &'static Registry, src: &str, name: &str, base: &Path, overrides: &[(String, Given)]) -> R<Doc> {
        Self::load_named(registry, src, name, Path::new(name), base, overrides).map_err(|(e, _)| e)
    }

    /// Load source text read from `path`, named `name` to the user; when
    /// it fails, with what it read before failing (`load_file_reads`).
    fn load_named(
        registry: &'static Registry,
        src: &str,
        name: &str,
        path: &Path,
        base: &Path,
        overrides: &[(String, Given)],
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let origin = Origin { src: src.to_string(), name: name.to_string(), path: path.to_path_buf(), base: base.to_path_buf() };
        Self::load_origin(registry, origin, overrides, None)
    }

    /// The file run loaded again, the arguments `target`'s instances pass
    /// (a top-level assembly of the file) built as coordinates (`coords`):
    /// what a package's answer differentiates placements along (a dynamics
    /// package built on the kernel's, say). Every other
    /// load builds them as written, so nothing else reads the coordinates'
    /// params, and no package adds what it adds once a file is read
    /// (Package::finish).
    pub fn with_coordinates(&self, target: &str) -> R<Doc> {
        self.reloader().ok_or_else(|| InkError::msg("only a document loaded from a source can be loaded again"))?.with_coordinates(target)
    }

    /// Whether it is a load `with_states` (Doc::with_states): the mechanism
    /// at any state, its moving coordinates' values, rates and accelerations
    /// params.
    pub fn states(&self) -> bool {
        self.states
    }

    /// How to load the file run again (`with_coordinates`), apart from the
    /// document: for a provider, which works without it.
    pub fn reloader(&self) -> Option<Reload> {
        Some(Reload { registry: self.registry, origin: self.source.clone()?, overrides: self.overrides.clone() })
    }

    fn load_origin(
        registry: &'static Registry,
        origin: Origin,
        overrides: &[(String, Given)],
        coordinates_of: Option<(String, bool)>,
    ) -> Result<Doc, (InkError, Vec<(std::path::PathBuf, Option<u128>)>)> {
        let mut d = Doc::empty(registry);
        d.overrides = overrides.to_vec();
        d.file = Rc::from(origin.name.as_str());
        (d.coordinates_of, d.states) = coordinates_of.map_or((None, false), |(t, s)| (Some(t), s));
        // Known while reading: a provider made then (a simulation) loads it
        // again.
        d.source = Some(origin.clone());
        let mut read = d.run_source_as(&origin.src, &origin.name, &origin.path, &origin.base, None);
        if read.is_ok() {
            d.rate_joints();
        }
        // What reads names after loading (a region on the command line)
        // reads the file's.
        d.scope = Scope { file: 1 };
        if read.is_ok() && d.coordinates_of.is_none() {
            read = registry.packages.iter().try_for_each(|p| {
                let mut r = Reader::new(&mut d);
                p.finish(&mut r)
            });
        }
        // A solve that does not solve with every other one: the error, once
        // all are read (Doc::solve_at_load); one depending on t, if what
        // rests on it was used at load time, where t = 0.
        let unsolved = d.unsolved.take().or_else(|| {
            let mut used: Vec<(usize, String)> = d.materialized.borrow().iter().map(|(&k, site)| (k, site.clone())).collect();
            used.sort();
            used.into_iter().find_map(|(k, site)| {
                let (_, e) = d.unsolved_t.iter().find(|(i, _)| d.solves_under(&[k]).contains(i))?;
                Some(InkError {
                    hint: Some(format!(
                        "'{}' is used at load time at {site}, where t = 0. Use it as a field (solved again at each moment), or make the solve close at t = 0",
                        d.params[k].name
                    )),
                    ..e.clone()
                })
            })
        });
        // The directions at the params' values (Doc::unit), once every
        // solve is read; those depending on t at each moment a view is
        // made at.
        let read = read.and_then(|()| if unsolved.is_none() { d.check_conditions(&d.load_ctx(), None, &[], true) } else { Ok(()) });
        let failed = match (read, unsolved) {
            (Ok(()), None) => return Ok(d),
            (Ok(()), Some(u)) => u,
            (Err(e), None) => e,
            (Err(mut e), Some(u)) => {
                let note = format!("before this, at {}:{}: {}. What reads its unknowns after it gets no solution", u.file, u.line, u.message);
                e.hint = Some(e.hint.map_or(note.clone(), |h| format!("{h}; {note}")));
                e
            }
        };
        Err((failed, d.files_read))
    }

    /// Context for load-time evaluation: parameters declared so far,
    /// overridden, with solved params re-solved for the overrides.
    fn load_ctx(&self) -> Ctx {
        let mut params = self.load_params();
        for i in (0..self.solves.len()).filter(|&i| self.provided_always(i)) {
            if self.provide_one(i, &mut params, 0.0).is_err() {
                break;
            }
        }
        Ctx::new(0.0, params)
    }

    /// A load-time calculation needs only the providers its inputs read,
    /// transitively. A geometry query unrelated to density must not measure
    /// every published-mass body read so far.
    fn load_ctx_for(&self, roots: &[NodeId]) -> Ctx {
        let mut params = self.load_params();
        for i in self.solves_read(roots) {
            if self.provide_one(i, &mut params, 0.0).is_err() {
                break;
            }
        }
        Ctx::new(0.0, params)
    }

    /// The params' load-time values before solving: defaults and overrides.
    fn load_params(&self) -> Vec<f64> {
        let mut params: Vec<f64> = self.params.iter().map(|p| p.default).collect();
        for (k, v) in &self.overrides {
            if let Some(i) = self.params.iter().position(|p| p.name == *k) {
                params[i] = v.value;
            }
        }
        params
    }

    /// What a file's run may make besides values and functions: params,
    /// bodies, assemblies, claims, extension forms and solves, or a span
    /// of t. A run that made none of them defined what it defines
    /// from its source alone.
    fn made(&self) -> (usize, usize, usize, usize, usize, usize, usize, usize) {
        (self.params.len(), self.bodies.len(), self.assemblies.len(), self.asserts.len(), self.ext.len(), self.solves.len(), self.arg_defs.len(), self.file_spans.len())
    }

    /// Imports library `canon` under namespace `ns` from its run under
    /// another (`libraries`): a file under `ns` that reads its names as that
    /// run's definitions, which are those a run under `ns` would make. A
    /// part family's file read by each of its parts so runs once (each part
    /// kept its own copy of every function of its family, 128 KB a part).
    /// False: no such run.
    fn library_under(&mut self, canon: &str, ns: Option<&str>, path: &Path) -> bool {
        let (Some(ns), Some((file, _))) = (ns, self.libraries.get(canon).cloned()) else { return false };
        let src = &self.files[file as usize];
        let mut f = SourceFile::new(src.name.clone(), canon.to_string(), Some(ns.to_string()), path.parent().map_or(src.dir.clone(), Path::to_path_buf));
        f.bindings = src.bindings.clone();
        // Defined there, as a run here would have (a part's own name
        // cannot shadow its family's).
        for (n, k) in &f.bindings {
            if let Some(&at) = self.origin.get(k) {
                self.origin.insert(format!("{ns}.{n}"), at);
            }
        }
        self.files.push(f);
        true
    }

    /// Run a file's forms (read from `path`, named `name` to the user); with
    /// a namespace, every top-level name it defines is `ns.name` in the
    /// document (Doc::public), and `name` as the file reads it.
    fn run_source_as(&mut self, src: &str, name: &str, path: &Path, base: &Path, ns: Option<&str>) -> R<()> {
        let id = if name == "std" {
            0
        } else {
            let path = path.canonicalize().map_or(path.display().to_string(), |p| p.display().to_string());
            let mut f = SourceFile::new(Rc::from(name), path, ns.map(str::to_string), base.to_path_buf());
            f.given = std::mem::take(&mut self.giving);
            self.files.push(f);
            self.files.len() as u32 - 1
        };
        self.read(name, (name != "std").then_some(path), src.as_bytes());
        let saved = std::mem::replace(&mut self.file, Rc::from(name));
        let saved_scope = std::mem::replace(&mut self.scope, Scope { file: id });
        self.running.push(id);
        let forms = parse(src).map_err(|e| InkError { message: e.msg, file: name.to_string(), line: e.line, col: e.col, hint: None, called_from: Vec::new() })?;
        // Where the file defines each value, to say so of a read before it.
        let ahead = forms
            .iter()
            .filter_map(|f| match f.list() {
                Some([h, n, ..]) if matches!(h.sym(), Some("def" | "fn" | "param")) => Some((n.sym()?.to_string(), f.line)),
                _ => None,
            })
            .collect();
        self.files[id as usize].ahead = ahead;
        self.files[id as usize].declared = top_names(&forms, self.registry);
        let saved_ns = std::mem::replace(&mut self.ns, ns.map(str::to_string));
        let r = self.run_forms(&forms, base);
        self.ns = saved_ns;

        self.file = saved;
        self.scope = saved_scope;
        self.running.pop();
        r
    }

    /// Read a file's top-level forms in order, each in the names defined
    /// before it. A `for` is expanded when it is reached, so its list reads
    /// the file's names too, and the forms it makes are read there, under
    /// what its list rested on (Doc::under), a `for` among them expanded as
    /// it is reached in turn.
    fn run_forms(&mut self, forms: &[Sexp], base: &Path) -> R<()> {
        for f in forms {
            if let Some((var, list, body)) = author::for_parts(f) {
                let (var, values) = self.at_form(f, |d| d.for_values(var, list))?;
                for (val, text) in values {
                    let made = author::for_made(body, &var, &val, &text);
                    self.run_forms(&made, base)?;
                }
                continue;
            }
            self.at_form(f, |d| d.top(f, base))?;
        }
        Ok(())
    }

    /// Run a top-level form, or a `for`'s list: what it makes a number at
    /// load time is used there (`form_at`: Doc::materialized), and the
    /// undefined cases those numbers pass through are reported, never
    /// silent (`load_eval`).
    fn at_form<T>(&mut self, f: &Sexp, run: impl FnOnce(&mut Self) -> R<T>) -> R<T> {
        let head: Vec<String> = f.list().unwrap_or_default().iter().take(2).filter_map(|x| x.sym().map(str::to_string)).collect();
        self.form_at = format!("{}:{} ({} ...)", self.file, f.line, head.join(" "));
        self.form_name = head.last().cloned().unwrap_or_default();
        let r = run(self)?;
        for what in self.undefined_at_load.take() {
            let w = format!("{}:{}: {what} in a number used at load time. An undefined case gives 0 (see `ink spec`)", self.file, f.line);
            if !self.load_warnings.contains(&w) {
                self.load_warnings.push(w);
            }
        }
        Ok(r)
    }

    /// `f` (no x, y, z) as a number at load time, its undefined cases kept
    /// for `at_form` to report.
    fn load_eval(&self, f: NodeId) -> f64 {
        let (tape, ctx) = (Tape::new(&self.g, f), self.load_ctx_for(&[f]));
        let cases = ink_core::eval::undefined_cases(&tape, &ctx, &[[0.0; 3]]);
        self.undefined_at_load.borrow_mut().extend(cases.iter().map(|u| u.what));
        eval_point(&tape, &ctx, [0.0; 3])
    }

    fn err(&self, at: &Sexp, msg: impl Into<String>) -> InkError {
        InkError { message: msg.into(), file: self.file.to_string(), line: at.line, col: at.col, hint: None, called_from: Vec::new() }
    }

    fn err_hint(&self, at: &Sexp, msg: impl Into<String>, hint: Option<String>) -> InkError {
        InkError { hint, ..self.err(at, msg) }
    }

    // ---------------------------------------------------------------- lookup

    /// Whether an assert is written in the file loaded, not in one it imports.
    /// The motion what is written in `file` (`AssertDef::file`,
    /// `Assembly::file`) is checked over: that file's `(t-span A B)`, an
    /// import's as much as the run file's; where it states none, the run
    /// file's. None: neither states one.
    pub fn span_of(&self, file: &str) -> Option<(f64, f64)> {
        self.file_spans.get(file).copied().or(self.t_span)
    }

    /// The hull of every span the model's files state (the file run's
    /// and its imports'): what `interfere` asks a pair reading t over.
    pub fn spans_hull(&self) -> Option<(f64, f64)> {
        self.file_spans.values().copied().chain(self.t_span).reduce(|(a, b), (c, d)| (a.min(c), b.max(d)))
    }

    pub fn is_own(&self, a: &AssertDef) -> bool {
        self.is_own_file(&a.file)
    }

    /// Whether `file` (where something is written, as `AssertDef::file`)
    /// is the file loaded, not one it imports.
    pub fn is_own_file(&self, file: &str) -> bool {
        file == &*self.file
    }

    /// The file a param, def or fn named `name` is written in. None: a
    /// name of std's, or none.
    pub fn written_in(&self, name: &str) -> Option<&str> {
        self.origin.get(name).map(|&(f, _)| &*self.files[f as usize].name)
    }

    /// The files the file loaded imports, directly or through others, in
    /// the order they were read; std is not one.
    pub fn imported_files(&self) -> Vec<&str> {
        let mut out: Vec<&str> = Vec::new();
        for f in self.files.iter().skip(2) {
            if !out.contains(&&*f.name) {
                out.push(&f.name);
            }
        }
        out
    }

    pub fn param_index(&self, name: &str) -> Option<usize> {
        self.param_keys.get(&*self.scoped(name)).copied()
    }

    pub fn body(&self, name: &str) -> Option<&Body> {
        self.bodies.iter().find(|b| b.name == name)
    }

    pub fn assembly(&self, name: &str) -> Option<&Assembly> {
        self.assemblies.iter().find(|a| a.name == name)
    }

    pub fn def_index_of(&self, name: &str) -> Option<usize> {
        self.def_index.get(&*self.scoped(name)).copied()
    }

    pub fn def(&self, name: &str) -> Option<&Value> {
        self.def_index_of(name).map(|i| &self.defs[i].1)
    }

    /// Whether param k is ranged where claims are answered: a tolerance,
    /// or a (ranged lo hi).
    pub fn has_range(&self, k: usize) -> bool {
        self.params[k].tol.is_some() || self.ranged.iter().any(|r| r.0 == k)
    }

    /// What to do where an answer needs one value of each of the ranged
    /// params `ks` it reads, over `spans` (the view's): a tolerance is
    /// checked at one size of it, each named with the ends of its band and
    /// the flags that pin it there (which end is worst is not proven, so
    /// each answer is that size's); a (ranged lo hi) is a range of the
    /// question, which no claim takes at one value, so its answer is
    /// looked at at one (eval, render, a viewer, ink check's conditions).
    pub fn at_one_value(&self, ks: &[usize], spans: &ink_core::eval::ParamSpans) -> String {
        let ranged: Vec<&str> = ks.iter().filter(|&&k| self.ranged.iter().any(|r| r.0 == k)).map(|&k| self.params[k].name.as_str()).collect();
        if !ranged.is_empty() {
            return format!(
                "'{}' is (ranged lo hi), a range of the question, so no claim takes it at one value. See one with ink eval or ink render --param NAME=V",
                ranged.join("', '")
            );
        }
        // An end as it is written: the band rounded outward is within a
        // float or two of it.
        let num = |x: f64| format!("{}", (x * 1e9).round() / 1e9);
        let each: Vec<String> = ks
            .iter()
            .map(|&k| {
                let name = &self.params[k].name;
                match spans.get(k as u32) {
                    Some((a, b)) if a < b => format!("--param {name}={a}:{a} or --param {name}={b}:{b} (its ends)", a = num(a), b = num(b)),
                    _ => format!("--param {name}=V:V"),
                }
            })
            .collect();
        match each.is_empty() {
            true => "check at fixed params with --nominal or --param NAME=V:V".into(),
            false => format!(
                "check it at one size of each: {}. Or use --nominal: each param at its value, no tolerance. Which size is worst is not proven: each answer is that size's",
                each.join("; ")
            ),
        }
    }

    /// The (ranged lo hi) named `name`, what --param NAME reaches: the def
    /// that made it, or its generated name ("t-slot-reaches 0"). Its names
    /// are unique, and none is a declared param's (a def's name is not; a
    /// generated one has a space), so load_params finds it by name too.
    pub fn ranged_named(&self, name: &str) -> Option<(usize, (f64, f64))> {
        self.ranged.iter().find(|r| self.params[r.0].name == name).copied()
    }

    /// The (ranged lo hi) a value is, if it is one itself.
    fn ranged_of(&self, v: &Value) -> Option<(usize, (f64, f64))> {
        let &Value::Field(f, _) = v else { return None };
        let graph::Node::Var(Var::Param(k)) = self.g.nodes[f as usize] else { return None };
        self.ranged.iter().find(|r| r.0 == k as usize).copied()
    }

    fn value_names(&self) -> Vec<String> {
        let mut v: Vec<String> = self.params.iter().map(|p| p.name.clone()).collect();
        v.extend(self.defs.iter().map(|d| d.0.clone()));
        v.extend(self.fns.keys().filter(|n| !n.contains(':')).cloned());
        v
    }

    /// A name a param is set by on the command line, `--param name=v`: a
    /// param's, written or generated (a `for`'s, a solve's unknown), or an
    /// import's namespace, which prefixes its params' names with `.`. It is
    /// one plain identifier (`param_word`), so every param's name, `ns.p`
    /// included, is set by --param as it is printed, unquoted, in any shell.
    pub(crate) fn check_param_name(&self, at: &Sexp, name: &str, what: &str) -> R<()> {
        if crate::command::param_word(name) {
            return Ok(());
        }
        Err(self.err_hint(
            at,
            format!(
                "{what} {name:?} is not a plain name. A param is set on the command line (--param name=v), so its name is letters, digits, '-' and '_', starting with a letter"
            ),
            Some("rename it, e.g. with '-' for a space".into()),
        ))
    }

    fn check_new_value_name(&self, at: &Sexp, name: &str) -> R<()> {
        if is_reserved(name) {
            return Err(self.err_hint(at, format!("'{name}' is a built-in name and cannot be redefined"), Some(reserved_hint())));
        }
        // std's names are stored under "std:<name>" (and under the name
        // until a file defines it): outside std a name is defined once, and
        // may shadow std's, which std itself keeps reading (Doc::scoped).
        let already = match self.scope.file {
            0 => self.files[0].bindings.contains_key(name).then(|| ("std".to_string(), 0)),
            _ => self
                .origin
                .get(&*self.scoped(name))
                .filter(|&&(f, _)| f == self.scope.file)
                .or_else(|| self.origin.get(&self.public(name)).filter(|_| !Self::std_name(name)))
                .map(|&(f, line)| (format!("{}:{line}", self.files[f as usize].name), f)),
        };
        if let Some((at_, f)) = already {
            // An imported param: given at its import.
            let hint = if f != self.scope.file && self.params.iter().any(|p| p.name == name) {
                format!("an importer gives an imported param its value at the import, (import \"{}\" ({name} V))", self.files[f as usize].name)
            } else {
                "names are immutable; pick a new name".into()
            };
            return Err(self.err_hint(at, format!("'{name}' is already defined, at {at_}"), Some(hint)));
        }
        Ok(())
    }

    /// Expression `e` as a number known at load time, with what says its
    /// quantity (units::Said); `none`: the error where it is no number.
    fn constant(&mut self, e: &Sexp, none: &str) -> R<(f64, Option<Rc<units::Said>>)> {
        let v = self.eval(e, &mut Vec::new())?;
        let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
        let f = self.as_field(v, e)?;
        let x = self.g.as_const(f).ok_or_else(|| self.err(e, none))?;
        Ok((x, said))
    }

    /// A param's `(tol lo hi)`: deviations from a value said by `said`;
    /// with the units its ends are said in, each where.
    fn tol(&mut self, t: &Sexp, said: Option<&units::Said>, usage: &str) -> R<((f64, f64), Vec<(Rc<units::Said>, u32, u32)>)> {
        let Some([_, lo, hi]) = t.list().filter(|l| l.first().and_then(|h| h.sym()) == Some("tol")) else {
            return Err(self.err_hint(t, "expected (tol lo hi)", Some(usage.into())));
        };
        let mut units = Vec::new();
        let mut num = |d: &mut Doc, e: &Sexp| -> R<f64> {
            let v = d.eval(e, &mut Vec::new())?;
            d.fits_change(&v, e, "a tolerance's deviation")?;
            if let Some(s) = said {
                d.fits(&v, e, s.unit.quantity, "a tolerance's deviation")?;
            }
            if let Value::Field(_, Some(s)) = &v {
                units.push((s.clone(), e.line, e.col));
            }
            let f = d.as_field(v, e)?;
            d.g.as_const(f).filter(|x| x.is_finite()).ok_or_else(|| d.err(e, "a tolerance's ends must be constant finite numbers"))
        };
        let (lo, hi) = (num(self, lo)?, num(self, hi)?);
        if lo > hi {
            return Err(self.err_hint(t, format!("tol needs lo <= hi, got {lo} > {hi}"), Some(usage.into())));
        }
        Ok(((lo, hi), units))
    }

    /// An import's values, `(name value [(tol lo hi)])` each, read in the
    /// importer's names, with the clause each is written in.
    fn import_values(&mut self, values: &[Sexp], usage: &str) -> R<Vec<(ImportValue, Sexp)>> {
        let mut given: Vec<(ImportValue, Sexp)> = Vec::new();
        for c in values {
            let Some([name, value, tol @ ..]) = c.list().filter(|l| l.len() <= 3) else {
                return Err(self.err_hint(c, "expected (name value [(tol lo hi)]): a value for a param of the imported file", Some(usage.into())));
            };
            let n = self.name_of(name)?;
            if given.iter().any(|(g, _)| g.name == n) {
                return Err(self.err(name, format!("'{n}' is given twice")));
            }
            let (x, said) = self.constant(value, "an import's value is a constant number, as a param default is")?;
            let (tol, mut units) = match tol.first() {
                Some(t) => self.tol(t, said.as_deref(), usage).map(|(b, u)| (Some(b), u))?,
                None => (None, Vec::new()),
            };
            units.extend(said.clone().map(|s| (s, value.line, value.col)));
            given.push((ImportValue { name: n, value: x, said, tol, file: self.file.clone(), line: value.line, units }, c.clone()));
        }
        Ok(given)
    }

    /// Import value `g` (clause `at`) sets a param written in file `i`
    /// (imported as `path`): refused where it names no `(param ...)` of
    /// that file, or, where the file was already run (first imported at
    /// `first`), where that run gave the param another value.
    fn import_sets(&self, i: usize, g: &ImportValue, at: &Sexp, path: &str, first: Option<&str>) -> R<()> {
        let (f, n) = (&self.files[i], &g.name);
        let own = |key: &String| self.origin.get(key).is_some_and(|o| o.0 as usize == i);
        let Some(key) = f.bindings.get(n) else {
            let params = f.bindings.iter().filter(|(_, k)| self.param_keys.contains_key(*k) && own(k)).map(|(n, _)| n.as_str());
            return Err(self.err_hint(at, format!("{} has no param '{n}'", f.name), closest(n, params)));
        };
        let not = |what: &str| Err(self.err(at, format!("'{n}' is {what}: an import sets params")));
        let k = match self.param_keys.get(key) {
            Some(_) if !own(key) => {
                let from = self.origin.get(key).map_or("another file", |o| &*self.files[o.0 as usize].name);
                return Err(self.err(at, format!("{n} is {from}'s: import {from} with it, before {}", f.name)));
            }
            Some(&k) => k,
            None if self.ranged_named(&f.ns.as_ref().map_or(n.clone(), |ns| format!("{ns}.{n}"))).is_some() => return not("a (ranged ...), a range of the question"),
            None if self.fns.contains_key(key) => return not("a fn, not a param"),
            None => return not("a def, not a param"),
        };
        let p = &self.params[k];
        let Some(first) = first else {
            // Not set by its param form: a package form's (a solve's unknown).
            return if self.set_at.contains_key(&k) { Ok(()) } else { not("made by a package form, not by (param ...)") };
        };
        if let Some(s) = self.param_said.get(&k) {
            g.fits(s.unit.quantity, &format!("{}'s {n}", f.name))?;
        }
        let asked = &at.list().expect("a clause")[1];
        let copy = format!("import a second copy under a namespace, (import \"{path}\" as NS ({n} {asked}))");
        let unit = self.param_said.get(&k).map_or("", |s| s.unit.quantity.unit());
        match self.set_at.get(&k) {
            Some(_) if p.default == g.value && g.tol.is_none_or(|t| p.tol == Some(t)) => Ok(()),
            Some((_, by)) => Err(self.err_hint(
                at,
                format!("{}'s {n} is already {}{unit}, set by {by}; this import asks {asked}", f.name, p.default),
                Some(format!("one load gives a file's param one value: {copy}")),
            )),
            None => Err(self.err_hint(
                at,
                format!("{} was read at {first} with {n} {}{unit}", f.name, p.default),
                Some(format!("give the value at {}'s first import ({first}), or {copy}", f.name)),
            )),
        }
    }

    fn std_name(name: &str) -> bool {
        static NAMES: std::sync::OnceLock<HashSet<String>> = std::sync::OnceLock::new();
        NAMES
            .get_or_init(|| {
                std_forms()
                    .iter()
                    .filter_map(|f| match f.list() {
                        Some([h, n, ..]) if matches!(h.sym(), Some("fn" | "def")) => n.sym().map(str::to_string),
                        _ => None,
                    })
                    .collect()
            })
            .contains(name)
    }

    /// Keep every file's definitions under their own keys, alongside the
    /// public names. Importing a file copies its visible bindings; later
    /// shadows in another file cannot change those keys. `name` as the
    /// file writes it; the public name is Doc::public's.
    fn defined(&mut self, name: &str, line: u32) -> String {
        let key = if self.scope.file == 0 { format!("std:{name}") } else { format!("file:{}:{name}", self.scope.file) };
        self.files[self.scope.file as usize].bindings.insert(name.to_string(), key.clone());
        if self.scope.file != 0 {
            let public = self.public(name);
            self.origin.insert(key.clone(), (self.scope.file, line));
            self.origin.insert(public.clone(), (self.scope.file, line));
            self.def_index.remove(&public);
            self.fns.remove(&public);
        }
        key
    }

    /// A name a top-level form of the file being read defines, as the
    /// document names it: in a file imported under a namespace, `ns.name`
    /// (`p.plate`, `p.q.plate` through a namespace inside it); else as written.
    fn public(&self, name: &str) -> String {
        match &self.files[self.scope.file as usize].ns {
            Some(ns) => format!("{ns}.{name}"),
            None => name.to_string(),
        }
    }

    /// A body or assembly, or a path from one (`arm/tip`), as the file
    /// being read writes it, as the document names it: its first segment
    /// read in the file's names (SourceFile::objects).
    fn resolve(&self, path: &str) -> String {
        let (first, rest) = path.split_once('/').map_or((path, None), |(a, b)| (a, Some(b)));
        match self.files[self.scope.file as usize].objects.get(first) {
            Some(n) => rest.map_or(n.clone(), |r| format!("{n}/{r}")),
            None => path.to_string(),
        }
    }

    /// A name (a body, an assembly, a path from one) written at `e`, as
    /// the document names it (Doc::resolve).
    fn path_of(&self, e: &Sexp) -> R<String> {
        Ok(self.resolve(&self.name_of(e)?))
    }

    fn scoped<'a>(&'a self, name: &'a str) -> std::borrow::Cow<'a, str> {
        if name.contains(':') || CONSTANTS.contains(&name) {
            return std::borrow::Cow::Borrowed(name);
        }
        if let Some(key) = self.files[self.scope.file as usize].bindings.get(name) { std::borrow::Cow::Borrowed(key) } else { std::borrow::Cow::Owned(format!("std:{name}")) }
    }

    /// A read of a name the file being read defines further on, not yet
    /// defined: an error, at a top-level form or in a function (a file's
    /// name means its own definition throughout, never std's before it).
    fn later(&self, s: &str, at: &Sexp) -> Option<InkError> {
        let file = &self.files[self.scope.file as usize];
        let line = file.ahead.get(s).filter(|_| self.scope.file != 0 && self.origin.get(&*self.scoped(s)).is_none_or(|&(f, _)| f != self.scope.file))?;
        Some(self.err_hint(at, format!("unknown name '{s}': it is defined at line {line}, after this"), Some("a file is read in order, so define a name before it is used".into())))
    }

    /// The directory of the file being read: the paths it names are
    /// relative to it (a function's, in its body).
    fn dir(&self) -> &Path {
        &self.files[self.scope.file as usize].dir
    }

    /// Why a name defined outside std is not visible here (Doc::scoped).
    fn hidden_hint(&self, name: &str) -> Option<String> {
        let &(file, line) = self.origin.get(name).filter(|_| self.scope.file != 0 && !self.files[self.scope.file as usize].bindings.contains_key(name))?;
        Some(format!(
            "'{name}' is defined at {}:{line}, in a file this one does not import. A file sees only its own names, its imports' and std's",
            self.files[file as usize].name
        ))
    }

    /// Where a form headed by a top-level form's name, written inside
    /// `inside`, belongs: (t-span 0 2) in an assembly.
    fn top_form_hint(&self, c: &Sexp, inside: &str) -> Option<String> {
        let head = c.list()?.first()?.sym()?;
        self.registry
            .top_forms()
            .split(", ")
            .any(|f| f == head)
            .then(|| format!("({head} ...) is a top-level form. Write it at the file's top level, not inside {inside} (`ink spec --section \"Top-level forms\"`)"))
    }

    fn check_new_object_name(&self, at: &Sexp, name: &str) -> R<()> {
        if name.contains('/') {
            return Err(self.err(at, format!("'{name}': body/assembly names cannot contain '/'")));
        }
        if self.body(name).is_some() || self.assembly(name).is_some() {
            return Err(self.err(at, format!("a body or assembly named '{name}' already exists")));
        }
        Ok(())
    }

    // ---------------------------------------------------------------- top level

    fn top(&mut self, form: &Sexp, base: &Path) -> R<()> {
        let items = form.list().ok_or_else(|| self.err(form, format!("top-level forms must be lists: one of {}", self.registry.top_forms())))?;
        let head = items.first().and_then(|h| h.sym()).ok_or_else(|| self.err(form, "expected a form name"))?;
        let args = &items[1..];
        // Any form but a claim may define what shape-of reads.
        if head != "assert" {
            self.shape_of.clear();
        }
        match head {
            "import" => {
                let usage = "usage: (import \"std\"), (import \"path.ink\") or (import \"path.ink\" as ns): with a namespace \
                             its names read ns.name; then (name value [(tol lo hi)]) ... sets the file's params";
                // A plain import inside a namespaced file joins its namespace,
                // so a part's family stays out of the importer's names.
                // `as`: the namespace its names read under here (`as`), and
                // under in the document (`ns`: this file's, then it).
                let (target, ns, as_, values) = match args {
                    [t, a, n, values @ ..] if a.sym() == Some("as") => {
                        let at = n;
                        let n = self.name_of(n)?;
                        self.check_param_name(at, &n, "the namespace")?;
                        (t, Some(self.public(&n)), Some(n), values)
                    }
                    [t, values @ ..] => (t, self.ns.clone(), None, values),
                    [] => return Err(self.err_hint(form, "bad import", Some(usage.into()))),
                };
                let SK::Str(p) = &target.k else { return Err(self.err(target, "import path must be a string")) };
                if p == "std" {
                    if args.len() > 1 {
                        return Err(self.err(target, "std is imported as it is, (import \"std\"): your names shadow its names"));
                    }
                    if self.imported.insert("std".into(), String::new()).is_none() {
                        self.run_source_as(STD_SRC, "std", Path::new("std"), base, None)?;
                    }
                } else {
                    let given = self.import_values(values, usage)?;
                    let (path, absent) = import_path(base, p);
                    self.files_read.extend(absent.into_iter().map(|f| (f, None)));
                    let canon = path.canonicalize().unwrap_or(path.clone()).display().to_string();
                    if let Some(i) = self.running.iter().position(|&r| self.files[r as usize].path == canon) {
                        let chain: Vec<&str> = self.running[i..].iter().map(|&r| &*self.files[r as usize].name).chain([&*self.files[self.running[i] as usize].name]).collect();
                        return Err(self.err_hint(
                            target,
                            format!("import cycle: {}", chain.join(" -> ")),
                            Some("a file cannot import itself, directly or through the files it imports".into()),
                        ));
                    }
                    let key = ns.as_ref().map_or(canon.clone(), |n| format!("{canon} as {n}"));
                    let first = self.imported.get(&key).cloned();
                    if first.is_none() {
                        self.imported.insert(key, format!("{}:{}", self.file, form.line));
                    }
                    if first.is_none() && !self.library_under(&canon, ns.as_deref(), &path) {
                        let src = std::fs::read_to_string(&path).inspect_err(|_| self.files_read.push((path.clone(), None))).map_err(|e| {
                            let dirs = import_dirs();
                            let hint = std::iter::once(path.clone()).chain(dirs.iter().map(|d| d.join(p))).find_map(|f| sibling(&f, p)).or_else(|| {
                                let dirs: Vec<String> = dirs.iter().map(|d| d.display().to_string()).collect();
                                let dirs = if dirs.is_empty() { "none: INK_PATH is unset".into() } else { dirs.join(", ") };
                                Path::new(p).is_relative().then(|| format!("an import is read beside its file, else from INK_PATH's directories ({dirs})"))
                            });
                            self.err_hint(target, format!("cannot read {}: {e}", shown(&path)), hint)
                        })?;
                        let dir = path.parent().unwrap_or(base).to_path_buf();
                        self.import_depth += 1;
                        let (before, varied) = (self.made(), self.varies);
                        self.giving = given.iter().map(|(v, _)| v.clone()).collect();
                        self.run_source_as(&src, &shown(&path), &path, &dir, ns.as_deref())?;
                        self.import_depth -= 1;
                        // A run that made only values and functions, of
                        // nothing that differs by namespace, is the same
                        // run under any namespace.
                        if let Some(n) = &ns
                            && self.made() == before
                            && self.varies == varied
                            && let Some(i) = self.files.iter().rposition(|f| f.path == canon && f.ns == ns)
                        {
                            self.libraries.entry(canon.clone()).or_insert((i as u32, n.clone()));
                        }
                    } else {
                        self.read_again(&path);
                    }
                    if !self.libraries.contains_key(&canon) {
                        self.varies += 1;
                    }
                    // The file sees the names of what it imports (Scope),
                    // under its namespace here (`as`).
                    if let Some(i) = self.files.iter().rposition(|f| f.path == canon && f.ns == ns) {
                        for (v, at) in &given {
                            self.import_sets(i, v, at, p, first.as_deref())?;
                        }
                        let named = |n: &String| as_.as_ref().map_or(n.clone(), |a| format!("{a}.{n}"));
                        let bindings: Vec<(String, String)> = self.files[i].bindings.iter().map(|(n, k)| (named(n), k.clone())).collect();
                        let objects: Vec<(String, String)> = self.files[i].objects.iter().map(|(n, o)| (named(n), o.clone())).collect();
                        let me = self.scope.file as usize;
                        for (name, key) in bindings {
                            let own = self.files[me].bindings.get(&name).and_then(|k| self.origin.get(k)).is_some_and(|&(f, _)| f == self.scope.file);
                            if !own {
                                self.files[me].bindings.insert(name, key);
                            }
                        }
                        for (name, o) in objects {
                            if self.files[me].objects.get(&name) != Some(&self.public(&name)) {
                                self.files[me].objects.insert(name, o);
                            }
                        }
                    }
                }
            }
            "param" => {
                let usage = "usage: (param name default-value [(tol lo hi)]). tol: deviations from the value, lo <= hi. (param d 8 (tol 0 0.022)) is 8.000..8.022";
                let (name, default, tol) = match args {
                    [name, default] => (name, default, None),
                    [name, default, tol] => (name, default, Some(tol)),
                    _ => return Err(self.err_hint(form, "bad param", Some(usage.into()))),
                };
                let n = self.name_of(name)?;
                self.check_param_name(name, &n, "the param")?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                let (mut v, mut said) = self.constant(default, "param default must be a constant number")?;
                let mut tol = tol.map(|t| self.tol(t, said.as_deref(), usage).map(|(b, _)| b)).transpose()?;
                // Its import's value (ImportValue), in the quantity written
                // here, else in its own.
                if let Some(g) = self.files[self.scope.file as usize].given.iter().find(|g| g.name == w).cloned() {
                    if let Some(s) = said.as_ref().or(g.said.as_ref()) {
                        g.fits(s.unit.quantity, &format!("{}'s {w}", self.file))?;
                    }
                    self.set_at.insert(self.params.len(), (v, format!("{}:{}", g.file, g.line)));
                    (v, said, tol) = (g.value, g.said.or(said), g.tol.or(tol));
                }
                let key = self.defined(&w, name.line);
                if let Some(s) = self.param_given(&n, said, form)? {
                    self.param_said.insert(self.params.len(), s);
                }
                self.param_keys.insert(key, self.params.len());
                let p = Param { name: n, default: v, tol };
                if let Some(w) = p.outside_band() {
                    let w = format!("{}:{}: {w}", self.file, form.line);
                    if !self.load_warnings.contains(&w) {
                        self.load_warnings.push(w);
                    }
                }
                self.params.push(p);
            }
            "def" => {
                if let Some(written) = fn_written_as_def(args) {
                    return Err(self.err_hint(form, "a function is defined by fn, not def", Some(format!("write {written}"))));
                }
                let [name, e] = args else { return Err(self.err(form, "usage: (def name expr)")) };
                let n = self.name_of(name)?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                let made = self.ranged.len();
                let mut v = self.eval(e, &mut Vec::new())?;
                // A (ranged lo hi) the def makes and is takes its name, the
                // one --param reaches it by (Doc::ranged_named).
                if let Some((k, _)) = self.ranged_of(&v).filter(|r| self.ranged[made..].contains(r)) {
                    self.params[k].name = n.clone();
                    if let Value::Field(f, said) = v {
                        let said = self.param_given(&n, said, e)?;
                        self.param_said.remove(&k);
                        self.param_said.extend(said.clone().map(|s| (k, s)));
                        v = Value::Field(f, said);
                    }
                }
                // Every def of a spatial field names a construction: its surface
                // points report this name (see eval::provenance). Constants stay
                // untagged so they keep folding.
                if let Value::Field(f, said) = &v
                    && self.g.depends_on(*f, |x| matches!(x, Var::X | Var::Y | Var::Z))
                {
                    v = Value::Field(self.g.tag(&n, *f), said.clone());
                }
                let key = self.defined(&w, name.line);
                self.def_index.insert(key, self.defs.len());
                if !self.origin.contains_key(&n) || self.scope.file != 0 {
                    self.def_index.insert(n.clone(), self.defs.len());
                }
                self.defs.push((n, v));
            }
            "fn" => {
                let [name, ps, body] = args else { return Err(self.err(form, "usage: (fn name (arg ...) body)")) };
                let n = self.name_of(name)?;
                self.check_new_value_name(name, &n)?;
                let (w, n) = (n.clone(), self.public(&n));
                let (params, rest) = self.arg_names(ps, "fn")?;
                let key = self.defined(&w, name.line);
                let def = Rc::new(FnDef { params, rest, body: body.clone(), file: self.file.clone(), scope: self.scope });
                self.fns.insert(key, def.clone());
                if !self.origin.contains_key(&n) || self.scope.file != 0 {
                    self.fns.insert(n, def);
                }
            }
            "body" => self.top_def(form, args, true)?,
            "assembly" => self.top_def(form, args, false)?,
            "assert" => {
                let [name, claim, options @ ..] = args else { return Err(self.err(form, "usage: (assert name claim [(res r)])")) };
                let n = self.public(&self.name_of(name)?);
                if self.asserts.iter().any(|a| a.name == n) {
                    return Err(self.err(name, format!("assert '{n}' already exists")));
                }
                let (c, bound) = self.claim(claim)?;
                let (res, tol) = self.assert_options(&c, options)?;
                let quantity = self.accuracy_quantity(&c);
                self.asserts.push(AssertDef {
                    name: n,
                    claim: c,
                    text: claim.to_string(),
                    bound,
                    res,
                    accuracy: tol,
                    quantity,
                    file: self.file.clone(),
                    line: form.line,
                    terms: Vec::new(),
                    fixed: Vec::new(),
                });
            }
            "t-span" => {
                let usage = "usage: (t-span A B): claims that read t are checked at every moment from A to B";
                let [a, b] = args else { return Err(self.err_hint(form, "wrong arguments", Some(usage.into()))) };
                let end = |d: &mut Doc, e: &Sexp| {
                    let f = d.field_of(e, Quantity::Time, "t-span", &mut Vec::new())?;
                    d.const_of(e, f)
                };
                let (a, b) = (end(self, a)?, end(self, b)?);
                if !(a < b && a.is_finite() && b.is_finite()) {
                    return Err(self.err_hint(form, "a time span needs finite A < B", Some(usage.into())));
                }
                // The file's motion, which its own claims and assemblies are
                // checked over, in an importer too (span_of); the run file's
                // is the default of what it states.
                if self.import_depth == 0 {
                    self.t_span = Some((a, b));
                }
                self.file_spans.insert(self.file.clone(), (a, b));
            }
            other => {
                if let Some((package, head)) = self.registry.package(other) {
                    return self.top_ext(package, head, form, args);
                }
                let all = self.registry.top_forms();
                return Err(self.err_hint(&items[0], format!("unknown top-level form '{other}'"), closest(other, all.split(", ")).or(Some(format!("top-level forms: {all}")))));
            }
        }
        Ok(())
    }

    /// `(body HEAD clauses...)` or `(assembly HEAD clauses...)`: HEAD is
    /// NAME, or `(NAME (P DEFAULT) ...)` for one taking arguments. Built
    /// here with its defaults; an instance passing other values builds it
    /// again (`instantiate`).
    fn top_def(&mut self, form: &Sexp, args: &[Sexp], body: bool) -> R<()> {
        let usage = if body {
            "usage: (body name (shape expr) (attr key value) (frame name F) ...), or with arguments (body (name (p default) ...) ...)"
        } else {
            "usage: (assembly name (instance name child [WHERE]) (frame name F) ...), or with arguments (assembly (name (p default) ...) ...)"
        };
        let Some((head, clauses)) = args.split_first() else { return Err(self.err(form, usage)) };
        let (name, ps) = match head.list() {
            Some([name, ps @ ..]) => (name, ps),
            Some([]) => return Err(self.err(head, usage)),
            None => (head, &[][..]),
        };
        let w = self.name_of(name)?;
        let n = self.public(&w);
        self.check_new_object_name(name, &n)?;
        self.files[self.scope.file as usize].objects.insert(w, n.clone());
        let mut params: Vec<(String, Sexp)> = Vec::new();
        for p in ps {
            let Some([pn, d]) = p.list() else { return Err(self.err_hint(p, "an argument is (name default)", Some(usage.into()))) };
            let pn = self.name_of(pn)?;
            if params.iter().any(|(q, _)| *q == pn) {
                return Err(self.err(p, format!("duplicate argument '{pn}'")));
            }
            params.push((pn, d.clone()));
        }
        let def = Rc::new(ArgDef { params, clauses: clauses.to_vec(), file: self.file.clone(), scope: self.scope });
        let (env, texts) = self.def_args(&def, &n, Vec::new(), form)?;
        let key = format!("{n}{env:?}");
        self.build(n.clone(), None, &def, env, texts, None, body, form)?;
        if !def.params.is_empty() {
            self.built.insert(key, n.clone());
            self.arg_defs.insert(n, def);
        }
        Ok(())
    }

    /// The values a definition's arguments take: those passed (as values,
    /// and as written), then the defaults, each read in the definition's
    /// file with the arguments before it.
    fn def_args(&mut self, def: &ArgDef, name: &str, args: Vec<(Value, String)>, at: &Sexp) -> R<(Env, Vec<(String, String)>)> {
        if args.len() > def.params.len() {
            let ps: Vec<&str> = def.params.iter().map(|p| p.0.as_str()).collect();
            return Err(self.err(at, format!("'{name}' takes {} argument(s) ({}), got {}", ps.len(), ps.join(" "), args.len())));
        }
        let mut env = Env::new();
        let mut texts = Vec::new();
        let mut args = args.into_iter();
        let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
        let r = (|| {
            for (p, d) in &def.params {
                let (v, text) = match args.next() {
                    Some(a) => a,
                    None => (self.eval(d, &mut env)?, d.to_string()),
                };
                env.push((p.clone(), v));
                texts.push((p.clone(), text));
            }
            Ok(())
        })();
        (self.file, self.scope) = saved;
        r.map(|_| (env, texts))
    }

    /// A body or an assembly named `n` from its clauses, its arguments
    /// taking `env`, in the names of the file it is written in; added to
    /// the document. `passed`: where the arguments were passed (the claims
    /// an assembly's joints add are there). A body's clauses read only its
    /// arguments: not the instances of an assembly it is built for.
    fn build(&mut self, n: String, of: Option<String>, def: &ArgDef, env: Env, params: Vec<(String, String)>, passed: Option<(Rc<str>, u32)>, body: bool, at: &Sexp) -> R<()> {
        // The file that defines it: where it is written, or where the
        // arguments it is built with were passed.
        let file = passed.as_ref().map_or(def.file.clone(), |p| p.0.clone());
        let (name, defined) = (n.clone(), of.is_none());
        let building = if body {
            None
        } else {
            let asm = Assembly {
                name: n.clone(),
                params: params.clone(),
                of: of.clone(),
                instances: Vec::new(),
                frames: Vec::new(),
                frame_on: Vec::new(),
                relations: Vec::new(),
                clauses: Vec::new(),
                file: file.clone(),
                mass: None,
            };
            Some(Building { asm, later: Vec::new(), passed, mass: Vec::new(), read: Vec::new() })
        };
        let parent = std::mem::replace(&mut self.building, building);
        let args = std::mem::replace(&mut self.args, env);
        let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
        let r = (|| {
            let clauses = self.expand_for(&def.clauses)?;
            if body {
                return self.body_clauses(n, of, params, &clauses, at).map(Some);
            }
            let later = clauses.iter().filter_map(|c| match c.list() {
                Some([h, i, ..]) if h.sym() == Some("instance") => i.sym().map(str::to_string),
                _ => None,
            });
            self.building.as_mut().unwrap().later = later.collect();
            self.assembly_clauses(&clauses).map(|_| None)
        })();
        (self.file, self.scope) = saved;
        self.args = args;
        let b = std::mem::replace(&mut self.building, parent);
        match r? {
            Some(body) => {
                let body = self.weighed_as_defined(body);
                self.bodies.push(Body { file, ..body })
            }
            None => {
                let b = b.unwrap();
                self.assemblies.push(b.asm);
                // Built again with other arguments, it places the bodies its
                // definition weighed.
                let get = |k: &str| b.mass.iter().find(|m| m.0 == k);
                if let (Some((_, mass, at)), true) = (get("mass"), defined) {
                    let saved = (std::mem::replace(&mut self.file, def.file.clone()), std::mem::replace(&mut self.scope, def.scope));
                    let r = self.assembly_mass(at, &name, mass, get("mass-params").map(|m| &m.1));
                    (self.file, self.scope) = saved;
                    r?;
                }
            }
        }
        Ok(())
    }

    /// A body from its clauses (expanded): one shape, attributes, frames.
    fn body_clauses(&mut self, n: String, of: Option<String>, params: Vec<(String, String)>, clauses: &[Sexp], at: &Sexp) -> R<Body> {
        let what = "body clauses are (shape expr), (attr key value), (frame name F) or (mark name (vec g r) ...)";
        let mut shape = None;
        let mut attrs: Vec<(String, Value)> = Vec::new();
        let mut frames = Vec::new();
        let mut marks = Vec::new();
        for c in clauses {
            let parts = c.list().ok_or_else(|| self.err(c, what))?;
            match parts.first().and_then(|h| h.sym()) {
                Some("shape") => {
                    let [_, e] = parts else { return Err(self.err(c, "usage: (shape expr)")) };
                    if shape.is_some() {
                        return Err(self.err(c, "body has more than one (shape ...)"));
                    }
                    shape = Some(self.field(e, &mut self.arg_env())?);
                }
                Some("attr") => {
                    let [_, k, v] = parts else { return Err(self.err(c, "usage: (attr key value)")) };
                    let k = self.name_of(k)?;
                    if attrs.iter().any(|(a, _)| *a == k) {
                        return Err(self.err(c, format!("duplicate attr '{k}'")));
                    }
                    let e = v;
                    let v = self.eval(e, &mut self.arg_env())?;
                    let q = match k.as_str() {
                        "density" => Some(Quantity::Density),
                        "mass" => Some(Quantity::Mass),
                        "spin-inertia" => Some(Quantity::Inertia),
                        _ => None,
                    };
                    if let Some(q) = q {
                        self.fits(&v, e, q, &format!("attr {k}"))?;
                    }
                    // A look is read once, here: clients take look::Look,
                    // never the layout.
                    let v = look::read(&mut self.g, &k, v).map_err(|e| self.err(c, format!("attr {k}: {e}")))?;
                    // A density is g/mm^3: a number in another unit is
                    // refused here, not read a million times off.
                    let density = match (k.as_str(), &v) {
                        ("density", Value::Field(f, _)) => materials::density_field_error(&self.g, *f),
                        ("material", Value::Data(d)) => materials::material_density_error(d.value()),
                        _ => None,
                    };
                    if let Some(e) = density {
                        return Err(self.err(c, format!("attr {k}: {e}")));
                    }
                    // A material is something physics can look up: a name
                    // in std/materials.ink (or data). Free text would make
                    // every study reading it fail far from here. Another
                    // spelling of a material is refused, naming it.
                    if let (true, Value::Str(m)) = (k == "material", &v)
                        && materials::lookup(m).is_none()
                    {
                        if materials::instead(m).is_some() {
                            let (what, hint) = materials::refused(m);
                            return Err(self.err_hint(c, format!("body '{}': {what}", of.as_deref().unwrap_or(&n)), Some(hint)));
                        }
                        let w = format!(
                            "{}:{}: body '{}': material '{m}' is not in std/materials.ink ({}). Name one of those, or give (data ...) with its keys. Free text goes in (attr material-note \"...\")",
                            self.file,
                            c.line,
                            of.as_deref().unwrap_or(&n),
                            materials::names().join(" ")
                        );
                        if !self.load_warnings.contains(&w) {
                            self.load_warnings.push(w);
                        }
                    }
                    attrs.push((k, v));
                }
                Some("frame") => frames.push(self.frame_clause(c, &frames)?),
                Some("mark") => marks.extend(self.mark_clause(c, parts)?),
                // A function's call is a shape written bare.
                head => {
                    let shape = head.filter(|h| self.fns.contains_key(&*self.scoped(h))).map(|h| format!("write (shape ({h} ...))"));
                    return Err(self.err_hint(c, what, self.top_form_hint(c, "a body").or(shape)));
                }
            }
        }
        let shape = shape.ok_or_else(|| self.err(at, format!("body '{n}' needs a (shape expr) clause")))?;
        // The body's name is a construction name too (provenance: faces of this body).
        let shape = self.g.tag(&n, shape);
        self.published_mass(at, &n, of.as_deref(), shape, &mut attrs)?;
        Ok(Body { name: n, params, of, shape, attrs: Attrs::new(attrs), frames, marks, file: self.file.clone() })
    }

    /// `(mark NAME CURVE ...)`: curves drawn on the part, each CURVE `(vec
    /// g r)`, where field g is 0 inside field r < 0, or a vec of them.
    fn mark_clause(&mut self, c: &Sexp, parts: &[Sexp]) -> R<Vec<(String, NodeId, NodeId)>> {
        let usage = "usage: (mark NAME CURVE ...): curves drawn on the part (a sheet's bend lines: (mark bend (sheet-bend-lines S))), \
                     each CURVE (vec g r), where field g is 0 inside field r < 0, or a vec of them (none: (vec)); `ink slice --out` writes them on layer NAME";
        let [_, n, curves @ ..] = parts else { return Err(self.err_hint(c, "bad mark", Some(usage.into()))) };
        let name = self.name_of(n)?;
        let mut out = Vec::new();
        fn add(v: &Value, out: &mut Vec<(NodeId, NodeId)>) -> bool {
            let Value::Vec(l) = v else { return false };
            match l.as_slice() {
                [Value::Field(g, _), Value::Field(r, _)] => {
                    out.push((*g, *r));
                    true
                }
                l => l.iter().all(|x| add(x, out)),
            }
        }
        for e in curves {
            let v = self.eval(e, &mut self.arg_env())?;
            let mut got = Vec::new();
            if !add(&v, &mut got) {
                return Err(self.err_hint(e, format!("a mark's curve is (vec g r), two fields, or a vec of them, not a {}", v.kind()), Some(usage.into())));
            }
            out.extend(got.into_iter().map(|(g, r)| (name.clone(), g, r)));
        }
        Ok(out)
    }

    /// What an instance named `iname` places for `(CHILD ARG ...)`: the
    /// body or assembly CHILD built with the arguments' values (once per
    /// distinct values), named after this instance's path.
    fn instantiate(&mut self, c: &Sexp, iname: &str, child: &str, args: &[Sexp]) -> R<String> {
        let body = self.body(child).is_some();
        let Some(def) = self.arg_defs.get(child).cloned() else {
            return Err(self.err_hint(
                c,
                format!("'{child}' is {} without arguments: it takes none", if body { "a body" } else { "an assembly" }),
                Some("write the child as a name, or give it arguments: (body (name (p default) ...) ...), (assembly (name (p default) ...) ...)".into()),
            ));
        };
        let mut vals = Vec::new();
        for a in args {
            let mut env = self.arg_env();
            let v = self.eval(a, &mut env)?;
            if let Value::Field(f, _) = v
                && self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z))
            {
                return Err(self.err(a, "a body's or an assembly's arguments are values of params and t, not fields of x y z"));
            }
            vals.push((v, a.to_string()));
        }
        // A load with coordinates: each number the target's own instances
        // pass is its value plus a param at 0 (Coord).
        let b = self.building.as_ref().unwrap();
        if self.coordinates_of.as_deref() == Some(b.asm.name.as_str()) && b.asm.of.is_none() && self.import_depth == 0 {
            for ((v, _), (p, _)) in vals.iter_mut().zip(&def.params) {
                let Value::Field(f, ref said) = *v else { continue };
                let said = said.clone();
                let name = format!("{iname}.{p}");
                let param = |d: &mut Doc, what: &str| {
                    d.params.push(Param { name: format!("{what}{name}"), default: 0.0, tol: None });
                    (d.params.len() - 1, d.g.var(Var::Param(d.params.len() as u32 - 1)))
                };
                // A value that is its coordinate's own (a simulated motion:
                // Reader::as_coordinate), passed itself, is that coordinate
                // (once): what else reads it moves along it too. An argument
                // that otherwise reads such a value (twice it, its negative,
                // a sum), or passes it again, is no coordinate of its own:
                // passed as it is, it moves along that one.
                let own = self.coordinate_params.iter().find(|(x, _)| *x == f).map(|&(_, k)| k);
                let reads_own = self.coordinate_params.iter().any(|&(_, k)| self.g.depends_on(f, |x| x == Var::Param(k as u32)));
                if reads_own && (own.is_none() || self.coords.iter().any(|c| Some(c.param) == own)) {
                    continue;
                }
                let (k, e) = match own {
                    Some(k) => (k, self.g.var(Var::Param(k as u32))),
                    None => param(self, "coordinate "),
                };
                if self.states && self.g.depends_on(f, |x| x == Var::T) {
                    // At t0: e, at speed rate and acceleration acc.
                    let ((kr, rate), (ka, acc)) = (param(self, "rate of coordinate "), param(self, "acceleration of coordinate "));
                    let t0 = match self.state_t {
                        Some(k) => self.g.var(Var::Param(k as u32)),
                        None => {
                            self.params.push(Param { name: "moment of the state".into(), default: 0.0, tol: None });
                            self.state_t = Some(self.params.len() - 1);
                            self.g.var(Var::Param(self.params.len() as u32 - 1))
                        }
                    };
                    let t = self.g.var(Var::T);
                    let dt = self.g.sub(t, t0);
                    let (half, dt2) = (self.g.constant(0.5), self.g.mul(dt, dt));
                    let (r, a) = (self.g.mul(rate, dt), self.g.mul(acc, dt2));
                    let a = self.g.mul(half, a);
                    let ra = self.g.add(r, a);
                    let passed = self.g.add(e, ra);
                    *v = Value::Field(passed, said);
                    self.coords.push(Coord { name, param: k, value: f, jet: Some((kr, ka)), passed });
                } else {
                    let passed = if own.is_some() { f } else { self.g.add(f, e) };
                    *v = Value::Field(passed, said);
                    self.coords.push(Coord { name, param: k, value: f, jet: None, passed });
                }
            }
        }
        let (env, texts) = self.def_args(&def, child, vals, c)?;
        // Equal values build equal bodies and assemblies (fields are
        // interned: equal expressions are one node).
        let key = format!("{child}{env:?}");
        if let Some(n) = self.built.get(&key) {
            return Ok(n.clone());
        }
        let b = self.building.as_ref().unwrap();
        let name = format!("{}/{iname}", b.asm.name);
        let passed = b.passed.clone().unwrap_or((self.file.clone(), c.line));
        self.calls.push((child.to_string(), self.file.clone(), c.line));
        let r = self.build(name.clone(), Some(child.to_string()), &def, env, texts, Some(passed), body, c);
        self.calls.pop();
        r.map_err(|mut e| {
            e.called_from.push((child.to_string(), self.file.to_string(), c.line, c.col));
            e
        })?;
        self.built.insert(key, name.clone());
        Ok(name)
    }

    /// The values of the arguments of the body or assembly being read.
    fn arg_env(&self) -> Env {
        self.args.clone()
    }

    fn assembly_clauses(&mut self, clauses: &[Sexp]) -> R<()> {
        let owned: Vec<String> = self.registry.packages.iter().flat_map(|p| p.clauses().iter().map(|h| format!(", ({h} ...)"))).collect();
        let what = format!("assembly clauses are (instance ...), (frame name F), (contact A B ZONE), (relate A B KIND){} or (attr mass M)", owned.concat());
        let what = what.as_str();
        for c in clauses {
            let parts = c.list().ok_or_else(|| self.err(c, what))?;
            match parts.first().and_then(|h| h.sym()) {
                Some("instance") => {
                    if parts.len() < 3 {
                        return Err(self.err_hint(c, "usage: (instance name child [WHERE [CFRAME] [MOTION] [(limits LO HI)] [(speed-limits LO HI)]])", Some(INSTANCE.into())));
                    }
                    let iname = self.name_of(&parts[1])?;
                    let asm = &self.building.as_ref().unwrap().asm;
                    if iname.contains('/') || asm.instances.iter().any(|i| i.name == iname) {
                        return Err(self.err(&parts[1], format!("instance name '{iname}' is invalid or duplicated")));
                    }
                    // (CHILD ARG ...): an assembly given arguments.
                    let (head, args) = match parts[2].list() {
                        Some([h, args @ ..]) => (h, args),
                        _ => (&parts[2], &[][..]),
                    };
                    let mut child = self.path_of(head)?;
                    if self.body(&child).is_none() && self.assembly(&child).is_none() {
                        let names = self
                            .bodies
                            .iter()
                            .filter(|b| b.of.is_none())
                            .map(|b| b.name.as_str())
                            .chain(self.assemblies.iter().filter(|a| a.of.is_none()).map(|a| a.name.as_str()));
                        return Err(self.err_hint(head, format!("no body or assembly named '{child}' (define it before the assembly)"), closest(&child, names)));
                    }
                    if parts[2].list().is_some() {
                        child = self.instantiate(&parts[2], &iname, &child, args)?;
                    }
                    let (at, joint, on, met) = self.instance_placement(c, &iname, &child, &parts[3..])?;
                    let b = self.building.as_mut().unwrap();
                    b.later.retain(|l| *l != iname);
                    b.asm.instances.push(Instance { name: iname, child, at, joint, on, met });
                }
                Some("frame") => {
                    let frames = self.building.as_ref().unwrap().asm.frames.clone();
                    self.building.as_mut().unwrap().read.clear();
                    let f = self.frame_clause(c, &frames)?;
                    let b = self.building.as_mut().unwrap();
                    b.asm.frame_on.push((f.0.clone(), std::mem::take(&mut b.read)));
                    b.asm.frames.push(f);
                }
                Some("contact") => {
                    let usage = "(contact A B ZONE): instances A and B may overlap inside ZONE, a field in A's own coordinates (a thread in a nut: (shape-of screw-body thread))";
                    let [_, a, b, zone] = parts else { return Err(self.err_hint(c, "bad contact", Some(usage.into()))) };
                    let (pa, pb) = self.relation_ends(a, b, usage)?;
                    let zone = Zone { a: self.field(zone, &mut self.arg_env())?, b: None };
                    self.building.as_mut().unwrap().asm.relations.push(Relation { a: pa, b: pb, zone: Some(zone), kind: None });
                }
                Some("relate") => self.relate(c, parts)?,
                Some(h) if let Some((package, head)) = self.registry.clause(h) => {
                    let mut r = Reader::in_env(self, self.arg_env());
                    let kept = package.clause(&mut r, c, head, &parts[1..])?;
                    let reads = r.into_reads();
                    let (file, line) = (self.file.clone(), c.line);
                    self.building.as_mut().unwrap().asm.clauses.push(OwnedClause { package: package.name(), head, kept, reads, file, line });
                }
                Some("attr") => {
                    let usage = "an assembly's attrs are its published (attr mass M) and (attr mass-params (read-json FILE)), FILE a JSON object of param values";
                    let [_, k, v] = parts else { return Err(self.err_hint(c, "usage: (attr key value)", Some(usage.into()))) };
                    let k = self.name_of(k)?;
                    if !matches!(k.as_str(), "mass" | "mass-params") {
                        return Err(self.err_hint(c, format!("an assembly has no attr '{k}'"), Some(usage.into())));
                    }
                    if self.building.as_ref().unwrap().mass.iter().any(|m| m.0 == k) {
                        return Err(self.err(c, format!("duplicate attr '{k}'")));
                    }
                    let e = v;
                    let v = self.eval(e, &mut self.arg_env())?;
                    if k == "mass" {
                        self.fits(&v, e, Quantity::Mass, "attr mass")?;
                    }
                    self.building.as_mut().unwrap().mass.push((k, v, c.clone()));
                }
                _ => return Err(self.err_hint(c, what, self.top_form_hint(c, "an assembly"))),
            }
        }
        Ok(())
    }

    /// The body or assembly a relation's end (an instance name, or a path
    /// into a sub-assembly, read so far) places.
    fn end_child(&self, path: &str) -> Result<String, String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or_default();
        let b = self.building.as_ref().unwrap();
        let mut child = b.asm.instances.iter().find(|i| i.name == first).map(|i| i.child.clone()).ok_or_else(|| format!("no instance '{first}'"))?;
        for seg in segs {
            child = self
                .assembly(&child)
                .and_then(|a| a.instances.iter().find(|i| i.name == seg))
                .map(|i| i.child.clone())
                .ok_or_else(|| format!("no instance '{seg}' in '{child}'"))?;
        }
        Ok(child)
    }

    /// The shape of a relation's end (an instance name, or a path into a
    /// sub-assembly, read so far) in its own coordinates: its body's, or
    /// the union of its assembly's bodies.
    fn end_shape(&mut self, path: &str) -> Result<NodeId, String> {
        let child = self.end_child(path)?;
        // The child by its name, which is a path where it was built with
        // arguments (`pair/p`), not an instance path to resolve.
        if let Some(b) = self.body(&child) {
            return Ok(b.shape);
        }
        let id = self.identity();
        let items = self.flatten(&child, &id, &child, &[])?;
        let shapes: Vec<NodeId> = items.iter().map(|i| self.placed_field(i)).collect();
        Ok(self.reduce(BinOp::Min, &shapes))
    }

    /// The ends of a contact or relation: instance names, or paths into
    /// sub-assemblies, read so far.
    fn relation_ends(&self, a: &Sexp, b: &Sexp, usage: &str) -> R<(String, String)> {
        let (pa, pb) = (self.name_of(a)?, self.name_of(b)?);
        let asm = &self.building.as_ref().unwrap().asm;
        for (who, s) in [(&pa, a), (&pb, b)] {
            self.check_contact_path(asm, who).map_err(|m| self.err_hint(s, m, Some(usage.into())))?;
        }
        Ok((pa, pb))
    }

    /// `(relate A B KIND)`: KIND, a function, called with A's and B's
    /// placements in this assembly's coordinates, gives a `(relation ...)`:
    /// its zone is intended overlap, as a contact's; each of its claims is
    /// the value claim `relate:ASSEMBLY/A~B/NAME`, stated where the
    /// assembly's arguments were passed, else here (as limits are). An end
    /// written `(PATH FRAME)` is read at FRAME, a frame in its own
    /// coordinates: the kind gets its placement composed with FRAME and its
    /// shape in FRAME's coordinates, and its zone there is put back in the
    /// end's own.
    fn relate(&mut self, c: &Sexp, parts: &[Sexp]) -> R<()> {
        let usage = "(relate A B KIND): KIND a function of A's and B's frames (in this assembly's coordinates) and their shapes (each in \
                     its own coordinates) giving a (relation ...), as std's kinds do, e.g. (relate pinion wheel (gear-mesh 15 45 1)); an \
                     end (PATH FRAME) is read at FRAME, a frame in its own coordinates (a gear whose axis is not its z)";
        let [_, a, b, kind] = parts else { return Err(self.err_hint(c, "bad relate", Some(usage.into()))) };
        // An end's FRAME is in its part's own coordinates: a frame the part
        // declares, by name, or one written there.
        let end = |me: &mut Self, s: &Sexp| -> R<(Sexp, Option<Xf>)> {
            let Some(l) = s.list() else { return Ok((s.clone(), None)) };
            let [p, f] = l else { return Err(me.err_hint(s, "a relation's end is a path or (PATH FRAME)", Some(usage.into()))) };
            if let Some(name) = f.sym() {
                let child = me.end_child(&me.name_of(p)?).map_err(|m| me.err(p, m))?;
                let frames = me.frames_of(&child);
                let Some((_, xf)) = frames.iter().find(|(n, _)| n == name) else { return Err(me.err(f, Self::no_frame(&child, name, &frames))) };
                return Ok((p.clone(), Some(*xf)));
            }
            me.own_frame = true;
            let xf = me.frame_value(f, &mut me.arg_env(), "a relation end's frame");
            me.own_frame = false;
            Ok((p.clone(), Some(xf?)))
        };
        let ((a, la), (b, lb)) = (end(self, a)?, end(self, b)?);
        let (pa, pb) = self.relation_ends(&a, &b, usage)?;
        let f = match self.eval(kind, &mut self.arg_env())? {
            Value::Fn(f) => f,
            other => return Err(self.err_hint(kind, format!("KIND must be a function of two frames and two shapes, not a {}", other.kind()), Some(usage.into()))),
        };
        let read = |me: &mut Self, p: &str, l: &Option<Xf>| -> R<(Value, NodeId)> {
            let (f, s) = (me.frame_of(p, None).map_err(|m| me.err(c, m))?, me.end_shape(p).map_err(|m| me.err(c, m))?);
            let Some(l) = l else { return Ok((Value::Frame(f), s)) };
            let (inv, at) = (me.inverse(l), me.compose(l, &f));
            let Value::Field(s, _) = me.to_world(&Value::Field(s, None), &inv) else { unreachable!("a placed field is a field") };
            Ok((Value::Frame(at), s))
        };
        let ((fa, sa), (fb, sb)) = (read(self, &pa, &la)?, read(self, &pb, &lb)?);
        let before: HashSet<usize> = self.materialized.borrow().keys().copied().collect();
        let r = match self.call(&f, vec![fa, fb, Value::Field(sa, None), Value::Field(sb, None)], kind)? {
            Value::Relation(r) => r,
            other => return Err(self.err_hint(kind, format!("KIND gave a {}, not a (relation ...)", other.kind()), Some(usage.into()))),
        };
        // A ranged param the kind fixed at load (by bounds of a part's shape)
        // is this relation's claims' to answer for, not the file's: a run
        // ranging it refuses them alone (AssertDef::fixed).
        let fixed: Vec<usize> = {
            let mut m = self.materialized.borrow_mut();
            let mut new: Vec<usize> = m.keys().copied().filter(|k| !before.contains(k) && self.has_range(*k)).collect();
            new.sort_unstable();
            new.iter().for_each(|k| _ = m.remove(k));
            new
        };
        let asm = self.building.as_ref().unwrap().asm.name.clone();
        let (file, line) = self.building.as_ref().unwrap().passed.clone().unwrap_or((self.file.clone(), c.line));
        let mut claims = Vec::new();
        for (n, v, lo, hi, tol) in &r.claims {
            let name = format!("relate:{asm}/{pa}~{pb}/{n}");
            if self.asserts.iter().any(|x| x.name == name) {
                return Err(self.err(c, format!("'{pa}' and '{pb}' are related twice with a claim '{n}'")));
            }
            let (l, h) = (self.g.sub(*v, *lo), self.g.sub(*hi, *v));
            let claim = Claim::Value(self.g.min(l, h), Cmp::Ge, 0.0);
            let text = format!("{} {pa}~{pb}: lo <= {n} <= hi, as value = min({n} - lo, hi - {n}) >= 0", r.kind);
            let terms = vec![(n.clone(), *v), ("lo".into(), *lo), ("hi".into(), *hi)];
            self.asserts.push(AssertDef {
                name: name.clone(),
                claim,
                text,
                bound: None,
                res: None,
                accuracy: *tol,
                quantity: None,
                file: file.clone(),
                line,
                terms,
                fixed: fixed.clone(),
            });
            claims.push(name);
        }
        let site = format!("{}:{}", self.file, c.line);
        let kind = Kind { name: r.kind.clone(), skipped: r.zone.as_ref().map(|z| z.0.clone()), claims, data: r.data.clone(), clauses: r.clauses.clone(), site };
        let back = |me: &mut Self, z: NodeId, l: &Option<Xf>| {
            let Some(l) = l else { return z };
            let Value::Field(z, _) = me.to_world(&Value::Field(z, None), l) else { unreachable!("a placed field is a field") };
            z
        };
        let zone = r.zone.as_ref().map(|z| z.1).map(|z| Zone { a: back(self, z.a, &la), b: z.b.map(|zb| back(self, zb, &lb)) });
        self.building.as_mut().unwrap().asm.relations.push(Relation { a: pa, b: pb, zone, kind: Some(Rc::new(kind)) });
        Ok(())
    }

    /// CLAUSES with each `(optional CLAUSE)` replaced by what CLAUSE expands
    /// to (a `for`'s clauses) where a package this ink loads owns it, and
    /// left out where none does (a mesh's coupling, a stepper's stiffness:
    /// the kernel stands without them). A kernel clause (KERNEL_CLAUSES, of
    /// either kind, `for`'s included) in it is an error: it is never optional. One no loaded
    /// package owns warns where some loaded package reads clauses of this
    /// kind (HEADS, a package's list of them), since there it can only be a
    /// misspelling or a package this ink lacks, and dropping it would drop
    /// that package's check; where none does (the kernel's ink), nothing
    /// could read it.
    fn optional_clauses(&mut self, clauses: &[Sexp], env: &mut Env, heads: fn(&dyn ext::Package) -> &'static [&'static str], usage: &str) -> R<Vec<Sexp>> {
        let mut out = Vec::with_capacity(clauses.len());
        for c in clauses {
            let inner = match c.list().unwrap_or(&[]) {
                [h, inner] if h.sym() == Some("optional") => inner,
                _ => {
                    out.push(c.clone());
                    continue;
                }
            };
            for x in self.expand_for_in(std::slice::from_ref(inner), env)? {
                let head = x.list().and_then(|l| l.first()).and_then(|h| h.sym());
                match head {
                    Some(h) if self.registry.packages.iter().any(|p| heads(*p).contains(&h)) => out.push(x.clone()),
                    Some(h) if !KERNEL_CLAUSES.contains(&h) && h != "optional" => {
                        let loaded: Vec<&str> = self.registry.packages.iter().filter(|p| !heads(**p).is_empty()).map(|p| p.name()).collect();
                        if !loaded.is_empty() {
                            let w = format!(
                                "{}:{}: (optional ({h} ...)) is left out: no package this ink loads ({}) reads '{h}'. Check its spelling",
                                self.file,
                                x.line,
                                loaded.join(", ")
                            );
                            if !self.load_warnings.contains(&w) {
                                self.load_warnings.push(w);
                            }
                        }
                    }
                    _ => return Err(self.err_hint(&x, "(optional CLAUSE) is for a package's clause, not the kernel's", Some(usage.into()))),
                }
            }
        }
        Ok(out)
    }

    /// `(relation NAME [(zone WHAT ZA [ZB])] (claim NAME VALUE LO HI) ...
    /// [(data KEY VALUE ...)] [CLAUSE ...])`: what a kind gives `relate`.
    fn eval_relation(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let usage = "(relation NAME [(zone WHAT ZA [ZB])] (claim NAME VALUE LO HI [(accuracy x)]) ... \
                     [(data KEY VALUE ...)] [CLAUSE ...] [(optional CLAUSE)]), the value of a kind for (relate A B KIND): overlap of A and B \
                     inside ZA (in A's coordinates) and ZB (in B's) is intended, WHAT saying what that leaves unchecked; each claim LO <= \
                     VALUE <= HI is a value claim (numbers, params, t, frames), allowing VALUE x past them with (accuracy x); data is shown with the \
                     pair; CLAUSE a package's, refused where no package this ink loads owns it, (optional CLAUSE) left out there";
        let Some((n, clauses)) = args.split_first() else { return Err(self.err_hint(e, "bad relation", Some(usage.into()))) };
        let kind = match &n.k {
            SK::Sym(s) | SK::Str(s) => s.clone(),
            _ => return Err(self.err_hint(n, "a relation's name is a name or a string", Some(usage.into()))),
        };
        let mut r = RelationValue { kind, zone: None, claims: Vec::new(), data: Vec::new(), clauses: Vec::new() };
        let mut data = false;
        let clauses = &self.optional_clauses(clauses, env, |p| p.relation_clauses(), usage)?;
        for c in clauses {
            let parts = c.list().unwrap_or(&[]);
            let head = parts.first().and_then(|h| h.sym());
            if let Some((package, h)) = head.and_then(|h| self.registry.relation_clause(h)) {
                if r.clauses.iter().any(|x| x.head == h) {
                    return Err(self.err(c, format!("one ({h} ...) per relation")));
                }
                let mut rd = Reader::in_env(self, env.clone());
                let kept = package.clause(&mut rd, c, h, &parts[1..])?;
                let reads = rd.into_reads();
                r.clauses.push(OwnedClause { package: package.name(), head: h, kept, reads, file: self.file.clone(), line: c.line });
                continue;
            }
            match (head, &parts[1.min(parts.len())..]) {
                (Some("zone"), [what, za, zb @ ..]) if zb.len() <= 1 && r.zone.is_none() => {
                    let Value::Str(what) = self.eval(what, env)? else { return Err(self.err_hint(what, "a zone's WHAT is a string", Some(usage.into()))) };
                    let a = self.field(za, env)?;
                    let b = zb.first().map(|zb| self.field(zb, env)).transpose()?;
                    r.zone = Some((what, Zone { a, b }));
                }
                (Some("claim"), [name, v, lo, hi, tol @ ..]) if tol.len() <= 1 => {
                    let name = self.name_of(name)?;
                    if r.claims.iter().any(|x| x.0 == name) {
                        return Err(self.err(c, format!("two claims named '{name}'")));
                    }
                    let mut f = |x: &Sexp| -> R<NodeId> {
                        let f = self.field(x, env)?;
                        if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                            return Err(self.err(x, "a relation's claims are about numbers (constants, params, t, frames), not fields of x y z"));
                        }
                        Ok(f)
                    };
                    let (v, lo, hi) = (f(v)?, f(lo)?, f(hi)?);
                    let tol = tol.first().map(|t| self.accuracy_clause(t, None, env)).transpose()?;
                    r.claims.push((name, v, lo, hi, tol));
                }
                (Some("data"), kvs) if kvs.len() % 2 == 0 && !data => {
                    data = true;
                    for kv in kvs.chunks(2) {
                        let k = self.name_of(&kv[0])?;
                        let v = self.eval(&kv[1], env)?;
                        r.data.push((k, v));
                    }
                }
                (Some(h), _) if !["zone", "claim", "data"].contains(&h) => {
                    return Err(self.err_hint(
                        c,
                        format!("bad relation clause: '{h}' is not the kernel's, nor from a package this ink loads. Check its name; a package's clause needs an ink that loads that package"),
                        Some(usage.into()),
                    ));
                }
                _ => return Err(self.err_hint(c, "bad relation clause (one zone, one data)", Some(usage.into()))),
            }
        }
        Ok(Value::Relation(Rc::new(r)))
    }

    /// An instance's placement from its clauses after the child: WHERE,
    /// then optionally CFRAME (a frame the child declares, met face to face
    /// with WHERE), MOTION and limits. Steps in order: without CFRAME
    /// (MOTION, WHERE); with it (C^-1, FLIP, MOTION, WHERE), FLIP the exact
    /// half-turn about x. Limits add the claim `limits:ASSEMBLY/NAME`.
    fn instance_placement(&mut self, c: &Sexp, iname: &str, child: &str, rest: &[Sexp]) -> R<(Xf, Option<Joint>, Vec<String>, Option<String>)> {
        let Some((w, rest)) = rest.split_first() else { return Ok((self.identity(), None, Vec::new(), None)) };
        self.building.as_mut().unwrap().read.clear();
        let anchor = self.frame_value(w, &mut self.arg_env(), "WHERE")?;
        let on = std::mem::take(&mut self.building.as_mut().unwrap().read);
        // A `for` among the clauses repeats them, its list read at load time
        // with the assembly's arguments: a clause one placement has and
        // another lacks (an actuator's limits only at a supply its record
        // publishes them for).
        let rest = &self.expand_for(rest)?;
        let rest = &self.optional_clauses(rest, &mut self.arg_env(), |p| p.joint_clauses(), INSTANCE)?;
        let (mut cframe, mut motion) = (None, None);
        // (limits LO HI) and (speed-limits LO HI): bounds on the drive and
        // on its rate, each with the accuracy it states, (accuracy x).
        const BOUNDS: [&str; 2] = ["limits", "speed-limits"];
        let mut bounds: [Option<(NodeId, NodeId, String, String, Option<&Sexp>)>; 2] = [None, None];
        let mut written: [Option<[(Value, &Sexp); 2]>; 2] = [None, None];
        // Packages' joint clauses, read once the motion is known.
        let mut owned: Vec<(&'static dyn ext::Package, &'static str, &Sexp)> = Vec::new();
        for x in rest {
            let head = x.list().and_then(|l| l.first()).and_then(|h| h.sym());
            match (x.sym(), head, x.list().unwrap_or(&[])) {
                (Some(n), ..) if cframe.is_none() && motion.is_none() && bounds.iter().all(Option::is_none) && owned.is_empty() => cframe = Some((n.to_string(), x)),
                (None, Some("turn" | "slide"), [m, a, v]) => {
                    if motion.is_some() {
                        return Err(self.err(x, "one motion per instance; put fixed offsets in WHERE, e.g. (at (move 0 0 -8) (frame-of m shaft))"));
                    }
                    let axis = match a.sym() {
                        Some("x") => 0,
                        Some("y") => 1,
                        Some("z") => 2,
                        _ => return Err(self.err(a, format!("the axis of {} must be x, y or z: an axis of the frame WHERE places the instance at", m.sym().unwrap()))),
                    };
                    let drive = self.field(v, &mut self.arg_env())?;
                    if self.g.depends_on(drive, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                        return Err(self.err(v, "placements cannot depend on x, y or z"));
                    }
                    motion = Some((head == Some("turn"), axis, drive, v.to_string()));
                }
                (None, Some(k @ ("limits" | "speed-limits")), [_, lo, hi, tol @ ..]) if tol.len() <= 1 => {
                    let slot = BOUNDS.iter().position(|b| *b == k).unwrap();
                    if bounds[slot].is_some() {
                        return Err(self.err(x, format!("one ({k} LO HI) per instance")));
                    }
                    let (vl, vh) = (self.eval(lo, &mut self.arg_env())?, self.eval(hi, &mut self.arg_env())?);
                    let (l, h) = (self.as_field(vl.clone(), lo)?, self.as_field(vh.clone(), hi)?);
                    let fixed = |v: Var| matches!(v, Var::X | Var::Y | Var::Z | Var::T);
                    if self.g.depends_on(l, fixed) || self.g.depends_on(h, fixed) {
                        return Err(self.err(x, format!("{k} are numbers or params, not fields of x y z or t")));
                    }
                    if let (Some(a), Some(b)) = (self.g.as_const(l), self.g.as_const(h))
                        && a > b
                    {
                        return Err(self.err(x, format!("{k} {a} > {b}: write ({k} LO HI) with LO <= HI")));
                    }
                    bounds[slot] = Some((l, h, lo.to_string(), hi.to_string(), tol.first()));
                    written[slot] = Some([(vl, lo), (vh, hi)]);
                }
                (None, Some(h), _) if let Some((package, h)) = self.registry.joint_clause(h) => {
                    if owned.iter().any(|o| o.1 == h) {
                        return Err(self.err(x, format!("one ({h} ...) per instance")));
                    }
                    owned.push((package, h, x));
                }
                (None, Some(h), _) if !["turn", "slide", "limits", "speed-limits"].contains(&h) => {
                    let msg = format!(
                        "bad instance clause: '{h}' is not the kernel's, nor from a package this ink loads. Check its name; a package's clause needs an ink that loads that package"
                    );
                    return Err(self.err_hint(x, msg, Some(INSTANCE.into())));
                }
                _ => return Err(self.err_hint(x, "bad instance clause", Some(INSTANCE.into()))),
            }
        }
        if motion.is_none() {
            if let Some(k) = bounds.iter().position(Option::is_some) {
                return Err(self.err(c, format!("({} LO HI) bounds a motion: add (turn AXIS deg) or (slide AXIS mm)", BOUNDS[k])));
            }
            if let Some((_, h, x)) = owned.first() {
                return Err(self.err(x, format!("({h} ...) is a joint's: add (turn AXIS deg) or (slide AXIS mm)")));
            }
        }
        // A slide's limits are lengths, its speed limits speeds; a turn's
        // speed limits are turning rates (its limits are angles: degrees).
        // Their accuracies are in their units.
        let (mut tols, mut quantities) = ([None; 2], [None; 2]);
        if let Some((turn, ..)) = motion {
            for (k, b) in written.iter().enumerate() {
                let q = match (k, turn) {
                    (0, true) => None,
                    (0, false) => Some(Quantity::Length),
                    (_, true) => Some(Quantity::AngularSpeed),
                    (_, false) => Some(Quantity::Speed),
                };
                quantities[k] = q;
                for (v, e) in b.iter().flatten() {
                    if let Some(q) = q {
                        self.fits(v, e, q, BOUNDS[k])?;
                    }
                }
                if let Some(t) = bounds[k].as_ref().and_then(|b| b.4) {
                    tols[k] = Some(self.accuracy_clause(t, q, &mut self.arg_env())?);
                }
            }
        }
        let mut at = self.identity();
        if let Some((n, x)) = &cframe {
            let frames = self.frames_of(child);
            let Some((_, f)) = frames.iter().find(|(m, _)| m == n) else { return Err(self.err(x, Self::no_frame(child, n, &frames))) };
            let f = *f;
            let inv = self.inverse(&f);
            let m = self.g.constant(-1.0);
            let mut flip = self.identity();
            (flip.r[1][1], flip.r[2][2]) = (m, m);
            at = self.compose(&inv, &flip);
        }
        let below = at;
        let joint = match motion {
            Some((turn, axis, drive, text)) => {
                let z = self.g.constant(0.0);
                let m = if turn {
                    let mut u = [z; 3];
                    u[axis] = self.g.constant(1.0);
                    self.rotation(u, drive, c)?
                } else {
                    let mut m = self.identity();
                    m.t[axis] = drive;
                    m
                };
                at = self.compose(&at, &m);
                // Each bound is a value claim, LO <= v <= HI, on the drive
                // or on its rate, stated where the motion is set: where
                // the arguments were passed, else this clause.
                let asm = self.building.as_ref().unwrap().asm.name.clone();
                let (file, line) = self.building.as_ref().unwrap().passed.clone().unwrap_or((self.file.clone(), c.line));
                // The rate as the drive reads so far: taken again once the
                // file is read (rate_joints).
                let rate = self.g.derivative(drive, Var::T);
                let what = [("limits", drive, text.clone(), "drive"), ("speed-limits", rate, format!("(diff {text} t)"), "rate")];
                for ((((lo, hi, lt, ht, _), tol), quantity), (kind, v, vt, term)) in
                    bounds.iter().zip(tols).zip(quantities).zip(what).filter_map(|(((b, t), q), w)| Some((((b.as_ref()?, t), q), w)))
                {
                    if kind == "speed-limits" {
                        self.rated.push((self.asserts.len(), drive, *lo, *hi));
                    }
                    let (a, b) = (self.g.sub(v, *lo), self.g.sub(*hi, v));
                    let claim = Claim::Value(self.g.min(a, b), Cmp::Ge, 0.0);
                    let text = format!("{lt} <= {vt} <= {ht}, as value = min({vt} - {lt}, {ht} - {vt}) >= 0");
                    let terms = vec![(term.into(), v), ("lo".into(), *lo), ("hi".into(), *hi)];
                    self.asserts.push(AssertDef {
                        name: format!("{kind}:{asm}/{iname}"),
                        claim,
                        text,
                        bound: None,
                        res: None,
                        accuracy: tol,
                        quantity,
                        file: file.clone(),
                        line,
                        terms,
                        fixed: Vec::new(),
                    });
                }
                let mut clauses = Vec::new();
                for (package, head, x) in owned {
                    let mut r = Reader::in_joint(self, self.arg_env(), turn);
                    let kept = package.clause(&mut r, x, head, &x.list().unwrap()[1..])?;
                    let reads = r.into_reads();
                    clauses.push(OwnedClause { package: package.name(), head, kept, reads, file: file.clone(), line });
                }
                let [limits, speed_limits] = bounds.map(|b| b.map(|l| (l.0, l.1)));
                Some(Joint { turn, axis, anchor, below, drive, text, limits, speed_limits, clauses })
            }
            None => None,
        };
        let met = cframe.map(|(n, _)| n);
        Ok((self.compose(&at, &anchor), joint, on, met))
    }

    /// `(frame NAME F)` in a body or an assembly: F in its coordinates.
    fn frame_clause(&mut self, c: &Sexp, frames: &[(String, Xf)]) -> R<(String, Xf)> {
        let [_, n, f] = c.list().unwrap_or(&[]) else { return Err(self.err(c, "usage: (frame name F), F a frame such as (at (move 0 0 5))")) };
        let n = self.name_of(n)?;
        if frames.iter().any(|(m, _)| *m == n) {
            return Err(self.err(c, format!("duplicate frame '{n}'")));
        }
        Ok((n, self.frame_value(f, &mut self.arg_env(), "F")?))
    }

    /// `e` as a placement (`what` names the slot in the error): a frame
    /// whose numbers read params, t and local names, not x, y or z, however
    /// it was made. Every placement is read here: an instance's WHERE, a
    /// frame clause, `place`.
    fn frame_value(&mut self, e: &Sexp, env: &mut Env, what: &str) -> R<Xf> {
        let how = "a placement is a frame: (at STEP ...) with steps (move dx dy dz) (rot-x deg) (rot-y deg) (rot-z deg) (rot ax ay az deg), or (frame-of inst name)";
        match self.eval(e, env).map_err(|err| self.misplaced(err, e, how))? {
            Value::Frame(f) => {
                let numbers: Vec<NodeId> = f.r.iter().flatten().chain(&f.t).copied().collect();
                if self.g.reachable(&numbers).iter().any(|&n| matches!(self.g.nodes[n as usize], graph::Node::Var(Var::X | Var::Y | Var::Z))) {
                    return Err(self.err(e, "placements cannot depend on x, y or z"));
                }
                Ok(f)
            }
            other => Err(self.err_hint(e, format!("{what} must be a frame, not a {}", other.kind()), Some("e.g. (at (move 0 0 5)) or (frame-of inst name)".into()))),
        }
    }

    /// `err` from evaluating `e`, written where a placement or a step of
    /// one goes: a call of std's (a shape's move or turn) or of no function
    /// that failed as written (`(move 0 0 5)`, `(rotate-z 90)`, `(rotz
    /// 90)`), with no hint of its own, says how one is written (`how`). A
    /// call of the model's own function keeps its error as it is.
    fn misplaced(&self, mut err: InkError, e: &Sexp, how: &str) -> InkError {
        let Some(head) = e.list().and_then(|l| l.first()) else { return err };
        let call = head.sym().is_some_and(|h| {
            let key = self.scoped(h);
            match self.fns.get(&*key) {
                Some(f) => &*f.file == "std",
                None => self.def(&key).is_none() && !SPECIAL.contains(&h) && builtin(h).is_none() && !COMPARISONS.contains(&h),
            }
        });
        if call && err.hint.is_none() && err.file == *self.file && [e, head].iter().any(|s| (s.line, s.col) == (err.line, err.col)) {
            err.hint = Some(how.into());
        }
        err
    }

    /// `(frame-of PATH [NAME])`: the frame NAME declared by what PATH ends
    /// at, placed; NAME omitted, the placement itself. In an assembly's
    /// clauses a path starting at one of its instances is in its
    /// coordinates; otherwise PATH is a body (its coordinates) or a path
    /// from a top-level assembly (world coordinates).
    fn frame_of(&mut self, path: &str, name: Option<&str>) -> Result<Xf, String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or_default();
        let mut xf = self.identity();
        let mut child = self.resolve(first);
        let mut on = false;
        if let Some(b) = &mut self.building {
            if let Some(i) = b.asm.instances.iter().find(|i| i.name == first) {
                if self.own_frame {
                    return Err(format!(
                        "a relation end's FRAME is in its part's own coordinates, but (frame-of {path} ...) is in assembly '{}''s. Name a frame the part declares, (PART NAME), or write one, (PART (at ...))",
                        b.asm.name
                    ));
                }
                (xf, child) = (i.at, i.child.clone());
                on = true;
            } else if b.later.iter().any(|l| l == first) {
                return Err(format!("'{first}' is placed later in '{}': move its clause above this one", b.asm.name));
            }
        }
        if self.body(&child).is_none() && self.assembly(&child).is_none() {
            let mut names: Vec<&str> = self.bodies.iter().map(|b| b.name.as_str()).chain(self.assemblies.iter().map(|a| a.name.as_str())).collect();
            if let Some(b) = &self.building {
                names.extend(b.asm.instances.iter().map(|i| i.name.as_str()));
            }
            let hint = closest(first, names.into_iter()).unwrap_or_default();
            return Err(format!("no instance, body or assembly named '{first}'. {hint}").trim().to_string());
        }
        for seg in segs {
            let Some(sub) = self.assembly(&child) else { return Err(format!("'{child}' is a body; '{path}' goes through it")) };
            let Some(i) = sub.instances.iter().find(|i| i.name == seg) else {
                let hint = choices(seg, "Its instances", sub.instances.iter().map(|i| i.name.as_str()));
                return Err(format!("'{child}' has no instance '{seg}'. {hint}"));
            };
            let (at, c) = (i.at, i.child.clone());
            xf = self.compose(&at, &xf);
            child = c;
        }
        // What the instance being built is placed on: the instance under
        // PATH the frame is of (an assembly's frame of one of its own
        // instances: a servo's output frame is its output's).
        let under: Vec<String> = name.and_then(|n| self.assembly(&child)?.frame_on.iter().find(|(f, _)| f == n).map(|(_, on)| on.clone())).unwrap_or_default();
        if on {
            let read = &mut self.building.as_mut().unwrap().read;
            match under.is_empty() {
                true => read.push(path.to_string()),
                false => read.extend(under.iter().map(|u| format!("{path}/{u}"))),
            }
        }
        let Some(name) = name else { return Ok(xf) };
        let frames = self.frames_of(&child);
        match frames.iter().find(|(n, _)| n == name) {
            Some((_, f)) => Ok(self.compose(f, &xf)),
            None => Err(Self::no_frame(&child, name, &frames)),
        }
    }

    /// The frames a body or an assembly declares.
    fn frames_of(&self, child: &str) -> Vec<(String, Xf)> {
        self.body(child).map(|b| b.frames.clone()).or(self.assembly(child).map(|a| a.frames.clone())).unwrap_or_default()
    }

    fn no_frame(child: &str, name: &str, frames: &[(String, Xf)]) -> String {
        if frames.is_empty() {
            return format!("'{child}' has no frame '{name}' (it declares none). Add (frame {name} F) to it, or place it by WHERE alone");
        }
        let names: Vec<&str> = frames.iter().map(|(n, _)| n.as_str()).collect();
        format!("'{child}' has no frame '{name}'. Its frames: {}", names.join(", "))
    }

    /// The inverse of a rigid placement: (R^T, -R^T t).
    pub fn inverse(&mut self, a: &Xf) -> Xf {
        let mut r = a.r;
        let mut t = a.t;
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = a.r[j][i];
            }
        }
        for i in 0..3 {
            let mut acc = self.g.constant(0.0);
            for k in 0..3 {
                let p = self.g.mul(r[i][k], a.t[k]);
                acc = self.g.add(acc, p);
            }
            t[i] = self.g.un(UnOp::Neg, acc);
        }
        Xf { r, t }
    }

    /// Contact ends are instance names, or paths into sub-assemblies
    /// (`inst/sub-inst/...`), relative to the assembly being defined.
    fn check_contact_path(&self, asm: &Assembly, path: &str) -> Result<(), String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or("");
        let Some(inst) = asm.instances.iter().find(|i| i.name == first) else {
            let hint = closest(first, asm.instances.iter().map(|i| i.name.as_str())).unwrap_or_default();
            return Err(format!("no instance '{first}' (read so far). {hint}").trim().to_string());
        };
        let mut child = inst.child.clone();
        for seg in segs {
            let sub = self.assembly(&child).ok_or_else(|| format!("'{child}' is a body; '{path}' goes too deep"))?;
            let i = sub
                .instances
                .iter()
                .find(|i| i.name == seg)
                .ok_or_else(|| format!("'{child}' has no instance '{seg}'. {}", choices(seg, "Its instances", sub.instances.iter().map(|i| i.name.as_str()))))?;
            child = i.child.clone();
        }
        Ok(())
    }

    /// Every relation of what `target` places (Doc::target): those
    /// declared in it and under it, and those an enclosing assembly
    /// declares between parts of it. Its ends are full paths (matching
    /// `Item::label`) and its zone in world coordinates: A's part placed
    /// where A is, B's where B is.
    pub fn relations_under(&mut self, target: &str) -> Vec<PlacedRelation> {
        match self.target(target) {
            Some(t) => self.relations_of(&t),
            None => Vec::new(),
        }
    }

    /// The relations of a resolved target (Doc::relations_under).
    pub fn relations_of(&mut self, t: &Target) -> Vec<PlacedRelation> {
        let mut out = Vec::new();
        // An end is the target's where it is at or under it, or above it
        // (an enclosing instance holds all of it).
        let within = |end: &str| end == t.path || end.starts_with(&format!("{}/", t.path)) || t.path.starts_with(&format!("{end}/"));
        for (path, asm, xf) in &t.above {
            let mut declared = Vec::new();
            self.relations_here(asm, path, xf, &mut declared);
            out.extend(declared.into_iter().filter(|r| within(&r.a) && within(&r.b)));
        }
        self.collect_relations(&t.def, &t.path, &t.at, &mut out);
        out
    }
    /// Every clause `package` owns written in an assembly placed under
    /// `target` (a top-level assembly, or a path into one), each with the
    /// path of its assembly's instance and where that is placed in the world.
    pub fn clauses_under(&mut self, target: &str, package: &str) -> Vec<PlacedClause> {
        let mut out = Vec::new();
        let Some((xf, asm_name)) = self.placement(target) else { return out };
        self.collect_clauses(&asm_name, target, &xf, package, &mut out);
        out
    }

    fn collect_clauses(&mut self, asm: &str, prefix: &str, xf: &Xf, package: &str, out: &mut Vec<PlacedClause>) {
        let Some(a) = self.assembly(asm) else { return };
        out.extend(a.clauses.iter().filter(|c| c.package == package).map(|c| PlacedClause { path: prefix.to_string(), at: *xf, clause: c.clone() }));
        let insts: Vec<(String, String, Xf)> = a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at)).collect();
        for (name, child, at) in insts {
            let xf = self.compose(&at, xf);
            self.collect_clauses(&child, &format!("{prefix}/{name}"), &xf, package, out);
        }
    }

    /// Where an instance path (`asm/inst/...`) is placed in its top
    /// assembly's coordinates, its instances composed from the top (as
    /// `items` places what is under it), and what it places. A top-level
    /// name: the identity and itself. None: no such path.
    pub fn placement(&mut self, path: &str) -> Option<(Xf, String)> {
        self.target(path).map(|t| (t.at, t.def))
    }

    /// The body or assembly `path` names (a definition name, or an
    /// instance path from a top-level assembly), where it is placed, the
    /// file that defines it and the assemblies enclosing it. None: no
    /// such body, assembly or path.
    pub fn target(&mut self, path: &str) -> Option<Target> {
        let mut segs = path.split('/');
        let mut def = segs.next().unwrap_or_default().to_string();
        let (mut at, mut above, mut here) = (self.identity(), Vec::new(), def.clone());
        for seg in segs {
            let (i, next) = self.assembly(&def).and_then(|a| a.instances.iter().find(|i| i.name == seg)).map(|i| (i.at, i.child.clone()))?;
            above.push((here.clone(), def, at));
            at = self.compose(&i, &at);
            here = format!("{here}/{seg}");
            def = next;
        }
        let file = self.assembly(&def).map(|a| a.file.clone()).or_else(|| self.body(&def).map(|b| b.file.clone()))?;
        Some(Target { path: path.to_string(), def, file, at, above })
    }

    fn collect_relations(&mut self, asm: &str, prefix: &str, xf: &Xf, out: &mut Vec<PlacedRelation>) {
        self.relations_here(asm, prefix, xf, out);
        let Some(a) = self.assembly(asm) else { return };
        let insts: Vec<(String, String, Xf)> = a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at)).collect();
        for (name, child, at) in insts {
            let xf = self.compose(&at, xf);
            self.collect_relations(&child, &format!("{prefix}/{name}"), &xf, out);
        }
    }

    /// The relations assembly `asm` itself declares, placed at `xf` under
    /// path `prefix`.
    fn relations_here(&mut self, asm: &str, prefix: &str, xf: &Xf, out: &mut Vec<PlacedRelation>) {
        let Some(a) = self.assembly(asm) else { return };
        let relations: Vec<(String, String, Option<Zone>, Option<Rc<Kind>>)> = a.relations.iter().map(|r| (r.a.clone(), r.b.clone(), r.zone, r.kind.clone())).collect();
        for (pa, pb, zone, kind) in relations {
            // A field in an end's coordinates, placed: its path's
            // instances, in order.
            let placed = |doc: &mut Self, path: &str, f: NodeId| {
                let (mut at, mut child) = (*xf, asm.to_string());
                for seg in path.split('/') {
                    let i = doc.assembly(&child).and_then(|a| a.instances.iter().find(|i| i.name == seg)).map(|i| (i.at, i.child.clone())).unwrap();
                    at = doc.compose(&i.0, &at);
                    child = i.1;
                }
                let Value::Field(z, _) = doc.to_world(&Value::Field(f, None), &at) else { unreachable!() };
                z
            };
            let placed_zone = zone.map(|z| {
                let w = placed(self, &pa, z.a);
                match z.b {
                    Some(b) => {
                        let wb = placed(self, &pb, b);
                        self.g.max(w, wb)
                    }
                    None => w,
                }
            });
            out.push(PlacedRelation { a: format!("{prefix}/{pa}"), b: format!("{prefix}/{pb}"), zone: placed_zone, local: zone, kind, at: *xf });
        }
    }

    /// The cell size of a `(res h)` clause (an assert's, a sampled field's,
    /// a study's): a number known at load time, above 0 and finite, as for
    /// `--res`.
    fn cell_size(&mut self, e: &Sexp, env: &mut Env) -> R<f64> {
        let f = self.field_of(e, Quantity::Length, "res", env)?;
        let h = self.const_of(e, f)?;
        self.check_cell(e, h)
    }

    /// A warning where `sin`, `cos` or `tan` (radians) is given what reads
    /// as degrees: a whole multiple of 15 from 15 up, alone or a factor of
    /// a product (`(sin 30)`, `(sin (* 360 t))`).
    fn degrees_warning(&mut self, name: &str, x: &Sexp) {
        let degrees = |e: &Sexp| matches!(e.k, SK::Num(v, None) if v.abs() >= 15.0 && v % 15.0 == 0.0);
        let seen = match x.head_args() {
            Some(("*", factors)) => factors.iter().find(|f| degrees(f)),
            _ => Some(x).filter(|x| degrees(x)),
        };
        if let Some(d) = seen {
            let w = format!("{}:{}: ({name} {x}) takes radians, but {d} looks like degrees. If so, write ({name} (deg->rad {x}))", self.file, x.line);
            if !self.load_warnings.contains(&w) {
                self.load_warnings.push(w);
            }
        }
    }

    /// Err where value `v` of expression `e`, in a slot that takes
    /// quantity `q` (`what`: the slot, for the message), was said in
    /// another quantity's unit (units::Said), or a vec of it has one that was.
    fn fits(&self, v: &Value, e: &Sexp, q: Quantity, what: &str) -> R<()> {
        match said(v).into_iter().find_map(|s| s.misfit(q, what)) {
            Some(m) => Err(self.err(e, m)),
            None => Ok(()),
        }
    }

    /// Value `v` of `e` is a change of a quantity (an accuracy, a
    /// tolerance's deviation): a unit with a zero of its own (degF, K)
    /// reads a point on its scale (1degF is -17.2 degC), not a change, so
    /// it is refused.
    fn fits_change(&self, v: &Value, e: &Sexp, what: &str) -> R<()> {
        let zero = said(v).into_iter().find_map(|s| match s.unit.scale {
            units::Scale::Offset { factor, .. } => Some((s, factor)),
            _ => None,
        });
        match zero {
            Some((s, k)) => Err(self.err_hint(
                e,
                format!("{what} is a change, and {} reads a point on its scale (written {})", s.unit.name, s.cite()),
                Some(format!("write the change in {}: one of 1{} is {}{}", s.unit.quantity.unit(), s.unit.name, (k * 1e6).round() / 1e6, s.unit.quantity.unit())),
            )),
            None => Ok(()),
        }
    }

    /// What param `name`'s value is said by: the unit `--param` gave it in
    /// (`overrides`), which must be the quantity its default (`said`, at
    /// `at`) was said in, else the default's.
    fn param_given(&self, name: &str, said: Option<Rc<units::Said>>, at: &Sexp) -> R<Option<Rc<units::Said>>> {
        let Some((_, Given { value, unit: Some(unit) })) = self.overrides.iter().find(|(k, _)| k == name) else { return Ok(said) };
        let (v, unit) = (*value, *unit);
        if let Some(s) = &said
            && s.unit.quantity != unit.quantity
        {
            return Err(self.err(at, format!("--param {name}: {name} is {}, written {}, but {} is {}", s.unit.quantity.name(), s.cite(), unit.name, unit.quantity.name())));
        }
        Ok(Some(Rc::new(units::Said { unit, written: units::to(v, unit), file: std::sync::Arc::from(format!("--param {name}")), line: 0 })))
    }

    /// Expression `e` as a field where a slot takes quantity `q` (`what`:
    /// the slot), its value checked as `fits` does.
    fn field_of(&mut self, e: &Sexp, q: Quantity, what: &str, env: &mut Env) -> R<NodeId> {
        let v = self.eval(e, env)?;
        self.fits(&v, e, q, what)?;
        self.as_field(v, e)
    }

    /// What a sum's value `(name ARGS)` is said by (units::Said): its
    /// terms' one quantity, an error where they say two. Other ops carry
    /// none (a product's units combine; a comparison is 1 or 0).
    fn sum_said(&self, name: &str, args: &[Value], at: &Sexp) -> R<Option<Rc<units::Said>>> {
        let terms = match name {
            "+" | "-" | "neg" | "abs" | "min" | "max" => args,
            "select" if args.len() == 3 => &args[1..],
            _ => return Ok(None),
        };
        self.one_said(name, terms, at)
    }

    /// What values that are one quantity (a sum's terms, a range's ends)
    /// are said by: the first said, an error at `at` where two say two
    /// quantities (`what`: the form, for the message).
    fn one_said(&self, what: &str, terms: &[Value], at: &Sexp) -> R<Option<Rc<units::Said>>> {
        let mut said = terms.iter().filter_map(|v| if let Value::Field(_, Some(s)) = v { Some(s) } else { None });
        let Some(first) = said.next() else { return Ok(None) };
        let mut kept = first;
        for s in said {
            if s.unit.quantity != first.unit.quantity {
                let m = format!("'{what}' takes values of one quantity, but {} is {} and {} is {}", first.cite(), first.unit.quantity.name(), s.cite(), s.unit.quantity.name());
                return Err(self.err(at, m));
            }
            // Said in the unit units::kept keeps (degF over degC).
            if units::kept(kept.unit, s.unit).is_ok_and(|u| !std::ptr::eq(u, kept.unit)) {
                kept = s;
            }
        }
        Ok(Some(kept.clone()))
    }

    fn check_cell(&self, e: &Sexp, h: f64) -> R<f64> {
        if !(h > 0.0 && h.is_finite()) {
            return Err(self.err(e, format!("(res h) expects a cell size above 0 and finite, got {h}")));
        }
        Ok(h)
    }

    /// An assert's options, the clauses after its claim: `(res r)`, the
    /// cell size the claim needs (a constant above 0 and finite, as for
    /// `--res`), and `(accuracy x)`, the accuracy it is checked to, in its
    /// units (`accuracy_clause`), for the claims that have one. Returns them.
    fn assert_options(&mut self, c: &Claim, options: &[Sexp]) -> R<(Option<f64>, Option<f64>)> {
        const USAGE: &str = "assert options: (res r), the cell size the claim needs; (accuracy x), the accuracy it is checked to";
        let mut res = None;
        let mut accuracy = None;
        for o in options {
            match o.list() {
                Some([head, _]) if head.sym() == Some("accuracy") => {
                    if accuracy.is_some() {
                        return Err(self.err(o, "assert has more than one (accuracy ...)"));
                    }
                    if c.accuracy() == Accuracy::Exact {
                        return Err(self.err_hint(
                            o,
                            "this claim is decided exactly, and has no accuracy to state",
                            Some("a point claim is decided by its sign there, and a range compares values exactly. Write any slack into its bounds".into()),
                        ));
                    }
                    let q = self.accuracy_quantity(c);
                    accuracy = Some(self.accuracy_clause(o, q, &mut Vec::new())?);
                }
                Some([head, r]) if head.sym() == Some("res") => {
                    if res.is_some() {
                        return Err(self.err(o, "assert has more than one (res ...)"));
                    }
                    if matches!(c, Claim::Metric { metric, .. } if metric.shape.is_none()) {
                        return Err(self.err_hint(
                            o,
                            "(res r) is the cell size of a geometry claim, and this claim is not checked on cells",
                            Some("a package form's metric is solved at its form's own (res ...)".into()),
                        ));
                    }
                    // A constant: a claim's res is not a design parameter.
                    let f = self.field_of(r, Quantity::Length, "res", &mut Vec::new())?;
                    let h = self.g.as_const(f).ok_or_else(|| self.err(r, "expected a constant number"))?;
                    res = Some(self.check_cell(r, h)?);
                }
                Some([head, _]) if head.sym() == Some("tol") => {
                    return Err(self.err_hint(o, "(tol ...) is a param's tolerance. A claim states its accuracy: (accuracy x)", Some(USAGE.into())));
                }
                _ => return Err(self.err_hint(o, format!("unknown assert option {o}"), Some(USAGE.into()))),
            }
        }
        Ok((res, accuracy))
    }

    /// What a claim's accuracy (and a comparison's value) is in, where
    /// known: a length for a claim about shapes or a thickness, a volume's
    /// volume, a metric's quantity. None: the model's own units (a value,
    /// an integral, a mean), or none at all (points and ranges).
    pub fn accuracy_quantity(&self, c: &Claim) -> Option<Quantity> {
        match c {
            _ if c.accuracy_is_length() => Some(Quantity::Length),
            Claim::Metric { metric, .. } => self.registry.metric_quantity(metric.metric),
            _ => None,
        }
    }

    /// An `(accuracy x)` clause (an assert's, a joint's limits', a
    /// relation's claim's): how far the claim's answer may be off, a
    /// constant above 0 and finite, in quantity `q` where the claim's units
    /// are known. `(tol ...)` is a part's tolerance, a param's band.
    fn accuracy_clause(&mut self, o: &Sexp, q: Option<Quantity>, env: &mut Env) -> R<f64> {
        const HINT: &str = "(accuracy x): how far the claim's answer may be off, in its units";
        let x = match o.list() {
            Some([head, x]) if head.sym() == Some("accuracy") => x,
            Some([head, _]) if head.sym() == Some("tol") => {
                return Err(self.err_hint(o, "(tol ...) is a param's tolerance. A claim states its accuracy instead", Some(HINT.into())));
            }
            _ => return Err(self.err_hint(o, format!("expected (accuracy x), not {o}"), Some(HINT.into()))),
        };
        let v = self.eval(x, env)?;
        if let Some(q) = q {
            self.fits(&v, x, q, "accuracy")?;
        }
        self.fits_change(&v, x, "an accuracy")?;
        let f = self.as_field(v, x)?;
        let t = self.g.as_const(f).ok_or_else(|| self.err(x, "expected a constant number"))?;
        if !(t > 0.0 && t.is_finite()) {
            return Err(self.err(x, format!("(accuracy x) expects an accuracy above 0 and finite, got {t}")));
        }
        Ok(t)
    }

    /// A claim, and its bound as a field of params (see `claim_at`).
    fn claim(&mut self, c: &Sexp) -> R<(Claim, Option<NodeId>)> {
        let claim = self.claim_parse(c);
        let bound = self.bound.take();
        Ok((claim?, bound))
    }

    fn claim_parse(&mut self, c: &Sexp) -> R<Claim> {
        let usage = format!(
            "claims: (empty f) (disjoint a b) (subset a b) (within f (x0 y0 z0) (x1 y1 z1)) \
             (clearance a b d) (touches a b) (range f lo hi region) (volume f <|<=|>|>= v) \
             (value expr < v) (overlaps a b) (integral target g <|<=|>|>= v) \
             (mean target h [(weight w)] <|<=|>|>= v), and a package's metric of a form or a shape, (metric form <|<=|>|>= v): {}",
            self.registry.metrics()
        );
        let (head, a) = c.head_args().ok_or_else(|| self.err_hint(c, "a claim must be a list", Some(usage.clone())))?;
        let env = &mut Vec::new();
        let bad = |d: &Doc| d.err_hint(c, format!("malformed claim '{head}'"), Some(usage.clone()));
        // Bounds are numbers: constants and params (at their current values),
        // like load magnitudes and clearance distances.
        const BOUND: &str = "a claim's bound must be a number (constants and params), not a field of x y z or t";
        Ok(match (head, a) {
            ("empty", [f]) => Claim::Empty(self.field(f, env)?),
            // A point overlaps a shape where it is inside it.
            ("overlaps", [p, q]) => match (self.shape_or_point(p, env)?, self.shape_or_point(q, env)?) {
                (Ok(p), Ok(q)) => Claim::Overlaps(p, q),
                (Err(at), Ok(f)) | (Ok(f), Err(at)) => Claim::Point { f, at, inside: true },
                (Err(_), Err(_)) => return Err(self.err(c, "overlaps takes a shape and a shape or a point, not two points")),
            },
            ("touches", [p, q]) => {
                let (p, q) = (self.field(p, env)?, self.field(q, env)?);
                Claim::Touches(p, q)
            }
            ("value", [f, op, v]) => {
                let f = self.field(f, env)?;
                if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                    return Err(self.err(c, "value claims take numbers (params, t), not fields of x y z"));
                }
                let cmp = self.cmp_of(op)?;
                let b = self.field(v, env)?;
                let b = self.claim_bound(b).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Value(f, cmp, b)
            }
            // A point (std's probe) and a shape: the point outside it, or
            // inside it.
            ("disjoint", [p, q]) => match (self.shape_or_point(p, env)?, self.shape_or_point(q, env)?) {
                (Ok(p), Ok(q)) => Claim::Disjoint(p, q),
                (Err(at), Ok(f)) | (Ok(f), Err(at)) => Claim::Point { f, at, inside: false },
                (Err(_), Err(_)) => return Err(self.err(c, "disjoint takes a shape and a shape or a point, not two points")),
            },
            ("subset", [p, q]) => match (self.shape_or_point(p, env)?, self.field(q, env)?) {
                (Ok(p), q) => Claim::Subset(p, q),
                (Err(at), f) => Claim::Point { f, at, inside: true },
            },
            ("within", [f, lo, hi]) => {
                let f = self.field(f, env)?;
                // Corners of constants and params: the box stays a field of
                // its params, so a ranged one is proven over its band.
                let (lo, hi) = (self.corner(lo)?, self.corner(hi)?);
                let b = self.box_nodes(lo, hi);
                Claim::Subset(f, b)
            }
            ("clearance", [p, q, d]) => {
                let (p, q) = (self.field(p, env)?, self.field(q, env)?);
                let d = self.field_of(d, Quantity::Length, "clearance", env)?;
                let d = self.claim_bound(d).ok_or_else(|| self.err(c, "clearance distance must be a number (constants and params), not a field"))?;
                Claim::Clearance(p, q, d)
            }
            ("range", [f, lo, hi, region]) => {
                let (f, lo, hi, region) = (self.field(f, env)?, self.field(lo, env)?, self.field(hi, env)?, self.field(region, env)?);
                let below = self.g.sub(f, lo);
                let above = self.g.sub(hi, f);
                let violation = self.g.min(below, above);
                Claim::Range(self.g.max(region, violation))
            }
            ("volume", [f, op, v]) => {
                let f = self.field(f, env)?;
                let cmp = self.cmp_of(op)?;
                let v = self.field_of(v, Quantity::Volume, "volume", env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                // The integral of 1 over the shape.
                let one = self.g.constant(1.0);
                Claim::Integral { parts: vec![(f, vec![one])], frames: vec![None], cmp, bound }
            }
            ("integral", [target, g, op, v]) => {
                let (parts, frames) = self.integrands(target, &[g], env)?;
                let cmp = self.cmp_of(op)?;
                let v = self.field(v, env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Integral { parts: parts.into_iter().map(|(s, mut g)| (s, g.remove(0))).collect(), frames, cmp, bound }
            }
            ("mean", [target, h, rest @ ..]) if matches!(rest.len(), 2 | 3) => {
                let (w, op, v) = match rest {
                    [op, v] => (None, op, v),
                    [w, op, v] => match w.list() {
                        Some([head, w]) if head.sym() == Some("weight") => (Some(w), op, v),
                        _ => {
                            return Err(self.err_hint(
                                w,
                                "expected (weight w)",
                                Some("(mean target h (weight w) < v): the mean of h weighted by w, e.g. (weight (attr density)) for the centre of mass".into()),
                            ));
                        }
                    },
                    _ => unreachable!(),
                };
                let one = Sexp { k: SK::Num(1.0, None), ..h.clone() };
                let (parts, frames) = self.integrands(target, &[w.unwrap_or(&one), h], env)?;
                if parts.iter().any(|(_, f)| f[0].len() != 1) {
                    return Err(self.err(c, "a mean's weight is one field"));
                }
                let cmp = self.cmp_of(op)?;
                let v = self.field(v, env)?;
                let bound = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Mean { parts: parts.into_iter().map(|(s, mut f)| (s, f[0][0], f.remove(1))).collect(), frames, cmp, bound }
            }
            // After the kernel's kinds, so a package cannot shadow one.
            (m, [of, op, v, rest @ ..]) if self.registry.metric(m).is_some() => {
                // (excluding KEY ...): what the metric leaves out that the
                // claim leaves out on purpose.
                let mut excluding = Vec::new();
                let known = self.registry.exclusions(m);
                let usage = || {
                    let list: Vec<String> = known.iter().map(|(k, about)| format!("{k}: {about}")).collect();
                    if list.is_empty() { format!("'{m}' leaves nothing out") } else { format!("'{m}' may leave out {}", list.join("; ")) }
                };
                for x in rest {
                    let Some([head, keys @ ..]) = x.list() else { return Err(self.err_hint(x, "expected (excluding KEY ...)", Some(usage()))) };
                    if head.sym() != Some("excluding") || keys.is_empty() {
                        return Err(self.err_hint(x, "expected (excluding KEY ...)", Some(usage())));
                    }
                    for k in keys {
                        match known.iter().find(|e| Some(e.0) == k.sym()) {
                            Some(e) => excluding.push(e.0),
                            None => return Err(self.err_hint(k, format!("'{m}' leaves out nothing called {}", k.sym().unwrap_or("that")), Some(usage()))),
                        }
                    }
                }
                let metric = self.form_metric(m, of, env, "assert")?;
                let cmp = self.cmp_of(op)?;
                let v = match self.registry.metric_quantity(m) {
                    Some(q) => self.field_of(v, q, m, env)?,
                    None => self.field(v, env)?,
                };
                let value = self.claim_bound(v).ok_or_else(|| self.err(c, BOUND))?;
                Claim::Metric { metric, cmp, value, excluding }
            }
            _ => return Err(bad(self)),
        })
    }

    /// A claim's bound: its value at load time; its field is returned with
    /// the claim (see `claim`), and is what claims are checked against
    /// (`claim_at`), so its params are not materialized, nor baked: an
    /// optimizer differentiates the bound's field with the claim.
    fn claim_bound(&mut self, f: NodeId) -> Option<f64> {
        self.bound = Some(f);
        if let Some(v) = self.g.as_const(f) {
            return Some(v);
        }
        if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
            return None;
        }
        Some(self.load_eval(f))
    }

    /// An integral's parts: for each body of `target` (a body, an assembly,
    /// or an instance path: its placed bodies, in world space), its shape
    /// and the fields of each of `exprs` (a field, or a vec of them), read
    /// with `(attr NAME)` meaning that body's attribute, placed with it. A
    /// target that names none of those is a field, one part with no
    /// attributes; one with no bodies (an empty assembly) is one part of
    /// an empty shape, its fields read once with no attributes: each
    /// integral is 0, and how many there are is the fields'. Every body
    /// gives each expression the same number of fields (one integral per
    /// field): an error names two that differ. With each part, its body's
    /// own frame, where it is a placed body.
    fn integrands(&mut self, target: &Sexp, exprs: &[&Sexp], env: &mut Env) -> R<(Vec<(NodeId, Vec<Vec<NodeId>>)>, Vec<Option<Placed>>)> {
        let path = target.sym().map(|s| self.resolve(s)).filter(|s| s.contains('/') || self.body(s).is_some() || self.assembly(s).is_some());
        let mut items: Vec<(NodeId, String, String, Result<Attrs, &'static str>, Option<Placed>)> = match path {
            Some(path) => {
                let items = self.items(&path).map_err(|m| self.err(target, m))?;
                items
                    .into_iter()
                    .map(|it| {
                        let who = self.named(&it);
                        let (shape, attrs) = (self.placed_field(&it), self.placed_attrs(&it));
                        let frame = it.instance.map(|(_, at)| Placed { shape: it.shape, at });
                        (shape, it.label, who, Ok(attrs), frame)
                    })
                    .collect()
            }
            None => vec![(self.field(target, env)?, target.to_string(), format!("'{target}'"), Err("is a field, with no attributes"), None)],
        };
        if items.is_empty() {
            let empty = self.g.constant(1.0);
            items.push((empty, target.to_string(), format!("'{target}'"), Err("has no bodies, so no attributes"), None));
        }
        let mut out: Vec<(NodeId, Vec<Vec<NodeId>>)> = Vec::new();
        let mut frames = Vec::new();
        let mut first = String::new();
        for (shape, label, who, attrs, frame) in items {
            let outer = self.integrand_of.replace((who, attrs));
            let fields = exprs.iter().map(|e| self.fields(e, env)).collect::<R<Vec<_>>>();
            self.integrand_of = outer;
            let fields = fields?;
            match out.first() {
                None => first = label,
                Some((_, those)) => {
                    let count = |n: usize| if n == 1 { "1 field".to_string() } else { format!("{n} fields") };
                    if let Some(((e, f), f0)) = exprs.iter().zip(&fields).zip(those).find(|((_, f), f0)| f.len() != f0.len()) {
                        return Err(self.err_hint(
                            e,
                            format!("{e} is {} for '{first}' but {} for '{label}'", count(f0.len()), count(f.len())),
                            Some("each field is one integral, summed over the bodies: every body must give the same number".into()),
                        ));
                    }
                }
            }
            out.push((shape, fields));
            frames.push(frame);
        }
        Ok((out, frames))
    }

    /// A field, or each of a vec of them (at least one).
    fn fields(&mut self, e: &Sexp, env: &mut Env) -> R<Vec<NodeId>> {
        match self.eval(e, env)? {
            Value::Field(f, _) => Ok(vec![f]),
            Value::Vec(v) if !v.is_empty() => v
                .into_iter()
                .map(|x| match x {
                    Value::Field(f, _) => Ok(f),
                    other => Err(self.err(e, format!("expected fields, got a {} in the vec: {e}", other.kind()))),
                })
                .collect(),
            other => Err(self.err(e, format!("expected a field or a vec of fields, got a {}: {e}", other.kind()))),
        }
    }

    /// The claim with its bound evaluated at `ctx` (solved params at their
    /// solved values, wherever the solve is in the file): the one path for
    /// asserts and optimization requirements.
    pub fn claim_at(&self, c: &Claim, bound: Option<NodeId>, ctx: &Ctx) -> Claim {
        let at = |v: f64| bound.map_or(v, |b| self.scalar(b, ctx));
        match c {
            Claim::Empty(f) => Claim::Empty(*f),
            Claim::Range(f) => Claim::Range(*f),
            Claim::Disjoint(p, q) => Claim::Disjoint(*p, *q),
            Claim::Subset(p, q) => Claim::Subset(*p, *q),
            Claim::Touches(p, q) => Claim::Touches(*p, *q),
            Claim::Point { f, at, inside } => Claim::Point { f: *f, at: *at, inside: *inside },
            Claim::Clearance(p, q, d) => Claim::Clearance(*p, *q, at(*d)),
            Claim::Metric { metric, cmp, value, excluding } => Claim::Metric { metric: metric.clone(), cmp: *cmp, value: at(*value), excluding: excluding.clone() },
            Claim::Value(f, cmp, v) => Claim::Value(*f, *cmp, at(*v)),
            Claim::Overlaps(p, q) => Claim::Overlaps(*p, *q),
            Claim::Integral { parts, frames, cmp, bound: v } => Claim::Integral { parts: parts.clone(), frames: frames.clone(), cmp: *cmp, bound: at(*v) },
            Claim::Mean { parts, frames, cmp, bound: v } => Claim::Mean { parts: parts.clone(), frames: frames.clone(), cmp: *cmp, bound: at(*v) },
        }
    }

    /// A construction name as written, matched to the names in the graph.
    /// Inside a namespaced file a name may be renamed `ns.name` (it is also
    /// a def there) while the tag that made it, in a family file, kept the
    /// plain name, or the other way round: the exact name wins, then the
    /// name without, then with, the current namespace.
    fn construction_names(&self, name: &str) -> Vec<String> {
        let mut v = vec![name.to_string()];
        if let Some(ns) = &self.ns {
            match name.strip_prefix(&format!("{ns}.")) {
                Some(plain) => v.push(plain.to_string()),
                None => v.push(format!("{ns}.{name}")),
            }
        }
        v
    }

    fn resolve_tag(&self, name: &str) -> Option<(u32, String)> {
        self.construction_names(name).into_iter().find_map(|n| self.g.tag_names.iter().position(|t| *t == n).map(|k| (k as u32, n)))
    }

    fn cmp_of(&self, op: &Sexp) -> R<Cmp> {
        Ok(match op.sym() {
            Some("<") => Cmp::Lt,
            Some("<=") => Cmp::Le,
            Some(">") => Cmp::Gt,
            Some(">=") => Cmp::Ge,
            _ => return Err(self.err(op, "comparison must be one of < <= > >=")),
        })
    }

    /// The params these expressions read.
    pub fn params_read(&self, nodes: &[NodeId]) -> Vec<usize> {
        self.g
            .reachable(nodes)
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Var(Var::Param(k)) => Some(k as usize),
                _ => None,
            })
            .collect()
    }

    /// Where param k's value first became a number at load time (a loop
    /// count, a sampled or measured field, a solve's input...): file:line
    /// (head name). None: it stays symbolic, so it can be ranged.
    pub fn materialized(&self, k: usize) -> Option<String> {
        self.materialized.borrow().get(&k).cloned()
    }

    /// Params baked into the graph at load time (see `bake`), sorted.
    pub fn baked_params(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.baked.borrow().iter().copied().collect();
        v.sort();
        v
    }

    /// What Rust builds the model from is noted where outside ink's range
    /// (Graph::note): an axis of 1.8e-162 normalized squares to 0.
    fn load_value(&self, f: NodeId) -> Option<f64> {
        let v = match self.g.as_const(f) {
            Some(v) => v,
            None if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) => return None,
            None => {
                self.bake(&[f]);
                self.load_eval(f)
            }
        };
        self.g.note(v, crate::graph::wild(v));
        Some(v)
    }

    pub fn identity(&mut self) -> Xf {
        let (o, z) = (self.g.constant(1.0), self.g.constant(0.0));
        Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: [z, z, z] }
    }

    /// `a` then `b`: p -> Rb (Ra p + ta) + tb.
    pub fn compose(&mut self, a: &Xf, b: &Xf) -> Xf {
        let mut r = a.r;
        let mut t = a.t;
        for i in 0..3 {
            for j in 0..3 {
                let mut acc = self.g.constant(0.0);
                for k in 0..3 {
                    let p = self.g.mul(b.r[i][k], a.r[k][j]);
                    acc = self.g.add(acc, p);
                }
                r[i][j] = acc;
            }
            let mut acc = b.t[i];
            for k in 0..3 {
                let p = self.g.mul(b.r[i][k], a.t[k]);
                acc = self.g.add(acc, p);
            }
            t[i] = acc;
        }
        Xf { r, t }
    }

    /// The unit vector along `v`: v / |v|. A direction must have a length,
    /// or the division by zero makes it 0 (a rotation about it a uniform
    /// scale). Where `v` reads no x, y, z (a placement's axis, a frame's
    /// direction), a zero length is refused: a constant one here, at
    /// `at`; otherwise its squared length is a condition (`conditions`),
    /// checked once the document is read at the params' values, at each
    /// moment a view is made at, and proven positive over every scope it
    /// reads (Doc::check_conditions). A `v` over space is a field: its
    /// zero is an undefined case, reported at samples like any division.
    fn unit(&mut self, v: &[NodeId], at: &Sexp, what: &str, hint: &str) -> R<Vec<NodeId>> {
        let g = &mut self.g;
        let sq: Vec<NodeId> = v.iter().map(|&a| g.un(UnOp::Square, a)).collect();
        let n2 = sq[1..].iter().fold(sq[0], |acc, &s| g.add(acc, s));
        if !self.g.depends_on(n2, |x| matches!(x, Var::X | Var::Y | Var::Z)) {
            match self.g.as_const(n2) {
                Some(c) if c > 0.0 && c.is_finite() => {}
                Some(c) => {
                    let why = if c > 0.0 { "is too long to square in doubles" } else { "has zero length" };
                    return Err(self.err_hint(at, format!("{what} {why}: it has no direction"), Some(hint.into())));
                }
                None if self.conditions.iter().any(|c| c.len2 == n2) => {}
                None => {
                    let c = Condition { len2: n2, what: what.to_string(), site: self.form_at.clone(), at: self.err_hint(at, "", Some(hint.into())) };
                    self.conditions.push(c);
                }
            }
        }
        let g = &mut self.g;
        let n = g.un(UnOp::Sqrt, n2);
        Ok(v.iter().map(|&a| g.bin(BinOp::Div, a, n)).collect())
    }

    /// Whether the directions (Doc::unit) have a length at ctx: at its
    /// moment and params' values, and, where they read them, for every
    /// value of its ranged params (ranged::value_span, as a value claim)
    /// and every moment of `motion` (value_span with t ranged over it). Err: the first that does not, naming it and where it is
    /// zero (a witness), or that its length could not be proven nonzero
    /// over the scope. `at_load`: only those not depending on t (load-time
    /// values are at t = 0; a moment is checked when a view is made at it).
    /// Solves that failed at ctx (`failed`) are skipped: what rests on them
    /// is refused already (ModelView::refusal), and over a motion so is
    /// what rests on a solve depending on t (solved again at each step).
    pub fn check_conditions(&self, ctx: &Ctx, motion: Option<(f64, f64)>, failed: &[usize], at_load: bool) -> R<()> {
        use ink_kernel::ranged::{Cmp, Comparison, ProofOptions, Status, value_span};
        for c in &self.conditions {
            let solves = self.solves_read(&[c.len2]);
            let solved_on_t = solves.iter().any(|&i| self.solves[i].on_t);
            let reads_t = solved_on_t || self.g.depends_on(c.len2, |v| matches!(v, Var::T));
            if solves.iter().any(|i| failed.contains(i)) || (at_load && reads_t) || (motion.is_some() && solved_on_t) {
                continue;
            }
            let ranged = self.g.depends_on(c.len2, |v| matches!(v, Var::Param(k) if ctx.p_spans.get(k).is_some()));
            let tape = Tape::new(&self.g, c.len2);
            let err = |message: String| InkError { message, ..c.at.clone() };
            let named = |t: f64, values: &[(usize, f64)]| -> String {
                let mut at: Vec<String> = values.iter().map(|&(k, v)| format!("{} = {v}", self.params[k].name)).collect();
                if reads_t {
                    at.insert(0, format!("t = {t}"));
                }
                if at.is_empty() { String::new() } else { format!(" at {}", at.join(", ")) }
            };
            let over = motion.filter(|_| reads_t);
            if let Some(span) = over {
                // t over the span, an axis of the one search as a ranged
                // param is (with the params ranged, where they are).
                let opts = ProofOptions { res: None, finest: 0.0, depth: ink_kernel::scope::PARAM_SPLITS, tol: 0.0 };
                let over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..ctx.clone() };
                let v = value_span(&tape, &over, [0.0; 3], Comparison::exact(Cmp::Gt, 0.0), &opts);
                match (v.status, v.witness) {
                    (Status::Proven, _) if v.range.1.is_finite() => {}
                    (Status::Disproven, Some((at, len2))) => {
                        let values: Vec<(usize, f64)> = at.params.iter().map(|&(k, x)| (k as usize, x)).collect();
                        return Err(err(format!(
                            "{} at {} has zero length{} (squared length {len2}), within t in [{}, {}]: it has no direction there",
                            c.what,
                            c.site,
                            named(at.t, &values),
                            span.0,
                            span.1
                        )));
                    }
                    _ => {
                        return Err(err(format!(
                            "{} at {} is not proven to have a length at every moment of t in [{}, {}]: a direction must have one at every moment",
                            c.what, c.site, span.0, span.1
                        )));
                    }
                }
            } else if ranged {
                let opts = ProofOptions { res: None, finest: 0.0, depth: ink_kernel::scope::PARAM_SPLITS, tol: 0.0 };
                let v = value_span(&tape, ctx, [0.0; 3], Comparison::exact(Cmp::Gt, 0.0), &opts);
                let scope: Vec<String> = ctx
                    .p_spans
                    .iter()
                    .filter(|s| self.params_read(&[c.len2]).contains(&(s.k as usize)))
                    .map(|s| format!("{} in [{}, {}]", self.params[s.k as usize].name, s.span.0, s.span.1))
                    .collect();
                let scope = scope.join(", ");
                match v.status {
                    Status::Proven if v.range.1.is_finite() => {}
                    Status::Disproven => {
                        let (at, len2) = v.witness.expect("a disproof has a witness");
                        let values: Vec<(usize, f64)> = at.params.iter().map(|&(k, x)| (k as usize, x)).collect();
                        return Err(err(format!(
                            "{} at {} has zero length{} (squared length {len2}), within {scope}: it has no direction there",
                            c.what,
                            c.site,
                            named(at.t, &values)
                        )));
                    }
                    // Where the bounds reach 0 but no value is found to
                    // (a zero within rounding), not proven, and refused.
                    _ => {
                        return Err(err(format!(
                            "{} at {} is not proven to have a length over {scope} (its squared length is bounded within [{}, {}] there): a direction must have one at every value",
                            c.what, c.site, v.range.0, v.range.1
                        )));
                    }
                }
            } else {
                let len2 = eval_point(&tape, ctx, [0.0; 3]);
                if !(len2 > 0.0 && len2.is_finite()) {
                    let values: Vec<(usize, f64)> = self.params_read(&[c.len2]).into_iter().map(|k| (k, ctx.params[k])).collect();
                    let why = if len2 > 0.0 { "is too long to square in doubles" } else { "has zero length" };
                    return Err(err(format!("{} at {} {why}{}: it has no direction", c.what, c.site, named(ctx.t, &values))));
                }
            }
        }
        Ok(())
    }

    fn rotation(&mut self, axis: [NodeId; 3], deg: NodeId, at: &Sexp) -> R<Xf> {
        let hint = "(rot ax ay az deg) turns about the axis (ax, ay, az), which needs a nonzero length. If it reads a param or t, keep it nonzero over their ranges";
        let u = self.unit(&axis, at, &format!("the axis of {at}"), hint)?;
        let g = &mut self.g;
        let [s, c] = turn_sin_cos(g, deg);
        let one = g.constant(1.0);
        let omc = g.sub(one, c);
        // Rodrigues: R = c I + s [u]x + (1 - c) u u^T. A diagonal entry is
        // written u_i^2 + c (1 - u_i^2), c read once: along a constant axis
        // it is exactly 1, where (1 - c) + c would read the angle (and t,
        // over a motion) in an entry that never changes.
        let mut r = [[0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let uu = g.mul(u[i], u[j]);
                let mut e;
                if i == j {
                    let rest = g.sub(one, uu);
                    let cr = g.mul(c, rest);
                    e = g.add(uu, cr);
                } else {
                    e = g.mul(omc, uu);
                    // [u]x[i][j] = -eps_ijk u_k
                    let kk = 3 - i - j;
                    let sign = if (i + 1) % 3 == j { -1.0 } else { 1.0 };
                    let sgn = g.constant(sign);
                    let su = g.mul(s, u[kk]);
                    let term = g.mul(sgn, su);
                    e = g.add(e, term);
                }
                r[i][j] = e;
            }
        }
        let z = g.constant(0.0);
        Ok(Xf { r, t: [z, z, z] })
    }

    /// `(at STEP ...)`: steps apply in order; a step is (move dx dy dz),
    /// (rot-x deg), (rot-y deg), (rot-z deg), (rot ax ay az deg) or any
    /// frame. A placement reads no x, y or z (frame_value).
    fn eval_at(&mut self, steps: &[Sexp], env: &mut Env) -> R<Xf> {
        let usage = "usage: (at step ...) with steps (move dx dy dz) (rot-x deg) (rot-y deg) (rot-z deg) (rot ax ay az deg) or a frame, applied in order";
        let mut xf = self.identity();
        for step in steps {
            let head = step.list().and_then(|l| l.first()).and_then(|h| h.sym()).unwrap_or("");
            let s = if matches!(head, "move" | "rot-x" | "rot-y" | "rot-z" | "rot") {
                let sp = step.list().unwrap();
                let ps = match head {
                    "move" => "dx dy dz",
                    "rot" => "ax ay az deg",
                    _ => "deg",
                };
                if sp.len() - 1 != ps.split(' ').count() {
                    return Err(self.err_hint(step, format!("the step '{head}' takes {} number(s) ({ps}), got {}", ps.split(' ').count(), sp.len() - 1), Some(usage.into())));
                }
                let mut nums = Vec::new();
                for a in &sp[1..] {
                    nums.push(self.field(a, env)?);
                }
                let (z, o) = (self.g.constant(0.0), self.g.constant(1.0));
                match (head, nums.as_slice()) {
                    ("move", [dx, dy, dz]) => {
                        let mut m = self.identity();
                        m.t = [*dx, *dy, *dz];
                        m
                    }
                    ("rot-x", [d]) => self.rotation([o, z, z], *d, step)?,
                    ("rot-y", [d]) => self.rotation([z, o, z], *d, step)?,
                    ("rot-z", [d]) => self.rotation([z, z, o], *d, step)?,
                    ("rot", [ax, ay, az, d]) => self.rotation([*ax, *ay, *az], *d, step)?,
                    _ => unreachable!("counted above"),
                }
            } else {
                match self.eval(step, env).map_err(|err| self.misplaced(err, step, usage))? {
                    Value::Frame(f) => f,
                    other => return Err(self.err_hint(step, format!("bad placement step: a {}, not a step or a frame", other.kind()), Some(usage.into()))),
                }
            };
            xf = self.compose(&xf, &s);
        }
        Ok(xf)
    }

    /// `(axes (origin ox oy oz) (x xx xy xz) (z zx zy zz))`: a frame as ink
    /// prints one, its clauses in any order. z keeps its direction (a
    /// face's normal); x is made square to it, its part along z dropped
    /// (refused where nothing is left: x along z); y is z cross x. A
    /// placement reads no x, y or z (frame_value).
    fn eval_axes(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Xf> {
        let usage = "usage: (axes (origin ox oy oz) (x xx xy xz) (z zx zy zz)): the frame at origin with those directions, as ink check prints one \
                     (each clause 3 numbers or a vec of 3); z keeps its direction, x is made square to it, y is z cross x";
        let mut got: [Option<[NodeId; 3]>; 3] = [None; 3];
        for a in args {
            let (k, nums) = match a.list() {
                Some([k, nums @ ..]) => (k.sym().unwrap_or(""), nums),
                _ => ("", &[][..]),
            };
            let Some(slot) = ["origin", "x", "z"].iter().position(|&n| n == k) else {
                let why = if k == "y" { "(axes ...) takes no y: it is z cross x".to_string() } else { format!("bad clause {a} of (axes ...)") };
                return Err(self.err_hint(a, why, Some(usage.into())));
            };
            let v = match nums {
                [v] => self.fields(v, env)?,
                _ => nums.iter().map(|n| self.field(n, env)).collect::<R<Vec<_>>>()?,
            };
            if got[slot].is_some() || v.len() != 3 {
                let why = if v.len() != 3 { format!("({k} ...) takes 3 numbers or a vec of 3") } else { format!("({k} ...) is given twice") };
                return Err(self.err_hint(a, why, Some(usage.into())));
            }
            got[slot] = Some([v[0], v[1], v[2]]);
        }
        let [Some(o), Some(x), Some(z)] = got else { return Err(self.err_hint(e, "(axes ...) needs origin, x and z", Some(usage.into()))) };
        let z = self.unit(&z, e, &format!("the z of {e}"), "z is the frame's normal, which needs a nonzero length. If it reads a param or t, keep it nonzero over their ranges")?;
        // x less its part along z.
        let g = &mut self.g;
        let along = (0..3).fold(g.constant(0.0), |acc, i| {
            let p = g.mul(x[i], z[i]);
            g.add(acc, p)
        });
        let x: Vec<NodeId> = (0..3)
            .map(|i| {
                let p = g.mul(along, z[i]);
                g.sub(x[i], p)
            })
            .collect();
        let x = self.unit(&x, e, &format!("the x of {e}, square to its z,"), "x must not be parallel to z: x is made square to z, and of an x along z nothing is left")?;
        let g = &mut self.g;
        let y: Vec<NodeId> = (0..3)
            .map(|i| {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                let (a, b) = (g.mul(z[j], x[k]), g.mul(z[k], x[j]));
                g.sub(a, b)
            })
            .collect();
        Ok(Xf { r: std::array::from_fn(|i| [x[i], y[i], z[i]]), t: o })
    }

    /// Express a local-space value in world space: substitute p_local = R^T (p - t).
    pub fn to_world(&mut self, v: &Value, xf: &Xf) -> Value {
        let map = self.inverse_map(xf);
        let w = self.placed_value(v, xf, &map);
        if let (Value::Field(local, _), Value::Field(world, _)) = (v, &w) {
            self.record_placed(*world, *local, xf);
        }
        w
    }

    /// Record that `world` is `local` placed at `at` (`placed`). Equal
    /// fields are one node, so an entry is one true decomposition of the
    /// field, not the item's: a reader takes the shape and the placement
    /// from the same entry. One whose placement changes nothing (world is
    /// local) says nothing and would hide one that does.
    fn record_placed(&mut self, world: NodeId, local: NodeId, at: &Xf) {
        if world != local {
            self.placed.insert(world, Placed { shape: local, at: *at });
        }
    }

    /// `v` put at `xf`: a field warped by its inverse (`map`), a point
    /// moved by it (R p + T), a vec each of its items.
    fn placed_value(&mut self, v: &Value, xf: &Xf, map: &[(Var, NodeId)]) -> Value {
        match v {
            Value::Point(p) => Value::Point(std::array::from_fn(|i| {
                let mut acc = xf.t[i];
                for j in 0..3 {
                    let m = self.g.mul(xf.r[i][j], p[j]);
                    acc = self.g.add(acc, m);
                }
                acc
            })),
            Value::Vec(items) => Value::Vec(items.iter().map(|i| self.placed_value(i, xf, map)).collect()),
            other => self.warp_value(other, map),
        }
    }

    pub(crate) fn inverse_map(&mut self, xf: &Xf) -> Vec<(Var, NodeId)> {
        let p = [self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z)];
        let d: Vec<NodeId> = (0..3).map(|j| self.g.sub(p[j], xf.t[j])).collect();
        let mut out = Vec::new();
        for (i, var) in [Var::X, Var::Y, Var::Z].into_iter().enumerate() {
            let mut acc = self.g.constant(0.0);
            for j in 0..3 {
                let m = self.g.mul(xf.r[j][i], d[j]);
                acc = self.g.add(acc, m);
            }
            out.push((var, acc));
        }
        out
    }

    /// The subexpressions of `f` named `name` (by a def or a tag), among
    /// the pieces its value is built of (Graph::pieces).
    fn named_in(&self, f: NodeId, name: &str) -> Vec<NodeId> {
        self.g
            .pieces(&[f])
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Tag(t, a) if self.g.tag_names[t as usize] == name => Some(a),
                _ => None,
            })
            .collect()
    }

    /// Every construction name inside `f`.
    fn named_in_all(&self, f: NodeId) -> Vec<&str> {
        let mut names: Vec<&str> = self
            .g
            .pieces(&[f])
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Tag(t, _) => Some(self.g.tag_names[t as usize].as_str()),
                _ => None,
            })
            .collect();
        names.sort_unstable();
        names.dedup();
        names
    }

    fn warp_value(&mut self, v: &Value, map: &[(Var, NodeId)]) -> Value {
        match v {
            Value::Field(f, said) => Value::Field(self.g.subst(*f, map), said.clone()),
            Value::Vec(items) => Value::Vec(items.iter().map(|i| self.warp_value(i, map)).collect()),
            other => other.clone(),
        }
    }

    // ---------------------------------------------------------------- targets

    /// Resolve a body, assembly, def, or instance path (`asm/inst/...`) to
    /// its items: each body's own field and its placement (Item).
    pub fn items(&mut self, target: &str) -> Result<Vec<Item>, String> {
        let parts: Vec<&str> = target.split('/').collect();
        if parts.len() == 1 {
            if let Some(b) = self.body(target) {
                let (name, shape, attrs) = (b.name.clone(), b.shape, b.attrs.clone());
                let id = self.identity();
                return Ok(vec![Item { label: name.clone(), shape, attrs, instance: Some((name, id)), links: Vec::new() }]);
            }
            if self.assembly(target).is_some() {
                let id = self.identity();
                return self.flatten(target, &id, target, &[]);
            }
            if let Some(v) = self.def(target) {
                return match v {
                    Value::Field(f, _) => Ok(vec![Item { label: target.into(), shape: *f, attrs: Attrs::default(), instance: None, links: Vec::new() }]),
                    other => Err(format!("def '{target}' is a {}, not a field", other.kind())),
                };
            }
            let names = self.bodies.iter().map(|b| b.name.as_str()).chain(self.assemblies.iter().map(|a| a.name.as_str()));
            let hint = closest(target, names).unwrap_or_default();
            return Err(format!("no body, assembly, or field def named '{target}'. {hint}").trim().to_string());
        }
        let mut asm_name = parts[0].to_string();
        if self.assembly(&asm_name).is_none() {
            return Err(format!("'{asm_name}' is not an assembly"));
        }
        let mut xf = self.identity();
        // The instances from the world down to the target: an item's links
        // start with them, as its world placement does (a placement above
        // the target moves and scales what is under it too).
        let mut above: Vec<Link> = Vec::new();
        for (k, seg) in parts[1..].iter().enumerate() {
            let asm = self.assembly(&asm_name).unwrap();
            let inst = asm
                .instances
                .iter()
                .find(|i| i.name == *seg)
                .ok_or_else(|| format!("assembly '{asm_name}' has no instance '{seg}'. {}", choices(seg, "Its instances", asm.instances.iter().map(|i| i.name.as_str()))))?;
            let (child, at, joint, on) = (inst.child.clone(), inst.at, inst.joint.clone(), inst.on.clone());
            xf = self.compose(&at, &xf);
            above.push(Link { path: parts[..k + 2].join("/"), at, joint, on });
            let last = k == parts.len() - 2;
            if let Some(b) = self.body(&child) {
                if !last {
                    return Err(format!("'{seg}' is a body, with no instances"));
                }
                return Ok(vec![Item { label: target.to_string(), shape: b.shape, attrs: b.attrs.clone(), instance: Some((child, xf)), links: above }]);
            }
            if last {
                return self.flatten(&child, &xf, target, &above);
            }
            asm_name = child;
        }
        unreachable!()
    }

    /// The items of assembly `asm` placed at `parent`, labelled under
    /// `prefix`, each with the instances above it (`links`, then its own).
    fn flatten(&mut self, asm: &str, parent: &Xf, prefix: &str, links: &[Link]) -> Result<Vec<Item>, String> {
        let a = self.assembly(asm).unwrap();
        let insts: Vec<(String, String, Xf, Option<Joint>, Vec<String>)> =
            a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at, i.joint.clone(), i.on.clone())).collect();
        let mut out = Vec::new();
        for (name, child, at, joint, on) in insts {
            let xf = self.compose(&at, parent);
            let label = format!("{prefix}/{name}");
            let mut path = links.to_vec();
            path.push(Link { path: label.clone(), at, joint, on });
            if let Some(b) = self.body(&child) {
                out.push(Item { label, shape: b.shape, attrs: b.attrs.clone(), instance: Some((child, xf)), links: path });
            } else {
                out.extend(self.flatten(&child, &xf, &label, &path)?);
            }
        }
        Ok(out)
    }

    /// An item's field where it is placed (`warped`), recorded as its
    /// body's field and placement (`Doc::placed`).
    /// What world field `w` was placed from, where it was (`Doc::placed`,
    /// or a field warped by a placement: a zone, an attribute): its field
    /// in its own frame and the placement.
    pub fn placement_of(&self, w: NodeId) -> Option<Placed> {
        self.placed
            .get(&w)
            .copied()
            .or_else(|| self.placed_fields.iter().filter(|(k, x)| **x == w && k.0 != w).map(|(k, _)| Placed { shape: k.0, at: k.1 }).min_by_key(|p| p.shape))
    }

    /// Every world field placed (`placement_of`), in node order.
    pub fn placements(&self) -> Vec<(NodeId, Placed)> {
        let mut out: Vec<(NodeId, Placed)> = self.placed.iter().map(|(w, p)| (*w, *p)).collect();
        out.extend(self.placed_fields.iter().filter(|((f, _), w)| *w != f && !self.placed.contains_key(w)).map(|(k, w)| (*w, Placed { shape: k.0, at: k.1 })));
        out.sort_unstable_by_key(|x| (x.0, x.1.shape));
        out.dedup_by_key(|x| x.0);
        out
    }

    pub fn placed_field(&mut self, it: &Item) -> NodeId {
        let Some((_, xf)) = &it.instance else { return it.shape };
        let f = self.warped(it.shape, xf);
        self.record_placed(f, it.shape, xf);
        f
    }

    /// A field of an item's frame (an attribute) where the item is placed:
    /// one reading none of x, y and z is itself.
    pub fn place(&mut self, f: NodeId, it: &Item) -> NodeId {
        match &it.instance {
            Some((_, xf)) if self.xyz_reads(f).contains(&true) => self.warped(f, xf),
            _ => f,
        }
    }

    /// `f` warped by placement `xf`, once per field and placement (by the
    /// identity: itself).
    pub fn warped(&mut self, f: NodeId, xf: &Xf) -> NodeId {
        if let Some(&w) = self.placed_fields.get(&(f, *xf)) {
            return w;
        }
        let map = self.inverse_map(xf);
        let identity = map.iter().all(|&(v, u)| matches!(self.g.nodes[u as usize], crate::graph::Node::Var(w) if w == v));
        let w = if identity { f } else { self.g.subst(f, &map) };
        self.placed_fields.insert((f, *xf), w);
        w
    }

    /// A value as an error message shows it: text as written, a number (a
    /// param's at its value) in its shortest form (24, 2.5), a vec's items, a
    /// record as JSON, anything else by kind.
    fn value_text(&self, v: &Value) -> String {
        match v {
            Value::Str(s) => s.clone(),
            Value::Field(f, _) => self.load_value(*f).map_or_else(|| "(a field of x, y, z or t)".into(), |x| format!("{x}")),
            Value::Vec(items) => format!("[{}]", items.iter().map(|x| self.value_text(x)).collect::<Vec<_>>().join(", ")),
            Value::Data(d) => d.value().to_string(),
            other => format!("(a {})", other.kind()),
        }
    }

    /// An item's attributes where it is placed: their fields `place`d.
    pub fn placed_attrs(&mut self, it: &Item) -> Attrs {
        fn placed(doc: &mut Doc, v: &Value, it: &Item) -> Value {
            match v {
                Value::Field(f, said) => Value::Field(doc.place(*f, it), said.clone()),
                Value::Vec(items) => Value::Vec(items.iter().map(|x| placed(doc, x, it)).collect()),
                other => other.clone(),
            }
        }
        Attrs::new(it.attrs.iter().map(|(k, v)| (k.clone(), placed(self, v, it))).collect())
    }

    /// Fields reading what the items' placed fields read (params, samples,
    /// t), none of them placed: each item's own field and, for each of x, y
    /// and z it reads, what its placement puts there.
    pub fn placed_reads(&mut self, items: &[Item]) -> Vec<NodeId> {
        let mut out = Vec::new();
        for it in items {
            out.push(it.shape);
            let Some((_, xf)) = &it.instance else { continue };
            let reads = self.xyz_reads(it.shape);
            if reads.contains(&true) {
                let map = self.inverse_map(xf);
                out.extend(map.into_iter().zip(reads).filter(|(_, r)| *r).map(|((_, u), _)| u));
            }
        }
        out
    }

    /// Which of x, y and z field `f` reads (once per field).
    fn xyz_reads(&mut self, f: NodeId) -> [bool; 3] {
        let g = &self.g;
        *self.xyz_read.entry(f).or_insert_with(|| {
            let r = g.reachable(&[f]);
            [Var::X, Var::Y, Var::Z].map(|v| r.iter().any(|&n| matches!(g.nodes[n as usize], crate::graph::Node::Var(w) if w == v)))
        })
    }

    /// The default target: `root` if it exists, else every body.
    /// The assemblies written in the file loaded that none of its other
    /// assemblies places, in the order written: the models it puts
    /// together (what `ink render` draws without a target or a root).
    pub fn top_assemblies(&self) -> Vec<&str> {
        let own: Vec<&Assembly> = self.assemblies.iter().filter(|a| a.of.is_none() && self.is_own_file(&a.file)).collect();
        own.iter().filter(|a| !own.iter().any(|b| b.instances.iter().any(|i| i.child == a.name))).map(|a| a.name.as_str()).collect()
    }

    pub fn default_items(&mut self) -> Result<Vec<Item>, String> {
        if self.body("root").is_some() || self.assembly("root").is_some() {
            return self.items("root");
        }
        if self.bodies.is_empty() {
            return Err("document has no bodies; add (body name (shape expr))".into());
        }
        let names: Vec<String> = self.bodies.iter().filter(|b| b.of.is_none()).map(|b| b.name.clone()).collect();
        let mut out = Vec::new();
        for n in names {
            out.extend(self.items(&n)?);
        }
        Ok(out)
    }

    // ---------------------------------------------------------------- helpers

    /// The window: where fields that are not bounded are examined (and
    /// viewed). Bounded shapes are analysed whole wherever they lie.
    pub fn domain(&self) -> (V3, V3) {
        DEFAULT_DOMAIN
    }

    /// The finest cell a proof refines to unless given a resolution: half
    /// the accuracy (ACCURACY of the window), the cell that proves touching
    /// surfaces to it: the same for every model, so a part's claims look as
    /// fine in an assembly as the assembly's own.
    pub fn finest(&self) -> f64 {
        let (lo, hi) = self.domain();
        ACCURACY / 2.0 * crate::query::max_extent(lo, hi)
    }

    /// Where analysis of field `f` looks (over `ctx.t_span`, if set, at
    /// every moment of it):
    /// a box proven to hold all of its material (`query::extent`, a search
    /// over all of space), or, if the field is not proven bounded, the
    /// window (`domain`), and then results that need all of the shape say so.
    /// A shape inside the window gets the window itself.
    pub fn reach(&self, f: NodeId, ctx: &Ctx) -> Reach {
        reach_in(&self.g, self.domain(), f, ctx)
    }

    /// Field that is negative strictly inside an axis-aligned box.
    pub fn box_field(&mut self, lo: V3, hi: V3) -> NodeId {
        let lo = lo.map(|c| self.g.constant(c));
        let hi = hi.map(|c| self.g.constant(c));
        self.box_nodes(lo, hi)
    }

    /// `box_field` with corners that are fields (of constants and params).
    fn box_nodes(&mut self, lo: [NodeId; 3], hi: [NodeId; 3]) -> NodeId {
        let vars = [Var::X, Var::Y, Var::Z];
        let mut terms = Vec::new();
        for i in 0..3 {
            let v = self.g.var(vars[i]);
            terms.push(self.g.sub(lo[i], v));
            terms.push(self.g.sub(v, hi[i]));
        }
        self.reduce(BinOp::Max, &terms)
    }

    /// Intersect a field with a box (keeps every query bounded and closed).
    pub fn clip(&mut self, f: NodeId, lo: V3, hi: V3) -> NodeId {
        let b = self.box_field(lo, hi);
        self.g.max(f, b)
    }

    /// `f` clipped to where analysis looks (see `reach`), and that region.
    /// Clipping a bounded field to its reach changes no sign: it has no
    /// material outside.
    pub fn clip_reach(&mut self, f: NodeId, ctx: &Ctx) -> (NodeId, Reach) {
        let r = self.reach(f, ctx);
        (self.clip(f, r.lo, r.hi), r)
    }

    /// A context at moment t with every solve solved there (`ctx_at`; Err:
    /// the first that fails, or a direction with no length there:
    /// Doc::check_conditions). Rates no shape reads are not needed there:
    /// what reads them refuses them itself (Doc::needed).
    pub fn ctx(&self, overrides: &[(String, Given)], t: f64) -> Result<Ctx, String> {
        let (ctx, failed) = self.ctx_at(overrides, t)?;
        match failed.iter().find(|(i, _)| self.needed(*i)) {
            Some((_, e)) => Err(e.message.clone()),
            None => {
                self.check_conditions(&ctx, None, &[], false).map_err(|e| e.message)?;
                Ok(ctx)
            }
        }
    }

    /// Evaluate a space-independent node to a number.
    pub fn scalar(&self, n: NodeId, ctx: &Ctx) -> f64 {
        eval_point(&Tape::new(&self.g, n), ctx, [0.0; 3])
    }

    /// The fields combined by `op` (min or max), as a balanced tree. Of
    /// none: the union is the empty shape (as `integrands` reads one), the
    /// intersection all of space.
    pub fn reduce(&mut self, op: BinOp, ids: &[NodeId]) -> NodeId {
        if ids.is_empty() {
            assert!(matches!(op, BinOp::Min | BinOp::Max), "{op:?} of no fields");
            return self.g.constant(if op == BinOp::Min { 1.0 } else { -1.0 });
        }
        if ids.len() == 1 {
            return ids[0];
        }
        let m = ids.len() / 2;
        let a = self.reduce(op, &ids[..m]);
        let b = self.reduce(op, &ids[m..]);
        self.g.bin(op, a, b)
    }

    /// A fn's or lambda's argument names. The last may be followed by
    /// `...`: it then takes the remaining arguments (none or more), as a vec.
    fn arg_names(&self, ps: &Sexp, what: &str) -> R<(Vec<String>, bool)> {
        let pl = ps.list().ok_or_else(|| self.err(ps, format!("{what} arguments must be a list of names")))?;
        let rest = pl.last().and_then(|p| p.sym()) == Some("...");
        let names = if rest { &pl[..pl.len() - 1] } else { pl };
        if rest && names.is_empty() {
            return Err(self.err(ps, "'...' follows the argument that takes the rest: (a rest ...)"));
        }
        let mut params: Vec<String> = Vec::new();
        for p in names {
            let pn = self.name_of(p)?;
            if pn == "..." {
                return Err(self.err(p, "'...' goes only at the end, after the argument that takes the rest: (a rest ...)"));
            }
            if params.contains(&pn) {
                return Err(self.err(p, format!("duplicate argument '{pn}'")));
            }
            params.push(pn);
        }
        Ok((params, rest))
    }

    /// A name: a bare word. A word quoted or in parentheses says to write it bare.
    fn name_of(&self, e: &Sexp) -> R<String> {
        e.sym().map(str::to_string).ok_or_else(|| {
            let word = match &e.k {
                SK::Str(s) => parse(s).ok().and_then(|f| f.first().and_then(|w| w.sym()).filter(|_| f.len() == 1).map(str::to_string)),
                SK::List(l) => l.first().and_then(|w| w.sym()).filter(|_| l.len() == 1).map(str::to_string),
                _ => None,
            };
            self.err_hint(e, format!("expected a name, got {e}"), word.map(|w| format!("a name is a bare word: {w}")))
        })
    }

    fn const_v3(&mut self, e: &Sexp) -> R<V3> {
        let parts = e.list().filter(|p| p.len() == 3).ok_or_else(|| self.err(e, "expected a point (x y z)"))?;
        let mut out = [0.0; 3];
        for i in 0..3 {
            let f = self.field(&parts[i], &mut Vec::new())?;
            out[i] = self.g.as_const(f).ok_or_else(|| self.err(&parts[i], "expected a constant number"))?;
        }
        Ok(out)
    }

    /// A point (x y z) of numbers: constants and params, kept as fields
    /// (not evaluated, so its params stay symbolic).
    /// A claim's operand: a shape (Ok), or a point (Err: `probe`).
    fn shape_or_point(&mut self, e: &Sexp, env: &mut Env) -> R<Result<NodeId, [NodeId; 3]>> {
        match self.eval(e, env)? {
            Value::Field(f, _) => Ok(Ok(f)),
            Value::Point(p) => Ok(Err(p)),
            other => Err(self.err(e, format!("expected a shape (a field) or a point (probe), got a {}: {e}", other.kind()))),
        }
    }

    fn corner(&mut self, e: &Sexp) -> R<[NodeId; 3]> {
        let parts = e.list().filter(|p| p.len() == 3).ok_or_else(|| self.err(e, "expected a point (x y z)"))?;
        let mut out = [0; 3];
        for i in 0..3 {
            let f = self.field(&parts[i], &mut Vec::new())?;
            if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
                return Err(self.err(&parts[i], "a corner is a number (constants and params), not a field of x y z or t"));
            }
            out[i] = f;
        }
        Ok(out)
    }

    fn field(&mut self, e: &Sexp, env: &mut Env) -> R<NodeId> {
        let v = self.eval(e, env)?;
        self.as_field(v, e)
    }

    /// Value `v` of `e` as a field.
    fn as_field(&self, v: Value, e: &Sexp) -> R<NodeId> {
        match v {
            Value::Field(f, _) => Ok(f),
            other => Err(self.err(e, format!("expected a field (number expression), got a {}: {e}", other.kind()))),
        }
    }

    // ---------------------------------------------------------------- expressions

    fn lookup(&mut self, s: &str, at: &Sexp, env: &Env) -> R<Value> {
        if let Some((_, v)) = env.iter().rev().find(|(k, _)| k == s) {
            return Ok(v.clone());
        }
        // Lexical: what is visible where it is read (Doc::scoped).
        if let Some(e) = self.later(s, at) {
            return Err(e);
        }
        let key = self.scoped(s);
        if let Some(i) = self.param_index(&key) {
            return Ok(Value::Field(self.g.var(Var::Param(i as u32)), self.param_said.get(&i).cloned()));
        }
        if let Some(v) = self.def(&key) {
            return Ok(v.clone());
        }
        let v = match s {
            "x" => self.g.var(Var::X),
            "y" => self.g.var(Var::Y),
            "z" => self.g.var(Var::Z),
            "t" => self.g.var(Var::T),
            "pi" => self.g.constant(std::f64::consts::PI),
            "inf" => self.g.constant(f64::INFINITY),
            _ => {
                if let Some(f) = self.fns.get(&*key) {
                    return Ok(Value::Fn(Callable::User(s.to_string(), f.clone())));
                }
                if let Some(hint) = self.hidden_hint(s) {
                    return Err(self.err_hint(at, format!("unknown name '{s}'"), Some(hint)));
                }
                let o = self.resolve(s);
                if let Some(kind) = self.body(&o).map(|_| "body").or(self.assembly(&o).map(|_| "assembly")) {
                    return Err(self.err_hint(
                        at,
                        format!("'{s}' is {} {kind}, not a value", if kind == "body" { "a" } else { "an" }),
                        Some(format!("its field is (shape-of {s})")),
                    ));
                }
                if builtin(s).is_some() || COMPARISONS.contains(&s) || FORM_FUNCTIONS.contains(&s) {
                    return Ok(Value::Fn(Callable::Builtin(s.to_string())));
                }
                if SPECIAL.contains(&s) {
                    return Err(self.err_hint(
                        at,
                        format!("'{s}' is a form, not a function: it cannot be passed as a value"),
                        Some(format!("wrap it in a function of its arguments, e.g. (lambda (a b) ({s} a b))")),
                    ));
                }
                let mut names = self.value_names();
                names.extend(env.iter().map(|(k, _)| k.clone()));
                names.extend(CONSTANTS.iter().map(|c| c.to_string()));
                // Declared as something else (a body or an assembly further on,
                // an instance) or an import's, it is not to be declared again.
                let instance = self.building.as_ref().is_some_and(|b| b.asm.instances.iter().any(|i| i.name == s) || b.later.iter().any(|l| l == s));
                let undeclared = !s.contains('.') && !instance && !self.files[self.scope.file as usize].declared.contains(s);
                let hint = closest(s, names.iter().map(|n| n.as_str())).or_else(|| undeclared.then(|| format!("declare it: (param {s} VALUE) or (def {s} EXPR)")));
                return Err(self.err_hint(at, format!("unknown name '{s}'"), hint));
            }
        };
        Ok(Value::Field(v, None))
    }

    fn callable(&mut self, e: &Sexp, env: &mut Env) -> R<Callable> {
        match self.eval(e, env)? {
            Value::Fn(f) => Ok(f),
            other => Err(self.err(e, format!("expected a function, got a {}", other.kind()))),
        }
    }

    fn list(&mut self, e: &Sexp, env: &mut Env) -> R<Vec<Value>> {
        match self.eval(e, env)? {
            Value::Vec(v) => Ok(v),
            other => Err(self.err(e, format!("expected a list (vec), got a {}", other.kind()))),
        }
    }

    fn const_num(&mut self, e: &Sexp, env: &mut Env) -> R<f64> {
        let f = self.field(e, env)?;
        self.const_of(e, f)
    }

    /// Field `f` of `e` as a number known at load time.
    fn const_of(&mut self, e: &Sexp, f: NodeId) -> R<f64> {
        self.load_number(e, f)?.ok_or_else(|| self.err_hint(e, "expected a number known at load time", Some("constants and params work here; x, y, z and t do not".into())))
    }

    /// The value of data `d`: numbers as constants, strings, arrays as
    /// vecs, an object as data (the place, not a copy).
    fn json_value(&mut self, d: &Data) -> Result<Value, String> {
        use serde_json::Value as J;
        let at = |v: &J, place: &dyn Fn() -> Data, me: &mut Self| -> Result<Value, String> {
            Ok(match v {
                J::Number(n) => Value::Field(me.g.constant(n.as_f64().ok_or("number out of range")?), None),
                J::Bool(b) => Value::Field(me.g.constant(if *b { 1.0 } else { 0.0 }), None),
                J::String(s) => Value::Str(s.clone()),
                J::Array(_) | J::Object(_) => return me.json_value(&place()),
                J::Null => return Err("null is not a value".into()),
            })
        };
        match d.value() {
            J::Object(_) => Ok(Value::Data(d.clone())),
            J::Array(a) => Ok(Value::Vec(a.iter().enumerate().map(|(i, v)| at(v, &|| d.nth(i), self)).collect::<Result<_, _>>()?)),
            v => at(v, &|| d.clone(), self),
        }
    }

    /// Whether a promise made here (its k and field, `roots`) reads an input
    /// of the call it is made in that varies over space. A promise in a
    /// function is about its body as a function of x, y, z with its inputs
    /// fixed (a sphere's radius): true of every fixed radius, it says
    /// nothing of a radius that reads x, whose slope adds to the field's
    /// (the derivative of r(x) - |p| has r'(x) in it). Inputs are the
    /// fields among what the call received (arguments and captured values,
    /// vecs included; t and params may vary, a promise holds for each of
    /// their values). Functions received are code, not inputs: a promise
    /// over what one returns is its author's.
    fn reads_moving_input(&self, roots: &[NodeId], env: &Env) -> bool {
        let n = self.inputs.last().map_or(0, |&n| n.min(env.len()));
        let mut inputs = Vec::new();
        env[..n].iter().for_each(|(_, v)| v.fields_into(&mut inputs, false));
        inputs.retain(|&i| self.g.as_const(i).is_none());
        if inputs.is_empty() {
            return false;
        }
        let read = self.g.reachable(roots);
        inputs.iter().any(|&i| read.binary_search(&i).is_ok() && self.g.depends_on(i, |v| matches!(v, Var::X | Var::Y | Var::Z)))
    }

    /// Call a function value with evaluated arguments.
    fn call(&mut self, f: &Callable, args: Vec<Value>, at: &Sexp) -> R<Value> {
        let (label, params, rest, body, file, mut env) = match f {
            Callable::Builtin(name) if FORM_FUNCTIONS.contains(&name.as_str()) => return self.call_form(name, args, at),
            Callable::Builtin(name) => return self.apply_builtin(name, args, at),
            Callable::User(name, d) => (name.clone(), d.params.clone(), d.rest, d.body.clone(), (d.file.clone(), d.scope), Vec::new()),
            Callable::Closure(c) => ("lambda".to_string(), c.params.clone(), c.rest, c.body.clone(), (c.file.clone(), c.scope), c.env.clone()),
        };
        let mut args = args;
        if rest && args.len() + 1 >= params.len() {
            let more = args.split_off(params.len() - 1);
            args.push(Value::Vec(more));
        }
        if args.len() != params.len() {
            let (n, ps) = if rest { (format!("{} or more", params.len() - 1), format!("{} ...", params.join(" "))) } else { (params.len().to_string(), params.join(" ")) };
            return Err(self.err(at, format!("'{label}' takes {n} argument(s) ({ps}), got {}", args.len())));
        }
        env.extend(params.into_iter().zip(args));
        if self.depth >= MAX_CALL_DEPTH {
            return Err(self.err(at, format!("call depth exceeded {MAX_CALL_DEPTH} in '{label}' (recursion is not allowed)")));
        }
        self.depth += 1;
        self.inputs.push(env.len());
        self.calls.push((label, self.file.clone(), at.line));
        // The body reads the names of where it was written (Scope).
        let saved = (std::mem::replace(&mut self.file, file.0), std::mem::replace(&mut self.scope, file.1));
        let r = self.eval(&body, &mut env);
        (self.file, self.scope) = saved;
        let (label, ..) = self.calls.pop().expect("pushed above");
        self.inputs.pop();
        self.depth -= 1;
        // An error inside the body: record this call, so it can be traced
        // back to the author's line.
        r.map_err(|mut e| {
            e.called_from.push((label, self.file.to_string(), at.line, at.col));
            e
        })
    }

    /// Call a form that is a function (FORM_FUNCTIONS) with evaluated
    /// arguments: the form, written at the call site, reading each argument
    /// by a name bound to its value, as a call of `(lambda (a ...) (form a
    /// ...))` would.
    fn call_form(&mut self, name: &str, args: Vec<Value>, at: &Sexp) -> R<Value> {
        let sx = |k| Sexp { k, line: at.line, col: at.col };
        let mut items = vec![sx(SK::Sym(name.to_string()))];
        let mut env: Env = Vec::with_capacity(args.len());
        for (i, v) in args.into_iter().enumerate() {
            let a = format!("arg-{}", i + 1);
            items.push(sx(SK::Sym(a.clone())));
            env.push((a, v));
        }
        self.inputs.push(env.len());
        let r = self.eval(&sx(SK::List(items)), &mut env);
        self.inputs.pop();
        r
    }

    /// Apply a kernel op or comparison. Vec arguments broadcast element-wise
    /// (scalars repeat), so vector maths reads naturally: (- p a), (* v 2).
    fn apply_builtin(&mut self, name: &str, args: Vec<Value>, at: &Sexp) -> R<Value> {
        if let Some(n) = args.iter().find_map(|v| if let Value::Vec(v) = v { Some(v.len()) } else { None }) {
            // Every vec the same length (an empty one too), and as many
            // arguments as the op takes.
            let takes = match builtin(name) {
                Some(Builtin::Un(_)) => Some(1),
                Some(Builtin::Bin(_)) => Some(2),
                Some(Builtin::Select) => Some(3),
                Some(Builtin::Minus) if !(1..=2).contains(&args.len()) => Some(2),
                _ if COMPARISONS.contains(&name) => Some(2),
                _ => None,
            };
            if let Some(k) = takes.filter(|&k| k != args.len()) {
                return Err(self.err(at, format!("'{name}' takes {k} argument(s), got {}", args.len())));
            }
            if let Some(m) = args.iter().find_map(|v| if let Value::Vec(v) = v { Some(v.len()).filter(|&m| m != n) } else { None }) {
                return Err(self.err(at, format!("'{name}' on vecs of lengths {n} and {m}")));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                let row = args
                    .iter()
                    .map(|v| match v {
                        Value::Vec(v) if v.len() == n => Ok(v[i].clone()),
                        Value::Vec(v) => Err(self.err(at, format!("'{name}' on vecs of lengths {n} and {}", v.len()))),
                        other => Ok(other.clone()),
                    })
                    .collect::<R<Vec<_>>>()?;
                out.push(self.apply_builtin(name, row, at)?);
            }
            return Ok(Value::Vec(out));
        }
        // A choice of frames is one of their numbers: a part that changes
        // hands (picked up at a moment) is placed by (select c F G).
        if let (Some(Builtin::Select), [Value::Field(c, _), Value::Frame(a), Value::Frame(b)]) = (builtin(name), args.as_slice()) {
            let mut pick = |p: NodeId, q: NodeId| self.g.select(*c, p, q);
            let r = std::array::from_fn(|i| std::array::from_fn(|k| pick(a.r[i][k], b.r[i][k])));
            return Ok(Value::Frame(Xf { r, t: std::array::from_fn(|i| pick(a.t[i], b.t[i])) }));
        }
        let mut fs = Vec::with_capacity(args.len());
        for v in &args {
            match v {
                Value::Field(f, _) => fs.push(*f),
                other => {
                    let frames = if name == "select" { ", or a number and two frames" } else { "" };
                    let got = match other {
                        Value::Str(s) => format!("the string \"{s}\""),
                        o => format!("a {}", o.kind()),
                    };
                    return Err(self.err(at, format!("'{name}' expects numbers or vecs{frames}, got {got}")));
                }
            }
        }
        let arity = |n: usize, d: &Doc| {
            if fs.len() != n { Err(d.err(at, format!("'{name}' takes {n} argument(s), got {}", fs.len()))) } else { Ok(()) }
        };
        if COMPARISONS.contains(&name) {
            arity(2, self)?;
            let (one, zero) = (self.g.constant(1.0), self.g.constant(0.0));
            let (a, b) = (fs[0], fs[1]);
            let (ab, ba) = (self.g.sub(a, b), self.g.sub(b, a));
            let id = match name {
                "<" => self.g.select(ab, one, zero),
                ">=" => self.g.select(ab, zero, one),
                ">" => self.g.select(ba, one, zero),
                "<=" => self.g.select(ba, zero, one),
                _ => {
                    let ne = self.g.select(ba, zero, one);
                    self.g.select(ab, zero, ne)
                }
            };
            return Ok(Value::Field(id, None));
        }
        let id = match builtin(name).expect("builtin") {
            Builtin::Un(op) => {
                arity(1, self)?;
                self.g.un(op, fs[0])
            }
            Builtin::Bin(op) => {
                arity(2, self)?;
                self.g.bin(op, fs[0], fs[1])
            }
            Builtin::Variadic(op) => {
                if fs.is_empty() {
                    return Err(self.err(at, format!("'{name}' needs at least one argument")));
                }
                self.reduce(op, &fs)
            }
            Builtin::Minus => match fs.as_slice() {
                [] => return Err(self.err(at, "'-' needs at least one argument")),
                [a] => self.g.neg(*a),
                // (- a b c ...) = a - b - c - ...
                [a, rest @ ..] => rest.iter().fold(*a, |acc, b| self.g.sub(acc, *b)),
            },
            Builtin::Select => {
                arity(3, self)?;
                self.g.select(fs[0], fs[1], fs[2])
            }
        };
        Ok(Value::Field(id, self.sum_said(name, &args, at)?))
    }

    fn eval(&mut self, e: &Sexp, env: &mut Env) -> R<Value> {
        let items = match &e.k {
            SK::Num(v, unit) => {
                let f = if self.scope == Scope::STD {
                    self.g.own(*v)
                } else {
                    if self.wild_written.is_none() && crate::graph::wild(*v) {
                        self.wild_written = Some((*v, format!("{}:{}", self.file, e.line)));
                    }
                    self.g.constant(*v)
                };
                let said = unit.as_ref().map(|w| match w {
                    sexp::Written::Here(unit, written) => Rc::new(units::Said { unit, written: *written, file: std::sync::Arc::from(&*self.file), line: e.line }),
                    sexp::Written::From(s) => Rc::new((**s).clone()),
                });
                return Ok(Value::Field(f, said));
            }
            SK::Str(s) => return Ok(Value::Str(s.clone())),
            SK::Sym(s) => return self.lookup(s, e, env),
            SK::List(items) => items,
        };
        let Some(head) = items.first() else { return Err(self.err(e, "empty list () is not an expression")) };
        let Some(name) = head.sym() else {
            // ((lambda (x) ...) 3), or any expression that yields a function.
            let Value::Fn(f) = self.eval(head, env)? else {
                return Err(self.err(head, format!("expected a function at the head of {e}")));
            };
            let vals = items[1..].iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.call(&f, vals, e);
        };
        let args = &items[1..];
        match name {
            "let" => {
                let [binds, body] = args else { return Err(self.err(e, "usage: (let ((name expr) ...) body)")) };
                let bl = binds.list().ok_or_else(|| self.err(binds, "let bindings must be a list of (name expr)"))?;
                let mark = env.len();
                for b in bl {
                    let [n, v] = b.list().unwrap_or(&[]) else { return Err(self.err(b, "let binding must be (name expr)")) };
                    let n = self.name_of(n)?;
                    if is_reserved(&n) {
                        return Err(self.err_hint(b, format!("'{n}' is a built-in name"), Some(reserved_hint())));
                    }
                    let v = self.eval(v, env)?;
                    env.push((n, v));
                }
                let r = self.eval(body, env);
                env.truncate(mark);
                return r;
            }
            "vec" => {
                let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
                return Ok(Value::Vec(vals));
            }
            "unit" => {
                let [v] = args else { return Err(self.err(e, "usage: (unit vec): the unit vector along a vec of fields")) };
                let Value::Vec(items) = self.eval(v, env)? else { return Err(self.err(v, "unit expects a vec")) };
                let fs = items.iter().map(|i| if let Value::Field(f, _) = i { Some(*f) } else { None }).collect::<Option<Vec<_>>>();
                let fs = fs.filter(|f| !f.is_empty()).ok_or_else(|| self.err(v, "unit expects a vec of fields (numbers)"))?;
                let hint = "(unit v) is v / |v|, which needs v nonzero (a path's tangent is zero between two equal points, and a sweep's first normal where `up` lies along its first segment); where it reads a param or t, keep it nonzero over their ranges";
                let u = self.unit(&fs, e, &format!("the vector of {e}"), hint)?;
                return Ok(Value::Vec(u.into_iter().map(|f| Value::Field(f, None)).collect()));
            }
            "lambda" | "map" | "reduce" | "fold" | "flatten" | "apply" | "range" | "len" | "concat" | "char-codes" | "slice" | "if" | "read-json" | "get" | "diff" => {
                return self.eval_author(name, e, args, env);
            }
            "nth" => {
                let [v, i] = args else { return Err(self.err(e, "usage: (nth vec index)")) };
                let Value::Vec(items) = self.eval(v, env)? else { return Err(self.err(v, "nth expects a vec")) };
                let i = self.whole(i, env, "an nth index")?;
                // Negative indices count from the end: -1 is the last element.
                let k = if i < 0.0 { items.len() as f64 + i } else { i };
                return (k >= 0.0)
                    .then(|| items.get(k as usize).cloned())
                    .flatten()
                    .ok_or_else(|| self.err(e, format!("index {i} out of range for vec of length {}", items.len())));
            }
            "data" => {
                let v = match args {
                    [single] => single.to_json(),
                    _ => serde_json::Value::Array(args.iter().map(|a| a.to_json()).collect()),
                };
                return Ok(Value::Data(Data::new(v)));
            }
            "warp" => {
                if args.len() != 4 && args.len() != 5 {
                    return Err(self.err(e, "usage: (warp expr fx fy fz [ft])"));
                }
                let v = self.eval(&args[0], env)?;
                let mut map = Vec::new();
                for (var, a) in [Var::X, Var::Y, Var::Z, Var::T].into_iter().zip(&args[1..]) {
                    map.push((var, self.field(a, env)?));
                }
                // A warp says where a shape's points are read from, not
                // where a point goes: that is its inverse, which rounding
                // can put anywhere (a matrix singular in real numbers whose
                // rounded determinant is not 0).
                if has_point(&v) {
                    return Err(self.err_hint(
                        e,
                        "a warp moves a shape, not a point (probe): it says where the shape is read from, and a point would need its inverse",
                        Some("put a point at a frame, which moves it forward: (place (probe x y z) (at (move dx dy dz) (rot-z deg)))".into()),
                    ));
                }
                return Ok(self.warp_value(&v, &map));
            }
            "probe" => {
                let [px, py, pz] = args else {
                    return Err(self.err_hint(e, "wrong arguments", Some("usage: (probe x y z): a point, for claims that it is in a shape or not".into())));
                };
                let mut at = [0; 3];
                for (k, a) in [px, py, pz].into_iter().enumerate() {
                    let f = self.field(a, env)?;
                    if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                        return Err(self.err(a, "a point's coordinates are numbers (constants and params), not fields of x y z"));
                    }
                    at[k] = f;
                }
                return Ok(Value::Point(at));
            }
            "tag" => {
                let usage = "usage: (tag name field [opaque]) names a construction for provenance (faces); the name is a symbol, a string, \
                             or an expression giving a string; opaque: points under it report this name only, not the names inside";
                let (n, f, opaque) = match args {
                    [n, f] => (n, f, false),
                    [n, f, o] if o.sym() == Some("opaque") => (n, f, true),
                    _ => return Err(self.err_hint(e, "bad tag", Some(usage.into()))),
                };
                let name = match &n.k {
                    SK::Sym(s) | SK::Str(s) => s.clone(),
                    SK::List(_) => match self.eval(n, env)? {
                        Value::Str(s) => s,
                        other => return Err(self.err_hint(n, format!("a tag name must be a string, not a {}", other.kind()), Some(usage.into()))),
                    },
                    _ => return Err(self.err_hint(n, "tag name must be a name or string", Some(usage.into()))),
                };
                let f = self.field(f, env)?;
                return Ok(Value::Field(if opaque { self.g.tag_opaque(&name, f) } else { self.g.tag(&name, f) }, None));
            }
            "lipschitz" => {
                let [k, f] = args else { return Err(self.err(e, "usage: (lipschitz k field) promises |grad field| <= k in x y z (exact distance: 1); proofs trust it")) };
                let (k, f) = (self.field(k, env)?, self.field(f, env)?);
                if self.reads_moving_input(&[k, f], env) {
                    return Ok(Value::Field(f, None));
                }
                return Ok(Value::Field(self.g.lipschitz(f, k), None));
            }
            "attr-of" => {
                let usage = "usage: (attr-of PATH NAME): the attribute NAME of the body placed at PATH (or the body PATH), as the body writes it: a record, string or number";
                let [p, k] = args else { return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))) };
                let (path, k) = (self.path_of(p)?, self.name_of(k)?);
                let items = self.items(&path).map_err(|m| self.err(p, m))?;
                let [it] = items.as_slice() else {
                    return Err(self.err_hint(p, format!("'{path}' places {} bodies; (attr-of PATH NAME) reads one", items.len()), Some(usage.into())));
                };
                let body = it.instance.as_ref().map_or(path.clone(), |(b, _)| b.clone());
                // Placed as the part is, as (attr NAME) reads it: a field attribute where the part lies.
                let attrs = self.placed_attrs(it);
                return match attrs.iter().find(|(n, _)| *n == k) {
                    Some((_, v)) => Ok(v.clone()),
                    None => {
                        let names: Vec<&str> = attrs.iter().map(|(n, _)| n.as_str()).collect();
                        Err(self.err_hint(p, format!("body '{body}' (at '{path}') has no attribute '{k}'"), Some(format!("its attributes: {}", names.join(", ")))))
                    }
                };
            }
            "attr" => {
                let usage = "(attr NAME) is an attribute of the body being integrated (its own, else its material's), placed with it: it is read in the fields of an \
                             integral or mean claim over a body, an assembly or an instance path, e.g. (mean root z (weight (attr density)) < 80)";
                let [k] = args else { return Err(self.err_hint(e, "usage: (attr NAME)", Some(usage.into()))) };
                let k = self.name_of(k)?;
                let Some((label, attrs)) = self.integrand_of.clone() else { return Err(self.err_hint(e, "(attr NAME) outside an integrand", Some(usage.into()))) };
                let attrs = match attrs {
                    Ok(attrs) => attrs,
                    Err(why) => return Err(self.err_hint(e, format!("the target {label} {why}"), Some(usage.into()))),
                };
                return match materials::required(&mut self.g, &label, &attrs, &k) {
                    Ok(f) => {
                        self.want(&[f]);
                        Ok(Value::Field(f, None))
                    }
                    Err(m) => Err(self.err_hint(e, m.message, m.hint)),
                };
            }
            "grid-bake" | "import-stl" | "redistance" | "bounds" | "profile-reach" => return self.eval_services(name, e, args, env),
            "material" | "by-body" => return self.eval_materials(name, e, args, env),
            "check" => return self.eval_check(e, args, env),
            "ranged" => {
                let usage = "usage: (ranged lo hi [at]): a new variable, which every claim reading it ranges over [lo, hi] (a range of the question: \
                             how far a tool is withdrawn, how far a view looks), and no run pins (--nominal, --param); at (lo if not given) where nothing ranges (eval, mesh, view)";
                if !(2..=3).contains(&args.len()) {
                    return Err(self.err_hint(e, "bad ranged", Some(usage.into())));
                }
                let ends = args.iter().map(|x| self.eval(x, env)).collect::<R<Vec<_>>>()?;
                let said = self.one_said("ranged", &ends, e)?;
                let end = |i: usize| -> R<f64> {
                    let f = self.as_field(ends[i].clone(), &args[i])?;
                    self.g.as_const(f).filter(|v| v.is_finite()).ok_or_else(|| self.err_hint(&args[i], "a range's ends are constant finite numbers", Some(usage.into())))
                };
                let (lo, hi) = (end(0)?, end(1)?);
                let at = if args.len() == 3 { end(2)? } else { lo };
                if !(lo <= at && at <= hi) {
                    return Err(self.err_hint(e, format!("(ranged {lo} {hi} {at}) needs lo <= at <= hi"), Some(usage.into())));
                }
                // Named by its form, "t-slot-reaches 0", the first number free:
                // no declared name has a space (a def that is one renames it).
                let base = &self.form_name;
                let name = (0..).map(|i| format!("{base} {i}")).find(|n| !self.ranged.iter().any(|r| self.params[r.0].name == *n)).unwrap();
                // A value given it by that name (--param "b 0=2mm") says
                // its unit as one given a declared param does.
                let said = self.param_given(&name, said, e)?;
                self.params.push(Param { name, default: at, tol: None });
                self.ranged.push((self.params.len() - 1, (lo, hi)));
                if let Some(s) = said.clone() {
                    self.param_said.insert(self.params.len() - 1, s);
                }
                return Ok(Value::Field(self.g.var(Var::Param(self.params.len() as u32 - 1)), said));
            }
            "error" => {
                let usage = "usage: (error \"text\" [VALUE ...]): a load error saying text and the values after it (strings as written, numbers, records), where it is evaluated (a library refusing an input no part has)";
                let [m, rest @ ..] = args else { return Err(self.err_hint(e, "bad error", Some(usage.into()))) };
                let Value::Str(mut m) = self.eval(m, env)? else { return Err(self.err_hint(e, "an error's text is a string", Some(usage.into()))) };
                for x in rest {
                    let v = self.eval(x, env)?;
                    m += &self.value_text(&v);
                }
                return Err(self.err(e, m));
            }
            "met-by" => {
                let usage = "usage: (met-by PATH): the name of the frame the instance PATH was placed by (its CFRAME, met face to face with WHERE)";
                let [p] = args else { return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))) };
                let path = self.path_of(p)?;
                // Resolved as attr-of resolves a path, with its errors.
                self.items(&path).map_err(|m| self.err(p, m))?;
                let (parent, name) = path.rsplit_once('/').ok_or_else(|| self.err_hint(p, format!("'{path}' is no instance path"), Some(usage.into())))?;
                let parent = self.placement(parent).map(|(_, child)| child);
                let met = parent.as_deref().and_then(|a| self.assembly(a)).and_then(|a| a.instances.iter().find(|i| i.name == name)).map(|i| i.met.clone());
                return match met {
                    Some(Some(n)) => Ok(Value::Str(n)),
                    Some(None) => Err(self.err_hint(p, format!("'{path}' is placed by no frame of its own (no CFRAME)"), Some(usage.into()))),
                    None => Err(self.err_hint(p, format!("'{path}' is no instance path"), Some(usage.into()))),
                };
            }
            "relation" => return self.eval_relation(e, args, env),
            "polygon" => return self.eval_polygon(e, args, env),
            "revolve-z" => return self.eval_revolve(e, args, env),
            "study-field" => return self.eval_result(e, args, env),
            "at" => return Ok(Value::Frame(self.eval_at(args, env)?)),
            "axes" => return Ok(Value::Frame(self.eval_axes(e, args, env)?)),
            "frame-of" => {
                let usage = "usage: (frame-of PATH [NAME]): the frame NAME declared by the body or assembly PATH ends at, placed; NAME omitted, PATH's placement";
                let (p, n) = match args {
                    [p] => (p, None),
                    [p, n] => (p, Some(self.name_of(n)?)),
                    _ => return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))),
                };
                let path = self.name_of(p)?;
                return self.frame_of(&path, n.as_deref()).map(Value::Frame).map_err(|m| self.err(e, m));
            }
            "place" => {
                let [v, f] = args else { return Err(self.err_hint(e, "wrong arguments", Some("usage: (place SHAPE F): SHAPE, written in F's coordinates, put at F".into()))) };
                let v = self.eval(v, env)?;
                let f = self.frame_value(f, env, "place's second argument")?;
                return Ok(self.to_world(&v, &f));
            }
            "shape-of" => {
                let usage = "usage: (shape-of name-or-path [construction | (except PATH ...)]): an instance's world field, that of a named construction inside it, or its parts but those under each PATH; a (for var list PATH ...) among them repeats its paths";
                // (except PATH ...): the parts under the path but those at or under each PATH,
                // a `for` among them repeating its paths for a list known at load time.
                let except = match args {
                    [_, n] => n.list().filter(|l| l.first().and_then(|w| w.sym()) == Some("except")).map(|l| l[1..].to_vec()),
                    _ => None,
                };
                let except = match except {
                    Some(left) => Some(self.expand_for_in(&left, env)?),
                    None => None,
                };
                let (p, part) = match args {
                    [p] => (p, None),
                    [p, _] if except.is_some() => (p, None),
                    [p, n] => (p, Some(self.name_of(n)?)),
                    _ => return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))),
                };
                let path = self.path_of(p)?;
                if part.is_none()
                    && except.is_none()
                    && let Some(&f) = self.shape_of.get(&path)
                {
                    return Ok(Value::Field(f, None));
                }
                // What it was made from (a body's shape, placements): what
                // reads it rests on that too.
                let mut items = self.items(&path).map_err(|m| self.err(p, m))?;
                if let Some(left) = &except {
                    if left.is_empty() {
                        return Err(self.err_hint(e, "(except) names no path", Some(usage.into())));
                    }
                    // At or under `out`: the label itself, or one below it past a '/'.
                    let under = |label: &str, out: &str| label == out || label.strip_prefix(out).is_some_and(|r| r.starts_with('/'));
                    let mut outs: Vec<String> = Vec::new();
                    for x in left {
                        let out = self.path_of(x)?;
                        if let Some(o) = outs.iter().find(|o| under(&out, o) || under(o, &out)) {
                            let what = if *o == out {
                                "is named twice".to_string()
                            } else if under(&out, o) {
                                format!("is under '{o}', already left out")
                            } else {
                                format!("holds '{o}', named too")
                            };
                            return Err(self.err_hint(x, format!("(except ...): '{out}' {what}"), Some("name each part, or the assembly holding it, once".into())));
                        }
                        if !items.iter().any(|it| under(&it.label, &out)) {
                            let hint = format!("a path from '{path}' down to a part it places, as {}", items.first().map_or("root/part", |it| it.label.as_str()));
                            return Err(self.err_hint(x, format!("'{out}' names no part of '{path}'"), Some(hint)));
                        }
                        outs.push(out);
                    }
                    items.retain(|it| !outs.iter().any(|o| under(&it.label, o)));
                    if items.is_empty() {
                        return Err(self.err_hint(
                            e,
                            format!("(except ...) leaves no part of '{path}'"),
                            Some("a claim against nothing holds whatever it asks: leave out less".into()),
                        ));
                    }
                }
                let mut shapes = Vec::new();
                for it in &items {
                    let Some(part) = &part else {
                        shapes.push(self.placed_field(it));
                        continue;
                    };
                    let Some((body, xf)) = it.instance.clone() else { return Err(self.err(p, format!("'{path}' is not a body or an instance"))) };
                    let local = self.body(&body).unwrap().shape;
                    let found = self.construction_names(part).iter().map(|n| self.named_in(local, n)).find(|f| !f.is_empty()).unwrap_or_default();
                    if found.is_empty() {
                        let names = self.named_in_all(local);
                        return Err(self.err_hint(p, format!("body '{body}' has no piece built as (tag {part} ...)"), Some(format!("its pieces: {}", names.join(", ")))));
                    }
                    // Every occurrence (a patterned def appears once per copy), placed like the instance.
                    let f = self.reduce(BinOp::Min, &found);
                    let map = self.inverse_map(&xf);
                    let world = self.g.subst(f, &map);
                    self.record_placed(world, f, &xf);
                    shapes.push(world);
                }
                // Of no items (an empty assembly), the empty shape.
                let f = self.reduce(BinOp::Min, &shapes);
                if part.is_none() && except.is_none() {
                    self.shape_of.insert(path, f);
                }
                return Ok(Value::Field(f, None));
            }
            _ => {}
        }
        // A local name bound to a function value (e.g. an argument of a lambda).
        if let Some((_, Value::Fn(f))) = env.iter().rev().find(|(k, _)| k == name) {
            let f = f.clone();
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.call(&f, vals, e);
        }
        if builtin(name).is_some() || COMPARISONS.contains(&name) {
            if let ("sin" | "cos" | "tan", [x]) = (name, args) {
                self.degrees_warning(name, x);
            }
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.apply_builtin(name, vals, e);
        }
        // A package's value (`(reaction PATH)`), read in this scope.
        if let Some((p, h)) = self.registry.value(name) {
            let mut r = Reader::in_env(self, env.clone());
            return p.value(&mut r, e, h, args);
        }
        // Only here is this head a user call. Resolve it in the file that
        // wrote it: syntax readers above (and placements, claims, packages)
        // must see their own words, even if a module defines one as a fn.
        let key = self.scoped(name);
        if let Some(mut err) = self.later(name, head) {
            // (def t-span (t-span 0 2)): the def's own name, read as a call.
            err.hint = self.top_form_hint(e, "an expression").or(err.hint);
            return Err(err);
        }
        if let Some(f) = self.fns.get(&*key).cloned() {
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.call(&Callable::User(name.to_string(), f), vals, e);
        }
        // A def whose value is a function is called as a local name holding one is.
        if let Some(Value::Fn(f)) = self.def(&key).cloned() {
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.call(&f, vals, e);
        }
        if let Some(v) = self.def(&key) {
            return Err(self.err(head, format!("'{name}' is a def holding a {}, not a function", v.kind())));
        }
        let mut names: Vec<String> = self.fns.keys().filter(|n| !n.contains(':')).cloned().collect();
        names.extend(SPECIAL.iter().map(|s| s.to_string()));
        names.extend(BUILTINS.iter().map(|(n, _)| n.to_string()));
        names.sort();
        let hint = if let Some(h) = self.hidden_hint(name) {
            Some(h)
        } else if self.registry.top_forms().split(", ").any(|f| f == name) {
            // A top-level form's head inside an expression.
            let here = match name {
                "fn" => "; a function value here is (lambda (args) body)",
                "def" | "param" => "; a local name here is (let ((name value)) body)",
                _ => "",
            };
            Some(format!("'{name}' is a top-level form, written only at a file's top level (`ink spec --section \"Top-level forms\"`){here}"))
        } else if let Some(axis) = name.strip_prefix("rot-").filter(|a| ["x", "y", "z"].contains(a)) {
            Some(format!("'{name}' is a placement step, written inside (at ...); a shape turns by std's (rotate-{axis} f deg)"))
        } else if !self.imported.contains_key("std") && in_std(name) {
            Some(format!("'{name}' is in the standard library: add (import \"std\")"))
        } else {
            unknown_function_hint(name, &names, self.imported.contains_key("std"))
        };
        Err(self.err_hint(head, format!("unknown function '{name}'"), hint))
    }
}

/// Doc::reach of `f` in graph `g` with window (`wlo`, `whi`).
pub fn reach_in(g: &ink_core::graph::Graph, window: (V3, V3), f: NodeId, ctx: &Ctx) -> Reach {
    reach_of(&std::sync::Arc::new(Tape::new(g, f)), window, ctx)
}

/// Halvings of a ranged param's span `extent_in_pieces` may make.
const REACH_PIECES: u32 = 6;

/// The extent of `tape` over ctx's scope as the union of its extents over
/// pieces of it: each ranged param it reads halved (the widest first, its
/// halves each found again) where a piece's is unbounded, up to `depth`
/// halvings. Over a whole turn a turned part's far points bound as any
/// direction, so no direction is shown outside (unbounded); over pieces of
/// the turn each is. The union holds every piece's material, so all of
/// it. Unbounded where a piece still is at the last halving, or no param
/// is read.
fn extent_in_pieces(tape: &std::sync::Arc<Tape>, ctx: &Ctx, window: (V3, V3), depth: u32) -> crate::query::Extent {
    use crate::query::Extent;
    let ks: Vec<u32> = ctx
        .p_spans
        .iter()
        .filter(|s| !s.derived && s.span.0 < s.span.1 && tape.ssa.iter().any(|i| matches!(i.op, ink_core::tape::Op::Var(Var::Param(k)) if k == s.k)))
        .map(|s| s.k)
        .collect();
    let Some(&k) = ks.iter().max_by(|&&a, &&b| {
        let w = |k: u32| ctx.p_spans.get(k).map_or(0.0, |s| s.1 - s.0);
        w(a).total_cmp(&w(b))
    }) else {
        return Extent::Unbounded;
    };
    if depth == 0 {
        return Extent::Unbounded;
    }
    let (lo, hi) = ctx.p_spans.get(k).expect("a ranged param");
    let m = ink_core::eval::span_mid((lo, hi));
    if !(lo < m && m < hi) {
        return Extent::Unbounded;
    }
    let mut out = Extent::Empty;
    for half in [(lo, m), (m, hi)] {
        let c = ctx.with_param(k, half);
        let e = match crate::query::extent(tape, &c, window) {
            Extent::Unbounded => extent_in_pieces(tape, &c, window, depth - 1),
            e => e,
        };
        out = match (out, e) {
            (Extent::Unbounded, _) | (_, Extent::Unbounded) => return Extent::Unbounded,
            (Extent::Empty, e) => e,
            (o, Extent::Empty) => o,
            (Extent::Bounded(a, b), Extent::Bounded(c, d)) => {
                let (lo, hi) = ink_core::v3::box_union((a, b), (c, d));
                Extent::Bounded(lo, hi)
            }
        };
    }
    out
}

/// Whether the first of the finest pieces `extent_in_pieces` may reach (each
/// param it reads at the first 2^-depth of its span) has a bounded extent:
/// where even that one has none, no union of pieces does, and the search
/// for one (each extent's own up to its largest budget, down every halving)
/// costs as much as the claim (a belt's shape driven along its loop).
fn finest_piece_bounded(tape: &std::sync::Arc<Tape>, ctx: &Ctx, window: (V3, V3), depth: u32) -> bool {
    let mut c = ctx.clone();
    let mut any = false;
    for s in ctx.p_spans.iter().filter(|s| !s.derived && s.span.0 < s.span.1 && tape.ssa.iter().any(|i| matches!(i.op, ink_core::tape::Op::Var(Var::Param(k)) if k == s.k))) {
        let (lo, hi) = s.span;
        c = c.with_param(s.k, (lo, lo + (hi - lo) * 2f64.powi(-(depth as i32))));
        any = true;
    }
    any && !matches!(crate::query::extent(tape, &c, window), crate::query::Extent::Unbounded)
}

/// Doc::reach of the field `tape` computes, with window (`wlo`, `whi`).
pub fn reach_of(tape: &std::sync::Arc<Tape>, (wlo, whi): (V3, V3), ctx: &Ctx) -> Reach {
    let extent = match crate::query::extent(tape, ctx, (wlo, whi)) {
        crate::query::Extent::Unbounded if finest_piece_bounded(tape, ctx, (wlo, whi), REACH_PIECES) => extent_in_pieces(tape, ctx, (wlo, whi), REACH_PIECES),
        e => e,
    };
    match extent {
        crate::query::Extent::Bounded(lo, hi) if (0..3).any(|k| lo[k] < wlo[k] || hi[k] > whi[k]) => Reach { lo, hi, whole: true },
        // No material anywhere, or all of it inside the window.
        crate::query::Extent::Empty | crate::query::Extent::Bounded(..) => Reach { lo: wlo, hi: whi, whole: true },
        crate::query::Extent::Unbounded => Reach { lo: wlo, hi: whi, whole: false },
    }
}

/// The kernel's own clauses of a relation (zone, claim, data) and of an
/// instance (its joint and limits): never `(optional ...)`, in either place.
const KERNEL_CLAUSES: &[&str] = &["zone", "claim", "data", "turn", "slide", "limits", "speed-limits"];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_load_time_scalar_runs_only_its_provider_dependencies() {
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        #[derive(Clone)]
        struct Counted {
            name: &'static str,
            owns: [usize; 1],
            read: usize,
            calls: Arc<AtomicUsize>,
        }
        impl Form for Counted {
            fn name(&self) -> &str {
                self.name
            }
            fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
                Some(Arc::new(self.clone()))
            }
        }
        impl ParamProvider for Counted {
            fn name(&self) -> &str {
                self.name
            }
            fn owns(&self) -> &[usize] {
                &self.owns
            }
            fn uses_t(&self) -> bool {
                true
            }
            fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                params[self.owns[0]] = 2.0 * params[self.read] + t;
                Ok(None)
            }
        }
        let mut d = Doc::load_str_in(&CORE, "(param p0 3) (param p1 3) (param p2 3) (param p3 3)", "test", Path::new("."), &[]).unwrap();
        let calls: Vec<_> = (0..3).map(|_| Arc::new(AtomicUsize::new(0))).collect();
        for (j, (name, read)) in [("first", 0), ("second", 1), ("unrelated", 0)].into_iter().enumerate() {
            let field = d.g.var(Var::Param(read as u32));
            let form = parse(&format!("(probe {name})")).unwrap().remove(0);
            d.keep_form("probe", "probe", &form, Arc::new(Counted { name, owns: [j + 1], read, calls: calls[j].clone() }), Reads { fields: vec![field], forms: vec![] }).unwrap();
        }
        calls.iter().for_each(|c| c.store(0, Ordering::Relaxed));
        d.provided.borrow_mut().clear();
        let second = d.g.var(Var::Param(2));
        assert_eq!(d.load_eval(second), 12.0);
        d.overrides = vec![("p0".into(), 5.0.into())];
        assert_eq!(d.load_eval(second), 20.0);
        assert_eq!(d.load_eval(second), 20.0); // unchanged inputs do no work
        let constant = d.g.constant(7.0);
        assert_eq!(d.load_eval(constant), 7.0);
        for (head, field) in [("bounds", "(max (- (abs x) 1) (- (abs y) 1) (- (abs z) 1))"), ("profile-reach", "(max (- (abs x) 1) (- (abs y) 1))")] {
            let f = parse(&format!("({head} {field})")).unwrap().remove(0);
            d.eval_services(head, &f, &f.list().unwrap()[1..], &mut Vec::new()).unwrap();
        }
        assert_eq!(calls.iter().map(|c| c.load(Ordering::Relaxed)).collect::<Vec<_>>(), [2, 2, 0]);
        // A solver's owned values are its warm-start inputs too.
        d.params[1].default = 8.0;
        assert_eq!(d.load_eval(second), 20.0);
        assert_eq!(calls[0].load(Ordering::Relaxed), 3);
        let at = d.ctx(&[], 1.0).unwrap();
        assert_eq!(at.params[2], 15.0);
        let counts: Vec<_> = calls.iter().map(|c| c.load(Ordering::Relaxed)).collect();
        assert_eq!(d.ctx(&[], 1.0).unwrap().params[2], 15.0);
        assert_eq!(calls.iter().map(|c| c.load(Ordering::Relaxed)).collect::<Vec<_>>(), counts);
        assert_eq!(d.ctx(&[("p0".into(), 7.0.into())], 1.0).unwrap().params[2], 31.0);
    }

    #[test]
    fn a_provider_on_demand_runs_only_where_its_params_are_read() {
        // A published mass's density: not measured at load, nor in a
        // context nothing reading it asks for (check, render, distance);
        // wanted (a mass weighed, a claim over a density), measured once:
        // what it set before is not an input to it.
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        #[derive(Clone)]
        struct Lazy {
            calls: Arc<AtomicUsize>,
        }
        impl Form for Lazy {
            fn name(&self) -> &str {
                "lazy"
            }
            fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
                Some(Arc::new(self.clone()))
            }
        }
        impl ParamProvider for Lazy {
            fn name(&self) -> &str {
                "lazy"
            }
            fn owns(&self) -> &[usize] {
                &[1]
            }
            fn uses_t(&self) -> bool {
                false
            }
            fn on_demand(&self) -> bool {
                true
            }
            fn warm_starts(&self) -> bool {
                false
            }
            fn provide(&self, params: &mut [f64], _: f64) -> Result<Option<String>, String> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                params[1] = 4.0;
                Ok(None)
            }
        }
        let mut d = Doc::load_str_in(&CORE, "(param p0 3) (param density 0)", "test", Path::new("."), &[]).unwrap();
        let calls = Arc::new(AtomicUsize::new(0));
        let form = parse("(probe lazy)").unwrap().remove(0);
        d.keep_form("probe", "probe", &form, Arc::new(Lazy { calls: calls.clone() }), Reads::default()).unwrap();
        assert_eq!((d.load_ctx().params[1], d.ctx(&[], 0.0).unwrap().params[1]), (0.0, 0.0));
        assert_eq!(calls.load(Ordering::Relaxed), 0, "provided at load or for a context nothing reading it asks for");
        let density = d.g.var(Var::Param(1));
        d.want(&[density]);
        let at: Vec<f64> = [0.0, 0.5].iter().map(|&t| d.ctx(&[], t).unwrap().params[1]).collect();
        let mut again = d.load_params();
        d.provide_one(0, &mut again, 0.0).unwrap();
        d.provide_one(0, &mut again, 0.0).unwrap();
        assert_eq!((at, again[1]), (vec![4.0, 4.0], 4.0));
        assert_eq!(calls.load(Ordering::Relaxed), 1, "measured more than once");
    }

    #[test]
    fn a_provider_that_does_not_read_t_runs_once_over_a_motion() {
        // A published mass's density, at each moment a GLB or a motion
        // claim samples: measured once, not at every moment.
        use std::sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        };
        #[derive(Clone)]
        struct Counted {
            uses_t: bool,
            owns: [usize; 1],
            calls: Arc<AtomicUsize>,
        }
        impl Form for Counted {
            fn name(&self) -> &str {
                "counted"
            }
            fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
                Some(Arc::new(self.clone()))
            }
        }
        impl ParamProvider for Counted {
            fn name(&self) -> &str {
                "counted"
            }
            fn owns(&self) -> &[usize] {
                &self.owns
            }
            fn uses_t(&self) -> bool {
                self.uses_t
            }
            fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
                self.calls.fetch_add(1, Ordering::Relaxed);
                params[self.owns[0]] = if self.uses_t { t } else { 4.0 };
                Ok(None)
            }
        }
        for uses_t in [false, true] {
            let mut d = Doc::load_str_in(&CORE, "(param p0 3)", "test", Path::new("."), &[]).unwrap();
            let calls = Arc::new(AtomicUsize::new(0));
            let form = parse("(probe counted)").unwrap().remove(0);
            d.keep_form("probe", "probe", &form, Arc::new(Counted { uses_t, owns: [0], calls: calls.clone() }), Reads::default()).unwrap();
            calls.store(0, Ordering::Relaxed);
            d.provided.borrow_mut().clear();
            let at: Vec<f64> = [0.0, 0.25, 0.5].iter().map(|&t| d.ctx(&[], t).unwrap().params[0]).collect();
            assert_eq!(at, if uses_t { [0.0, 0.25, 0.5] } else { [4.0; 3] });
            assert_eq!(calls.load(Ordering::Relaxed), if uses_t { 3 } else { 1 }, "uses_t {uses_t}");
        }
    }

    /// A directory of a test's own (`ink-TAG-PID` in the system's temporary
    /// directory), removed when dropped, a test that panics too.
    struct Scratch(std::path::PathBuf);

    impl Scratch {
        fn new(tag: &str) -> Scratch {
            let dir = std::env::temp_dir().join(format!("ink-{tag}-{}", std::process::id()));
            _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            Scratch(dir)
        }
    }

    impl std::ops::Deref for Scratch {
        type Target = Path;
        fn deref(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            _ = std::fs::remove_dir_all(&self.0);
        }
    }

    /// The language alone, with no packages.
    static CORE: Registry = Registry::new(&[]);

    #[test]
    fn a_json_file_is_read_once_and_shared() {
        // A table read by several files is one document in memory: every
        // data value read from it, and `get` into it, is a place in it, not
        // a copy (each catalog part held a whole copy of its family's
        // table: 2.4 MB a part).
        let dir = Scratch::new("json-shared");
        std::fs::write(dir.join("table.json"), r#"{"sizes": {"a/b": {"d": 3, "tags": [{"k": 1}]}}}"#).unwrap();
        let reader = r#"(def table (read-json "table.json")) (def row (get (get table "sizes") "a/b")) (def d (get row "d")) (def first (nth (get row "tags") 0))"#;
        std::fs::write(dir.join("a.ink"), reader).unwrap();
        std::fs::write(dir.join("b.ink"), reader).unwrap();
        let d = Doc::load_str_in(&CORE, r#"(import "a.ink" as p) (import "b.ink" as q)"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let data = |n: &str| match d.def(n) {
            Some(Value::Data(x)) => x.clone(),
            v => panic!("{n}: {v:?}"),
        };
        let all = [data("p.table"), data("q.table"), data("p.row"), data("q.row"), data("q.first")];
        assert!(all.iter().all(|x| Rc::ptr_eq(&x.root, &all[0].root)));
        assert_eq!(all[3].value()["d"], 3);
        assert_eq!(all[4].value()["k"], 1);
        assert!(matches!(d.def("q.d"), Some(Value::Field(f, _)) if d.g.as_const(*f) == Some(3.0)));
    }

    #[test]
    fn get_takes_a_default_where_the_key_is_missing() {
        // A record's optional key (a belt's working tension, per its width
        // or whole) is read with a default, evaluated only where the key is
        // missing; without one, the error lists the keys there are.
        let dir = Scratch::new("json-default");
        std::fs::write(dir.join("belt.json"), r#"{"whole": 600, "at": "x"}"#).unwrap();
        let d = Doc::load_str_in(&CORE, r#"(def b (read-json "belt.json")) (def w (get b "whole" (get b "per-mm"))) (def p (get b "per-mm" 7))"#, "t", &dir, &[])
            .unwrap_or_else(|e| panic!("{e:?}"));
        let num = |n: &str| match d.def(n) {
            Some(Value::Field(f, _)) => d.g.as_const(*f),
            v => panic!("{n}: {v:?}"),
        };
        assert_eq!((num("w"), num("p")), (Some(600.0), Some(7.0)));
        let e = Doc::load_str_in(&CORE, r#"(def w (get (read-json "belt.json") "per-mm"))"#, "t", &dir, &[]).err().expect("a missing key");
        assert!(e.message == "no key 'per-mm'" && e.hint.as_deref() == Some("its keys: whole, at"), "{e:?}");
    }

    #[test]
    fn a_function_reads_paths_beside_its_own_file() {
        // A path a function reads is relative to the file it is written
        // in, wherever it is called from: a caller with a table.json of its
        // own still gets the library's (docs/IR.md, names are lexical).
        let dir = Scratch::new("fn-path");
        std::fs::create_dir_all(dir.join("lib")).unwrap();
        std::fs::write(dir.join("lib/lookup.ink"), r#"(fn lookup (k) (get (read-json "table.json") k))"#).unwrap();
        std::fs::write(dir.join("lib/table.json"), r#"{"member": 7}"#).unwrap();
        std::fs::write(dir.join("table.json"), r#"{"member": 99}"#).unwrap();
        let d = Doc::load_str_in(&CORE, r#"(import "lib/lookup.ink" as l) (def v (l.lookup "member"))"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
        assert!(matches!(d.def("v"), Some(Value::Field(f, _)) if d.g.as_const(*f) == Some(7.0)));
    }

    #[test]
    fn a_turn_keeps_its_axis_exact() {
        // A turn about z leaves z alone: its matrix's z row and column are
        // the constants 0 and 1, so a frame on a spinning part reads t only
        // where it moves. A z entry of (1 - c) + c would read t, and a
        // claim on a spinning gear's axis or centre (a relation's) would be
        // stepped through every turn of it instead of checked once.
        let src = "(body b (shape (- x 1))) (assembly root (instance a b (at (move 0 0 5)) (turn z (* 360 t))))";
        let d = Doc::load_str_in(&CORE, src, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let at = d.assembly("root").unwrap().instances[0].at;
        let c = |e: NodeId| d.g.as_const(e);
        assert_eq!([at.r[2][2], at.r[0][2], at.r[1][2], at.r[2][0], at.r[2][1]].map(c), [Some(1.0), Some(0.0), Some(0.0), Some(0.0), Some(0.0)]);
        assert_eq!(at.t.map(c), [Some(0.0), Some(0.0), Some(5.0)]);
    }

    #[test]
    fn a_family_file_runs_once_for_all_its_parts() {
        // Each part imports its family plainly, joining the part's
        // namespace. A family of values and functions runs once: under the
        // second part's namespace its names are bound to what the first
        // run defined (each part kept its own copy of every function of its
        // family). A family with a param runs for each part:
        // a param is set by its name, the part's.
        let dir = Scratch::new("family-once");
        std::fs::write(dir.join("family.ink"), "(fn twice (v) (* 2 v)) (def k (twice 4))").unwrap();
        std::fs::write(dir.join("sized.ink"), "(param size 3) (fn sized () size)").unwrap();
        for (i, n) in [(1, 1), (2, 5)] {
            std::fs::write(dir.join(format!("part{i}.ink")), format!("(import \"family.ink\") (import \"sized.ink\") (def v (+ (twice {n}) (sized) k))")).unwrap();
        }
        let d = Doc::load_str_in(&CORE, r#"(import "part1.ink" as p) (import "part2.ink" as q)"#, "t", &dir, &[]).unwrap_or_else(|e| panic!("{e:?}"));
        let v = |n: &str| match d.def(n) {
            Some(Value::Field(f, _)) => d.g.as_const(*f),
            v => panic!("{n}: {v:?}"),
        };
        assert_eq!((v("p.v"), v("q.v"), v("q.k")), (None, None, Some(8.0)));
        assert_eq!(d.scoped("q.twice"), d.scoped("p.twice"), "the family ran again");
        assert_eq!(d.fns.values().filter(|f| f.params == ["v"]).map(Rc::as_ptr).collect::<HashSet<_>>().len(), 1);
        let names: Vec<&str> = d.params.iter().map(|p| p.name.as_str()).collect();
        assert_eq!(names, ["p.size", "q.size"]);
        let ctx = |q: f64| Ctx::new(0.0, vec![3.0, q]);
        let at = |n: &str, c: &Ctx| match d.def(n) {
            Some(Value::Field(f, _)) => ink_core::eval::eval_point(&ink_core::tape::Tape::new(&d.g, *f), c, [0.0; 3]),
            v => panic!("{n}: {v:?}"),
        };
        assert_eq!((at("p.v", &ctx(3.0)), at("q.v", &ctx(4.0))), (13.0, 22.0));
        // A part's own name cannot shadow its family's, as when it ran.
        std::fs::write(dir.join("part3.ink"), "(import \"family.ink\") (def k 1)").unwrap();
        let e = Doc::load_str_in(&CORE, r#"(import "part1.ink" as p) (import "part3.ink" as r)"#, "t", &dir, &[]).err().expect("k defined twice");
        assert!(e.message.contains("already defined"), "{e:?}");
    }

    #[test]
    fn grids_over_package_results_key_by_their_values() {
        // A grid sampled from a package's result is keyed by what it reads
        // (no package is named): the result's values are part of the field's
        // fingerprint, so a solver that gives other values (here: a probe
        // package's version) misses, and one giving the same values hits.
        use std::sync::atomic::{AtomicU32, Ordering};
        static VERSION: AtomicU32 = AtomicU32::new(1);
        struct Probe;
        struct ProbeForm(String);
        struct Solved(f32);
        impl Package for Probe {
            fn name(&self) -> &'static str {
                "probe"
            }
            fn forms(&self) -> &'static [&'static str] {
                &["probe"]
            }
            fn read(&self, r: &mut Reader, _form: &Sexp, _head: &str, args: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
                Ok(std::sync::Arc::new(ProbeForm(r.name(&args[0])?)))
            }
        }
        impl Form for ProbeForm {
            fn name(&self) -> &str {
                &self.0
            }
            fn results(&self) -> Vec<ResultSpec> {
                vec![ResultSpec::scalar("v", "")]
            }
            fn run(&self, _doc: &mut Doc, _ctx: &Ctx, _opts: &RunOpts) -> R<Box<dyn Outcome>> {
                Ok(Box::new(Solved(VERSION.load(Ordering::SeqCst) as f32)))
            }
        }
        impl Outcome for Solved {
            fn field(&self, _: &str) -> Option<crate::grid::GridData> {
                Some(sampled_result(vec![self.0; 8], [2, 2, 2], [-5.0; 3], 10.0, 0.0))
            }
        }
        static PROBES: Registry = Registry::new(&[&Probe]);
        let key = |version: u32| {
            VERSION.store(version, Ordering::SeqCst);
            let src = "(probe p)\n(def f (max (- (abs x) 5) (study-field p v)))";
            let d = Doc::load_str_in(&PROBES, src, "test", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}"));
            let Some(Value::Field(f, _)) = d.def("f") else { panic!("f is not a field") };
            d.grid_key("grid-bake", Some(1.0), None, Some(*f), &[], &d.load_ctx().params, &[])
        };
        let (a, b) = (key(1), key(2));
        assert_ne!(a, b);
        assert_eq!(key(1), a);
    }

    /// A `for` among an instance's clauses reads its assembly's arguments:
    /// a drive is limited at one supply and not at another.
    #[test]
    fn an_instance_clause_for_reads_the_assembly_arguments() {
        let d = Doc::load_str_in(
            &CORE,
            "(import \"std\") (body a (shape (box 1 1 1)))
            (assembly (drive (angle 0) (volts 24))
              (instance r a (at) (turn z angle) (for lo (if (< volts 20) (vec) (vec -5)) (limits lo 5))))
            (assembly root (instance low (drive (* 10 t) 12)) (instance at (drive (* 10 t))))",
            "t",
            Path::new("."),
            &[],
        )
        .unwrap();
        let limited = |path: &str| d.assemblies.iter().find(|a| a.name == path).unwrap().instances[0].joint.as_ref().unwrap().limits.is_some();
        assert!(!limited("root/low") && limited("root/at"));
    }

    /// A package's clauses (a joint's effort limits and stiffness, a
    /// relation's coupling) are refused by an ink that loads no package
    /// owning them, saying so; `(optional CLAUSE)` is left out there, so a
    /// gear mesh's claims stand without its coupling.
    #[test]
    fn a_package_clause_is_refused_where_no_package_owns_it() {
        let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\") (body a (shape (box 1 1 1)))\n{src}"), "t", Path::new("."), &[]);
        for (src, head) in [
            ("(assembly root (instance r a (at) (turn z 0) (effort-limits -5 5)))", "'effort-limits'"),
            ("(assembly root (instance r a (at) (slide x 0) (stiffness 10)))", "'stiffness'"),
            ("(assembly root (instance p a) (instance q a (at (move 9 0 0))) (relate p q (belt-drive 5 5 10)))", "'couple'"),
        ] {
            let e = load(src).err().unwrap_or_else(|| panic!("{src} loads"));
            assert!(e.message.contains(head) && e.message.contains("a package this ink loads"), "{src}: {}", e.message);
        }
        let d = load(
            "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
               (relate p q (lambda (fa fb sa sb) (relation r (optional (couple (vec) 0)) (data n 1)))))",
        )
        .unwrap_or_else(|e| panic!("{e:?}"));
        let k = d.assemblies[0].relations[0].kind.clone().unwrap();
        assert!(k.clauses.is_empty() && k.data.len() == 1);
        // A joint's, as a catalog stepper writes them: the kernel places it
        // without them; a kernel clause is never optional.
        let d = load(
            "(assembly root (instance r a (at) (turn z 0) (optional (stiffness 10))
               (for c (vec 0 1) (optional (effort-limits -5 5)))))",
        )
        .unwrap_or_else(|e| panic!("{e:?}"));
        let j = d.assemblies[0].instances[0].joint.as_ref().unwrap();
        assert!(j.clauses.is_empty(), "{:?}", j.clauses.len());
        let e = load("(assembly root (instance r a (at) (optional (turn z 0))))").err().unwrap();
        assert!(e.message.contains("for a package's clause"), "{}", e.message);
        // What it expands to is checked: a kernel clause in a `for` is
        // still the kernel's, in an instance and in a relation.
        for src in [
            "(assembly root (instance r a (at) (turn z 5) (optional (for k (vec 0) (limits -1 1)))))",
            "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
               (relate p q (lambda (fa fb sa sb) (relation r (optional (for k (vec 0) (limits -1 1)))))))",
            "(assembly root (instance p a) (instance q a (at (move 9 0 0)))
               (relate p q (lambda (fa fb sa sb) (relation r (optional (for k (vec 0) (data n 1)))))))",
        ] {
            let e = load(src).err().unwrap_or_else(|| panic!("{src} loads"));
            assert!(e.message.contains("for a package's clause"), "{src}: {}", e.message);
        }
        // Where a loaded package reads joint or relation clauses, one none
        // reads (a misspelt stiffness) warns: dropping it there would drop
        // that package's check. Where none does (the kernel's ink, whose
        // packages read neither), the kernel stands without it, silently.
        struct Other(bool);
        impl Package for Other {
            fn name(&self) -> &'static str {
                if self.0 { "joints" } else { "forms-only" }
            }
            fn forms(&self) -> &'static [&'static str] {
                &[]
            }
            fn read(&self, _: &mut Reader, f: &Sexp, _: &str, _: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
                unreachable!("{f}")
            }
            fn joint_clauses(&self) -> &'static [&'static str] {
                if self.0 { &["stiffness"] } else { &[] }
            }
            fn relation_clauses(&self) -> &'static [&'static str] {
                if self.0 { &["couple"] } else { &[] }
            }
        }
        static OTHER: Registry = Registry::new(&[&Other(true)]);
        static FORMS_ONLY: Registry = Registry::new(&[&Other(false)]);
        let typo = "(import \"std\") (body a (shape (box 1 1 1)))
            (assembly root (instance p a) (instance q a (at (move 9 0 0)) (turn z 0) (optional (stifness 10)))
              (relate p q (lambda (fa fb sa sb) (relation r (optional (cuple (vec) 0))))))";
        let warned = |reg: &'static Registry| Doc::load_str_in(reg, typo, "t", Path::new("."), &[]).unwrap_or_else(|e| panic!("{e:?}")).load_warnings;
        let w = warned(&OTHER);
        assert!(w.iter().any(|w| w.contains("'stifness'")) && w.iter().any(|w| w.contains("'cuple'")), "{w:?}");
        assert!(warned(&CORE).is_empty() && warned(&FORMS_ONLY).is_empty());
    }

    /// The papercuts agents hit, each as
    /// they wrote it: the error says what to write.
    #[test]
    fn errors_say_what_to_write() {
        let e = Doc::load_str_in(&CORE, "(def a (sqr x))", "t", Path::new("."), &[]).err().unwrap();
        assert!(e.message.contains("unknown function 'sqr'"));
        assert!(e.hint.unwrap().contains("sqrt"));
        let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
        let hint = |src: &str| load(src).err().unwrap().hint.unwrap_or_default();
        // A CAD word gets std's name for it, not a near spelling ('hole'
        // is not 'fold'); a misspelt part of a name finds the name.
        let h = hint("(def a (hole 1 2))");
        assert!(h.contains("'subtract'") && h.contains("'cylinder-z'") && !h.contains("fold"), "{h}");
        assert!(hint("(def a (fillet 1 2))").contains("'fillet-union'"));
        assert!(hint("(def a (extrude 1 2))").contains("'extrude-z'"));
        assert!(hint("(def a (cylnder 1 2))").contains("'cylinder-z'"));
        assert_eq!(hint("(def a (foo 1 2))"), "");
        // Rounding a number is floor's; a top-level form inside an
        // expression says what is written there ('fn', not 'ln').
        assert!(hint("(def a (round 1.5))").contains("(floor (+ v 0.5))"));
        let h = hint("(def f (map (fn (a) (* 2 a)) (vec 1 2)))");
        assert!(h.contains("top-level form") && h.contains("(lambda (args) body)") && !h.contains("'ln'"), "{h}");
        // A function written as a def, three ways.
        assert_eq!(hint("(def (hole z) (move (cylinder-z 3 40) 0 6 z))"), "write (fn hole (z) body)");
        assert_eq!(hint("(def peg (x) (move (sphere 2) x 25 0))"), "write (fn peg (x) body)");
        assert_eq!(hint("(def slot (fn (x d) (move (box d 2 2) x 0 0)))"), "write (fn slot (x d) body)");
        // A def holding a function is called by name, as a local one is.
        let d = load("(def mark (lambda (z) (move (cylinder-z 1 1) 0 0 z)))\n(def marks (min (mark 10) (mark 20)))").unwrap_or_else(|e| panic!("{e:?}"));
        assert!(matches!(d.def("marks"), Some(Value::Field(..))));
        assert!(load("(def w 5)\n(def a (w 3))").err().unwrap().message.contains("'w' is a def holding a field, not a function"));
        // A reserved name lists them all.
        for n in ["t", "len", "ln", "floor", "fold", "at"] {
            let h = hint(&format!("(param {n} 3)"));
            assert!(h.starts_with("the built-in names: x y z t pi inf;") && h.split([' ', ';']).any(|w| w == n), "{n}: {h}");
        }
        // A name quoted or in parentheses.
        assert_eq!(hint("(body \"pen-cup\" (shape (box 1 1 1)))"), "a name is a bare word: pen-cup");
        assert_eq!(hint("(body b (shape (box 1 1 1)))\n(assembly root (for i (range 2) (instance (motor-$i) b)))"), "a name is a bare word: motor-0");
        // A placement step written as a shape's turn.
        assert!(hint("(def b (rot-y 90 (cylinder-z 8 55)))").contains("(rotate-y f deg)"));
        // A shape's move or turn written as a placement or a step of one
        // a step's misspelling,
        // a step nested in a step; a function
        // misspelt there keeps its own hint.
        let plate = "(body plate (shape (box 40 40 5)))\n(fn seat (i) (at (move i 0 5)))\n";
        for (src, says) in [
            ("(body axle (shape (sphere 2)) (frame left-axle (move -24 0 0)))", "(at STEP ...)"),
            ("(assembly root (instance p plate (move 0 0 5)))", "(at STEP ...)"),
            ("(assembly root (instance p plate (at (move 0 0 5) (rotate-z 90))))", "(rot-z deg)"),
            ("(assembly root (instance p plate (at (move 0 0 5) (rotz 90))))", "(rot-z deg)"),
            ("(assembly root (instance p plate (sead 1)))", "'read-json', 'seat'?"),
        ] {
            assert!(hint(&format!("{plate}{src}")).contains(says), "{src}: {}", hint(&format!("{plate}{src}")));
        }
        // The model's own function making a frame, miscounted, is told so
        // only.
        for src in ["(assembly root (instance p plate (seat)))", "(assembly root (instance p plate (at (seat 2 3))))"] {
            let e = load(&format!("{plate}{src}")).err().unwrap();
            assert!(e.message.starts_with("'seat' takes 1 argument(s) (i)") && e.hint.is_none(), "{src}: {e:?}");
        }
        // A shape written bare in a body.
        assert_eq!(hint("(body b (box 1 1 1))"), "write (shape (box ...))");
        // A one-letter name undeclared is no typo of t.
        assert_eq!(hint("(body plate (shape (box w 40 3)))"), "declare it: (param w VALUE) or (def w EXPR)");
        // A name declared as a body or an assembly further on, or an
        // instance, is not to be declared again.
        for src in
            ["(body b (shape (box root 1 1)))\n(assembly root (instance x b))", "(body p (shape (box 4 4 1)))\n(assembly root (instance a p) (instance b p (at (move a 0 0))))"]
        {
            assert!(!hint(src).contains("declare it"), "{src}: {}", hint(src));
        }
        let e = load(&format!("{plate}(assembly root (instance r plate (at (move 0 0 14 (move 0 0 4) (move 0 0 8)))))")).err().unwrap();
        assert!(e.message == "the step 'move' takes 3 number(s) (dx dy dz), got 5" && e.col == 38, "{e:?}");
    }
    /// std's sheet refuses at load what it cannot build, saying what to
    /// change: the sizes it cannot take, crossing bends, a bend written
    /// the wrong way round, a frame in a bend's zone.
    #[test]
    fn sheet_refuses_what_it_cannot_build() {
        let load = |src: &str| Doc::load_str_in(&CORE, &format!("(import \"std\")\n{src}"), "t", Path::new("."), &[]);
        let sheet = |th: &str, k: &str, bends: &str| {
            format!("(def S (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40))) (vec) {th} {k} (vec {bends})))\n(def f (sheet-folded S))")
        };
        let says = |src: &str| load(src).err().map(|e| e.message).unwrap_or_default();
        let l = "(sheet-bend (vec 30 40) (vec 30 0) 90 2)";
        assert!(load(&sheet("2", "0.4", l)).is_ok());
        assert!(says(&sheet("0", "0.4", l)).contains("thickness th must be above 0"));
        assert!(says(&sheet("2", "-1.5", l)).contains("K-factor k") && says(&sheet("2", "1.2", l)).contains("0 to 1"));
        assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 40) (vec 30 0) 90 0)")).contains("inside radius must be above 0"));
        assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 40) (vec 30 0) 190 2)")).contains("at most 180"));
        // bend 2 runs across bend 1's zone; its midpoint lies outside it
        let cross = format!("{l} (sheet-bend (vec 20 35) (vec 60 35) 90 2)");
        assert!(says(&sheet("2", "0.4", &cross)).contains("crosses another bend's zone"));
        // a lip on the upright written b -> a: its left faces the base
        let lip = format!("{l} (sheet-bend (vec 50 0) (vec 50 40) 90 2)");
        assert!(says(&sheet("2", "0.4", &lip)).contains("swap a and b"));
        assert!(says(&format!("{}\n(def fr (sheet-frame S 30 20 1))", sheet("2", "0.4", l))).contains("in a bend's zone"));
        assert!(says(&sheet("2", "0.4", "(sheet-bend (vec 30 20) (vec 30 20) 90 2)")).contains("a and b must differ"));
        // two tabs on one edge and a bend beyond both: which flange is it on?
        let two = "(sheet-bend (vec 30 10) (vec 30 0) 90 2) (sheet-bend (vec 30 40) (vec 30 30) 45 2) (sheet-bend (vec 45 40) (vec 45 0) 90 2)";
        assert!(says(&sheet("2", "0.4", two)).contains("a bend's line lies beyond two bends neither of which lies beyond the other"));
    }

    /// A number's unit goes with its value as it is computed, and the slot
    /// that uses the value checks it: through lets, calls, the branch an
    /// `if` takes, functions passed and bound, and sums; not through a
    /// branch not taken, nor a product (whose units combine).
    #[test]
    fn a_value_carries_the_unit_it_was_said_in() {
        static CORE: Registry = Registry::new(&[]);
        let fns = "(fn twice (p) (+ p p)) (fn choose (v) (if 1 v 1psi)) (fn other (v) (if 0 v 1psi)) (fn pass (f v) (f v)) (fn id (v) v)";
        let ids: String = (1..=17).map(|i| format!("(fn id{i} (v) ({} v))", if i == 1 { "id".into() } else { format!("id{}", i - 1) })).collect();
        let load = |mass: &str| {
            let src = format!("{fns} {ids} (param pp 50psi) (def r (ranged 1psi 2psi)) (body b (shape (- (max (abs x) (abs y) (abs z)) 1)) (attr mass {mass}))");
            Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).err().map(|e| e.message)
        };
        let nested = |x: &str| (0..23).fold(x.to_string(), |a, _| format!("(twice {a})"));
        for mass in [
            "(let ((p 50psi)) p)".to_string(),
            nested("1psi"),
            "(other 5g)".into(),
            "(pass id 1psi)".into(),
            "(let ((g id)) (g 1psi))".into(),
            "(id17 1psi)".into(),
            "(nth (map id (vec 1psi)) 0)".into(),
            "(reduce + (vec 1psi 2psi))".into(),
            "(- (max pp 0) 1)".into(),
            "r".into(),
        ] {
            let e = load(&mass).unwrap_or_else(|| panic!("{mass}: accepted"));
            assert!(e.starts_with("attr mass takes a mass") && e.contains("psi (t:1)"), "{mass}: {e}");
        }
        for mass in ["(choose 5g)", "(* 2 1psi)", &nested("1g")] {
            assert_eq!(load(mass), None, "{mass}");
        }
        // A sum's terms are one quantity.
        let e = load("(+ 1g 1psi)").unwrap();
        assert!(e.starts_with("'+' takes values of one quantity, but 1g (t:1) is a mass and 1psi (t:1) is a pressure"), "{e}");
    }

    /// A temperature in degF or K reads a point on its scale, so a change
    /// (a tolerance, an accuracy) in one is refused in a sum too, whichever
    /// term it is.
    #[test]
    fn a_sum_keeps_a_unit_with_a_zero_of_its_own() {
        for dev in ["(+ 1degC 33.8degF)", "(+ 33.8degF 1degC)", "(+ 1 1K)"] {
            let src = format!("(param temp 20degC (tol -1 {dev}))");
            let e = Doc::load_str_in(&CORE, &src, "t", Path::new("."), &[]).err().unwrap_or_else(|| panic!("{dev}: accepted"));
            assert!(e.message.starts_with("a tolerance's deviation is a change"), "{dev}: {}", e.message);
        }
        // A deviation with a unit is its default's quantity; a bare term
        // in a sum is in the language's unit.
        let e = Doc::load_str_in(&CORE, "(param dim 1mm (tol -1g 1g))", "t", Path::new("."), &[]).err().unwrap();
        assert!(e.message.starts_with("a tolerance's deviation takes a length"), "{}", e.message);
        assert!(Doc::load_str_in(&CORE, "(param dim 1mm (tol (- 0.1 1um) (+ 0.1 1um)))", "t", Path::new("."), &[]).is_ok());
    }

    /// A `for` makes numbers of its list's values: each keeps the unit its
    /// value was said in.
    #[test]
    fn a_for_keeps_its_values_units() {
        let load = |list: &str| Doc::load_str_in(&CORE, &format!("(for q {list} (body b (shape (- x 1)) (attr mass q)))"), "t", Path::new("."), &[]).err().map(|e| e.message);
        for list in ["(vec 50psi)", "(vec (+ 50psi 1))", "(let ((p 50psi)) (vec p))"] {
            let e = load(list).unwrap_or_else(|| panic!("{list}: accepted"));
            assert!(e.starts_with("attr mass takes a mass") && e.contains("written 50psi (t:1)"), "{list}: {e}");
        }
        assert_eq!(load("(vec 5g)"), None);
    }

    /// A --param value in a unit is checked as its default would be:
    /// against the default's quantity, and where the param is used.
    #[test]
    fn a_param_given_in_a_unit_is_said_in_it() {
        let dir = Scratch::new("param-units");
        let load = |default: &str, given: &str| {
            let f = dir.join("m.ink");
            std::fs::write(&f, format!("(param m {default})\n(body b (shape (- x 1)) (attr mass m))")).unwrap();
            let a = crate::command::Args::parse(&["eval".into(), f.display().to_string(), "--param".into(), format!("m={given}")]);
            crate::command::param_args(&a).and_then(|p| Doc::load_file_in(&CORE, &f, &p.values)).err().map(|e| e.message)
        };
        let says = |default: &str, given: &str, text: &str| {
            let e = load(default, given).unwrap_or_else(|| panic!("{default} {given}: accepted"));
            assert!(e.contains(text), "{default} {given}: {e}");
        };
        says("5g", "50psi", "--param m: m is a mass, written 5g (");
        says("5g", "50psi", "m.ink:1), but psi is a pressure");
        says("5", "50psi", "attr mass takes a mass");
        says("5", "50psi", "written 50psi (--param m)");
        says("5", "1psi:2psi", "attr mass takes a mass");
        says("5", "1g:2psi", "--param m: a range's ends are one quantity, but g is a mass and psi is a pressure");
        for given in ["50g", "50", "1g:2g", "1mg:2"] {
            assert_eq!(load("5g", given), None, "{given}");
        }
    }

    /// A value given a param keeps its unit wherever it travels: to a
    /// range by the name its form gives it, through a reload that changes
    /// the value (an optimizer's design), and across a range's ends (degF
    /// over degC, as in a sum, whichever end it is).
    #[test]
    fn a_given_unit_travels_with_its_value() {
        let dir = Scratch::new("given-travels");
        let run = |name: &str, src: &str, given: &str| {
            let f = dir.join(name);
            std::fs::write(&f, src).unwrap();
            let a = crate::command::Args::parse(&["eval".into(), f.display().to_string(), "--param".into(), given.into()]);
            crate::command::param_args(&a).and_then(|p| crate::view::ModelView::load_in(&CORE, &f, &p.values, 0.0))
        };
        // An unnamed range, by its form's name.
        let size = "(body b (shape (- (max (abs x) (abs y) (abs z)) (ranged 1mm 5mm))))";
        let e = run("size.ink", size, "b 0=500psi").err().expect("a pressure as a length");
        assert!(e.message.contains("--param b 0: b 0 is a length, written 1mm") && e.message.contains("but psi is a pressure"), "{}", e.message);
        assert!(run("size.ink", size, "b 0=0.1in").is_ok());
        // A reload changing the value keeps its unit: m is a pressure.
        let opt = "(param m 0)\n(def v (if (> m 1) m 2mm3))\n(body b (shape (- (max (abs x) (abs y) (abs z)) 0.5)))\n(assert a (volume (shape-of b) < v))";
        let v = run("opt.ink", opt, "m=50psi").unwrap_or_else(|e| panic!("{e:?}"));
        let e = v.reload(&[("m".into(), 3.0.into())]).err().expect("a pressure as a volume");
        assert!(e.message.starts_with("volume takes a volume") && e.message.contains("(--param m)"), "{}", e.message);
        assert!(run("opt.ink", opt, "m=50").unwrap().reload(&[("m".into(), 3.0.into())]).is_ok());
        // A range's ends: degF or K is kept over degC, either way round.
        let eps = "(param eps 1) (for v (vec eps) (assert a (value 1 > 0) (accuracy v)))";
        for given in ["eps=1degC:274.15K", "eps=274.15K:1degC"] {
            let e = run("eps.ink", eps, given).err().unwrap_or_else(|| panic!("{given}: accepted"));
            assert!(e.message.starts_with("an accuracy is a change, and K reads a point on its scale"), "{given}: {}", e.message);
        }
    }
}
