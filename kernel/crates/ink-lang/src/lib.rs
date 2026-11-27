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
//!
//! This file holds the document's types; what a Doc does is in a file each:
//! loading (load.rs), names (lookup.rs, words.rs), top-level forms
//! (top.rs), bodies and assemblies (build.rs, place.rs, transform.rs),
//! units (said.rs), claims (assert.rs), params (params.rs), targets
//! (targets.rs), reach (reach.rs) and expressions (expr.rs).

// The crates it builds on, at the paths the language's code uses.
use ink_core::{V3, eval, graph, tape};
mod error;
pub use error::{InkError, Open, Refused};
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
use std::path::{Path, PathBuf};
use std::rc::Rc;

pub mod answer;
mod assert;
mod author;
mod bridge;
mod build;
/// Conditions builders state on what they build (`check`).
mod carried;
pub mod check;
/// Checking a document's claims on a model view.
pub mod claims;
pub mod command;
pub mod cost;
mod estimate;
mod expr;
mod ext;
pub mod family;
pub mod follow;
mod font;
mod frame;
mod inexact;
pub mod interfere;
mod load;
pub mod look;
mod lookup;
mod magnitude;
pub use magnitude::Magnitude;
pub mod materials;
pub mod motion;
pub mod pairs;
mod params;
mod place;
mod profile;
mod rates;
mod reach;
mod reader;
mod records;
mod results;
mod said;
mod services;
pub mod sexp;
mod support;
mod targets;
mod top;
mod transform;
pub mod units;
pub use units::Given;
pub mod validity;
/// A loaded model and what analysis has learned of it.
pub mod view;
mod words;

pub use self::ext::*;
pub use self::follow::Follow;
pub use self::inexact::*;
use self::load::{canonical, index_read, not_read, sibling};
pub use self::load::{import_dirs, import_path, shown};
pub use self::reach::{Reach, plane_reach, reach_in, reach_of};
pub use self::reader::*;
pub use self::results::sampled_result;
use self::words::{
    BUILTINS, Builtin, CAD_WORDS, COMPARISONS, CONSTANTS, FORM_FUNCTIONS, SPECIAL, builtin, comment_hint, fn_written_as_def, has_point, in_std, is_reserved, near, reserved_hint,
    std_forms, unknown_function_hint,
};
pub use self::words::{builtin_words, choices, closest};

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

#[derive(Clone, Debug)]
pub enum Value {
    /// A number expression, and the number written with a unit its
    /// quantity was said by, if any (units::Said).
    Field(NodeId, Option<Rc<units::Said>>),
    /// A region whose exact membership is between two fields. Only an
    /// emptiness proof may consume the enclosure.
    Region {
        certain: NodeId,
        possible: NodeId,
    },
    /// A point (`probe`): three numbers, constants and params. Claims read
    /// it as a point (its sign in a shape, exactly); `place` puts it at a
    /// frame (R p + T, forward); it is no field, and no warp moves it.
    Point([NodeId; 3]),
    Str(String),
    /// A list, shared: a copy of the value is a copy of the pointer.
    Vec(Rc<Vec<Value>>),
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
        let owned = self.clauses.iter().flat_map(|c| c.reads.fields.iter().chain(&c.reads.own_clock).copied());
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
    assembly: Option<Rc<str>>,
}

/// Where a function reads the names it does not bind (lexical scope): the
/// file it is written in (an index into Doc::files; 0 is std). It sees the
/// names of that file, of the files it imports (SourceFile::bindings), and std's:
/// never one of the file calling it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Scope {
    pub file: u32,
}

impl Scope {
    pub const STD: Scope = Scope { file: 0 };
}

/// A lexical binding is a reference to the declaration itself. A source
/// symbol (including one containing a colon) can only reach one through
/// its file's bindings or std's.
#[derive(Clone)]
struct Binding {
    file: u32,
    line: u32,
    value: BoundName,
}

#[derive(Clone)]
enum BoundName {
    Param(usize),
    Def(usize),
    Text { param: usize, def: usize },
    Function(Rc<FnDef>),
}

impl Binding {
    fn def(&self) -> Option<usize> {
        match self.value {
            BoundName::Def(i) | BoundName::Text { def: i, .. } => Some(i),
            _ => None,
        }
    }
}

/// A file run: its name, its canonical path and the namespace it is
/// imported under (one file under two is run twice), its directory (the
/// paths it reads are relative to it), and the names it reads, as written
/// in it: its own and those of the files it imports, transitively (an
/// import `as p` adds them as `p.name`). `bindings`: each value's and
/// function's declaration; `objects`: each body's and assembly's name in the
/// document (`p.plate` for `plate` in a file imported as p: Doc::resolve).
/// `ahead`: where it defines each value, and `declared` every name its
/// forms declare (top_names: bodies and assemblies too). `given`: the
/// values its import gave its params (ImportValue). `own`: where it
/// defines a built-in's name (line, col: its definition is read from
/// there on), and `built_in_read` the line it first read each built-in at
/// (Doc::shadows).
struct SourceFile {
    span: Option<(f64, f64)>,
    name: Rc<str>,
    path: PathBuf,
    ns: Option<String>,
    dir: std::path::PathBuf,
    bindings: HashMap<String, Binding>,
    objects: HashMap<String, String>,
    ahead: HashMap<String, u32>,
    declared: HashSet<String>,
    given: Vec<ImportValue>,
    own: HashMap<String, (u32, u32)>,
    built_in_read: HashMap<String, u32>,
    /// The member of its family it was run as (family::member_of); None:
    /// the file itself (its `members`' first).
    member: Option<family::Member>,
    /// Its `(members ...)` read: one to a file.
    members_read: bool,
}

impl SourceFile {
    fn new(name: Rc<str>, path: PathBuf, ns: Option<String>, dir: std::path::PathBuf) -> SourceFile {
        SourceFile {
            span: None,
            name,
            path,
            ns,
            dir,
            bindings: HashMap::new(),
            objects: HashMap::new(),
            ahead: HashMap::new(),
            declared: HashSet::new(),
            given: Vec::new(),
            own: HashMap::new(),
            built_in_read: HashMap::new(),
            member: None,
            members_read: false,
        }
    }
}

#[derive(Clone)]
struct TolEnd {
    value: f64,
    said: Option<Rc<units::Said>>,
    line: u32,
    col: u32,
}

impl TolEnd {
    fn check(&self, rule: &units::UnitRequirement, name: &str) -> Result<(), String> {
        rule.check(name, &format!("tol bound {}", self.value), self.said.as_ref().is_some_and(|s| s.explicit_when.is_none()), self.value == 0.0)
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
    /// A text param's value, where it gives text (then `value` is NaN).
    text: Option<Rc<str>>,
    said: Option<Rc<units::Said>>,
    tol: Option<(f64, f64)>,
    file: Rc<str>,
    line: u32,
    col: u32,
    ends: Vec<TolEnd>,
    requirement: Option<units::UnitRequirement>,
}

impl ImportValue {
    fn check_units(&self, rule: &units::UnitRequirement, name: &str) -> R<()> {
        let error = |message, line, col| InkError { message, file: self.file.to_string(), line, col, hint: None, called_from: Vec::new(), open: None };
        rule.check(name, &format!("the import value {}", self.value), self.said.as_ref().is_some_and(|s| s.explicit_when.is_none()), self.value == 0.0)
            .map_err(|m| error(m, self.line, self.col))?;
        for end in &self.ends {
            end.check(rule, name).map_err(|m| error(m, end.line, end.col))?;
        }
        Ok(())
    }

    /// Refused, at the import, where its value or band is said in a unit
    /// of another quantity than `q`, its param's (`what`).
    fn fits(&self, q: Quantity, what: &str) -> R<()> {
        let mut units = self.said.iter().map(|s| (s, self.line, self.col)).chain(self.ends.iter().filter_map(|e| Some((e.said.as_ref()?, e.line, e.col))));
        match units.find_map(|(s, line, col)| Some((s.misfit(q, what)?, line, col))) {
            Some((message, line, col)) => Err(InkError { message, file: self.file.to_string(), line, col, hint: None, called_from: Vec::new(), open: None }),
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
    /// The string each of a vec of vecs starts with (the names of a part's `derives`).
    pub fn strings_first(list: &Value) -> Vec<String> {
        let Value::Vec(items) = list else { return Vec::new() };
        items.iter().filter_map(|item| if let Value::Vec(parts) = item { if let Some(Value::Str(s)) = parts.first() { Some(s.to_string()) } else { None } } else { None }).collect()
    }

    /// A list of values.
    pub fn list(items: Vec<Value>) -> Value {
        Value::Vec(Rc::new(items))
    }

    /// What it is, for messages ("field", "frame").
    pub fn kind(&self) -> &'static str {
        match self {
            Value::Field(..) => "field",
            Value::Region { .. } => "region bounded for empty",
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
    pub fn fields_into(&self, out: &mut Vec<NodeId>, closures: bool) {
        match self {
            Value::Field(f, _) => out.push(*f),
            Value::Region { certain, possible } => out.extend([certain, possible]),
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

/// A text param, `(param name "Ada" [(chars S)] [(len lo hi)])`: a string
/// read while the file loads (what a customer types), set by `--param
/// name=Grace` or an import, every value checked against what it allows.
#[derive(Clone, Debug)]
pub struct TextParam {
    pub name: String,
    pub default: String,
    /// The value this load read: the default, an import's or `--param`'s.
    pub value: String,
    /// The characters it may hold, where it says (`(chars S)`).
    pub chars: Option<String>,
    /// How many characters it may hold, where it says (`(len lo hi)`).
    pub len: Option<(usize, usize)>,
    /// Where an import set it ("b.ink:4").
    pub set_at: Option<String>,
}

impl TextParam {
    /// Why `v` is not a value it allows (`by`: who gave it), else None.
    /// Why `v` is not a text it allows, else None.
    pub fn refuses(&self, v: &str) -> Option<Refusal> {
        let n = v.chars().count();
        if let Some((lo, hi)) = self.len
            && !(lo..=hi).contains(&n)
        {
            return Some(Refusal::Length { n, lo, hi });
        }
        let allowed = self.chars.as_deref()?;
        let (at, c) = v.chars().enumerate().find(|(_, c)| !allowed.contains(*c))?;
        Some(Refusal::Char { at, c })
    }

    /// `refuses`, said: who gave `v` (`by`), and why the param refuses it.
    fn refused(&self, v: &str, by: &str) -> Option<String> {
        self.refuses(v).map(|r| format!("{by}: {}", r.said(self, v)))
    }
}

/// Why a text param does not allow a text (TextParam::refuses).
#[derive(Clone, Debug, PartialEq)]
pub enum Refusal {
    /// It holds `n` characters, not `lo` to `hi`.
    Length { n: usize, lo: usize, hi: usize },
    /// Its character `c`, the `at`-th (from 0), is not one it allows.
    Char { at: usize, c: char },
}

impl Refusal {
    /// What param `p` says of text `v` it refuses so.
    pub fn said(&self, p: &TextParam, v: &str) -> String {
        let v = serde_json::json!(v);
        match *self {
            Refusal::Length { n, lo, hi } => format!("{v} holds {n} characters; {} holds {lo} to {hi}", p.name),
            Refusal::Char { at, c } => format!(
                "{v} holds {} (character {}), which {} does not allow; it allows {}",
                serde_json::json!(c.to_string()),
                at + 1,
                p.name,
                serde_json::json!(p.chars.as_deref().unwrap_or_default())
            ),
        }
    }
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
            "param {name} is {v} with its band {a}..{b}, which excludes {v} (a basic size with a fit's deviations): require and interfere cover the band, but the commands that work at one value (distance, mass, eval, render, mesh) measure {v}, a size no real part has, and say so (outside_band). \
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

/// Every turn uses the same degree reduction and exact quarter turns.
pub(crate) fn turn_sin_cos(g: &mut Graph, deg: NodeId) -> [NodeId; 2] {
    let (angle, quadrant) = match g.as_const(deg).filter(|d| d.is_finite()) {
        Some(d) => {
            let r = d % 360.0;
            let q = (r / 90.0).round();
            // The remainder is exact, and subtraction from the nearest
            // quadrant is exact (within a factor of two, or subtracting zero).
            (g.constant(r - 90.0 * q), q.rem_euclid(4.0) as u8)
        }
        None => (deg, 0),
    };
    let [s, c] = g.turn_components(angle);
    match quadrant {
        0 => [s, c],
        1 => [c, g.neg(s)],
        2 => [g.neg(s), g.neg(c)],
        _ => [g.neg(c), s],
    }
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
    /// What the design states about the part placed here, `(fact WHAT
    /// [FRAME] VALUE)`: data, read by nothing until a `require` asks.
    pub facts: Vec<Stated>,
}

/// `(fact WHAT [WHERE] VALUE)` in an instance: a number the design states
/// about the part placed there, or with WHERE about a part under it
/// (`under`, the instances down to it) or one of its frames (`frame`), its
/// value read as any number is (params, units), and where it is written.
#[derive(Clone, Debug)]
pub struct Stated {
    pub what: String,
    pub under: Vec<String>,
    pub frame: Option<String>,
    pub value: NodeId,
    /// The language's unit of the quantity its value is, where ink kept
    /// one for it by its own rules (units::Said: a number written with its
    /// unit, a name, a sum of such, scaling by a bare number). None: ink
    /// has none for it (a product, a bare number, a sum with a bare term, a
    /// temperature that was worked on).
    pub unit: Option<&'static str>,
    pub file: Rc<str>,
    pub line: u32,
}

/// A part's formula for a fact (Doc::formula): the facts it reads, each with the variable that stands for it and
/// the unit it is read in ("": a pure number); its value and unit; and where it holds, with why not elsewhere.
#[derive(Clone, Debug)]
pub struct Formula {
    pub what: String,
    pub reads: Vec<(String, usize, &'static str)>,
    pub value: NodeId,
    pub unit: Option<&'static str>,
    pub only: Option<(NodeId, String)>,
    pub file: Rc<str>,
    pub line: u32,
}

/// `(require NAME ASKED)`: a question the design asks. `(within-ratings
/// PATH)`: that every rating of every part at or under PATH holds
/// (`within`). A claim, with its clauses: that claim, which is the
/// `AssertDef` of this name (`within` is None). Data: nothing
/// answers it at load.
#[derive(Clone, Debug)]
pub struct Require {
    pub name: String,
    pub within: Option<String>,
    /// `(within-ratings PATH (excluding P ...))`: the parts at or under
    /// PATH the design says, on purpose, the question is not about.
    pub excluding: Vec<String>,
    pub file: Rc<str>,
    pub line: u32,
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
    /// Scalar arguments passed to direct instances, before a state load
    /// replaces them with independent coordinates.
    coordinate_inputs: Vec<NodeId>,
    /// Declared frames, in the assembly's coordinates.
    pub frames: Vec<(String, Xf)>,
    /// The instances (paths in it) each frame is of (Instance::on); and,
    /// for a frame written as exactly the language's `(frame-of PATH ITS)`,
    /// which of that instance's child's frames ITS is (its place in the
    /// child's `frames`): the two are one place.
    pub frame_on: Vec<(String, Vec<String>, Option<u32>)>,
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
    /// The run of that file in which the claim or clause is stated.
    pub scope: Scope,
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
    passed: Option<(Rc<str>, u32, Scope)>,
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
    /// Emptiness of a region enclosed between an inner and an outer solid.
    EmptyBounds { certain: NodeId, possible: NodeId },
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
            Claim::Empty(_) | Claim::EmptyBounds { .. } | Claim::Disjoint(..) | Claim::Subset(..) | Claim::Clearance(..) | Claim::Overlaps(..) | Claim::Touches(..) => {
                Accuracy::Length
            }
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
            Claim::EmptyBounds { certain, possible } => vec![*certain, *possible],
            Claim::Clearance(p, q, _) | Claim::Disjoint(p, q) | Claim::Subset(p, q) | Claim::Overlaps(p, q) | Claim::Touches(p, q) => vec![*p, *q],
            Claim::Point { f, at, .. } => std::iter::once(*f).chain(at.iter().copied()).collect(),
            Claim::Integral { parts, frames, .. } => parts.iter().flat_map(|(s, gs)| std::iter::once(*s).chain(gs.iter().copied())).chain(framed(frames)).collect(),
            Claim::Mean { parts, frames, .. } => parts.iter().flat_map(|(s, w, hs)| [*s, *w].into_iter().chain(hs.iter().copied())).chain(framed(frames)).collect(),
            Claim::Metric { metric, .. } => metric.shape.into_iter().collect(),
        }
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
    /// The run of that file in which the claim or clause is stated.
    pub scope: Scope,
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
    /// Text params, `(param name "Ada" ...)`: load-time strings, each
    /// bound as a def of its value (TextParam).
    pub texts: Vec<TextParam>,
    /// At a load given declared cases (Given::case: a combination ink
    /// assert asks), a text param's refusal of the text it took there: the
    /// param, the text and why. The load stopped there (run_forms), so the
    /// document holds only what came before: no combination to ask.
    pub refused: Vec<(TextParam, String, Refusal)>,
    /// Each param that names cases, `(cases v ...)`, and their values:
    /// what `ink require` checks its claims at besides the param's value.
    pub cases: Vec<(String, Vec<Given>)>,
    /// The fonts read and glyphs built (font.rs).
    fonts: font::Fonts,
    pub defs: Vec<(String, Value)>,
    pub fns: HashMap<String, Rc<FnDef>>,
    pub bodies: Vec<Body>,
    pub assemblies: Vec<Assembly>,
    /// Bodies and assemblies with arguments, by name: what an instance
    /// passing arguments builds again.
    pub arg_defs: HashMap<String, Rc<ArgDef>>,
    /// Instantiations by definition and argument values. An assembly also
    /// owns claims in the scope that passed those arguments.
    built: HashMap<(String, String, Option<Scope>), String>,
    /// The values of the arguments of the body or assembly whose clauses
    /// are being read: what they read the arguments' names as.
    args: Env,
    pub asserts: Vec<AssertDef>,
    /// The questions the file asks with `require`, in order.
    pub requires: Vec<Require>,
    // For a generated placement claim, where its underlying limit was declared.
    claim_declarations: HashMap<String, String>,
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
    /// The run file's `(t-span A B)`: the motion its claims that read t
    /// are checked over by default.
    pub t_span: Option<(f64, f64)>,
    /// World fields of placed bodies (an instance, a construction in one,
    /// `place`), with the body's own field and the placement: a body that
    /// moves rigidly, which motion checks bound by how far its points move.
    /// Keyed by field, so one of the decompositions of equal fields: an
    /// item's own body and placement are its instance's and links'.
    pub placed: HashMap<NodeId, Placed>,
    /// Each placed body's world field, by the path of the instance placing it (an Item's label): a
    /// failure's witness names the bodies there (claims::bodies_at).
    pub placed_labels: HashMap<NodeId, String>,
    /// A file runs once in a namespace. Repeated imports check the directory
    /// its relative paths resolve from and the values of any supplied params.
    imported: HashMap<(PathBuf, Option<String>), (usize, String)>,
    std_imported: bool,
    /// Each param an import set (ImportValue), by index: its written
    /// default and where the import is ("b.ink:4").
    pub set_at: HashMap<usize, (f64, String)>,
    /// Only syntax is shared across file runs, by exact source text.
    /// Function bodies retain their shared subtrees after the load.
    parsed: HashMap<String, std::sync::Arc<[Sexp]>>,
    /// The JSON files read (`read-json`), by path: read once, shared.
    json: HashMap<std::path::PathBuf, Data>,
    /// The fingerprint of each file's content the load read (`Doc::read`):
    /// see `sources`.
    sources_read: Vec<u128>,
    /// The files on disk the load depended on: see `files_read`.
    files_read: Vec<(std::path::PathBuf, Option<u128>)>,
    /// Params a solve sets (by index): numerical at their values (a
    /// numerical method's answer), certified where a view derives them over
    /// ranged params (ModelView::set_param_spans).
    pub solved_params: HashSet<usize>,
    /// Sampled data found at load time from other sampled data (a box of
    /// it: measured_box), each with what it was found from, which an
    /// answer reading it reads too (Doc::inexact).
    found_from: Vec<(std::sync::Arc<crate::grid::Source>, Vec<std::sync::Arc<crate::grid::Source>>)>,
    /// The variables `(ranged lo hi)` made (params at lo, no tolerance),
    /// each with [lo, hi]: a range of the question, which every claim
    /// reading one ranges over and no run pins (ModelView::declared_spans).
    pub ranged: Vec<(usize, (f64, f64))>,
    /// Params a construction sets at 0 that stand for how far the shape it
    /// built may be from the one written (polygon's when it is not built
    /// exactly), each with that distance (`allowed`): a band on what a claim
    /// reads, as a solve's ParamProvider::error is.
    pub allowances: Vec<(usize, f64)>,
    /// Free depths made by overhangs. A claim of emptiness reading one
    /// must exclude the defect outright, however thin the roof is.
    support_depths: HashSet<u32>,
    /// Witness and coverage fields for a complete, unclipped overhang mask.
    print_bounds: HashMap<NodeId, support::Overhang>,
    /// Each value claim's bound as written, by its value's node and the bound's float: the node whose real
    /// value it is (Graph::real_value), which a float alone cannot find (an inexact constant has its own).
    pub value_bounds: HashMap<(NodeId, u64), NodeId>,
    /// Each number param's values written in digits (its default, its cases), by their float: the exact
    /// fraction each writes (Doc::exact_value).
    param_exact: HashMap<usize, Vec<(u64, ink_core::exact::Ratio)>>,
    /// Each solved param's rate along each input asked for (Doc::rate):
    /// the param its rates' solve sets, or None where it does not vary.
    /// Params a derivative held fixed along an input (no solve set them
    /// when it was taken), with where: a solve read later may not make one
    /// vary along it (Doc::varies_held).
    held: HashMap<(usize, Var), String>,
    /// Joints' speed-limit claims (index in `asserts`), each with its
    /// drive and limits: their rates are taken once the file is read
    /// (Doc::rate_joints), through every solve their drives read, however
    /// the file orders them (a loop closed over the assembly's frames
    /// comes after it).
    rated: Vec<(usize, NodeId, NodeId, NodeId, bool)>,
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
    magnitudes: Vec<magnitude::Guard>,
    /// A claim command installs its scope after loading. Its midpoint
    /// cannot stand in for the range when validating material fields.
    claim_materials: bool,
    /// The bound of the claim being parsed (see `claim`).
    bound: Option<NodeId>,
    /// Warnings found while loading (undefined cases in constants).
    pub load_warnings: Vec<String>,
    // Optional limit declarations by assembly, expanded to paths only for the placed root.
    dropped_limits: HashMap<String, Vec<(String, UncheckedClaim)>>,
    /// What each param's quantity was said by, by index, where its value
    /// was (units::Said): a slot taking a known quantity checks it.
    param_said: HashMap<usize, Rc<units::Said>>,
    param_requirements: HashMap<String, units::UnitRequirement>,
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
    /// Each field a `warp` reads elsewhere, with its copy (magnitude.rs: a copy is checked as its original).
    warped: Vec<(NodeId, NodeId)>,
    /// Each field a `warp` of x, y and z alone made (not of t), with the field it reads and where: the
    /// warped field is that field at the three coordinates (carried::rigidly reads a move or a turn
    /// back as a placement).
    pub(crate) warps: HashMap<NodeId, (NodeId, [NodeId; 3])>,
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
    /// Where each public value or function name was last defined.
    origin: HashMap<String, (u32, u32)>,
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
    /// Each ranged param a load-time box or reach was found over the
    /// whole band of (`load_scope`), with the band and where: a run may
    /// range one within the band, not beyond it (ModelView::scope_spans).
    covered: Vec<(usize, (f64, f64), String)>,
    /// The span of t each load-time box or reach of a field reading t
    /// holds over (`load_scope`), and where: a moment or motion beyond one
    /// is refused (Doc::covers_t).
    load_t: Vec<((f64, f64), String)>,
    /// Undefined cases the current form's load-time numbers passed through
    /// (`load_eval`), for `at_form` to report.
    undefined_at_load: std::cell::RefCell<Vec<&'static str>>,
    /// `(shape-of path)`'s fields by path, while no top-level form but
    /// claims has run since: what is defined does not change between them,
    /// and a whole assembly's is a walk of every item (an assembly's claims
    /// read it again and again: the RC car's 100 tool-access claims, 3 s of
    /// its 3.4 s load).
    shape_of: HashMap<(Scope, String), NodeId>,
    /// The top-level form being run: file:line (head name).
    form_at: String,
    /// The first number written outside ink's range and where (file:line):
    /// what the range's unknown cites (validity::numbers).
    pub(crate) wild_written: Option<(f64, String)>,
    /// Its name (the def's, the claim's), which a `ranged` it makes is named by.
    form_name: String,
    /// The assembly being read, if any.
    building: Option<Building>,
    /// The assembly whose local instance names are in lexical scope. A
    /// top-level function has none, even when called to place an instance.
    local_assembly: Option<Rc<str>>,
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
    /// The member the next file run is (family::member_of): an import's,
    /// or the file loaded's.
    member: Option<family::Member>,
}

/// The source a document was loaded from (Doc::load_named's inputs).
#[derive(Clone)]
struct Origin {
    src: String,
    name: String,
    path: std::path::PathBuf,
    base: std::path::PathBuf,
    /// The member `path` names, its family's file run as it (family::member_of).
    member: Option<family::Member>,
    /// The directories its imports were read from before `$INK_PATH`'s
    /// (Doc::load_file_importing).
    import_dirs: Vec<std::path::PathBuf>,
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
        Doc::load_origin(self.registry, self.origin.clone(), &self.overrides, Some((target.to_string(), false)), false).map_err(|(e, _)| e)
    }

    /// Doc::with_states.
    pub fn with_states(&self, target: &str) -> R<Doc> {
        Doc::load_origin(self.registry, self.origin.clone(), &self.overrides, Some((target.to_string(), true)), false).map_err(|(e, _)| e)
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

/// A direction that is not a constant (Doc::unit): its squared length
/// or, for a single axial component, its absolute value,
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

#[cfg(test)]
mod tests;
