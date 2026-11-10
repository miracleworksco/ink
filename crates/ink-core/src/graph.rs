//! Hash-consed expression DAG with constant folding.
//!
//! Node ids are topologically ordered: every node's children have smaller ids.
//! Equal math is stored once.

// These maps index graph-owned keys; node ids follow insertion order,
// and traversals sort their ids before rebuilding. Hashing changes neither.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

pub type NodeId = u32;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Var {
    X,
    Y,
    Z,
    T,
    Param(u32),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum UnOp {
    Neg,
    Abs,
    Sqrt,
    Square,
    Exp,
    Ln,
    Sin,
    Cos,
    Tan,
    Floor,
    /// Degrees in radians: a pi / 180, pi exact. A turn's angle is read so
    /// (its sine and cosine), and a whole turn of degrees is then exactly
    /// 2 pi: the turns either side of a wrap at 360 agree exactly
    /// (eval::scope_run's jumps).
    Rad,
    /// Radians in degrees: a 180 / pi, pi exact (atan2's branches, a whole
    /// turn apart, are exactly 360 apart).
    Deg,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BinOp {
    Add,
    Sub,
    Mul,
    Div,
    Min,
    Max,
    Atan2,
    Pow,
    /// Euclidean max: sqrt(a^2 + b^2) where both are positive, else max(a, b).
    /// The exact distance to the intersection of two exact fields with
    /// orthogonal gradients (box, cylinder, extrusion). Gradient bounds know
    /// its branches never add up, so it stays 1-Lipschitz where the sum
    /// formula min(max(a, b), 0) + len(max(a, 0), max(b, 0)) looks 2.
    Emax,
}

impl UnOp {
    /// Every operation, each at its `index` (a new one must be listed).
    pub const ALL: [UnOp; 12] = [UnOp::Neg, UnOp::Abs, UnOp::Sqrt, UnOp::Square, UnOp::Exp, UnOp::Ln, UnOp::Sin, UnOp::Cos, UnOp::Tan, UnOp::Floor, UnOp::Rad, UnOp::Deg];

    pub fn index(self) -> usize {
        match self {
            UnOp::Neg => 0,
            UnOp::Abs => 1,
            UnOp::Sqrt => 2,
            UnOp::Square => 3,
            UnOp::Exp => 4,
            UnOp::Ln => 5,
            UnOp::Sin => 6,
            UnOp::Cos => 7,
            UnOp::Tan => 8,
            UnOp::Floor => 9,
            UnOp::Rad => 10,
            UnOp::Deg => 11,
        }
    }
}

impl BinOp {
    /// Every operation, each at its `index` (a new one must be listed).
    pub const ALL: [BinOp; 9] = [BinOp::Add, BinOp::Sub, BinOp::Mul, BinOp::Div, BinOp::Min, BinOp::Max, BinOp::Atan2, BinOp::Pow, BinOp::Emax];

    pub fn index(self) -> usize {
        match self {
            BinOp::Add => 0,
            BinOp::Sub => 1,
            BinOp::Mul => 2,
            BinOp::Div => 3,
            BinOp::Min => 4,
            BinOp::Max => 5,
            BinOp::Atan2 => 6,
            BinOp::Pow => 7,
            BinOp::Emax => 8,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum Node {
    Const(f64),
    Var(Var),
    Un(UnOp, NodeId),
    Bin(BinOp, NodeId, NodeId),
    /// `Select(c, a, b)`: `a` where `c < 0`, else `b`.
    Select(NodeId, NodeId, NodeId),
    /// Identity carrying a name (index into `Graph::tag_names`): marks where a
    /// named construction enters the math, for provenance (named faces).
    Tag(u32, NodeId),
    /// A sampled field (index into `Graph::grids`): component `mask` of its
    /// trilinear interpolant at coordinates (x, y, z), clamped into the grid
    /// (see grid::GridData::component: the value for mask 0, else the
    /// derivative along the mask's axes). A read of the field (Node::Read)
    /// is it plus the outside continuation, in plain ops.
    Grid(u32, u8, [NodeId; 3]),
    /// `Read(g, path, a, at)`: sampled field g (index into `Graph::grids`)
    /// read at the coordinates of the Frame `at`, differentiated along
    /// `path` (index into `Graph::paths`: the grid's axes, in the order
    /// differentiated; none for the value). In value it is `a`: the same
    /// derivative of the samples' interpolant and continuation there
    /// (`Graph::read` makes it from the coordinates).
    Read(u32, u32, NodeId, NodeId),
    /// `Lipschitz(a, k, frame)`: `a` in value, plus the promise that `a` is a
    /// k-Lipschitz function of the frame's coordinates (exact distance
    /// fields: k = 1). Warps rewrite the frame like any expression, so the
    /// promise follows the shape; gradient bounds use k times the frame's
    /// stretch.
    Lipschitz(NodeId, NodeId, NodeId),
    /// The coordinates (u, v, w) of a Lipschitz promise. Value 0.
    Frame([NodeId; 3]),
    /// `Deriv(a, f, axis)`: `a` in value, which is f's derivative along x,
    /// y or z (axis 0, 1, 2) as `Graph::diff` writes it. Ranges then also
    /// hold it within f's gradient bounds (interval::G's partials), which
    /// keep relations the formula's own ranges lose: the slope of
    /// sqrt(u1^2 + u2^2 + ...) along x is at most sqrt(sum (du_i/dx)^2),
    /// across the distance's skeleton too, where the formula is a quotient
    /// of two ranges through 0. Along t (axis 3) it only marks a as read
    /// from f's derivative, which the formula drops where f jumps: a claim
    /// on it over a motion is unknown where f may jump (ink-lang's
    /// motion::ranged). A warp keeps the mark on the warped f.
    Deriv(NodeId, NodeId, u8),
    /// `Real(a)`: `a`, a constants-and-params value, read as a real number:
    /// what reads it is arithmetic in real numbers, bounded outward, not a
    /// constants-and-params value taken in floats (Tape::fixed). A rewrite
    /// puts it where spatial arithmetic's operands became constants and
    /// params (Graph::rebuild): `x - y` along x = y is 0, and 1e16 plus it
    /// plus 0.5 still means 1e16 + 0.5, not its float.
    Real(NodeId),
    /// `Within(lo, hi, v)`: a number known only to lie within the constants
    /// lo..hi (a circle's radius, say, made from its data with a certified
    /// bound on the rounding): ranges bound what reads it over all of
    /// lo..hi, so an answer holds for whichever value it is, and a point's
    /// value is the constant v, a float within them. Never fixed. The u64
    /// is its identity (Graph::within): one number for one identity.
    Within(NodeId, NodeId, NodeId, u64),
}

impl Node {
    pub fn children(&self) -> ([NodeId; 3], usize) {
        match *self {
            Node::Const(_) | Node::Var(_) => ([0; 3], 0),
            Node::Un(_, a) => ([a, 0, 0], 1),
            Node::Bin(_, a, b) => ([a, b, 0], 2),
            Node::Select(c, a, b) | Node::Lipschitz(c, a, b) | Node::Within(c, a, b, _) => ([c, a, b], 3),
            Node::Tag(_, a) | Node::Real(a) => ([a, 0, 0], 1),
            Node::Read(_, _, a, at) => ([a, at, 0], 2),
            Node::Deriv(a, f, _) => ([a, f, 0], 2),
            Node::Grid(_, _, c) | Node::Frame(c) => (c, 3),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Key {
    Const(u64),
    Var(Var),
    Un(UnOp, NodeId),
    Bin(BinOp, NodeId, NodeId),
    Select(NodeId, NodeId, NodeId),
    Tag(u32, NodeId),
    Grid(u32, u8, [NodeId; 3]),
    Read(u32, u32, NodeId, NodeId),
    Lipschitz(NodeId, NodeId, NodeId),
    Frame([NodeId; 3]),
    Deriv(NodeId, NodeId, u8),
    Real(NodeId),
    /// By its identity alone (see Graph::within).
    Within(u64),
}

#[derive(Default)]
pub struct Graph {
    pub nodes: Vec<Node>,
    index: HashMap<Key, NodeId>,
    /// Per node, a hash of its expression (not of where in the graph it
    /// is, nor of the names in it), which orders the operands of
    /// commutative operations: so an expression has one form however the
    /// graph holding it was built (Tape::new), and its tape one
    /// fingerprint.
    shapes: Vec<u128>,
    /// Per node: its value is made of constants and params alone, as
    /// Tape::fixed (an operation, select or name of those), so it is taken
    /// in floats.
    fixed: Vec<bool>,
    pub tag_names: Vec<String>,
    /// Each param's name, by number (the model names them: `name_param`):
    /// what orders operands alike but for the params they read
    /// (Tape::new), which their numbers, where the model declared them,
    /// cannot.
    pub param_names: Vec<String>,
    /// Names provenance stops at (a fillet: its surface is the fillet's,
    /// not that of the faces it blends).
    pub opaque: Vec<u32>,
    pub grids: Vec<std::sync::Arc<crate::grid::GridData>>,
    /// The derivatives sampled fields are read along (Node::Read): grid
    /// axes (0 x, 1 y, 2 z), in the order differentiated.
    pub paths: Vec<Vec<u8>>,
    path_index: HashMap<Vec<u8>, u32>,
    /// Per grid and path: the read at the coordinates (x, y, z), which
    /// `Graph::read` puts its coordinates into.
    read_forms: HashMap<(u32, u32), NodeId>,
    /// The first number outside ink's range (docs/CORE.md Answers) the
    /// model is built from (`note`), its bits (0: none).
    wild: std::sync::atomic::AtomicU64,
    /// std's own constants outside the range (`own`), and folds of them.
    owned: HashSet<NodeId>,
}

/// Whether `v` is outside ink's range: finite, not 0, and smaller than
/// 1e-100 or larger than 1e100 in size.
pub fn wild(v: f64) -> bool {
    v != 0.0 && v.is_finite() && !(1e-100..=1e100).contains(&v.abs())
}

/// Whether floats lost `r`, computed from `args` (by an operation that
/// `scales`: a product, a quotient, a power, exp): a subnormal, 0 from
/// finite numbers that are not, or an infinity from them (an overflow; ln 0
/// is exactly -inf).
pub fn lost(r: f64, args: &[f64], scales: bool) -> bool {
    let real = args.iter().all(|&a| a != 0.0 && a.is_finite());
    (r != 0.0 && r.abs() < f64::MIN_POSITIVE) || (real && (r.is_infinite() || (scales && r == 0.0)))
}

impl Graph {
    pub fn new() -> Self {
        Self::default()
    }

    fn intern(&mut self, n: Node) -> NodeId {
        let key = match n {
            Node::Const(v) => Key::Const(if v == 0.0 { 0 } else { v.to_bits() }),
            Node::Var(v) => Key::Var(v),
            Node::Un(o, a) => Key::Un(o, a),
            Node::Bin(o, a, b) => Key::Bin(o, a, b),
            Node::Select(c, a, b) => Key::Select(c, a, b),
            Node::Tag(t, a) => Key::Tag(t, a),
            Node::Grid(g, m, c) => Key::Grid(g, m, c),
            Node::Read(g, d, a, at) => Key::Read(g, d, a, at),
            Node::Lipschitz(a, k, f) => Key::Lipschitz(a, k, f),
            Node::Frame(c) => Key::Frame(c),
            Node::Deriv(a, f, k) => Key::Deriv(a, f, k),
            Node::Real(a) => Key::Real(a),
            Node::Within(.., id) => Key::Within(id),
        };
        if let Some(&id) = self.index.get(&key) {
            return id;
        }
        let id = self.nodes.len() as NodeId;
        let shape = self.shape_of(n);
        let fixed = match n {
            Node::Const(_) | Node::Var(Var::Param(_)) => true,
            Node::Un(..) | Node::Bin(..) | Node::Select(..) | Node::Tag(..) => {
                let (ch, k) = n.children();
                ch[..k].iter().all(|&c| self.fixed[c as usize])
            }
            _ => false,
        };
        self.nodes.push(n);
        self.shapes.push(shape);
        self.fixed.push(fixed);
        self.index.insert(key, id);
        id
    }

    /// The hash of `n`'s expression (`shapes`).
    fn shape_of(&self, n: Node) -> u128 {
        let h = |a: NodeId| self.shapes[a as usize];
        let mut f = Fingerprint::default();
        match n {
            Node::Const(v) => f.word(1).f64(v),
            // Not its number: a param's number is where the model declared
            // or imported it, so an expression's form (and its tape's
            // fingerprint) is the same in any model that imports it.
            Node::Var(Var::Param(_)) => f.word(2),
            Node::Var(v) => f.word(3).word(match v {
                Var::X => 0,
                Var::Y => 1,
                Var::Z => 2,
                _ => 3,
            }),
            Node::Un(op, a) => f.word(4).word(op as u64).key(h(a)),
            Node::Bin(op, a, b) => f.word(5).word(op as u64).key(h(a)).key(h(b)),
            Node::Select(c, a, b) => f.word(6).key(h(c)).key(h(a)).key(h(b)),
            // Names are transparent: renaming a part reorders nothing.
            Node::Tag(_, a) => return h(a),
            Node::Grid(g, m, c) => {
                self.grids[g as usize].fingerprint(f.word(8).word(m as u64));
                c.iter().for_each(|n| _ = f.key(h(*n)));
                &mut f
            }
            Node::Read(g, path, a, at) => {
                self.grids[g as usize].fingerprint(f.word(11));
                self.paths[path as usize].iter().for_each(|&k| _ = f.word(k as u64));
                f.key(h(a)).key(h(at))
            }
            Node::Lipschitz(a, k, fr) => f.word(9).key(h(a)).key(h(k)).key(h(fr)),
            Node::Frame(c) => f.word(10).key(h(c[0])).key(h(c[1])).key(h(c[2])),
            Node::Deriv(a, fd, k) => f.word(12).word(k as u64).key(h(a)).key(h(fd)),
            Node::Real(a) => f.word(13).key(h(a)),
            Node::Within(l, hi, v, id) => f.word(14).word(id).key(h(l)).key(h(hi)).key(h(v)),
        };
        f.value()
    }

    /// A constant: a number the model is built from (written, read from a
    /// file, made in Rust), noted where outside the range.
    pub fn constant(&mut self, v: f64) -> NodeId {
        self.note(v, wild(v));
        self.number(v)
    }

    fn number(&mut self, v: f64) -> NodeId {
        self.intern(Node::Const(if v == 0.0 { 0.0 } else { v }))
    }

    /// Note `v` as the model's first number outside the range, if `outside`
    /// (a 0 floats lost as -0, whose bits are not 0's).
    pub fn note(&self, v: f64, outside: bool) {
        use std::sync::atomic::Ordering::Relaxed;
        if outside {
            let _ = self.wild.compare_exchange(0, if v == 0.0 { (-0.0f64).to_bits() } else { v.to_bits() }, Relaxed, Relaxed);
        }
    }

    /// The first number noted outside the range.
    pub fn wild(&self) -> Option<f64> {
        let b = self.wild.load(std::sync::atomic::Ordering::Relaxed);
        (b != 0).then(|| f64::from_bits(b))
    }

    /// `constant` for std's own implementation (norm2-up's 2^1000): not
    /// the model's number, so not noted.
    pub fn own(&mut self, v: f64) -> NodeId {
        let id = self.number(v);
        if wild(v) {
            self.owned.insert(id);
        }
        id
    }

    /// A fold `v` of the constants `args` (`vals`): noted where floats lost
    /// it (`lost`), unless all of them are std's own.
    fn folded(&mut self, v: f64, args: &[NodeId], vals: &[f64], scales: bool) -> NodeId {
        if args.iter().all(|a| self.owned.contains(a)) {
            return self.own(v);
        }
        self.note(v, lost(v, vals, scales));
        self.number(v)
    }

    pub fn var(&mut self, v: Var) -> NodeId {
        self.intern(Node::Var(v))
    }

    /// `a` read as a real number (Node::Real), where it is a constants-and-
    /// params value; anything else already is one.
    pub fn real(&mut self, a: NodeId) -> NodeId {
        if self.fixed[a as usize] { self.intern(Node::Real(a)) } else { a }
    }

    /// A number known only within lo..hi, its point value v (Node::Within);
    /// v itself where they are equal. `id` is its identity: a hash of what
    /// defines the exact number (its construction's exact inputs and which
    /// of its numbers it is), never of its bounds (two exact numbers can
    /// round to the same ones) nor of when it was made. One identity is one
    /// number: one node, so x - x folds to 0, c c is a square, and keys
    /// and structural comparisons read the identity; two identities are
    /// two numbers, whatever their bounds.
    pub fn within(&mut self, lo: f64, hi: f64, v: f64, id: u64) -> NodeId {
        assert!(lo <= v && v <= hi, "within: {lo} <= {v} <= {hi}");
        let (l, h, c) = (self.constant(lo), self.constant(hi), self.constant(v));
        if lo == hi { c } else { self.intern(Node::Within(l, h, c, id)) }
    }

    /// An operation that read a spatial value, now on constants and params
    /// alone (a rewrite's: `Graph::rebuild`): folded where its float is its
    /// real value (`exact_un`), else kept, on its operand read as a real
    /// number, so its range is rounded outward like any spatial value's.
    fn un_real(&mut self, op: UnOp, a: NodeId) -> NodeId {
        if let Some(x) = self.as_const(a)
            && let Some(v) = exact_un(op, x)
        {
            return self.folded(v, &[a], &[x], matches!(op, UnOp::Square | UnOp::Exp));
        }
        let a = self.real(a);
        self.un(op, a)
    }

    /// `un_real` for two operands (`exact_bin`).
    fn bin_real(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        if let (Some(x), Some(y)) = (self.as_const(a), self.as_const(b))
            && let Some(v) = exact_bin(op, x, y)
        {
            return self.folded(v, &[a, b], &[x, y], matches!(op, BinOp::Mul | BinOp::Div | BinOp::Pow | BinOp::Emax));
        }
        let a = self.real(a);
        self.bin(op, a, b)
    }

    pub fn as_const(&self, id: NodeId) -> Option<f64> {
        match self.nodes[id as usize] {
            Node::Const(v) => Some(v),
            _ => None,
        }
    }

    /// Building, not rewriting: a field rewritten (variables put in,
    /// nodes replaced, names removed) is remade by `rebuild`, whose rule
    /// keeps spatial arithmetic real where its operands became constants.
    ///
    /// Constants fold, except into an undefined case: that stays an
    /// operation, so it is reported where a value depends on it (an
    /// evaluation's active path, `eval::undefined_cases`) like any other,
    /// and not where the value is discarded (multiplied by a constant 0).
    pub fn un(&mut self, op: UnOp, a: NodeId) -> NodeId {
        if let Some(v) = self.as_const(a)
            && undefined_un(op, v).is_none()
        {
            return self.folded(un_f64(op, v), &[a], &[v], matches!(op, UnOp::Square | UnOp::Exp));
        }
        match (op, self.nodes[a as usize]) {
            (UnOp::Neg, Node::Un(UnOp::Neg, inner)) => return inner,
            (UnOp::Abs, Node::Un(UnOp::Abs, _)) => return a,
            _ => {}
        }
        self.intern(Node::Un(op, a))
    }

    pub fn bin(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        use BinOp::*;
        let (ca, cb) = (self.as_const(a), self.as_const(b));
        if let (Some(x), Some(y)) = (ca, cb) {
            return match undefined_bin(op, x, y) {
                None => self.folded(bin_f64(op, x, y), &[a, b], &[x, y], matches!(op, Mul | Div | Pow | Emax)),
                Some(_) => self.intern(Node::Bin(op, a, b)),
            };
        }
        match op {
            Add => {
                if ca == Some(0.0) {
                    return b;
                }
                if cb == Some(0.0) {
                    return a;
                }
            }
            Sub => {
                if cb == Some(0.0) {
                    return a;
                }
                if ca == Some(0.0) {
                    return self.un(UnOp::Neg, b);
                }
                if a == b {
                    return self.constant(0.0);
                }
            }
            Mul => {
                if ca == Some(1.0) {
                    return b;
                }
                if cb == Some(1.0) {
                    return a;
                }
                if ca == Some(0.0) || cb == Some(0.0) {
                    return self.constant(0.0);
                }
                if ca == Some(-1.0) {
                    return self.un(UnOp::Neg, b);
                }
                if cb == Some(-1.0) {
                    return self.un(UnOp::Neg, a);
                }
            }
            Div => {
                if cb == Some(1.0) {
                    return a;
                }
            }
            Min | Max => {
                if a == b {
                    return a;
                }
            }
            Pow => {
                if cb == Some(1.0) {
                    return a;
                }
                if cb == Some(2.0) {
                    return self.un(UnOp::Square, a);
                }
            }
            Atan2 | Emax => {}
        }
        // Commutative: operands in the order of their expressions' hashes
        // (which read no param's number), then of their nodes (`tied`).
        let order = |n: NodeId| (self.shapes[n as usize], n);
        let (a, b) = if matches!(op, Add | Mul | Min | Max | Emax) && order(b) < order(a) { (b, a) } else { (a, b) };
        self.intern(Node::Bin(op, a, b))
    }

    /// Whether `n` is a commutative operation whose operands' expressions
    /// are alike up to the names in them (params, tags), so ordered by
    /// node: Tape::new orders them by those names in the question.
    pub(crate) fn tied(&self, n: NodeId) -> bool {
        matches!(self.nodes[n as usize], Node::Bin(BinOp::Add | BinOp::Mul | BinOp::Min | BinOp::Max | BinOp::Emax, a, b) if a != b && self.shapes[a as usize] == self.shapes[b as usize])
    }

    pub fn select(&mut self, c: NodeId, a: NodeId, b: NodeId) -> NodeId {
        if let Some(v) = self.as_const(c) {
            return if v < 0.0 { a } else { b };
        }
        if a == b {
            return a;
        }
        self.intern(Node::Select(c, a, b))
    }

    /// Name a subexpression (identity in value).
    pub fn tag(&mut self, name: &str, a: NodeId) -> NodeId {
        let t = match self.tag_names.iter().position(|n| n == name) {
            Some(i) => i as u32,
            None => {
                self.tag_names.push(name.to_string());
                (self.tag_names.len() - 1) as u32
            }
        };
        self.intern(Node::Tag(t, a))
    }

    /// The hash of n's expression (`shapes`).
    pub(crate) fn shape(&self, n: NodeId) -> u128 {
        self.shapes[n as usize]
    }

    /// Name param k (`param_names`).
    pub fn name_param(&mut self, k: u32, name: &str) {
        let k = k as usize;
        if self.param_names.len() <= k {
            self.param_names.resize(k + 1, String::new());
        }
        self.param_names[k] = name.to_string();
    }

    /// A name that provenance stops at: points under it report this name
    /// only, not the names inside.
    pub fn tag_opaque(&mut self, name: &str, a: NodeId) -> NodeId {
        let n = self.tag(name, a);
        if let Node::Tag(t, _) = self.nodes[n as usize]
            && !self.opaque.contains(&t)
        {
            self.opaque.push(t);
        }
        n
    }

    /// A sampled field read at the given coordinates (Node::Read): in value,
    /// its trilinear interpolant at the coordinates clamped into the grid,
    /// plus the distance to the grid box outside it unless the grid holds
    /// its boundary value there (`GridData::holds`) (in plain ops, so the
    /// core rules give its derivatives and ranges; the distance's squares
    /// are max(x - hi, 0)^2 + max(lo - x, 0)^2, exactly 0 inside).
    pub fn grid(&mut self, data: std::sync::Arc<crate::grid::GridData>, at: [NodeId; 3]) -> NodeId {
        let id = match self.grids.iter().position(|g| std::sync::Arc::ptr_eq(g, &data)) {
            Some(i) => i,
            None => {
                self.grids.push(data.clone());
                self.grids.len() - 1
            }
        } as u32;
        let value = self.path(Vec::new());
        self.read(id, value, at)
    }

    /// The index of a derivative path (`Graph::paths`).
    fn path(&mut self, p: Vec<u8>) -> u32 {
        if let Some(&i) = self.path_index.get(&p) {
            return i;
        }
        self.paths.push(p.clone());
        self.path_index.insert(p, (self.paths.len() - 1) as u32);
        (self.paths.len() - 1) as u32
    }

    /// Grid g read at `at` along derivative `path` (Node::Read): its value
    /// is made from the coordinates alone, so equal reads are equal nodes.
    fn read(&mut self, g: u32, path: u32, at: [NodeId; 3]) -> NodeId {
        let form = self.read_form(g, path);
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let map: HashMap<NodeId, NodeId> = (0..3).map(|k| (xyz[k], at[k])).collect();
        let a = self.rebuild(form, &map);
        let frame = self.intern(Node::Frame(at));
        self.intern(Node::Read(g, path, a, frame))
    }

    /// The form of grid g's read along `path` (read_form), made with every
    /// read of it (Graph::read).
    pub fn read_form_of(&self, g: u32, path: u32) -> NodeId {
        self.read_forms[&(g, path)]
    }

    /// The samples' read of grid g along `path` at (x, y, z): the value
    /// (interpolant at the clamped point, plus the distance to the grid
    /// box unless the grid holds), and each further axis its derivative
    /// along it.
    fn read_form(&mut self, g: u32, path: u32) -> NodeId {
        if let Some(&f) = self.read_forms.get(&(g, path)) {
            return f;
        }
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let f = match self.paths[path as usize].clone().split_last() {
            Some((&axis, before)) => {
                let before = self.path(before.to_vec());
                let f = self.read_form(g, before);
                self.diff(f, [Var::X, Var::Y, Var::Z][axis as usize])
            }
            None => {
                let (lo, hi) = self.grids[g as usize].corners();
                let holds = self.grids[g as usize].holds();
                let zero = self.constant(0.0);
                let mut q = xyz;
                let mut out = zero;
                for k in 0..3 {
                    let (l, h) = (self.constant(lo[k]), self.constant(hi[k]));
                    let below_hi = self.min(xyz[k], h);
                    q[k] = self.max(below_hi, l);
                    if holds {
                        continue;
                    }
                    let up = self.sub(xyz[k], h);
                    let up = self.max(up, zero);
                    let dn = self.sub(l, xyz[k]);
                    let dn = self.max(dn, zero);
                    let (up2, dn2) = (self.un(UnOp::Square, up), self.un(UnOp::Square, dn));
                    let both = self.add(up2, dn2);
                    out = self.add(out, both);
                }
                let raw = self.intern(Node::Grid(g, 0, q));
                if holds {
                    raw
                } else {
                    let dist = self.un(UnOp::Sqrt, out);
                    self.add(raw, dist)
                }
            }
        };
        self.read_forms.insert((g, path), f);
        f
    }

    /// Promise that `a` is k-Lipschitz in x, y, z (identity in value). A k
    /// that is not positive promises nothing.
    pub fn lipschitz(&mut self, a: NodeId, k: NodeId) -> NodeId {
        if self.as_const(a).is_some() || self.as_const(k).is_some_and(|k| !(k > 0.0)) {
            return a;
        }
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let frame = self.intern(Node::Frame(xyz));
        self.intern(Node::Lipschitz(a, k, frame))
    }

    // Convenience builders.
    pub fn add(&mut self, a: NodeId, b: NodeId) -> NodeId {
        self.bin(BinOp::Add, a, b)
    }
    pub fn sub(&mut self, a: NodeId, b: NodeId) -> NodeId {
        self.bin(BinOp::Sub, a, b)
    }
    pub fn mul(&mut self, a: NodeId, b: NodeId) -> NodeId {
        self.bin(BinOp::Mul, a, b)
    }
    pub fn max(&mut self, a: NodeId, b: NodeId) -> NodeId {
        self.bin(BinOp::Max, a, b)
    }
    pub fn min(&mut self, a: NodeId, b: NodeId) -> NodeId {
        self.bin(BinOp::Min, a, b)
    }
    pub fn neg(&mut self, a: NodeId) -> NodeId {
        self.un(UnOp::Neg, a)
    }

    /// All nodes reachable from `roots`, sorted ascending (a valid evaluation order).
    pub fn reachable(&self, roots: &[NodeId]) -> Vec<NodeId> {
        self.reach(roots, true)
    }

    /// The nodes `roots`' values are built of: `reachable`, but not through
    /// what only decides between them (a select's condition) or what a
    /// derivative is taken of. A name found here is a piece of the shape;
    /// one read only to decide (std's fillets test a ball against each
    /// operand, moved) is not.
    pub fn pieces(&self, roots: &[NodeId]) -> Vec<NodeId> {
        let mut seen: HashSet<NodeId> = HashSet::default();
        let mut stack: Vec<NodeId> = roots.to_vec();
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            let (ch, k) = match self.nodes[n as usize] {
                Node::Select(_, a, b) => ([a, b, 0], 2),
                Node::Deriv(..) => ([0; 3], 0),
                node => node.children(),
            };
            stack.extend(ch[..k].iter().filter(|c| !seen.contains(c)));
        }
        let mut out: Vec<NodeId> = seen.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// `reachable`, through a read's value (Node::Read's `a`) or not: what
    /// remakes reads from their coordinates (rewriting, differentiating,
    /// fingerprints) needs only those.
    fn reach(&self, roots: &[NodeId], values: bool) -> Vec<NodeId> {
        // Cost proportional to the subgraph, not the whole graph (graphs grow
        // large when big assemblies are flattened; most queries touch little).
        let mut seen: HashSet<NodeId> = HashSet::default();
        let mut stack: Vec<NodeId> = roots.to_vec();
        let mut out = Vec::new();
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            out.push(n);
            let (ch, k) = match self.nodes[n as usize] {
                Node::Read(_, _, _, at) if !values => ([at, 0, 0], 1),
                node => node.children(),
            };
            for &c in &ch[..k] {
                if !seen.contains(&c) {
                    stack.push(c);
                }
            }
        }
        out.sort_unstable();
        out
    }

    /// Substitute variables (this is `warp`). Lipschitz promises follow
    /// where they still hold (`rebuild`).
    pub fn subst(&mut self, root: NodeId, map: &[(Var, NodeId)]) -> NodeId {
        let nodes: HashMap<NodeId, NodeId> = map.iter().map(|(v, n)| (self.intern(Node::Var(*v)), *n)).collect();
        self.rebuild(root, &nodes)
    }

    /// `root` with the nodes `map` names replaced, rebuilt above them (a
    /// pair's placements relative to a moving frame: ink-lang's
    /// motion::moving_frame).
    /// Lipschitz promises follow where they still hold (`rebuild`).
    pub fn replace(&mut self, root: NodeId, map: &std::collections::HashMap<NodeId, NodeId>) -> NodeId {
        self.rebuild(root, map)
    }

    /// `root` without its names: each tag replaced by the field it names,
    /// every other node kept (the same values, derivatives and promises),
    /// for a reader of values alone, which then sees one field however it
    /// was named (a package's solve that reads no face region's provenance).
    pub fn untagged(&mut self, root: NodeId) -> NodeId {
        self.rebuild_as(root, &HashMap::default(), true, None)
    }

    /// `root` with `map`'s nodes replaced. A promise holds in its frame's
    /// coordinates for every value of t and the params, so it survives
    /// variables replaced: x, y or z by anything (the frame is rewritten
    /// with the field), t or a param by an expression of them. It is
    /// dropped where it reads a variable of t or a param replaced by an
    /// expression of x, y or z (a slope d/dt times dt/dx it does not
    /// cover), or any other node replaced (it held of the field as
    /// written, not of one with a part swapped: x (cos^2 p + sin^2 p) is
    /// 1-Lipschitz, x (1000 + sin^2 p) is not).
    fn rebuild<S: std::hash::BuildHasher>(&mut self, root: NodeId, map: &std::collections::HashMap<NodeId, NodeId, S>) -> NodeId {
        self.rebuild_as(root, map, false, None)
    }

    /// `root` with each of `maxes` (`max` nodes, with the operands under
    /// them grouped) written as the `emax` of each group's max: the same
    /// sign, another value where groups are positive (distance::exact,
    /// which groups operands over square directions). A promise above one
    /// is dropped: it held of the max.
    pub fn with_emax(&mut self, root: NodeId, maxes: &HashMap<NodeId, Vec<Vec<NodeId>>>) -> NodeId {
        self.rebuild_as(root, &HashMap::default(), false, Some(maxes))
    }

    /// `rebuild`, and with `untag` each name replaced by what it names
    /// (`untagged`): the one rewrite every other is, so each keeps spatial
    /// arithmetic real where its operands became constants (an untagged
    /// p - tag(p) is 0, and 1e16 + 0.5 above it still its real value).
    /// Untagging changes no value, so a derivative's mark stays.
    fn rebuild_as<S: std::hash::BuildHasher>(
        &mut self,
        root: NodeId,
        map: &std::collections::HashMap<NodeId, NodeId, S>,
        untag: bool,
        emax: Option<&HashMap<NodeId, Vec<Vec<NodeId>>>>,
    ) -> NodeId {
        let ids = self.reach(&[root], false);
        let mut memo: HashMap<NodeId, NodeId> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        // The nodes reading a replacement that breaks promises.
        let mut breaks: HashSet<NodeId> = HashSet::default();
        for &id in &ids {
            if let Some(&to) = map.get(&id) {
                let kept = match self.nodes[id as usize] {
                    Node::Var(Var::X | Var::Y | Var::Z) => true,
                    Node::Var(_) => !self.depends_on(to, |v| matches!(v, Var::X | Var::Y | Var::Z)),
                    _ => false,
                };
                if !kept {
                    breaks.insert(id);
                }
                // A coordinate replaced by params (a moving point's fixed
                // one) is still a coordinate: arithmetic on it is real.
                let to = if matches!(self.nodes[id as usize], Node::Var(Var::X | Var::Y | Var::Z)) && self.as_const(to).is_none() { self.real(to) } else { to };
                memo.insert(id, to);
                continue;
            }
            let (ch, k) = match self.nodes[id as usize] {
                Node::Read(_, _, _, at) => ([at, 0, 0], 1),
                node => node.children(),
            };
            if ch[..k].iter().any(|c| breaks.contains(c)) {
                breaks.insert(id);
            }
            let new = match self.nodes[id as usize] {
                Node::Const(_) | Node::Var(_) => id,
                // An operation that read a spatial value (x, a moment's t)
                // keeps its meaning in real numbers where its operands
                // became constants and params (`un_real`): x - y along
                // x = y is 0, and 1e16 + 0.5 above it is not its float.
                // One of constants and params as written stays one, taken
                // in floats as a def is.
                Node::Un(op, a) => {
                    let a = memo[&a];
                    if !self.fixed[id as usize] && self.fixed[a as usize] { self.un_real(op, a) } else { self.un(op, a) }
                }
                Node::Bin(..) if let Some(groups) = emax.and_then(|m| m.get(&id)) => {
                    breaks.insert(id);
                    let maxes: Vec<NodeId> = groups.iter().map(|grp| grp[1..].iter().fold(memo[&grp[0]], |acc, o| self.bin(BinOp::Max, acc, memo[o]))).collect();
                    maxes[1..].iter().fold(maxes[0], |acc, &m| self.bin(BinOp::Emax, acc, m))
                }
                Node::Bin(op, a, b) => {
                    let (a, b) = (memo[&a], memo[&b]);
                    if !self.fixed[id as usize] && self.fixed[a as usize] && self.fixed[b as usize] { self.bin_real(op, a, b) } else { self.bin(op, a, b) }
                }
                Node::Real(a) => {
                    let a = memo[&a];
                    self.real(a)
                }
                Node::Within(..) => id,
                Node::Select(c, a, b) => {
                    let (c, a, b) = (memo[&c], memo[&a], memo[&b]);
                    self.select(c, a, b)
                }
                Node::Tag(_, a) if untag => memo[&a],
                Node::Tag(t, a) => {
                    let a = memo[&a];
                    self.intern(Node::Tag(t, a))
                }
                Node::Grid(g, m, c) => {
                    let c = c.map(|n| memo[&n]);
                    self.intern(Node::Grid(g, m, c))
                }
                // Read again at the rewritten coordinates (its value is
                // theirs: a map cannot reach into it).
                Node::Read(g, path, _, at) => {
                    let Node::Frame(c) = self.nodes[at as usize] else { unreachable!("a read's coordinates are a frame") };
                    let c = c.map(|n| memo[&n]);
                    self.read(g, path, c)
                }
                Node::Lipschitz(a, k, f) => {
                    let (a, k, f) = (memo[&a], memo[&k], memo[&f]);
                    if self.as_const(a).is_some() || breaks.contains(&id) { a } else { self.intern(Node::Lipschitz(a, k, f)) }
                }
                Node::Frame(c) => {
                    let c = c.map(|n| memo[&n]);
                    self.intern(Node::Frame(c))
                }
                // Rewritten, a is no longer the derivative of f along the
                // axis (a warp's chain rule changes it): plain a. Along t
                // the mark stays on the rewritten f, which is what jumps
                // (its speed is taken from f, not from a).
                Node::Deriv(a0, f0, k) => {
                    let (a, f) = (memo[&a0], memo[&f0]);
                    if a == a0 && f == f0 {
                        id
                    } else if untag {
                        self.intern(Node::Deriv(a, f, k))
                    } else if k == 3 && self.depends_on(f, |v| v == Var::T) {
                        self.intern(Node::Deriv(a, f, 3))
                    } else {
                        a
                    }
                }
            };
            memo.insert(id, new);
        }
        memo[&root]
    }

    /// Derivative of min/max: the branch taken; at a tie, the mean of both,
    /// da/2 + db/2 (as the evaluators do): halving each first, since the sum
    /// of constant sides can overflow where their mean does not (the mean of
    /// 1e308 and 1.5e308). One constant slope on both sides is that slope
    /// where its halves are exact; below 2^-1021 they round (the mean of
    /// 5e-324 and 5e-324 is 0), so it is the mean there too.
    fn d_branch(&mut self, op: BinOp, a: NodeId, b: NodeId, da: NodeId, db: NodeId) -> NodeId {
        if da == db && self.as_const(da).is_some_and(|v| 0.5 * v + 0.5 * v == v) {
            return da;
        }
        let (lo, hi) = if op == BinOp::Min { (da, db) } else { (db, da) };
        let c = self.sub(a, b);
        let nc = self.neg(c);
        let half = self.constant(0.5);
        let (ha, hb) = (self.mul(half, da), self.mul(half, db));
        let mean = self.add(ha, hb);
        let tie_or_hi = self.select(nc, hi, mean);
        self.select(c, lo, tie_or_hi)
    }

    /// Symbolic derivative of `root` with respect to a variable, as a new
    /// expression. min/max/abs/select differentiate the branch taken; a
    /// sampled field's read, by the chain rule through its coordinates,
    /// with the reads one axis further along its path (Node::Read: about
    /// its source, unknown). Each op's derivative is written as the
    /// evaluators compute it (eval::eval_grads, eval_dual): its partials,
    /// from its operands' values (dbin, dun; a read's: the reads one axis
    /// further), times its operands' derivatives, summed in order. So the
    /// derivative's value is theirs, bit for bit, at every scale.
    pub fn diff(&mut self, root: NodeId, var: Var) -> NodeId {
        let ids = self.reach(&[root], false);
        let mut d: HashMap<NodeId, NodeId> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        let zero = self.constant(0.0);
        let one = self.constant(1.0);
        for &id in &ids {
            let node = self.nodes[id as usize];
            let dn = match node {
                Node::Const(_) => zero,
                Node::Var(v) => {
                    if v == var {
                        one
                    } else {
                        zero
                    }
                }
                Node::Un(op, a) => {
                    let da = d[&a];
                    if da == zero {
                        zero
                    } else {
                        match op {
                            UnOp::Neg => self.neg(da),
                            UnOp::Abs => {
                                // At 0 the average of both sides: 0.
                                let n = self.neg(da);
                                let na = self.neg(a);
                                let pos = self.select(na, da, zero);
                                self.select(a, n, pos)
                            }
                            UnOp::Sqrt => {
                                // (0.5 / sqrt a) da where a > 0, else 0.
                                let half = self.constant(0.5);
                                let s = self.un(UnOp::Sqrt, a);
                                let f = self.bin(BinOp::Div, half, s);
                                let q = self.mul(f, da);
                                let na = self.neg(a);
                                self.select(na, q, zero)
                            }
                            UnOp::Square => {
                                let two = self.constant(2.0);
                                let t = self.mul(two, a);
                                self.mul(t, da)
                            }
                            UnOp::Exp => self.mul(id, da),
                            UnOp::Ln => {
                                // ln of a non-positive is a constant (-inf, or 0).
                                let f = self.bin(BinOp::Div, one, a);
                                let q = self.mul(f, da);
                                let na = self.neg(a);
                                self.select(na, q, zero)
                            }
                            UnOp::Sin => {
                                let c = self.un(UnOp::Cos, a);
                                self.mul(c, da)
                            }
                            UnOp::Cos => {
                                let s = self.un(UnOp::Sin, a);
                                let n = self.neg(s);
                                self.mul(n, da)
                            }
                            UnOp::Tan => {
                                let t2 = self.un(UnOp::Square, id);
                                let f = self.add(one, t2);
                                self.mul(f, da)
                            }
                            UnOp::Floor => zero,
                            // Linear: the same exact scaling of the
                            // derivative, read as a real number where it is
                            // constants and params (never folded into a
                            // rounded coefficient).
                            UnOp::Rad | UnOp::Deg => {
                                let r = self.real(da);
                                self.un(op, r)
                            }
                        }
                    }
                }
                Node::Bin(op, a, b) => {
                    let (da, db) = (d[&a], d[&b]);
                    match op {
                        BinOp::Add => self.add(da, db),
                        BinOp::Sub => self.sub(da, db),
                        BinOp::Mul => {
                            let x = self.mul(da, b);
                            let y = self.mul(a, db);
                            self.add(x, y)
                        }
                        BinOp::Div => {
                            // (1 / b) da - ((a / b) / b) db (b^2 underflows to
                            // 0 for tiny b; x / 0 = 0 at b = 0).
                            let fa = self.bin(BinOp::Div, one, b);
                            let x = self.mul(fa, da);
                            let q = self.bin(BinOp::Div, a, b);
                            let qb = self.bin(BinOp::Div, q, b);
                            let y = self.mul(qb, db);
                            let r = self.sub(x, y);
                            let nb = self.un(UnOp::Abs, b);
                            let nb = self.neg(nb);
                            self.select(nb, r, zero)
                        }
                        BinOp::Min | BinOp::Max => self.d_branch(op, a, b, da, db),
                        BinOp::Emax => {
                            // Both positive: (a / emax) da + (b / emax) db (dividing
                            // first, as the evaluators do: a da can overflow), each
                            // weight at most 1, as it is in real numbers (so bounds
                            // of a slope through it see that); else as max.
                            let e = self.bin(BinOp::Emax, a, b);
                            let (ae, be) = (self.bin(BinOp::Div, a, e), self.bin(BinOp::Div, b, e));
                            let (ae, be) = (self.min(ae, one), self.min(be, one));
                            let x = self.mul(ae, da);
                            let y = self.mul(be, db);
                            let round = self.add(x, y);
                            let sharp = self.d_branch(BinOp::Max, a, b, da, db);
                            let both = self.min(a, b);
                            let nb = self.neg(both);
                            self.select(nb, round, sharp)
                        }
                        BinOp::Atan2 => {
                            // atan2(a, b): (b / r2) da + (-a / r2) db, with
                            // r2 = a^2 + b^2 > 0.
                            let a2 = self.un(UnOp::Square, a);
                            let b2 = self.un(UnOp::Square, b);
                            let den = self.add(a2, b2);
                            let fa = self.bin(BinOp::Div, b, den);
                            let na = self.neg(a);
                            let fb = self.bin(BinOp::Div, na, den);
                            let x = self.mul(fa, da);
                            let y = self.mul(fb, db);
                            let q = self.add(x, y);
                            let nd = self.neg(den);
                            self.select(nd, q, zero)
                        }
                        BinOp::Pow => {
                            // Where a < 0 and b is not a whole number the value is
                            // the constant 0: so is its derivative.
                            let fl = self.un(UnOp::Floor, b);
                            let fr = self.sub(b, fl);
                            let fr = self.un(UnOp::Abs, fr);
                            let nfr = self.neg(fr);
                            let zero_here = self.max(a, nfr);
                            if db == zero {
                                // b a^(b-1) da
                                let bm1 = self.sub(b, one);
                                let p = self.bin(BinOp::Pow, a, bm1);
                                let t = self.mul(b, p);
                                let r = self.mul(t, da);
                                self.select(zero_here, zero, r)
                            } else {
                                // b a^(b-1) da + a^b ln(a) db, the second term
                                // only where a > 0 (elsewhere a^b does not vary
                                // smoothly with b), as the evaluators.
                                let bm1 = self.sub(b, one);
                                let p = self.bin(BinOp::Pow, a, bm1);
                                let fa = self.mul(b, p);
                                let x = self.mul(fa, da);
                                let la = self.un(UnOp::Ln, a);
                                let vl = self.mul(id, la);
                                let na = self.neg(a);
                                let fb = self.select(na, vl, zero);
                                let y = self.mul(fb, db);
                                let r = self.add(x, y);
                                self.select(zero_here, zero, r)
                            }
                        }
                    }
                }
                Node::Select(c, a, b) => {
                    let (da, db) = (d[&a], d[&b]);
                    self.select(c, da, db)
                }
                Node::Tag(_, a) | Node::Lipschitz(a, _, _) | Node::Deriv(a, _, _) | Node::Real(a) => d[&a],
                Node::Frame(_) | Node::Within(..) => zero,
                // The interpolant: chain rule through the coordinates, with
                // the component with one more axis (0 along an axis it
                // already has, inside a cell).
                Node::Grid(g, m, c) => {
                    let mut acc = zero;
                    for k in 0..3 {
                        let dc = d[&c[k]];
                        if dc == zero || m >> k & 1 == 1 {
                            continue;
                        }
                        let slope = self.intern(Node::Grid(g, m | 1 << k, c));
                        let term = self.mul(slope, dc);
                        acc = self.add(acc, term);
                    }
                    acc
                }
                // A read: chain rule through its coordinates, with the read
                // one axis further along its path (whatever its value, the
                // samples', simplifies to: 0 along an axis twice in a cell,
                // where the source's is not). Where the coordinates do not
                // move with `var`, neither does the source: 0.
                Node::Read(g, path, _, at) => {
                    let Node::Frame(c) = self.nodes[at as usize] else { unreachable!("a read's coordinates are a frame") };
                    let mut acc = zero;
                    for k in 0..3 {
                        let dc = d[&c[k]];
                        if dc == zero {
                            continue;
                        }
                        let mut along = self.paths[path as usize].clone();
                        along.push(k as u8);
                        let along = self.path(along);
                        let slope = self.read(g, along, c);
                        let term = self.mul(slope, dc);
                        acc = self.add(acc, term);
                    }
                    acc
                }
            };
            d.insert(id, dn);
        }
        d[&root]
    }

    /// `diff`, and where `var` is x, y or z, marked as the derivative of
    /// `root` (Node::Deriv), so ranges also hold it within root's gradient
    /// bounds; where it is t and root reads t, marked so even where it is
    /// constant, so a claim on a speed checks the motion for jumps. What
    /// the language's `diff` and `grad` make.
    pub fn derivative(&mut self, root: NodeId, var: Var) -> NodeId {
        let a = self.diff(root, var);
        let axis = match var {
            Var::X => 0,
            Var::Y => 1,
            Var::Z => 2,
            // Through marks: a speed whose formula is constant still marks
            // its own drive, so an acceleration checks both.
            Var::T if self.depends_on(root, |v| v == Var::T) => 3,
            _ => return a,
        };
        if axis < 3 && self.as_const(a).is_some() {
            return a;
        }
        self.intern(Node::Deriv(a, root, axis))
    }

    /// `derivative` of `root` along t, marked as the speed of each of
    /// `read` instead of root's: what reads root (a turn's sine and cosine
    /// of its angle), each continuous where root jumps by what none of them
    /// sees (a whole turn), so its speed is the formula's there.
    pub fn speed_read(&mut self, root: NodeId, read: &[NodeId]) -> NodeId {
        let a = self.diff(root, Var::T);
        if !self.depends_on(root, |v| v == Var::T) {
            return a;
        }
        read.iter().fold(a, |m, &f| self.intern(Node::Deriv(m, f, 3)))
    }

    /// Content hash of the expressions at `roots`, with parameters replaced by
    /// their values: equal hashes mean equal functions of x, y, z, t (up to
    /// hash collisions), wherever and however they were built. Stable across
    /// runs and machines; for caching results. Construction names (tags)
    /// count only `with_names`: they matter to provenance (named faces), not
    /// to values.
    pub fn fingerprint(&self, roots: &[NodeId], params: &[f64], with_names: bool, fp: &mut Fingerprint) {
        let h = self.fingerprints(roots, params, with_names);
        for r in roots {
            fp.key(h[r]);
        }
    }

    /// `fingerprint`'s content hash of each node the `roots` reach, for a
    /// key with structure: a record writes each of its fields' hashes in
    /// the slot the field fills.
    pub fn fingerprints(&self, roots: &[NodeId], params: &[f64], with_names: bool) -> std::collections::HashMap<NodeId, u128> {
        let ids = self.reach(roots, false);
        let mut h: HashMap<NodeId, u128> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        for &id in &ids {
            if let Node::Tag(_, a) = self.nodes[id as usize]
                && !with_names
            {
                h.insert(id, h[&a]); // names are transparent
                continue;
            }
            let mut f = Fingerprint::default();
            match self.nodes[id as usize] {
                Node::Const(v) => f.word(1).f64(v),
                Node::Var(Var::Param(k)) => f.word(2).f64(params.get(k as usize).copied().unwrap_or(f64::NAN)),
                Node::Var(v) => f.word(3).word(match v {
                    Var::X => 0,
                    Var::Y => 1,
                    Var::Z => 2,
                    _ => 3,
                }),
                Node::Un(op, a) => f.word(4).word(op as u64).key(h[&a]),
                Node::Bin(op, a, b) => f.word(5).word(op as u64).key(h[&a]).key(h[&b]),
                Node::Select(c, a, b) => f.word(6).key(h[&c]).key(h[&a]).key(h[&b]),
                // Whether provenance stops at the name counts too: it changes
                // which names a face reports.
                Node::Tag(t, a) => f.word(7).str(&self.tag_names[t as usize]).word(self.opaque.contains(&t) as u64).key(h[&a]),
                Node::Grid(g, m, c) => {
                    self.grids[g as usize].fingerprint(f.word(8).word(m as u64));
                    c.iter().for_each(|n| _ = f.key(h[n]));
                    &mut f
                }
                // Its value is made from these (Graph::read).
                Node::Read(g, path, _, at) => {
                    self.grids[g as usize].fingerprint(f.word(11));
                    f.word(self.paths[path as usize].len() as u64);
                    self.paths[path as usize].iter().for_each(|&k| _ = f.word(k as u64));
                    f.key(h[&at])
                }
                Node::Lipschitz(a, k, fr) => f.word(9).key(h[&a]).key(h[&k]).key(h[&fr]),
                Node::Frame(c) => f.word(10).key(h[&c[0]]).key(h[&c[1]]).key(h[&c[2]]),
                Node::Deriv(a, fd, k) => f.word(12).word(k as u64).key(h[&a]).key(h[&fd]),
                Node::Real(a) => f.word(13).key(h[&a]),
                Node::Within(l, hi, v, id) => f.word(14).word(id).key(h[&l]).key(h[&hi]).key(h[&v]),
            };
            h.insert(id, f.value());
        }
        h.into_iter().collect()
    }

    /// The sources of the sampled fields the fields at `roots` read (what
    /// each was sampled from, and how far it may be from it), each once, in
    /// the order the grids were made.
    pub fn sources(&self, roots: &[NodeId]) -> Vec<std::sync::Arc<crate::grid::Source>> {
        let mut read = vec![false; self.grids.len()];
        for n in self.reachable(roots) {
            if let Node::Grid(g, ..) | Node::Read(g, ..) = self.nodes[n as usize] {
                read[g as usize] = true;
            }
        }
        let mut out: Vec<std::sync::Arc<crate::grid::Source>> = Vec::new();
        for (g, _) in read.iter().enumerate().filter(|r| *r.1) {
            if let Some(s) = self.grids[g].source()
                && !out.iter().any(|o| std::sync::Arc::ptr_eq(o, s))
            {
                out.push(s.clone());
            }
        }
        out
    }

    pub fn depends_on(&self, root: NodeId, pred: impl Fn(Var) -> bool) -> bool {
        self.reachable(&[root]).iter().any(|&n| matches!(self.nodes[n as usize], Node::Var(v) if pred(v)))
    }

    /// Whether `root`'s value reads a variable `pred` takes: `depends_on`,
    /// but not through what only bounds it (a Lipschitz promise's frame and
    /// constant, a derivative's mark), so a profile in x and y with a
    /// distance promise does not read z.
    pub fn value_depends_on(&self, root: NodeId, pred: impl Fn(Var) -> bool) -> bool {
        let mut seen = HashSet::default();
        let mut stack = vec![root];
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            let (ch, k) = match self.nodes[n as usize] {
                Node::Var(v) if pred(v) => return true,
                Node::Lipschitz(a, ..) | Node::Deriv(a, ..) => ([a, 0, 0], 1),
                node => node.children(),
            };
            stack.extend(&ch[..k]);
        }
        false
    }
}

// ---- scalar semantics (shared by folding and every evaluator) ----
//
// The core's value rules (docs/CORE.md): the IEEE-754 result, with one zero
// (-0 reads as 0), x / 0 = 0, and anything else undefined (NaN) = 0. So
// sqrt and ln of a negative are 0, pow of a negative base with an exponent
// that is not a whole number is 0; ln 0 is -inf. conformance/ tests every
// evaluator against an independent reference.

/// Where each operation is undefined (docs/IR.md, Fields, Value rules: its
/// value is 0, and it is reported): the one table of operator domains, read
/// over operands' ranges [lo, hi] (a point for a value). The case an op
/// may meet there, if any. Every check reads it: the evaluators' report and
/// constant folding (at points: `undefined_un`, `undefined_bin`), and
/// whether a promise's constants are real numbers (eval::promise_valid).
/// sqrt of a negative is no such case: it is sqrt(max(a, 0)), continuous.
pub fn undefined_un_over(op: UnOp, [lo, hi]: [f64; 2]) -> Option<&'static str> {
    match op {
        UnOp::Ln if lo < 0.0 => Some("ln of a negative"),
        UnOp::Sin | UnOp::Cos | UnOp::Tan if lo.is_infinite() || hi.is_infinite() => Some("sin, cos or tan of infinity"),
        UnOp::Neg | UnOp::Abs | UnOp::Square | UnOp::Exp | UnOp::Ln | UnOp::Sqrt | UnOp::Sin | UnOp::Cos | UnOp::Tan | UnOp::Floor | UnOp::Rad | UnOp::Deg => None,
    }
}

/// See `undefined_un_over`. A NaN bound may be anything.
pub fn undefined_bin_over(op: BinOp, a: [f64; 2], b: [f64; 2]) -> Option<&'static str> {
    let nan = a.iter().chain(&b).any(|x| x.is_nan());
    let (zero, infinite) = (|r: [f64; 2]| r[0] <= 0.0 && r[1] >= 0.0, |r: [f64; 2]| r[0] == f64::NEG_INFINITY || r[1] == f64::INFINITY);
    match op {
        BinOp::Div if zero(b) => Some("division by zero"),
        BinOp::Div if nan || (infinite(a) && infinite(b)) => Some("infinity divided by infinity"),
        BinOp::Add if nan || (a[1] == f64::INFINITY && b[0] == f64::NEG_INFINITY) || (a[0] == f64::NEG_INFINITY && b[1] == f64::INFINITY) => Some("infinity minus infinity"),
        BinOp::Sub if nan || (a[1] == f64::INFINITY && b[1] == f64::INFINITY) || (a[0] == f64::NEG_INFINITY && b[0] == f64::NEG_INFINITY) => Some("infinity minus infinity"),
        BinOp::Mul if nan || (zero(a) && infinite(b)) || (zero(b) && infinite(a)) => Some("zero times infinity"),
        BinOp::Pow if a[0] < 0.0 && !(b[0] == b[1] && b[0].fract() == 0.0) => Some("a negative base to a power that is not a whole number"),
        BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Div | BinOp::Pow | BinOp::Min | BinOp::Max | BinOp::Atan2 | BinOp::Emax => None,
    }
}

/// The undefined case an op meets at these arguments (`undefined_un_over`
/// at a point).
pub fn undefined_un(op: UnOp, a: f64) -> Option<&'static str> {
    undefined_un_over(op, [a, a])
}

/// See `undefined_un`.
pub fn undefined_bin(op: BinOp, a: f64, b: f64) -> Option<&'static str> {
    undefined_bin_over(op, [a, a], [b, b])
}

/// The core's value for an IEEE result: NaN (undefined) is 0. Library
/// functions (exp, ln, sin, cos, tan, atan, atan2, pow, hypot) are the
/// `libm` crate's, plain f64 code with no CPU dispatch, not the system's:
/// a model means the same, bit for bit, on every machine, so its answers
/// do too (and the result cache needs no maths library in its keys).
#[inline]
pub fn defined(v: f64) -> f64 {
    if v.is_nan() { 0.0 } else { v }
}

/// `a / b` with `x / 0 = 0`.
#[inline]
pub fn div_f64(a: f64, b: f64) -> f64 {
    if b == 0.0 { 0.0 } else { defined(a / b) }
}

#[inline]
pub fn pow_f64(a: f64, b: f64) -> f64 {
    // + 0.0: one zero, so pow(-0, -1) is +inf like pow(0, -1).
    if a < 0.0 && b.fract() != 0.0 { 0.0 } else { libm::pow(a + 0.0, b) }
}

#[inline]
pub fn ln_f64(a: f64) -> f64 {
    if a > 0.0 {
        libm::log(a)
    } else if a == 0.0 {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

#[inline]
pub fn atan2_f64(y: f64, x: f64) -> f64 {
    // + 0.0: one zero, so the negative x axis is +pi from either side of 0.
    libm::atan2(y + 0.0, x + 0.0)
}

#[inline]
pub fn un_f64(op: UnOp, a: f64) -> f64 {
    match op {
        UnOp::Neg => -a,
        UnOp::Abs => a.abs(),
        UnOp::Sqrt => a.max(0.0).sqrt(),
        UnOp::Square => a * a,
        UnOp::Exp => libm::exp(a),
        UnOp::Ln => ln_f64(a),
        UnOp::Sin => defined(libm::sin(a)),
        UnOp::Cos => defined(libm::cos(a)),
        UnOp::Tan => defined(libm::tan(a)),
        UnOp::Floor => a.floor(),
        UnOp::Rad => a * (std::f64::consts::PI / 180.0),
        UnOp::Deg => a * (180.0 / std::f64::consts::PI),
    }
}

#[inline]
pub fn bin_f64(op: BinOp, a: f64, b: f64) -> f64 {
    match op {
        BinOp::Add => defined(a + b),
        BinOp::Sub => defined(a - b),
        BinOp::Mul => defined(a * b),
        BinOp::Div => div_f64(a, b),
        BinOp::Min => a.min(b),
        BinOp::Max => a.max(b),
        BinOp::Atan2 => atan2_f64(a, b),
        BinOp::Pow => pow_f64(a, b),
        BinOp::Emax => emax_f64(a, b),
    }
}

#[inline]
pub fn emax_f64(a: f64, b: f64) -> f64 {
    if a > 0.0 && b > 0.0 { libm::hypot(a, b) } else { a.max(b) }
}

/// Below this, a product's or quotient's rounding error may itself round
/// away, so its float is not shown exact.
const EXACT_TINY: f64 = 1e-270;

/// `op` of `a` where its float is its real value, shown so (an exact
/// square root, a square without rounding, sin 0), else None: what a
/// rewrite may fold an operation that read a spatial value into
/// (`Graph::un_real`).
pub fn exact_un(op: UnOp, a: f64) -> Option<f64> {
    if !a.is_finite() || undefined_un(op, a).is_some() {
        return None;
    }
    let v = un_f64(op, a);
    let exact = match op {
        UnOp::Neg | UnOp::Abs | UnOp::Floor => true,
        UnOp::Square => a == 0.0 || (v.is_finite() && v >= EXACT_TINY && a.mul_add(a, -v) == 0.0),
        UnOp::Sqrt => a <= 0.0 || (a >= EXACT_TINY && v.mul_add(v, -a) == 0.0),
        UnOp::Exp | UnOp::Cos => a == 0.0,
        UnOp::Sin | UnOp::Tan => a == 0.0,
        UnOp::Ln => a == 1.0,
        UnOp::Rad | UnOp::Deg => a == 0.0,
    };
    exact.then_some(v)
}

/// `exact_un` for two operands: a sum or difference with no rounding
/// error, a product or quotient that is exact, min and max.
pub fn exact_bin(op: BinOp, a: f64, b: f64) -> Option<f64> {
    if !a.is_finite() || !b.is_finite() || undefined_bin(op, a, b).is_some() {
        return None;
    }
    let v = bin_f64(op, a, b);
    if !v.is_finite() {
        return None;
    }
    // The sum's rounding error (Knuth's two-sum), which is a float.
    let two_sum = |a: f64, b: f64| {
        let s = a + b;
        let bb = s - a;
        (a - (s - bb)) + (b - bb)
    };
    let exact = match op {
        BinOp::Add => two_sum(a, b) == 0.0,
        BinOp::Sub => two_sum(a, -b) == 0.0,
        BinOp::Mul => a == 0.0 || b == 0.0 || (v.abs() >= EXACT_TINY && a.mul_add(b, -v) == 0.0),
        BinOp::Div => a == 0.0 || (v.abs() >= EXACT_TINY && a.abs() >= EXACT_TINY && v.mul_add(b, -a) == 0.0),
        BinOp::Min | BinOp::Max => true,
        BinOp::Emax => !(a > 0.0 && b > 0.0),
        BinOp::Atan2 | BinOp::Pow => false,
    };
    exact.then_some(v)
}

// Derivatives: each formula evaluated under the value rules, on the branch
// the value took; at an exact tie between branches (min, max, abs at 0) the
// average of both sides, so the result does not depend on argument order.

/// Value and derivative d/da.
#[inline]
pub fn dun(op: UnOp, a: f64) -> (f64, f64) {
    match op {
        UnOp::Neg => (-a, -1.0),
        UnOp::Abs => (
            a.abs(),
            if a < 0.0 {
                -1.0
            } else if a > 0.0 {
                1.0
            } else {
                0.0
            },
        ),
        UnOp::Sqrt => {
            let v = a.max(0.0).sqrt();
            (v, if v > 0.0 { 0.5 / v } else { 0.0 })
        }
        UnOp::Square => (a * a, 2.0 * a),
        UnOp::Exp => {
            let e = libm::exp(a);
            (e, e)
        }
        UnOp::Ln => (ln_f64(a), if a > 0.0 { 1.0 / a } else { 0.0 }),
        UnOp::Sin => (defined(libm::sin(a)), defined(libm::cos(a))),
        UnOp::Cos => (defined(libm::cos(a)), defined(-libm::sin(a))),
        UnOp::Tan => {
            let t = defined(libm::tan(a));
            (t, 1.0 + t * t)
        }
        UnOp::Floor => (a.floor(), 0.0),
        UnOp::Rad => {
            let k = std::f64::consts::PI / 180.0;
            (a * k, k)
        }
        UnOp::Deg => {
            let k = 180.0 / std::f64::consts::PI;
            (a * k, k)
        }
    }
}

/// Value and partials d/da, d/db.
#[inline]
pub fn dbin(op: BinOp, a: f64, b: f64) -> (f64, f64, f64) {
    match op {
        BinOp::Add => (defined(a + b), 1.0, 1.0),
        BinOp::Sub => (defined(a - b), 1.0, -1.0),
        BinOp::Mul => (defined(a * b), b, a),
        BinOp::Div => {
            if b == 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                // -(a / b) / b, not -a / b^2: b^2 underflows to 0 for tiny b.
                (defined(a / b), 1.0 / b, defined(-(a / b) / b))
            }
        }
        // At exact ties the average of the two branch derivatives: surface
        // points often land exactly on ties (e.g. min(d, 0) on a box face),
        // and a sharp edge gets the direction between its two faces.
        BinOp::Min => {
            if a == b {
                (a, 0.5, 0.5)
            } else if a < b {
                (a, 1.0, 0.0)
            } else {
                (b, 0.0, 1.0)
            }
        }
        BinOp::Max => {
            if a == b {
                (a, 0.5, 0.5)
            } else if a > b {
                (a, 1.0, 0.0)
            } else {
                (b, 0.0, 1.0)
            }
        }
        BinOp::Atan2 => {
            let r2 = a * a + b * b;
            if r2 > 0.0 { (atan2_f64(a, b), defined(b / r2), defined(-a / r2)) } else { (atan2_f64(a, b), 0.0, 0.0) }
        }
        BinOp::Pow => {
            let v = pow_f64(a, b);
            // 0 where the value is the constant 0 (a negative base, an
            // exponent that is not a whole number); b - 1 can round to a
            // whole number, so test b itself.
            let fa = if b == 0.0 || (a < 0.0 && b.fract() != 0.0) { 0.0 } else { defined(b * pow_f64(a, b - 1.0)) };
            let fb = if a > 0.0 { defined(v * libm::log(a)) } else { 0.0 };
            (v, fa, fb)
        }
        BinOp::Emax => {
            if a > 0.0 && b > 0.0 {
                // Each weight at most 1 (hypot may round below a or b).
                let e = libm::hypot(a, b);
                (e, defined(a / e).min(1.0), defined(b / e).min(1.0))
            } else {
                dbin(BinOp::Max, a, b)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_number_known_within_bounds_is_read_over_them() {
        // Ranges read all of its bounds, through arithmetic on constants
        // too (never folded to its value); a point reads its value; it is
        // constant in space.
        use crate::eval::{Ctx, eval_interval, eval_point};
        use crate::interval::I;
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let k = g.within(1.0, 2.0, 1.5, 1);
        let two = g.constant(2.0);
        let k2 = g.mul(k, two);
        let f = g.sub(x, k2);
        let t = crate::tape::Tape::new(&g, f);
        let ctx = Ctx::new(0.0, vec![]);
        let p = I::point(0.0);
        assert_eq!(eval_interval(&t, &ctx, [p, p, p], None), I::new(-4.0, -2.0));
        assert_eq!(eval_point(&t, &ctx, [0.0; 3]), -3.0);
        let d = g.diff(f, Var::X);
        assert_eq!(g.as_const(d), Some(1.0));
        assert_eq!(g.within(1.0, 1.0, 1.0, 2), g.constant(1.0));
    }

    #[test]
    fn numbers_known_within_the_same_bounds_are_two_numbers() {
        // Two exact numbers can round to the same bounds: two identities,
        // their difference bounded, never folded to 0; one identity is one
        // number, less itself 0.
        use crate::eval::{Ctx, eval_interval};
        use crate::interval::I;
        let mut g = Graph::new();
        let (a, b) = (g.within(1.0, 2.0, 1.5, 1), g.within(1.0, 2.0, 1.5, 2));
        assert_eq!(g.within(1.0, 2.0, 1.5, 1), a);
        assert_ne!(a, b);
        let d = g.sub(a, b);
        assert_eq!(g.as_const(d), None);
        let t = crate::tape::Tape::new(&g, d);
        let p = I::point(0.0);
        assert_eq!(eval_interval(&t, &Ctx::new(0.0, vec![]), [p, p, p], None), I::new(-1.0, 1.0));
        let z = g.sub(a, a);
        assert_eq!(g.as_const(z), Some(0.0));
    }

    #[test]
    fn sealed_names_change_the_fingerprint() {
        // A sealed name (provenance stops there: a fillet's face reports the
        // fillet) means different face names, so results keyed by content
        // must not be shared with the same expression under a plain name.
        let print = |sealed: bool| {
            let mut g = Graph::new();
            let x = g.var(Var::X);
            let r = if sealed { g.tag_opaque("f", x) } else { g.tag("f", x) };
            let mut fp = Fingerprint::default();
            g.fingerprint(&[r], &[], true, &mut fp);
            fp.value()
        };
        assert_ne!(print(true), print(false));
    }

    /// Two operands of a sum that differ only in which param they read
    /// are ordered by where the question first reads each, not by where
    /// the graph made them: a model that declares or imports its params in
    /// another order asks the same question.
    #[test]
    fn a_sum_of_two_params_has_one_key_whatever_their_order() {
        let key = |q_first: bool| {
            let mut g = Graph::new();
            let (p, q) = if q_first {
                let q = g.var(Var::Param(1));
                (g.var(Var::Param(0)), q)
            } else {
                let p = g.var(Var::Param(0));
                (p, g.var(Var::Param(1)))
            };
            let (sp, sq) = if q_first {
                let sq = g.un(UnOp::Sin, q);
                (g.un(UnOp::Sin, p), sq)
            } else {
                let sp = g.un(UnOp::Sin, p);
                (sp, g.un(UnOp::Sin, q))
            };
            let sum = g.add(sp, sq);
            let root = g.sub(p, sum);
            let tape = crate::tape::Tape::new(&g, root);
            let ctx = crate::eval::Ctx::new(0.0, vec![0.25, 0.5]);
            let names = crate::tape::Names::of(&[&tape], &ctx);
            let mut fp = Fingerprint::default();
            tape.fingerprint(&mut fp, &names);
            fp.value()
        };
        assert_eq!(key(false), key(true));
    }

    /// Neither read before, two params' terms are ordered by the params'
    /// names, not by which the model declared or the graph made first: the
    /// question, its values included, has one key.
    #[test]
    fn a_sum_of_two_unread_params_has_one_key_whatever_their_order() {
        let key = |q_first: bool| {
            let mut g = Graph::new();
            // q_first: q is param 0, declared and built first.
            let (kp, kq) = if q_first { (1, 0) } else { (0, 1) };
            g.name_param(kp, "p");
            g.name_param(kq, "q");
            let (sp, sq) = if q_first {
                let q = g.var(Var::Param(kq));
                let sq = g.un(UnOp::Sin, q);
                let p = g.var(Var::Param(kp));
                (g.un(UnOp::Sin, p), sq)
            } else {
                let p = g.var(Var::Param(kp));
                let sp = g.un(UnOp::Sin, p);
                let q = g.var(Var::Param(kq));
                (sp, g.un(UnOp::Sin, q))
            };
            let root = g.add(sp, sq);
            let tape = crate::tape::Tape::new(&g, root);
            let mut values = vec![0.0; 2];
            (values[kp as usize], values[kq as usize]) = (0.25, 0.5);
            let ctx = crate::eval::Ctx::new(0.0, values);
            let names = crate::tape::Names::of(&[&tape], &ctx);
            let mut fp = Fingerprint::default();
            tape.fingerprint(&mut fp, &names);
            assert!(ctx.fingerprint(&mut fp, &names));
            fp.value()
        };
        assert_eq!(key(false), key(true));
    }

    /// Terms alike but for their params are ordered by their own terms in
    /// order, inner ties first: x - (sin (p + zz) + sin (q + r)) has one
    /// key whichever of p and zz the graph made first (an unrelated
    /// `(def unused zz)` makes zz first).
    #[test]
    fn nested_ties_are_ordered_from_the_inside_out() {
        let key = |zz_first: bool| {
            let mut g = Graph::new();
            for (k, n) in ["p", "zz", "q", "r"].iter().enumerate() {
                g.name_param(k as u32, n);
            }
            let (p, zz) = if zz_first {
                let zz = g.var(Var::Param(1));
                (g.var(Var::Param(0)), zz)
            } else {
                let p = g.var(Var::Param(0));
                (p, g.var(Var::Param(1)))
            };
            let (q, r) = (g.var(Var::Param(2)), g.var(Var::Param(3)));
            let pz = g.add(p, zz);
            let qr = g.add(q, r);
            let (a, b) = (g.un(UnOp::Sin, pz), g.un(UnOp::Sin, qr));
            let sum = g.add(a, b);
            let x = g.var(Var::X);
            let root = g.sub(x, sum);
            let tape = crate::tape::Tape::new(&g, root);
            let ctx = crate::eval::Ctx::new(0.0, vec![0.1, 0.2, 0.3, 0.4]);
            let names = crate::tape::Names::of(&[&tape], &ctx);
            let mut fp = Fingerprint::default();
            tape.fingerprint(&mut fp, &names);
            assert!(ctx.fingerprint(&mut fp, &names));
            fp.value()
        };
        assert_eq!(key(false), key(true));
    }

    /// Library functions give musl's bits, which no CPU dispatch changes.
    /// glibc's differ in the last bit at these inputs (by CPU: sin at p,
    /// where 1e16 (sin p + 0.3775604968083285) is 0 or -0.56, so a model
    /// would be proven empty on one machine and disproven on another, and
    /// the result cache would hand one's answer to the other).
    #[test]
    fn library_functions_are_the_same_on_every_machine() {
        let p = -21.60398819379739;
        assert_eq!(un_f64(UnOp::Sin, p), -0.3775604968083285);
        assert_eq!(un_f64(UnOp::Sin, -1499.931286), 0.9839855162426463);
        assert_eq!(un_f64(UnOp::Cos, -1499.877923), -0.23047820418116294);
        assert_eq!(un_f64(UnOp::Tan, -1499.88304), -4.320450143326324);
        assert_eq!(bin_f64(BinOp::Atan2, -1499.548242, 1.0), -1.5701294593852633);
        assert_eq!(un_f64(UnOp::Exp, -2.142849832857143), 0.11732002370048947);
        assert_eq!(un_f64(UnOp::Ln, 1499.959064), 7.313193096051238);
        // The derivatives' values, the constants a graph folds and the
        // range rules' ends are the same functions.
        assert_eq!(dun(UnOp::Sin, -1499.931286).0, 0.9839855162426463);
        let mut g = Graph::default();
        let c = g.constant(-1499.931286);
        let s = g.un(UnOp::Sin, c);
        assert_eq!(g.as_const(s), Some(0.9839855162426463));
        let r = crate::interval::un(UnOp::Exp, crate::interval::I::point(-2.142849832857143));
        assert_eq!(r, crate::interval::widen(crate::interval::I::point(0.11732002370048947), 4));
    }

    #[test]
    fn a_rewrite_keeps_a_promise_only_where_it_still_holds() {
        // (lipschitz 1 (- x p)) holds for every p; with p warped to 1000 x
        // it is -999 x: kept, the promise would prove a clearance of 2
        // across a 0.9 gap.
        use crate::eval::{Ctx, eval_interval_grad};
        use crate::interval::{I, lipschitz};
        let slope = |g: &Graph, f: NodeId| {
            let t = crate::tape::Tape::new(g, f);
            lipschitz(&eval_interval_grad(&t, &Ctx::new(0.0, vec![0.5]), [I::new(-1.0, 1.0); 3]))
        };
        let mut g = Graph::new();
        let (x, p, t) = (g.var(Var::X), g.var(Var::Param(0)), g.var(Var::T));
        let (one, k) = (g.constant(1.0), g.constant(1000.0));
        let d = g.sub(x, p);
        let f = g.lipschitz(d, one);
        let kx = g.mul(k, x);
        let warped = g.subst(f, &[(Var::Param(0), kx)]);
        let replaced = g.replace(f, &[(p, kx)].into());
        assert!(slope(&g, warped) >= 999.0 && slope(&g, replaced) >= 999.0);
        // x (cos^2 p + sin^2 p) is x; with cos^2 p replaced by 1000, a
        // non-spatial node, it is not 1-Lipschitz.
        let (c, s) = (g.un(UnOp::Cos, p), g.un(UnOp::Sin, p));
        let (c2, s2) = (g.un(UnOp::Square, c), g.un(UnOp::Square, s));
        let sum = g.add(c2, s2);
        let e = g.mul(x, sum);
        let f2 = g.lipschitz(e, one);
        let swapped = g.replace(f2, &[(c2, k)].into());
        assert!(slope(&g, swapped) >= 1000.0);
        // Kept: x warped (the frame stretches with it), p to an expression
        // of t.
        let two = g.constant(2.0);
        let (x2, t2) = (g.mul(two, x), g.mul(two, t));
        for w in [g.subst(f, &[(Var::X, x2)]), g.subst(f, &[(Var::Param(0), t2)]), g.replace(f, &[(x, x2)].into())] {
            assert!(matches!(g.nodes[w as usize], Node::Lipschitz(..)));
        }
    }
}

/// The splitmix64 finaliser: a bijective mix of 64 bits.
pub const fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

/// A bijective mix of 128 bits, each output bit depending on every input
/// bit: spreads FNV's weak low bits.
const fn mix128(z: u128) -> u128 {
    let (lo, hi) = (z as u64, (z >> 64) as u64);
    let hi = mix64(hi ^ mix64(lo));
    let lo = mix64(lo ^ hi);
    (hi as u128) << 64 | lo as u128
}

/// A deterministic 128-bit hash (FNV-1a over bytes, then a final mix), for
/// content keys that must be stable across runs: at 128 bits, two contents
/// never share a key by chance.
#[derive(Clone, Copy, Debug)]
pub struct Fingerprint(u128);

impl Default for Fingerprint {
    fn default() -> Self {
        Fingerprint(0x6c62272e07bb014262b821756295c58d)
    }
}

impl Fingerprint {
    pub fn word(&mut self, w: u64) -> &mut Self {
        self.raw(&w.to_le_bytes())
    }
    /// Another hash (a node's, a key).
    pub fn key(&mut self, k: u128) -> &mut Self {
        self.raw(&k.to_le_bytes())
    }
    pub fn f64(&mut self, v: f64) -> &mut Self {
        self.word(v.to_bits())
    }
    pub fn str(&mut self, s: &str) -> &mut Self {
        self.word(s.len() as u64).raw(s.as_bytes())
    }
    /// Bytes as they are: the caller frames them (a length before them).
    pub fn raw(&mut self, bytes: &[u8]) -> &mut Self {
        for b in bytes {
            self.0 = (self.0 ^ *b as u128).wrapping_mul(0x0000000001000000000000000000013b);
        }
        self
    }
    pub fn value(&self) -> u128 {
        mix128(self.0)
    }
}
