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
    /// interpolant at coordinates (x, y, z), clamped into the grid
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
    /// Exact ball-clearance corner test (`corners[ci]` = `(a, b)`), radius `r`,
    /// active where `guard < 0`, evaluated at `at` (a `Frame([fx, fy, fz])`).
    /// Holds `a` and `b` by reference rather than inlining gradients and warped
    /// copies, so nested fillets do not multiply the tape.
    CornerReach(u32, NodeId, NodeId, NodeId),
}

impl Node {
    pub fn children(&self) -> ([NodeId; 3], usize) {
        match *self {
            Node::Const(_) | Node::Var(_) => ([0; 3], 0),
            Node::Un(_, a) => ([a, 0, 0], 1),
            Node::Bin(_, a, b) => ([a, b, 0], 2),
            Node::Select(c, a, b) | Node::Lipschitz(c, a, b) | Node::Within(c, a, b, _) | Node::CornerReach(_, c, a, b) => ([c, a, b], 3),
            Node::Tag(_, a) | Node::Real(a) => ([a, 0, 0], 1),
            Node::Read(_, _, a, at) => ([a, at, 0], 2),
            Node::Deriv(a, f, _) => ([a, f, 0], 2),
            Node::Grid(_, _, c) | Node::Frame(c) => (c, 3),
        }
    }
}
