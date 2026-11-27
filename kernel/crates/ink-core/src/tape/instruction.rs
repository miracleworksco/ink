//! Tape opcodes, their operands, and the cost of evaluating them.

use crate::graph::{BinOp, UnOp, Var};

#[derive(Clone, Copy, Debug)]
pub enum Op {
    Const(f64),
    Var(Var),
    Un(UnOp),
    Bin(BinOp),
    /// a = condition, b = value where a < 0, c = value otherwise.
    Select,
    /// Identity carrying a provenance name (see graph::Node::Tag).
    Tag(u32),
    /// Component `mask` of sampled field `Tape::grids[i]` at coordinates
    /// (a, b, c) (see graph::Node::Grid).
    Grid(u32, u8),
    /// A sampled field read at the frame b (see graph::Node::Read): a. c is
    /// not an argument: it indexes
    /// `Tape::forms`, the tape of the read's form over x, y, z, whose
    /// gradient at the frame's coordinates is the read's partials (the
    /// reads one axis further), which the point derivatives chain (the
    /// frame's coordinates stay live until the read for them).
    Read,
    /// a, promised k-Lipschitz (k = b) in the frame c (see graph::Node::Lipschitz).
    Lipschitz,
    /// The coordinates (a, b, c) of a Lipschitz promise. Value 0.
    Frame,
    /// a, the derivative of b along an axis (see graph::Node::Deriv): its
    /// range is also within b's gradient bounds along x, y or z (0, 1, 2);
    /// along t (3), a mark only.
    Deriv(u8),
    /// a, read as a real number (see graph::Node::Real): never fixed.
    Real,
    /// A number known only within the constants a..b, c at a point, of
    /// that identity (see graph::Node::Within): never fixed.
    Within(u64),
    /// Exact ball-clearance corner test (`forms[fa]`, `forms[fb]`): `a` = radius `r`,
    /// `b` = frame `at` (`Op::Frame`), `c` = `guard` (active where `< 0`, else 0).
    CornerReach(u32, u32),
}

#[derive(Clone, Copy, Debug)]
pub struct Ins {
    pub op: Op,
    pub out: u32,
    pub a: u32,
    pub b: u32,
    pub c: u32,
}

/// How a tape is evaluated, for what that costs (`Op::work`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pass {
    /// At points, in lanes (eval::eval_points).
    Point,
    /// Its value and gradient at points (eval::eval_grads).
    Gradient,
    /// Its range over a box (eval::eval_interval).
    Box,
}

impl Op {
    /// The machine instructions evaluating it retires, by `Pass`: its own,
    /// beyond its operands', counted on a core of a Zen 5 with AVX-512 (a
    /// point's share of the lanes it runs in). Arithmetic costs about one
    /// at a point, a transcendental 15 to 270 times that, a sampled
    /// field's read about 800; over a box, 60 to 1,000.
    pub fn work(self, pass: Pass) -> usize {
        let [point, gradient, range] = match self {
            Op::Const(_) | Op::Var(_) => [1, 50, 60],
            Op::Bin(BinOp::Mul) => [1, 117, 300],
            Op::Bin(BinOp::Div) => [2, 130, 255],
            Op::Un(UnOp::Floor) => [1, 65, 275],
            Op::Un(UnOp::Rad | UnOp::Deg) => [15, 65, 350],
            Op::Bin(BinOp::Emax) => [70, 175, 520],
            Op::Un(UnOp::Exp | UnOp::Ln) => [78, 130, 375],
            Op::Un(UnOp::Sin | UnOp::Cos) => [93, 218, 490],
            Op::Un(UnOp::Tan) => [121, 171, 490],
            Op::Bin(BinOp::Atan2) => [120, 235, 690],
            Op::Bin(BinOp::Pow) => [273, 775, 980],
            Op::Grid(..) | Op::Read => [780, 8300, 6300],
            Op::CornerReach(..) => [100, 150, 300],
            _ => [1, 70, 120],
        };
        match pass {
            Pass::Point => point,
            Pass::Gradient => gradient,
            Pass::Box => range,
        }
    }
}

impl Ins {
    /// Argument slots (or ssa indices) and how many there are.
    pub fn args(&self) -> ([u32; 3], usize) {
        match self.op {
            Op::Const(_) | Op::Var(_) => ([0; 3], 0),
            Op::Un(_) | Op::Tag(_) | Op::Real => ([self.a, 0, 0], 1),
            Op::Bin(_) | Op::Read | Op::Deriv(_) => ([self.a, self.b, 0], 2),
            Op::Select | Op::Grid(..) | Op::Lipschitz | Op::Frame | Op::Within(_) | Op::CornerReach(..) => ([self.a, self.b, self.c], 3),
        }
    }
}
