//! Hash-consed expression DAG with constant folding.
//!
//! Node ids are topologically ordered: every node's children have smaller ids.
//! Equal math is stored once.

mod build;
mod differentiate;
mod exact;
mod fingerprint;
mod grid;
mod motion;
mod node;
mod rewrite;
mod scalar;
mod traverse;

#[cfg(test)]
mod tests;

pub use fingerprint::{Fingerprint, mix64};
pub use motion::{DriveRead, MotionRead};
pub use node::{BinOp, Node, NodeId, UnOp, Var};
pub use scalar::{
    atan2_f64, bin_f64, dbin, defined, div_f64, dun, emax_f64, exact_bin, exact_un, ln_f64, lost, pow_f64, un_f64, undefined_bin, undefined_bin_over, undefined_un,
    undefined_un_over, wild,
};

use build::Key;
// These maps index graph-owned keys; node ids follow insertion order,
// and traversals sort their ids before rebuilding. Hashing changes neither.
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

#[derive(Default)]
pub struct Graph {
    pub nodes: Vec<Node>,
    index: HashMap<Key, NodeId>,
    domains: motion::Domains,
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
    /// Per node: an operation in it has operands it gives 0 for, having no
    /// result there (a quotient, a power, ln; and sqrt, which is 0 below
    /// 0). While `strict`, a number made of constants and params that is
    /// one keeps its form (`hides`): no fold or identity turns it into a
    /// plain 0, so what reads the number can ask whether it is one
    /// (`Graph::may_lack_a_result`).
    partial: Vec<bool>,
    /// Per node: it reads a marker (`marker`): a param whose being read
    /// is what matters, not its number (an estimate's, 0 where it applies).
    /// No fold or identity drops one (`hides`): 0 times it and it less
    /// itself keep their form, so what reads them still reads the marker.
    marked: Vec<bool>,
    /// The params that are markers (`marker`).
    markers: HashSet<u32>,
    /// Set while a number is built that must not hide an operation with
    /// no result (a stated fact's): 0 times such a number, such a number
    /// less itself and the sqrt of a negative constant are then kept as
    /// written, where otherwise they are 0 (a zero-length segment's 0 / 0
    /// times its length is discarded: a shape's rule, kept).
    pub strict: bool,
    /// Real expressions behind precomputed constants, for symbolic identities.
    constructions: HashMap<NodeId, Node>,
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
    /// Operand pairs `(a, b)` for `Node::CornerReach(ci, r, guard, at)`.
    pub corners: Vec<(NodeId, NodeId)>,
    corner_index: HashMap<(NodeId, NodeId), u32>,
    /// The first number outside ink's range (docs/CORE.md Answers) the
    /// model is built from (`note`), its bits (0: none).
    wild: std::sync::atomic::AtomicU64,
    /// std's own constants outside the range (`own`), and folds of them.
    owned: HashSet<NodeId>,
    /// Constants of the arithmetic on a term that read x, y, z or t, which
    /// cancelled (`zero_of`) or folded exactly with one (`exact_number`):
    /// each a constant to every reader, kept apart from the same number
    /// written, so a fold of one is exact or a band, never a rounding.
    exact: HashSet<NodeId>,
    /// The real value of each constant that is not its float (`real_constant`): a literal's exact fraction,
    /// pi's interval, what folding them gave. Each such constant is a node of its own, never one shared
    /// with the same float made otherwise. Only a value claim on constants reads it.
    reals: HashMap<NodeId, crate::exact::Exact>,
}
