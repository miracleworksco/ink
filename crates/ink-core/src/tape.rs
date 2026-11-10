//! Flat instruction tapes compiled from the DAG.
//!
//! A tape keeps two views of the same program:
//! - `ssa`: one value per instruction, arguments are instruction indices.
//!   Used for simplification.
//! - `code`: the same instructions with register-allocated slots.
//!   Used by the evaluators, so scratch memory stays small.
//!
//! The last instruction is the output.

use crate::graph::{BinOp, Fingerprint, Graph, Node, NodeId, UnOp, Var};
use crate::interval::I;

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
}

#[derive(Clone, Copy, Debug)]
pub struct Ins {
    pub op: Op,
    pub out: u32,
    pub a: u32,
    pub b: u32,
    pub c: u32,
}

impl Ins {
    /// Argument slots (or ssa indices) and how many there are.
    pub fn args(&self) -> ([u32; 3], usize) {
        match self.op {
            Op::Const(_) | Op::Var(_) => ([0; 3], 0),
            Op::Un(_) | Op::Tag(_) | Op::Real => ([self.a, 0, 0], 1),
            Op::Bin(_) | Op::Read | Op::Deriv(_) => ([self.a, self.b, 0], 2),
            Op::Select | Op::Grid(..) | Op::Lipschitz | Op::Frame | Op::Within(_) => ([self.a, self.b, self.c], 3),
        }
    }
}

/// The first instruction from `j` on that every later use of j's value
/// passes through (no other value depending on j is still to be read),
/// and from which the output is `monotone`. `dep` is scratch, all false.
fn join_of(ssa: &[Ins], last: &[u32], monotone: &[bool], dep: &mut [bool], j: usize) -> usize {
    let n = ssa.len();
    if monotone[j] {
        return j;
    }
    dep[j] = true;
    // Values that depend on j and are still to be read.
    let mut live = 1usize;
    let mut at = n - 1;
    for i in j + 1..n {
        let (args, k) = Tape::reads(ssa, i);
        let mut d = false;
        for (m, &a) in args[..k].iter().enumerate() {
            let a = a as usize;
            if dep[a] && !args[..m].contains(&args[m]) {
                d = true;
                if last[a] as usize == i {
                    live -= 1;
                }
            }
        }
        if d {
            dep[i] = true;
            live += 1;
            if live == 1 && monotone[i] {
                at = i;
                break;
            }
        }
    }
    dep[j..=at].fill(false);
    at
}

/// Tape::forked_steps over instructions `from..to`.
fn forked_steps(jumps: &[u32], rejoin: &[u32], from: usize, to: usize) -> u64 {
    let (mut steps, mut i) = (0u64, from);
    for &j in &jumps[jumps.partition_point(|&j| (j as usize) < from)..] {
        let (j, r) = (j as usize, rejoin[j as usize] as usize);
        if j >= to {
            break;
        }
        if j >= i && r > j && r < to {
            steps = steps.saturating_add((j + 1 - i) as u64).saturating_add(forked_steps(jumps, rejoin, j + 1, r + 1).saturating_mul(2));
            i = r + 1;
        }
    }
    steps.saturating_add((to - i) as u64)
}

/// Whether instruction `u`'s bound is its operand `a`'s ends mapped, with
/// its other operands fixed (so a bound of it from the hull of two bounds
/// of `a` is the hull of the bounds from each): sums, min, max, emax, a
/// select's branches, a monotone function, a product with a constant; `a`
/// read once.
fn ends_mapped(ssa: &[Ins], u: usize, a: u32) -> bool {
    let ins = &ssa[u];
    let (args, k) = ins.args();
    if args[..k].iter().filter(|&&x| x == a).count() != 1 {
        return false;
    }
    let constant = |x: u32| matches!(ssa[x as usize].op, Op::Const(c) if !c.is_nan());
    match ins.op {
        Op::Bin(BinOp::Add | BinOp::Sub | BinOp::Min | BinOp::Max | BinOp::Emax) => true,
        Op::Bin(BinOp::Mul) => constant(if ins.a == a { ins.b } else { ins.a }),
        Op::Bin(BinOp::Div) => ins.a == a && constant(ins.b),
        Op::Un(UnOp::Neg | UnOp::Sqrt | UnOp::Exp | UnOp::Floor) | Op::Tag(_) | Op::Real => true,
        Op::Select => ins.a != a,
        Op::Lipschitz => ins.a == a,
        _ => false,
    }
}

/// Choice codes recorded by interval evaluation for min/max/select.
pub const BOTH: u8 = 0;
pub const FIRST: u8 = 1;
pub const SECOND: u8 = 2;

#[derive(Clone, Debug)]
pub struct Tape {
    pub ssa: Vec<Ins>,
    pub code: Vec<Ins>,
    pub nslots: usize,
    /// Sampled fields referenced by Op::Grid (the graph's, shared).
    pub grids: Vec<std::sync::Arc<crate::grid::GridData>>,
    /// Tag names provenance stops at (graph::Graph::opaque).
    pub opaque: Vec<u32>,
    /// Per instruction: its value depends only on constants and params (in
    /// the model, not just in a simplified tape, where a picked branch can
    /// leave spatial arithmetic with constant operands). Params are fixed for
    /// a query, so interval evaluation takes these as exact points computed
    /// in floats, as load-time defs are; unless one reads a ranged param
    /// (eval::ParamSpans), whose value it means in real numbers.
    pub fixed: Vec<bool>,
    /// Per fixed instruction: bit k % 64 of each param k it reads (in the
    /// model, as `fixed`): a prefilter, where a context's ranged params'
    /// bits are not set it reads none of them (two params 64 apart share a
    /// bit).
    pub pmask: Vec<u64>,
    /// Per fixed instruction that reads params, which ones, exactly (in the
    /// model, as `fixed`): (instruction, start, len) by instruction, the
    /// params sorted in `plist[start..start + len]` (`params_of`), a list
    /// shared where it is an argument's. Read only where `pmask` meets a
    /// context's ranged params: a ranged param this reads takes the
    /// instruction off the float path, one sharing only its bit does not.
    pub pread: Vec<(u32, u32, u32)>,
    pub plist: Vec<u32>,
    /// Per instruction: whether a derivative's range reads its gradient
    /// bounds (it is in the cone of a field an Op::Deriv bounds by). Empty
    /// where none does.
    pub slopes: Vec<bool>,
    /// The tapes of the read forms Op::Read names (graph::Graph::read_form).
    pub forms: Vec<std::sync::Arc<Tape>>,
    /// The instructions that can jump (floor, atan2, a select whose
    /// condition reads no space), in order, which interval evaluation may
    /// take as a choice of two sides (eval::forked). Derived from `ssa`.
    pub jumps: Vec<u32>,
    /// Where jumps join (Tape::rejoin), derived from `ssa` when a jump is
    /// first straddled.
    joins: std::sync::OnceLock<Joins>,
    /// Each promise's operands (Tape::promise) by its instruction, each
    /// found when first read.
    promises: std::sync::OnceLock<Box<[(u32, std::sync::OnceLock<Promise>)]>>,
    /// Per instruction, bit k of each spatial axis k (x, y, z) its value
    /// reads (Tape::axes_read), derived from `ssa` when first read.
    axes: std::sync::OnceLock<Vec<u8>>,
    /// Where pieces may meet (crate::sign::levels), derived from `ssa` when
    /// first read.
    pub(crate) levels: std::sync::OnceLock<Vec<(u32, f64)>>,
    /// Per min or max of a and a + c (or a - c), through names, in order:
    /// the instruction, the slot holding c there (kept until it), and
    /// whether the second operand less the first is -c (`Tape::apart`).
    /// Derived from `ssa`.
    pub apart: Vec<(u32, u32, bool)>,
}

/// What interval evaluation reads again of a promise, Lipschitz instruction
/// `i` of the SSA form (eval::promised, eval::promise_valid): the
/// instructions its field, its constant and its frame read, in order,
/// renumbered among themselves; which are a select's condition's (`cond`)
/// and which the field's and the frame's coordinates' (`at`, bounded at its
/// centre); and k's and the coordinates' alone (`over`, bounded over the
/// box), renumbered among themselves.
#[derive(Clone, Debug)]
pub struct Promise {
    /// Each instruction's index in the tape, and it with its arguments
    /// renumbered.
    pub ins: Vec<(u32, Ins)>,
    pub cond: Vec<bool>,
    pub at: Vec<bool>,
    pub over: Vec<(u32, Ins)>,
    /// The promise, its arguments renumbered; where the field, the
    /// constant and the frame's coordinates are among them.
    pub promise: Ins,
    pub f: u32,
    pub frame: [u32; 3],
    /// Where k and the coordinates are in `over`.
    pub over_k: u32,
    pub over_frame: [u32; 3],
}

impl Promise {
    fn of(ssa: &[Ins], i: usize) -> Promise {
        let (l, fr) = (ssa[i], ssa[ssa[i].c as usize]);
        let frame = [fr.a, fr.b, fr.c];
        let cone = |roots: &[u32]| {
            let mut seen = std::collections::HashSet::new();
            let mut stack = roots.to_vec();
            while let Some(j) = stack.pop() {
                if seen.insert(j) {
                    let (args, k) = ssa[j as usize].args();
                    stack.extend_from_slice(&args[..k]);
                }
            }
            seen
        };
        let all = cone(&[l.a, l.b, l.c]);
        let mut order: Vec<u32> = all.iter().copied().collect();
        order.sort_unstable();
        let conditions: Vec<u32> = order.iter().filter(|&&j| matches!(ssa[j as usize].op, Op::Select)).map(|&j| ssa[j as usize].a).collect();
        let (cond, over, at) = (cone(&conditions), cone(&[l.b, frame[0], frame[1], frame[2]]), cone(&[l.a, frame[0], frame[1], frame[2]]));
        let renumber = |order: &[u32]| {
            let pos: std::collections::HashMap<u32, u32> = order.iter().enumerate().map(|(k, &j)| (j, k as u32)).collect();
            let ins = order
                .iter()
                .map(|&j| {
                    let mut c = ssa[j as usize];
                    let (_, n) = c.args();
                    for a in [&mut c.a, &mut c.b, &mut c.c].into_iter().take(n) {
                        *a = pos[&*a];
                    }
                    (j, c)
                })
                .collect::<Vec<_>>();
            (ins, pos)
        };
        let (ins, pos) = renumber(&order);
        let mut over_order: Vec<u32> = over.iter().copied().collect();
        over_order.sort_unstable();
        let (over_ins, over_pos) = renumber(&over_order);
        let marks = |set: &std::collections::HashSet<u32>| order.iter().map(|j| set.contains(j)).collect();
        Promise {
            ins,
            cond: marks(&cond),
            at: marks(&at),
            over: over_ins,
            promise: Ins { a: pos[&l.a], b: pos[&l.b], c: pos[&l.c], ..l },
            f: pos[&l.a],
            frame: frame.map(|c| pos[&c]),
            over_k: over_pos[&l.b],
            over_frame: frame.map(|c| over_pos[&c]),
        }
    }
}

/// Where each of a tape's jumps joins (Tape::rejoin). An instruction is
/// monotone where its value reaches the output only through operations
/// that each read it once and whose bounds are their operand's ends mapped
/// (sums, min, max, a monotone function): a bound of the output from the
/// hull of two bounds of it is then the hull of the bounds from each, so
/// nothing is lost by joining them there (eval::forked).
#[derive(Clone, Debug)]
struct Joins {
    /// Per jump (Tape::jumps): Tape::rejoin's instruction; 0 for others.
    rejoin: Vec<u32>,
    /// Tape::forked_steps.
    forked: u64,
}

/// What a region holds besides what a mesh draws there (Tape::sheet): one
/// sheet, nothing else; what is thinner than this (along the direction
/// checked); or unknown.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Sheet {
    One,
    Thin(f64),
    Unknown,
}

/// `n`'s operands in the order its tape visits them: a commutative
/// operation's whose expressions differ only in their names
/// (Graph::tied) by `named_shape`, so the tape, and the question's
/// numbering of its params, do not depend on which the graph made first,
/// nor on where the model declared them. Else as the node holds them.
fn ordered(g: &Graph, n: NodeId, memo: &mut std::collections::HashMap<NodeId, u128>) -> [NodeId; 3] {
    let (ch, _) = g.nodes[n as usize].children();
    if g.tied(n) && named_shape(g, ch[1], memo) < named_shape(g, ch[0], memo) { [ch[1], ch[0], 0] } else { ch }
}

/// A hash of `n`'s expression with the names in it (a param's
/// `Graph::param_names`, or its number where it has none; a tag's), a
/// commutative operation's operands taken in either order alike: found
/// from the inside out, so it is the same however the graph was built.
/// Memoised in `memo`.
fn named_shape(g: &Graph, n: NodeId, memo: &mut std::collections::HashMap<NodeId, u128>) -> u128 {
    let mut stack = vec![(n, false)];
    while let Some((m, operands_done)) = stack.pop() {
        if memo.contains_key(&m) {
            continue;
        }
        let node = g.nodes[m as usize];
        let (cs, k) = node.children();
        if !operands_done {
            stack.push((m, true));
            stack.extend(cs[..k].iter().filter(|c| !memo.contains_key(c)).map(|&c| (c, false)));
            continue;
        }
        let h = |c: NodeId| memo[&c];
        let mut f = Fingerprint::default();
        match node {
            Node::Var(Var::Param(p)) => match g.param_names.get(p as usize).filter(|s| !s.is_empty()) {
                Some(name) => f.word(1).str(name),
                None => f.word(2).word(p as u64),
            },
            Node::Tag(t, a) => f.word(3).str(&g.tag_names[t as usize]).key(h(a)),
            Node::Bin(op @ (BinOp::Add | BinOp::Mul | BinOp::Min | BinOp::Max | BinOp::Emax), a, b) => {
                let (x, y) = (h(a).min(h(b)), h(a).max(h(b)));
                f.word(4).word(op as u64).key(x).key(y)
            }
            _ => {
                f.word(5).key(g.shape(m));
                cs[..k].iter().for_each(|&c| _ = f.key(h(c)));
                &mut f
            }
        };
        memo.insert(m, f.value());
    }
    memo[&n]
}

impl Tape {
    /// The tape of `root`: its instructions in one order, operands before
    /// their consumers, depth first from the root, each node's operands in
    /// turn (a commutative operation's alike up to names, by those names:
    /// `ordered`); with the grids it reads (renumbered in that order) and the
    /// names provenance stops at among its own. So the tape of an
    /// expression is the same however the graph holding it was built, and
    /// whatever else it holds (Tape::fingerprint).
    pub fn new(g: &Graph, root: NodeId) -> Tape {
        let mut idx: std::collections::HashMap<NodeId, u32> = std::collections::HashMap::new();
        let mut ids: Vec<NodeId> = Vec::new();
        // Each tied operand's named_shape, found once.
        let mut memo: std::collections::HashMap<NodeId, u128> = std::collections::HashMap::new();
        // Commutative operations whose operands are visited swapped.
        let mut swapped: std::collections::HashSet<NodeId> = std::collections::HashSet::new();
        // Per node on the stack: its operands in the order visited, and how
        // many have been.
        let mut stack: Vec<(NodeId, [NodeId; 3], usize)> = vec![(root, ordered(g, root, &mut memo), 0)];
        while let Some((n, ch, k)) = stack.last_mut() {
            let count = g.nodes[*n as usize].children().1;
            if *k < count {
                let c = ch[*k];
                *k += 1;
                if !idx.contains_key(&c) {
                    stack.push((c, ordered(g, c, &mut memo), 0));
                }
            } else {
                let (n, ch) = (*n, *ch);
                stack.pop();
                if ch != g.nodes[n as usize].children().0 {
                    swapped.insert(n);
                }
                if let std::collections::hash_map::Entry::Vacant(e) = idx.entry(n) {
                    e.insert(ids.len() as u32);
                    ids.push(n);
                }
            }
        }
        let mut ssa = Vec::with_capacity(ids.len());
        let mut forms: Vec<std::sync::Arc<Tape>> = Vec::new();
        let mut form_index: std::collections::HashMap<NodeId, u32> = std::collections::HashMap::new();
        let mut grids: Vec<u32> = Vec::new();
        let mut grid = |gi: u32| match grids.iter().position(|&x| x == gi) {
            Some(k) => k as u32,
            None => {
                grids.push(gi);
                (grids.len() - 1) as u32
            }
        };
        let mut tags: Vec<u32> = Vec::new();
        for (i, &id) in ids.iter().enumerate() {
            let m = |n: NodeId| idx[&n];
            let out = i as u32;
            let ins = match g.nodes[id as usize] {
                Node::Const(v) => Ins { op: Op::Const(v), out, a: 0, b: 0, c: 0 },
                Node::Var(v) => Ins { op: Op::Var(v), out, a: 0, b: 0, c: 0 },
                Node::Un(o, a) => Ins { op: Op::Un(o), out, a: m(a), b: 0, c: 0 },
                // A commutative operation's operands in the order visited.
                Node::Bin(o, a, b) if swapped.contains(&id) => Ins { op: Op::Bin(o), out, a: m(b), b: m(a), c: 0 },
                Node::Bin(o, a, b) => Ins { op: Op::Bin(o), out, a: m(a), b: m(b), c: 0 },
                Node::Select(c, a, b) => Ins { op: Op::Select, out, a: m(c), b: m(a), c: m(b) },
                Node::Tag(t, a) => {
                    tags.push(t);
                    Ins { op: Op::Tag(t), out, a: m(a), b: 0, c: 0 }
                }
                Node::Grid(gi, mask, c) => Ins { op: Op::Grid(grid(gi), mask), out, a: m(c[0]), b: m(c[1]), c: m(c[2]) },
                Node::Read(gi, path, a, at) => {
                    let form = g.read_form_of(gi, path);
                    let f = *form_index.entry(form).or_insert_with(|| {
                        forms.push(std::sync::Arc::new(Tape::new(g, form)));
                        (forms.len() - 1) as u32
                    });
                    Ins { op: Op::Read, out, a: m(a), b: m(at), c: f }
                }
                Node::Lipschitz(a, k, f) => Ins { op: Op::Lipschitz, out, a: m(a), b: m(k), c: m(f) },
                Node::Frame(c) => Ins { op: Op::Frame, out, a: m(c[0]), b: m(c[1]), c: m(c[2]) },
                Node::Deriv(a, f, k) => Ins { op: Op::Deriv(k), out, a: m(a), b: m(f), c: 0 },
                Node::Real(a) => Ins { op: Op::Real, out, a: m(a), b: 0, c: 0 },
                Node::Within(l, h, v, id) => Ins { op: Op::Within(id), out, a: m(l), b: m(h), c: m(v) },
            };
            ssa.push(ins);
        }
        let mut t = Self::from_ssa(ssa);
        t.forms = forms;
        t.grids = grids.iter().map(|&gi| g.grids[gi as usize].clone()).collect();
        t.opaque = g.opaque.iter().copied().filter(|o| tags.contains(o)).collect();
        t
    }

    /// What instruction i of `ssa` needs live: its arguments, and a read
    /// its frame's coordinates too (Op::Read).
    pub fn reads(ssa: &[Ins], i: usize) -> ([u32; 5], usize) {
        let (args, k) = ssa[i].args();
        let mut u = [0; 5];
        u[..k].copy_from_slice(&args[..k]);
        if let Op::Read = ssa[i].op {
            let (c, _) = ssa[ssa[i].b as usize].args();
            u[2..].copy_from_slice(&c);
            return (u, 5);
        }
        (u, k)
    }

    /// Of a min or max instruction `i`, whose operands are a and a + c or
    /// a - c, through names (either way round): c, and whether the second
    /// operand less the first is -c. Their order is c's sign, so bounds of
    /// c decide it where bounds of the two, which each read a, cannot (a
    /// smooth blend below the min of what it blends).
    pub fn apart_of(ssa: &[Ins], i: usize) -> Option<(u32, bool)> {
        if !matches!(ssa[i].op, Op::Bin(BinOp::Min | BinOp::Max)) {
            return None;
        }
        let strip = |mut j: u32| {
            while let Op::Tag(_) = ssa[j as usize].op {
                j = ssa[j as usize].a;
            }
            j
        };
        let (a, b) = (strip(ssa[i].a), strip(ssa[i].b));
        // y as x + c (false) or x - c (true): c and the sign.
        let off = |y: u32, x: u32| {
            let ins = ssa[y as usize];
            match ins.op {
                Op::Bin(BinOp::Sub) if strip(ins.a) == x => Some((ins.b, true)),
                Op::Bin(BinOp::Add) if strip(ins.a) == x => Some((ins.b, false)),
                Op::Bin(BinOp::Add) if strip(ins.b) == x => Some((ins.a, false)),
                _ => None,
            }
        };
        if a == b {
            return None;
        }
        off(b, a).or_else(|| off(a, b).map(|(c, neg)| (c, !neg)))
    }

    fn from_ssa(ssa: Vec<Ins>) -> Tape {
        let n = ssa.len();
        // A min or max of a and a + c reads c too (Tape::apart): its slot
        // is kept until there.
        let uses = |i: usize| {
            let (mut u, k) = Self::reads(&ssa, i);
            match Self::apart_of(&ssa, i) {
                Some((c, _)) => {
                    u[k] = c;
                    (u, k + 1)
                }
                None => (u, k),
            }
        };
        let mut last = vec![0usize; n];
        for i in 0..n {
            let (args, k) = uses(i);
            for &a in &args[..k] {
                last[a as usize] = i;
            }
        }
        last[n - 1] = usize::MAX;

        let mut slot = vec![0u32; n];
        let mut free: Vec<u32> = Vec::new();
        let mut nslots = 0u32;
        let mut code = Vec::with_capacity(n);
        for (i, ins) in ssa.iter().enumerate() {
            let (_, k) = ins.args();
            let mut c = *ins;
            if k > 0 {
                c.a = slot[ins.a as usize];
            }
            if k > 1 {
                c.b = slot[ins.b as usize];
            }
            if k > 2 {
                c.c = slot[ins.c as usize];
            }
            // Free argument slots whose last use is here; the output may reuse one.
            let (args, k) = uses(i);
            for j in 0..k {
                let a = args[j] as usize;
                if last[a] == i && !args[..j].contains(&args[j]) {
                    free.push(slot[a]);
                }
            }
            let s = free.pop().unwrap_or_else(|| {
                nslots += 1;
                nslots - 1
            });
            slot[i] = s;
            c.out = s;
            code.push(c);
        }
        let apart = (0..n).filter_map(|i| Self::apart_of(&ssa, i).map(|(c, neg)| (i as u32, slot[c as usize], neg))).collect();
        let (mut fixed, mut pmask) = (vec![false; n], vec![0u64; n]);
        for (i, ins) in ssa.iter().enumerate() {
            let (args, k) = ins.args();
            fixed[i] = match ins.op {
                Op::Const(_) | Op::Var(Var::Param(_)) => true,
                Op::Un(_) | Op::Bin(_) | Op::Select | Op::Tag(_) => args[..k].iter().all(|&a| fixed[a as usize]),
                Op::Var(_) | Op::Grid(..) | Op::Read | Op::Lipschitz | Op::Frame | Op::Deriv(_) | Op::Real | Op::Within(_) => false,
            };
            pmask[i] = match ins.op {
                Op::Var(Var::Param(p)) => 1 << (p % 64),
                _ if fixed[i] => args[..k].iter().fold(0, |m, &a| m | pmask[a as usize]),
                _ => 0,
            };
        }
        // The params each fixed instruction reads: a param its own; else
        // its arguments', an argument's list where that holds them all.
        let (mut range, mut plist) = (vec![(0u32, 0u32); n], Vec::new());
        for (i, ins) in ssa.iter().enumerate().filter(|&(i, _)| fixed[i] && pmask[i] != 0) {
            let (args, k) = ins.args();
            let mut u: Vec<u32> = match ins.op {
                Op::Var(Var::Param(p)) => vec![p],
                _ => args[..k].iter().flat_map(|&a| plist[range[a as usize].0 as usize..(range[a as usize].0 + range[a as usize].1) as usize].to_vec()).collect(),
            };
            u.sort_unstable();
            u.dedup();
            range[i] = match args[..k].iter().map(|&a| range[a as usize]).find(|r| r.1 as usize == u.len()) {
                Some(r) => r,
                None => {
                    let start = plist.len() as u32;
                    plist.extend(&u);
                    (start, u.len() as u32)
                }
            };
        }
        let pread = (0..n).filter(|&i| range[i].1 > 0).map(|i| (i as u32, range[i].0, range[i].1)).collect();
        let mut slopes = vec![false; n];
        for i in (0..n).rev() {
            if let Op::Deriv(k) = ssa[i].op
                && k < 3
            {
                slopes[ssa[i].b as usize] = true;
            }
            if slopes[i] {
                let (args, k) = ssa[i].args();
                args[..k].iter().for_each(|&a| slopes[a as usize] = true);
            }
        }
        if !slopes.contains(&true) {
            slopes = Vec::new();
        }
        // Selects whose conditions read no x, y or z: those may switch at
        // a moment of the scope (eval::moment), where scope_run forks them.
        let mut space = vec![false; n];
        for (i, ins) in ssa.iter().enumerate() {
            let (args, k) = Self::reads(&ssa, i);
            space[i] = matches!(ins.op, Op::Var(Var::X | Var::Y | Var::Z)) || args[..k].iter().any(|&a| space[a as usize]);
        }
        let jumps = (0..n as u32)
            .filter(|&i| match ssa[i as usize].op {
                Op::Un(UnOp::Floor) | Op::Bin(BinOp::Atan2) => true,
                Op::Select => !space[ssa[i as usize].a as usize],
                _ => false,
            })
            .collect();
        Tape {
            ssa,
            code,
            nslots: nslots as usize,
            grids: Vec::new(),
            opaque: Vec::new(),
            fixed,
            pmask,
            pread,
            plist,
            slopes,
            forms: Vec::new(),
            jumps,
            joins: Default::default(),
            promises: Default::default(),
            axes: Default::default(),
            levels: Default::default(),
            apart,
        }
    }

    /// Of min or max instruction `i` whose operands are a and a + c or
    /// a - c (`apart_of`): the slot holding c, and whether the second less
    /// the first is -c.
    pub fn apart(&self, i: usize) -> Option<(u32, bool)> {
        if self.apart.is_empty() {
            return None;
        }
        self.apart.binary_search_by_key(&(i as u32), |x| x.0).ok().map(|k| (self.apart[k].1, self.apart[k].2))
    }

    /// Which of x, y and z the field's value reads, not through what only
    /// bounds it (a Lipschitz promise's constant and frame, a derivative's
    /// mark), as graph::Graph::value_depends_on: along an axis it reads
    /// not, the field is the same all along (a 2D profile's z).
    pub fn space_read(&self) -> [bool; 3] {
        let r = self.axes()[self.ssa.len() - 1];
        std::array::from_fn(|k| r & 1 << k != 0)
    }

    /// The space axes each value reads (as space_read), a bit each.
    fn axes(&self) -> Vec<u8> {
        let mut reads = vec![0u8; self.ssa.len()];
        for i in 0..self.ssa.len() {
            reads[i] = match self.ssa[i].op {
                Op::Var(v) => [Var::X, Var::Y, Var::Z].iter().position(|&w| w == v).map_or(0, |a| 1 << a),
                Op::Lipschitz | Op::Deriv(_) => reads[self.ssa[i].a as usize],
                _ => {
                    let (args, m) = Self::reads(&self.ssa, i);
                    args[..m].iter().fold(0, |r, &a| r | reads[a as usize])
                }
            };
        }
        reads
    }

    fn joins(&self) -> &Joins {
        self.joins.get_or_init(|| {
            let (ssa, n) = (&self.ssa, self.ssa.len());
            let mut last = vec![0u32; n];
            let mut count = vec![0u32; n];
            for i in 0..n {
                let (args, k) = Self::reads(ssa, i);
                for &a in &args[..k] {
                    last[a as usize] = i as u32;
                    count[a as usize] += 1;
                }
            }
            last[n - 1] = u32::MAX;
            // Backwards from the output: read once (so by its last reader),
            // by an operation monotone in that operand, itself monotone.
            let mut monotone = vec![false; n];
            monotone[n - 1] = true;
            for i in (0..n.saturating_sub(1)).rev() {
                let u = last[i] as usize;
                monotone[i] = count[i] == 1 && monotone[u] && ends_mapped(ssa, u, i as u32);
            }
            let (mut rejoin, mut dep) = (vec![0u32; n], vec![false; n]);
            for &j in &self.jumps {
                rejoin[j as usize] = join_of(ssa, &last, &monotone, &mut dep, j as usize) as u32;
            }
            let forked = forked_steps(&self.jumps, &rejoin, 0, n);
            Joins { rejoin, forked }
        })
    }

    /// Lipschitz instruction `i`'s operands as interval evaluation reads
    /// them again (Promise), found once.
    pub fn promise(&self, i: usize) -> &Promise {
        let all = self.promises.get_or_init(|| (0..self.ssa.len() as u32).filter(|&j| matches!(self.ssa[j as usize].op, Op::Lipschitz)).map(|j| (j, Default::default())).collect());
        let k = all.binary_search_by_key(&(i as u32), |p| p.0).expect("a promise");
        all[k].1.get_or_init(|| Promise::of(&self.ssa, i))
    }

    /// Where the bounds of a jump at instruction `j` (Tape::jumps;
    /// eval::forked) can be joined again: the first instruction from j on
    /// through which every later use of j's value passes, and from which the
    /// output is monotone (Joins). The output at the latest.
    pub fn rejoin(&self, j: usize) -> usize {
        self.joins().rejoin[j] as usize
    }

    /// The steps a run over the tape takes where every jump is straddled
    /// and forked (eval::forked: each window bounded twice, the forks
    /// inside it each time), saturating.
    pub fn forked_steps(&self) -> u64 {
        self.joins().forked
    }

    /// Bit k of each spatial axis k (x, y, z) instruction `i`'s value reads.
    pub fn axes_read(&self, i: usize) -> u8 {
        self.axes.get_or_init(|| {
            let mut m: Vec<u8> = Vec::with_capacity(self.ssa.len());
            for (j, ins) in self.ssa.iter().enumerate() {
                let own = match ins.op {
                    Op::Var(Var::X) => 1,
                    Op::Var(Var::Y) => 2,
                    Op::Var(Var::Z) => 4,
                    _ => 0,
                };
                let (args, k) = Self::reads(&self.ssa, j);
                m.push(args[..k].iter().fold(own, |a, &x| a | m[x as usize]));
            }
            m
        })[i]
    }

    /// Hash all of the tape into `fp`: every instruction of both views,
    /// the grids and forms it reads, and what each instruction is marked
    /// with, its tags and params by their numbers in the question (`Names`).
    /// Equal hashes mean the same tape, instruction for instruction, up to
    /// the model's numbering of names (and hash collisions), so the same
    /// answer to any query on it in the same context: for keys of cached
    /// answers (ink_kernel::cache).
    pub fn fingerprint(&self, fp: &mut Fingerprint, names: &Names) {
        // Every field, named: one added must be hashed too.
        // `jumps`, `joins`, `promises`, `axes`, `levels` and `apart` are derived from `ssa`; `pmask`, `pread` and
        // `plist` enter as what they decide (below).
        let Tape { ssa, code, nslots, grids, opaque, fixed, pmask: _, pread: _, plist: _, slopes, forms, jumps: _, joins: _, promises: _, axes: _, levels: _, apart: _ } = self;
        for list in [ssa, code] {
            fp.word(list.len() as u64);
            for &Ins { op, out, a, b, c } in list {
                match op {
                    Op::Const(v) => fp.word(0).f64(v),
                    Op::Var(v) => fp.word(1).word(match v {
                        Var::X => 0,
                        Var::Y => 1,
                        Var::Z => 2,
                        Var::T => 3,
                        Var::Param(k) => 4 + names.param(k) as u64,
                    }),
                    Op::Un(u) => fp.word(2).word(u as u64),
                    Op::Bin(b) => fp.word(3).word(b as u64),
                    Op::Select => fp.word(4),
                    Op::Tag(t) => fp.word(5).word(names.tag(t) as u64),
                    Op::Grid(g, m) => fp.word(6).word(g as u64).word(m as u64),
                    Op::Read => fp.word(7),
                    Op::Lipschitz => fp.word(8),
                    Op::Frame => fp.word(9),
                    Op::Deriv(k) => fp.word(10).word(k as u64),
                    Op::Real => fp.word(11),
                    Op::Within(id) => fp.word(12).word(id),
                };
                fp.word(out as u64).word(a as u64).word(b as u64).word(c as u64);
            }
        }
        fp.word(*nslots as u64).word(grids.len() as u64);
        grids.iter().for_each(|g| g.fingerprint(fp));
        fp.word(opaque.len() as u64);
        opaque.iter().for_each(|&t| _ = fp.word(names.tag(t) as u64));
        fp.word(fixed.len() as u64);
        fixed.iter().for_each(|f| _ = fp.word(*f as u64));
        // Which fixed instructions read a ranged param, exactly
        // (Tape::pread): those it takes in real numbers.
        for i in (0..fixed.len()).filter(|&i| fixed[i]) {
            fp.word(self.params_of(i).iter().any(|k| names.ranged.binary_search(k).is_ok()) as u64);
        }
        fp.word(slopes.len() as u64);
        slopes.iter().for_each(|s| _ = fp.word(*s as u64));
        fp.word(forms.len() as u64);
        forms.iter().for_each(|f| f.fingerprint(fp, names));
    }

    /// The params fixed instruction `i` reads, exactly, sorted (`pread`):
    /// none for any other.
    pub fn params_of(&self, i: usize) -> &[u32] {
        match self.pread.binary_search_by_key(&(i as u32), |e| e.0) {
            Ok(j) => {
                let (_, start, len) = self.pread[j];
                &self.plist[start as usize..(start + len) as usize]
            }
            Err(_) => &[],
        }
    }

    /// `pread` and `plist` of a tape whose instruction `n` is this one's
    /// `old` (`kept`, by n), each list copied once.
    fn reread(&self, kept: impl Iterator<Item = (u32, usize)>) -> (Vec<(u32, u32, u32)>, Vec<u32>) {
        let (mut pread, mut plist, mut at) = (Vec::new(), Vec::new(), std::collections::HashMap::new());
        for (n, old) in kept {
            if let Ok(j) = self.pread.binary_search_by_key(&(old as u32), |e| e.0) {
                let (_, start, len) = self.pread[j];
                let to = *at.entry(start).or_insert_with(|| {
                    plist.extend_from_slice(&self.plist[start as usize..(start + len) as usize]);
                    (plist.len() - len as usize) as u32
                });
                pread.push((n, to, len));
            }
        }
        (pread, plist)
    }

    /// The params the tape's value depends on (its forms' included), each
    /// once: those its instructions read as evaluation reads them
    /// (`pread`), which a simplified tape keeps where it folded a param
    /// away (a ranged one still takes the arithmetic it fed in real
    /// numbers).
    pub fn params_read(&self) -> Vec<u32> {
        let mut out = self.plist.clone();
        out.extend(self.forms.iter().flat_map(|f| f.params_read()));
        out.sort_unstable();
        out.dedup();
        out
    }

    pub fn len(&self) -> usize {
        self.ssa.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ssa.is_empty()
    }

    pub fn out_slot(&self) -> usize {
        self.code.last().unwrap().out as usize
    }

    /// The tape computing instruction `i` (an ssa index) alone: its cone of
    /// instructions, in time proportional to that cone.
    pub fn cone(&self, i: usize) -> Tape {
        self.cone_max(&[i])
    }

    /// The tape computing the max of instructions `roots` (ssa indices;
    /// one: that instruction alone, `cone`): their cones, then the max of
    /// their values, in time proportional to those cones.
    pub fn cone_max(&self, roots: &[usize]) -> Tape {
        self.cone_terms(&roots.iter().map(|&r| (r, true)).collect::<Vec<_>>())
    }

    /// The tape computing the max of terms `roots` (ssa indices, each
    /// negated where its flag is false: a term of `conjuncts`), as
    /// `cone_max`.
    pub fn cone_terms(&self, roots: &[(usize, bool)]) -> Tape {
        let mut seen: std::collections::HashSet<u32> = roots.iter().map(|&(r, _)| r as u32).collect();
        let mut stack: Vec<u32> = seen.iter().copied().collect();
        while let Some(j) = stack.pop() {
            let (args, k) = self.ssa[j as usize].args();
            for &a in &args[..k] {
                if seen.insert(a) {
                    stack.push(a);
                }
            }
        }
        let mut ids: Vec<u32> = seen.into_iter().collect();
        ids.sort_unstable();
        let remap: std::collections::HashMap<u32, u32> = ids.iter().enumerate().map(|(n, &j)| (j, n as u32)).collect();
        let mut ssa: Vec<Ins> = ids
            .iter()
            .enumerate()
            .map(|(n, &j)| {
                let mut ins = self.ssa[j as usize];
                let (_, k) = ins.args();
                if k > 0 {
                    ins.a = remap[&ins.a];
                }
                if k > 1 {
                    ins.b = remap[&ins.b];
                }
                if k > 2 {
                    ins.c = remap[&ins.c];
                }
                ins.out = n as u32;
                ins
            })
            .collect();
        // Each root negated where asked, and after the first joined by a
        // max, the last the output (one root as it is: its own, ids' last,
        // as it reads every other).
        let term = |ssa: &mut Vec<Ins>, (r, up): (usize, bool)| {
            let r = remap[&(r as u32)];
            if up {
                return r;
            }
            let n = ssa.len() as u32;
            ssa.push(Ins { op: Op::Un(UnOp::Neg), out: n, a: r, b: 0, c: 0 });
            n
        };
        let mut out = term(&mut ssa, roots[0]);
        for &r in &roots[1..] {
            let b = term(&mut ssa, r);
            let n = ssa.len() as u32;
            ssa.push(Ins { op: Op::Bin(BinOp::Max), out: n, a: out, b, c: 0 });
            out = n;
        }
        // As this tape has them (a simplified tape's are not recomputable);
        // a joining max or negation is not fixed (bounded as any value is:
        // soundly).
        let joined = ssa.len() - ids.len();
        let mut t = Tape::from_ssa(ssa);
        t.fixed = ids.iter().map(|&j| self.fixed[j as usize]).chain(std::iter::repeat_n(false, joined)).collect();
        t.pmask = ids.iter().map(|&j| self.pmask[j as usize]).chain(std::iter::repeat_n(0, joined)).collect();
        (t.pread, t.plist) = self.reread(ids.iter().enumerate().map(|(n, &j)| (n as u32, j as usize)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        t.forms = self.forms.clone();
        t
    }

    /// The tape with each instruction `i` (an ssa index) of `values` the
    /// constant given (a floor's value over a region it is constant in).
    pub fn with_consts(&self, values: &[(usize, f64)]) -> Tape {
        if values.is_empty() {
            return self.clone();
        }
        let mut ssa = self.ssa.clone();
        for &(i, v) in values {
            ssa[i] = Ins { op: Op::Const(v), out: i as u32, a: 0, b: 0, c: 0 };
        }
        let mut t = Tape::from_ssa(ssa);
        (t.fixed, t.pmask) = (self.fixed.clone(), self.pmask.clone());
        for &(i, _) in values {
            (t.fixed[i], t.pmask[i]) = (true, 0);
        }
        (t.pread, t.plist) = self.reread((0..self.ssa.len()).filter(|i| !values.iter().any(|v| v.0 == *i)).map(|i| (i as u32, i)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        t.forms = self.forms.clone();
        t
    }

    /// Whether select `i` is continuous where it switches, by its form, in
    /// the region `choices` were found for: its condition is x - y (or
    /// y - x), x and y continuous there (`continuous`: else the condition
    /// may jump past 0, never meeting x = y), and its branches one
    /// expression but for x and y swapped (g(x) and g(y): select(x - y, x,
    /// y) is a min), so where the condition reaches 0, x = y, they are one
    /// value.
    pub fn continuous_select(&self, i: usize, choices: &[u8]) -> bool {
        let ssa = &self.ssa;
        let strip = |mut j: usize| {
            while let Op::Tag(_) = ssa[j].op {
                j = ssa[j].a as usize;
            }
            j
        };
        let c = strip(ssa[i].a as usize);
        let (x, y) = match ssa[c].op {
            Op::Bin(BinOp::Sub) => (strip(ssa[c].a as usize), strip(ssa[c].b as usize)),
            _ => return false,
        };
        if !self.continuous(x, choices) || !self.continuous(y, choices) {
            return false;
        }
        // p and q one expression but for x and y swapped (memoised: the
        // tape is a DAG).
        fn alike(ssa: &[Ins], strip: &dyn Fn(usize) -> usize, x: usize, y: usize, p: usize, q: usize, seen: &mut std::collections::HashMap<(usize, usize), bool>) -> bool {
            let (p, q) = (strip(p), strip(q));
            if p == q || (p, q) == (x, y) || (p, q) == (y, x) {
                return true;
            }
            if let Some(&r) = seen.get(&(p, q)) {
                return r;
            }
            seen.insert((p, q), false);
            let (a, b) = (&ssa[p], &ssa[q]);
            let same_op = match (a.op, b.op) {
                (Op::Const(u), Op::Const(v)) => u.to_bits() == v.to_bits(),
                // One number by identity alone, never by its bounds.
                (Op::Within(u), Op::Within(v)) => u == v,
                (u, v) => std::mem::discriminant(&u) == std::mem::discriminant(&v) && format!("{u:?}") == format!("{v:?}"),
            };
            let r = same_op && {
                let ((pa, k), (qa, _)) = (a.args(), b.args());
                let mut go = |i: usize, j: usize| alike(ssa, strip, x, y, pa[i] as usize, qa[j] as usize, seen);
                // A sum, product, min or max in either order.
                let commutes = matches!(a.op, Op::Bin(BinOp::Add | BinOp::Mul | BinOp::Min | BinOp::Max | BinOp::Emax));
                (0..k).all(|j| go(j, j)) || commutes && go(0, 1) && go(1, 0)
            };
            seen.insert((p, q), r);
            r
        }
        alike(ssa, &strip, x, y, ssa[i].b as usize, ssa[i].c as usize, &mut std::collections::HashMap::new())
    }

    /// Whether value `j` is continuous in the region `choices` were found
    /// for, by what it is built of: variables, constants, sums, products,
    /// negations, magnitudes, squares, square roots, sines, cosines,
    /// exponentials, degrees and radians, mins and maxes, and selects (and
    /// mins and maxes) decided there, through the branch they take. A
    /// floor, quotient, select that may switch, or anything else is not
    /// taken as continuous.
    pub fn continuous(&self, j: usize, choices: &[u8]) -> bool {
        let n = self.ssa.len();
        let mut seen = vec![false; n];
        let mut stack = vec![j];
        while let Some(i) = stack.pop() {
            if std::mem::replace(&mut seen[i], true) {
                continue;
            }
            let ins = &self.ssa[i];
            let next: &[u32] = match (ins.op, choices[i]) {
                // A number known within bounds is one constant (its
                // identity's), whatever its bounds' width.
                (Op::Var(_) | Op::Const(_) | Op::Frame | Op::Within(_), _) => &[],
                (Op::Tag(_) | Op::Lipschitz | Op::Real, _) => &[ins.a],
                (Op::Un(UnOp::Neg | UnOp::Abs | UnOp::Sqrt | UnOp::Square | UnOp::Exp | UnOp::Sin | UnOp::Cos | UnOp::Rad | UnOp::Deg), _) => &[ins.a],
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => &[ins.a],
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => &[ins.b],
                (Op::Select, FIRST) => &[ins.b],
                (Op::Select, SECOND) => &[ins.c],
                (Op::Bin(BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Min | BinOp::Max | BinOp::Emax), _) => &[ins.a, ins.b],
                _ => return false,
            };
            stack.extend(next.iter().map(|&a| a as usize));
        }
        true
    }

    /// Whether a select the output reads (through the branches `choices`
    /// pick) may switch in the region they were found for, between branches
    /// that are not one field: a seam, which bounds cannot show continuous
    /// (a sweep's mitre between two segments' fields).
    pub fn switches(&self, choices: &[u8]) -> bool {
        let n = self.ssa.len();
        let strip = |mut j: usize| {
            while let Op::Tag(_) = self.ssa[j].op {
                j = self.ssa[j].a as usize;
            }
            j
        };
        let mut live = vec![false; n];
        live[n - 1] = true;
        for i in (0..n).rev() {
            if !live[i] {
                continue;
            }
            let ins = &self.ssa[i];
            let picked = match (ins.op, choices[i]) {
                (Op::Select, BOTH) if strip(ins.b as usize) != strip(ins.c as usize) && !self.continuous_select(i, choices) => return true,
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => Some(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => Some(ins.b),
                (Op::Select, FIRST) => Some(ins.b),
                (Op::Select, SECOND) => Some(ins.c),
                _ => None,
            };
            match picked {
                Some(a) => live[a as usize] = true,
                None => {
                    let (args, k) = ins.args();
                    for &a in &args[..k] {
                        live[a as usize] = true;
                    }
                }
            }
        }
        false
    }

    /// Whether the output's 0 over a region has no volume (holds no ball),
    /// from `grads` found over it (eval_interval_grad_every's). A value's 0
    /// has none where its bounds exclude 0, or its rate along an axis keeps
    /// one sign throughout (a field 0 on a ball has rate 0 there); a min's,
    /// max's or select's 0 is among its branches' that may be picked (where two
    /// may tie, both), a
    /// product's its factors', so it has none where theirs have none (each
    /// such field is piecewise continuous: a 0 with no volume is one of
    /// measure 0, and a union of those holds no ball); a negation's, an
    /// absolute value's or a square's is its argument's.
    pub fn zeros_thin(&self, grads: &[crate::interval::G]) -> bool {
        self.thin_zeros(grads)[self.ssa.len() - 1]
    }

    /// `zeros_thin` of every value.
    fn thin_zeros(&self, grads: &[crate::interval::G]) -> Vec<bool> {
        let n = self.ssa.len();
        let mut thin: Vec<bool> = Vec::with_capacity(n);
        for (i, ins) in self.ssa.iter().enumerate() {
            let g = &grads[i];
            let t = |j: u32| thin[j as usize];
            let v = |j: u32| grads[j as usize].v;
            let by_bounds = g.v.lo > 0.0 || g.v.hi < 0.0 || g.d.iter().any(|d| d.lo > 0.0 || d.hi < 0.0);
            let by_args = match ins.op {
                Op::Tag(_) | Op::Lipschitz | Op::Real | Op::Un(UnOp::Neg | UnOp::Abs | UnOp::Square) => t(ins.a),
                Op::Bin(BinOp::Mul) => t(ins.a) && t(ins.b),
                // A branch may be picked unless the other is below it (above,
                // for a max) throughout: values that may tie both may be.
                Op::Bin(BinOp::Min) => (v(ins.a).lo > v(ins.b).hi || t(ins.a)) && (v(ins.b).lo > v(ins.a).hi || t(ins.b)),
                Op::Bin(BinOp::Max | BinOp::Emax) => (v(ins.a).hi < v(ins.b).lo || t(ins.a)) && (v(ins.b).hi < v(ins.a).lo || t(ins.b)),
                Op::Select => (v(ins.a).none_inside() || t(ins.b)) && (v(ins.a).all_inside() || t(ins.c)),
                _ => false,
            };
            thin.push(by_bounds || by_args);
        }
        thin
    }

    /// Whether every 0 of the output over a region has points above 0
    /// beside it (it is no highest point of the field there), from `grads`
    /// found over it (eval_interval_grad_every's). Each value's 0s have
    /// points below 0 beside them, and points above, where its bounds
    /// exclude 0 or its rate along an axis keeps one sign throughout (at a
    /// point inside the region); a negation swaps the two, an absolute
    /// value's or a square's 0s have points above where its argument's
    /// have either, a product with a value of one sign throughout keeps
    /// (or swaps) its factor's; a min's 0 is a 0 of a branch that may be
    /// picked, with points below beside it where that branch's have, and
    /// points above where that branch's have and the other is above 0
    /// throughout; a max the reverse. Each such field is piecewise
    /// continuous.
    pub fn zeros_rise(&self, grads: &[crate::interval::G]) -> bool {
        let mut beside: Vec<(bool, bool)> = Vec::with_capacity(self.ssa.len());
        for (i, ins) in self.ssa.iter().enumerate() {
            let g = &grads[i];
            let s = |j: u32| beside[j as usize];
            let v = |j: u32| grads[j as usize].v;
            let swap = |(lo, hi): (bool, bool)| (hi, lo);
            let by_bounds = g.v.lo > 0.0 || g.v.hi < 0.0 || g.d.iter().any(|d| d.lo > 0.0 || d.hi < 0.0);
            let by_args = match ins.op {
                Op::Tag(_) | Op::Lipschitz | Op::Real => s(ins.a),
                Op::Un(UnOp::Neg) => swap(s(ins.a)),
                Op::Un(UnOp::Abs | UnOp::Square) => (false, s(ins.a).0 || s(ins.a).1),
                Op::Bin(BinOp::Mul) => match (v(ins.a), v(ins.b)) {
                    (_, w) if w.lo > 0.0 => s(ins.a),
                    (_, w) if w.hi < 0.0 => swap(s(ins.a)),
                    (w, _) if w.lo > 0.0 => s(ins.b),
                    (w, _) if w.hi < 0.0 => swap(s(ins.b)),
                    _ => (false, false),
                },
                // A branch is never picked where the other is below it
                // (above, for a max) throughout.
                Op::Bin(BinOp::Min) => {
                    let (na, nb) = (v(ins.a).lo > v(ins.b).hi, v(ins.b).lo > v(ins.a).hi);
                    ((na || s(ins.a).0) && (nb || s(ins.b).0), (na || s(ins.a).1 && v(ins.b).lo > 0.0) && (nb || s(ins.b).1 && v(ins.a).lo > 0.0))
                }
                Op::Bin(BinOp::Max | BinOp::Emax) => {
                    let (na, nb) = (v(ins.a).hi < v(ins.b).lo, v(ins.b).hi < v(ins.a).lo);
                    ((na || s(ins.a).0 && v(ins.b).hi < 0.0) && (nb || s(ins.b).0 && v(ins.a).hi < 0.0), (na || s(ins.a).1) && (nb || s(ins.b).1))
                }
                _ => (false, false),
            };
            beside.push((by_bounds || by_args.0, by_bounds || by_args.1));
        }
        beside.last().is_some_and(|b| b.1)
    }

    /// Whether the output changes sign at most once along unit direction `d`
    /// on every line through a region (from below 0 to not, or the
    /// reverse, the same way on every line), from `choices` and `grads`
    /// found over it (eval_interval's and eval_interval_grad_every's). Then
    /// its zero set there, but for 0s of no volume inside the material
    /// (drawn as material, crate::sign::drawn), is one sheet: no piece of it
    /// lies wholly inside.
    /// Each value is so where its rate along d, carried by the chain rule (a
    /// min's, max's or select's undecided the hull of its branches'),
    /// excludes 0; where it keeps one sign throughout (at most 0 with a 0
    /// of no volume counts as below it: pieces meeting flush, drawn as
    /// material); and where its sign is that of a value that is so (or its
    /// negation's, the other way): a min or max (or emax, whose sign is
    /// its max's) of two that change sign the same way (where one is not
    /// below 0 is a union of such, and an intersection), or of one and one
    /// that keeps the sign that leaves the other's (a max with a value
    /// below 0 throughout, a min with one not below 0), a product with a
    /// value above 0 throughout, and |x| - y where x + y is above 0
    /// throughout (its sign x - y's), or x - y below 0 (its sign -x - y's),
    /// either's rate excluding 0.
    pub fn monotone(&self, choices: &[u8], grads: &[crate::interval::G], d: crate::V3) -> bool {
        matches!(self.sheet(choices, grads, None, d), Sheet::One)
    }

    /// `monotone`, or else, where the output's sign set is that of a select
    /// switching between two branches each changing sign once the same way
    /// (a seam that bounds cannot show continuous, `switches`), how thick
    /// what it may hide is: the field lies between its branches' min and
    /// max, so along each line along d every sign change is within
    /// sup |b - c| / (their least rate) of the others, sup |b - c| bounded
    /// by its value at the region's centre (`centre`: the values there,
    /// eval_interval_grad_every over that point) and its gradient's bounds
    /// over the region's half extents (`centre`'s second). Carried through
    /// what keeps the sign set (as in `monotone`), two such summed.
    pub fn sheet(&self, choices: &[u8], grads: &[crate::interval::G], centre: Option<(&[crate::interval::G], crate::V3)>, d: crate::V3) -> Sheet {
        self.sheet_way(choices, grads, centre, d).0
    }

    /// `sheet`, and the way the output changes sign along d where it does
    /// (1: below 0, then not; -1: the reverse; 2: it keeps one sign).
    pub fn sheet_way(&self, choices: &[u8], grads: &[crate::interval::G], centre: Option<(&[crate::interval::G], crate::V3)>, d: crate::V3) -> (Sheet, i8) {
        use crate::interval::{I, add_o, div_o, hull, mul_o, sub_o};
        let n = self.ssa.len();
        let mut r: Vec<I> = Vec::with_capacity(n);
        // Each value's way of changing sign along d: 1 (below 0, then
        // not), -1 (not, then below), 0 (unknown); and how thick what it
        // may hide is (0: nothing).
        let mut way: Vec<i8> = Vec::with_capacity(n);
        let mut thick: Vec<f64> = Vec::with_capacity(n);
        let v = |j: u32| grads[j as usize].v;
        // Inside throughout (or at most 0, its 0 of no volume: pieces
        // meeting flush, drawn as material, crate::sign::solid), or
        // nowhere (crate::sign).
        let thin = self.thin_zeros(grads);
        let strip = |mut j: usize| {
            while let Op::Tag(_) = self.ssa[j].op {
                j = self.ssa[j].a as usize;
            }
            j
        };
        let below = |j: u32| v(j).all_inside() || v(j).hi <= 0.0 && thin[j as usize];
        let not_below = |j: u32| v(j).none_inside();
        for i in 0..n {
            let ins = &self.ssa[i];
            let at = |j: u32| r[j as usize];
            let rate = match (ins.op, choices[i]) {
                (Op::Var(Var::X), _) => I::point(d[0]),
                (Op::Var(Var::Y), _) => I::point(d[1]),
                (Op::Var(Var::Z), _) => I::point(d[2]),
                (Op::Var(_) | Op::Const(_) | Op::Within(_), _) => I::point(0.0),
                (Op::Tag(_) | Op::Lipschitz | Op::Real, _) => at(ins.a),
                (Op::Bin(BinOp::Add), _) => add_o(at(ins.a), at(ins.b)),
                (Op::Bin(BinOp::Sub), _) => sub_o(at(ins.a), at(ins.b)),
                (Op::Un(UnOp::Neg), _) => I::new(-at(ins.a).hi, -at(ins.a).lo),
                (Op::Bin(BinOp::Mul), _) => add_o(mul_o(v(ins.a), at(ins.b)), mul_o(v(ins.b), at(ins.a))),
                (Op::Un(UnOp::Square), _) => mul_o(mul_o(I::point(2.0), v(ins.a)), at(ins.a)),
                (Op::Un(UnOp::Sqrt), _) if v(ins.a).lo > 0.0 => {
                    let (lo, hi) = (v(ins.a).lo, v(ins.a).hi);
                    let k = I::new(0.5 / hi.sqrt() * (1.0 - 1e-15), 0.5 / lo.sqrt() * (1.0 + 1e-15));
                    mul_o(k, at(ins.a))
                }
                (Op::Un(UnOp::Abs), _) if v(ins.a).lo > 0.0 => at(ins.a),
                (Op::Un(UnOp::Abs), _) if v(ins.a).hi < 0.0 => I::new(-at(ins.a).hi, -at(ins.a).lo),
                // Its argument's rate, either sign.
                (Op::Un(UnOp::Abs), _) => hull(at(ins.a), I::new(-at(ins.a).hi, -at(ins.a).lo)),
                // (a' b - a b') / b^2, where b is not 0.
                (Op::Bin(BinOp::Div), _) if v(ins.b).lo > 0.0 || v(ins.b).hi < 0.0 => {
                    div_o(sub_o(mul_o(at(ins.a), v(ins.b)), mul_o(v(ins.a), at(ins.b))), mul_o(v(ins.b), v(ins.b)))
                }
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => at(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => at(ins.b),
                (Op::Bin(BinOp::Min | BinOp::Max), BOTH) => hull(at(ins.a), at(ins.b)),
                // max(a, b) where either is at most 0, else hypot(a, b),
                // whose rate is its sides' weighted by a / hypot and
                // b / hypot (squares summing to 1): their hull, up to
                // sqrt 2 times.
                (Op::Bin(BinOp::Emax), BOTH) => {
                    let h = hull(at(ins.a), at(ins.b));
                    let s = std::f64::consts::SQRT_2 * (1.0 + 1e-15);
                    I::new(h.lo.min(s * h.lo), h.hi.max(s * h.hi))
                }
                (Op::Select, FIRST) => at(ins.b),
                (Op::Select, SECOND) => at(ins.c),
                // A select that may switch along d jumps where it does,
                // unless its condition is constant along d (it does not
                // switch on a line along d) or its branches are one field.
                (Op::Select, BOTH) => {
                    let strip = |mut j: usize| {
                        while let Op::Tag(_) = self.ssa[j].op {
                            j = self.ssa[j].a as usize;
                        }
                        j
                    };
                    let still = at(ins.a) == I::point(0.0);
                    match still || strip(ins.b as usize) == strip(ins.c as usize) || self.continuous_select(i, choices) {
                        true => hull(at(ins.b), at(ins.c)),
                        false => I::new(f64::NEG_INFINITY, f64::INFINITY),
                    }
                }
                // Anything else: its gradient's bounds, along d.
                _ => {
                    let g = grads[i].d;
                    add_o(add_o(mul_o(g[0], I::point(d[0])), mul_o(g[1], I::point(d[1]))), mul_o(g[2], I::point(d[2])))
                }
            };
            // Both its chain rule's bounds and its gradient's, along d, hold:
            // where they meet.
            let g = grads[i].d;
            let by_grad = add_o(add_o(mul_o(g[0], I::point(d[0])), mul_o(g[1], I::point(d[1]))), mul_o(g[2], I::point(d[2])));
            let rate = match (rate.lo.max(by_grad.lo), rate.hi.min(by_grad.hi)) {
                (lo, hi) if lo <= hi => I::new(lo, hi),
                _ => rate,
            };
            r.push(rate);
            let w = |j: u32| way[j as usize];
            let by_rate = if rate.lo > 0.0 {
                1
            } else if rate.hi < 0.0 {
                -1
            } else {
                0
            };
            // 2: keeps one sign throughout (changes as either way does).
            let alike = |a: i8, b: i8| match (a, b) {
                (2, x) | (x, 2) => x,
                (a, b) if a == b => a,
                _ => 0,
            };
            let t = |j: u32| thick[j as usize];
            let both = |a: u32, b: u32| (alike(w(a), w(b)), t(a) + t(b));
            let by_sign: (i8, f64) = match (ins.op, choices[i]) {
                (Op::Tag(_) | Op::Lipschitz | Op::Real, _) => (w(ins.a), t(ins.a)),
                // A min, max or positive product maps which points are
                // inside exactly; a negation swaps inside and not only away
                // from 0 (a 0 is not inside either way), so it changes sign
                // the other way only where its argument's 0 has no volume,
                // and keeps one sign only where its own bounds show it.
                (Op::Un(UnOp::Neg), _) => match w(ins.a) {
                    x @ (1 | -1) if thin[ins.a as usize] => (-x, t(ins.a)),
                    _ => (0, 0.0),
                },
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => (w(ins.a), t(ins.a)),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => (w(ins.b), t(ins.b)),
                (Op::Select, FIRST) => (w(ins.b), t(ins.b)),
                (Op::Select, SECOND) => (w(ins.c), t(ins.c)),
                (Op::Bin(BinOp::Max | BinOp::Emax), BOTH) if below(ins.a) => (w(ins.b), t(ins.b)),
                (Op::Bin(BinOp::Max | BinOp::Emax), BOTH) if below(ins.b) => (w(ins.a), t(ins.a)),
                (Op::Bin(BinOp::Min), BOTH) if not_below(ins.a) => (w(ins.b), t(ins.b)),
                (Op::Bin(BinOp::Min), BOTH) if not_below(ins.b) => (w(ins.a), t(ins.a)),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), BOTH) => both(ins.a, ins.b),
                (Op::Bin(BinOp::Mul), _) if v(ins.b).lo > 0.0 => (w(ins.a), t(ins.a)),
                (Op::Bin(BinOp::Mul), _) if v(ins.a).lo > 0.0 => (w(ins.b), t(ins.b)),
                // |x| - y is below 0 where both x - y and -x - y are (a max
                // of the two): where one of them is below 0 throughout, as
                // the other is, so as its rate shows (a wall's face, in a
                // box reaching its middle but not its other face).
                (Op::Bin(BinOp::Sub), _) if matches!(self.ssa[strip(ins.a as usize)].op, Op::Un(UnOp::Abs)) => {
                    let (x, y) = (self.ssa[strip(ins.a as usize)].a, ins.b);
                    let rate = if add_o(v(x), v(y)).lo > 0.0 {
                        sub_o(r[x as usize], r[y as usize])
                    } else if sub_o(v(x), v(y)).hi < 0.0 {
                        let s = add_o(r[x as usize], r[y as usize]);
                        I::new(-s.hi, -s.lo)
                    } else {
                        I::point(0.0)
                    };
                    (i8::from(rate.lo > 0.0) - i8::from(rate.hi < 0.0), 0.0)
                }
                // A seam between branches that each change sign once, the
                // same way: what it may hide is within their sheets' reach.
                (Op::Select, BOTH) => match centre {
                    Some((at, half)) => {
                        let (rb, rc) = (r[ins.b as usize], r[ins.c as usize]);
                        let least = if rb.lo > 0.0 && rc.lo > 0.0 {
                            rb.lo.min(rc.lo)
                        } else if rb.hi < 0.0 && rc.hi < 0.0 {
                            (-rb.hi).min(-rc.hi)
                        } else {
                            0.0
                        };
                        let (gb, gc) = (grads[ins.b as usize].d, grads[ins.c as usize].d);
                        let gap = (at[ins.b as usize].v.hi - at[ins.c as usize].v.lo).abs().max((at[ins.c as usize].v.hi - at[ins.b as usize].v.lo).abs())
                            + (0..3).map(|k| crate::interval::mag(sub_o(gb[k], gc[k])) * half[k]).sum::<f64>();
                        match least > 0.0 && gap.is_finite() {
                            true => (if rb.lo > 0.0 { 1 } else { -1 }, gap / least * (1.0 + 1e-12)),
                            false => (0, 0.0),
                        }
                    }
                    None => (0, 0.0),
                },
                _ => (0, 0.0),
            };
            let kept = below(i as u32) || not_below(i as u32);
            let (wi, ti) = if kept {
                (2, 0.0)
            } else if by_rate != 0 {
                (by_rate, 0.0)
            } else {
                by_sign
            };
            way.push(wi);
            thick.push(ti);
        }
        let held = match (way[n - 1], thick[n - 1]) {
            (0, _) => Sheet::Unknown,
            (_, 0.0) => Sheet::One,
            (_, t) => Sheet::Thin(t),
        };
        (held, way[n - 1])
    }

    /// Whether the output's gradient may turn abruptly in a region, from
    /// `choices` and `grads` found over it (eval_interval's and
    /// eval_interval_grad_every's): through the branches they pick, a min,
    /// max, emax or select it reads is undecided between branches whose
    /// gradients there point apart throughout (interval::apart; two that
    /// may point alike, f and 2 f, meet in no crease), or an abs of a value
    /// that may be 0 there (its sides' gradients opposite), or a floor that
    /// may step there.
    pub fn may_crease(&self, choices: &[u8], grads: &[crate::interval::G]) -> bool {
        use crate::interval::apart;
        let neg = |d: [crate::interval::I; 3]| d.map(|i| crate::interval::I::new(-i.hi, -i.lo));
        let n = self.ssa.len();
        let mut live = vec![false; n];
        live[n - 1] = true;
        for i in (0..n).rev() {
            if !live[i] {
                continue;
            }
            let ins = &self.ssa[i];
            let g = |j: u32| grads[j as usize].d;
            let creased = match (ins.op, choices[i]) {
                (Op::Bin(BinOp::Min | BinOp::Max), BOTH) => apart(&g(ins.a), &g(ins.b)),
                // An emax is hypot (smooth) where both sides are above 0,
                // and turns into its max smoothly: it creases only where
                // both may be at most 0, as a max.
                (Op::Bin(BinOp::Emax), BOTH) => grads[ins.a as usize].v.lo <= 0.0 && grads[ins.b as usize].v.lo <= 0.0 && apart(&g(ins.a), &g(ins.b)),
                (Op::Select, BOTH) => apart(&g(ins.b), &g(ins.c)),
                (Op::Un(UnOp::Abs), _) => {
                    let v = grads[ins.a as usize].v;
                    v.lo <= 0.0 && v.hi >= 0.0 && apart(&g(ins.a), &neg(g(ins.a)))
                }
                (Op::Un(UnOp::Floor), _) => {
                    let v = grads[ins.a as usize].v;
                    v.lo.floor() != v.hi.floor()
                }
                _ => false,
            };
            if creased {
                return true;
            }
            let picked = match (ins.op, choices[i]) {
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => Some(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => Some(ins.b),
                (Op::Select, FIRST) => Some(ins.b),
                (Op::Select, SECOND) => Some(ins.c),
                _ => None,
            };
            match picked {
                Some(a) => live[a as usize] = true,
                None => {
                    let (args, k) = ins.args();
                    for &a in &args[..k] {
                        live[a as usize] = true;
                    }
                }
            }
        }
        false
    }

    /// Choices that keep the output's sign, not its value, over a region
    /// where instruction i lies in `vals[i]` (sound bounds): at a min, max
    /// or emax whose sign reaches the output's (through mins, maxes, emaxes,
    /// whose sign is their max's, negations, tags and products with a
    /// factor above 0 throughout: each maps where a value is below, at or
    /// above 0 to where the next is), the branch above 0 throughout dropped
    /// from a min and the one below 0 throughout from a max, where the other
    /// is not so too: the other branch's sign is then the min's or max's (a
    /// slab's faces past a profile's walls drop out of an extrusion inside
    /// its span; inside both, neither is, so how deep a point lies is still
    /// read along both). The tape simplified by them has the output's sign
    /// at every point of the region (below, at or above 0), not its value.
    /// `choices` starts with none made; true where it makes one.
    pub fn sign_choices(&self, vals: &[I], choices: &mut [u8]) -> bool {
        let n = self.ssa.len();
        // Each value: 0 not reached yet (dead, if never), 1 its sign
        // reaching the output's, 2 not.
        let mut way = vec![0u8; n];
        way[n - 1] = 1;
        let mut any = false;
        for i in (0..n).rev() {
            let (ins, w) = (&self.ssa[i], way[i]);
            if w == 0 {
                continue;
            }
            let (a, b) = (vals[ins.a as usize], vals[ins.b as usize]);
            let pick = match ins.op {
                _ if w != 1 => None,
                Op::Bin(BinOp::Min) if b.lo > 0.0 && !(a.lo > 0.0) => Some(ins.a),
                Op::Bin(BinOp::Min) if a.lo > 0.0 && !(b.lo > 0.0) => Some(ins.b),
                Op::Bin(BinOp::Max | BinOp::Emax) if b.hi < 0.0 && !(a.hi < 0.0) => Some(ins.a),
                Op::Bin(BinOp::Max | BinOp::Emax) if a.hi < 0.0 && !(b.hi < 0.0) => Some(ins.b),
                _ => None,
            };
            if let Some(k) = pick {
                choices[i] = if k == ins.a { FIRST } else { SECOND };
                any = true;
                way[k as usize] = way[k as usize].max(1);
                continue;
            }
            let (args, k) = Self::reads(&self.ssa, i);
            for &x in &args[..k] {
                let keeps = w == 1
                    && match ins.op {
                        Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) | Op::Un(UnOp::Neg) | Op::Tag(_) => true,
                        Op::Bin(BinOp::Mul) => ins.a != ins.b && vals[if x == ins.a { ins.b } else { ins.a } as usize].lo > 0.0,
                        _ => false,
                    };
                way[x as usize] = way[x as usize].max(if keeps { 1 } else { 2 });
            }
        }
        any
    }

    /// The values each point where the output is at most 0 has at most 0
    /// (`true`) or at least 0 (`false`), and the space axes each reads (a
    /// bit each): the output's terms joined by maxes (an emax is at most 0
    /// where its max is) and, negated, by mins. A term reading one axis
    /// alone bounds where along it the output can be at most 0 (a slab's
    /// span). Each value is visited once each way, so a term shared by
    /// many is listed once, and the walk is as long as the tape.
    pub fn conjuncts(&self) -> Vec<(usize, bool, u8)> {
        let n = self.ssa.len();
        let reads = self.axes();
        let mut seen = vec![[false; 2]; n];
        let (mut out, mut stack) = (Vec::new(), vec![(n - 1, true)]);
        while let Some((i, below)) = stack.pop() {
            if std::mem::replace(&mut seen[i][below as usize], true) {
                continue;
            }
            let ins = &self.ssa[i];
            match (ins.op, below) {
                (Op::Bin(BinOp::Max | BinOp::Emax), true) | (Op::Bin(BinOp::Min), false) => stack.extend([(ins.a as usize, below), (ins.b as usize, below)]),
                (Op::Un(UnOp::Neg), _) => stack.push((ins.a as usize, !below)),
                (Op::Tag(_), _) => stack.push((ins.a as usize, below)),
                _ => out.push((i, below, reads[i])),
            }
        }
        out
    }

    /// The operand min/max/select instruction `i` reduces to by `choices`
    /// (interval evaluation proved the other cannot win), if one.
    fn pick(&self, choices: &[u8], i: usize) -> Option<u32> {
        let ins = &self.ssa[i];
        match (ins.op, choices[i]) {
            (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => Some(ins.a),
            (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => Some(ins.b),
            (Op::Select, FIRST) => Some(ins.b),
            (Op::Select, SECOND) => Some(ins.c),
            _ => None,
        }
    }

    /// Per instruction, whether the output reads its value, through the
    /// branches `choices` leaves (the rest, wherever `choices` were found,
    /// change nothing the output is).
    pub fn live(&self, choices: &[u8]) -> Vec<bool> {
        let n = self.ssa.len();
        let mut live = vec![false; n];
        live[n - 1] = true;
        for i in (0..n).rev() {
            if !live[i] {
                continue;
            }
            if let Some(a) = self.pick(choices, i) {
                live[a as usize] = true;
            } else {
                let (args, k) = self.ssa[i].args();
                for &a in &args[..k] {
                    live[a as usize] = true;
                }
            }
        }
        live
    }

    /// Drop branches of min/max/select that interval evaluation proved can't win.
    pub fn simplify(&self, choices: &[u8]) -> Tape {
        let n = self.ssa.len();
        let pick = |i: usize| self.pick(choices, i);
        let live = self.live(choices);

        let mut remap = vec![u32::MAX; n];
        let mut out: Vec<Ins> = Vec::with_capacity(n);
        for i in 0..n {
            if !live[i] {
                continue;
            }
            if let Some(a) = pick(i) {
                remap[i] = remap[a as usize];
                continue;
            }
            let mut ins = self.ssa[i];
            let (_, k) = ins.args();
            if k > 0 {
                ins.a = remap[ins.a as usize];
            }
            if k > 1 {
                ins.b = remap[ins.b as usize];
            }
            if k > 2 {
                ins.c = remap[ins.c as usize];
            }
            ins.out = out.len() as u32;
            remap[i] = ins.out;
            out.push(ins);
        }
        // If the root was aliased, everything after its target is dead.
        // Held in its own length, not its parent's: the many simplified
        // tapes of an octree::Local would hold their parents' many times.
        out.truncate(remap[n - 1] as usize + 1);
        out.shrink_to_fit();
        let mut t = Tape::from_ssa(out);
        // What is fixed stays as it was: arithmetic on a branch picked for
        // the region still means real arithmetic, not a constants-and-params
        // value in floats (select(x, c, 2)^2 for x < 0 is c^2 in reals); and
        // so do the params it reads (a branch picked over a param's span
        // keeps reading the param).
        let kept = || (0..n).filter(|&i| live[i] && pick(i).is_none() && (remap[i] as usize) < t.fixed.len());
        let restored: Vec<usize> = kept().collect();
        for &i in &restored {
            t.fixed[remap[i] as usize] = self.fixed[i];
            t.pmask[remap[i] as usize] = self.pmask[i];
        }
        (t.pread, t.plist) = self.reread(restored.into_iter().map(|i| (remap[i], i)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        t.forms = self.forms.clone();
        t
    }
}

/// A question's names: the tags and params its tapes read, numbered as they
/// first appear in them (instructions, then names provenance stops at, then
/// forms), then in the equations of the derivations they read
/// (Ctx::read_by: whether t is read counts those too), then the context's
/// other ranged params. The model's own numbers
/// say where it declared or imported each, so a cached answer is keyed by
/// these (`Tape::fingerprint`, `Ctx::fingerprint`): an unchanged part keeps
/// its key in any model that imports it, wherever it does. An answer that
/// names a tag or param by the model's number also hashes those (`raw`).
#[derive(Clone, Debug, Default)]
pub struct Names {
    tags: Vec<u32>,
    params: Vec<u32>,
    /// The context's ranged params, the model's numbers, sorted.
    pub ranged: Vec<u32>,
    /// Whether a tape reads t: else the moment is not in the question (a
    /// part that does not move has one answer over a motion).
    pub reads_t: bool,
}

impl Names {
    pub fn of(tapes: &[&Tape], ctx: &crate::eval::Ctx) -> Names {
        fn add(list: &mut Vec<u32>, v: u32) {
            if !list.contains(&v) {
                list.push(v);
            }
        }
        fn walk(t: &Tape, n: &mut Names) {
            for i in &t.ssa {
                match i.op {
                    Op::Tag(g) => add(&mut n.tags, g),
                    Op::Var(Var::Param(k)) => add(&mut n.params, k),
                    Op::Var(Var::T) => n.reads_t = true,
                    _ => {}
                }
            }
            t.opaque.iter().for_each(|&g| add(&mut n.tags, g));
            t.forms.iter().for_each(|f| walk(f, n));
        }
        let mut n = Names { ranged: ctx.p_spans.iter().map(|s| s.k).collect(), ..Names::default() };
        n.ranged.sort_unstable();
        tapes.iter().for_each(|t| walk(t, &mut n));
        ctx.read_by(tapes).1.iter().flat_map(|d| d.tapes()).for_each(|t| walk(t, &mut n));
        ctx.p_spans.iter().for_each(|s| add(&mut n.params, s.k));
        n
    }

    /// Tag t's number in the question.
    pub fn tag(&self, t: u32) -> u32 {
        self.tags.iter().position(|&g| g == t).map_or(u32::MAX, |i| i as u32)
    }

    /// Param k's number in the question.
    pub fn param(&self, k: u32) -> u32 {
        self.params.iter().position(|&p| p == k).map_or(u32::MAX, |i| i as u32)
    }

    /// The params in the question's order (the model's numbers).
    pub fn params(&self) -> &[u32] {
        &self.params
    }

    /// The model's numbers of the tags and params, in the question's order.
    pub fn raw(&self, fp: &mut Fingerprint) {
        fp.word(self.tags.len() as u64).word(self.params.len() as u64);
        self.tags.iter().chain(&self.params).for_each(|&v| _ = fp.word(v as u64));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::{Ctx, eval_interval_clamped, eval_point};
    use crate::interval::I;

    /// A term shared by many is one conjunct, whatever the paths to it: a
    /// sphere's field under twenty emaxes of itself with itself is one term
    /// (walked as a tree, 2^20 of them; thirty levels, 16 GiB).
    #[test]
    fn conjuncts_are_as_many_as_the_instructions_at_most() {
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let (xx, yy, zz) = (g.un(UnOp::Square, x), g.un(UnOp::Square, y), g.un(UnOp::Square, z));
        let xy = g.add(xx, yy);
        let r2 = g.add(xy, zz);
        let r = g.un(UnOp::Sqrt, r2);
        let one = g.constant(1.0);
        let mut f = g.sub(r, one);
        for _ in 0..20 {
            f = g.bin(BinOp::Emax, f, f);
        }
        let t = Tape::new(&g, f);
        let terms = t.conjuncts();
        assert!(terms.len() <= t.len() && terms.len() == 1, "{} terms of a tape {} long", terms.len(), t.len());
    }

    /// A min of a and a - c, c > 0 (a blend below what it rounds, named),
    /// is the blend wherever c's bounds are above 0, though a's and the
    /// blend's bounds overlap (each reads a, here x^2 over a box about
    /// 0).
    #[test]
    fn a_min_of_a_and_a_less_c_is_ordered_by_c() {
        let mut g = Graph::new();
        let (x, y) = (g.var(Var::X), g.var(Var::Y));
        let a = g.un(UnOp::Square, x);
        let yy = g.un(UnOp::Square, y);
        let one = g.constant(1.0);
        let c = g.add(yy, one);
        let less = g.sub(a, c);
        let blend = g.tag("blend", less);
        for f in [g.min(a, blend), g.min(blend, a)] {
            let t = Tape::new(&g, f);
            let mut ch = vec![BOTH; t.len()];
            let b = [I::new(-1.0, 1.0), I::new(-0.5, 0.5), I::point(0.0)];
            let r = crate::eval::eval_interval(&t, &Ctx::new(0.0, vec![]), b, Some(&mut ch));
            // The blend's own bound: x^2 - (y^2 + 1) over the box.
            assert!(r.lo >= -1.25 - 1e-12 && r.hi <= 0.0, "{r:?}");
            let i = t.len() - 1;
            let picked = match ch[i] {
                FIRST => t.ssa[i].a,
                SECOND => t.ssa[i].b,
                _ => panic!("undecided: {r:?}"),
            };
            assert!(matches!(t.ssa[picked as usize].op, Op::Tag(_)), "{:?}", t.ssa[i]);
        }
    }

    /// A tape simplified by its sign choices over a box has its sign (below,
    /// at or above 0) at every point of the box: through mins, maxes, an
    /// emax, a negation and a product with a positive factor, and where a
    /// branch reads 0 throughout a slab (max(0, floor(z + 0.5)), which no
    /// choice may drop as below 0).
    #[test]
    fn sign_choices_keep_the_sign_everywhere() {
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let (quarter, half, zero, two) = (g.constant(0.25), g.constant(0.5), g.constant(0.0), g.constant(2.0));
        let zh = g.add(z, half);
        let step = g.un(UnOp::Floor, zh);
        let plateau = g.max(zero, step);
        let xq = g.sub(x, quarter);
        let left = g.min(xq, plateau);
        let yy = g.mul(y, y);
        let ring = g.sub(yy, quarter);
        let slab = g.sub(z, half);
        let cap = g.bin(BinOp::Emax, ring, slab);
        let out = g.neg(cap);
        let scaled = g.mul(two, out);
        let f = g.max(left, scaled);
        let t = Tape::new(&g, f);
        let ctx = Ctx::new(0.0, vec![]);
        let class = |v: f64| {
            if v < 0.0 {
                -1
            } else if v > 0.0 {
                1
            } else {
                0
            }
        };
        let mut made = 0;
        for i in 0..64 {
            let at = |k: usize| -1.0 + 0.5 * ((i >> (2 * k)) & 3) as f64;
            let b: [I; 3] = std::array::from_fn(|k| I::new(at(k), at(k) + 0.5));
            let mut vals = Vec::new();
            eval_interval_clamped(&t, &ctx, b, &[], &mut vals).unwrap();
            let mut ch = vec![BOTH; t.len()];
            if !t.sign_choices(&vals, &mut ch) {
                continue;
            }
            made += 1;
            let s = t.simplify(&ch);
            assert!(s.len() < t.len());
            for q in 0..216 {
                let p: crate::V3 = std::array::from_fn(|k| b[k].lo + (b[k].hi - b[k].lo) * ((q / [1, 6, 36][k]) % 6) as f64 / 5.0);
                assert_eq!(class(eval_point(&t, &ctx, p)), class(eval_point(&s, &ctx, p)), "at {p:?} in {b:?}");
            }
        }
        assert!(made > 8, "choices made in {made} boxes");
    }
}
