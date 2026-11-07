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
            Op::Select | Op::Grid(..) | Op::Lipschitz | Op::Frame => ([self.a, self.b, self.c], 3),
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
    /// The instructions that can jump (floor, atan2), in order, which
    /// interval evaluation may take as a choice of two sides
    /// (eval::forked). Derived from `ssa`.
    pub jumps: Vec<u32>,
    /// Where jumps join (Tape::rejoin), derived from `ssa` when a jump is
    /// first straddled.
    joins: std::sync::OnceLock<Joins>,
    /// Each promise's operands (Tape::promise) by its instruction, each
    /// found when first read.
    promises: std::sync::OnceLock<Box<[(u32, std::sync::OnceLock<Promise>)]>>,
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
}

impl Tape {
    /// The tape of `root`: its instructions in one order, operands before
    /// their consumers, depth first from the root, each node's operands in
    /// turn; with the grids it reads (renumbered in that order) and the
    /// names provenance stops at among its own. So the tape of an
    /// expression is the same however the graph holding it was built, and
    /// whatever else it holds (Tape::fingerprint).
    pub fn new(g: &Graph, root: NodeId) -> Tape {
        let mut idx: std::collections::HashMap<NodeId, u32> = std::collections::HashMap::new();
        let mut ids: Vec<NodeId> = Vec::new();
        let mut stack: Vec<(NodeId, usize)> = vec![(root, 0)];
        while let Some((n, k)) = stack.last_mut() {
            let (ch, count) = g.nodes[*n as usize].children();
            if *k < count {
                let c = ch[*k];
                *k += 1;
                if !idx.contains_key(&c) {
                    stack.push((c, 0));
                }
            } else {
                let n = *n;
                stack.pop();
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

    fn from_ssa(ssa: Vec<Ins>) -> Tape {
        let n = ssa.len();
        let uses = |i: usize| Self::reads(&ssa, i);
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
        let (mut fixed, mut pmask) = (vec![false; n], vec![0u64; n]);
        for (i, ins) in ssa.iter().enumerate() {
            let (args, k) = ins.args();
            fixed[i] = match ins.op {
                Op::Const(_) | Op::Var(Var::Param(_)) => true,
                Op::Un(_) | Op::Bin(_) | Op::Select | Op::Tag(_) => args[..k].iter().all(|&a| fixed[a as usize]),
                Op::Var(_) | Op::Grid(..) | Op::Read | Op::Lipschitz | Op::Frame | Op::Deriv(_) | Op::Real => false,
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
        let jumps = (0..n as u32).filter(|&i| matches!(ssa[i as usize].op, Op::Un(UnOp::Floor) | Op::Bin(BinOp::Atan2))).collect();
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
        }
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
            Joins { rejoin }
        })
    }

    /// Lipschitz instruction `i`'s operands as interval evaluation reads
    /// them again (Promise), found once.
    pub fn promise(&self, i: usize) -> &Promise {
        let all = self.promises.get_or_init(|| (0..self.ssa.len() as u32).filter(|&j| matches!(self.ssa[j as usize].op, Op::Lipschitz)).map(|j| (j, Default::default())).collect());
        let k = all.binary_search_by_key(&(i as u32), |p| p.0).expect("a promise");
        all[k].1.get_or_init(|| Promise::of(&self.ssa, i))
    }

    /// Where the bounds of a jump at instruction `j` (a floor or atan2;
    /// eval::forked) can be joined again: the first instruction from j on
    /// through which every later use of j's value passes, and from which the
    /// output is monotone (Joins). The output at the latest.
    pub fn rejoin(&self, j: usize) -> usize {
        self.joins().rejoin[j] as usize
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
        // `jumps`, `joins` and `promises` are derived from `ssa`; `pmask`, `pread` and
        // `plist` enter as what they decide (below).
        let Tape { ssa, code, nslots, grids, opaque, fixed, pmask: _, pread: _, plist: _, slopes, forms, jumps: _, joins: _, promises: _ } = self;
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

    /// The params the tape reads (its forms' included), each once.
    pub fn params_read(&self) -> Vec<u32> {
        let mut out: Vec<u32> = self.ssa.iter().filter_map(|i| if let Op::Var(Var::Param(k)) = i.op { Some(k) } else { None }).collect();
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
        let mut seen: std::collections::HashSet<u32> = roots.iter().map(|&r| r as u32).collect();
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
        // Each root after the first joined by a max, the last the output
        // (one root: its own, ids' last, as it reads every other).
        let mut out = remap[&(roots[0] as u32)];
        for &r in &roots[1..] {
            let n = ssa.len() as u32;
            ssa.push(Ins { op: Op::Bin(BinOp::Max), out: n, a: out, b: remap[&(r as u32)], c: 0 });
            out = n;
        }
        let mut t = Tape::from_ssa(ssa);
        // As this tape has them (a simplified tape's are not recomputable);
        // a joining max is not fixed (bounded as any value is: soundly).
        let joined = roots.len() - 1;
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
                (Op::Select, BOTH) if strip(ins.b as usize) != strip(ins.c as usize) => return true,
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

    /// Whether the output changes sign at most once along unit direction `d`
    /// on every line through a region (from below 0 to not, or the
    /// reverse, the same way on every line), from `choices` and `grads`
    /// found over it (eval_interval's and eval_interval_grad_every's). Then
    /// its zero set there is one sheet: no piece of it lies wholly inside.
    /// Each value is so where its rate along d, carried by the chain rule (a
    /// min's, max's or select's undecided the hull of its branches'),
    /// excludes 0; where it keeps one sign throughout; and where its sign
    /// is that of a value that is so: a min or max (or emax, whose sign is
    /// its max's) of two that change sign the same way (where one is not
    /// below 0 is a union of such, and an intersection), or of one and one
    /// that keeps the sign that leaves the other's (a max with a value
    /// below 0 throughout, a min with one not below 0), and a product with
    /// a value above 0 throughout.
    pub fn monotone(&self, choices: &[u8], grads: &[crate::interval::G], d: crate::V3) -> bool {
        use crate::interval::{I, add_o, hull, mul_o, sub_o};
        let n = self.ssa.len();
        let mut r: Vec<I> = Vec::with_capacity(n);
        // Each value's way of changing sign along d: 1 (below 0, then
        // not), -1 (not, then below), 0 (unknown).
        let mut way: Vec<i8> = Vec::with_capacity(n);
        let v = |j: u32| grads[j as usize].v;
        // Inside throughout, or nowhere (crate::sign).
        let below = |j: u32| v(j).all_inside();
        let not_below = |j: u32| v(j).none_inside();
        for i in 0..n {
            let ins = &self.ssa[i];
            let at = |j: u32| r[j as usize];
            let rate = match (ins.op, choices[i]) {
                (Op::Var(Var::X), _) => I::point(d[0]),
                (Op::Var(Var::Y), _) => I::point(d[1]),
                (Op::Var(Var::Z), _) => I::point(d[2]),
                (Op::Var(_) | Op::Const(_), _) => I::point(0.0),
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
                    match still || strip(ins.b as usize) == strip(ins.c as usize) {
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
            let by_sign = match (ins.op, choices[i]) {
                (Op::Tag(_) | Op::Lipschitz | Op::Real, _) => w(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => w(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => w(ins.b),
                (Op::Select, FIRST) => w(ins.b),
                (Op::Select, SECOND) => w(ins.c),
                (Op::Bin(BinOp::Max | BinOp::Emax), BOTH) if below(ins.a) => w(ins.b),
                (Op::Bin(BinOp::Max | BinOp::Emax), BOTH) if below(ins.b) => w(ins.a),
                (Op::Bin(BinOp::Min), BOTH) if not_below(ins.a) => w(ins.b),
                (Op::Bin(BinOp::Min), BOTH) if not_below(ins.b) => w(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), BOTH) => alike(w(ins.a), w(ins.b)),
                (Op::Bin(BinOp::Mul), _) if v(ins.b).lo > 0.0 => w(ins.a),
                (Op::Bin(BinOp::Mul), _) if v(ins.a).lo > 0.0 => w(ins.b),
                _ => 0,
            };
            let kept = below(i as u32) || not_below(i as u32);
            way.push(if kept {
                2
            } else if by_rate != 0 {
                by_rate
            } else {
                by_sign
            });
        }
        way[n - 1] != 0
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

    /// Drop branches of min/max/select that interval evaluation proved can't win.
    pub fn simplify(&self, choices: &[u8]) -> Tape {
        let n = self.ssa.len();
        let pick = |i: usize| -> Option<u32> {
            let ins = &self.ssa[i];
            match (ins.op, choices[i]) {
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => Some(ins.a),
                (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => Some(ins.b),
                (Op::Select, FIRST) => Some(ins.b),
                (Op::Select, SECOND) => Some(ins.c),
                _ => None,
            }
        };

        let mut live = vec![false; n];
        live[n - 1] = true;
        for i in (0..n).rev() {
            if !live[i] {
                continue;
            }
            if let Some(a) = pick(i) {
                live[a as usize] = true;
            } else {
                let (args, k) = self.ssa[i].args();
                for &a in &args[..k] {
                    live[a as usize] = true;
                }
            }
        }

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
/// forms), then the context's other ranged params. The model's own numbers
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
