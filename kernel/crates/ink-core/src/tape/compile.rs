//! Graph lowering, register allocation, and metadata computed from SSA.

use super::{Ins, Op, Tape};
use crate::graph::{BinOp, Fingerprint, Graph, Node, NodeId, UnOp, Var};

static BUILT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

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
                Node::CornerReach(ci, r, guard, at) => {
                    let (a, b) = g.corners[ci as usize];
                    let mut form_of = |node: NodeId| {
                        *form_index.entry(node).or_insert_with(|| {
                            forms.push(std::sync::Arc::new(Tape::new(g, node)));
                            (forms.len() - 1) as u32
                        })
                    };
                    let fa = form_of(a);
                    let fb = form_of(b);
                    Ins { op: Op::CornerReach(fa, fb), out, a: m(r), b: m(at), c: m(guard) }
                }
            };
            ssa.push(ins);
        }
        let mut t = Self::from_ssa(ssa);
        (t.fixed, t.pmask, t.pread, t.plist) = Self::fixed_of(&t.ssa);
        t.forms = forms;
        t.grids = grids.iter().map(|&gi| g.grids[gi as usize].clone()).collect();
        t.opaque = g.opaque.iter().copied().filter(|o| tags.contains(o)).collect();
        t
    }

    /// What instruction i of `ssa` needs live: its arguments, and a read or
    /// corner test its frame's coordinates too (Op::Read, Op::CornerReach).
    ///
    /// The array is built whole, from values: written a slice at a time,
    /// then returned, it was read back wider than written, which stalled
    /// (a tenth of the Falcon swerve check's time, in tape compilation).
    #[inline(always)]
    pub fn reads(ssa: &[Ins], i: usize) -> ([u32; 6], usize) {
        // args() leaves 0 past its count, as the array's tail is.
        let ([a, b, c], k) = ssa[i].args();
        match ssa[i].op {
            Op::Read | Op::CornerReach(..) => {
                let ([x, y, z], _) = ssa[ssa[i].b as usize].args();
                if k == 2 { ([a, b, x, y, z, 0], 5) } else { ([a, b, c, x, y, z], 6) }
            }
            _ => ([a, b, c, 0, 0, 0], k),
        }
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

    /// Of a subtraction `i` that is a remainder, a - p floor(a / p) (std's
    /// mod, through names): p. In real numbers its value lies between 0
    /// and p, which bounds of it, reading a twice, lose where the floor
    /// takes several values.
    pub fn remainder_of(ssa: &[Ins], i: usize) -> Option<u32> {
        let strip = |mut j: u32| {
            while let Op::Tag(_) = ssa[j as usize].op {
                j = ssa[j as usize].a;
            }
            j
        };
        let ins = ssa[i];
        let Op::Bin(BinOp::Sub) = ins.op else { return None };
        let m = ssa[strip(ins.b) as usize];
        let Op::Bin(BinOp::Mul) = m.op else { return None };
        [(m.a, m.b), (m.b, m.a)].into_iter().find_map(|(p, f)| {
            let f = ssa[strip(f) as usize];
            let Op::Un(UnOp::Floor) = f.op else { return None };
            let q = ssa[strip(f.a) as usize];
            (matches!(q.op, Op::Bin(BinOp::Div)) && strip(q.a) == strip(ins.a) && strip(q.b) == strip(p)).then_some(p)
        })
    }

    /// Tapes built in this process (`Tape::new`, a cone, a simplified
    /// tape): what a query's repeated derivations cost (`kept`).
    pub fn built() -> usize {
        BUILT.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// What instruction i of `ssa` reads as interval evaluation runs it:
    /// `reads`, c of a min or max of a and a + c (Tape::apart), whose
    /// bound reads c's slot to order the two, and p of a remainder
    /// (Tape::remainder), whose bound reads p's.
    pub(super) fn uses(ssa: &[Ins], i: usize) -> ([u32; 6], usize) {
        let (mut u, k) = Self::reads(ssa, i);
        match Self::apart_of(ssa, i).map(|(c, _)| c).or_else(|| Self::remainder_of(ssa, i)) {
            Some(c) => {
                u[k] = c;
                (u, k + 1)
            }
            None => (u, k),
        }
    }

    /// The tape of `ssa`: its slots, and what follows from its
    /// instructions alone. What is fixed (`fixed`, `pmask`, `pread`,
    /// `plist`) is the model's: the caller's to set (`fixed_of`, for a tape
    /// of the model; a tape derived from another keeps the other's).
    pub(super) fn from_ssa(ssa: Vec<Ins>) -> Tape {
        BUILT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let n = ssa.len();
        // One walk finds each value's last use as `uses` reads it (each
        // slot is kept until it), each min or max of a and a + c with its
        // c and each remainder with its p (by instruction, their slots set
        // once known), and the instructions that can jump: floors, atan2s
        // and selects whose conditions read no x, y or z (those may switch
        // at a moment of the scope, eval::moment, where scope_run forks
        // them).
        let (mut last, mut space) = (vec![0u32; n], vec![false; n]);
        let (mut apart, mut remainders, mut jumps) = (Vec::new(), Vec::new(), Vec::new());
        let mut derivs = false;
        for (i, ins) in ssa.iter().enumerate() {
            let (mut args, mut k) = Self::reads(&ssa, i);
            space[i] = matches!(ins.op, Op::Var(Var::X | Var::Y | Var::Z)) || args[..k].iter().any(|&a| space[a as usize]);
            match ins.op {
                Op::Un(UnOp::Floor) | Op::Bin(BinOp::Atan2) => jumps.push(i as u32),
                Op::Select if !space[ins.a as usize] => jumps.push(i as u32),
                Op::Deriv(d) => derivs |= d < 3,
                _ => {}
            }
            let order = match Self::apart_of(&ssa, i) {
                Some((c, neg)) => {
                    apart.push((i as u32, c, neg));
                    Some(c)
                }
                None => Self::remainder_of(&ssa, i).inspect(|&p| remainders.push((i as u32, p))),
            };
            if let Some(c) = order {
                args[k] = c;
                k += 1;
            }
            for &a in &args[..k] {
                last[a as usize] = i as u32;
            }
        }
        last[n - 1] = u32::MAX;

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
            let (args, k) = Self::uses(&ssa, i);
            for j in 0..k {
                let a = args[j] as usize;
                if last[a] == i as u32 && !args[..j].contains(&args[j]) {
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
        for e in &mut apart {
            e.1 = slot[e.1 as usize];
        }
        for e in &mut remainders {
            e.1 = slot[e.1 as usize];
        }
        let mut slopes = Vec::new();
        if derivs {
            slopes = vec![false; n];
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
        }
        Tape {
            ssa,
            code,
            nslots: nslots as usize,
            grids: Vec::new(),
            opaque: Vec::new(),
            fixed: Vec::new(),
            pmask: Vec::new(),
            pread: Vec::new(),
            plist: Vec::new(),
            slopes,
            forms: Vec::new(),
            jumps,
            joins: Default::default(),
            promises: Default::default(),
            axes: Default::default(),
            levels: Default::default(),
            apart,
            remainders,
        }
    }

    /// Of a tape of the model, `ssa`: per instruction, whether it is fixed
    /// (Tape::fixed), the bits of the params it reads (`pmask`), and which
    /// those are (`pread`, `plist`).
    fn fixed_of(ssa: &[Ins]) -> (Vec<bool>, Vec<u64>, Vec<(u32, u32, u32)>, Vec<u32>) {
        let n = ssa.len();
        let (mut fixed, mut pmask) = (vec![false; n], vec![0u64; n]);
        for (i, ins) in ssa.iter().enumerate() {
            let (args, k) = ins.args();
            fixed[i] = match ins.op {
                Op::Const(_) | Op::Var(Var::Param(_)) => true,
                Op::Un(_) | Op::Bin(_) | Op::Select | Op::Tag(_) => args[..k].iter().all(|&a| fixed[a as usize]),
                Op::Var(_) | Op::Grid(..) | Op::Read | Op::Lipschitz | Op::Frame | Op::Deriv(_) | Op::Real | Op::Within(_) | Op::CornerReach(..) => false,
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
        (fixed, pmask, pread, plist)
    }
}
