use super::{Graph, Node, NodeId, Var};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

impl Graph {
    /// All nodes reachable from `roots`, sorted ascending (a valid evaluation order).
    pub fn reachable(&self, roots: &[NodeId]) -> Vec<NodeId> {
        self.reach(roots, true, true)
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
                Node::Deriv(..) | Node::CornerReach(..) => ([0; 3], 0),
                node => node.children(),
            };
            stack.extend(ch[..k].iter().filter(|c| !seen.contains(c)));
        }
        let mut out: Vec<NodeId> = seen.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// `reachable`, through a read's value (Node::Read's `a`) or not, and
    /// through a corner test's operands (`corners`) or not.
    pub(super) fn reach(&self, roots: &[NodeId], values: bool, corners: bool) -> Vec<NodeId> {
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
            let node = self.nodes[n as usize];
            let (ch, k) = match node {
                Node::Read(_, _, _, at) if !values => ([at, 0, 0], 1),
                node => node.children(),
            };
            for &c in &ch[..k] {
                if !seen.contains(&c) {
                    stack.push(c);
                }
            }
            if corners && let Node::CornerReach(ci, ..) = node {
                let (a, b) = self.corners[ci as usize];
                for c in [a, b] {
                    if !seen.contains(&c) {
                        stack.push(c);
                    }
                }
            }
        }
        out.sort_unstable();
        out
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

    /// Whether `root` reads a variable `pred` takes, through everything it
    /// is built of. A corner test's operands are read at its point (`at`),
    /// not at x y z: their coordinates count only as the point's do, so a
    /// fillet warped to a fixed point reads no x y z.
    pub fn depends_on(&self, root: NodeId, pred: impl Fn(Var) -> bool) -> bool {
        self.reads_dyn(root, &pred, true)
    }

    /// Whether `root`'s value reads a variable `pred` takes: `depends_on`,
    /// but not through what only bounds it (a Lipschitz promise's frame and
    /// constant, a derivative's mark), so a profile in x and y with a
    /// distance promise does not read z.
    pub fn value_depends_on(&self, root: NodeId, pred: impl Fn(Var) -> bool) -> bool {
        self.value_depends_on_dyn(root, &pred)
    }

    fn value_depends_on_dyn(&self, root: NodeId, pred: &dyn Fn(Var) -> bool) -> bool {
        self.reads_dyn(root, pred, false)
    }

    /// `depends_on` (`bounds`: through a Lipschitz promise's frame and
    /// constant and a derivative's mark) or `value_depends_on` (not).
    fn reads_dyn(&self, root: NodeId, pred: &dyn Fn(Var) -> bool, bounds: bool) -> bool {
        // Memoized post-order pass: bit 0/1/2 = reads X/Y/Z; bit 3 = reads a
        // non-spatial var matching `pred`. Avoids exponential re-traversal
        // through nested Node::CornerReach sub-graphs.
        let mut memo: HashMap<NodeId, u8> = HashMap::default();
        let mut stack = vec![(root, false)];
        while let Some((n, done)) = stack.pop() {
            if memo.contains_key(&n) {
                continue;
            }
            if done {
                let m = match self.nodes[n as usize] {
                    Node::Var(Var::X) => 1,
                    Node::Var(Var::Y) => 2,
                    Node::Var(Var::Z) => 4,
                    Node::Var(v) => u8::from(pred(v)) << 3,
                    Node::Lipschitz(a, ..) | Node::Deriv(a, ..) if !bounds => memo[&a],
                    Node::CornerReach(ci, r, g, at) => {
                        let (a, b) = self.corners[ci as usize];
                        let Node::Frame(fc) = self.nodes[at as usize] else { unreachable!("a corner's coordinates are a frame") };
                        let sub = memo[&a] | memo[&b];
                        let mut out = memo[&r] | memo[&g] | (sub & 8);
                        for (j, &c) in fc.iter().enumerate() {
                            if sub & (1 << j) != 0 {
                                out |= memo[&c];
                            }
                        }
                        out
                    }
                    node => {
                        let (ch, k) = node.children();
                        ch[..k].iter().fold(0u8, |acc, &c| acc | memo[&c])
                    }
                };
                memo.insert(n, m);
            } else {
                stack.push((n, true));
                match self.nodes[n as usize] {
                    Node::Lipschitz(a, ..) | Node::Deriv(a, ..) if !bounds => {
                        if !memo.contains_key(&a) {
                            stack.push((a, false));
                        }
                    }
                    Node::CornerReach(ci, r, g, at) => {
                        let (a, b) = self.corners[ci as usize];
                        let Node::Frame(fc) = self.nodes[at as usize] else { unreachable!("a corner's coordinates are a frame") };
                        for c in [a, b, r, g, fc[0], fc[1], fc[2]] {
                            if !memo.contains_key(&c) {
                                stack.push((c, false));
                            }
                        }
                    }
                    node => {
                        let (ch, k) = node.children();
                        for &c in &ch[..k] {
                            if !memo.contains_key(&c) {
                                stack.push((c, false));
                            }
                        }
                    }
                }
            }
        }
        let m = memo[&root];
        (m & 1 != 0 && pred(Var::X)) || (m & 2 != 0 && pred(Var::Y)) || (m & 4 != 0 && pred(Var::Z)) || (m & 8 != 0)
    }
}
