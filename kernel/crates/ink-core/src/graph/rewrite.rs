use super::{BinOp, Graph, Node, NodeId, UnOp, Var};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};

impl Graph {
    /// Substitute variables (this is `warp`). Lipschitz promises follow
    /// where they still hold (`rebuild`).
    pub fn subst(&mut self, root: NodeId, map: &[(Var, NodeId)]) -> NodeId {
        let nodes: HashMap<NodeId, NodeId> = map.iter().map(|(v, n)| (self.intern(Node::Var(*v)), *n)).collect();
        let value = self.rebuild(root, &nodes);
        self.mapped_domain(root, value, map)
    }

    /// `root` with the nodes `map` names replaced, rebuilt above them (a
    /// pair's placements relative to a moving frame: ink-lang's
    /// motion::moving_frame).
    /// Lipschitz promises follow where they still hold (`rebuild`).
    pub fn replace(&mut self, root: NodeId, map: &std::collections::HashMap<NodeId, NodeId>) -> NodeId {
        let value = self.rebuild(root, map);
        self.keeping_inputs(value, &[value, root])
    }

    /// `replace`, where each node is replaced by one equal to it in real numbers at every point and
    /// every value of t and the params (a sum with what cancels in it cancelled): the field is the
    /// one written, so every promise above holds of it still.
    pub fn replace_equal(&mut self, root: NodeId, map: &std::collections::HashMap<NodeId, NodeId>) -> NodeId {
        let value = self.rebuild_as(root, map, false, None, true);
        self.keeping_inputs(value, &[value, root])
    }

    /// `root` without its names: each tag replaced by the field it names,
    /// every other node kept (the same values, derivatives and promises),
    /// for a reader of values alone, which then sees one field however it
    /// was named (a package's solve that reads no face region's provenance).
    pub fn untagged(&mut self, root: NodeId) -> NodeId {
        let value = self.rebuild_as(root, &HashMap::default(), true, None, false);
        self.keeping_inputs(value, &[root])
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
    pub(super) fn rebuild<S: std::hash::BuildHasher>(&mut self, root: NodeId, map: &std::collections::HashMap<NodeId, NodeId, S>) -> NodeId {
        self.rebuild_as(root, map, false, None, false)
    }

    /// `root` with each of `maxes` (`max` nodes, with the operands under
    /// them grouped) written as the `emax` of each group's max: the same
    /// sign, another value where groups are positive (distance::exact,
    /// which groups operands over square directions). A promise above one
    /// is dropped: it held of the max.
    pub fn with_emax(&mut self, root: NodeId, maxes: &HashMap<NodeId, Vec<Vec<NodeId>>>) -> NodeId {
        let value = self.rebuild_as(root, &HashMap::default(), false, Some(maxes), false);
        self.keeping_inputs(value, &[value, root])
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
        equal: bool,
    ) -> NodeId {
        let is_spatial_var = |id: NodeId| matches!(self.nodes[id as usize], Node::Var(Var::X | Var::Y | Var::Z));
        let spatial_only = !untag && emax.is_none() && map.keys().all(|&from| is_spatial_var(from));
        let has_spatial = map.keys().any(|&from| is_spatial_var(from));
        let non_spatial_map: Option<HashMap<NodeId, NodeId>> =
            (!spatial_only && has_spatial).then(|| map.iter().filter(|&(&k, _)| !is_spatial_var(k)).map(|(&k, &v)| (k, v)).collect());
        let ids = self.reach(&[root], false, !spatial_only && !has_spatial);
        let mut memo: HashMap<NodeId, NodeId> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        // The nodes reading a replacement that breaks promises.
        let mut breaks: HashSet<NodeId> = HashSet::default();
        // A coordinate replaced by what jumps (a select, a floor: a warp that folds): a promise whose
        // k reads x, y or z is not kept above it, since its k is read at the warp's image and a way
        // between two points' images need not lie in it (Promise::leaps, for a tape's own frame).
        let jumps = map.iter().any(|(&from, &to)| {
            matches!(self.nodes[from as usize], Node::Var(Var::X | Var::Y | Var::Z))
                && self.reachable(&[to]).iter().any(|&n| matches!(self.nodes[n as usize], Node::Select(..) | Node::Un(UnOp::Floor, _)))
        });
        for &id in &ids {
            if let Some(&to) = map.get(&id).or_else(|| self.domains.original.get(&id).and_then(|n| map.get(n))) {
                let kept = match self.nodes[id as usize] {
                    Node::Var(Var::X | Var::Y | Var::Z) => true,
                    Node::Var(_) => !self.depends_on(to, |v| matches!(v, Var::X | Var::Y | Var::Z)),
                    _ => false,
                };
                if !kept && !equal {
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
                Node::CornerReach(ci, r, g, at) => {
                    let (r, g, at) = (memo[&r], memo[&g], memo[&at]);
                    let ci = if spatial_only {
                        ci
                    } else if let Some(nsm) = &non_spatial_map {
                        let (a, b) = self.corners[ci as usize];
                        let a = self.rebuild_as(a, nsm, untag, emax, equal);
                        let b = self.rebuild_as(b, nsm, untag, emax, equal);
                        let (a, b) = (self.numeric(a), self.numeric(b));
                        self.intern_corner(a, b)
                    } else {
                        let (a, b) = self.corners[ci as usize];
                        let (a, b) = (self.numeric(memo[&a]), self.numeric(memo[&b]));
                        self.intern_corner(a, b)
                    };
                    self.intern(Node::CornerReach(ci, r, g, at))
                }
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
                    let leaps = jumps && self.depends_on(k, |v| matches!(v, Var::X | Var::Y | Var::Z));
                    let (a, k, f) = (memo[&a], memo[&k], memo[&f]);
                    if self.as_const(a).is_some() || breaks.contains(&id) || leaps { a } else { self.intern(Node::Lipschitz(a, k, f)) }
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
                    } else if k == 3 && self.depends_on(f, |v| matches!(v, Var::T | Var::Param(_))) {
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
}
