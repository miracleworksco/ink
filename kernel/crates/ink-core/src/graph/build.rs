use super::{BinOp, Graph, Node, NodeId, UnOp, Var, bin_f64, exact_bin, exact_un, lost, un_f64, undefined_bin, undefined_un, wild};

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(super) enum Key {
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
    CornerReach(u32, NodeId, NodeId, NodeId),
    /// A constant a dropped term made (Graph::exact_number), by its bits.
    Exact(u64),
    /// A constant whose real value is not its float (Graph::real_constant), by a hash of its content.
    Inexact(u64),
}

impl Graph {
    pub fn new() -> Self {
        Self::default()
    }

    pub(super) fn intern(&mut self, n: Node) -> NodeId {
        // An inexact constant (`real_constant`) is its own node only in numbers of constants and params
        // (a param's arithmetic on a written 16.7, which a value claim reads in real numbers, Doc's
        // exact_value): as anything else's operand it is the shared constant of its float, so every
        // other node is the one it would be without it. Either way its float is the same, and so is
        // the expression's form (`shape_of`).
        let (ch, k) = n.children();
        let fixed = matches!(n, Node::Un(..) | Node::Bin(..) | Node::Select(..) | Node::Tag(..)) && ch[..k].iter().all(|&c| self.fixed[c as usize]);
        let n = if !fixed && ch[..k].iter().any(|c| self.reals.contains_key(c)) {
            let mut plain = |c: NodeId| match self.reals.contains_key(&c) {
                true => {
                    let v = self.constant_value(c).expect("a constant");
                    self.number(v)
                }
                false => c,
            };
            match n {
                Node::Const(_) | Node::Var(_) => n,
                Node::Un(op, a) => Node::Un(op, plain(a)),
                Node::Bin(op, a, b) => Node::Bin(op, plain(a), plain(b)),
                Node::Select(c, a, b) => Node::Select(plain(c), plain(a), plain(b)),
                Node::Tag(t, a) => Node::Tag(t, plain(a)),
                Node::Grid(g, m, at) => Node::Grid(g, m, at.map(&mut plain)),
                Node::Read(g, p, a, at) => Node::Read(g, p, plain(a), plain(at)),
                Node::Lipschitz(a, k, f) => Node::Lipschitz(plain(a), plain(k), plain(f)),
                Node::Frame(at) => Node::Frame(at.map(&mut plain)),
                Node::Deriv(a, f, k) => Node::Deriv(plain(a), plain(f), k),
                Node::Real(a) => Node::Real(plain(a)),
                Node::Within(l, h, v, id) => Node::Within(plain(l), plain(h), plain(v), id),
                Node::CornerReach(ci, r, g, at) => Node::CornerReach(ci, plain(r), plain(g), plain(at)),
            }
        } else {
            n
        };
        let (ch, k) = n.children();
        if !self.domains.original.is_empty() && ch[..k].iter().any(|c| self.domains.original.contains_key(c)) {
            // Aliases distinguish design reads, never mathematical variables.
            // Every numerical node and tape still has one copy of each input.
            let value = |c| self.numeric(c);
            let numeric = match n {
                Node::Const(_) | Node::Var(_) => n,
                Node::Un(op, a) => Node::Un(op, value(a)),
                Node::Bin(op, a, b) => Node::Bin(op, value(a), value(b)),
                Node::Select(c, a, b) => Node::Select(value(c), value(a), value(b)),
                Node::Tag(t, a) => Node::Tag(t, value(a)),
                Node::Grid(g, m, at) => Node::Grid(g, m, at.map(value)),
                Node::Read(g, p, a, at) => Node::Read(g, p, value(a), value(at)),
                Node::Lipschitz(a, k, f) => Node::Lipschitz(value(a), value(k), value(f)),
                Node::Frame(at) => Node::Frame(at.map(value)),
                Node::Deriv(a, f, k) => Node::Deriv(value(a), value(f), k),
                Node::Real(a) => Node::Real(value(a)),
                Node::Within(l, h, v, id) => Node::Within(value(l), value(h), value(v), id),
                Node::CornerReach(ci, r, g, at) => Node::CornerReach(ci, value(r), value(g), value(at)),
            };
            let domain = self.domains.node(n);
            let value = self.intern_value(numeric);
            return self.with_domain(value, domain);
        }
        self.intern_value(n)
    }

    pub(super) fn numeric(&self, n: NodeId) -> NodeId {
        self.domains.original.get(&n).copied().unwrap_or(n)
    }

    fn intern_value(&mut self, n: Node) -> NodeId {
        let key = Self::key_of(n);
        self.intern_keyed(n, key)
    }

    pub(super) fn key_of(n: Node) -> Key {
        match n {
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
            Node::CornerReach(ci, r, g, at) => Key::CornerReach(ci, r, g, at),
        }
    }

    pub(super) fn intern_keyed(&mut self, n: Node, key: Key) -> NodeId {
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
        let partial = matches!(n, Node::Bin(BinOp::Div | BinOp::Pow, ..) | Node::Un(UnOp::Ln | UnOp::Sqrt, _)) || {
            let (ch, k) = n.children();
            ch[..k].iter().any(|&c| self.partial[c as usize])
        };
        self.partial.push(partial);
        let marked = match n {
            Node::Var(Var::Param(k)) => self.markers.contains(&k),
            // Its fields are kept beside it (`corners`), not among its children.
            Node::CornerReach(ci, r, g, at) => {
                let (a, b) = self.corners[ci as usize];
                [a, b, r, g, at].iter().any(|&c| self.marked[c as usize])
            }
            _ => {
                let (ch, k) = n.children();
                ch[..k].iter().any(|&c| self.marked[c as usize])
            }
        };
        self.marked.push(marked);
        let domain = self.domains.node(n);
        self.domains.of.push(domain);
        self.nodes.push(n);
        self.shapes.push(shape);
        self.fixed.push(fixed);
        self.index.insert(key, id);
        id
    }

    /// A constant: a number the model is built from (written, read from a
    /// file, made in Rust), noted where outside the range.
    pub fn constant(&mut self, v: f64) -> NodeId {
        self.note(v, wild(v));
        self.number(v)
    }

    pub(super) fn number(&mut self, v: f64) -> NodeId {
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
    pub(super) fn folded(&mut self, v: f64, args: &[NodeId], vals: &[f64], scales: bool) -> NodeId {
        if args.iter().all(|a| self.owned.contains(a)) {
            return self.own(v);
        }
        self.note(v, lost(v, vals, scales));
        self.number(v)
    }

    pub fn var_node(&self, v: Var) -> Option<NodeId> {
        self.index.get(&Key::Var(v)).copied()
    }

    pub fn var(&mut self, v: Var) -> NodeId {
        self.intern(Node::Var(v))
    }

    /// Param `k`, a marker (Graph::marked): what reads it is never folded
    /// so as not to (an estimate's guard: a claim reading it is to know).
    /// Made a marker before anything reads it, so every node that does is
    /// marked.
    pub fn marker(&mut self, k: u32) -> NodeId {
        debug_assert!(!self.index.contains_key(&Key::Var(Var::Param(k))), "param {k} was read before it was a marker");
        self.markers.insert(k);
        self.var(Var::Param(k))
    }

    /// A literal or an expression of fixed scalar inputs, evaluated in floats
    /// unless one of those inputs is ranged. Real constructions keep bounds.
    pub fn fixed_scalar(&self, id: NodeId) -> bool {
        self.fixed[id as usize]
    }

    pub fn as_const(&self, id: NodeId) -> Option<f64> {
        if self.domains.of[id as usize] != 0 {
            return None;
        }
        self.constant_value(id)
    }

    /// The folded numerical coefficient, independent of checks required to
    /// use it. For inspecting a Jacobian's sparsity, for example. Use
    /// `as_const` when replacing a design expression with a number.
    pub fn constant_value(&self, id: NodeId) -> Option<f64> {
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
        let value = self.un_value(op, self.numeric(a));
        self.keeping_inputs(value, &[a])
    }

    /// The 0 that dropping `of` leaves (0 times it, it less itself): where one reads x, y, z or
    /// t, an exact 0 (`exact_number`), so what folds with it is exact or a band: 1e16 + (x - x +
    /// 0.5) less 1e16 is 0.5, not 1e16 less 1e16 in floats, and so of t. Of constants and params
    /// as written, the constant 0.
    fn zero_of(&mut self, of: &[NodeId]) -> NodeId {
        match of.iter().all(|&n| self.fixed[n as usize]) {
            true => self.constant(0.0),
            false => self.exact_number(0.0),
        }
    }

    /// Whether dropping `n` from an expression (0 times it, it less itself)
    /// would hide a marker it reads (`marker`), or a number that may have
    /// no result: a constants-and-params value with such an operation in it.
    pub(super) fn hides(&self, n: NodeId) -> bool {
        self.marked[n as usize] || (self.strict && self.fixed[n as usize] && self.partial[n as usize])
    }

    /// Whether a number made of constants and params has an operation in
    /// it that gives 0 where it has no result (`partial`): ask each such
    /// operation over its operands' ranges (`undefined_un_over`,
    /// `undefined_bin_over`; a sqrt, whether its operand is below 0).
    pub fn may_lack_a_result(&self, n: NodeId) -> bool {
        self.partial[n as usize]
    }

    pub(super) fn un_value(&mut self, op: UnOp, a: NodeId) -> NodeId {
        if let Some(v) = self.constant_value(a)
            && undefined_un(op, v).is_none()
            && !(self.strict && op == UnOp::Sqrt && v < 0.0)
        {
            if self.is_exact(a) {
                if let Some(x) = exact_un(op, v) {
                    return self.exact_number(x);
                }
                let band = crate::interval::un(op, crate::interval::I::point(v));
                if let Some(n) = self.exact_band(band, un_f64(op, v), Node::Un(op, a)) {
                    return n;
                }
                let a = self.real(a);
                return self.intern(Node::Un(op, a));
            }
            let r = self.folded(un_f64(op, v), &[a], &[v], matches!(op, UnOp::Square | UnOp::Exp));
            return self.fold_real(r, Ok(op), &[a]).unwrap_or(r);
        }
        match (op, self.nodes[a as usize]) {
            (UnOp::Neg, Node::Un(UnOp::Neg, inner)) => return inner,
            (UnOp::Abs, Node::Un(UnOp::Abs, _)) => return a,
            _ => {}
        }
        self.intern(Node::Un(op, a))
    }

    pub fn bin(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        let value = self.bin_value(op, self.numeric(a), self.numeric(b));
        self.keeping_inputs(value, &[a, b])
    }

    pub(super) fn bin_value(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        use BinOp::*;
        let (ca, cb) = (self.constant_value(a), self.constant_value(b));
        let (ea, eb) = (self.is_exact(a), self.is_exact(b));
        if let (Some(x), Some(y)) = (ca, cb) {
            if (ea || eb) && undefined_bin(op, x, y).is_none() {
                // Arithmetic on a term of x, y, z or t: exact, or the band its rounding leaves.
                if let Some(v) = exact_bin(op, x, y) {
                    return self.exact_number(v);
                }
                let (band, _) = crate::interval::bin(op, crate::interval::I::point(x), crate::interval::I::point(y));
                if let Some(n) = self.exact_band(band, bin_f64(op, x, y), Node::Bin(op, a, b)) {
                    return n;
                }
                let a = self.real(a);
                return self.intern(Node::Bin(op, a, b));
            }
            return match undefined_bin(op, x, y) {
                None => {
                    let r = self.folded(bin_f64(op, x, y), &[a, b], &[x, y], matches!(op, Mul | Div | Pow | Emax));
                    self.fold_real(r, Err(op), &[a, b]).unwrap_or(r)
                }
                Some(_) => self.intern(Node::Bin(op, a, b)),
            };
        }
        // An exact 0 adds nothing and keeps what it meets real (a param beside it is read as a real
        // number); times what has a result it is itself. An exact constant beside a param reads it
        // as a real number, so floats never decide their sum.
        let (za, zb) = (ea && ca == Some(0.0), eb && cb == Some(0.0));
        match op {
            Add | Sub if zb => return self.real(a),
            Add if za => return self.real(b),
            Sub if za => {
                let b = self.real(b);
                return self.un(UnOp::Neg, b);
            }
            Mul if (za && !self.hides(b)) || (zb && !self.hides(a)) => return if za { a } else { b },
            _ => {}
        }
        let (a, b) = match (ea, eb) {
            (true, false) if self.fixed[b as usize] => (a, self.real(b)),
            (false, true) if self.fixed[a as usize] => (self.real(a), b),
            _ => (a, b),
        };
        // An inexact constant is no identity: times a written 1.00000000000000001 is not the number itself.
        let plain = |n: NodeId| self.constant_value(n).filter(|_| !self.reals.contains_key(&n));
        let (ca, cb) = (plain(a), plain(b));
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
                if a == b && !self.hides(a) {
                    return self.zero_of(&[a]);
                }
            }
            Mul => {
                if ca == Some(1.0) {
                    return b;
                }
                if cb == Some(1.0) {
                    return a;
                }
                if ca == Some(0.0) && !self.hides(b) {
                    return self.zero_of(&[b]);
                }
                if cb == Some(0.0) && !self.hides(a) {
                    return self.zero_of(&[a]);
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
        if let Some(v) = self.constant_value(c) {
            let chosen = if v < 0.0 { a } else { b };
            return self.keeping_inputs(chosen, &[c, chosen]);
        }
        let value = self.select_value(self.numeric(c), self.numeric(a), self.numeric(b));
        self.keeping_inputs(value, &[c, a, b])
    }

    fn select_value(&mut self, c: NodeId, a: NodeId, b: NodeId) -> NodeId {
        // Either way the same, unless the condition reads a marker (`hides`).
        if a == b && !self.hides(c) {
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

    /// Promise that `a` is k-Lipschitz in x, y, z (identity in value). A k
    /// that is not positive promises nothing.
    pub fn lipschitz(&mut self, a: NodeId, k: NodeId) -> NodeId {
        let value = self.lipschitz_value(a, k);
        self.keeping_inputs(value, &[a, k])
    }

    fn lipschitz_value(&mut self, a: NodeId, k: NodeId) -> NodeId {
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

    /// Exact ball-clearance corner test (`Node::CornerReach`): evaluates `a` and `b`
    /// at the centres of radius-`r` balls along each other's gradients, active where
    /// `guard < 0` (0 otherwise), without copying `a` and `b` into the tape.
    pub fn corner_reach(&mut self, a: NodeId, b: NodeId, r: NodeId, guard: NodeId) -> NodeId {
        let (an, bn) = (self.numeric(a), self.numeric(b));
        let ci = self.intern_corner(an, bn);
        let coords = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let at = self.intern(Node::Frame(coords));
        let node = self.intern(Node::CornerReach(ci, r, guard, at));
        self.keeping_inputs(node, &[node, a, b])
    }

    pub(super) fn intern_corner(&mut self, a: NodeId, b: NodeId) -> u32 {
        if let Some(&ci) = self.corner_index.get(&(a, b)) {
            return ci;
        }
        let ci = self.corners.len() as u32;
        self.corners.push((a, b));
        self.corner_index.insert((a, b), ci);
        ci
    }
}
