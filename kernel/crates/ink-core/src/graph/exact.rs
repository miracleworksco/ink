use super::{BinOp, Fingerprint, Graph, Key, Node, NodeId, UnOp, bin_f64, exact_bin, exact_un, un_f64};

impl Graph {
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

    /// Evaluate a constant turn once. Keep the real angle for identities
    /// such as sin(a)^2 + cos(a)^2, and its enclosure for numeric evaluation.
    pub fn turn_components(&mut self, angle: NodeId) -> [NodeId; 2] {
        use crate::interval::{I, un};
        let Some(degrees) = self.as_const(angle).filter(|d| d.is_finite() && *d != 0.0) else {
            let radians = self.un_real(UnOp::Rad, angle);
            return [self.un_real(UnOp::Sin, radians), self.un_real(UnOp::Cos, radians)];
        };
        // Opposite turns share the same cosine and opposite sines.
        let magnitude = degrees.abs();
        let positive = self.constant(magnitude);
        let radians = self.un_real(UnOp::Rad, positive);
        let angle_bound = un(UnOp::Rad, I::point(magnitude));
        let [sin, cos] = [UnOp::Sin, UnOp::Cos].map(|op| {
            // The 30-60-90 and 45-45-90 triangles give tighter (sometimes
            // exact) coefficients for the turns common in circular patterns.
            let (bound, value) = match (magnitude, op) {
                (30.0, UnOp::Sin) => (I::point(0.5), 0.5),
                (30.0, UnOp::Cos) => (un(UnOp::Sqrt, I::point(0.75)), 0.75f64.sqrt()),
                (45.0, _) => (un(UnOp::Sqrt, I::point(0.5)), 0.5f64.sqrt()),
                _ => (un(op, angle_bound), un_f64(op, un_f64(UnOp::Rad, magnitude))),
            };
            let id = Fingerprint::default().str("constant turn").word(op as u64).f64(magnitude).value() as u64;
            // Keep even the singleton's trig identity for symbolic proofs.
            let (l, h, v) = (self.constant(bound.lo), self.constant(bound.hi), self.constant(value));
            let n = self.intern(Node::Within(l, h, v, id));
            self.constructions.insert(n, Node::Un(op, radians));
            n
        });
        [if degrees < 0.0 { self.neg(sin) } else { sin }, cos]
    }

    /// A stored, spatially constant coefficient: its point value for
    /// drawing and its full enclosure for proofs.
    pub fn enclosed_constant(&self, n: NodeId) -> Option<(f64, crate::interval::I)> {
        use crate::interval::I;
        match self.nodes[n as usize] {
            Node::Const(v) => Some((v, I::point(v))),
            Node::Within(l, h, v, _) => Some((self.as_const(v)?, I::new(self.as_const(l)?, self.as_const(h)?))),
            Node::Un(UnOp::Neg, a) => self.enclosed_constant(a).map(|(v, b)| (-v, I::new(-b.hi, -b.lo))),
            _ => None,
        }
    }

    /// Real arithmetic in a frame. Keep a product or sum of enclosed
    /// constants constant too, with its symbolic form for rigid identities:
    /// even where its bound is one float (sin 30 sin 30 is 0.25), so two
    /// composed turns keep sin^2 + cos^2 = 1.
    pub fn frame_bin(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        let n = self.bin_real(op, a, b);
        // Already one (-1 times -sin 30 is sin 30): kept, its form with it.
        let enclosed = |m: NodeId| match self.nodes[m as usize] {
            Node::Un(UnOp::Neg, i) => matches!(self.nodes[i as usize], Node::Within(..)),
            m => matches!(m, Node::Within(..)),
        };
        if n == a || n == b || self.as_const(n).is_some() || enclosed(n) {
            return n;
        }
        let (Some((x, u)), Some((y, v))) = (self.enclosed_constant(a), self.enclosed_constant(b)) else { return n };
        let value = bin_f64(op, x, y);
        let bound = crate::interval::bin(op, u, v).0;
        if !(bound.lo.is_finite() && bound.hi.is_finite() && bound.contains(value)) {
            return n;
        }
        let id = Fingerprint::default().str("constant frame coefficient").key(self.shape(n)).value() as u64;
        let made = self.nodes[n as usize];
        // A singleton made of a turn's coefficients keeps its construction
        // too (as `turn_components`' own do); of plain bands, it is its float.
        let symbolic = |m: NodeId| {
            let m = match self.nodes[m as usize] {
                Node::Un(UnOp::Neg, i) => i,
                _ => m,
            };
            self.constructions.contains_key(&self.numeric(m))
        };
        let coefficient = match bound.lo == bound.hi && (symbolic(a) || symbolic(b)) {
            true => {
                let (l, h, c) = (self.constant(value), self.constant(value), self.constant(value));
                self.intern(Node::Within(l, h, c, id))
            }
            false => self.within(bound.lo, bound.hi, value, id),
        };
        if self.constant_value(coefficient).is_none() {
            self.constructions.insert(coefficient, made);
        }
        self.with_reads_of(coefficient, n)
    }

    /// The expression a rigid-frame proof may use for this coefficient.
    pub fn construction(&self, n: NodeId) -> Node {
        self.constructions.get(&self.numeric(n)).copied().unwrap_or(self.nodes[n as usize])
    }

    /// An operation that read a spatial value, now on constants and params
    /// alone (a rewrite's: `Graph::rebuild`): folded where its float is its
    /// real value (`exact_un`), else kept, on its operand read as a real
    /// number, so its range is rounded outward like any spatial value's.
    pub fn un_real(&mut self, op: UnOp, a: NodeId) -> NodeId {
        let value = self.un_real_value(op, self.numeric(a));
        self.keeping_inputs(value, &[a])
    }

    fn un_real_value(&mut self, op: UnOp, a: NodeId) -> NodeId {
        if let Some(x) = self.as_const(a)
            && let Some(v) = exact_un(op, x)
            && !(self.strict && op == UnOp::Sqrt && x < 0.0)
        {
            return self.folded(v, &[a], &[x], matches!(op, UnOp::Square | UnOp::Exp));
        }
        let a = self.real(a);
        self.un_value(op, a)
    }

    /// `un_real` for two operands (`exact_bin`).
    pub fn bin_real(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        let value = self.bin_real_value(op, self.numeric(a), self.numeric(b));
        self.keeping_inputs(value, &[a, b])
    }

    fn bin_real_value(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        if let (Some(x), Some(y)) = (self.as_const(a), self.as_const(b))
            && let Some(v) = exact_bin(op, x, y)
        {
            return self.folded(v, &[a, b], &[x, y], matches!(op, BinOp::Mul | BinOp::Div | BinOp::Pow | BinOp::Emax));
        }
        if !self.fixed[a as usize] || !self.fixed[b as usize] {
            return self.bin_value(op, a, b);
        }
        // Leave a literal identity (0 times a parameter, a parameter over
        // 1) available to the simplifier. Equal operands remain equal.
        if a == b {
            let a = self.real(a);
            return self.bin_value(op, a, a);
        }
        let (a, b) = if self.as_const(a).is_some() && self.as_const(b).is_none() { (a, self.real(b)) } else { (self.real(a), b) };
        self.bin_value(op, a, b)
    }

    /// The constant `v`, made by arithmetic on a term that read x, y, z or t (`exact`): a
    /// constant to every reader, apart from `v` written.
    pub(super) fn exact_number(&mut self, v: f64) -> NodeId {
        let v = if v == 0.0 { 0.0 } else { v };
        let id = self.intern_keyed(Node::Const(v), Key::Exact(v.to_bits()));
        self.exact.insert(id);
        id
    }

    /// The constant `plain` (the shared node of its float `v`, made as its caller makes one) whose real
    /// value is `e`: `plain` itself where `e` is exactly `v`, else a node of its own, keyed by a hash of
    /// its content (its float, fraction and interval: the same number in two models is one key, whatever
    /// order they were built in). Two numbers whose hashes collide are kept apart: the second takes the
    /// next key.
    pub fn real_constant(&mut self, plain: NodeId, e: crate::exact::Exact) -> NodeId {
        use std::hash::{Hash, Hasher};
        let v = self.constant_value(plain).expect("a constant");
        if (e.band.lo == v && e.band.hi == v) || (e.ratio.is_some() && e.ratio == crate::exact::Ratio::of_f64(v)) {
            return plain;
        }
        let owned = self.owned.contains(&plain);
        let mut h = std::collections::hash_map::DefaultHasher::new();
        (v.to_bits(), e.ratio.map(|r| (r.n, r.d)), e.band.lo.to_bits(), e.band.hi.to_bits()).hash(&mut h);
        let mut key = h.finish();
        loop {
            match self.index.get(&Key::Inexact(key)).copied() {
                Some(n) if self.constant_value(n).map(f64::to_bits) == Some(v.to_bits()) && self.reals.get(&n) == Some(&e) => return n,
                Some(_) => key = key.wrapping_add(1),
                None => {
                    let n = self.intern_keyed(Node::Const(if v == 0.0 { 0.0 } else { v }), Key::Inexact(key));
                    self.reals.insert(n, e);
                    if owned {
                        self.owned.insert(n);
                    }
                    return n;
                }
            }
        }
    }

    /// Constant `n` in real numbers: its own (`real_constant`), else its float. None: no constant.
    pub fn real_value(&self, n: NodeId) -> Option<crate::exact::Exact> {
        let v = self.constant_value(n)?;
        Some(match self.reals.get(&n) {
            Some(e) => *e,
            None => crate::exact::Exact { ratio: crate::exact::Ratio::of_f64(v), band: crate::interval::I::point(v) },
        })
    }

    /// A fold `v` of constants `args` by `op` that rounds (`real_constant`): its real value, by their
    /// fractions where the operation keeps them, else over their intervals (a square root of 2, a sine).
    /// None: the fold is exactly its float, the shared constant.
    pub(super) fn fold_real(&mut self, plain: NodeId, op: Result<UnOp, BinOp>, args: &[NodeId]) -> Option<NodeId> {
        let v = self.constant_value(plain)?;
        use crate::exact::Exact;
        let (a, b) = (self.real_value(args[0])?, args.get(1).and_then(|&b| self.real_value(b)));
        let ratio = match (op, b) {
            (Err(BinOp::Add), Some(b)) => a.ratio.zip(b.ratio).and_then(|(x, y)| x.plus(y)),
            (Err(BinOp::Sub), Some(b)) => a.ratio.zip(b.ratio).and_then(|(x, y)| x.minus(y)),
            (Err(BinOp::Mul), Some(b)) => a.ratio.zip(b.ratio).and_then(|(x, y)| x.times(y)),
            (Err(BinOp::Div), Some(b)) => a.ratio.zip(b.ratio).and_then(|(x, y)| x.over(y)),
            (Ok(UnOp::Neg), None) => a.ratio.map(|x| x.negate()),
            _ => None,
        };
        let e = match (ratio, op, b) {
            (Some(x), ..) => Exact::of(x),
            (None, Ok(u), None) => Exact { ratio: None, band: crate::interval::un(u, a.band) },
            (None, Err(o), Some(b)) => Exact { ratio: None, band: crate::interval::bin(o, a.band, b.band).0 },
            _ => return None,
        };
        if (e.band.lo == v && e.band.hi == v) || (e.ratio.is_some() && e.ratio == crate::exact::Ratio::of_f64(v)) {
            return None;
        }
        Some(self.real_constant(plain, e))
    }

    /// Whether `n` is a constant arithmetic on x, y, z or t made (`exact_number`).
    pub fn is_exact(&self, n: NodeId) -> bool {
        self.exact.contains(&n)
    }

    /// An operation's float `v` on an exact constant where it rounds: within `b` (its operands'
    /// bounds, rounded outward), a number of its own identity, so what reads it is bounded over
    /// all of the band; where that has no finite bound, kept as an operation on a real number.
    pub(super) fn exact_band(&mut self, b: crate::interval::I, v: f64, made: Node) -> Option<NodeId> {
        if !(b.lo.is_finite() && b.hi.is_finite() && b.lo <= v && v <= b.hi) {
            return None;
        }
        use std::hash::{Hash, Hasher};
        let mut h = std::collections::hash_map::DefaultHasher::new();
        (Self::key_of(made), v.to_bits()).hash(&mut h);
        Some(self.within(b.lo, b.hi, v, h.finish() | 1 << 63))
    }
}
