use super::{BinOp, Graph, Node, NodeId, UnOp, Var};
use rustc_hash::FxHashMap as HashMap;

impl Graph {
    /// Derivative of min/max: the branch taken; at a tie, the mean of both,
    /// da/2 + db/2 (as the evaluators do): halving each first, since the sum
    /// of constant sides can overflow where their mean does not (the mean of
    /// 1e308 and 1.5e308). One constant slope on both sides is that slope
    /// where its halves are exact; below 2^-1021 they round (the mean of
    /// 5e-324 and 5e-324 is 0), so it is the mean there too.
    fn d_branch(&mut self, op: BinOp, a: NodeId, b: NodeId, da: NodeId, db: NodeId) -> NodeId {
        if da == db && self.as_const(da).is_some_and(|v| 0.5 * v + 0.5 * v == v) {
            return da;
        }
        let (lo, hi) = if op == BinOp::Min { (da, db) } else { (db, da) };
        let c = self.bin_real(BinOp::Sub, a, b);
        let nc = self.neg(c);
        let half = self.constant(0.5);
        let (ha, hb) = (self.d_mul(half, da), self.d_mul(half, db));
        let mean = self.bin_real(BinOp::Add, ha, hb);
        let tie_or_hi = self.select(nc, hi, mean);
        self.select(c, lo, tie_or_hi)
    }

    /// A derivative's term `l * r`, one of them a derivative: where that is the constant 0 (what does
    /// not vary along the axis), the term is that 0, not a dropped value (`zero_of`): a slope of 0 is
    /// exact, and reads nothing more than the derivative it came from.
    fn d_mul(&mut self, l: NodeId, r: NodeId) -> NodeId {
        let zero = |g: &Self, n: NodeId, o: NodeId| g.as_const(n) == Some(0.0) && !g.hides(o);
        if zero(self, l, r) {
            return l;
        }
        if zero(self, r, l) {
            return r;
        }
        self.bin_real(BinOp::Mul, l, r)
    }

    /// Symbolic derivative of `root` with respect to a variable, as a new
    /// expression. min/max/abs/select differentiate the branch taken; a
    /// sampled field's read, by the chain rule through its coordinates,
    /// with the reads one axis further along its path (Node::Read: about
    /// its source, unknown). Each op's derivative is written as the
    /// evaluators compute it (eval::eval_grads, eval_dual): its partials,
    /// from its operands' values (dbin, dun; a read's: the reads one axis
    /// further), times its operands' derivatives, summed in order. So the
    /// derivative's value is theirs, bit for bit, at every scale.
    pub fn diff(&mut self, root: NodeId, var: Var) -> NodeId {
        let ids = self.reach(&[root], false, false);
        let mut d: HashMap<NodeId, NodeId> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        let zero = self.constant(0.0);
        let one = self.constant(1.0);
        for &id in &ids {
            let node = self.nodes[id as usize];
            let dn = match node {
                Node::Const(_) => zero,
                Node::Var(v) => {
                    if v == var {
                        one
                    } else {
                        zero
                    }
                }
                Node::Un(op, a) => {
                    let da = d[&a];
                    if da == zero {
                        zero
                    } else {
                        match op {
                            UnOp::Neg => self.neg(da),
                            UnOp::Abs => {
                                // At 0 the average of both sides: 0.
                                let n = self.neg(da);
                                let na = self.neg(a);
                                let pos = self.select(na, da, zero);
                                self.select(a, n, pos)
                            }
                            UnOp::Sqrt => {
                                // (0.5 / sqrt a) da where a > 0, else 0.
                                let half = self.constant(0.5);
                                let s = self.un_real(UnOp::Sqrt, a);
                                let f = self.bin_real(BinOp::Div, half, s);
                                let q = self.d_mul(f, da);
                                let na = self.neg(a);
                                self.select(na, q, zero)
                            }
                            UnOp::Square => {
                                let two = self.constant(2.0);
                                let t = self.bin_real(BinOp::Mul, two, a);
                                self.d_mul(t, da)
                            }
                            UnOp::Exp => self.d_mul(id, da),
                            UnOp::Ln => {
                                // ln of a non-positive is a constant (-inf, or 0).
                                let f = self.bin_real(BinOp::Div, one, a);
                                let q = self.d_mul(f, da);
                                let na = self.neg(a);
                                self.select(na, q, zero)
                            }
                            UnOp::Sin => {
                                let c = self.un_real(UnOp::Cos, a);
                                self.d_mul(c, da)
                            }
                            UnOp::Cos => {
                                let s = self.un_real(UnOp::Sin, a);
                                let n = self.neg(s);
                                self.d_mul(n, da)
                            }
                            UnOp::Tan => {
                                let t2 = self.un_real(UnOp::Square, id);
                                let f = self.bin_real(BinOp::Add, one, t2);
                                self.d_mul(f, da)
                            }
                            UnOp::Floor => zero,
                            // Linear: the same exact scaling of the
                            // derivative, read as a real number where it is
                            // constants and params (never folded into a
                            // rounded coefficient).
                            UnOp::Rad | UnOp::Deg => {
                                let r = self.real(da);
                                self.un_real(op, r)
                            }
                        }
                    }
                }
                Node::Bin(op, a, b) => {
                    let (da, db) = (d[&a], d[&b]);
                    match op {
                        BinOp::Add => self.bin_real(BinOp::Add, da, db),
                        BinOp::Sub => self.bin_real(BinOp::Sub, da, db),
                        BinOp::Mul => {
                            let x = self.d_mul(da, b);
                            let y = self.d_mul(a, db);
                            self.bin_real(BinOp::Add, x, y)
                        }
                        BinOp::Div => {
                            // (1 / b) da - ((a / b) / b) db (b^2 underflows to
                            // 0 for tiny b; x / 0 = 0 at b = 0).
                            let fa = self.bin_real(BinOp::Div, one, b);
                            let x = self.d_mul(fa, da);
                            let q = self.bin_real(BinOp::Div, a, b);
                            let qb = self.bin_real(BinOp::Div, q, b);
                            let y = self.d_mul(qb, db);
                            let r = self.bin_real(BinOp::Sub, x, y);
                            let nb = self.un_real(UnOp::Abs, b);
                            let nb = self.neg(nb);
                            self.select(nb, r, zero)
                        }
                        BinOp::Min | BinOp::Max => self.d_branch(op, a, b, da, db),
                        BinOp::Emax => {
                            // Both positive: (a / emax) da + (b / emax) db (dividing
                            // first, as the evaluators do: a da can overflow), each
                            // weight at most 1, as it is in real numbers (so bounds
                            // of a slope through it see that); else as max.
                            let e = self.bin_real(BinOp::Emax, a, b);
                            let (ae, be) = (self.bin_real(BinOp::Div, a, e), self.bin_real(BinOp::Div, b, e));
                            let (ae, be) = (self.min(ae, one), self.min(be, one));
                            let x = self.d_mul(ae, da);
                            let y = self.d_mul(be, db);
                            let round = self.bin_real(BinOp::Add, x, y);
                            let sharp = self.d_branch(BinOp::Max, a, b, da, db);
                            let both = self.min(a, b);
                            let nb = self.neg(both);
                            self.select(nb, round, sharp)
                        }
                        BinOp::Atan2 => {
                            // atan2(a, b): (b / r2) da + (-a / r2) db, with
                            // r2 = a^2 + b^2 > 0.
                            let a2 = self.un_real(UnOp::Square, a);
                            let b2 = self.un_real(UnOp::Square, b);
                            let den = self.bin_real(BinOp::Add, a2, b2);
                            let fa = self.bin_real(BinOp::Div, b, den);
                            let na = self.neg(a);
                            let fb = self.bin_real(BinOp::Div, na, den);
                            let x = self.d_mul(fa, da);
                            let y = self.d_mul(fb, db);
                            let q = self.bin_real(BinOp::Add, x, y);
                            let nd = self.neg(den);
                            self.select(nd, q, zero)
                        }
                        BinOp::Pow => {
                            // Where a < 0 and b is not a whole number the value is
                            // the constant 0: so is its derivative.
                            let fl = self.un_real(UnOp::Floor, b);
                            let fr = self.bin_real(BinOp::Sub, b, fl);
                            let fr = self.un_real(UnOp::Abs, fr);
                            let nfr = self.neg(fr);
                            let zero_here = self.max(a, nfr);
                            if db == zero {
                                // b a^(b-1) da
                                let bm1 = self.bin_real(BinOp::Sub, b, one);
                                let p = self.bin_real(BinOp::Pow, a, bm1);
                                let t = self.bin_real(BinOp::Mul, b, p);
                                let r = self.d_mul(t, da);
                                self.select(zero_here, zero, r)
                            } else {
                                // b a^(b-1) da + a^b ln(a) db, the second term
                                // only where a > 0 (elsewhere a^b does not vary
                                // smoothly with b), as the evaluators.
                                let bm1 = self.bin_real(BinOp::Sub, b, one);
                                let p = self.bin_real(BinOp::Pow, a, bm1);
                                let fa = self.bin_real(BinOp::Mul, b, p);
                                let x = self.d_mul(fa, da);
                                let la = self.un_real(UnOp::Ln, a);
                                let vl = self.bin_real(BinOp::Mul, id, la);
                                let na = self.neg(a);
                                let fb = self.select(na, vl, zero);
                                let y = self.d_mul(fb, db);
                                let r = self.bin_real(BinOp::Add, x, y);
                                self.select(zero_here, zero, r)
                            }
                        }
                    }
                }
                Node::Select(c, a, b) => {
                    let (da, db) = (d[&a], d[&b]);
                    self.select(c, da, db)
                }
                Node::Tag(_, a) | Node::Lipschitz(a, _, _) | Node::Deriv(a, _, _) | Node::Real(a) => d[&a],
                // Frame, Within, and CornerReach carry no symbolic derivative:
                // CornerReach is a condition-only corner predicate whose sign
                // drives select in fillet-union / chamfer-union.
                Node::Frame(_) | Node::Within(..) | Node::CornerReach(..) => zero,
                // The interpolant: chain rule through the coordinates.
                Node::Grid(g, m, c) => {
                    let mut acc = zero;
                    for k in 0..3 {
                        let dc = d[&c[k]];
                        if dc == zero {
                            continue;
                        }
                        let Some(dm) = self.grids[g as usize].derivative(m, k) else { continue };
                        let slope = self.intern(Node::Grid(g, dm, c));
                        let term = self.bin_real(BinOp::Mul, slope, dc);
                        acc = self.bin_real(BinOp::Add, acc, term);
                    }
                    acc
                }
                // A read: chain rule through its coordinates, with the read
                // one axis further along its path. The samples' derivative
                // can simplify where the source's does not. Where the
                // coordinates do not move with `var`, neither does the source.
                Node::Read(g, path, _, at) => {
                    let Node::Frame(c) = self.nodes[at as usize] else { unreachable!("a read's coordinates are a frame") };
                    let mut acc = zero;
                    for k in 0..3 {
                        let dc = d[&c[k]];
                        if dc == zero {
                            continue;
                        }
                        let mut along = self.paths[path as usize].clone();
                        along.push(k as u8);
                        let along = self.path(along);
                        let slope = self.read(g, along, c);
                        let term = self.bin_real(BinOp::Mul, slope, dc);
                        acc = self.bin_real(BinOp::Add, acc, term);
                    }
                    acc
                }
            };
            d.insert(id, dn);
        }
        self.derivative_domain(d[&root], root, var)
    }

    /// `diff`, and where `var` is x, y or z, marked as the derivative of
    /// `root` (Node::Deriv), so ranges also hold it within root's gradient
    /// bounds; where it is t and root reads t, marked so even where it is
    /// constant, so a claim on a speed checks the motion for jumps. What
    /// the language's `diff` and `grad` make.
    pub fn derivative(&mut self, root: NodeId, var: Var) -> NodeId {
        let a = self.diff(root, var);
        let axis = match var {
            Var::X => 0,
            Var::Y => 1,
            Var::Z => 2,
            // Through marks: a speed whose formula is constant still marks
            // its own drive, so an acceleration checks both.
            Var::T if self.depends_on(root, |v| v == Var::T) => 3,
            _ => return a,
        };
        if axis < 3 && self.as_const(a).is_some() {
            return a;
        }
        self.intern(Node::Deriv(a, root, axis))
    }

    /// `derivative` of `root` along t, marked as the speed of each of
    /// `read` instead of root's: what reads root (a turn's sine and cosine
    /// of its angle), each continuous where root jumps by what none of them
    /// sees (a whole turn), so its speed is the formula's there.
    pub fn speed_read(&mut self, root: NodeId, read: &[NodeId]) -> NodeId {
        let a = self.diff(root, Var::T);
        if !self.depends_on(root, |v| v == Var::T) {
            // The rate folds to zero, but evaluating the supplied reads
            // can still require another drive to be defined.
            return read.iter().fold(a, |n, &f| self.with_reads_of(n, f));
        }
        self.mark_speed(a, read)
    }

    /// Keep the drives whose continuity a time derivative requires, even
    /// when its formula simplifies to a constant.
    pub fn mark_speed(&mut self, rate: NodeId, drives: &[NodeId]) -> NodeId {
        drives.iter().fold(rate, |m, &f| self.intern(Node::Deriv(m, f, 3)))
    }
}
