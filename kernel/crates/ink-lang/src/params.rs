//! Params: adding them, the ones a field reads, and their values (exact,
//! at load, over t).

use super::*;

impl Doc {
    /// The params these expressions read.
    /// Field `f` in exact fractions, its params at `params`: where every constant it reads has one, each
    /// param is at a value written in digits (its default or a case), and every operation keeps them
    /// (sums, differences, products, quotients, floors, signs, minima and maxima, choices). None
    /// elsewhere: a param computed through something inexact or set by a solve, a sine, a moment.
    pub fn exact_value(&self, f: NodeId, params: &[f64]) -> Option<ink_core::exact::Ratio> {
        use crate::graph::{BinOp, Node, UnOp};
        use ink_core::exact::Ratio;
        let mut nodes = self.g.reachable(&[f]);
        nodes.sort_unstable();
        let mut at: HashMap<NodeId, Ratio> = HashMap::new();
        for n in nodes {
            let get = |m: NodeId| at.get(&m).copied();
            let r = match self.g.nodes[n as usize] {
                Node::Const(_) => self.g.real_value(n).and_then(|e| e.ratio),
                Node::Var(Var::Param(k)) if !self.solved_params.contains(&(k as usize)) => {
                    let x = params.get(k as usize)?.to_bits();
                    self.param_exact.get(&(k as usize)).and_then(|w| w.iter().find(|(b, _)| *b == x)).map(|w| w.1)
                }
                Node::Tag(_, a) => get(a),
                Node::Un(op, a) => get(a).and_then(|a| match op {
                    UnOp::Neg => Some(a.negate()),
                    UnOp::Abs => Some(if a.n < 0 { a.negate() } else { a }),
                    UnOp::Square => a.times(a),
                    UnOp::Floor => Ratio::new(a.n.div_euclid(a.d), 1),
                    _ => None,
                }),
                Node::Bin(op, a, b) => get(a).zip(get(b)).and_then(|(a, b)| match op {
                    BinOp::Add => a.plus(b),
                    BinOp::Sub => a.minus(b),
                    BinOp::Mul => a.times(b),
                    BinOp::Div if b.n != 0 => a.over(b),
                    BinOp::Min => Some(if a.cmp_to(b)? <= 0 { a } else { b }),
                    BinOp::Max => Some(if a.cmp_to(b)? >= 0 { a } else { b }),
                    _ => None,
                }),
                Node::Select(c, a, b) => get(c).and_then(|c| if c.n < 0 { get(a) } else { get(b) }),
                _ => None,
            };
            match r {
                Some(r) => at.insert(n, r),
                None if n == f => return None,
                None => continue,
            };
        }
        at.get(&f).copied()
    }

    pub fn params_read(&self, nodes: &[NodeId]) -> Vec<usize> {
        self.g
            .reachable(nodes)
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Var(Var::Param(k)) => Some(k as usize),
                _ => None,
            })
            .collect()
    }

    /// Declare a param, named in the graph too (Graph::param_names); its
    /// number.
    /// The params the model's files set (each with its index): not the
    /// allowances its constructions add (`allowances`), which no one writes
    /// or sets, so no listing of its params shows.
    pub fn written_params(&self) -> impl Iterator<Item = (usize, &Param)> {
        self.params.iter().enumerate().filter(|(k, _)| !self.allowances.iter().any(|a| a.0 == *k))
    }

    /// Field `f` moved by a param at 0 that the shape's construction may be
    /// off by up to `width` (`allowances`), named `what` and a number.
    pub(crate) fn allowed(&mut self, f: NodeId, width: f64, what: &str) -> NodeId {
        let name = (0..).map(|i| format!("{what} {i}")).find(|n| !self.params.iter().any(|p| p.name == *n)).expect("a free name");
        let k = self.add_param(Param { name, default: 0.0, tol: None });
        self.allowances.push((k, width));
        let u = self.g.var(Var::Param(k as u32));
        self.g.add(f, u)
    }

    pub fn add_param(&mut self, p: Param) -> usize {
        let k = self.params.len();
        self.g.name_param(k as u32, &p.name);
        self.params.push(p);
        k
    }

    /// Err where moment `t`, or `motion`, is beyond the span a load-time
    /// box or reach of a field reading t holds over (`load_t`).
    pub fn covers_t(&self, t: f64, motion: Option<(f64, f64)>) -> Result<(), InkError> {
        let (lo, hi) = motion.unwrap_or((t, t));
        let Some(((a, b), site)) = self.load_t.iter().find(|((a, b), _)| !(*a <= lo && hi <= *b)) else { return Ok(()) };
        let asked = if lo == hi { format!("t = {t}") } else { format!("the motion {lo} to {hi}") };
        Err(InkError {
            hint: Some("state the motion, (t-span A B), before what is bounded at load, covering every moment asked".into()),
            ..InkError::msg(format!("{asked} is beyond t in [{a}, {b}], the moments a box found at load at {site} holds over"))
        })
    }

    /// The ranged params a load-time box or reach was found over the band
    /// of, with it and where (`covered`).
    pub fn covered(&self, k: usize) -> Option<((f64, f64), &str)> {
        self.covered.iter().find(|c| c.0 == k).map(|c| (c.1, c.2.as_str()))
    }

    /// Where param k's value first became a number at load time (a loop
    /// count, a sampled or measured field, a solve's input...): file:line
    /// (head name). None: it stays symbolic, so it can be ranged.
    pub fn materialized(&self, k: usize) -> Option<String> {
        self.materialized.borrow().get(&k).cloned()
    }

    /// What Rust builds the model from is noted where outside ink's range
    /// (Graph::note): an axis of 1.8e-162 normalized squares to 0.
    pub(super) fn load_value(&self, f: NodeId) -> Option<f64> {
        if motion::needs_speed(self, f) {
            return None;
        }
        let v = match self.g.as_const(f) {
            Some(v) => v,
            None if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) => return None,
            None => {
                self.bake(&[f]);
                self.load_eval(f)
            }
        };
        self.g.note(v, crate::graph::wild(v));
        Some(v)
    }
}
