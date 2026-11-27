//! Continuity, zero sets, and gradient changes over a region.

use super::{BOTH, FIRST, Ins, Op, SECOND, Tape};
use crate::graph::{BinOp, UnOp, Var};

/// What a region holds besides what a mesh draws there (Tape::sheet): one
/// sheet, nothing else; what is thinner than this (along the direction
/// checked); or unknown.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Sheet {
    One,
    Thin(f64),
    Unknown,
}

impl Tape {
    /// Whether select `i` is continuous where it switches, by its form, in
    /// the region `choices` were found for: its condition is x - y (or
    /// y - x), x and y continuous there (`continuous`: else the condition
    /// may jump past 0, never meeting x = y), and its branches one
    /// expression but for x and y swapped (g(x) and g(y): select(x - y, x,
    /// y) is a min), so where the condition reaches 0, x = y, they are one
    /// value.
    pub fn continuous_select(&self, i: usize, choices: &[u8]) -> bool {
        let ssa = &self.ssa;
        let strip = |j: usize| self.through(choices, j);
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
        self.seams(choices, &self.live(choices)).contains(&true)
    }

    /// Per instruction, whether it is a select the output reads (`live`)
    /// that may switch between branches that are not one field (`switches`).
    fn seams(&self, choices: &[u8], live: &[bool]) -> Vec<bool> {
        let strip = |j: usize| self.through(choices, j);
        (0..self.ssa.len())
            .map(|i| {
                let ins = &self.ssa[i];
                live[i] && matches!((ins.op, choices[i]), (Op::Select, BOTH)) && strip(ins.b as usize) != strip(ins.c as usize) && !self.continuous_select(i, choices)
            })
            .collect()
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
        self.sheets(choices, grads).way(centre, d)
    }

    /// What `sheet_way` reads of a region whatever the direction, from
    /// `choices` and `grads` found over it: found once, walked along any
    /// number of directions (`Sheets::way`).
    pub fn sheets<'a>(&'a self, choices: &'a [u8], grads: &'a [crate::interval::G]) -> Sheets<'a> {
        let live = self.live(choices);
        Sheets { tape: self, choices, grads, seams: self.seams(choices, &live), thin: self.thin_zeros(grads), live }
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
        self.turns(choices, grads, &crate::interval::apart)
    }

    /// Whether `may_crease` may hold for some gradients: a min, max, emax
    /// or select it reads is undecided, an abs's value may be 0 or a
    /// floor's step, from `grads`' values alone. Where none is, the
    /// output's gradient turns smoothly whatever its gradients.
    pub fn may_turn(&self, choices: &[u8], grads: &[crate::interval::G]) -> bool {
        self.turns(choices, grads, &|_, _| true)
    }

    /// `may_crease` with `apart` telling where two branches' gradients
    /// point apart.
    fn turns(&self, choices: &[u8], grads: &[crate::interval::G], apart: &dyn Fn(&[crate::interval::I; 3], &[crate::interval::I; 3]) -> bool) -> bool {
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
}

/// A region's account for `Tape::sheet_way` that does not depend on the
/// direction (`Tape::sheets`): the values the output reads through the
/// branches `choices` pick (`live`), its selects that may switch between
/// branches that are not one field (`seams`), and the values whose 0s
/// have no volume there (`thin`).
pub struct Sheets<'a> {
    tape: &'a Tape,
    choices: &'a [u8],
    grads: &'a [crate::interval::G],
    live: Vec<bool>,
    seams: Vec<bool>,
    thin: Vec<bool>,
}

impl Sheets<'_> {
    /// Whether a select may switch there between branches that are not
    /// one field (`Tape::switches`).
    pub fn switches(&self) -> bool {
        self.seams.contains(&true)
    }

    /// The instructions a walk along a direction runs: those the output
    /// reads.
    pub fn reads(&self) -> usize {
        self.live.iter().filter(|&&l| l).count()
    }

    /// `Tape::sheet`.
    pub fn sheet(&self, centre: Option<(&[crate::interval::G], crate::V3)>, d: crate::V3) -> Sheet {
        self.way(centre, d).0
    }

    /// `Tape::sheet_way`: each value the output reads, in order (the
    /// others are never read).
    pub fn way(&self, centre: Option<(&[crate::interval::G], crate::V3)>, d: crate::V3) -> (Sheet, i8) {
        use crate::interval::{I, add_o, div_o, hull, mul_o, sub_o};
        let (tape, choices, grads, thin) = (self.tape, self.choices, self.grads, &self.thin);
        let n = tape.ssa.len();
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
        let strip = |j: usize| tape.through(choices, j);
        let below = |j: u32| v(j).all_inside() || v(j).hi <= 0.0 && thin[j as usize];
        let not_below = |j: u32| v(j).none_inside();
        for i in 0..n {
            if !self.live[i] {
                r.push(I::point(0.0));
                way.push(0);
                thick.push(0.0);
                continue;
            }
            let ins = &tape.ssa[i];
            let at = |j: u32| r[j as usize];
            let by_gradient = || {
                let g = grads[i].d;
                add_o(add_o(mul_o(g[0], I::point(d[0])), mul_o(g[1], I::point(d[1]))), mul_o(g[2], I::point(d[2])))
            };
            let mut used_gradient = false;
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
                    let still = at(ins.a) == I::point(0.0);
                    match still || !self.seams[i] {
                        true => hull(at(ins.b), at(ins.c)),
                        false => I::new(f64::NEG_INFINITY, f64::INFINITY),
                    }
                }
                // This is already the directional gradient bound.
                _ => {
                    used_gradient = true;
                    by_gradient()
                }
            };
            // Intersect a chain-rule bound with the directional gradient.
            // A singleton (in particular an enclosed constant's zero rate)
            // cannot be tightened. Keep the chain rule if the bounds miss.
            let rate = if used_gradient || rate.lo == rate.hi {
                rate
            } else {
                let grad = by_gradient();
                let (lo, hi) = (rate.lo.max(grad.lo), rate.hi.min(grad.hi));
                if lo <= hi { I::new(lo, hi) } else { rate }
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
                (Op::Bin(BinOp::Sub), _) if matches!(tape.ssa[strip(ins.a as usize)].op, Op::Un(UnOp::Abs)) => {
                    let (x, y) = (tape.ssa[strip(ins.a as usize)].a, ins.b);
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
}
