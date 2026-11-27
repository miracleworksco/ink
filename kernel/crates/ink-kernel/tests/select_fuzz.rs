//! Interval bounds over a box and a step of t, where selects on t fork
//! (and selects on one condition fork together), against the field
//! sampled at points and moments: random fields of selects sharing and
//! nearly sharing conditions, nested in each other, through placements
//! turning with t, time remapped, min, max, abs and floor (of values
//! exact in floats: bounds hold the field in real numbers, and a floor of
//! a rounded value steps where the real one does not), along t or a
//! ranged param. Every sampled
//! value must lie in the bound over the step, within the first-order
//! bound from the step's centre (rates and jump), and, where nothing may
//! jump, move no faster than the rate. A field `still` must not move.
//!
//! 1,200 fields here (seconds); INK_FUZZ_CASES=n runs n, from
//! INK_FUZZ_SEED, and INK_FUZZ_OUT=file writes each bound's width (to
//! compare two builds' cost).

use ink_core::eval::{Ctx, Inputs, eval_interval_rates_jump, eval_interval_scope, eval_interval_scope_every, eval_points, span_mid};
use ink_core::graph::{Graph, Node, NodeId, UnOp, Var};
use ink_core::interval::{I, mag};
use ink_core::tape::Tape;
use std::fmt::Write as _;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
    fn pick<T: Copy>(&mut self, xs: &[T]) -> T {
        xs[self.below(xs.len())]
    }
}

/// A random field: its graph, root, and the moments its conditions switch.
struct Gen {
    g: Graph,
    rng: Rng,
    conds: Vec<NodeId>,
    pool: Vec<NodeId>,
    seams: Vec<f64>,
    /// The input the motion is along: t, or a ranged param.
    tv: Var,
}

impl Gen {
    fn c(&mut self, v: f64) -> NodeId {
        self.g.constant(v)
    }
    fn k(&mut self, xs: &[f64]) -> NodeId {
        let v = self.rng.pick(xs);
        self.c(v)
    }
    fn t(&mut self) -> NodeId {
        self.g.var(self.tv)
    }
    /// Conditions switching at `s`, written several ways (the same node
    /// reused most often), and one switching a hair later.
    fn conditions(&mut self, s: f64) {
        let t = self.t();
        let cs = self.c(s);
        let a = self.g.sub(t, cs);
        self.conds.extend([a, a, a]);
        let b = self.g.sub(cs, t);
        self.conds.push(b);
        let two = self.c(2.0);
        let t2 = self.g.mul(two, t);
        let s2 = self.c(2.0 * s);
        let c2 = self.g.sub(t2, s2);
        self.conds.push(c2);
        let three = self.c(3.0);
        let c3 = self.g.mul(a, three);
        self.conds.push(c3);
        let near = self.c(s + 1.0 / 1024.0);
        let c4 = self.g.sub(t, near);
        self.conds.push(c4);
        // A condition through a floor of t: not a Moment, bounded across.
        if self.rng.chance(0.3) {
            let four = self.c(4.0);
            let ft = self.g.mul(four, t);
            let f = self.g.un(UnOp::Floor, ft);
            let k = self.c((4.0 * s).floor() + 0.5);
            let c5 = self.g.sub(f, k);
            self.conds.push(c5);
        }
    }
    fn leaf(&mut self) -> NodeId {
        match self.rng.below(8) {
            0 => self.g.var(Var::X),
            1 => self.g.var(Var::Y),
            2 => self.g.var(Var::Z),
            3 | 4 => self.t(),
            5 => {
                let v = (self.rng.below(13) as f64 - 6.0) / 2.0;
                self.c(v)
            }
            _ if !self.pool.is_empty() => {
                let k = self.rng.below(self.pool.len());
                self.pool[k]
            }
            _ => self.t(),
        }
    }
    fn expr(&mut self, depth: usize) -> NodeId {
        if depth == 0 || self.rng.chance(0.2) {
            return self.leaf();
        }
        let n = match self.rng.below(16) {
            0..=4 => {
                let c = self.rng.pick(&self.conds.clone());
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.select(c, a, b)
            }
            5 => {
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.add(a, b)
            }
            6 => {
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.sub(a, b)
            }
            7 => {
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.mul(a, b)
            }
            8 => {
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.min(a, b)
            }
            9 => {
                let (a, b) = (self.expr(depth - 1), self.expr(depth - 1));
                self.g.max(a, b)
            }
            10 => {
                let a = self.expr(depth - 1);
                let op = self.rng.pick(&[UnOp::Neg, UnOp::Abs, UnOp::Sin, UnOp::Cos, UnOp::Square]);
                self.g.un(op, a)
            }
            11 => {
                // A turn at a rate, or a time remap.
                let t = self.t();
                let k = self.k(&[1.0, 2.0, 3.0, 7.5, -4.0]);
                let a = self.g.mul(k, t);
                if self.rng.chance(0.5) { self.g.un(UnOp::Sin, a) } else { self.g.un(UnOp::Cos, a) }
            }
            12 => {
                // A floor of what is exact in floats (bounds hold the
                // field in real numbers; a floor of a rounded value can
                // step where the real one does not).
                let a = match self.rng.below(4) {
                    0 => self.g.var(Var::X),
                    1 => self.g.var(Var::Z),
                    _ => self.t(),
                };
                let k = self.k(&[1.0, 2.0, 4.0]);
                let m = self.g.mul(k, a);
                self.g.un(UnOp::Floor, m)
            }
            13 => {
                // A select inside a time remap: t := a t + b.
                let e = self.expr(depth - 1);
                let t = self.t();
                let a = self.k(&[1.0, 2.0, 0.5, -1.0]);
                let b = self.k(&[0.0, 0.25, 0.1, 1.0]);
                let at = self.g.mul(a, t);
                let r = self.g.add(at, b);
                self.g.subst(e, &[(self.tv, r)])
            }
            14 => {
                // Moved along x by a select on t: a placement's entry.
                let e = self.expr(depth - 1);
                let c = self.rng.pick(&self.conds.clone());
                let (p, q) = (self.expr(depth - 1), self.expr(depth - 1));
                let off = self.g.select(c, p, q);
                let x = self.g.var(Var::X);
                let xm = self.g.sub(x, off);
                self.g.subst(e, &[(Var::X, xm)])
            }
            _ => {
                let a = self.expr(depth - 1);
                let s = self.g.un(UnOp::Square, a);
                let one = self.c(1.0);
                let s1 = self.g.add(s, one);
                self.g.un(UnOp::Sqrt, s1)
            }
        };
        self.pool.push(n);
        n
    }
    /// A part handed between placements by selects on two conditions,
    /// turning in the middle one: each entry of its placement a select.
    fn handed(&mut self) -> NodeId {
        let (c1, c2) = (self.rng.pick(&self.conds.clone()), self.rng.pick(&self.conds.clone()));
        let turns = 1 + self.rng.below(3);
        let mut r: [[NodeId; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| self.g.constant(if i == j { 1.0 } else { 0.0 })));
        for k in 0..turns {
            let t = self.t();
            let w = self.k(&[11.0, 12.0, 3.0, -5.0]);
            let a = self.g.mul(w, t);
            let (cs, sn) = (self.g.un(UnOp::Cos, a), self.g.un(UnOp::Sin, a));
            let nsn = self.g.neg(sn);
            let (zero, one) = (self.c(0.0), self.c(1.0));
            let m = if k % 2 == 0 { [[one, zero, zero], [zero, cs, nsn], [zero, sn, cs]] } else { [[cs, nsn, zero], [sn, cs, zero], [zero, zero, one]] };
            r = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    let mut s = self.c(0.0);
                    for (l, row) in m.iter().enumerate() {
                        let p = self.g.mul(r[i][l], row[j]);
                        s = self.g.add(s, p);
                    }
                    s
                })
            });
        }
        let id: [[NodeId; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| self.g.constant(if i == j { 1.0 } else { 0.0 })));
        let tr = [[10.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0]].map(|v| v.map(|x| self.g.constant(x)));
        let rot: [[NodeId; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let inner = self.g.select(c2, r[i][j], id[i][j]);
                self.g.select(c1, id[i][j], inner)
            })
        });
        let tt: [NodeId; 3] = std::array::from_fn(|i| {
            let inner = self.g.select(c2, tr[1][i], tr[2][i]);
            self.g.select(c1, tr[0][i], inner)
        });
        let p = [Var::X, Var::Y, Var::Z].map(|v| self.g.var(v));
        let d: [NodeId; 3] = std::array::from_fn(|i| self.g.sub(p[i], tt[i]));
        // q = R^T (p - T)
        let q: [NodeId; 3] = std::array::from_fn(|i| {
            let mut s = self.c(0.0);
            for (j, dj) in d.iter().enumerate() {
                let m = self.g.mul(rot[j][i], *dj);
                s = self.g.add(s, m);
            }
            s
        });
        let ball = |gn: &mut Gen, c: [f64; 3], rad: f64| {
            let mut s = gn.c(0.0);
            for (i, qi) in q.iter().enumerate() {
                let ci = gn.c(c[i]);
                let u = gn.g.sub(*qi, ci);
                let sq = gn.g.un(UnOp::Square, u);
                s = gn.g.add(s, sq);
            }
            let r = gn.g.un(UnOp::Sqrt, s);
            let rad = gn.c(rad);
            gn.g.sub(r, rad)
        };
        let (a, b) = (ball(self, [0.0; 3], 2.0), ball(self, [3.0, 0.0, 0.0], 1.0));
        let part = self.g.min(a, b);
        if self.rng.chance(0.5) {
            let x = self.g.var(Var::X);
            let k = self.k(&[-20.0, 5.0, 9.0]);
            let wall = self.g.sub(k, x);
            let wall = self.g.neg(wall);
            self.g.min(part, wall)
        } else {
            part
        }
    }
}

fn generate(seed: u64, tv: Var) -> (Graph, NodeId, Vec<f64>) {
    let mut gn = Gen { g: Graph::new(), rng: Rng(seed), conds: Vec::new(), pool: Vec::new(), seams: Vec::new(), tv };
    let n = 1 + gn.rng.below(3);
    for _ in 0..n {
        let s = gn.rng.pick(&[0.25, 0.3, 0.5, 0.6, 0.75, 0.375, 0.125]);
        gn.seams.push(s);
        gn.conditions(s);
    }
    let root = match gn.rng.below(4) {
        0 => gn.handed(),
        1 => {
            let h = gn.handed();
            let e = gn.expr(3);
            let m = gn.g.min(h, e);
            if gn.rng.chance(0.5) {
                let t = gn.t();
                let two = gn.c(2.0);
                let tt = gn.g.mul(two, t);
                gn.g.subst(m, &[(gn.tv, tt)])
            } else {
                m
            }
        }
        _ => {
            let d = 3 + gn.rng.below(4);
            let a = gn.expr(d);
            let b = gn.expr(d);
            match gn.rng.below(3) {
                0 => gn.g.min(a, b),
                1 => gn.g.add(a, b),
                _ => gn.g.max(a, b),
            }
        }
    };
    // Speeds read through the derivative along t.
    let root = if tv == Var::T && gn.rng.chance(0.15) { gn.g.derivative(root, Var::T) } else { root };
    let seams = gn.seams.clone();
    (gn.g, root, seams)
}

fn show(g: &Graph, n: NodeId, depth: usize, out: &mut String) {
    if depth > 40 {
        out.push_str("...");
        return;
    }
    match g.nodes[n as usize] {
        Node::Const(v) => _ = write!(out, "{v}"),
        Node::Var(v) => _ = write!(out, "{v:?}"),
        Node::Un(op, a) => {
            _ = write!(out, "({op:?} ");
            show(g, a, depth + 1, out);
            out.push(')');
        }
        Node::Bin(op, a, b) => {
            _ = write!(out, "({op:?} ");
            show(g, a, depth + 1, out);
            out.push(' ');
            show(g, b, depth + 1, out);
            out.push(')');
        }
        Node::Select(c, a, b) => {
            out.push_str("(select ");
            show(g, c, depth + 1, out);
            out.push(' ');
            show(g, a, depth + 1, out);
            out.push(' ');
            show(g, b, depth + 1, out);
            out.push(')');
        }
        Node::Deriv(a, f, axis) => {
            out.push_str("(deriv ");
            show(g, a, depth + 1, out);
            out.push(' ');
            show(g, f, depth + 1, out);
            _ = write!(out, " {axis})");
        }
        other => _ = write!(out, "{other:?}"),
    }
}

/// How far a float evaluation may stray from the real number a bound holds.
fn slack(x: f64) -> f64 {
    1e-9 * (1.0 + x.abs())
}

#[derive(Default)]
struct Tally {
    cases: usize,
    checks: usize,
    forked: usize,
    failures: Vec<String>,
    widths: String,
    still: usize,
    still_not_rate0: usize,
    rate0_not_still: usize,
    skipped: usize,
}

fn spans(rng: &mut Rng, seams: &[f64]) -> Vec<(f64, f64)> {
    let mut out = vec![(0.0, 1.0)];
    for _ in 0..5 {
        let s = rng.pick(seams);
        let w = rng.pick(&[0.3, 0.05, 1e-3, 1e-7]);
        out.push(match rng.below(4) {
            0 => (s - w, s + w),
            1 => (s, s + w),
            2 => (s - w, s),
            _ => {
                let a = rng.unit();
                (a, (a + w).min(1.0).max(a + 1e-9))
            }
        });
    }
    out
}

fn boxes(rng: &mut Rng) -> Vec<[I; 3]> {
    let mut out = Vec::new();
    for _ in 0..4 {
        let c: [f64; 3] = std::array::from_fn(|_| rng.pick(&[0.0, 10.0, 9.0, -1.0, 3.0]) + (rng.unit() - 0.5) * rng.pick(&[0.0, 2.0, 8.0]));
        let h = rng.pick(&[0.0, 0.01, 0.3, 2.0]);
        out.push(c.map(|x| I::new(x - h, x + h)));
    }
    out
}

fn points(rng: &mut Rng, b: [I; 3]) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for k in 0..8 {
        out.push(std::array::from_fn(|d| if k >> d & 1 == 0 { b[d].lo } else { b[d].hi }));
    }
    out.push(b.map(|x| x.mid()));
    for _ in 0..12 {
        out.push(std::array::from_fn(|d| b[d].lo + (b[d].hi - b[d].lo) * rng.unit()));
    }
    out.dedup();
    out
}

fn moments(rng: &mut Rng, s: (f64, f64), seams: &[f64]) -> Vec<f64> {
    let mut out = vec![s.0, s.1, span_mid(s)];
    for &x in seams {
        for y in [x, x.next_down(), x.next_up(), x + 1.0 / 1024.0] {
            if s.0 <= y && y <= s.1 {
                out.push(y);
            }
        }
    }
    for k in 0..24 {
        out.push(s.0 + (s.1 - s.0) * (k as f64 + rng.unit()) / 24.0);
    }
    out.retain(|&t| s.0 <= t && t <= s.1);
    out.sort_by(f64::total_cmp);
    out.dedup();
    out
}

/// The context at moment `t` along `tv`.
fn at(tv: Var, t: f64) -> Ctx {
    if tv == Var::T { Ctx::new(t, vec![]) } else { Ctx::new(0.0, vec![t]) }
}

/// Values of `tape` at `pts` at each of `ts` along `tv`: [t][point].
fn sample(tv: Var, tape: &Tape, ts: &[f64], pts: &[[f64; 3]]) -> Vec<Vec<f64>> {
    let mut out = Vec::new();
    ts.iter()
        .map(|&t| {
            eval_points(tape, &at(tv, t), pts, &mut out);
            out.clone()
        })
        .collect()
}

/// The checks of one value's bound `v`, centre bound `m`, rate along t
/// `d` and jump `j` over box and step against its samples.
fn check(what: &str, v: I, m: I, d: I, j: f64, s: (f64, f64), ts: &[f64], vals: &[Vec<f64>], fail: &mut Vec<String>) {
    let tc = span_mid(s);
    let r = mag(d);
    for (ti, &t) in ts.iter().enumerate() {
        for (pi, &x) in vals[ti].iter().enumerate() {
            if !x.is_finite() {
                continue;
            }
            if !(v.lo - slack(x) <= x && x <= v.hi + slack(x)) {
                fail.push(format!("{what}: value {x} at t {t} point#{pi} outside {v:?}"));
                return;
            }
            let e = r * (t - tc).abs() + j;
            if e.is_finite() && !(m.lo - e - slack(x) <= x && x <= m.hi + e + slack(x)) {
                fail.push(format!("{what}: value {x} at t {t} point#{pi} outside centre {m:?} +- {e} (rate {d:?}, jump {j})"));
                return;
            }
        }
    }
    // Nothing may jump: no faster than the rate between any two moments.
    if j == 0.0 && r.is_finite() {
        for ti in 1..ts.len() {
            for pi in 0..vals[ti].len() {
                let (a, b) = (vals[ti - 1][pi], vals[ti][pi]);
                let dt = ts[ti] - ts[ti - 1];
                if a.is_finite() && b.is_finite() && (b - a).abs() > r * dt * (1.0 + 1e-9) + slack(a) + slack(b) {
                    fail.push(format!("{what}: moves {} in {dt} at t {} point#{pi}: faster than rate {d:?}", b - a, ts[ti]));
                    return;
                }
            }
        }
    }
}

fn run_case(seed: u64, tally: &mut Tally, write_widths: bool) {
    // One field in four moves along a ranged param, not t.
    let tv = if seed % 4 == 3 { Var::Param(0) } else { Var::T };
    let (g, root, seams) = generate(seed, tv);
    let tape = Tape::new(&g, root);
    // Forking every jump inside the last doubles the work per jump: a
    // build forking selects on one condition apart (as before they rode
    // one fork) takes minutes on some fields. Those are left out (of both
    // builds, so they compare).
    let most: u64 = std::env::var("INK_FUZZ_STEPS").ok().and_then(|s| s.parse().ok()).unwrap_or(1 << 20);
    if tape.forked_steps() > most {
        tally.skipped += 1;
        return;
    }
    let mut rng = Rng(seed ^ 0xABCD_EF01);
    tally.cases += 1;
    let inputs = if tv == Var::T { Inputs { t: true, params: Vec::new() } } else { Inputs { t: false, params: vec![0] } };
    let n = tape.ssa.len();
    let mut fail = Vec::new();
    for s in spans(&mut rng, &seams) {
        let ctx = if tv == Var::T {
            Ctx { t: span_mid(s), t_span: Some(s), ..Ctx::new(s.0, vec![]) }
        } else {
            Ctx { p_spans: ink_core::eval::ParamSpans::new([(0, s)]), ..Ctx::new(0.0, vec![span_mid(s)]) }
        };
        let ts = moments(&mut rng, s, &seams);
        for b in boxes(&mut rng) {
            let pts = points(&mut rng, b);
            let vals = sample(tv, &tape, &ts, &pts);
            let before = ink_core::eval::forked_count();
            let v = eval_interval_scope(&tape, &ctx, b, None);
            check("scope", v, I::ENTIRE, I::ENTIRE, f64::INFINITY, s, &ts, &vals, &mut fail);
            let (v2, d, j) = eval_interval_rates_jump(&tape, &ctx, b, &inputs);
            let d = d.first().copied().unwrap_or(I::point(0.0));
            check("rates", v2, v2, d, j, s, &ts, &vals, &mut fail);
            let sc = eval_interval_scope_every(&tape, &ctx, &inputs, b, true);
            check("every-out", sc.v[n - 1], sc.m[n - 1], sc.d[n - 1], sc.j[n - 1], s, &ts, &vals, &mut fail);
            // Some inner values, each against its own cone.
            for _ in 0..6 {
                let k = rng.below(n);
                if matches!(tape.ssa[k].op, ink_core::tape::Op::Const(_)) {
                    continue;
                }
                let cone = tape.cone(k);
                let cv = sample(tv, &cone, &ts, &pts);
                check(&format!("every[{k}]"), sc.v[k], sc.m[k], sc.d[k], sc.j[k], s, &ts, &cv, &mut fail);
                tally.checks += 1;
            }
            tally.forked += (ink_core::eval::forked_count() - before) as usize;
            tally.checks += 3;
            if write_widths {
                _ = writeln!(tally.widths, "{seed} {:?} {:?} {} {} {} {}", s, b.map(|x| (x.lo, x.hi)), v.hi - v.lo, v2.hi - v2.lo, mag(d), j);
            }
            if !fail.is_empty() {
                break;
            }
        }
        if !fail.is_empty() {
            let mut e = String::new();
            show(&g, root, 0, &mut e);
            tally.failures.push(format!("seed {seed} span {s:?}: {}\n  field: {e}", fail.join("; ")));
            return;
        }
    }
    // STILL-BEGIN
    // Whether it moves at all, over all of space and the whole span.
    if tv != Var::T {
        return;
    }
    let span = (0.0, 1.0);
    let ctx = Ctx::new(0.0, vec![]);
    let still = ink_kernel::motion::still(&tape, &ctx, [I::ENTIRE; 3], span);
    let rate0 = ink_kernel::motion::rate(&tape, &ctx, [I::ENTIRE; 3], span) == Some(0.0);
    tally.still += usize::from(still);
    tally.still_not_rate0 += usize::from(still && !rate0);
    tally.rate0_not_still += usize::from(rate0 && !still);
    if still || rate0 {
        let b = [I::new(-20.0, 20.0); 3];
        let pts = points(&mut rng, b);
        let ts = moments(&mut rng, span, &seams);
        let vals = sample(tv, &tape, &ts, &pts);
        for pi in 0..pts.len() {
            if ts.iter().enumerate().any(|(ti, _)| vals[ti][pi].to_bits() != vals[0][pi].to_bits() && (vals[ti][pi] - vals[0][pi]).abs() > slack(vals[0][pi])) {
                let mut e = String::new();
                show(&g, root, 0, &mut e);
                tally.failures.push(format!("seed {seed}: said still (still {still}, rate 0 {rate0}) but moves at point {:?}\n  field: {e}", pts[pi]));
                break;
            }
        }
    }
    // STILL-END
}

#[test]
fn interval_bounds_over_steps_of_t_hold_every_sampled_value() {
    let cases: u64 = std::env::var("INK_FUZZ_CASES").ok().and_then(|s| s.parse().ok()).unwrap_or(1200);
    let seed0: u64 = std::env::var("INK_FUZZ_SEED").ok().and_then(|s| s.parse().ok()).unwrap_or(1);
    let out = std::env::var("INK_FUZZ_OUT").ok();
    use rayon::prelude::*;
    let tallies: Vec<Tally> = (seed0..seed0 + cases)
        .into_par_iter()
        .fold(Tally::default, |mut t, seed| {
            let start = std::time::Instant::now();
            run_case(seed, &mut t, out.is_some());
            let took = start.elapsed().as_secs_f64();
            if took > 2.0 {
                eprintln!("slow: seed {seed} took {took:.1} s");
            }
            t
        })
        .collect();
    let mut all = Tally::default();
    for t in tallies {
        all.cases += t.cases;
        all.checks += t.checks;
        all.forked += t.forked;
        all.still += t.still;
        all.still_not_rate0 += t.still_not_rate0;
        all.rate0_not_still += t.rate0_not_still;
        all.skipped += t.skipped;
        all.failures.extend(t.failures);
        all.widths.push_str(&t.widths);
    }
    if let Some(f) = out {
        std::fs::write(f, &all.widths).unwrap();
    }
    eprintln!(
        "{} fields ({} more skipped, forking too much), {} bounds checked, {} forks; still {} (still but rate not 0: {}; rate 0 but not still: {}); {} failing",
        all.cases,
        all.skipped,
        all.checks,
        all.forked,
        all.still,
        all.still_not_rate0,
        all.rate0_not_still,
        all.failures.len()
    );
    for f in all.failures.iter().take(20) {
        eprintln!("{f}");
    }
    assert!(all.failures.is_empty(), "{} fields fail", all.failures.len());
}
