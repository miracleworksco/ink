//! Jump rejoin points and how the output moves with an instruction's value.

use super::{Ins, Op, Tape};
use crate::graph::{BinOp, UnOp, Var};

/// Where each of a tape's jumps joins (Tape::rejoin). An instruction is
/// monotone where its value reaches the output only through operations
/// that each read it once and whose bounds are their operand's ends mapped
/// (sums, min, max, a monotone function), the output moving the same way
/// with it through every read: a bound of the output from the hull of two
/// bounds of it is then the hull of the bounds from each (each end of the
/// output's bound is the end of every read that moves it, all from the
/// same end of the hull), so nothing is lost by joining them there
/// (eval::forked). A value read again through a negation is not: a
/// screw's thread, and the same thread cut from its nut.
#[derive(Clone, Debug)]
pub(super) struct Joins {
    /// Per jump (Tape::jumps): Tape::rejoin's instruction; 0 for others.
    rejoin: Vec<u32>,
    /// Tape::forked_steps.
    forked: u64,
}

/// The first instruction from `j` on that every later use of j's value
/// passes through (no other value depending on j is still to be read, as
/// interval evaluation reads them: Tape::uses), and from which the output
/// is `monotone`. `dep` is scratch, all false.
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
        let (args, k) = Tape::uses(ssa, i);
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

/// Where each of `jumps` joins (Tape::rejoin; 0 for other instructions),
/// given the instructions from which the output is monotone.
pub(super) fn rejoins(ssa: &[Ins], jumps: &[u32], monotone: &[bool]) -> Vec<u32> {
    let n = ssa.len();
    let mut last = vec![0u32; n];
    for i in 0..n {
        let (args, k) = Tape::uses(ssa, i);
        for &a in &args[..k] {
            last[a as usize] = i as u32;
        }
    }
    last[n - 1] = u32::MAX;
    let (mut rejoin, mut dep) = (vec![0u32; n], vec![false; n]);
    for &j in jumps {
        rejoin[j as usize] = join_of(ssa, &last, monotone, &mut dep, j as usize) as u32;
    }
    rejoin
}

/// Tape::forked_steps over instructions `from..to`.
fn forked_steps(jumps: &[u32], rejoin: &[u32], from: usize, to: usize) -> u64 {
    let (mut steps, mut i) = (0u64, from);
    for &j in &jumps[jumps.partition_point(|&j| (j as usize) < from)..] {
        let (j, r) = (j as usize, rejoin[j as usize] as usize);
        if j >= to {
            break;
        }
        if j >= i && r > j && r < to {
            steps = steps.saturating_add((j + 1 - i) as u64).saturating_add(forked_steps(jumps, rejoin, j + 1, r + 1).saturating_mul(2));
            i = r + 1;
        }
    }
    steps.saturating_add((to - i) as u64)
}

/// Whether instruction `u`'s bound is its operand `a`'s ends mapped, with
/// its other operands fixed (so a bound of it from the hull of two bounds
/// of `a` is the hull of the bounds from each): sums, min, max, emax, a
/// select's branches, or its condition where its branches are numbers (one
/// where the condition is below 0 over the bound, the other at or above,
/// else both), a monotone function, a product with a constant; `a` read
/// once.
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
        Op::Select => ins.a != a || (constant(ins.b) && constant(ins.c)),
        Op::Lipschitz => ins.a == a,
        _ => false,
    }
}

/// Which way instruction `u`'s value moves with its operand `a`, the
/// others held: 1 with it, -1 against it, 0 where its bound is not `a`'s
/// ends mapped (`ends_mapped`) or it moves both ways.
fn way_of(ssa: &[Ins], u: usize, a: u32) -> i8 {
    if !ends_mapped(ssa, u, a) {
        return 0;
    }
    let ins = &ssa[u];
    let value = |x: u32| match ssa[x as usize].op {
        Op::Const(c) => c,
        _ => 0.0,
    };
    let sign = |x: u32| match value(x) {
        c if c > 0.0 => 1,
        c if c < 0.0 => -1,
        _ => 0,
    };
    match ins.op {
        // From the first branch below 0 to the second at or above it.
        Op::Select if ins.a == a => match value(ins.c) - value(ins.b) {
            d if d > 0.0 => 1,
            d if d < 0.0 => -1,
            _ => 0,
        },
        Op::Bin(BinOp::Sub) if ins.b == a => -1,
        Op::Un(UnOp::Neg) => -1,
        Op::Bin(BinOp::Mul) => sign(if ins.a == a { ins.b } else { ins.a }),
        Op::Bin(BinOp::Div) => sign(ins.b),
        _ => 1,
    }
}

impl Tape {
    /// Backwards from the output, which way it moves with each value: the
    /// same way through every read, each by an operation whose bound is
    /// that operand's ends mapped (`ends_mapped`), or not monotone (0); none
    /// where the output does not read it. A reader's way is final before its
    /// operands'.
    fn ways(&self) -> Vec<Option<i8>> {
        let (ssa, n) = (&self.ssa, self.ssa.len());
        let mut way = vec![None::<i8>; n];
        if n > 0 {
            way[n - 1] = Some(1);
        }
        for u in (0..n).rev() {
            let w = way[u].unwrap_or(0);
            let (args, k) = Self::uses(ssa, u);
            for &a in &args[..k] {
                let by = if w == 0 { 0 } else { w * way_of(ssa, u, a) };
                let a = a as usize;
                way[a] = Some(match way[a] {
                    Some(x) if x != by => 0,
                    _ => by,
                });
            }
        }
        way
    }

    /// Which way the output moves with param k, at any values of the rest
    /// (`ways`): 1 with it, -1 against it, 0 either; None where it does not
    /// read it.
    pub fn param_way(&self, k: u32) -> Option<i8> {
        let ways = self.ways();
        self.ssa.iter().zip(&ways).filter(|(i, _)| matches!(i.op, Op::Var(Var::Param(p)) if p == k)).fold(None, |acc, (_, &w)| match (acc, w.unwrap_or(0)) {
            (None, w) => Some(w),
            (Some(a), w) if a == w => Some(a),
            _ => Some(0),
        })
    }

    fn joins(&self) -> &Joins {
        self.joins.get_or_init(|| {
            let ssa = &self.ssa;
            let monotone: Vec<bool> = self.ways().iter().map(|w| w.is_some_and(|w| w != 0)).collect();
            let rejoin = rejoins(ssa, &self.jumps, &monotone);
            let forked = forked_steps(&self.jumps, &rejoin, 0, ssa.len());
            Joins { rejoin, forked }
        })
    }

    /// Where the bounds of a jump at instruction `j` (Tape::jumps;
    /// eval::forked) can be joined again: the first instruction from j on
    /// through which every later use of j's value passes, and from which the
    /// output is monotone (Joins). The output at the latest.
    pub fn rejoin(&self, j: usize) -> usize {
        self.joins().rejoin[j] as usize
    }

    /// The steps a run over the tape takes where every jump is straddled
    /// and forked (eval::forked: each window bounded twice, the forks
    /// inside it each time), saturating.
    pub fn forked_steps(&self) -> u64 {
        self.joins().forked
    }
}
