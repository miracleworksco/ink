//! Operands and constant checks for Lipschitz promises.

use super::{Ins, Op, Tape};
use crate::graph::{UnOp, Var};

#[derive(Clone, Debug)]
pub(super) struct Promises {
    operands: Box<[(u32, std::sync::OnceLock<Promise>)]>,
    always_real: Vec<bool>,
}

impl Promises {
    fn of(ssa: &[Ins]) -> Self {
        let (mut spatial, mut always_real) = (vec![false; ssa.len()], vec![false; ssa.len()]);
        let mut operands = Vec::new();
        for (i, ins) in ssa.iter().enumerate() {
            let (args, n) = ins.args();
            let space = |a: u32| spatial[a as usize];
            let always_spatial = match ins.op {
                Op::Var(Var::X | Var::Y | Var::Z) => true,
                Op::Select => space(ins.a) || (space(ins.b) && space(ins.c)),
                _ => args[..n].iter().copied().any(space),
            };
            spatial[i] = always_spatial;
            always_real[i] = match ins.op {
                Op::Const(v) => v.is_finite(),
                Op::Var(v) => matches!(v, Var::X | Var::Y | Var::Z),
                _ => {
                    // A spatial operation cannot introduce an unreal
                    // constant. These other operators keep finite values
                    // finite even when every argument is spatially constant.
                    let finite =
                        matches!(ins.op, Op::Select | Op::Tag(_) | Op::Lipschitz | Op::Read | Op::Deriv(_) | Op::Real | Op::Frame | Op::Within(_) | Op::Un(UnOp::Neg | UnOp::Abs));
                    (always_spatial || finite) && args[..n].iter().all(|&a| always_real[a as usize])
                }
            };
            if matches!(ins.op, Op::Lipschitz) {
                operands.push((i as u32, Default::default()));
            }
        }
        Self { operands: operands.into_boxed_slice(), always_real }
    }
}

/// What interval evaluation reads again of a promise, Lipschitz instruction
/// `i` of the SSA form (eval::promise_valid): the instructions its field,
/// its constant and its frame read, by their indices in the SSA form, in
/// order; which are a select's condition's (`cond`). Indices, not copies:
/// a field of many promises over one shared frame (a wheel's 424, each
/// reading the same 11,000 instructions of coordinates) would hold each
/// copy once per promise, and every tape simplified over a box its own.
#[derive(Clone, Debug)]
pub struct Promise {
    pub ins: Vec<u32>,
    pub cond: Vec<bool>,
    /// Whether k reads x, y or z and a coordinate of the frame passes through a jump (a select, a
    /// floor): the frame's coordinates over a box need not be all that lies between the centre's and
    /// a point's, so k's bound there (read at the image) need not bound the slope along the way.
    pub leaps: bool,
}

impl Promise {
    fn of(ssa: &[Ins], i: usize) -> Promise {
        PROMISE_SCRATCH.with_borrow_mut(|scratch| scratch.build(ssa, i))
    }
}

/// Sparse work buffers for extracting a promise's operands. Only positions
/// touched by the cone are cleared between promises.
#[derive(Default)]
struct PromiseScratch {
    positions: Vec<u32>,
    stack: Vec<u32>,
    order: Vec<u32>,
    uses: Vec<u8>,
    /// A build began and did not end (a panic unwound through it): the
    /// buffers hold its marks, and the next build clears them first.
    dirty: bool,
}

impl PromiseScratch {
    fn build(&mut self, ssa: &[Ins], i: usize) -> Promise {
        let (l, fr) = (ssa[i], ssa[ssa[i].c as usize]);
        let frame = [fr.a, fr.b, fr.c];
        // Discover the operands once. Their order in the SSA is already a
        // topological order, so all three uses of the cone can be marked
        // together in a backwards pass over the compact instructions.
        const CONDITION: u8 = 1;
        const OVER: u8 = 2;
        let PromiseScratch { positions, stack, order, uses, dirty } = self;
        if std::mem::replace(dirty, true) {
            positions.fill(u32::MAX);
            stack.clear();
            order.clear();
            uses.clear();
        }
        if positions.len() < ssa.len() {
            positions.resize(ssa.len(), u32::MAX);
        }
        stack.extend([l.a, l.b, l.c]);
        while let Some(j) = stack.pop() {
            if std::mem::replace(&mut positions[j as usize], 0) == u32::MAX {
                order.push(j);
                let (args, k) = ssa[j as usize].args();
                stack.extend_from_slice(&args[..k]);
            }
        }
        order.sort_unstable();
        for (k, &j) in order.iter().enumerate() {
            positions[j as usize] = k as u32;
        }
        let pos = |j: u32| positions[j as usize];
        let ins: Vec<_> = order.iter().map(|&j| (j, renumber(ssa[j as usize], pos))).collect();
        let promise = renumber(l, pos);
        let local_frame = frame.map(pos);
        uses.resize(ins.len(), 0);
        uses[promise.b as usize] |= OVER;
        for j in local_frame {
            uses[j as usize] |= OVER;
        }
        for (n, &(_, i)) in ins.iter().enumerate().rev() {
            let read = uses[n];
            let (args, k) = i.args();
            for &arg in &args[..k] {
                uses[arg as usize] |= read;
            }
            if matches!(i.op, Op::Select) {
                uses[i.a as usize] |= CONDITION;
            }
        }
        let cond = uses.iter().map(|&u| u & CONDITION != 0).collect();
        // Whether k reads x, y or z, and a frame coordinate a jump, over the cone of just k and the
        // coordinates (each instruction's arguments come before it).
        let (mut spatial, mut jumps) = (vec![false; ins.len()], vec![false; ins.len()]);
        for (n, &(_, i)) in ins.iter().enumerate() {
            if uses[n] & OVER == 0 {
                continue;
            }
            let (args, k) = i.args();
            spatial[n] = matches!(i.op, Op::Var(Var::X | Var::Y | Var::Z)) || args[..k].iter().any(|&a| spatial[a as usize]);
            jumps[n] = matches!(i.op, Op::Select | Op::Un(UnOp::Floor)) || args[..k].iter().any(|&a| jumps[a as usize]);
        }
        let leaps = spatial[promise.b as usize] && local_frame.iter().any(|&j| jumps[j as usize]);
        for j in order.drain(..) {
            positions[j as usize] = u32::MAX;
        }
        uses.clear();
        *dirty = false;
        Promise { ins: ins.into_iter().map(|(j, _)| j).collect(), cond, leaps }
    }
}

thread_local! {
    static PROMISE_SCRATCH: std::cell::RefCell<PromiseScratch> = Default::default();
}

fn renumber(mut ins: Ins, pos: impl Fn(u32) -> u32) -> Ins {
    let (_, n) = ins.args();
    for a in [&mut ins.a, &mut ins.b, &mut ins.c].into_iter().take(n) {
        *a = pos(*a);
    }
    ins
}

impl Tape {
    /// Lipschitz instruction `i`'s operands as interval evaluation reads
    /// them again (Promise), found once.
    pub fn promise(&self, i: usize) -> &Promise {
        let all = &self.promises.get_or_init(|| Promises::of(&self.ssa)).operands;
        let k = all.binary_search_by_key(&(i as u32), |p| p.0).expect("a promise");
        all[k].1.get_or_init(|| Promise::of(&self.ssa, i))
    }

    /// The constants of this promise are real in every region and context.
    /// When the tape cannot show that, the evaluator checks them per box.
    pub(crate) fn promise_always_real(&self, i: usize) -> bool {
        self.promises.get_or_init(|| Promises::of(&self.ssa)).always_real[i]
    }
}
