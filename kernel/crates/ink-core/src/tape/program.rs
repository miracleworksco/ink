//! Tape storage and access to instruction and input metadata.

use super::joins::Joins;
use super::promise::Promises;
use super::{Ins, Op, Pass};
use crate::graph::Var;

#[derive(Clone, Debug)]
pub struct Tape {
    pub ssa: Vec<Ins>,
    pub code: Vec<Ins>,
    pub nslots: usize,
    /// Sampled fields referenced by Op::Grid (the graph's, shared).
    pub grids: Vec<std::sync::Arc<crate::grid::GridData>>,
    /// Tag names provenance stops at (graph::Graph::opaque).
    pub opaque: Vec<u32>,
    /// Per instruction: its value depends only on constants and params (in
    /// the model, not just in a simplified tape, where a picked branch can
    /// leave spatial arithmetic with constant operands). Params are fixed for
    /// a query, so interval evaluation takes these as exact points computed
    /// in floats, as load-time defs are; unless one reads a ranged param
    /// (eval::ParamSpans), whose value it means in real numbers.
    pub fixed: Vec<bool>,
    /// Per fixed instruction: bit k % 64 of each param k it reads (in the
    /// model, as `fixed`): a prefilter, where a context's ranged params'
    /// bits are not set it reads none of them (two params 64 apart share a
    /// bit).
    pub pmask: Vec<u64>,
    /// Per fixed instruction that reads params, which ones, exactly (in the
    /// model, as `fixed`): (instruction, start, len) by instruction, the
    /// params sorted in `plist[start..start + len]` (`params_of`), a list
    /// shared where it is an argument's. Read only where `pmask` meets a
    /// context's ranged params: a ranged param this reads takes the
    /// instruction off the float path, one sharing only its bit does not.
    pub pread: Vec<(u32, u32, u32)>,
    pub plist: Vec<u32>,
    /// Per instruction: whether a derivative's range reads its gradient
    /// bounds (it is in the cone of a field an Op::Deriv bounds by). Empty
    /// where none does.
    pub slopes: Vec<bool>,
    /// The tapes of the read forms Op::Read names (graph::Graph::read_form).
    pub forms: Vec<std::sync::Arc<Tape>>,
    /// The instructions that can jump (floor, atan2, a select whose
    /// condition reads no space), in order, which interval evaluation may
    /// take as a choice of two sides (eval::forked). Derived from `ssa`.
    pub jumps: Vec<u32>,
    /// Where jumps join (Tape::rejoin), derived from `ssa` when a jump is
    /// first straddled.
    pub(super) joins: std::sync::OnceLock<Joins>,
    /// Each promise's operands (Tape::promise) by its instruction, each
    /// found when first read.
    pub(super) promises: std::sync::OnceLock<Promises>,
    /// Per instruction, bit k of each spatial axis k (x, y, z) its value
    /// reads (Tape::axes_read), derived from `ssa` when first read.
    pub(super) axes: std::sync::OnceLock<Vec<u8>>,
    /// Where pieces may meet (crate::sign::levels), derived from `ssa` when
    /// first read.
    pub(crate) levels: std::sync::OnceLock<Vec<(u32, f64)>>,
    /// Per min or max of a and a + c (or a - c), through names, in order:
    /// the instruction, the slot holding c there (kept until it), and
    /// whether the second operand less the first is -c (`Tape::apart`).
    /// Derived from `ssa`.
    pub apart: Vec<(u32, u32, bool)>,
    /// Per remainder a - p floor(a / p) (std's mod), in order: the
    /// instruction and the slot holding p there (kept until it)
    /// (`Tape::remainder`). Derived from `ssa`.
    pub remainders: Vec<(u32, u32)>,
}

impl Tape {
    /// Of min or max instruction `i` whose operands are a and a + c or
    /// a - c (`apart_of`): the slot holding c, and whether the second less
    /// the first is -c.
    pub fn apart(&self, i: usize) -> Option<(u32, bool)> {
        if self.apart.is_empty() {
            return None;
        }
        self.apart.binary_search_by_key(&(i as u32), |x| x.0).ok().map(|k| (self.apart[k].1, self.apart[k].2))
    }

    /// Of remainder instruction `i` (`remainder_of`): the slot holding p.
    pub fn remainder(&self, i: usize) -> Option<u32> {
        if self.remainders.is_empty() {
            return None;
        }
        self.remainders.binary_search_by_key(&(i as u32), |x| x.0).ok().map(|k| self.remainders[k].1)
    }

    /// Which of x, y and z the field's value reads, not through what only
    /// bounds it (a Lipschitz promise's constant and frame, a derivative's
    /// mark), as graph::Graph::value_depends_on: along an axis it reads
    /// not, the field is the same all along (a 2D profile's z).
    pub fn space_read(&self) -> [bool; 3] {
        let r = self.axes()[self.ssa.len() - 1];
        std::array::from_fn(|k| r & 1 << k != 0)
    }

    /// The space axes each value reads (as space_read), a bit each.
    pub(super) fn axes(&self) -> Vec<u8> {
        let mut reads = vec![0u8; self.ssa.len()];
        for i in 0..self.ssa.len() {
            reads[i] = match self.ssa[i].op {
                Op::Var(v) => [Var::X, Var::Y, Var::Z].iter().position(|&w| w == v).map_or(0, |a| 1 << a),
                Op::Lipschitz | Op::Deriv(_) => reads[self.ssa[i].a as usize],
                Op::CornerReach(fa, fb) => {
                    let (ra, rb) = (self.forms[fa as usize].space_read(), self.forms[fb as usize].space_read());
                    let fr = &self.ssa[self.ssa[i].b as usize];
                    let fc = [fr.a, fr.b, fr.c];
                    let mut m = reads[self.ssa[i].a as usize] | reads[self.ssa[i].c as usize];
                    for k in 0..3 {
                        if ra[k] || rb[k] {
                            m |= reads[fc[k] as usize];
                        }
                    }
                    m
                }
                _ => {
                    let (args, m) = Self::reads(&self.ssa, i);
                    args[..m].iter().fold(0, |r, &a| r | reads[a as usize])
                }
            };
        }
        reads
    }

    /// Bit k of each spatial axis k (x, y, z) instruction `i`'s value reads.
    pub fn axes_read(&self, i: usize) -> u8 {
        self.axes.get_or_init(|| {
            let mut m: Vec<u8> = Vec::with_capacity(self.ssa.len());
            for (j, ins) in self.ssa.iter().enumerate() {
                let own = match ins.op {
                    Op::Var(Var::X) => 1,
                    Op::Var(Var::Y) => 2,
                    Op::Var(Var::Z) => 4,
                    _ => 0,
                };
                let (args, k) = Self::reads(&self.ssa, j);
                m.push(args[..k].iter().fold(own, |a, &x| a | m[x as usize]));
            }
            m
        })[i]
    }

    /// The params fixed instruction `i` reads, exactly, sorted (`pread`):
    /// none for any other.
    pub fn params_of(&self, i: usize) -> &[u32] {
        match self.pread.binary_search_by_key(&(i as u32), |e| e.0) {
            Ok(j) => {
                let (_, start, len) = self.pread[j];
                &self.plist[start as usize..(start + len) as usize]
            }
            Err(_) => &[],
        }
    }

    /// The params the tape's value depends on (its forms' included), each
    /// once: those its instructions read as evaluation reads them
    /// (`pread`), which a simplified tape keeps where it folded a param
    /// away (a ranged one still takes the arithmetic it fed in real
    /// numbers).
    pub fn params_read(&self) -> Vec<u32> {
        let mut out = self.plist.clone();
        out.extend(self.forms.iter().flat_map(|f| f.params_read()));
        out.sort_unstable();
        out.dedup();
        out
    }

    pub fn len(&self) -> usize {
        self.ssa.len()
    }

    /// The machine instructions a pass of it retires (`Op::work`).
    pub fn work(&self, pass: Pass) -> usize {
        self.ssa.iter().map(|i| i.op.work(pass)).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.ssa.is_empty()
    }

    pub fn out_slot(&self) -> usize {
        self.code.last().unwrap().out as usize
    }

    /// The bytes the tape holds, about: its lists, not what is derived
    /// from them when first read.
    pub fn bytes(&self) -> usize {
        use std::mem::size_of_val;
        let Tape { ssa, code, nslots: _, grids, opaque, fixed, pmask, pread, plist, slopes, forms, jumps, joins: _, promises: _, axes: _, levels: _, apart, remainders } = self;
        size_of_val(&**ssa)
            + size_of_val(&**code)
            + size_of_val(&**grids)
            + size_of_val(&**opaque)
            + size_of_val(&**fixed)
            + size_of_val(&**pmask)
            + size_of_val(&**pread)
            + size_of_val(&**plist)
            + size_of_val(&**slopes)
            + size_of_val(&**forms)
            + size_of_val(&**jumps)
            + size_of_val(&**apart)
            + size_of_val(&**remainders)
    }
}
