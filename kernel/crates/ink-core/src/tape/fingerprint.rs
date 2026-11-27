//! Tape identity for cached answers and for sharing derived tapes.

use super::{Ins, Op, Tape};
use crate::graph::{Fingerprint, Var};

/// A question's names: the tags and params its tapes read, numbered as they
/// first appear in them (instructions, then names provenance stops at, then
/// forms), then in the equations of the derivations they read
/// (Ctx::read_by: whether t is read counts those too), then the context's
/// other ranged params. The model's own numbers
/// say where it declared or imported each, so a cached answer is keyed by
/// these (`Tape::fingerprint`, `Ctx::fingerprint`): an unchanged part keeps
/// its key in any model that imports it, wherever it does. An answer that
/// names a tag or param by the model's number also hashes those (`raw`).
#[derive(Clone, Debug, Default)]
pub struct Names {
    tags: Vec<u32>,
    params: Vec<u32>,
    /// The context's ranged params, the model's numbers, sorted.
    pub ranged: Vec<u32>,
    /// Whether a tape reads t: else the moment is not in the question (a
    /// part that does not move has one answer over a motion).
    pub reads_t: bool,
}

impl Names {
    pub fn of(tapes: &[&Tape], ctx: &crate::eval::Ctx) -> Names {
        fn add(list: &mut Vec<u32>, v: u32) {
            if !list.contains(&v) {
                list.push(v);
            }
        }
        fn walk(t: &Tape, n: &mut Names) {
            for i in &t.ssa {
                match i.op {
                    Op::Tag(g) => add(&mut n.tags, g),
                    Op::Var(Var::Param(k)) => add(&mut n.params, k),
                    Op::Var(Var::T) => n.reads_t = true,
                    _ => {}
                }
            }
            t.opaque.iter().for_each(|&g| add(&mut n.tags, g));
            t.forms.iter().for_each(|f| walk(f, n));
        }
        let mut n = Names { ranged: ctx.p_spans.iter().map(|s| s.k).collect(), ..Names::default() };
        n.ranged.sort_unstable();
        tapes.iter().for_each(|t| walk(t, &mut n));
        ctx.read_by(tapes).1.iter().flat_map(|d| d.tapes()).for_each(|t| walk(t, &mut n));
        ctx.p_spans.iter().for_each(|s| add(&mut n.params, s.k));
        n
    }

    /// Tag t's number in the question.
    pub fn tag(&self, t: u32) -> u32 {
        self.tags.iter().position(|&g| g == t).map_or(u32::MAX, |i| i as u32)
    }

    /// Param k's number in the question.
    pub fn param(&self, k: u32) -> u32 {
        self.params.iter().position(|&p| p == k).map_or(u32::MAX, |i| i as u32)
    }

    /// The params in the question's order (the model's numbers).
    pub fn params(&self) -> &[u32] {
        &self.params
    }

    /// The model's numbers of the tags and params, in the question's order.
    pub fn raw(&self, fp: &mut Fingerprint) {
        fp.word(self.tags.len() as u64).word(self.params.len() as u64);
        self.tags.iter().chain(&self.params).for_each(|&v| _ = fp.word(v as u64));
    }
}

impl Tape {
    /// Hash all of the tape into `fp`: every instruction of both views,
    /// the grids and forms it reads, and what each instruction is marked
    /// with, its tags and params by their numbers in the question (`Names`).
    /// Equal hashes mean the same tape, instruction for instruction, up to
    /// the model's numbering of names (and hash collisions), so the same
    /// answer to any query on it in the same context: for keys of cached
    /// answers (ink_kernel::cache).
    pub fn fingerprint(&self, fp: &mut Fingerprint, names: &Names) {
        // Every field, named: one added must be hashed too.
        // `jumps`, `joins`, `promises`, `axes`, `levels`, `apart` and `remainders` are derived from `ssa`; `pmask`, `pread` and
        // `plist` enter as what they decide (below).
        let Tape {
            ssa,
            code,
            nslots,
            grids,
            opaque,
            fixed,
            pmask: _,
            pread: _,
            plist: _,
            slopes,
            forms,
            jumps: _,
            joins: _,
            promises: _,
            axes: _,
            levels: _,
            apart: _,
            remainders: _,
        } = self;
        for list in [ssa, code] {
            fp.word(list.len() as u64);
            for &Ins { op, out, a, b, c } in list {
                match op {
                    Op::Const(v) => fp.word(0).f64(v),
                    Op::Var(v) => fp.word(1).word(match v {
                        Var::X => 0,
                        Var::Y => 1,
                        Var::Z => 2,
                        Var::T => 3,
                        Var::Param(k) => 4 + names.param(k) as u64,
                    }),
                    Op::Un(u) => fp.word(2).word(u as u64),
                    Op::Bin(b) => fp.word(3).word(b as u64),
                    Op::Select => fp.word(4),
                    Op::Tag(t) => fp.word(5).word(names.tag(t) as u64),
                    Op::Grid(g, m) => fp.word(6).word(g as u64).word(m as u64),
                    Op::Read => fp.word(7),
                    Op::Lipschitz => fp.word(8),
                    Op::Frame => fp.word(9),
                    Op::Deriv(k) => fp.word(10).word(k as u64),
                    Op::Real => fp.word(11),
                    Op::Within(id) => fp.word(12).word(id),
                    Op::CornerReach(fa, fb) => fp.word(13).word(fa as u64).word(fb as u64),
                };
                fp.word(out as u64).word(a as u64).word(b as u64).word(c as u64);
            }
        }
        fp.word(*nslots as u64).word(grids.len() as u64);
        grids.iter().for_each(|g| g.fingerprint(fp));
        fp.word(opaque.len() as u64);
        opaque.iter().for_each(|&t| _ = fp.word(names.tag(t) as u64));
        fp.word(fixed.len() as u64);
        fixed.iter().for_each(|f| _ = fp.word(*f as u64));
        // Which fixed instructions read a ranged param, exactly
        // (Tape::pread): those it takes in real numbers.
        for i in (0..fixed.len()).filter(|&i| fixed[i]) {
            fp.word(self.params_of(i).iter().any(|k| names.ranged.binary_search(k).is_ok()) as u64);
        }
        fp.word(slopes.len() as u64);
        slopes.iter().for_each(|s| _ = fp.word(*s as u64));
        fp.word(forms.len() as u64);
        forms.iter().for_each(|f| f.fingerprint(fp, names));
    }

    /// An instruction's op as two words: its kind and what it carries (a
    /// constant to the bit).
    fn op_words(op: Op) -> (u64, u64) {
        match op {
            Op::Const(v) => (0, v.to_bits()),
            Op::Var(v) => (
                1,
                match v {
                    Var::X => 0,
                    Var::Y => 1,
                    Var::Z => 2,
                    Var::T => 3,
                    Var::Param(k) => 4 + k as u64,
                },
            ),
            Op::Un(u) => (2, u as u64),
            Op::Bin(b) => (3, b as u64),
            Op::Select => (4, 0),
            Op::Tag(t) => (5, t as u64),
            Op::Grid(g, m) => (6, (g as u64) << 8 | m as u64),
            Op::Read => (7, 0),
            Op::Lipschitz => (8, 0),
            Op::Frame => (9, 0),
            Op::Deriv(d) => (10, d as u64),
            Op::Real => (11, 0),
            Op::Within(id) => (12, id),
            Op::CornerReach(fa, fb) => (13, (fa as u64) << 32 | fb as u64),
        }
    }

    /// A hash of the lists `same` compares, for this process alone (its
    /// grids and forms by address).
    pub fn content(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut h = rustc_hash::FxHasher::default();
        let Tape { ssa, grids, opaque, fixed, pmask, pread, plist, slopes, forms, .. } = self;
        h.write_usize(ssa.len());
        for &Ins { op, out, a, b, c } in ssa {
            let (k, v) = Tape::op_words(op);
            [k, v, out as u64, a as u64, b as u64, c as u64].into_iter().for_each(|w| h.write_u64(w));
        }
        grids.iter().for_each(|g| h.write_usize(std::sync::Arc::as_ptr(g) as usize));
        forms.iter().for_each(|f| h.write_usize(std::sync::Arc::as_ptr(f) as usize));
        (opaque, fixed, pmask, pread, plist, slopes).hash(&mut h);
        h.finish()
    }

    /// Whether two tapes are one field, evaluated alike: every list equal,
    /// each constant to the bit, the same grids and forms.
    pub fn same(&self, other: &Tape) -> bool {
        // Every field, named: one added must be compared too. `code`,
        // `nslots`, `jumps`, `apart`, `remainders` and what is derived when
        // first read follow from `ssa`.
        let Tape {
            ssa,
            code: _,
            nslots: _,
            grids,
            opaque,
            fixed,
            pmask,
            pread,
            plist,
            slopes,
            forms,
            jumps: _,
            joins: _,
            promises: _,
            axes: _,
            levels: _,
            apart: _,
            remainders: _,
        } = self;
        let ins = |i: &Ins| (Tape::op_words(i.op), i.out, i.a, i.b, i.c);
        fn arcs<T>(a: &[std::sync::Arc<T>], b: &[std::sync::Arc<T>]) -> bool {
            a.len() == b.len() && a.iter().zip(b).all(|(x, y)| std::sync::Arc::ptr_eq(x, y))
        }
        ssa.len() == other.ssa.len()
            && ssa.iter().zip(&other.ssa).all(|(x, y)| ins(x) == ins(y))
            && arcs(grids, &other.grids)
            && arcs(forms, &other.forms)
            && (opaque, fixed, pmask, pread, plist, slopes) == (&other.opaque, &other.fixed, &other.pmask, &other.pread, &other.plist, &other.slopes)
    }
}
