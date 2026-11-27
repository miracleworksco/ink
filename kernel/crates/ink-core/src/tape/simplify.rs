//! Derived tapes and branch choices for restricting evaluation.

use super::{Ins, Op, Tape, kept};
use crate::graph::{BinOp, UnOp};
use crate::interval::I;

/// The most a box's own tape may be of its parent's for the box to keep it (`Tape::narrowed_much`);
/// longer, it shares the parent's. From one measurement: an AndyMark wheel's mass, whose tapes
/// simplifying shrinks little, held 39 GB keeping every box's own tape (at most all of the parent's),
/// 24.9, 21.8 and 7.5 GB at most three quarters, a half and a quarter (438, 672 and 879 s), and took
/// over 20 minutes always sharing; the catalog's sample took 704, 683 and 701 s at all, a half and a
/// quarter (the same answers), and 1,392 s always sharing. A quarter is the most that kept that mass
/// under 16 GB.
pub const OWN_TAPE_AT_MOST: f64 = 0.25;

/// Choice codes recorded by interval evaluation for min/max/select.
pub const BOTH: u8 = 0;
pub const FIRST: u8 = 1;
pub const SECOND: u8 = 2;

impl Tape {
    /// Choices that keep the output's sign, not its value, over a region
    /// where instruction i lies in `vals[i]` (sound bounds): at a min, max
    /// or emax whose sign reaches the output's (through mins, maxes, emaxes,
    /// whose sign is their max's, negations, tags and products with a
    /// factor above 0 throughout: each maps where a value is below, at or
    /// above 0 to where the next is), the branch above 0 throughout dropped
    /// from a min and the one below 0 throughout from a max, where the other
    /// is not so too: the other branch's sign is then the min's or max's (a
    /// slab's faces past a profile's walls drop out of an extrusion inside
    /// its span; inside both, neither is, so how deep a point lies is still
    /// read along both). The tape simplified by them has the output's sign
    /// at every point of the region (below, at or above 0), not its value.
    /// `choices` starts with none made; true where it makes one.
    pub fn sign_choices(&self, vals: &[I], choices: &mut [u8]) -> bool {
        let n = self.ssa.len();
        // Each value: 0 not reached yet (dead, if never), 1 its sign
        // reaching the output's, 2 not.
        let mut way = vec![0u8; n];
        way[n - 1] = 1;
        let mut any = false;
        for i in (0..n).rev() {
            let (ins, w) = (&self.ssa[i], way[i]);
            if w == 0 {
                continue;
            }
            let (a, b) = (vals[ins.a as usize], vals[ins.b as usize]);
            let pick = match ins.op {
                _ if w != 1 => None,
                Op::Bin(BinOp::Min) if b.lo > 0.0 && !(a.lo > 0.0) => Some(ins.a),
                Op::Bin(BinOp::Min) if a.lo > 0.0 && !(b.lo > 0.0) => Some(ins.b),
                Op::Bin(BinOp::Max | BinOp::Emax) if b.hi < 0.0 && !(a.hi < 0.0) => Some(ins.a),
                Op::Bin(BinOp::Max | BinOp::Emax) if a.hi < 0.0 && !(b.hi < 0.0) => Some(ins.b),
                _ => None,
            };
            if let Some(k) = pick {
                choices[i] = if k == ins.a { FIRST } else { SECOND };
                any = true;
                way[k as usize] = way[k as usize].max(1);
                continue;
            }
            let (args, k) = Self::reads(&self.ssa, i);
            for &x in &args[..k] {
                let keeps = w == 1
                    && match ins.op {
                        Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax) | Op::Un(UnOp::Neg) | Op::Tag(_) => true,
                        Op::Bin(BinOp::Mul) => ins.a != ins.b && vals[if x == ins.a { ins.b } else { ins.a } as usize].lo > 0.0,
                        _ => false,
                    };
                way[x as usize] = way[x as usize].max(if keeps { 1 } else { 2 });
            }
        }
        any
    }

    /// The values each point where the output is at most 0 has at most 0
    /// (`true`) or at least 0 (`false`), and the space axes each reads (a
    /// bit each): the output's terms joined by maxes (an emax is at most 0
    /// where its max is) and, negated, by mins. A term reading one axis
    /// alone bounds where along it the output can be at most 0 (a slab's
    /// span). Each value is visited once each way, so a term shared by
    /// many is listed once, and the walk is as long as the tape.
    pub fn conjuncts(&self) -> Vec<(usize, bool, u8)> {
        let n = self.ssa.len();
        let reads = self.axes();
        let mut seen = vec![[false; 2]; n];
        let (mut out, mut stack) = (Vec::new(), vec![(n - 1, true)]);
        while let Some((i, below)) = stack.pop() {
            if std::mem::replace(&mut seen[i][below as usize], true) {
                continue;
            }
            let ins = &self.ssa[i];
            match (ins.op, below) {
                (Op::Bin(BinOp::Max | BinOp::Emax), true) | (Op::Bin(BinOp::Min), false) => stack.extend([(ins.a as usize, below), (ins.b as usize, below)]),
                (Op::Un(UnOp::Neg), _) => stack.push((ins.a as usize, !below)),
                (Op::Tag(_), _) => stack.push((ins.a as usize, below)),
                _ => out.push((i, below, reads[i])),
            }
        }
        out
    }

    /// The operand min/max/select instruction `i` reduces to by `choices`
    /// (interval evaluation proved the other cannot win), if one.
    /// Value `j` as the tape simplified by `choices` reads it: past its
    /// names (tags), and through each min, max and select they decide to
    /// the branch it takes (`simplify` aliases it so).
    pub fn through(&self, choices: &[u8], mut j: usize) -> usize {
        loop {
            j = match (self.ssa[j].op, self.pick(choices, j)) {
                (Op::Tag(_), _) => self.ssa[j].a as usize,
                (_, Some(a)) => a as usize,
                _ => return j,
            };
        }
    }

    /// Value `j` as the tape simplified by `choices` holds it: through each
    /// min, max and select they decide to the branch it takes.
    pub fn picked(&self, choices: &[u8], mut j: usize) -> usize {
        while let Some(a) = self.pick(choices, j) {
            j = a as usize;
        }
        j
    }

    fn pick(&self, choices: &[u8], i: usize) -> Option<u32> {
        let ins = &self.ssa[i];
        match (ins.op, choices[i]) {
            (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), FIRST) => Some(ins.a),
            (Op::Bin(BinOp::Min | BinOp::Max | BinOp::Emax), SECOND) => Some(ins.b),
            (Op::Select, FIRST) => Some(ins.b),
            (Op::Select, SECOND) => Some(ins.c),
            _ => None,
        }
    }

    /// Per instruction, whether the output reads its value, through the
    /// branches `choices` leaves (the rest, wherever `choices` were found,
    /// change nothing the output is).
    pub fn live(&self, choices: &[u8]) -> Vec<bool> {
        let n = self.ssa.len();
        let mut live = vec![false; n];
        live[n - 1] = true;
        for i in (0..n).rev() {
            if !live[i] {
                continue;
            }
            if let Some(a) = self.pick(choices, i) {
                live[a as usize] = true;
            } else {
                let (args, k) = self.ssa[i].args();
                for &a in &args[..k] {
                    live[a as usize] = true;
                }
            }
        }
        live
    }

    /// The tape over a region where interval evaluation recorded `choices`
    /// on this one: itself where they pick no branch, else simplified
    /// (`simplify`), kept (`kept`): a region's parts that choose as it did
    /// share it.
    pub fn narrowed(self: &std::sync::Arc<Self>, choices: &[u8]) -> std::sync::Arc<Tape> {
        if choices.iter().all(|&c| c == BOTH) {
            return self.clone();
        }
        kept::Kept::global().get(self, kept::Ask::Narrowed, choices, || self.simplify(choices))
    }

    /// `narrowed`, where that is at most OWN_TAPE_AT_MOST of this tape; else this tape, shared. A search that keeps a tape for each box it has yet to test holds one per box: where
    /// simplifying over a box takes little away (a field whose parts all reach into it), a near copy.
    pub fn narrowed_much(self: &std::sync::Arc<Self>, choices: &[u8]) -> std::sync::Arc<Tape> {
        let sub = self.narrowed(choices);
        if std::sync::Arc::ptr_eq(&sub, self) || sub.len() as f64 > OWN_TAPE_AT_MOST * self.len() as f64 { self.clone() } else { sub }
    }

    /// `cone`, kept (`kept`).
    pub fn kept_cone(self: &std::sync::Arc<Self>, i: usize) -> std::sync::Arc<Tape> {
        kept::Kept::global().get(self, kept::Ask::Cone, &(i as u32).to_le_bytes(), || self.cone(i))
    }

    /// Drop branches of min/max/select that interval evaluation proved can't win.
    pub fn simplify(&self, choices: &[u8]) -> Tape {
        let n = self.ssa.len();
        let pick = |i: usize| self.pick(choices, i);
        let live = self.live(choices);

        let mut remap = vec![u32::MAX; n];
        let mut out: Vec<Ins> = Vec::with_capacity(n);
        for i in 0..n {
            if !live[i] {
                continue;
            }
            if let Some(a) = pick(i) {
                remap[i] = remap[a as usize];
                continue;
            }
            let mut ins = self.ssa[i];
            let (_, k) = ins.args();
            if k > 0 {
                ins.a = remap[ins.a as usize];
            }
            if k > 1 {
                ins.b = remap[ins.b as usize];
            }
            if k > 2 {
                ins.c = remap[ins.c as usize];
            }
            ins.out = out.len() as u32;
            remap[i] = ins.out;
            out.push(ins);
        }
        // If the root was aliased, everything after its target is dead.
        // Held in its own length, not its parent's: the many simplified
        // tapes of an octree::Local would hold their parents' many times.
        out.truncate(remap[n - 1] as usize + 1);
        out.shrink_to_fit();
        let mut t = Tape::from_ssa(out);
        // What is fixed stays as it was: arithmetic on a branch picked for
        // the region still means real arithmetic, not a constants-and-params
        // value in floats (select(x, c, 2)^2 for x < 0 is c^2 in reals); and
        // so do the params it reads (a branch picked over a param's span
        // keeps reading the param).
        let kept = || (0..n).filter(|&i| live[i] && pick(i).is_none() && (remap[i] as usize) < t.fixed.len());
        let restored: Vec<usize> = kept().collect();
        for &i in &restored {
            t.fixed[remap[i] as usize] = self.fixed[i];
            t.pmask[remap[i] as usize] = self.pmask[i];
        }
        (t.pread, t.plist) = self.reread(restored.into_iter().map(|i| (remap[i], i)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        if t.ssa.iter().any(|ins| matches!(ins.op, Op::Read | Op::CornerReach(..))) {
            t.forms = self.forms.clone();
        }
        t
    }

    /// The tape computing instruction `i` (an ssa index) alone: its cone of
    /// instructions, in time proportional to that cone.
    pub fn cone(&self, i: usize) -> Tape {
        self.cone_max(&[i])
    }

    /// The tape computing the max of instructions `roots` (ssa indices;
    /// one: that instruction alone, `cone`): their cones, then the max of
    /// their values, in time proportional to those cones.
    pub fn cone_max(&self, roots: &[usize]) -> Tape {
        self.cone_terms(&roots.iter().map(|&r| (r, true)).collect::<Vec<_>>())
    }

    /// The tape computing the max of terms `roots` (ssa indices, each
    /// negated where its flag is false: a term of `conjuncts`), as
    /// `cone_max`.
    pub fn cone_terms(&self, roots: &[(usize, bool)]) -> Tape {
        let mut seen: std::collections::HashSet<u32> = roots.iter().map(|&(r, _)| r as u32).collect();
        let mut stack: Vec<u32> = seen.iter().copied().collect();
        while let Some(j) = stack.pop() {
            let (args, k) = self.ssa[j as usize].args();
            for &a in &args[..k] {
                if seen.insert(a) {
                    stack.push(a);
                }
            }
        }
        let mut ids: Vec<u32> = seen.into_iter().collect();
        ids.sort_unstable();
        let remap: std::collections::HashMap<u32, u32> = ids.iter().enumerate().map(|(n, &j)| (j, n as u32)).collect();
        let mut ssa: Vec<Ins> = ids
            .iter()
            .enumerate()
            .map(|(n, &j)| {
                let mut ins = self.ssa[j as usize];
                let (_, k) = ins.args();
                if k > 0 {
                    ins.a = remap[&ins.a];
                }
                if k > 1 {
                    ins.b = remap[&ins.b];
                }
                if k > 2 {
                    ins.c = remap[&ins.c];
                }
                ins.out = n as u32;
                ins
            })
            .collect();
        // Each root negated where asked, and after the first joined by a
        // max, the last the output (one root as it is: its own, ids' last,
        // as it reads every other).
        let term = |ssa: &mut Vec<Ins>, (r, up): (usize, bool)| {
            let r = remap[&(r as u32)];
            if up {
                return r;
            }
            let n = ssa.len() as u32;
            ssa.push(Ins { op: Op::Un(UnOp::Neg), out: n, a: r, b: 0, c: 0 });
            n
        };
        let mut out = term(&mut ssa, roots[0]);
        for &r in &roots[1..] {
            let b = term(&mut ssa, r);
            let n = ssa.len() as u32;
            ssa.push(Ins { op: Op::Bin(BinOp::Max), out: n, a: out, b, c: 0 });
            out = n;
        }
        // As this tape has them (a simplified tape's are not recomputable);
        // a joining max or negation is not fixed (bounded as any value is:
        // soundly).
        let joined = ssa.len() - ids.len();
        let mut t = Tape::from_ssa(ssa);
        t.fixed = ids.iter().map(|&j| self.fixed[j as usize]).chain(std::iter::repeat_n(false, joined)).collect();
        t.pmask = ids.iter().map(|&j| self.pmask[j as usize]).chain(std::iter::repeat_n(0, joined)).collect();
        (t.pread, t.plist) = self.reread(ids.iter().enumerate().map(|(n, &j)| (n as u32, j as usize)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        if t.ssa.iter().any(|ins| matches!(ins.op, Op::Read | Op::CornerReach(..))) {
            t.forms = self.forms.clone();
        }
        t
    }

    /// The tape with each instruction `i` (an ssa index) of `values` the
    /// constant given (a floor's value over a region it is constant in).
    pub fn with_consts(&self, values: &[(usize, f64)]) -> Tape {
        if values.is_empty() {
            return self.clone();
        }
        let mut ssa = self.ssa.clone();
        for &(i, v) in values {
            ssa[i] = Ins { op: Op::Const(v), out: i as u32, a: 0, b: 0, c: 0 };
        }
        let mut t = Tape::from_ssa(ssa);
        (t.fixed, t.pmask) = (self.fixed.clone(), self.pmask.clone());
        for &(i, _) in values {
            (t.fixed[i], t.pmask[i]) = (true, 0);
        }
        (t.pread, t.plist) = self.reread((0..self.ssa.len()).filter(|i| !values.iter().any(|v| v.0 == *i)).map(|i| (i as u32, i)));
        t.grids = self.grids.clone();
        t.opaque = self.opaque.clone();
        if t.ssa.iter().any(|ins| matches!(ins.op, Op::Read | Op::CornerReach(..))) {
            t.forms = self.forms.clone();
        }
        t
    }

    /// `pread` and `plist` of a tape whose instruction `n` is this one's
    /// `old` (`kept`, by n), each list copied once.
    fn reread(&self, kept: impl Iterator<Item = (u32, usize)>) -> (Vec<(u32, u32, u32)>, Vec<u32>) {
        let (mut pread, mut plist, mut at) = (Vec::new(), Vec::new(), std::collections::HashMap::new());
        for (n, old) in kept {
            if let Ok(j) = self.pread.binary_search_by_key(&(old as u32), |e| e.0) {
                let (_, start, len) = self.pread[j];
                let to = *at.entry(start).or_insert_with(|| {
                    plist.extend_from_slice(&self.plist[start as usize..(start + len) as usize]);
                    (plist.len() - len as usize) as u32
                });
                pread.push((n, to, len));
            }
        }
        (pread, plist)
    }

    /// The machine instructions simplifying it over a box into `simple`
    /// retires beside the box's pass (the choices it reads, and
    /// `simplify`): 200 an instruction it reads and 600 one it keeps,
    /// counted as `Op::work` is.
    pub fn simplify_work(&self, simple: &Tape) -> usize {
        200 * self.len() + 600 * simple.len()
    }
}
