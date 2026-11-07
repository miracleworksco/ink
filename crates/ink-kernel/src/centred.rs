//! The clearance proofs' first-order bounds: a lower bound of a field over
//! a box from centred forms of the terms of its max, alone and in pairs,
//! where plain interval bounds lose that two facing surfaces' slopes cancel.
//! Two faces across a thin gap are proven with cells about sqrt(gap x
//! radius of curvature) across, and two flat faces at any box size. The
//! clearance, emptiness and overlap searches (query) are their users.

use ink_core::eval::{Ctx, Inputs, MAX_INPUTS, eval_interval_grad_every, eval_interval_scope_every};
use ink_core::graph::{BinOp, UnOp};
use ink_core::interval::{self, G, I, add_o, add_up, mag, mul_o, mul_up, norm, sub_o, sum_up};
use ink_core::tape::{Op, Tape};

/// A lower bound of the tape over the box from first-order centred forms of
/// the terms of its max: the best of each term alone
/// (term_lower) and of each pair (max_lower_centred_paired),
/// and exactly of a term and its own negation (Same::opposed: a part
/// outside a zone that is its own feature, at any box size).
/// See `max_terms` for what counts as a term. None when the max has too
/// many terms to try. Over the scope's inputs (a time span, ranged
/// params: Inputs::of), the terms at the centre are bounded first order in
/// them, and a pair also from both terms at the inputs' middles and the
/// pair's rates along them (two faces that one param, or a motion, moves
/// together: their rates cancel as their slopes do).
pub fn centred_lower(tape: &Tape, ctx: &Ctx, b: [I; 3]) -> Option<f64> {
    centred_lower_max(&[(tape, 0.0)], ctx, b)
}

/// `centred_lower` of max_k (f_k - s_k) for fields f_k (tapes) and amounts
/// s_k: the terms of every f_k's max, each lowered by its s_k, are terms of
/// the whole, so a pair can take one term from each (two surfaces facing
/// each other across a gap: a clearance's shape and the other grown by it).
/// The box may be unbounded: the centre is then a finite point of it, and an
/// axis along which a term (or a pair's combination) has an exactly zero
/// slope adds nothing however far the box reaches, so the max of parallel
/// flat faces (two slabs) is bounded over all of space.
pub fn centred_lower_max(parts: &[(&Tape, f64)], ctx: &Ctx, b: [I; 3]) -> Option<f64> {
    centred_bound(parts, ctx, b, None, true, 0.0).map(|(lo, _)| lo)
}

/// Whether `centred_lower_max` shows max_k (f_k - s_k) >= `level` over the
/// box: it stops at the first bound that does, and tries none where the
/// field at the box's centre is below `level` (no lower bound over the box
/// can pass a value in it: a cell inside the shape, as an overlap's search
/// visits everywhere).
pub fn centred_at_least(parts: &[(&Tape, f64)], ctx: &Ctx, b: [I; 3], level: f64) -> bool {
    centred_bound(parts, ctx, b, Some(level), true, 0.0).is_some_and(|(lo, _)| lo >= level)
}

/// How the first-order forms settled a box (`centred_settles`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Settles {
    No,
    /// The tape is at least the level all over it (centred_at_least).
    AtLeast,
    /// Where it is below the level, every point lies within the tolerance
    /// of a point where it is not.
    Covered,
}

/// Whether the tape is at least 0 over the box (`centred_at_least`, tried
/// only `bounded`), or else its solid (where it is below 0) in the box
/// holds no point deeper than `tol`, by a pair of terms of its max that
/// face each other along an axis (two faces that touch, or overlap by less
/// than the tolerance, over any area). Their first-order bound (exact, of a
/// term and its own negation: Same) is l a + (1 - l) b >= L over the box,
/// so at a point of the solid, where both are below 0, a > L / l and b > L
/// / (1 - l); where a rises along the axis at least m_a and b falls at
/// least m_b, a is 0 within -L / (l m_a) one way and b within -L / ((1 -
/// l) m_b) the other, and the max, at least each, is not below 0 there.
/// Where both are within the tolerance and half the box's side, one way
/// stays in the box (where the slopes hold) from every point of it. Without
/// the choices (a min a term), whose ways would each exit elsewhere.
pub fn centred_settles(tape: &Tape, ctx: &Ctx, b: [I; 3], tol: f64, bounded: bool) -> Settles {
    match centred_bound(&[(tape, 0.0)], ctx, b, Some(0.0), bounded, tol) {
        // Its least bound not inside: no point is.
        Some((lo, _)) if !ink_core::sign::inside(lo) => Settles::AtLeast,
        Some((_, true)) => Settles::Covered,
        _ => Settles::No,
    }
}

/// Terms a first-order bound of a max tries, each alone and each pair (in
/// the largest of the ways through its choices).
const MAX_TERMS: usize = 16;

/// Ways through a max's choices a first-order bound tries.
const MAX_CHOICES: usize = 4;

/// The ways of taking one member of each choice, given how many members
/// each has: if at least one and at most MAX_CHOICES. Counted without
/// wrapping (16 choices of 16 make 2^64 ways, which wraps to none, and an
/// enumeration of none checks nothing).
fn ways(members: impl IntoIterator<Item = usize>) -> Option<usize> {
    members.into_iter().try_fold(1usize, usize::checked_mul).filter(|w| (1..=MAX_CHOICES).contains(w))
}

/// `centred_lower_max` (tried where `bound`), stopping once the bound
/// reaches `enough`; and with a tolerance `tol` above 0, whether a facing
/// pair covers the box below `enough` to it (centred_settles).
fn centred_bound(parts: &[(&Tape, f64)], ctx: &Ctx, b: [I; 3], enough: Option<f64>, bound: bool, tol: f64) -> Option<(f64, bool)> {
    if parts.iter().any(|&(tape, _)| tape.ssa.is_empty()) {
        return None;
    }
    // Seen through offsets, negations and choices a max has no fewer terms
    // than without (each term becomes one or more): too many without, and
    // nothing need be evaluated.
    let plain = parts.iter().map(|&(tape, _)| max_terms(tape, &[], false, false, false, I::point(0.0)).map(|(t, ..)| t.len())).sum::<Option<usize>>();
    if plain.is_none_or(|n| n > MAX_TERMS) {
        return None;
    }
    // A centre inside the box (halved first, so it cannot overflow; on an
    // unbounded axis, a finite end, else 0), and the distance from it to the
    // box's far side on each axis, rounded up: the centre may round to one
    // end of a narrow side.
    let mid = b.map(|x| {
        let m = (0.5 * x.lo + 0.5 * x.hi).max(x.lo).min(x.hi);
        if m.is_finite() {
            m
        } else if x.lo.is_finite() {
            x.lo
        } else if x.hi.is_finite() {
            x.hi
        } else {
            0.0
        }
    });
    let hw: [f64; 3] = std::array::from_fn(|k| (mid[k] - b[k].lo).max(b[k].hi - mid[k]).next_up());
    let c = mid.map(I::point);
    let at_cs: Vec<Vec<interval::G>> = parts.iter().map(|&(tape, _)| eval_interval_grad_every(tape, ctx, c)).collect();
    // The field at the centre (over the span: at every moment of it) is at
    // most this, and no lower bound over the box passes it. A cover needs
    // the centre within about the tolerance of the surface (to first order
    // by the slope there: else none is tried, for speed alone).
    let centre = parts.iter().zip(&at_cs).map(|(&(_, s), g)| (g[g.len() - 1].v.hi - s).next_up()).fold(f64::NEG_INFINITY, f64::max);
    let bound = bound && enough.is_none_or(|level| centre >= level);
    let cover = tol > 0.0
        && enough.is_some_and(|level| {
            let slope = at_cs.iter().map(|g| interval::norm(g[g.len() - 1].d.map(interval::mag))).fold(0.0, f64::max);
            !(level - centre > 2.0 * tol * slope)
        });
    if !bound && !cover {
        return Some((f64::NEG_INFINITY, false));
    }
    // The inputs the terms are first order in (one list for every part, so
    // a pair can take a term of each), and how far each reaches from its
    // middle.
    let inputs = Inputs::of(parts.iter().map(|&(tape, _)| tape), ctx);
    let hu = inputs.half_widths(ctx);
    let mut evals = Vec::with_capacity(parts.len());
    for (&(tape, s), mut at_c) in parts.iter().zip(at_cs) {
        let over = eval_interval_grad_every(tape, ctx, b);
        // The values at the centre first order in the inputs (only the
        // terms' values at the centre are used below), and each one's value
        // at the inputs' middles and rates along them, for the pairs.
        let scoped = (!inputs.is_empty()).then(|| eval_interval_scope_every(tape, ctx, &inputs, c, true));
        if let Some(sc) = &scoped {
            for (g, v) in at_c.iter_mut().zip(&sc.v) {
                g.v = *v;
            }
        }
        evals.push((tape, s, at_c, over, scoped));
    }
    // The terms, seen through shifts, negations and negated mins, with the
    // choices a min among them makes, where they fit; else as the terms
    // were seen before those (no choices).
    let n = inputs.len();
    let as_g = |t: &[Term], k: usize| -> Vec<TermBounds> {
        let (_, s, at_c, over, scoped) = &evals[k];
        t.iter()
            .map(|tm| {
                let signed = |v: I, neg: bool| if neg { I::new(-v.hi, -v.lo) } else { v };
                let term = |g: &interval::G| {
                    let (v, d) = if tm.neg { (signed(g.v, true), g.d.map(|x| signed(x, true))) } else { (g.v, g.d) };
                    interval::G { v: add_o(v, tm.shift), d, l: g.l, s2: None, m: interval::NO_GRAM }
                };
                // At the inputs' middles, and its rates: the term's, the
                // part's amount, and each value its shift adds.
                let inputs = scoped.as_ref().filter(|sc| tm.by.is_some() && !sc.d.is_empty()).map(|sc| {
                    let (mut m, mut r) = (add_o(signed(sc.m[tm.i], tm.neg), I::point(-s)), [I::point(0.0); MAX_INPUTS]);
                    for (j, r) in r.iter_mut().enumerate().take(n) {
                        *r = signed(sc.rate(tm.i, n, j), tm.neg);
                    }
                    for &(k, neg) in tm.by.as_ref().map_or(&[][..], |b| b.as_slice()) {
                        m = add_o(m, signed(sc.m[k as usize], neg));
                        for (j, r) in r.iter_mut().enumerate().take(n) {
                            *r = add_o(*r, signed(sc.rate(k as usize, n, j), neg));
                        }
                    }
                    Box::new((m, r))
                });
                TermBounds { at_c: term(&at_c[tm.i]), over: term(&over[tm.i]), inputs, same: Same::of(k, tm) }
            })
            .collect()
    };
    // Those that fit: no more than MAX_TERMS terms in the largest of the
    // ways through the choices, and no more than MAX_CHOICES ways.
    let collect = |flat: bool, choose: bool| {
        let (mut terms, mut choices, mut split) = (Vec::new(), Vec::new(), Vec::new());
        for (k, (tape, s, _, over, _)) in evals.iter().enumerate() {
            let (t, c, sp) = max_terms(tape, over, flat, choose, !inputs.is_empty(), I::point(-s))?;
            terms.extend(as_g(&t, k));
            choices.extend(c.into_iter().map(|alts| alts.into_iter().map(|a| as_g(&a, k)).collect::<Vec<_>>()));
            split.extend(sp.iter().map(|tm| Same::of(k, tm)));
        }
        let largest = terms.len() + choices.iter().map(|a: &Vec<Vec<_>>| a.iter().map(Vec::len).max().unwrap_or(0)).sum::<usize>();
        (largest <= MAX_TERMS && ways(choices.iter().map(Vec::len)).is_some()).then_some((terms, choices, split))
    };
    let (terms, choices, split) = collect(true, true).or_else(|| collect(true, false)).or_else(|| collect(false, false))?;
    if terms.is_empty() && choices.is_empty() {
        return None;
    }
    let level = enough.unwrap_or(0.0);
    // The best bound of the max of these terms: each alone, and each pair
    // (or the first that is enough).
    let enough = enough.unwrap_or(f64::INFINITY);
    let pairs = |t: &[TermBounds]| {
        let mut best = f64::NEG_INFINITY;
        for (k, a) in t.iter().enumerate() {
            best = best.max(term_lower(a.at_c, a.over, hw));
            for x in &split {
                best = best.max(a.same.opposed(x));
            }
            for b in &t[k + 1..] {
                if best >= enough {
                    return best;
                }
                best = best.max(a.same.opposed(&b.same));
                best = best.max(max_lower_centred_paired(a.at_c, b.at_c, a.over, b.over, hw));
                if let (Some(x), Some(y)) = (&a.inputs, &b.inputs) {
                    best = best.max(max_lower_centred_inputs((x.0, &x.1[..n]), (y.0, &y.1[..n]), &a.over, &b.over, hw, &hu[..n]));
                }
            }
        }
        best
    };
    // max(T, min(u, v)) = min(max(T, u), max(T, v)): the least over every
    // way of taking one member of each choice (once one is not enough, no
    // other can make up for it). There is at least one way (with no
    // choices, one: the terms), so the bound is never the starting +inf.
    let mut worst = if bound { f64::INFINITY } else { f64::NEG_INFINITY };
    for w in (0..ways(choices.iter().map(Vec::len))?).filter(|_| bound) {
        let mut t = terms.clone();
        let mut rest = w;
        for c in &choices {
            t.extend_from_slice(&c[rest % c.len()]);
            rest /= c.len();
        }
        worst = worst.min(pairs(&t));
        if worst < enough && enough.is_finite() {
            break;
        }
    }
    if worst >= enough || !cover {
        return Some((worst, false));
    }
    // A facing pair among the terms without choices (a min a term).
    let plain = match choices.is_empty() {
        true => terms,
        false => collect(true, false).or_else(|| collect(false, false)).map_or_else(Vec::new, |(t, ..)| t),
    };
    let half: [f64; 3] = b.map(|x| (0.5 * (x.hi - x.lo)).next_down());
    let covers = plain.iter().enumerate().any(|(k, a)| plain[k + 1..].iter().any(|c| facing_cover(a, c, hw, half, level, tol)));
    Some((worst, covers))
}

/// Whether terms `a` and `c` cover the box's solid below `level` to `tol`
/// along an axis (centred_settles): their pair's bound L (below the level)
/// and weight l, and along some axis one rising at least m and the other
/// falling at least m', each reaching the level from any point of the
/// solid within (level - L) over its weight times its slope, no more than
/// the tolerance and `half` the box's side there.
fn facing_cover(a: &TermBounds, c: &TermBounds, hw: [f64; 3], half: [f64; 3], level: f64, tol: f64) -> bool {
    let exact = a.same.opposed(&c.same);
    let (lower, l) = match exact > f64::NEG_INFINITY {
        true => (exact, 0.5),
        false => max_lower_centred_paired_weighted(a.at_c, c.at_c, a.over, c.over, hw),
    };
    let gap = (level - lower).next_up();
    if !(gap > 0.0 && gap.is_finite() && l > 0.0 && l < 1.0) {
        return false;
    }
    // Rounded so that each reach is at least the exact one; a product that
    // rounds to 0 or below (subnormal) bounds no distance, so reaches none.
    let reach = |w: f64, slope: f64| {
        let d = (w * slope).next_down().next_down();
        if d > 0.0 { (gap / d).next_up() } else { f64::INFINITY }
    };
    (0..3).any(|k| {
        [(a, c, l), (c, a, 1.0 - l)].into_iter().any(|(up, down, w)| {
            let (rise, fall) = (up.over.d[k].lo, -down.over.d[k].hi);
            rise > 0.0 && fall > 0.0 && {
                let most = tol.min(half[k]);
                reach(w, rise) <= most && reach(1.0 - w, fall) <= most
            }
        })
    })
}

/// A term's bounds for the first-order forms: at the box's centre (over the
/// whole scope beyond space, first order in its inputs where it has any)
/// and over the box (spatial slopes over box and scope); and where the
/// inputs are taken (centred_bound), its value at the centre of both and
/// its rates along the inputs over the scope at the box's centre.
#[derive(Clone)]
struct TermBounds {
    at_c: interval::G,
    over: interval::G,
    inputs: Option<Box<(I, [I; MAX_INPUTS])>>,
    same: Same,
}

/// Which value a term (or a max split into its terms) is: its part, its
/// instruction and whether it is negated, with its shift.
#[derive(Clone, Copy)]
struct Same {
    part: usize,
    i: usize,
    neg: bool,
    shift: I,
}

impl Same {
    fn of(part: usize, tm: &Term) -> Same {
        Same { part, i: tm.i, neg: tm.neg, shift: tm.shift }
    }

    /// A lower bound of the max of this term and another: where the other
    /// is this one's own value negated (a zone that is a part's own feature,
    /// overlap outside it asked as max(part, -zone)), their mean is the mean
    /// of their shifts at every point, however the value varies over the box
    /// (bounds of the two taken apart lose that they cancel); else none.
    fn opposed(&self, o: &Same) -> f64 {
        if !(self.part == o.part && self.i == o.i && self.neg != o.neg) {
            return f64::NEG_INFINITY;
        }
        // Halving is exact but for a subnormal sum.
        let (sum, half) = (add_o(self.shift, o.shift).lo, 0.5 * add_o(self.shift, o.shift).lo);
        if 2.0 * half == sum { half } else { half.next_down() }
    }
}

/// Values a term's shift may add up (its `by`) and still be bounded first
/// order in the inputs.
const SHIFTS: usize = 4;

/// A term of a max: an instruction, whether it is negated, and a shift added
/// to it (its bound), with the values the shift adds up beyond the part's
/// own amount (each with whether it is subtracted: a shift that reads a
/// ranged param moves with it, which its bound alone forgets); None where
/// there are more than SHIFTS of them.
#[derive(Clone, Copy)]
struct Term {
    i: usize,
    neg: bool,
    shift: I,
    by: Option<Shifts>,
}

/// Up to SHIFTS values a shift adds up, each with whether it is subtracted.
#[derive(Clone, Copy)]
struct Shifts {
    at: [(u32, bool); SHIFTS],
    n: u8,
}

impl Shifts {
    const NONE: Shifts = Shifts { at: [(0, false); SHIFTS], n: 0 };

    /// These and one more, if they fit.
    fn with(mut self, x: (u32, bool)) -> Option<Shifts> {
        let n = self.n as usize;
        (n < SHIFTS).then(|| {
            self.at[n] = x;
            self.n += 1;
            self
        })
    }

    fn as_slice(&self) -> &[(u32, bool)] {
        &self.at[..self.n as usize]
    }
}

impl Term {
    /// Instruction `i` (negated with `neg`) with the shift `shift`, which
    /// adds the value `add` (and its sign) to this one's.
    fn then(&self, i: usize, neg: bool, shift: I, add: Option<(u32, bool)>) -> Term {
        let by = match add {
            None => self.by,
            Some(x) => self.by.and_then(|b| b.with(x)),
        };
        Term { i, neg, shift, by }
    }
}

/// The terms of the max a tape's root takes, for first-order lower bounds:
/// each an instruction, whether it is negated, and a shift added to it, so
/// the root is at least the max of (+-value + shift) over them. Through
/// nested max and emax (which is at least the max of its terms), tags, and
/// promises: a promised term (lipschitz) is a term itself, so its promise
/// bounds it whole even where its formula's own slope is unbounded (a
/// thread's atan2 cut); a promised max is also split into its terms. With
/// `flat`, also through a shift by a value constant in space over the box
/// (`over`: its slope is exactly 0; a constant, a param, a time), which
/// shifts each term (max(a, b) - c = max(a - c, b - c): a clearance's
/// offset), and a negation, which turns a min into a max of negated terms
/// (-min(a, b) = max(-a, -b): a union subtracted). With `choose`, the
/// choices too: a min among the terms (under a negation, a max: an
/// intersection subtracted, as a pocket's floor and wall) is at least the
/// least of its members, so the root is at least the least, over its
/// members, of the max with that member's terms in its place; each choice
/// lists its members' terms (a min within one is a term). Without, a min is
/// a term. `shift` is the root's own. None where there are more than
/// MAX_TERMS terms (or in a member of a choice), a choice has more than
/// MAX_CHOICES members, or the choices more than MAX_CHOICES ways: the walk
/// stops there, so a max whose terms share subexpressions (walked once for
/// each way to them: max(f, f + 1) nested doubles them each time) costs no
/// more than its first MAX_TERMS terms. With `by`, each term also records
/// the values its shift adds up (Term::by: over ranged params). Also the
/// maxes it split (outside the choices): the root is at least each of
/// those too, which a term negating one's own value cancels (Same).
fn max_terms(tape: &Tape, over: &[interval::G], flat: bool, choose: bool, by: bool, shift: I) -> Option<(Vec<Term>, Vec<Vec<Vec<Term>>>, Vec<Term>)> {
    let Some(root) = tape.ssa.len().checked_sub(1) else { return Some((Vec::new(), Vec::new(), Vec::new())) };
    let (mut choices, mut split) = (Vec::new(), Vec::new());
    let start = Term { i: root, neg: false, shift, by: by.then_some(Shifts::NONE) };
    let terms = terms_from(tape, over, flat, start, choose.then_some(&mut choices), &mut split)?;
    Some((terms, choices, split))
}

/// `max_terms` from one start; with `choices`, a min among the terms adds a
/// choice there, else it is a term.
fn terms_from(tape: &Tape, over: &[interval::G], flat: bool, start: Term, mut choices: Option<&mut Vec<Vec<Vec<Term>>>>, split: &mut Vec<Term>) -> Option<Vec<Term>> {
    let ssa = &tape.ssa;
    let untag = |mut j: usize| {
        while let Op::Tag(_) = ssa[j].op {
            j = ssa[j].a as usize;
        }
        j
    };
    let level = |j: usize| over[j].l == 0.0 && over[j].d.iter().all(|x| x.lo == 0.0 && x.hi == 0.0);
    let signed = |v: I, neg: bool| if neg { I::new(-v.hi, -v.lo) } else { v };
    // A max of terms, under a negation a min (its negation is a max).
    let splits = |op: Op, neg: bool| match op {
        Op::Bin(BinOp::Max | BinOp::Emax) => !neg,
        Op::Bin(BinOp::Min) => neg && flat,
        _ => false,
    };
    // A min of terms, under a negation a max (not emax: it is at least the
    // max, which bounds its negation the wrong way).
    let chooses = |op: Op, neg: bool| {
        flat && match op {
            Op::Bin(BinOp::Min) => !neg,
            Op::Bin(BinOp::Max) => neg,
            _ => false,
        }
    };
    let mut terms = Vec::new();
    let mut stack = vec![start];
    while let Some(tm) = stack.pop() {
        let (i, neg, s) = (tm.i, tm.neg, tm.shift);
        let ins = &ssa[i];
        let (a, b) = (ins.a as usize, ins.b as usize);
        match ins.op {
            Op::Tag(_) => stack.push(tm.then(a, neg, s, None)),
            op if splits(op, neg) => {
                split.push(tm);
                stack.extend([tm.then(a, neg, s, None), tm.then(b, neg, s, None)]);
            }
            Op::Un(UnOp::Neg) if flat => stack.push(tm.then(a, !neg, s, None)),
            // +-(a + k) + s = +-a + (s +- k), and so on.
            Op::Bin(BinOp::Add) if flat && level(b) => stack.push(tm.then(a, neg, add_o(s, signed(over[b].v, neg)), Some((b as u32, neg)))),
            Op::Bin(BinOp::Add) if flat && level(a) => stack.push(tm.then(b, neg, add_o(s, signed(over[a].v, neg)), Some((a as u32, neg)))),
            Op::Bin(BinOp::Sub) if flat && level(b) => stack.push(tm.then(a, neg, sub_o(s, signed(over[b].v, neg)), Some((b as u32, !neg)))),
            Op::Bin(BinOp::Sub) if flat && level(a) => stack.push(tm.then(b, !neg, add_o(s, signed(over[a].v, neg)), Some((a as u32, neg)))),
            Op::Lipschitz => {
                terms.push(tm);
                // Its formula as a term too would add nothing: the same
                // value and partials, without the promise.
                let a = untag(a);
                if matches!(ssa[a].op, Op::Lipschitz) || splits(ssa[a].op, neg) {
                    stack.push(tm.then(a, neg, s, None));
                }
            }
            op if chooses(op, neg) && choices.is_some() => {
                // Its members, through nested mins (maxes) and tags.
                let (mut members, mut more) = (Vec::new(), vec![a, b]);
                while let Some(j) = more.pop() {
                    let j = untag(j);
                    if matches!((ssa[j].op, op), (Op::Bin(x), Op::Bin(y)) if x == y) {
                        more.extend([ssa[j].a as usize, ssa[j].b as usize]);
                    } else {
                        members.push(j);
                        if members.len() > MAX_CHOICES {
                            return None;
                        }
                    }
                }
                let alts = members.into_iter().map(|m| terms_from(tape, over, flat, tm.then(m, neg, s, None), None, &mut Vec::new())).collect::<Option<_>>()?;
                let choices = choices.as_deref_mut().expect("checked");
                choices.push(alts);
                ways(choices.iter().map(Vec::len))?;
            }
            _ => terms.push(tm),
        }
        if terms.len() > MAX_TERMS {
            return None;
        }
    }
    Some(terms)
}

/// A lower bound of max(a, b) over a box within `hw` of a centre on each axis,
/// from a's and b's values at the centre (`ac`, `bc`) and their gradient
/// bounds over the box (`a`, `b`). max(a, b) >= l a + (1 - l) b for any l in
/// [0, 1], and the mean value theorem bounds that combination to first order.
/// Where two surfaces face each other across a thin gap, their gradients
/// nearly cancel in the combination, so this proves the gap with cells of
/// about sqrt(gap x radius of curvature), where plain bounds need cells as
/// small as the gap itself. Two flat faces' exactly opposite slopes cancel
/// exactly (l = 1/2), so their gap is proven over any box, even an
/// unbounded one. The bound is taken at the one l that best cancels the
/// gradients (where each term's own bound is taken apart: l = 0 and 1).
pub fn max_lower_centred_paired(ac: G, bc: G, a: G, b: G, hw: [f64; 3]) -> f64 {
    max_lower_centred_paired_weighted(ac, bc, a, b, hw).0
}

/// `max_lower_centred_paired`, and the weight l its bound of l a + (1 - l)
/// b is taken at (1 - l exact: the weights pair_lower takes).
pub fn max_lower_centred_paired_weighted(ac: G, bc: G, a: G, b: G, hw: [f64; 3]) -> (f64, f64) {
    let l = pair_weight(&a, &b);
    let v = pair_lower(ac, bc, a, b, hw, l);
    (if v.is_nan() { f64::NEG_INFINITY } else { v }, 1.0 - (1.0 - l))
}

/// `max_lower_centred_paired` over a box and a scope of inputs beyond space
/// (time, ranged params), centred in both: from a's and b's values at the
/// centre of both (`am.0`, `bm.0`), their rates along each input over the
/// scope at the box's centre (`am.1`, `bm.1`), their spatial partials over
/// the box and scope (`a`, `b`), and how far each input reaches from its
/// centre (`hu`). From the centre (c, u_c) to a point (x, u), first along
/// the inputs at c, then in space at u: the combination changes by at most
/// its rates' bound times |u - u_c| on the first leg and its slopes' bound
/// times |x - c| on the second (each leg inside the box and scope, where
/// those bounds hold). Two faces that one param moves together (two parts
/// touching where the param places both) have opposite rates as well as
/// opposite slopes, and both cancel; bounded at the centre over the whole
/// scope, as max_lower_centred_paired is, each face is blurred by its whole
/// motion. A Lipschitz promise bounds the slopes only (they are a's and
/// b's); the rates are their formulas'.
pub fn max_lower_centred_inputs(am: (I, &[I]), bm: (I, &[I]), a: &G, b: &G, hw: [f64; 3], hu: &[f64]) -> f64 {
    let l = pair_weight(a, b);
    // Weights that sum to exactly 1: 1 - (1 - l) is exact (Sterbenz).
    let m = 1.0 - l;
    let (li, mi) = (I::point(1.0 - m), I::point(m));
    let comb = |x: I, y: I| add_o(mul_o(li, x), mul_o(mi, y));
    let space = sum_up((0..3).map(|k| spread(mag(comb(a.d[k], b.d[k])), hw[k])));
    let inputs = sum_up((0..hu.len()).map(|j| spread(mag(comb(am.1[j], bm.1[j])), hu[j])));
    let v = (comb(am.0, bm.0).lo - add_up(space, inputs)).next_down();
    if v.is_nan() { f64::NEG_INFINITY } else { v }
}

/// The weight l of l a + (1 - l) b that best cancels the gradients (their
/// midpoints) at the centre. Any l in [0, 1] gives a valid bound; this one
/// is scaled so that it stays finite for slopes near either end of f64's
/// range (halved before the sum, and divided by the largest: the gradients
/// 1e308 and -1e308 give 1/2, where their sum's square would overflow to
/// a NaN weight and no bound).
fn pair_weight(a: &G, b: &G) -> f64 {
    let mid = |g: &G| g.d.map(|x| x.mid());
    let (ga, gb) = (mid(a), mid(b));
    let s = ga.iter().chain(&gb).fold(0.0f64, |m, v| m.max(v.abs()));
    if !(s > 0.0 && s.is_finite()) {
        return 0.5;
    }
    let (ga, gb) = (ga.map(|v| v / s), gb.map(|v| v / s));
    let diff: Vec<f64> = (0..3).map(|k| ga[k] - gb[k]).collect();
    let dd: f64 = diff.iter().map(|x| x * x).sum();
    if dd > 0.0 { (-(0..3).map(|k| diff[k] * gb[k]).sum::<f64>() / dd).clamp(0.0, 1.0) } else { 0.5 }
}

/// The first-order lower bound of l a + (1 - l) b over the box (NaN where
/// it is not defined).
fn pair_lower(ac: G, bc: G, a: G, b: G, hw: [f64; 3], l: f64) -> f64 {
    // Weights that sum to exactly 1: 1 - (1 - l) is exact (Sterbenz).
    let m = 1.0 - l;
    let (li, mi) = (I::point(1.0 - m), I::point(m));
    let comb = |x: I, y: I| add_o(mul_o(li, x), mul_o(mi, y));
    let spread = sum_up((0..3).map(|k| spread(mag(comb(a.d[k], b.d[k])), hw[k])));
    (comb(ac.v, bc.v).lo - spread).next_down()
}

/// How much a term whose slope along an axis is at most `slope` can change
/// within `half` of the centre along it: their product rounded up, and 0
/// for an exactly zero slope however far the box reaches.
fn spread(slope: f64, half: f64) -> f64 {
    if slope == 0.0 { 0.0 } else { mul_up(slope, half) }
}

/// A lower bound of a over a box within `hw` of a centre on each axis, from its
/// value there and its gradient bounds over the box (mean value theorem).
/// Tight where a's formula has dependencies that inflate interval values (a
/// branch cut, a rotation) but a finite gradient bound.
pub fn term_lower(ac: G, a: G, hw: [f64; 3]) -> f64 {
    let by_partials = sum_up((0..3).map(|k| spread(mag(a.d[k]), hw[k])));
    let by_scalar = spread(a.l, norm(hw));
    let v = (ac.v.lo - by_partials.min(by_scalar)).next_down();
    if v.is_nan() { f64::NEG_INFINITY } else { v }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_core::eval::ParamSpans;
    use ink_core::graph::{Graph, Var, bin_f64, un_f64};
    use ink_core::interval::NO_GRAM;

    /// eval_interval_scope_every (scope_run over the SSA form) keeps every value for its caller (the
    /// centred bound reads the terms of a max): where a jump is a choice,
    /// every value the fork sets bounds both sides, the jump's own and its
    /// window's, not only the one it joins at. atan2(p, -1) is near pi for
    /// p >= 0 and near -pi below, floor(p) 0 or -1: the last side alone
    /// holds one of them. In max(a, sin(a)) the jump itself is a max term
    /// (the centred bound read 3.04 against the field's -0.0995 at p =
    /// -0.1, and -5e-324 against -0.84 for the floor); in max(sin(atan2(p,
    /// -1) / 2), x) the half angle is inside the window, which joins at the
    /// sine. Each value's bound must hold it at every p of the span, and
    /// the centred lower bound must not pass the field.
    #[test]
    fn every_kept_value_of_a_fork_bounds_both_sides() {
        let mut cases = Vec::new();
        for jump in [0, 1] {
            let mut g = Graph::new();
            let (p, m1) = (g.var(Var::Param(0)), g.constant(-1.0));
            let a = if jump == 0 { g.bin(BinOp::Atan2, p, m1) } else { g.un(UnOp::Floor, p) };
            let s = g.un(UnOp::Sin, a);
            let f = g.bin(BinOp::Max, a, s);
            cases.push(Tape::new(&g, f));
        }
        let mut g = Graph::new();
        let (p, m1, half, x) = (g.var(Var::Param(0)), g.constant(-1.0), g.constant(0.5), g.var(Var::X));
        let a = g.bin(BinOp::Atan2, p, m1);
        let h = g.bin(BinOp::Mul, half, a);
        let s = g.un(UnOp::Sin, h);
        let f = g.bin(BinOp::Max, s, x);
        cases.push(Tape::new(&g, f));
        let ctx = Ctx { p_spans: ParamSpans::new([(0, (-0.1, 0.1))]), ..Ctx::new(0.0, vec![0.0]) };
        let inputs = Inputs { t: false, params: vec![0] };
        let b = [I::point(-5.0); 3];
        for tape in &cases {
            let sc = eval_interval_scope_every(tape, &ctx, &inputs, b, false);
            let lower = centred_lower(tape, &ctx, b).unwrap_or(f64::NEG_INFINITY);
            for k in 0..=20 {
                let p = -0.1 + 0.01 * k as f64;
                // Every SSA value at p, in floats.
                let mut v = vec![0.0f64; tape.ssa.len()];
                for (n, ins) in tape.ssa.iter().enumerate() {
                    v[n] = match ins.op {
                        Op::Const(c) => c,
                        Op::Var(Var::Param(_)) => p,
                        Op::Var(_) => -5.0,
                        Op::Un(op) => un_f64(op, v[ins.a as usize]),
                        Op::Bin(op) => bin_f64(op, v[ins.a as usize], v[ins.b as usize]),
                        _ => unreachable!(),
                    };
                }
                for (n, &x) in v.iter().enumerate() {
                    assert!(sc.v[n].contains(x), "value {n} ({:?}) at p {p} = {x} outside {:?}", tape.ssa[n].op, sc.v[n]);
                }
                assert!(lower <= v[v.len() - 1], "centred lower {lower} above the field {} at p {p}", v[v.len() - 1]);
            }
        }
    }

    #[test]
    fn opposite_flat_faces_bound_their_max_over_unbounded_boxes() {
        // max(z - 1, 1.04 - z) >= 0.02 everywhere: the pair's slopes cancel
        // exactly at l = 1/2, and a zero slope adds nothing however far the
        // box reaches. Slopes that do not cancel exactly give no bound there.
        let z = |s: f64| G { v: I::ENTIRE, d: [I::point(0.0), I::point(0.0), I::point(s)], l: 1.0, s2: None, m: NO_GRAM };
        let at = |v: f64, g: G| G { v: I::point(v), ..g };
        let hw = [f64::INFINITY; 3];
        let (a, b) = (z(1.0), z(-1.0));
        // (1.04 - 1) / 2 of the floats: 0.02 and 18e-18, exactly.
        let lo = max_lower_centred_paired(at(-1.0, a), at(1.04, b), a, b, hw);
        assert!(lo > 0.0199 && lo <= 0.5 * 1.04 - 0.5, "{lo}");
        let b2 = z(-1.0 - f64::EPSILON);
        assert_eq!(max_lower_centred_paired(at(-1.0, a), at(1.04, b2), a, b2, hw), f64::NEG_INFINITY);
        let w = G { d: [I::point(0.0), I::point(0.0), I::new(0.9, 1.0)], ..a };
        assert_eq!(max_lower_centred_paired(at(-1.0, w), at(1.04, b), w, b, hw), f64::NEG_INFINITY);
        assert_eq!(term_lower(at(1.0, a), a, hw), f64::NEG_INFINITY);
        let flat = z(0.0);
        assert_eq!(term_lower(at(1.0, flat), G { l: 0.0, ..flat }, hw), 1.0f64.next_down());
    }

    #[test]
    fn pair_bounds_hold_their_strength_at_the_ends_of_the_range() {
        // max(1 + 1e308 z, 1 - 1e308 z) >= 1: the slopes cancel at l = 1/2.
        // The weight holds where midpoints would overflow (1e308 slopes) and
        // where squares underflow (subnormal slopes).
        let z = |s: f64| G { v: I::ENTIRE, d: [I::point(0.0), I::point(0.0), I::point(s)], l: s.abs(), s2: None, m: NO_GRAM };
        let at = |v: f64, g: G| G { v: I::point(v), ..g };
        for (s, t) in [(1e308, -1e308), (1.5e308, -1e308), (5e-324, -5e-324), (1e-300, -3e-300)] {
            let (a, b) = (z(s), z(t));
            assert!((0.0..=1.0).contains(&pair_weight(&a, &b)));
            let lo = max_lower_centred_paired(at(1.0, a), at(1.0, b), a, b, [1.0; 3]);
            let want = if s == -t { 1.0f64.next_down() } else { f64::NEG_INFINITY };
            assert!(lo >= want && lo < 1.0, "{s}, {t}: {lo}");
        }
        assert_eq!(pair_weight(&z(1e308), &z(-1e308)), 0.5);
        assert_eq!(I::new(1e308, 1.5e308).mid(), 1.25e308);
        assert_eq!(I::new(5e-324, 5e-324).mid(), 5e-324);
    }
}
