use super::*;

/// A value at a point (no ranged param) against its comparison: it holds,
/// or fails (an undefined value, NaN, too); a value within rounding of
/// the bound moved by the accuracy is unknown.
pub(super) fn point_status(c: Comparison, x: f64) -> Status {
    let p = I::point(x);
    if c.holds_on(p) {
        Status::Holds
    } else if c.fails_on(p) || x.is_nan() {
        Status::Fails
    } else {
        Status::Unknown
    }
}

/// `(subset a b)`'s search field, max(a, -b): what its witnesses read, and
/// its question with a and b themselves (`query::Except`: no point inside a
/// and not inside b, which counts where b reads 0).
pub(crate) fn subset_field(g: &mut ink_core::graph::Graph, a: u32, b: u32) -> u32 {
    let nb = g.neg(b);
    g.max(a, nb)
}

/// What a subset's search covers: a's reach, where its points lie (the
/// search field reads 0, not below it, where b does, so that field's own
/// reach may leave out a's material). Its proof and its range check read
/// this one.
pub(crate) fn subset_reach(a: u32, _b: u32) -> u32 {
    a
}

/// The solid of field `f` holds no point deeper than `tol` (an emptiness
/// claim, or an overlap's, a clearance of 0's); with `except` (a, b), f
/// being max(a, -b), no point inside a and not inside b (a subset's).
pub(super) fn shallow(v: &mut ModelView, f: u32, tol: f64, res: Option<f64>) -> Check {
    shallow_except(v, f, None, tol, res)
}

pub(super) fn shallow_except(v: &mut ModelView, f: u32, except: Option<(u32, u32)>, tol: f64, res: Option<f64>) -> Check {
    let (p, whole, tape) = v.prove_except(f, except, res, tol);
    let points = if p.exact { Vec::new() } else { Held::scoped(&p, &tape, &v.ctx) };
    let held = Some(Held::confirmed_at(p.exact, points, &tape, &v.ctx));
    Check { held, ..emptied(p, whole) }
}

/// Only exact coverage proves support. A point outside that coverage
/// starts a separate search for a certified layer of air.
pub(super) fn check_overhang(v: &mut ModelView, region: crate::support::Overhang, res: Option<f64>, tol: f64) -> Check {
    let mut options = v.proof_options(res, tol);
    options.tol = 0.0;
    let (mut proof, whole) = v.prove_with(region.possible, &options);
    if whole && proof.exact() {
        return Check { held: Some(Held::Outright), ..emptied(proof, true) };
    }
    if let Some(gap) = region.witness(v, &proof) {
        return checked(Status::Fails, gap);
    }
    let inner = shallow(v, region.certain, tol, res);
    if inner.status == Status::Fails {
        return inner;
    }
    match &mut proof.verdict {
        Verdict::Disproven { witness: [x, y, z], .. } => {
            checked(Status::Unknown, Evidence::Unproven(format!("a continuous support path is unproved near ({x}, {y}, {z}) mm; a gap there is unconfirmed")))
        }
        Verdict::Unknown { attained, .. } => {
            *attained = None;
            Check { note: Some("a continuous support path is unproved in the undecided cells"), ..emptied(proof, whole) }
        }
        _ => checked(Status::Unknown, Evidence::Unexamined),
    }
}

/// An emptiness claim's answer from its proof, all of its field examined
/// or not (`whole`: else it holds in the window only).
pub(crate) fn emptied(p: Proof, whole: bool) -> Check {
    in_window(Status::of(p.status()), whole, Evidence::Empty(p))
}

/// A touches claim's note where its parts come no nearer each other than
/// its accuracy.
pub(crate) const NOT_TOUCHING: &str = "they come no nearer each other than the accuracy: not touching";

/// `(touches a b)`: they overlap no deeper than the accuracy `tol`
/// (`shallow`), and come within it of each other at every value of the
/// ranged params (`contact`).
pub(super) fn touches(v: &mut ModelView, a: u32, b: u32, tol: f64, res: Option<f64>) -> Check {
    let f = v.doc.g.max(a, b);
    let over = shallow(v, f, tol, res);
    let overlap = Box::new(over.evidence.clone());
    if over.status != Status::Holds {
        return Check { evidence: Evidence::Touches { overlap, gap: Box::new(Evidence::Unexamined) }, ..over };
    }
    let gap = contact(v, a, b, tol, res);
    let mut c = Check { note: gap.note, ..checked(gap.status, Evidence::Touches { overlap, gap: Box::new(gap.evidence) }) };
    if c.status == Status::Fails {
        c.note = Some(NOT_TOUCHING);
    }
    if gap.unbounded {
        (c.unbounded, c.in_window) = (true, gap.in_window);
    }
    c
}

/// A point of b within d of one of a, at every value of the ranged params
/// a and b read: proven by a pair of points closer than d (the exact
/// clearance search's witness, which disproves a clearance of d), and
/// disproven where that clearance is proven. A pair found at one value of
/// ranged params says nothing of the others: there it is unknown.
pub(crate) fn contact(v: &mut ModelView, a: u32, b: u32, d: f64, res: Option<f64>) -> Check {
    let finest = if d > 0.0 { 0.5 * d } else { v.finest() };
    let read: Vec<u32> = ranged_reads(v, &[a, b]).into_iter().map(|r| r.0).collect();
    let flip = |c: Check| {
        let status = match c.status {
            Status::Fails => Status::Holds,
            Status::Holds => Status::Fails,
            s => s,
        };
        Check { status, ..c }
    };
    let mut c = flip(clearance(v, a, b, d, res, finest));
    // Over ranged params the pair found is at one value of them, and says
    // nothing of the rest: apart by d at every value fails.
    if c.status == Status::Holds && !read.is_empty() {
        (c.status, c.note) = (Status::Unknown, Some(AT_ONE_VALUE));
    }
    c
}

/// What an overlap or a contact over ranged params says found at one value
/// of them.
pub(crate) const AT_ONE_VALUE: &str = "over ranged params, a point or pair found shows it at one value of them, not every value. Only apart at every value fails";

/// A point inside the shape `f` (`inside`) or not: its field's bounds
/// there, rounded outward (real numbers), against 0; over the ranged
/// params it reads, for every value of them (value_span, halving pieces).
/// The coordinates are numbers, at their values (constants and params);
/// coordinates that read a ranged param or t move the point with it: the
/// field read there is then checked as a value of the params (at a moment,
/// t one value: its outward bounds there decide, never the field at the
/// coordinates' floats, which may round the point elsewhere). A coordinate not
/// finite somewhere in the scope (one overflowed) is no point to be inside
/// or outside of: unknown.
pub(super) fn point(v: &mut ModelView, f: u32, at: [u32; 3], inside: bool) -> Check {
    let cmp = if inside { Cmp::Lt } else { Cmp::Ge };
    let p = at.map(|n| v.doc.scalar(n, &v.ctx));
    // Each coordinate's bounds over the scope, before the field reads it:
    // a moving one's values are not p's, and one the field does not read
    // still says where the point is.
    let coords = at.map(|n| ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, n), &v.ctx, [I::point(0.0); 3], None));
    if coords.iter().any(|c| !(c.lo.is_finite() && c.hi.is_finite())) {
        return checked(Status::Unknown, Evidence::Point { at: p, bounds: (f64::NEG_INFINITY, f64::INFINITY) });
    }
    let reads_t = |v: &ModelView, ns: &[u32]| ns.iter().any(|&n| v.doc.g.depends_on(n, |x| x == Var::T));
    // Substitute a coordinate whenever its rounded point cannot stand
    // for its real value, including a constant spatial construction.
    let moves: Vec<(Var, u32)> = [Var::X, Var::Y, Var::Z]
        .into_iter()
        .zip(at)
        .enumerate()
        .filter_map(|(i, (axis, n))| (coords[i] != I::point(p[i]) || !ranged_reads(v, &[n]).is_empty() || reads_t(v, &[n])).then_some((axis, n)))
        .collect();
    let moving = !moves.is_empty();
    let p = std::array::from_fn(|i| if coords[i].contains(p[i]) { p[i] } else { coords[i].mid() });
    // Searched over a scope: a ranged param read, or t's span where it reads t.
    let ranged = !ranged_reads(v, &at).is_empty() || !ranged_reads(v, &[f]).is_empty() || (v.ctx.t_span.is_some() && (reads_t(v, &at) || reads_t(v, &[f])));
    let tape = if moving {
        let w = v.doc.g.subst(f, &moves);
        Tape::new(&v.doc.g, w)
    } else {
        Tape::new(&v.doc.g, f)
    };
    // The field reads the coordinates that do not move at p, where the
    // point is at the params' values (the evidence's).
    let iv = ink_core::eval::eval_interval_scope(&tape, &v.ctx, p.map(I::point), None);
    let evidence = Evidence::Point { at: p, bounds: (iv.lo, iv.hi) };
    // On a 0 the shape is read about the point (sign::at): flush faces inside
    // the material are in it, a surface is not.
    if !ranged && iv == I::point(0.0) {
        return match ink_core::sign::at(&tape, &v.ctx, p, 0.0) {
            Some(i) => checked(if i == inside { Status::Holds } else { Status::Fails }, evidence),
            None => checked(Status::Unknown, evidence),
        };
    }
    let c = Comparison::exact(cmp, 0.0);
    if c.holds_on(iv) {
        return checked(Status::Holds, evidence);
    }
    if !ranged && c.fails_on(iv) {
        return checked(Status::Fails, evidence);
    }
    if !ranged {
        return checked(Status::Unknown, evidence);
    }
    let result = ranged::value_span(&tape, &v.ctx, p, c, &v.proof_options(None, 0.0));
    checked(Status::of(result.status), Evidence::ValueSpan(result))
}

/// Several comparisons, all claimed: it fails where one does, holds where
/// all do.
pub(super) fn all_of(each: impl Iterator<Item = Status>) -> Status {
    each.fold(Status::Holds, |all, s| match (all, s) {
        (Status::Fails, _) | (_, Status::Fails) => Status::Fails,
        (Status::Holds, Status::Holds) => Status::Holds,
        _ => Status::Unknown,
    })
}

/// No point of b within d of a, exactly, with boxes no finer than
/// `finest` (or a res given where finer): the one search
/// (query::clearance), covering one shape proven bounded with boxes and
/// growing the other by d around each, either way round, since what it
/// needs, the grown shape's slope within d of the covered, may be bounded
/// for one and not the other (a polygon's sign has no slope bound at its
/// surface; a field promised 1.41 grows by more than one exact near its
/// faces). Only the covered shape within d of the other's box matters.
/// With neither proven bounded, all of space is tested as one box (two
/// parallel slabs settle there); else only the window is searched, and
/// only a witness pair (a violation) is an answer.
pub(super) fn clearance(v: &mut ModelView, pa: u32, pb: u32, d: f64, res: Option<f64>, finest: f64) -> Check {
    // A shape proven empty. One not proven bounded and empty in the window
    // says nothing about the rest of it: a bounded partner gives the
    // search its region.
    if [pa, pb].into_iter().any(|f| v.whole(f) && v.bounds(f).is_none()) {
        return checked(Status::Holds, Evidence::Vacuous("one of the shapes is empty"));
    }
    let finest = res.map_or(finest, |r| r.min(finest));
    let (ta, tb) = (Arc::new(Tape::new(&v.doc.g, pa)), Arc::new(Tape::new(&v.doc.g, pb)));
    let whole = v.whole(pa) || v.whole(pb);
    if !whole {
        let all = query::Way { grown: ta.clone(), covered: tb.clone(), lo: [f64::NEG_INFINITY; 3], hi: [f64::INFINITY; 3], swapped: false };
        if let (verdict @ Verdict::Proven, pair, h) = query::clearance(&[all], &v.ctx, d, finest, v.search_boxes) {
            return checked(Status::Holds, Evidence::ClearanceSearch { verdict, pair, res: h });
        }
    }
    let Some(ways) = clearance_ways(v, (pa, &ta), (pb, &tb), d) else {
        // No point of the covered shape is within d of the other.
        return checked(Status::Holds, Evidence::Vacuous("the shapes' boxes are farther apart than the clearance"));
    };
    if ways.is_empty() {
        return Check { unbounded: true, ..checked(Status::Unknown, Evidence::Unexamined) };
    }
    let (verdict, pair, h) = query::clearance(&ways, &v.ctx, d, finest, v.search_boxes);
    let status = match verdict {
        Verdict::Proven => Status::Holds,
        Verdict::Disproven { .. } => Status::Fails,
        Verdict::Unknown { .. } => Status::Unknown,
    };
    let mut c = checked(status, Evidence::ClearanceSearch { verdict, pair, res: h });
    if !whole && status != Status::Fails {
        if status == Status::Holds {
            (c.status, c.in_window) = (Status::Unknown, Some(Status::Holds));
        }
        c.unbounded = true;
    }
    c
}

/// The ways round a clearance d's search (`clearance`, query::Way): each
/// shape proven bounded covered with boxes (with neither, the second, in
/// the window), within d of the other's box where that is proven bounded.
/// None: the boxes are farther apart than d.
fn clearance_ways(v: &mut ModelView, (pa, ta): (u32, &Arc<Tape>), (pb, tb): (u32, &Arc<Tape>), d: f64) -> Option<Vec<query::Way>> {
    let whole = v.whole(pa) || v.whole(pb);
    let mut ways = Vec::new();
    for (cov, grown, swapped) in [(pb, pa, false), (pa, pb, true)] {
        if !(v.whole(cov) || !whole && !swapped) {
            continue;
        }
        // Neither proven bounded, and the covered one empty in the window:
        // nothing to search.
        let Some((mut lo, mut hi, _)) = v.bounds(cov) else { continue };
        if v.whole(grown) {
            let (glo, ghi, _) = v.bounds(grown).expect("not proven empty");
            for k in 0..3 {
                lo[k] = lo[k].max((glo[k] - d).next_down());
                hi[k] = hi[k].min((ghi[k] + d).next_up());
            }
        }
        if (0..3).any(|k| lo[k] > hi[k]) {
            return None;
        }
        let grown = if swapped { tb.clone() } else { ta.clone() };
        ways.push(query::Way { grown, covered: v.clipped(cov), lo, hi, swapped });
    }
    Some(ways)
}

/// The minimum distance between two solids as a bracket (`distance`).
pub struct Distance {
    /// Proven: every point of each at least this from the other.
    pub lower: f64,
    /// A pair of points confirmed inside a and b, and their distance
    /// rounded up: an upper bound. None: no pair was confirmed.
    pub upper: Option<(f64, V3, V3)>,
    /// The search ended short of the accuracy (query::Ended::Finest or
    /// Limit) with d, the clearance it was proving, this: the bracket is as
    /// tight as the search goes at its budget.
    pub undecided: Option<f64>,
    /// Boxes the search tested.
    pub boxes: usize,
    /// The call's time limit stopped the search (ink_kernel::stop).
    pub stopped: bool,
}

/// The minimum distance between a and b, bracketed to the model's
/// accuracy by one search (query::min_distance): a clearance's, covering
/// each shape proven bounded with boxes (`clearance_ways`), its d shrinking
/// to each pair closer than it that the search finds, confirmed inside
/// both, less the accuracy, from the farthest the targets' boxes are apart
/// (where it finds none, again from a pair proven inside them).
/// The lower bound is what it proved: nothing where neither shape is proven
/// bounded (only the window was searched) or the targets read inputs that
/// are not exact (samples, solved values not certified), on which nothing
/// is proven. Each bound is sound whatever the fields: a field that is not
/// a distance (a solid cut by a plane) never stands in for one.
pub fn distance(v: &mut ModelView, fa: u32, fb: u32, opts: ClaimOptions) -> Result<Distance, InkError> {
    let (Some((alo, ahi, _)), Some((blo, bhi, _))) = (v.bounds(fa), v.bounds(fb)) else { return Err(InkError::msg("a target is empty")) };
    let far = ink_core::interval::norm(std::array::from_fn(|k| (ahi[k] - blo[k]).max(bhi[k] - alo[k]).next_up()));
    let (ta, tb) = (Arc::new(Tape::new(&v.doc.g, fa)), Arc::new(Tape::new(&v.doc.g, fb)));
    let tol = v.accuracy();
    let finest = opts.res.map_or(v.finest(), |r| r.min(v.finest())).min(if tol > 0.0 { 0.5 * tol } else { f64::INFINITY });
    let ways = clearance_ways(v, (fa, &ta), (fb, &tb), far).unwrap_or_default();
    let mut found = query::min_distance(&ways, &v.ctx, far, None, tol, finest, query::BOXES);
    // Undecided with no pair found: the search's seeks walk along a
    // field's slope, which a field flat away from its solid (a union with
    // a constant) does not have there. Again from a pair of points proven
    // inside them (query::inside_pair). Proven with none, a target has no
    // point inside (its field is never below 0), and there is no pair to
    // find.
    if found.upper.is_none()
        && matches!(found.ended, Some(query::Ended::Finest | query::Ended::Limit))
        && let Some(seed) = query::inside_pair(&ta, &tb, &v.ctx, (alo, ahi), (blo, bhi), finest)
    {
        let ways = clearance_ways(v, (fa, &ta), (fb, &tb), seed.0).unwrap_or_default();
        let again = query::min_distance(&ways, &v.ctx, far, Some(seed), tol, finest, query::BOXES);
        found = query::MinDistance { lower: again.lower.max(found.lower), tested: found.tested + again.tested, ..again };
    }
    let proven = (v.whole(fa) || v.whole(fb)) && !v.doc.inexact(&[fa, fb], &v.certified()).any();
    let stopped = found.ended == Some(query::Ended::Stopped);
    let undecided = matches!(found.ended, Some(query::Ended::Finest | query::Ended::Limit)).then(|| found.upper.map_or(far, |u| u.0) - tol);
    Ok(Distance { lower: if proven { found.lower } else { 0.0 }, upper: found.upper, undecided, boxes: found.tested, stopped })
}
