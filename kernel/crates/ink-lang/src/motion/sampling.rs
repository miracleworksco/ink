use super::*;

/// First piece of the span whose rate is still unbounded after subdivision.
pub(super) fn may_jump(span: (f64, f64), mut bounded: impl FnMut((f64, f64)) -> bool) -> Option<(f64, f64)> {
    let least = (span.1 - span.0) / MAX_STEPS as f64;
    let mut todo = vec![span];
    while let Some(s) = todo.pop() {
        if bounded(s) {
            continue;
        }
        let m = ink_core::eval::span_mid(s);
        if s.1 - s.0 <= least || !(s.0 < m && m < s.1) {
            return Some(s);
        }
        todo.extend([(m, s.1), (s.0, m)]);
    }
    None
}

/// A witness's value of param k (t ranged: its moment), where it has one.
pub(super) fn witness_t(e: &Evidence, k: u32) -> Option<f64> {
    use ink_kernel::query::Verdict;
    let of = |ps: &[(u32, f64)]| ps.iter().find(|p| p.0 == k).map(|p| p.1);
    match e {
        Evidence::Empty(p) | Evidence::Overlaps(p, _) => match &p.verdict {
            Verdict::Disproven { params, .. } => of(params),
            _ => None,
        },
        Evidence::ValueSpan(v) => v.witness.as_ref().and_then(|(at, _)| of(&at.params)),
        Evidence::ClearanceSearch { verdict: ink_kernel::query::Verdict::Disproven { params, .. }, .. } => of(params),
        Evidence::Touches { overlap, gap } => witness_t(overlap, k).or_else(|| witness_t(gap, k)),
        _ => None,
    }
}

/// The claim with each node it reads mapped by `f`, for the kinds a band
/// ranges (None: a package's metric, an integral or a mean: measured).
pub(super) fn with_nodes(c: &Claim, f: &mut impl FnMut(u32) -> u32) -> Option<Claim> {
    Some(match *c {
        Claim::Empty(a) => Claim::Empty(f(a)),
        Claim::EmptyBounds { certain, possible } => Claim::EmptyBounds { certain: f(certain), possible: f(possible) },
        Claim::Range(a) => Claim::Range(f(a)),
        Claim::Disjoint(a, b) => Claim::Disjoint(f(a), f(b)),
        Claim::Subset(a, b) => Claim::Subset(f(a), f(b)),
        Claim::Overlaps(a, b) => Claim::Overlaps(f(a), f(b)),
        Claim::Touches(a, b) => Claim::Touches(f(a), f(b)),
        Claim::Point { f: a, at, inside } => Claim::Point { f: f(a), at: at.map(&mut *f), inside },
        Claim::Clearance(a, b, d) => Claim::Clearance(f(a), f(b), d),
        Claim::Value(a, cmp, x) => Claim::Value(f(a), cmp, x),
        Claim::Metric { .. } | Claim::Integral { .. } | Claim::Mean { .. } => return None,
    })
}

/// The claim checked at claims::SAMPLES moments of `span`, evenly spaced,
/// its ends included (where a motion starts and stops: a part landing on
/// its stop at the last moment), in order, with its bound evaluated
/// there: an estimate (a failure at a moment is still a witness: the
/// first stops it). A moment what it rests on has no value at (a solve
/// that does not solve there, a simulation stopped at a dead point, a
/// rate not defined where the motion starts) is undecided, with why.
/// Where a provider of what it reads chooses the steps
/// (ParamProvider::steps), at those.
pub(super) fn measured(v: &mut ModelView, claim: Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Check {
    use ink_kernel::ranged::Status as S;
    let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
    let chosen: Vec<Result<crate::Steps, crate::Refused>> =
        v.doc.solves_read(&roots).into_iter().filter_map(|i| v.doc.ext[v.doc.solves[i].ext].form.provider(&v.doc)?.steps(&v.ctx.params, span)).collect();
    match chosen.into_iter().collect::<Result<Vec<_>, _>>().map(merged).transpose() {
        Some(Ok(steps)) => return measured_at(v, &claim, bound, &roots, steps, span, opts),
        Some(Err(why)) => {
            let m = Motion { undecided: vec![span], ..Motion::none(S::Unknown, span) };
            let c = finish(Some(checked(Status::Unknown, Evidence::Unsolved(InkError { open: why.open, ..InkError::msg(why.why) }))), m, false, None);
            return Check { motion: c.motion.map(|s| Stepped { basis: Some("not stepped: the steps rest on what has no value (unavailable)".into()), ..s }), ..c };
        }
        None => {}
    }
    let mut moments = claims::samples(span);
    moments.dedup();
    let mut m = Motion::none(S::Proven, span);
    let mut deciding: Option<Check> = None;
    for &t in &moments {
        if ink_kernel::stop::reached() {
            past_limit(&mut m, (t, t));
            continue;
        }
        let refused = match v.set_moment(t) {
            Err(e) => Some(e),
            Ok(()) => v.refusal(&roots),
        };
        let solved = refused.is_none();
        m.steps += 1;
        let c = match refused {
            Some(e) => checked(Status::Unknown, Evidence::Unsolved(e)),
            None => match validate_at(v, &claim, bound, opts) {
                Some(o) => o,
                None => claims::decide(v, v.doc.claim_at(&claim, bound, &v.ctx), opts).unwrap_or_else(|e| checked(Status::Unknown, Evidence::Unsolved(e))),
            },
        };
        match c.status {
            Status::Fails => {
                (m.status, m.t, deciding) = (S::Disproven, Some(t), Some(certified(v, &claim, bound, c, opts)));
                break;
            }
            Status::Holds => {
                deciding.get_or_insert(c);
            }
            Status::Unknown => {
                m.undecided.push((t, t));
                m.status = S::Unknown;
                m.stopped += usize::from(solved && ink_kernel::stop::reached());
                deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(c));
            }
        }
    }
    // There is no interior sample between adjacent floats. Passing just
    // the endpoints cannot even estimate what happens inside this span.
    if span.0 < span.1 && moments.len() <= 2 && m.status == S::Proven {
        m.status = S::Unknown;
        m.undecided.push(span);
        deciding = None;
    }
    let step = (moments.len() > 1).then(|| (span.1 - span.0) / (moments.len() - 1) as f64);
    finish(deciding, m, false, step)
}

/// `claim` (with its `bound`, a field) judged at the view's moment
/// (claims::validate) where it, its bound included, reads a solve on t,
/// before anything is put in its place (the bound's number there:
/// Doc::claim_at): over the motion it was judged before it was asked
/// (claims::ask_in) with what such a solve sets ranged where it could be;
/// where not (one only measured), it is judged here, at each moment
/// measured, as it is asked.
fn validate_at(v: &mut ModelView, claim: &Claim, bound: Option<u32>, opts: ClaimOptions) -> Option<Check> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    if !v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t) {
        return None;
    }
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let at = Ctx { t_span: None, ..v.ctx.clone() };
    claims::validate(v, claim, &roots, tol, opts.res, at, None)
}

/// A moment (or step) of a measured claim not checked: the call's time
/// limit had passed (ink_kernel::stop). A moment checked and left unknown
/// once it has is the limit's too (its searches stopped); one a solve
/// refuses is not.
fn past_limit(m: &mut Motion, s: (f64, f64)) {
    m.undecided.push(s);
    m.status = ink_kernel::ranged::Status::Unknown;
    m.stopped += 1;
}

/// Why a step a model chose is not covered: `not_followed`'s reason wherever a provider's steps are read.
pub const MODEL_GAP: &str = "what the model works out changes faster than its steps follow, or jumps";

/// What is said where moments read over a motion do not cover it: in `gaps` (steps of the motion, in
/// order; those that meet said as one), for `why`. One wording for every answer from moments (a claim
/// measured, a clash, a study loaded from the model, a rating). None: no gap.
pub fn not_followed(why: &str, gaps: &[(f64, f64)]) -> Option<String> {
    let mut joined: Vec<(f64, f64)> = Vec::new();
    for &s in gaps {
        match joined.last_mut() {
            Some(l) if l.1 >= s.0 => l.1 = l.1.max(s.1),
            _ => joined.push(s),
        }
    }
    let gaps = &joined[..];
    let (first, last) = (gaps.first()?, gaps.last()?);
    // To the digits that tell the two ends apart.
    let sig = |x: f64| if crate::claims::sig(first.0) == crate::claims::sig(last.1) { format!("{x}") } else { crate::claims::sig(x) };
    let over = if gaps.len() == 1 { format!("t {} to {}", sig(first.0), sig(first.1)) } else { format!("{} steps from t {} to {}", gaps.len(), sig(first.0), sig(last.1)) };
    Some(format!(
        "not settled between the moments read, over {over}: {why}. A yes at the moments says nothing of what lies between them there: ask a shorter span (--t-span A:B), or slow what \
         changes"
    ))
}

/// The steps several providers chose over one span, each a claim reads
/// (an effort and a reaction, two reactions): every end any chose, a step
/// unchecked or loose where it lies in one some provider left so, their
/// bases each once. None: none chose.
pub(crate) fn merged(all: Vec<crate::Steps>) -> Option<crate::Steps> {
    let mut it = all.into_iter();
    let first = it.next()?;
    let rest: Vec<crate::Steps> = it.collect();
    if rest.is_empty() {
        return Some(first);
    }
    let every: Vec<&crate::Steps> = std::iter::once(&first).chain(&rest).collect();
    let mut ends: Vec<f64> = every.iter().flat_map(|s| s.ends.iter().copied()).collect();
    ends.sort_by(f64::total_cmp);
    ends.dedup();
    let within = |s: (f64, f64), of: &[(f64, f64)]| of.iter().any(|o| o.0 <= s.0 && s.1 <= o.1);
    let steps: Vec<(f64, f64)> = ends.windows(2).map(|w| (w[0], w[1])).collect();
    let unchecked = steps.iter().copied().filter(|&s| every.iter().any(|x| within(s, &x.unchecked))).collect();
    let loose = steps.iter().copied().filter(|&s| every.iter().any(|x| within(s, &x.loose))).collect();
    let mut bases: Vec<String> = Vec::new();
    for b in every.iter().map(|s| s.basis.clone()) {
        if !bases.contains(&b) {
            bases.push(b);
        }
    }
    Some(crate::Steps { ends, unchecked, loose, basis: bases.join("; and at the steps of what else it reads: ") })
}

/// The claim measured at the span's ends and the middles of the steps a
/// provider chose (`steps`), in order: failing at one, it fails, and its
/// witness is the failing moment of least margin (a value claim's value
/// furthest past its bound: what a design fixed against it must clear;
/// the first where there is no margin to compare);
/// a moment it rests on a solve with no value at, undecided, with the
/// solve's error; the steps the provider could not choose as it says,
/// unchecked or loose. Holding, its evidence is the moment where it came
/// nearest to failing (a value claim's value of least margin: its greatest
/// under an upper bound, its least over a lower).
fn measured_at(v: &mut ModelView, claim: &Claim, bound: Option<u32>, roots: &[u32], steps: crate::Steps, span: (f64, f64), opts: ClaimOptions) -> Check {
    let mut m = Motion { unchecked: steps.unchecked.clone(), loose: steps.loose.clone(), ..Motion::none(ink_kernel::ranged::Status::Proven, span) };
    let skipped = |s: (f64, f64)| steps.unchecked.contains(&s) || steps.loose.contains(&s);
    let (mut deciding, mut least): (Option<Check>, f64) = (None, f64::INFINITY);
    // Its worst failing moment so far: its check, the moment, and its margin.
    let mut worst: Option<(Check, f64, f64)> = None;
    let margin = |c: &Check| {
        let x = match &c.evidence {
            Evidence::Value(x) => Some(*x),
            Evidence::ValueSpan(vs) => match c.status {
                Status::Fails => vs.witness.as_ref().map(|(_, x)| *x),
                Status::Holds => Some(match claim {
                    Claim::Value(_, Cmp::Lt | Cmp::Le, _) => vs.range.1,
                    _ => vs.range.0,
                }),
                Status::Unknown => None,
            },
            _ => None,
        };
        match (x, claim) {
            (Some(x), Claim::Value(_, Cmp::Lt | Cmp::Le, bound)) => bound - x,
            (Some(x), Claim::Value(_, _, bound)) => x - bound,
            (Some(x), _) => x,
            _ => f64::INFINITY,
        }
    };
    // What it reads of t itself, not through a provider: the provider's steps follow what the provider
    // rests on, not that, so it is shown over each step (`shown_over`) or the step is not settled. A
    // claim about shapes reads where they are, and a part a provider moves (a simulated coordinate,
    // read about each moment as the provider's own local motion) may pass through another between two
    // of its steps: every step of it is shown so.
    let reads_t = match claim {
        Claim::Value(..) => reads_t_itself(v, roots),
        _ => roots.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)),
    };
    let mut own: Vec<(f64, f64)> = Vec::new();
    let (mut asked, mut shown) = (0, 0);
    let middles = steps.ends.windows(2).map(|w| (w[0], w[1])).filter(|&s| !skipped(s)).map(|s| (ink_core::eval::span_mid(s), s));
    let moments: Vec<(f64, (f64, f64))> = std::iter::once((span.0, (span.0, span.0))).chain(middles).chain([(span.1, (span.1, span.1))]).collect();
    // What was left undecided before it failed (what it says of its steps then).
    let mut before = 0;
    for (t, s) in moments {
        if ink_kernel::stop::reached() {
            if worst.is_some() {
                break;
            }
            past_limit(&mut m, s);
            continue;
        }
        m.steps += 1;
        let refused = match v.set_moment(t) {
            Err(e) => Some(e),
            Ok(()) => v.refusal(roots),
        };
        if let Some(e) = refused {
            m.undecided.push(s);
            deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(checked(Status::Unknown, Evidence::Unsolved(e))));
            continue;
        }
        if let Some(o) = validate_at(v, claim, bound, opts) {
            m.undecided.push(s);
            deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(o));
            continue;
        }
        let c = v.doc.claim_at(claim, bound, &v.ctx);
        let decided = match banded(v, &c, roots, t, opts) {
            Some(c) => Ok(c),
            None => claims::decide(v, c, opts),
        };
        let Ok(c) = decided else {
            m.undecided.push(s);
            continue;
        };
        // Failed already: only a moment failing further is looked for.
        if worst.is_some() && c.status != Status::Fails {
            continue;
        }
        match c.status {
            Status::Fails => {
                let x = margin(&c);
                if worst.is_none() {
                    before = m.undecided.len();
                }
                if worst.as_ref().is_none_or(|w| x < w.2) {
                    worst = Some((c, t, x));
                }
            }
            Status::Unknown => {
                m.undecided.push(s);
                m.stopped += usize::from(ink_kernel::stop::reached());
                deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(c));
            }
            Status::Holds => {
                if reads_t && s.0 < s.1 {
                    // With none of the first BETWEEN_FIRST shown, the rest would not all be: not asked.
                    let over = match asked >= BETWEEN_FIRST && shown == 0 {
                        true => Over::Not,
                        false => shown_over(v, claim, bound, s, opts),
                    };
                    asked += 1;
                    match over {
                        Over::Shown => shown += 1,
                        Over::Not => own.push(s),
                        Over::Fails(c) => {
                            (m.status, m.t) = (ink_kernel::ranged::Status::Disproven, c.ranged.iter().find_map(|r| witness_t(&c.evidence, r.0)).or(Some(t)));
                            return Check { motion: Some(Stepped { basis: Some(steps.basis.clone()), ..Stepped::of(m, false, None) }), ..*c };
                        }
                    }
                }
                let x = margin(&c);
                if deciding.as_ref().is_none_or(|d| d.status == Status::Holds && x < least) {
                    (deciding, least) = (Some(c), x);
                }
            }
        }
    }
    if let Some((c, t, _)) = worst {
        m.undecided.truncate(before);
        (m.status, m.t) = (ink_kernel::ranged::Status::Disproven, Some(t));
        let c = certified(v, claim, bound, c, opts);
        return Check { motion: Some(Stepped { basis: Some(steps.basis.clone()), ..Stepped::of(m, false, None) }), ..c };
    }
    let between = not_followed("besides what its steps follow, it reads t itself, and that is not shown to keep it over the step", &own);
    m.loose.extend(own);
    if !(m.unchecked.is_empty() && m.undecided.is_empty() && m.loose.is_empty()) {
        m.status = ink_kernel::ranged::Status::Unknown;
    }
    let c = finish(deciding, m, false, None);
    Check { motion: c.motion.map(|s| Stepped { basis: Some(steps.basis), between: between.filter(|_| c.status == Status::Unknown), ..s }), ..c }
}
