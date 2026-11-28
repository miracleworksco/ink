use super::*;

/// The claim at every moment of `span` (`bound`: its bound, as a field, if
/// it reads params): with t ranged (`ranged`), else measured.
/// The view is left at its moment.
pub(crate) fn over_span(v: &mut ModelView, claim: Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    let t0 = v.ctx.t;
    let bound = bound.filter(|&b| claims::reads_t(v, &[b]));
    Asked::Done(()).then(v, move |v, ()| {
        // A body that only moves rigidly, where what is asked is the same of it moved so or carried by its
        // placement: asked once in its own frame.
        let rigid = crate::carried::over_span(v, &claim, bound, span, opts);
        v.set_moment(t0)?;
        if let Some(c) = rigid? {
            return Ok(Asked::Done(c));
        }
        let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
        let solved = v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t);
        // The sampled answer first: cheap, always there, an estimate (a
        // failure at a moment is a witness, guaranteed: nothing to prove,
        // but where it reads a solve on t, found there to its solver's
        // accuracy: a witness with t ranged is confirmed over its root's
        // enclosure).
        let sampled = measured(v, claim.clone(), bound, span, opts);
        v.set_moment(t0)?;
        if sampled.status == Status::Fails && sampled.accuracy.is_some() && !solved {
            return Ok(Asked::Done(sampled));
        }
        let keep = |sampled: Check, why: String, answer: Option<Check>| Check { motion: sampled.motion.map(|s| Stepped { left: Some(why), left_answer: answer.map(Box::new), ..s }), ..sampled };
        // Beside other searches (interfere's moving pairs): the sampled
        // answer is the pair's, and t ranged is tried once every pair has
        // one, with the time left (`Upgrade`), replacing it where it decides.
        // Where the sampled answer would stand as a yes, what lies between its moments is asked
        // (`between`), once every pair has an answer.
        let later = move |v: &mut ModelView, asked: &Claim, sampled: &Check, limit: bool| -> Result<Option<Check>, InkError> {
            if sampled.status != Status::Holds || sampled.motion.as_ref().is_none_or(|s| s.basis.is_some()) {
                return Ok(None);
            }
            if ink_kernel::stop::reached() {
                return Ok(None);
            }
            let (why, answer) = sampled.motion.as_ref().map_or((String::new(), None), |s| (s.left.clone().unwrap_or_default(), s.left_answer.as_deref().cloned()));
            let alone = ClaimOptions { beside: false, ..opts };
            let c = between(v, asked, bound, span, alone, sampled.clone(), why, answer, limit);
            v.set_moment(t0)?;
            let c = c?;
            // Past the call's time limit, an unfinished upgrade leaves the
            // earlier sampled answer standing instead of overwriting it with
            // unknown.
            if c.status == Status::Unknown && ink_kernel::stop::reached() {
                return Ok(None);
            }
            // The sampled answer, as it was, where every step held; one decided with t ranged (its
            // period asked in pieces) is finished as any such answer is.
            let sampled = c.status == Status::Holds && !c.motion.as_ref().is_some_and(|s| s.ranged);
            Ok(Some(if sampled { c } else { claims::finish_later(v, asked, bound, c, opts) }))
        };
        if opts.beside && v.doc.claim_motion_inputs(&claim, None).iter().any(|&f| own_field_moves(v, f)) {
            let kept = keep(sampled, "a shape's own field reads t (a belt's teeth running round, a chain's links). Beside other pairs, finding its region over the span is a search of its own, so t ranged is tried once every pair has an answer".into(), None);
            let (asked, first) = (claim.clone(), kept.clone());
            // Once every pair has an answer, t ranged is tried of it alone, as a claim's is.
            v.upgrades.push(Upgrade {
                start: Box::new(move |v| {
                    let alone = ClaimOptions { beside: false, ..opts };
                    let r = ranged(v, &asked, bound, span, alone);
                    v.set_moment(t0)?;
                    Ok(Later::Now(match r? {
                        Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                        // Its steps are asked whatever the whole search left: a belt's are short.
                        Ranged::Left(..) => later(v, &asked, &first, false)?,
                        Ranged::Wait(..) => unreachable!("only beside other searches is t ranged searched apart"),
                    }))
                }),
            });
            return Ok(Asked::Done(kept));
        }
        if opts.beside {
            let asked = claim.clone();
            let kept = keep(sampled, "t ranged is tried once every pair has an answer, with the time left: it did not decide it".into(), None);
            let first = kept.clone();
            v.upgrades.push(Upgrade {
                start: Box::new(move |v| {
                    let r = ranged(v, &asked, bound, span, opts);
                    v.set_moment(t0)?;
                    Ok(match r? {
                        Ranged::Decided(c) => Later::Now(Some(claims::finish_later(v, &asked, bound, c, opts))),
                        Ranged::Left(_, a) => Later::Now(later(v, &asked, &first, a.is_some_and(|a| at_limit(&a.evidence)))?),
                        Ranged::Wait(search, whole, read, world, moments, charged) => {
                            let (again, first) = (asked.clone(), first.clone());
                            let asked_in = (search.0.clone(), search.1.clone());
                            let then = move |v: &mut ModelView, p: Found| -> Result<Option<Check>, InkError> {
                                // Over a drive's value only a hold is an
                                // answer: else the sampled one stands.
                                if !moments {
                                    return Ok((p.status() == ink_kernel::ranged::Status::Proven).then(|| {
                                        let held = held_beside(v, &p, &asked_in, None).map(|h| match h.charged(charged) {
                                            claims::Held::Outright => claims::Held::Outright,
                                            claims::Held::ToAccuracy { .. } => claims::Held::ToAccuracy { near: None },
                                            h @ claims::Held::SetAside { .. } => h,
                                        });
                                        let c = Check { ranged: read, held, note: Some(DRIVE), ..claims::in_window(Status::Holds, whole, Evidence::Empty(p)) };
                                        match decided(c, span, false) {
                                            Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                                            _ => None,
                                        }
                                    }).flatten());
                                }
                                let held = held_beside(v, &p, &asked_in, world.as_ref()).map(|h| h.charged(charged));
                                // A witness in the frame shown none in the
                                // world: the sampled answer stands.
                                let e = match &world {
                                    Some(w) if p.status() != ink_kernel::ranged::Status::Unknown => match w.confirmed(v, Evidence::Empty(p.clone())) {
                                        Some(e) => e,
                                        None => return Ok(None),
                                    },
                                    _ => Evidence::Empty(p.clone()),
                                };
                                let c = Check { ranged: read, held, ..claims::in_window(Status::of(p.status()), whole, e) };
                                Ok(match decided(c, span, false) {
                                    Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                                    _ => None,
                                })
                            };
                            Later::Wait(
                                search,
                                Box::new(move |v: &mut ModelView, p: Found| {
                                    let limit = at_limit(&Evidence::Empty(p.clone()));
                                    match then(v, p)? {
                                        Some(c) => Ok(Some(c)),
                                        None => later(v, &again, &first, limit),
                                    }
                                }),
                            )
                        }
                    })
                }),
            });
            return Ok(Asked::Done(kept));
        }
        // Then t ranged, with the time left, to prove it or find a witness;
        // undecided, the sampled answer stands, saying why.
        let r = ranged(v, &claim, bound, span, opts);
        v.set_moment(t0)?;
        Ok(Asked::Done(match r? {
            Ranged::Decided(c) => c,
            // Read at steps a provider chose, what lies between them is the provider's to answer for
            // (a step it does not follow is loose or unchecked: `measured_at`).
            Ranged::Left(why, answer) if sampled.status == Status::Holds && sampled.motion.as_ref().is_some_and(|s| s.basis.is_none()) => {
                let limit = answer.as_ref().is_some_and(|a| at_limit(&a.evidence));
                let c = between(v, &claim, bound, span, opts, sampled, why, answer, limit);
                v.set_moment(t0)?;
                c?
            }
            Ranged::Left(why, answer) => keep(sampled, why, answer),
            Ranged::Wait(..) => unreachable!("only beside other searches is t ranged searched apart"),
        }))
    })
}

/// How many times a claim that holds at its moments is asked with t ranged over a piece of the span, to
/// cover what lies between them (`between`).
const BETWEEN_ASKS: usize = 64;

/// The pieces one period of a drive is asked in where its whole search ended at its limit of boxes.
const PERIOD_PIECES: usize = 16;

/// A claim over one period of its drive (`window`, the values of t's param k), asked in PERIOD_PIECES
/// pieces, each a search of its own share of the boxes one search may test: a mesh's teeth through one
/// pitch are a few boxes deep in each sixteenth and millions in the whole. Holds where every piece
/// does (outright where every one does), fails with the first that fails, and is the first undecided
/// one's otherwise. The view's scope is left over the window.
pub(super) fn in_pieces(v: &mut ModelView, asked: &Claim, opts: ClaimOptions, k: u32, window: (f64, f64)) -> Result<Check, InkError> {
    let boxes = std::mem::replace(&mut v.search_boxes, ink_kernel::query::BOXES / PERIOD_PIECES);
    let at = |i: usize| if i == PERIOD_PIECES { window.1 } else { window.0 + (window.1 - window.0) * i as f64 / PERIOD_PIECES as f64 };
    let (mut held, mut outright) = (None, true);
    let mut ended = None;
    for i in 0..PERIOD_PIECES {
        let piece = (at(i), at(i + 1));
        v.ctx.p_spans.narrow(k, piece);
        v.ctx.params[k as usize] = ink_core::eval::span_mid(piece);
        v.ctx.derive();
        match claims::decide(v, asked.clone(), opts) {
            Ok(c) if c.status == Status::Holds => {
                outright &= matches!(c.held, Some(claims::Held::Outright));
                held.get_or_insert(c);
            }
            other => {
                ended = Some(other);
                break;
            }
        }
    }
    v.search_boxes = boxes;
    v.ctx.p_spans.narrow(k, window);
    v.ctx.params[k as usize] = ink_core::eval::span_mid(window);
    v.ctx.derive();
    match ended {
        Some(other) => other,
        None => {
            let c = held.expect("a piece");
            Ok(Check { held: c.held.map(|h| if outright { h } else { claims::Held::ToAccuracy { near: None } }), ..c })
        }
    }
}

/// How far, in its own sizes, a body may travel within a step that is asked: a search of BETWEEN_BOXES
/// boxes halves t a few times, not the eight a whole search does (`FAR`).
const BETWEEN_FAR: f64 = 16.0;

/// Whether a search ended at its limit of boxes, or would have passed it at its rate: one given fewer
/// boxes ends sooner.
fn at_limit(e: &Evidence) -> bool {
    use ink_kernel::query::{Ended, Verdict};
    let limit = |v: &Verdict| matches!(v, Verdict::Unknown { ended: Ended::Limit, .. });
    match e {
        Evidence::Empty(p) | Evidence::Overlaps(p, _) => limit(&p.verdict),
        Evidence::ClearanceSearch { verdict, .. } => limit(verdict),
        Evidence::Touches { overlap, gap } => at_limit(overlap) || at_limit(gap),
        _ => false,
    }
}

/// The steps asked before any is shown: with none of them shown, no more are asked.
pub(super) const BETWEEN_FIRST: usize = 4;

/// The boxes each such ask may test: its search ends undecided past them, and the step is halved or
/// left. With BETWEEN_ASKS it is all the work a claim spends on what lies between its moments, the
/// same on every machine.
pub(super) const BETWEEN_BOXES: usize = 1 << 13;

/// A claim that holds at each moment read over `span` (`sampled`) and that t ranged over the whole span
/// left undecided (`why`, `answer`), asked of what lies between the moments: with t ranged over each step
/// between two of them, a step left undecided halved and asked again, to BETWEEN_ASKS asks. Over a short
/// step a part travels little, so the question is whether the gap the claim has at the moments is more
/// than what it reads can move within the step. Every step held: the sampled answer stands, an estimate
/// (each step's proof is to its own accuracy, and they are not made one). A witness in a step: that
/// answer, a no. A step left: not settled over it, saying so (`not_followed`).
#[allow(clippy::too_many_arguments)]
fn between(
    v: &mut ModelView,
    claim: &Claim,
    bound: Option<u32>,
    span: (f64, f64),
    opts: ClaimOptions,
    sampled: Check,
    why: String,
    answer: Option<Check>,
    limit: bool,
) -> Result<Check, InkError> {
    let mut moments = claims::samples(span);
    moments.dedup();
    let mut todo: std::collections::VecDeque<(f64, f64)> = moments.windows(2).map(|w| (w[0], w[1])).collect();
    let (mut asks, mut gaps, mut last, mut stopped) = (0, Vec::new(), None, 0);
    // The search over the whole span ended at its limit of boxes (or would have passed it): one over a
    // step, given BETWEEN_BOXES, ends sooner. None is asked. First, where the motion repeats, its one
    // period is asked in pieces, each a search of its own (`in_pieces`): that answer is the claim's
    // over the span, whatever the whole search left. Past the call's time limit nothing more is
    // asked: every step is left, the limit's.
    if !ink_kernel::stop::reached()
        && let Ranged::Decided(c) = ranged(v, claim, bound, span, ClaimOptions { in_pieces: true, ..opts })?
    {
        return Ok(c);
    }
    if ink_kernel::stop::reached() {
        stopped += 1;
        gaps.extend(todo.drain(..));
    } else if limit {
        gaps.extend(todo.drain(..));
    }
    let boxes = std::mem::replace(&mut v.search_boxes, BETWEEN_BOXES);
    let first = todo.len();
    let (mut held, mut left) = (0, 0);
    while let Some(s) = todo.pop_front() {
        // Each step left is halved, two asks more: once those pass what is allowed not every step
        // can be shown, and with none of the first BETWEEN_FIRST shown the rest would not all be. No
        // more is asked.
        if first + 2 * left > BETWEEN_ASKS || (asks == BETWEEN_FIRST.min(first) && held == 0) {
            gaps.push(s);
            gaps.extend(todo.drain(..));
            break;
        }
        if asks >= BETWEEN_ASKS || ink_kernel::stop::reached() {
            stopped += usize::from(asks < BETWEEN_ASKS);
            gaps.push(s);
            continue;
        }
        asks += 1;
        // A body that crosses many of its own sizes within the step is not followed by a search of
        // BETWEEN_BOXES boxes: not asked, as one it left.
        let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
        let asked = match sweeps(v, &roots, s, BETWEEN_FAR) {
            Some((travel, own)) => Ok(Ranged::Left(
                format!("a body it reads travels {travel:.1} within the step, over {BETWEEN_FAR} times its size ({own:.2}): further than the search a step is given follows"),
                None,
            )),
            // The claim as it is: its parts were asked over the whole span (`ranged`), and a step is no
            // more of them.
            None => ranged_whole(v, claim, bound, s, opts),
        };
        if asked.is_err() {
            v.search_boxes = boxes;
        }
        match asked? {
            Ranged::Decided(c) if c.status == Status::Holds => held += 1,
            // A witness in a step is one of the span asked.
            Ranged::Decided(c) => {
                v.search_boxes = boxes;
                return Ok(Check { motion: c.motion.map(|s| Stepped { motion: Motion { span, ..s.motion }, ..s }), ..c });
            }
            Ranged::Left(w, a) => {
                left += 1;
                let mid = ink_core::eval::span_mid(s);
                match s.0 < mid && mid < s.1 {
                    true => todo.extend([(s.0, mid), (mid, s.1)]),
                    false => gaps.push(s),
                }
                last = Some((w, a));
            }
            Ranged::Wait(..) => unreachable!("only beside other searches is t ranged searched apart"),
        }
    }
    v.search_boxes = boxes;
    gaps.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut joined: Vec<(f64, f64)> = Vec::new();
    for s in gaps {
        match joined.last_mut() {
            Some(l) if l.1 >= s.0 => l.1 = l.1.max(s.1),
            _ => joined.push(s),
        }
    }
    let reason = match (stopped, limit) {
        (0, false) => "what it reads may move further within a step than the room the claim has at the moments",
        (0, true) => "the search over the whole span ended at its limit of boxes, and a search over a step between two moments is given far fewer",
        _ => "the call's time limit passed before every step between them was asked",
    };
    let Some(gap) = not_followed(reason, &joined) else {
        let held = format!("{why}. Between the moments, t ranged over each of {asks} steps holds");
        return Ok(Check { motion: sampled.motion.map(|s| Stepped { left: Some(held), left_answer: answer.map(Box::new), ..s }), ..sampled });
    };
    // What the last step left says why where it is a contact, or was not searched (an estimate's
    // condition, a claim no range shows): a search's own account of itself names its usual budget,
    // not this one.
    let (w, a) = last.unwrap_or((why, answer));
    let contact = a.as_ref().is_some_and(|a| a.note == Some(claims::AT_ONE_VALUE));
    let (w, a) = match (contact || limit || a.is_none(), stopped) {
        (true, _) | (_, 1..) => (w, a),
        _ => (format!("each step was asked with t ranged over it, its search given {BETWEEN_BOXES} boxes, and {} of them left it undecided", joined.len()), None),
    };
    let motion = sampled.motion.map(|s| Stepped {
        motion: Motion { status: ink_kernel::ranged::Status::Unknown, loose: joined, stopped: s.motion.stopped + stopped, ..s.motion },
        left: Some(w),
        left_answer: a.map(Box::new),
        between: Some(gap),
        ..s
    });
    Ok(Check { status: Status::Unknown, accuracy: None, motion, ..sampled })
}

/// Whether every drive whose speed `claim` reads is bounded on moment
/// `t`'s two-sided neighbourhood, as the view asked at `t` alone finds it
/// (`ink ... --t`): where one is not, a value there is no witness.
pub(crate) fn bounded_at(v: &mut ModelView, claim: &Claim, roots: &[u32], t: f64) -> bool {
    let ctx = v.ctx.clone();
    v.at_point(ctx, |v| {
        v.motion = None;
        v.set_moment(t).is_ok() && jumps(v, claim, roots).is_none()
    })
}

/// A derivative's formula applies on either side of a jump. Before using it
/// over time, bound each marked drive's rate where the claim reads it.
/// Some: the claim's unknown answer, and the span a suspected jump was
/// narrowed to, if it was (a rate with no bound at all narrows nothing).
pub(crate) fn jumps(v: &mut ModelView, claim: &Claim, roots: &[u32]) -> Option<(Check, Option<(f64, f64)>)> {
    if !v.doc.g.has_speed_reads() {
        return None;
    }
    let span = v.motion.unwrap_or((v.ctx.t, v.ctx.t));
    let origin = [interval::I::point(0.0); 3];
    let space = [interval::I::ENTIRE; 3];
    let mut todo = Vec::new();
    match claim {
        Claim::Point { f, at, .. } => {
            if !v.doc.g.motion_reads(*f).is_empty() {
                let (mut ctx, _, covers) = claims::judged_scope(v, at, Some(span));
                ctx.t_span = Some(span);
                let b = if covers { at.map(|n| ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, n), &ctx, origin, None)) } else { space };
                todo.push((*f, b));
            }
            todo.extend(roots[1..].iter().map(|&f| (f, origin)));
        }
        Claim::Value(..) => todo.extend(roots.iter().map(|&f| (f, origin))),
        _ => todo.extend(roots.iter().map(|&f| (f, space))),
    }
    for (read, b, integrated) in required_drives(&v.doc, todo) {
        let from = if integrated { v.doc.t_span.map_or(span.0, |s| s.0.min(span.0)) } else { span.0 };
        if let Some((s, evidence, narrowed)) = drive_jump(v, read, b, (from, span.1)) {
            let mut c = checked(Status::Unknown, evidence);
            if let Some(span) = v.motion {
                // A failed input before this question invalidates the state
                // already at its start.
                let s = if s.0 < span.0 { span } else { s };
                let m = Motion { undecided: vec![s], t: Some(ink_core::eval::span_mid(s)), ..Motion::none(ink_kernel::ranged::Status::Unknown, span) };
                c.motion = Some(Stepped { ranged: true, ..Stepped::of(m, true, None) });
            }
            return Some((c, narrowed.then_some(s)));
        }
    }
    None
}
