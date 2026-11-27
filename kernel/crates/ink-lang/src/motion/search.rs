use super::*;

/// What `ranged` settled: an answer over every moment, or why the sampled
/// one stands.
pub(super) enum Ranged {
    Decided(Check),
    /// Why, and the search's own unknown answer where it searched.
    Left(String, Option<Check>),
    /// A disjoint pair's emptiness to search beside others, whether its
    /// shapes are proven bounded, the ranged params it reads, the frame it
    /// is searched in (`Frame`), whether t's param is a moment (not a
    /// drive's angle: `Drive`, where only a hold is an answer), and what
    /// the shapes asked may be off the world's (Held::charged).
    Wait(Search, bool, Vec<(u32, (f64, f64))>, Option<Frame>, bool, f64),
}

/// The answer with t ranged as param k: holding or failing, the claim's
/// over every moment (its witness's moment, param k's value); undecided,
/// the sampled answer stands (`own`: the time share stopped it).
pub(super) fn decided(c: Check, span: (f64, f64), own: bool) -> Ranged {
    match c.status {
        Status::Unknown => Ranged::Left(
            match own {
                true => "with t ranged over the span the search was not decided within half the call's time left".into(),
                false => "with t ranged over the span the search left it undecided".into(),
            },
            Some(c),
        ),
        status => {
            let k = c.ranged.iter().map(|r| r.0).max();
            let t = k.and_then(|k| witness_t(&c.evidence, k));
            let m = Motion { t, ..Motion::none(if status == Status::Holds { ink_kernel::ranged::Status::Proven } else { ink_kernel::ranged::Status::Disproven }, span) };
            Ranged::Decided(Check { motion: Some(Stepped { ranged: true, ..Stepped::of(m, true, None) }), ..c })
        }
    }
}

/// The claim with t a ranged param over `span` (the view's param "t",
/// made once): every node of it, and its bound, with t replaced by that
/// param (Graph::subst), decided as any claim over ranged params is
/// (claims::decide), the param's value a witness's moment. Only t's span is
/// added to the view's scope: what the view knows of fields that do not read it
/// holds. What a solve on t sets is derived from that param
/// (solved_over_t). Not ranged: a claim of a kind no band ranges (a
/// package's metric, an integral, a thickness, cavities), one reading a
/// solve on t not followed over the span, one whose bound moves.
/// Asked alone (require), it stops at half the call's time left (stop::within);
/// beside other searches (interfere), at the call's limit.
pub(super) fn ranged(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Ranged, InkError> {
    let whole = ranged_whole(v, claim, bound, span, opts)?;
    if matches!(whole, Ranged::Decided(_)) {
        return Ok(whole);
    }
    if !opts.in_pieces {
        return Ok(whole);
    }
    // Left so: of the part of a shape that repeats about its own axis, and of what is added to it,
    // each asked alone. Every one clear, the shape is (it is within them: `repeating_parts`). Once,
    // where the claim's period is asked in pieces (`between`): its parts have the same one.
    let Some(parts) = repeating_parts(v, claim) else { return Ok(whole) };
    let mut held: Option<Check> = None;
    for part in &parts {
        match ranged_whole(v, part, bound, span, opts)? {
            Ranged::Decided(c) if c.status == Status::Holds => held = Some(c),
            _ => return Ok(whole),
        }
    }
    Ok(match held {
        Some(c) => Ranged::Decided(Check { held: c.held.map(|_| claims::Held::ToAccuracy { near: None }), note: Some(PARTS), ..c }),
        None => whole,
    })
}

const PARTS: &str = "asked of the part of a shape that repeats about its own axis (a wheel's teeth and flanges, without its bore's flat or its set screw's hole), and of what is added to it: the shape lies within them, so clear of them is clear of it";

/// A claim that two shapes are apart (disjoint, a clearance) as the same claims of one shape's
/// parts, where that shape is placed and repeats about its own z but for small parts
/// (symmetry::nearly_about): of its core, which repeats, and of what is added to it. What is cut
/// from the shape is left in, so the shape lies within the two and clear of them is clear of it
/// (not the other way: a claim these fail is not shown to fail). The core repeats every tooth where
/// the whole wheel, with a flat on its bore, repeats only every turn. None: neither shape is so.
fn repeating_parts(v: &mut ModelView, claim: &Claim) -> Option<Vec<Claim>> {
    use ink_core::symmetry::{Order, nearly_about, turn_symmetry};
    let (a, b, gap) = match *claim {
        Claim::Disjoint(a, b) => (a, b, None),
        Claim::Clearance(a, b, d) => (a, b, Some(d)),
        _ => return None,
    };
    let again = |a: u32, b: u32| match gap {
        None => Claim::Disjoint(a, b),
        Some(d) => Claim::Clearance(a, b, d),
    };
    for (side, other, first) in [(b, a, false), (a, b, true)] {
        let reached: std::collections::HashSet<u32> = v.doc.g.reachable(&[side]).into_iter().collect();
        let placed: Vec<(u32, crate::Placed)> = v.doc.placements().into_iter().filter(|(w, _)| reached.contains(w)).collect();
        // The one shape placed: a named one is placed again without its name (Doc::placements),
        // inside itself.
        let outer: Vec<(u32, crate::Placed)> = placed.iter().filter(|(w, _)| !placed.iter().any(|(o, _)| o != w && v.doc.g.reachable(&[*o]).contains(w))).copied().collect();
        let [(w, p)] = outer[..] else { continue };
        if turn_symmetry(&v.doc.g, p.shape).order != Order::Of(1) {
            continue;
        }
        let Some(near) = nearly_about(&mut v.doc.g, p.shape, 2) else { continue };
        if near.at != [0.0; 3] || near.symmetry.order == Order::Of(1) {
            continue;
        }
        let of = |v: &mut ModelView, part: u32| {
            let crate::Value::Field(world, _) = v.doc.to_world(&crate::Value::Field(part, None), &p.at) else { return side };
            v.doc.g.replace(side, &[(w, world)].into())
        };
        let mut parts = vec![of(v, near.core)];
        parts.extend(near.added.map(|p| of(v, p)));
        return Some(parts.into_iter().map(|part| if first { again(part, other) } else { again(other, part) }).collect());
    }
    None
}

/// `ranged` of the claim as it is.
pub(super) fn ranged_whole(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Ranged, InkError> {
    // Through its drive where it has one (`Drive`); not decided there
    // but by a hold, asked over t. A witness in the frame that is none in
    // the world: asked there.
    let (mut framed, mut through) = (true, true);
    loop {
        match ranged_in(v, claim, bound, span, opts, framed, through)? {
            Some(r) => return Ok(r),
            None if through => through = false,
            None if framed => framed = false,
            None => unreachable!("a world search over t needs no confirming"),
        }
    }
}

/// `ranged`, in the first shape's frame where `framed` and it has one
/// (`moving_frame`), through its drive where `through` and it has one
/// (`drive`). None: decided there by a witness that is not shown one in
/// the world (`Frame::confirmed`), or through the drive other than by a
/// hold.
#[allow(clippy::too_many_arguments)]
fn ranged_in(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions, framed: bool, through: bool) -> Result<Option<Ranged>, InkError> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    // The whole scope first (`scope_over`), kept until the answer is given
    // in the world: every rewrite of the claim (its frame, simplified over
    // the scope) and every witness confirmed reads a solved input as
    // ranged, never as the number it is at the view's moment.
    let (scope, k) = match scope_over(v, &roots, span) {
        Ok(s) => s,
        Err(why) => return Ok(Some(Ranged::Left(why, None))),
    };
    let (ctx, motion) = (std::mem::replace(&mut v.ctx, scope), v.motion.take());
    let r = scoped(v, claim, bound, span, opts, (framed, through), k, &roots);
    (v.ctx, v.motion) = (ctx, motion);
    r
}

/// A claim's scope over a motion (`roots` its fields): the view's context
/// with t's param k (made once: ModelView::t_param) ranged over `span`, and
/// nothing on t a form sets read (solved_over_t).
/// It is what the claim is validated in (claims::ask_in), searched in and
/// confirmed in, in that order. Err: why a solve on t it reads is not
/// derived (it is then measured, moment by moment).
pub(crate) fn scope_over(v: &mut ModelView, roots: &[u32], span: (f64, f64)) -> Result<(Ctx, u32), String> {
    let k = v.t_param();
    solved_over_t(v, roots, span)?;
    let mut c = v.ctx.clone();
    c.params[k as usize] = ink_core::eval::span_mid(span);
    c.p_spans = c.p_spans.with_band(k, span);
    c.derive();
    Ok((c, k))
}

/// `ranged_in` in its scope (the view's context, t's param k ranged over
/// `span` and what solves on t set derived from it): in its frame, then
/// through its drive, then searched, then confirmed in the world.
#[allow(clippy::too_many_arguments)]
fn scoped(
    v: &mut ModelView,
    claim: &Claim,
    bound: Option<u32>,
    span: (f64, f64),
    opts: ClaimOptions,
    (framed, through): (bool, bool),
    k: u32,
    roots: &[u32],
) -> Result<Option<Ranged>, InkError> {
    if let Some(why) = sweeps_far(v, roots, span) {
        return Ok(Some(Ranged::Left(why, None)));
    }
    let t = v.doc.g.var(Var::Param(k));
    // Asked in pieces only of a period: where no drive shows one there is nothing more to ask.
    if opts.in_pieces && !through {
        return Ok(Some(Ranged::Left("its motion shows no period to ask in pieces".into(), None)));
    }
    // Searched where one of the two shapes rests (`frame`): the same claim
    // at every moment, its points moved rigidly.
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let frame = if framed { moving_frame(v, claim, span, Some(tol)) } else { None };
    let charge = frame.as_ref().map_or(0.0, |f| f.2);
    let tn = v.doc.g.var(Var::T);
    // Through the drive all its reads of t pass through, where one is
    // shown periodic and the claim in its frame still moves: t's param is
    // then the drive's value.
    let drive = match through {
        true => {
            let moves = roots.iter().any(|&n| {
                let n = frame.as_ref().map_or(n, |m| v.doc.g.replace(n, &m.0));
                v.doc.g.depends_on(n, |x| x == Var::T)
            });
            moves.then(|| drive(v, roots, span, k)).flatten()
        }
        false => None,
    };
    if through && drive.is_none() {
        return Ok(None);
    }
    let by: std::collections::HashMap<u32, u32> = [(drive.as_ref().map_or(tn, |d| d.node), t)].into();
    // Its witnesses are moments where t is the drive and nothing is
    // charged (`Drive::slack`); else only a hold is an answer.
    let moments = drive.as_ref().is_none_or(|d| d.node == tn && d.slack == 0.0);
    let world = frame.as_ref().filter(|_| moments).map(|f| Frame::new(v, claim, f.1, k, span, opts.accuracy, charge));
    let mut on = |n: u32| {
        let n = match &frame {
            Some(m) => v.doc.g.replace(n, &m.0),
            None => n,
        };
        v.doc.g.replace(n, &by)
    };
    let Some(asked) = with_nodes(claim, &mut on) else {
        return Ok(Some(Ranged::Left("a claim of this kind is not ranged".into(), None)));
    };
    let bound = bound.map(&mut on);
    let asked = match (asked, bound) {
        // A bound that moves: the value checked is the expression minus
        // the bound, against 0.
        (Claim::Value(f, cmp, _), Some(b)) => Claim::Value(v.doc.g.sub(f, b), cmp, 0.0),
        (_, Some(_)) => return Ok(Some(Ranged::Left("its bound moves (reads t), and is not ranged".into(), None))),
        (c, None) => c,
    };
    // What rides together, cancelled: a sum in which the motion drops out in real numbers is the sum
    // left (rigid::cancelled), the same field.
    // Only where the claim is seen from a frame that moves: in the world nothing rides.
    let asked = match frame.is_some() {
        false => asked,
        true => {
            let read = v.doc.claim_roots(&asked);
            let same = ink_kernel::rigid::cancelled(&mut v.doc.g, &read, &v.ctx, &|x| x == Var::T || x == Var::Param(k));
            match same.is_empty() {
                true => asked,
                false => with_nodes(&asked, &mut |n| v.doc.g.replace_equal(n, &same)).expect("a claim of its kind"),
            }
        }
    };
    let read = v.doc.claim_roots(&asked);
    if drive.is_some() && v.doc.g.reachable(&read).contains(&tn) {
        return Ok(None);
    }
    // What the shapes asked may be off the world's, charged to the claim:
    // the drive's slack (Drive::slack) and the pivots' drift
    // (`moving_frame`). An emptiness decided to its accuracy less it, a
    // clearance over its distance plus it, an overlap's witness deeper by
    // it at every moment (not through a drive: its witness is no moment).
    let slack = drive.as_ref().map_or(0.0, |d| d.slack);
    let charged = interval::add_up(charge, slack);
    let (asked, opts) = match asked {
        _ if charged == 0.0 => (asked, opts),
        // The outer region of a bridge needs exact emptiness in world coordinates.
        Claim::EmptyBounds { .. } => return Ok(None),
        Claim::Empty(..) | Claim::Disjoint(..) | Claim::Subset(..) => {
            let less = -interval::add_up(-tol, charged);
            if !(less > 0.0) {
                return Ok(None);
            }
            (asked, ClaimOptions { accuracy: Some(less), ..opts })
        }
        Claim::Clearance(a, b, gap) => (Claim::Clearance(a, b, interval::add_up(gap, charged)), opts),
        Claim::Overlaps(..) if slack == 0.0 => (asked, ClaimOptions { accuracy: Some(interval::add_up(tol, charged)), ..opts }),
        _ => return Ok(None),
    };
    match &drive {
        // One period of the drive's values stands for every moment.
        Some(d) => {
            v.ctx.p_spans.narrow(k, d.window);
            v.ctx.params[k as usize] = ink_core::eval::span_mid(d.window);
        }
        None => {
            // Over one period of a motion longer than one
            // (ink_core::period): the rest repeats it.
            let tapes: Vec<Tape> = read.iter().map(|&n| Tape::new(&v.doc.g, n)).collect();
            if let Some(one) = ink_core::period::quotient(&tapes.iter().collect::<Vec<_>>(), &v.ctx) {
                v.ctx = one;
            }
        }
    }
    // t's param is reported over the span (through the drive its window
    // is the drive's values, not moments).
    let at_span = |mut r: Vec<(u32, (f64, f64))>| {
        r.iter_mut().filter(|x| x.0 == k && !moments).for_each(|x| x.1 = span);
        r
    };
    match (opts.beside, &asked) {
        (true, Claim::Disjoint(a, b)) => {
            let tol = claims::claim_accuracy(v, &asked, opts.accuracy);
            let f = v.doc.g.max(*a, *b);
            let tape = v.clipped(f);
            let crate::Reach { lo, hi, whole } = v.reach(f);
            Ok(Some(Ranged::Wait(Search(tape, v.ctx.clone(), (lo, hi), v.proof_options(opts.res, tol)), whole, at_span(v.ranged_reads(&read)), world, moments, charged)))
        }
        _ => {
            let (c, own) = ink_kernel::stop::within(if opts.beside { 1.0 } else { 0.5 }, || match (&drive, opts.in_pieces) {
                (Some(d), true) => in_pieces(v, &asked, opts, k, d.window),
                _ => claims::decide(v, asked, opts),
            });
            // t's param among them where the world reads it, though the
            // frame cancels it (a witness, a near point, is at a moment).
            let mut ranged = at_span(v.ranged_reads(&read));
            if (world.is_some() || !moments) && !ranged.iter().any(|r| r.0 == k) {
                ranged.push((k, span));
            }
            let c = c?;
            // An overlap fails only proven apart outright, which shapes off
            // the world's do not show of it: asked there. The answer is to
            // the claim's own accuracy, and holds outright only of shapes
            // asked as they are (Held::charged).
            if charged > 0.0 && matches!(claim, Claim::Overlaps(..)) && c.status == Status::Fails {
                return Ok(None);
            }
            let c = Check { accuracy: c.accuracy.map(|_| tol), held: c.held.map(|h| h.charged(charged)), ..c };
            // Through the drive, a hold is the answer (its near point is
            // at a drive's value, not a moment: none given); else asked
            // over t.
            let c = match (moments, c.status) {
                (true, _) => c,
                (false, Status::Holds) => Check {
                    held: c.held.map(|h| match h {
                        claims::Held::Outright => claims::Held::Outright,
                        claims::Held::ToAccuracy { .. } => claims::Held::ToAccuracy { near: None },
                        h @ claims::Held::SetAside { .. } => h,
                    }),
                    note: c.note.or(Some(DRIVE)),
                    ..c
                },
                (false, Status::Fails) => return Ok(None),
                (false, _) => c,
            };
            // Its witnesses in the world, each shown one there (undecided,
            // there are none).
            let c = match &world {
                Some(w) if c.status != Status::Unknown => match w.confirmed(v, c.evidence.clone()) {
                    Some(evidence) => Check { evidence, held: c.held.clone().map(|h| w.held(v, h)), ..c },
                    None => return Ok(None),
                },
                _ => c,
            };
            Ok(Some(decided(Check { ranged, ..c }, span, own)))
        }
    }
}

/// Whether what `roots` read is bounded over a span of t: Err where one
/// reads a param a form on t sets, which is measured at each moment (an
/// effort at its steps, a simulation's state) and has no bound over a
/// span, with which; one that rests on a ranged param says that first.
fn solved_over_t(v: &ModelView, roots: &[u32], span: (f64, f64)) -> Result<(), String> {
    let follow = v.doc.follow();
    // A condition on t that an estimate applies under (estimate-when) is shown over the whole span or
    // the claim is not ranged over it: its marker is then nothing over the span, as at a moment.
    let mut guarded = Ok(());
    let unbounded = |i: &usize| {
        let s = &v.doc.solves[*i];
        let Some(guard) = crate::estimate::Estimate::of(&v.doc.ext[s.ext]).filter(|_| s.on_t) else { return s.on_t };
        let mut over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..v.ctx.clone() };
        if let (Err(why), true) = (guard.provide_over(&mut over), guarded.is_ok()) {
            guarded = Err(crate::estimate::not_over(&why, span));
        }
        false
    };
    let found = v.doc.solves_read(roots).into_iter().find(unbounded);
    let Some(i) = found else { return guarded };
    if let Some(&r) = follow.rests_on(i).iter().find(|&&r| v.ctx.p_spans.get(r as u32).is_some()) {
        return Err(format!("it reads {}, which rests on ranged '{}' and is followed along t at one value of it", follow.followed(i).2, v.doc.params[r].name));
    }
    Err(format!("it reads {}, which rests on a value measured at moments along t (an effort at its steps, a simulation's state), not bounded", follow.followed(i).2))
}

/// How a search beside others (`Search`'s tape and context: `asked`)
/// held, as `shallow` says it: its near point found in the scope, moved
/// to the world where it was asked in a frame. None: not proven.
pub(super) fn held_beside(v: &ModelView, p: &Found, (tape, ctx): &(std::sync::Arc<Tape>, Ctx), world: Option<&Frame>) -> Option<claims::Held> {
    if p.status() != ink_kernel::ranged::Status::Proven {
        return None;
    }
    let points = if p.exact { Vec::new() } else { claims::Held::scoped(p, tape, ctx) };
    let h = claims::Held::confirmed_at(p.exact, points, tape, ctx);
    Some(match world {
        Some(w) => w.held(v, h),
        None => h,
    })
}
