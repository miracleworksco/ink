use super::*;

/// Whether `roots` read t other than as a provider's own moment: t less a param a solve on t sets is
/// time since the moment it was provided at (a simulated value's local model about it), which that
/// provider's steps follow.
pub(super) fn reads_t_itself(v: &ModelView, roots: &[u32]) -> bool {
    let g = &v.doc.g;
    let moment = |n: u32| matches!(g.nodes[n as usize], Node::Var(Var::Param(k)) if v.doc.solves.iter().any(|s| s.on_t && s.owns.contains(&(k as usize))));
    let (mut seen, mut stack) = (std::collections::HashSet::new(), roots.to_vec());
    while let Some(n) = stack.pop() {
        if !seen.insert(n) {
            continue;
        }
        match g.nodes[n as usize] {
            Node::Var(Var::T) => return true,
            Node::Bin(BinOp::Sub, a, b) if matches!(g.nodes[a as usize], Node::Var(Var::T)) && moment(b) => {}
            node => {
                let (children, len) = node.children();
                stack.extend(&children[..len]);
            }
        }
    }
    false
}

/// What `shown_over` finds of a step.
pub(super) enum Over {
    Shown,
    Not,
    Fails(Box<Check>),
}

/// Whether a claim that holds at the middle of `s`, a step a provider chose (the view at that
/// moment), holds at every moment of it with what the provider set held as it is there. A value
/// claim: its value, less its bound, keeps its side over the step by interval bounds on t. A claim
/// about shapes: asked with t ranged over the step (a search of BETWEEN_BOXES boxes), where what a
/// provider sets about its moment (a simulated coordinate's local motion) moves with t as it says; a
/// witness there is a failure, to what the provider's own account of the step is good for.
pub(super) fn shown_over(v: &mut ModelView, claim: &Claim, bound: Option<u32>, s: (f64, f64), opts: ClaimOptions) -> Over {
    let Claim::Value(f, cmp, limit) = *claim else {
        let k = v.t_param();
        let tk = v.doc.g.var(Var::Param(k));
        let g = &mut v.doc.g;
        let Some(asked) = with_nodes(claim, &mut |n| g.subst(n, &[(Var::T, tk)])) else { return Over::Not };
        let (ctx, motion, boxes) = (v.ctx.clone(), v.motion.take(), std::mem::replace(&mut v.search_boxes, BETWEEN_BOXES));
        v.ctx.params[k as usize] = ink_core::eval::span_mid(s);
        v.ctx.p_spans = v.ctx.p_spans.clone().with_band(k, s);
        let r = claims::decide(v, asked, opts);
        (v.ctx, v.motion, v.search_boxes) = (ctx, motion, boxes);
        return match r {
            Ok(c) if c.status == Status::Holds => Over::Shown,
            Ok(c) if c.status == Status::Fails && c.accuracy.is_some() => Over::Fails(Box::new(Check { accuracy: None, ..c })),
            _ => Over::Not,
        };
    };
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let (node, limit) = match bound {
        Some(b) => (v.doc.g.sub(f, b), 0.0),
        None => (f, limit),
    };
    let over = Ctx { t_span: Some(s), ..v.ctx.clone() };
    let iv = ink_core::eval::eval_interval(&Tape::new(&v.doc.g, node), &over, [interval::I::point(0.0); 3], None);
    match (ink_kernel::ranged::Comparison { cmp, bound: limit, tol }).holds_on(iv) {
        true => Over::Shown,
        false => Over::Not,
    }
}

/// A failure found at a sampled moment (the view at it), a guarantee only
/// where the same claim, built from its own nodes with t the param `ranged`
/// asks over the span, registered at the moment's one value (every other
/// ranged input as it is), fails too: then with that answer's witness. A
/// sample reads t as a number, and a check's fast path may read what its
/// nodes give there in floats ((1e16 + t) - 1e16 is t, 0 in floats; a
/// probe moved by it frozen where it rounds): asked so, its points and
/// values are the original expressions' in real numbers, bounded outward
/// (Frame::ctx's range before its point). Otherwise an estimate. A claim
/// whose nodes read no t froze nothing: its sample is the question.
pub(super) fn certified(v: &mut ModelView, claim: &Claim, bound: Option<u32>, c: Check, opts: ClaimOptions) -> Check {
    let roots = v.doc.claim_motion_inputs(claim, bound);
    if c.accuracy.is_none() || !roots.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)) {
        return c;
    }
    let estimate = |c: Check| Check { accuracy: None, ..c };
    let k = v.t_param();
    let tk = v.doc.g.var(Var::Param(k));
    let g = &mut v.doc.g;
    let mut on = |n: u32| g.subst(n, &[(Var::T, tk)]);
    let Some(asked) = with_nodes(claim, &mut on) else { return estimate(c) };
    let asked = match (asked, bound.map(&mut on)) {
        (Claim::Value(f, cmp, _), Some(b)) => Claim::Value(v.doc.g.sub(f, b), cmp, 0.0),
        (_, Some(_)) => return estimate(c),
        (a, None) => a,
    };
    let (ctx, motion) = (v.ctx.clone(), v.motion.take());
    let t = v.ctx.t;
    v.ctx.params[k as usize] = t;
    v.ctx.p_spans = v.ctx.p_spans.clone().with_band(k, (t, t));
    let r = claims::decide(v, asked, opts);
    (v.ctx, v.motion) = (ctx, motion);
    match r {
        Ok(r) if r.status == Status::Fails && r.accuracy.is_some() => Check { evidence: r.evidence, ..c },
        _ => estimate(c),
    }
}

/// A claim at moment t reading params whose provider reads them over the
/// step holding t within a band (ParamProvider::band): a round reading as
/// the band's centre, the real values what it set. A value claim holds
/// where it holds for every value in the band (interval bounds; its value
/// then the band's worst), fails where it fails at the real values, and is
/// undecided between, saying so; another claim fails where it fails at the
/// real values, else is undecided. A band of no finite width decides
/// nothing by holding. None: no band.
pub(super) fn banded(v: &mut ModelView, claim: &Claim, roots: &[u32], t: f64, opts: ClaimOptions) -> Option<Check> {
    let mut band: Vec<(usize, f64, f64)> = Vec::new();
    for i in v.doc.solves_read(roots) {
        let s = &v.doc.solves[i];
        let Some(p) = v.doc.ext[s.ext].form.provider(&v.doc) else { continue };
        if let Some((w, centre)) = p.band(&v.ctx.params, t) {
            band.extend(s.owns.iter().zip(w.into_iter().zip(centre)).map(|(&k, (w, c))| (k, w, c)));
        }
    }
    let read = v.doc.params_read(roots);
    band.retain(|b| read.contains(&b.0) && b.1 != 0.0);
    if band.is_empty() {
        return None;
    }
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let x = match claim {
        Claim::Value(f, cmp, bound) => {
            let comparison = ink_kernel::ranged::Comparison { cmp: *cmp, bound: *bound, tol };
            let tape = Tape::new(&v.doc.g, *f);
            if band.iter().all(|b| b.1.is_finite()) {
                let spans: Vec<(u32, (f64, f64))> = band.iter().map(|&(k, w, c)| (k as u32, (c - w, c + w))).collect();
                let mut centre = v.ctx.clone();
                for &(k, _, c) in &band {
                    centre.params[k] = c;
                }
                let over = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&spans), ..centre };
                let iv = ink_core::eval::eval_interval(&tape, &over, [interval::I::point(0.0); 3], None);
                if comparison.holds_on(iv) {
                    let worst = if matches!(cmp, crate::Cmp::Lt | crate::Cmp::Le) { iv.hi } else { iv.lo };
                    return Some(Check { accuracy: None, ..checked(Status::Holds, Evidence::Value(worst)) });
                }
            }
            let x = ink_core::eval::eval_point(&tape, &v.ctx, [0.0; 3]);
            if comparison.fails_on(interval::I::point(x)) || x.is_nan() {
                return Some(Check { accuracy: None, ..checked(Status::Fails, Evidence::Value(x)) });
            }
            Some(x)
        }
        _ => {
            let at = v.doc.claim_at(claim, None, &v.ctx);
            match claims::decide(v, at, opts) {
                Ok(c) if c.status == Status::Fails => return Some(c),
                _ => None,
            }
        }
    };
    let widths: Vec<String> = band.iter().map(|&(k, w, _)| format!("{} within {w:.6}", v.doc.params[k].name)).collect();
    let why = format!(
        "what it reads is set at the moment's real value, and over its step read with parts that spin at rest, round, the real value within a band of \
         that ({}, in their units: motion.basis): neither holding over that band nor failing at this moment's real value{}",
        widths.join(", "),
        x.map_or(String::new(), |x| format!(" ({x:.6})"))
    );
    Some(checked(Status::Unknown, Evidence::Unsolved(InkError::msg(why))))
}

/// The stepped answer: its status from the steps, with the evidence of the
/// step that decided it (a disproof, the first undecided, else the
/// tightest proven). Measured between steps (not `bound`), it is an
/// estimate, except a failure: the failing step's own answer, a witness.
pub(super) fn finish(deciding: Option<Check>, m: Motion, bound: bool, step: Option<f64>) -> Check {
    let status = Status::of(m.status);
    let c = deciding.unwrap_or_else(|| checked(status, Evidence::Stepped));
    let accuracy = c.accuracy.filter(|_| bound || status == Status::Fails);
    Check { status, accuracy, in_window: c.in_window.filter(|_| status == Status::Unknown), motion: Some(Stepped::of(m, bound, step)), ..c }
}

/// The rigid group item `it` rides over the motion `span`: how many of its
/// links, from the top, reach the last whose placement moves (0: none
/// does, the world's). Items of one group (the same instance at that link)
/// keep their placements relative to each other at every moment: nothing
/// below it moves, and every placement down to it is proven an isometry
/// over the span (rigid::Isometry: a rotation exactly, its translation
/// finite), so what they share moves them both rigidly. Sharing a moving link is not
/// enough: a link whose matrix scales with t (a turn about an axis given
/// too small to normalize exactly) stretches the distance between its
/// parts as it turns, and one placement above it that is no rotation
/// stretches the turning distances by direction, and one whose
/// translation is infinite at a moment (a move by ln((t - 0.5)^2)) places
/// nothing there. None: its own field moves (it deforms), or a placement
/// down to its group is not so proven.
pub fn rigid_group(v: &ModelView, it: &crate::Item, span: (f64, f64)) -> Option<usize> {
    if claims::moves(v, &[it.shape], span) {
        return None;
    }
    let moves = |l: &crate::Link| claims::moves(v, &l.at.r.iter().flatten().chain(&l.at.t).copied().collect::<Vec<u32>>(), span);
    let k = it.links.iter().rposition(moves).map_or(0, |i| i + 1);
    let over = Ctx { t_span: Some(span), ..v.ctx.clone() };
    it.links[..k].iter().all(|l| Isometry::over(&v.doc.g, &l.at.r, &l.at.t, &over).is_some()).then_some(k)
}
