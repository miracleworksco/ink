use super::*;

/// Check one claim (`bound`: its bound, as a field, if it reads params). A
/// claim that reads inexact inputs (Doc::inexact: sampled fields, a solved
/// value not certified as a root) is an estimate: whatever proves it proves
/// it about those samples or numbers. So is one checked by a measuring
/// method (rays, a package's solve), whatever it reads.
///
/// Every claim is checked here, at a moment, over ranged params or over a
/// motion, and here alone is its accuracy settled: it is asked once, to
/// its own (a shape no thicker than twice it holds to it: nothing violates
/// it by more). The answer's accuracy is the one asked, where its method
/// guarantees it. The view is scoped as every claim's question first
/// (ModelView::scope_question).
pub fn check(v: &mut ModelView, claim: Claim, bound: Option<u32>, opts: ClaimOptions) -> Result<Check, InkError> {
    ask(v, claim, bound, opts)?.wait(v)
}

/// `check`, its searches over a motion left to the caller to run
/// (motion::Asked).
pub fn ask(v: &mut ModelView, claim: Claim, bound: Option<u32>, opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    ask_bounded(v, claim, bound, bound, opts)
}

/// `ask`, with the bound as written (`original_bound`) beside the one
/// checked: an assert whose ranged bound moved into its claim (`asserted`)
/// is asked here too, so no claim skips the support-depth rule.
pub(super) fn ask_bounded(v: &mut ModelView, claim: Claim, bound: Option<u32>, original_bound: Option<u32>, opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    v.scope_question()?;
    // A thin skin can be unsupported over its entire width. Whether the
    // region is asked empty directly, clipped, or compared with another
    // shape, only exact emptiness proves it. Accuracy still sets refinement.
    let shape_test = matches!(
        claim,
        Claim::Empty(_) | Claim::EmptyBounds { .. } | Claim::Disjoint(..) | Claim::Subset(..) | Claim::Overlaps(..) | Claim::Touches(..) | Claim::Clearance(_, _, 0.0)
    );
    if shape_test
        && !v.doc.support_depths.is_empty()
        && v.doc.claim_roots(&claim).iter().any(|&f| v.doc.g.depends_on(f, |var| matches!(var, Var::Param(k) if v.doc.support_depths.contains(&k))))
    {
        let tol = claim_accuracy(v, &claim, opts.accuracy);
        let exact = ClaimOptions { accuracy: Some(0.0), res: opts.res.or((tol > 0.0).then_some(tol / 2.0)), ..opts };
        return ask_in(v, claim, bound, original_bound, exact)?.then(v, |_, mut c| {
            // The depth of a cell in the search mask says nothing about how
            // far a physical roof projects. Preserve the witness or cells.
            if let Evidence::Empty(Proof { verdict: Verdict::Unknown { attained, .. }, .. }) = &mut c.evidence {
                *attained = None;
            }
            Ok(Asked::Done(c))
        });
    }
    ask_in(v, claim, bound, original_bound, opts)
}

/// A ranged bound may already be subtracted into the value that will be
/// checked. Keep the original bound among the inputs to validate even when
/// that subtraction cancels (a speed compared to itself at a jump).
fn ask_in(v: &mut ModelView, claim: Claim, bound: Option<u32>, original_bound: Option<u32>, opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    // A real construction on the right of a comparison needs the same
    // enclosure as one on the left. Preserve it through motion checks too.
    let real_bound = bound.filter(|&b| !v.doc.g.fixed_scalar(b));
    let (claim, bound, compared) = match (&claim, real_bound) {
        (Claim::Value(f, cmp, _), Some(b)) => (Claim::Value(v.doc.g.bin_real(ink_core::graph::BinOp::Sub, *f, b), *cmp, 0.0), None, true),
        _ => (claim, bound, false),
    };
    let mut roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(original_bound).collect();
    // A cancelled speed can still depend on an estimated provider. Preserve
    // its accuracy and validity scope when the numerical root is a literal.
    let mut parameters: Vec<_> = roots
        .iter()
        .filter(|&&f| crate::motion::needs_speed(&v.doc, f))
        .flat_map(|&f| v.doc.g.motion_reads(f))
        .filter_map(|r| match *r {
            ink_core::graph::MotionRead::Param(k) => v.doc.g.var_node(Var::Param(k)),
            _ => None,
        })
        .collect();
    parameters.sort_unstable();
    parameters.dedup();
    roots.extend(parameters);
    // Some package claims run their own model. Demand published masses
    // in this view too, before snapshotting their bounds for validation.
    let mut masses = Vec::new();
    for i in v.doc.solves_read(&roots) {
        let solve = &v.doc.solves[i];
        if v.doc.ext[solve.ext].package == "mass" {
            masses.push(v.doc.g.var(Var::Param(solve.owns[0] as u32)));
        }
    }
    if !masses.is_empty()
        && let Err(e) = v.want(&masses)
    {
        return Ok(Asked::Done(checked(Status::Unknown, Evidence::Unsolved(e))));
    }
    let tol = claim_accuracy(v, &claim, opts.accuracy);
    let inexact = v.doc.inexact(&roots, &v.certified());
    // A constant-looking result can still read an undefined speed (max of
    // a limit and a step's zero branch derivative). Check the drive before
    // deciding that a single moment or a rigid frame answers the claim.
    let inputs = v.doc.claim_motion_inputs(&claim, original_bound);
    let jump = crate::motion::jumps(v, &claim, &inputs);
    let motion = v.motion.filter(|&s| moves(v, &inputs, s));
    let moment = v.ctx.t;
    let opts = ClaimOptions { accuracy: Some(tol), ..opts };
    if jump.is_none() && motion.is_none() && real_bound.is_none() && planar(v, &claim, &roots) {
        return v.in_plane(|v| {
            let (scope, span, _) = judged_scope(v, &roots, v.motion);
            let c = match validate(v, &claim, &roots, tol, opts.res, scope, span) {
                Some(c) => c,
                None => decide(v, claim.clone(), opts)?,
            };
            Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact)))
        });
    }
    // Judged before anything answers it (`validate`), over the scope it is
    // asked over (`judged_scope`).
    let (scope, span, _) = judged_scope(v, &roots, v.motion);
    if !compared && let Some(b) = real_bound {
        let iv = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, b), &scope, [I::point(0.0); 3], None);
        if iv != I::point(v.doc.scalar(b, &v.ctx)) {
            let c = Check { note: Some("the comparison threshold is known only within rounding bounds"), ..checked(Status::Unknown, Evidence::Unexamined) };
            return Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact)));
        }
    }
    if let Some(c) = validate(v, &claim, &roots, tol, opts.res, scope, span) {
        return Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact)));
    }
    // Where a drive may jump, its continuity is unproved, so the claim
    // cannot hold over the span; but a witness that breaks it at a moment
    // needs no continuity. Search as usual and keep only a failure: any
    // other answer is the jump's unknown. A failure is a witness only at a
    // moment where the formula's rate is the drive's: not inside a span a
    // suspected jump was narrowed to, and, unless it is measured (an
    // estimate, which states its own error), at a moment where every
    // drive it reads is bounded (motion::bounded_at). Else unknown.
    let c = match motion {
        // Over a motion, how it held is said where t ranged decided it
        // (its near point at a moment), never of a moment measured.
        Some(span) => crate::motion::over_span(v, claim.clone(), bound, span, opts)?.then(v, |_, c| {
            let ranged = c.motion.as_ref().is_some_and(|m| m.ranged);
            Ok(Asked::Done(Check { held: c.held.clone().filter(|_| ranged), ..c }))
        })?,
        None => Asked::Done(decide(v, claim.clone(), opts)?),
    };
    c.then(v, move |v, mut c| {
        if compared {
            c.note = c.note.or(Some("the value checked is the expression minus its enclosed threshold, against zero"));
        }
        let Some((unknown, narrowed)) = jump else {
            return Ok(Asked::Done(finished(v, &claim, &roots, c, tol, inexact)));
        };
        let at = witness_moment(&c, moment);
        if c.status == Status::Fails && !narrowed.is_some_and(|(a, b)| at.is_none_or(|t| a <= t && t <= b)) {
            let done = finished(v, &claim, &roots, c, tol, inexact.clone());
            if done.inexact.any() || at.is_some_and(|t| crate::motion::bounded_at(v, &claim, &inputs, t)) {
                return Ok(Asked::Done(done));
            }
        }
        Ok(Asked::Done(finished(v, &claim, &roots, unknown, tol, inexact)))
    })
}

/// The moment a failure's witness is at: the one its search over t found,
/// else the view's (`moment`), where it was decided at one; None where it
/// was asked once for every moment.
fn witness_moment(c: &Check, moment: f64) -> Option<f64> {
    let at = match &c.evidence {
        Evidence::ValueSpan(ValueSpan { witness: Some((at, _)), .. }) => Some(at.t),
        _ => None,
    };
    match &c.motion {
        Some(m) => m.motion.t.or(at.filter(|_| m.ranged)),
        None => Some(at.unwrap_or(moment)),
    }
}

/// Whether `claim` (reading `roots`) is asked in the plane: a question the
/// same at every height, of shapes that read no z (2D profiles, prisms),
/// is a search over x and y alone (ModelView::in_plane).
fn planar(v: &ModelView, claim: &Claim, roots: &[u32]) -> bool {
    let same = match claim {
        Claim::Empty(_)
        | Claim::EmptyBounds { .. }
        | Claim::Range(_)
        | Claim::Clearance(..)
        | Claim::Disjoint(..)
        | Claim::Subset(..)
        | Claim::Touches(..)
        | Claim::Overlaps(..) => true,
        Claim::Metric { metric, .. } => metric.shape.is_some() && v.doc.registry.plane_metric(metric.metric),
        _ => false,
    };
    same && roots.iter().all(|&r| !v.doc.g.value_depends_on(r, |x| x == Var::Z))
}

/// The one judgement made of a claim before any path answers it: ink's
/// range (validity::outside) and the directions it reads resting on a
/// solve on t (Doc::check_conditions), over `scope` (`judged_scope`'s, or
/// one moment measured: motion::validate_at), over `motion` where it spans
/// one; and that it reads no samples a claim is not decided on
/// (grid::Source::undecided: a result whose band its form did not measure,
/// as a metric's claim without its band is unknown). Some: unknown, saying
/// why, the claim's answer.
pub(crate) fn validate(v: &mut ModelView, claim: &Claim, roots: &[u32], tol: f64, res: Option<f64>, scope: Ctx, motion: Option<(f64, f64)>) -> Option<Check> {
    if let Err(e) = v.doc.check_conditions(&scope, motion, &v.unsolved_solves(), false, Some(roots)) {
        return Some(Check { accuracy: None, ..checked(Status::Unknown, Evidence::Unsolved(e)) });
    }
    // A mass command may show an unsettled estimate. A verdict needs the
    // measured density to be plausible throughout its enclosure.
    for i in v.doc.solves_read(roots) {
        let solve = &v.doc.solves[i];
        if v.doc.ext[solve.ext].package != "mass" {
            continue;
        }
        let k = solve.owns[0];
        let field = v.doc.g.var(Var::Param(k as u32));
        // Over a mass's tolerance its density is derived over it (its band in the scope, ModelView's spans).
        let bounds = match scope.p_spans.get(k as u32) {
            Some((lo, hi)) => Some(I::new(lo, hi)),
            None => v.doc.provided_bounds(field, &scope),
        };
        if !bounds.is_some_and(|b| crate::materials::DENSITY.covers(b)) {
            let range = bounds.map_or_else(|| "unavailable".into(), |b| format!("{}..{}", b.lo, b.hi));
            let name = &v.doc.params[k].name;
            // The value itself in range, only its bounds reaching out: no mistake of the mass's unit, the
            // shape's volume is measured too loosely to show it (a thread, a fine lattice).
            let advice = if scope.p_spans.get(k as u32).is_some_and(|(lo, _)| lo <= 0.0) {
                "Over the published mass's tolerance its low end leaves these bodies nothing: the bodies with a density of their own may weigh as much as it. \
                 Narrow the tolerance, or check the known bodies' densities"
            } else if crate::materials::DENSITY.covers(I::point(scope.params[k])) {
                "The density itself is in range: the mass is plausible, and ink's bounds on this shape's volume (a thread, fine features) are too loose to show it. \
                 A claim is decided only on bounds inside the range, so this one is not: no unit is wrong"
            } else {
                "Check the volume and the published mass unit"
            };
            let e = InkError::msg(format!("{name} = {} g/mm³ has measurement bounds {range}; expected {}. {advice}", scope.params[k], crate::materials::DENSITY.expected()));
            return Some(Check { accuracy: None, ..checked(Status::Unknown, Evidence::Unsolved(e)) });
        }
    }
    if let Some((what, why)) = v.doc.inexact(roots, &[]).samples.iter().find_map(|s| Some((s.what.clone(), s.undecided.clone()?))) {
        let e = InkError::msg(format!("it reads {what}, on which no claim is decided: {why}"));
        return Some(Check { accuracy: None, ..checked(Status::Unknown, Evidence::Unsolved(e)) });
    }
    let why = crate::validity::outside(v, claim, roots, tol, res, scope)?;
    Some(Check { accuracy: None, ..checked(Status::Unknown, Evidence::Outside(why)) })
}

/// The scope a claim reading `roots` is judged over before it is asked
/// over `motion` (`validate`), the motion it spans, and whether it covers
/// every moment. Over a motion, where it reads a solve on t, with what
/// those set ranged (motion::scope_over); where one is not derived (only
/// measured), the view's moment alone, not covering the span: each moment
/// measured is judged then (motion::validate_at), and nothing answers the
/// claim for the whole span at once. Else the view's
/// (validity::scope).
pub(crate) fn judged_scope(v: &mut ModelView, roots: &[u32], motion: Option<(f64, f64)>) -> (Ctx, Option<(f64, f64)>, bool) {
    match motion.filter(|_| v.doc.solves_read(roots).iter().any(|&i| v.doc.solves[i].on_t)) {
        Some(span) => match crate::motion::scope_over(v, roots, span) {
            Ok((c, _)) => (Ctx { t_span: Some(span), ..c }, Some(span), true),
            Err(_) => (Ctx { t_span: None, ..v.ctx.clone() }, None, false),
        },
        None => (crate::validity::scope(v), motion, true),
    }
}

/// `ask_in`'s answer, as it is decided: `answer`, with the ranged params
/// it reads, t among them where it was ranged over a motion
/// (motion::ranged).
fn finished(v: &ModelView, claim: &Claim, roots: &[u32], c: Check, tol: f64, inexact: crate::Inexact) -> Check {
    let mut ranged = ranged_reads(v, roots);
    for &r in &c.ranged {
        if !ranged.contains(&r) {
            ranged.push(r);
        }
    }
    answer(c, tol, inexact, v.doc.accuracy_quantity(claim), ranged)
}

/// A claim's answer over a motion decided after it was first answered
/// (motion::Upgrade), finished as `ask` finishes one.
pub(crate) fn finish_later(v: &ModelView, claim: &Claim, bound: Option<u32>, c: Check, opts: ClaimOptions) -> Check {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    let tol = claim_accuracy(v, claim, opts.accuracy);
    let inexact = v.doc.inexact(&roots, &v.certified());
    finished(v, claim, &roots, c, tol, inexact)
}

/// A claim's answer as it is decided (`decide`, or a still pair of `ink
/// interfere`'s proven beside the others: interfere::prove_pairs),
/// finished: settled at the accuracy asked, `tol`; an estimate where it
/// reads inexact inputs; in the claim's `quantity`, with the ranged params
/// it reads, at a moment or over a motion alike (a witness gives each
/// one's value).
pub(crate) fn answer(c: Check, tol: f64, inexact: crate::Inexact, quantity: Option<crate::units::Quantity>, ranged: Vec<(u32, (f64, f64))>) -> Check {
    let c = if inexact.any() { c.estimate() } else { c };
    Check { inexact, quantity, ranged, ..c.settle(tol) }
}

/// The accuracy a claim is checked to: the one given, else its default:
/// for a claim about shapes its parts' (a length: `finest`), for a comparison
/// (a value, a volume, an integral, a mean, a metric, a thickness) 0, and
/// a point or a range is decided exactly (Claim::accuracy).
pub(crate) fn claim_accuracy(v: &ModelView, claim: &Claim, given: Option<f64>) -> f64 {
    match claim.accuracy() {
        crate::Accuracy::Length => given.unwrap_or_else(|| parts_accuracy(v, claim).0),
        crate::Accuracy::Value if v.doc.shape_accuracy(claim) => given.unwrap_or_else(|| parts_accuracy(v, claim).0),
        crate::Accuracy::Value => given.unwrap_or(0.0),
        crate::Accuracy::Exact => 0.0,
    }
}

/// A claim about shapes' accuracy where none is given: `finest` of the
/// parts its shapes are made of (Doc::parts_of).
fn parts_accuracy(v: &ModelView, claim: &Claim) -> (f64, Vec<(String, f64)>) {
    if v.doc.bodies.iter().all(|b| v.doc.tolerance(&b.attrs).is_none()) {
        return (v.accuracy(), Vec::new());
    }
    finest(v.accuracy(), v.doc.claim_roots(claim).into_iter().flat_map(|f| v.doc.parts_of(f)))
}

/// The accuracy a claim about parts is checked to where none is given:
/// the finest of its parts' (Doc::parts_of), each its declared tolerance,
/// `(attr tolerance x)`, or else `model`'s (ModelView::accuracy), as is
/// material of no part. A part made to 0.2 mm against one made to 0.01 is
/// checked to 0.01; two made to 0.2, to 0.2. With the tolerances declared,
/// where one of them is that accuracy (none where the model's is finer).
pub fn finest(model: f64, parts: impl IntoIterator<Item = Option<(String, f64)>>) -> (f64, Vec<(String, f64)>) {
    let (mut declared, mut undeclared) = (Vec::new(), false);
    for p in parts {
        match p {
            Some(p) if !declared.contains(&p) => declared.push(p),
            Some(_) => {}
            None => undeclared = true,
        }
    }
    let tol = declared.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    match tol <= model || !undeclared {
        true if !declared.is_empty() => (tol, declared),
        _ => (model, Vec::new()),
    }
}

/// Where a claim's accuracy came from, where it is not its default.
#[derive(Clone, Debug, PartialEq)]
pub enum AccuracyBy {
    /// The claim's own `(accuracy x)`.
    Claim,
    /// The run's (`ClaimOptions::accuracy`, the CLI's `--accuracy`).
    Run,
    /// Its parts' declared tolerances (`finest`), each part's name and its.
    Tolerance(Vec<(String, f64)>),
}

/// What gives claim `c` its accuracy: its own (`own`), the run's (`run`,
/// for the claims about shapes) or its parts' tolerances (`finest`); None,
/// its default.
pub fn accuracy_by(v: &ModelView, c: &Claim, own: Option<f64>, run: Option<f64>) -> Option<AccuracyBy> {
    let length = c.accuracy_is_length() || v.doc.shape_accuracy(c);
    match (own, run) {
        (Some(_), _) => Some(AccuracyBy::Claim),
        _ if !length => None,
        (_, Some(_)) => Some(AccuracyBy::Run),
        _ => Some(parts_accuracy(v, c).1).filter(|parts| !parts.is_empty()).map(AccuracyBy::Tolerance),
    }
}

/// Whether the point `p` is inside the solid of `tape` at every value of the ranged params `spans`
/// (in `ctx`): by interval bounds of the field at the point, the widest range halved where they do not
/// show it, to STAYS readings. A point that is inside at one value and deep enough stays inside while
/// what moves about it moves less than its depth, which short ranges show.
pub(super) fn stays_inside(tape: &Tape, ctx: &Ctx, p: V3, spans: &[(u32, (f64, f64))]) -> bool {
    const STAYS: usize = 4096;
    let width = |b: &[(u32, (f64, f64))], i: usize| (b[i].1.1 - b[i].1.0) / (spans[i].1.1 - spans[i].1.0).max(f64::MIN_POSITIVE);
    // The independent ones are halved; what is derived from them follows (Ctx::derive).
    let spans: Vec<(u32, (f64, f64))> = spans.iter().copied().filter(|&(k, _)| !ctx.p_spans.is_derived(k)).collect();
    let (mut todo, mut read) = (vec![spans.clone()], 0);
    while let Some(b) = todo.pop() {
        read += 1;
        if read > STAYS {
            return false;
        }
        let mut c = ctx.clone();
        for &(k, s) in &b {
            c.params[k as usize] = ink_core::eval::span_mid(s);
            c.p_spans.narrow(k, s);
        }
        c.derive();
        if eval_interval(tape, &c, p.map(I::point), None).hi < 0.0 {
            continue;
        }
        let Some(i) = (0..b.len()).max_by(|&i, &j| width(&b, i).total_cmp(&width(&b, j))) else { return false };
        let (lo, hi) = b[i].1;
        let mid = ink_core::eval::span_mid((lo, hi));
        if !(lo < mid && mid < hi) {
            return false;
        }
        for half in [(lo, mid), (mid, hi)] {
            let mut next = b.clone();
            next[i].1 = half;
            todo.push(next);
        }
    }
    true
}
