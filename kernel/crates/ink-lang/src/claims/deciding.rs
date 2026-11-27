use super::*;

/// `check` at the view's moment, before the inexact inputs it reads are
/// taken into account.
pub(crate) fn decide(v: &mut ModelView, claim: Claim, opts: ClaimOptions) -> Result<Check, InkError> {
    let res = opts.res;
    let roots = v.doc.claim_roots(&claim);
    // The claim's accuracy: a shape claim's a length, a comparison's in
    // its own units (`claim_accuracy`).
    let tol = claim_accuracy(v, &claim, opts.accuracy);
    if let Some(c) = estimated(v, &claim, &roots, tol, opts) {
        return Ok(Check { ranged: ranged_reads(v, &roots), ..c });
    }
    let mut c = match claim {
        // Values in a region are compared exactly: their accuracy is 0.
        Claim::Empty(f) => shallow(v, f, tol, res),
        Claim::EmptyBounds { certain, possible } if let Some(region) = v.doc.print_bounds.values().find(|r| r.certain == certain && r.possible == possible).copied() => {
            check_overhang(v, region, res, tol)
        }
        Claim::EmptyBounds { certain, possible } => {
            let inner = shallow(v, certain, tol, res);
            if inner.status == Status::Fails {
                inner
            } else {
                // The candidate layer is constructed by the checker. Its thickness
                // is not the depth of the physical overhang: only exact emptiness
                // can discharge it. The requested accuracy still bounds refinement.
                let mut options = v.proof_options(res, tol);
                options.tol = 0.0;
                let (proof, whole) = v.prove_with(possible, &options);
                let exact = proof.exact;
                let outer = emptied(proof, whole);
                if exact && outer.status == Status::Holds {
                    Check { held: Some(Held::Outright), ..outer }
                } else {
                    // A partial depth bound on this artificial layer is
                    // not an attained accuracy for the physical ceiling, and
                    // a point of it no failure: what the proof reached is
                    // why it is not settled.
                    Check {
                        note: Some(
                            "a ceiling may lack support: no row along x or y was shown to bridge it, but a bridge may run another way. (unbridged-along part z0 a h l w dx dy) checks along the direction its bridges run",
                        ),
                        ..checked(Status::Unknown, outer.evidence)
                    }
                }
            }
        }
        Claim::Range(f) => Check { held: None, ..shallow(v, f, tol, res) },
        // A clearance d: no point of one within d of the other, exactly,
        // with cells down to half the accuracy, or a sixteenth of d where
        // finer (a clearance sets its own scale): a gap within about that
        // of d is unknown, never proven, so the accuracy can only leave an
        // answer undecided. A clearance of 0 is the separation question,
        // disjointness, to the accuracy.
        Claim::Clearance(a, b, d) if d > 0.0 => clearance(v, a, b, d, res, (if tol > 0.0 { 0.5 * tol } else { v.finest() }).min(d / 16.0)),
        Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) => {
            let f = v.doc.g.max(a, b);
            shallow(v, f, tol, res)
        }
        Claim::Subset(a, b) => {
            let f = subset_field(&mut v.doc.g, a, b);
            shallow_except(v, f, Some((a, b)), tol, res)
        }
        Claim::Touches(a, b) => touches(v, a, b, tol, res),
        Claim::Point { f, at, inside } => point(v, f, at, inside),
        Claim::Integral { parts, frames, cmp, bound } => {
            let c = Comparison { cmp, bound, tol };
            let parts: Vec<Part> = parts.iter().enumerate().map(|(k, (s, gs))| (*s, gs.iter().map(|&g| Some(g)).collect(), frames.get(k).copied().flatten())).collect();
            let each = |i: &[f64]| -> Vec<Status> { i.iter().map(|&x| point_status(c, x)).collect() };
            let (s, sampled) = estimate_scope(v, &parts, res, |i| all_of(each(i).into_iter()));
            let integrals: Vec<(f64, Status)> = s.sums.iter().copied().zip(each(&s.sums)).collect();
            // Decided by its proven bounds: proven, but where ranged params were read at sampled values.
            let sampled_scope = sampled.is_some();
            let c = checked(s.status, Evidence::Integral { integrals, res: s.res, sampled });
            Check { unbounded: !s.whole, ..if sampled_scope { c.estimate() } else { c } }
        }
        Claim::Mean { parts, frames, cmp, bound } => mean(v, &parts, &frames, Comparison { cmp, bound, tol }, res),
        Claim::Metric { metric, cmp, value: bound, excluding } => {
            // A package solves at one point of the scope: where its form
            // reads a ranged axis (itself, or through a form it reads, as a
            // thermal study reads its heat study), its answer says nothing
            // about the range.
            if scope_reads(v, &roots) {
                let names = &v.doc.params;
                let num = |x: f64| (x * 1e9).round() / 1e9;
                let read: Vec<String> = ranged_reads(v, &roots).iter().map(|&(k, (a, b))| format!("{} ({}..{})", names[k as usize].name, num(a), num(b))).collect();
                let e = InkError {
                    hint: Some(at_one_value(v, &roots)),
                    ..InkError::msg(format!("a package's metric is solved at one value of each param, but this form reads ranged {}. No solve covers the range", read.join(", ")))
                };
                checked(Status::Unknown, Evidence::Unsolved(e))
            } else if metric.shape.is_some_and(|f| !v.whole(f)) {
                // Of a shape seen only in the window, the window's faces
                // would be its walls.
                Check { unbounded: true, ..checked(Status::Unknown, Evidence::Unexamined) }
            } else {
                // A count of voids is to a length on its shape, exact on
                // its value.
                let (tol, accuracy) = if metric.shape.is_some() && v.doc.registry.shape_accuracy(metric.metric) { (0.0, Some(tol)) } else { (tol, None) };
                metric_claim(v, &metric, Comparison { cmp, bound, tol }, if metric.shape.is_some() { res } else { None }, accuracy, &excluding)?
            }
        }
        // A constant against a constant: decided in real numbers, by their exact fractions where both have one
        // (12 in to the foot), else over the bands their roundings leave, and not settled within rounding of the
        // limit (sin 30 deg against 1/2). Its float is the value given. A constant whose band has no bound (a
        // fold of one through a pole) is read as its float, as any other claim reads it.
        Claim::Value(f, cmp, bound) if v.doc.g.real_value(f).is_some_and(|e| e.band.lo.is_finite() && e.band.hi.is_finite()) && ranged_reads(v, &[f]).is_empty() => {
            let (g, x) = (&v.doc.g, v.doc.g.constant_value(f).expect("a constant"));
            let (a, b) = (
                g.real_value(f).expect("a constant"),
                v.doc
                    .value_bounds
                    .get(&(f, bound.to_bits()))
                    .and_then(|&n| g.real_value(n))
                    .filter(|e| e.band.lo.is_finite() && e.band.hi.is_finite())
                    .unwrap_or(ink_core::exact::Exact { ratio: ink_core::exact::Ratio::of_f64(bound), band: I::point(bound) }),
            );
            let c = Comparison { cmp, bound: 0.0, tol };
            let status = match a.ratio.zip(b.ratio).and_then(|(p, q)| p.minus(q)) {
                Some(d) => point_status(c, d.n.signum() as f64),
                None => {
                    let d = ink_core::interval::bin(ink_core::graph::BinOp::Sub, a.band, b.band).0;
                    if c.holds_on(d) {
                        Status::Holds
                    } else if c.fails_on(d) {
                        Status::Fails
                    } else {
                        Status::Unknown
                    }
                }
            };
            checked(status, Evidence::Value(x))
        }
        Claim::Value(f, cmp, bound) if ranged_reads(v, &[f]).is_empty() => {
            let x = v.doc.scalar(f, &v.ctx);
            let c = Comparison { cmp, bound, tol };
            // A fixed scalar has a point bound. A real construction (a
            // derivative or a spatial expression evaluated at a station)
            // retains its enclosure even when no variable remains.
            // Each param it reads one value, but the arithmetic on it in real numbers: rounded outward, as
            // a range's is, not taken in floats as a def's value is (1e16 + 1 less 1e16 is 1, not 0).
            // A param written in digits is the fraction they write, as a literal is: where every number the
            // claim reads is exact and its arithmetic keeps fractions (334 mm against 20 links of 16.7),
            // the claim is decided on them.
            let exact = |n: ink_core::graph::NodeId| v.doc.exact_value(n, &v.ctx.params);
            let limit = v.doc.value_bounds.get(&(f, bound.to_bits())).map_or_else(|| ink_core::exact::Ratio::of_f64(bound), |&n| exact(n));
            if let Some(d) = exact(f).zip(limit).and_then(|(a, b)| a.minus(b)) {
                checked(point_status(Comparison { cmp, bound: 0.0, tol }, d.n as f64 / d.d as f64), Evidence::Value(x))
            } else {
                let ctx = v.doc.params_read(&[f]).into_iter().fold(v.ctx.clone(), |c, k| c.with_band(k as u32, (v.ctx.params[k], v.ctx.params[k])));
                let iv = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, f), &ctx, [I::point(0.0); 3], None);
                let value = if iv.contains(x) { x } else { iv.mid() };
                if c.holds_on(iv) {
                    checked(Status::Holds, Evidence::Value(value))
                } else if c.fails_on(iv) {
                    checked(Status::Fails, Evidence::Value(value))
                } else if v.doc.g.value_depends_on(f, |n| n == Var::T) && point_status(c, x) == Status::Holds {
                    Check { accuracy: None, ..checked(Status::Holds, Evidence::Value(x)) }
                } else {
                    checked(Status::Unknown, Evidence::Value(value))
                }
            }
        }
        Claim::Value(f, cmp, bound) => {
            // For every value of each ranged param: interval bounds over
            // pieces, halving undecided ones; a point violating it deeper
            // than the accuracy disproves.
            let tape = Tape::new(&v.doc.g, f);
            let result = ranged::value_span(&tape, &v.ctx, [0.0; 3], Comparison { cmp, bound, tol }, &v.proof_options(res, 0.0));
            checked(Status::of(result.status), Evidence::ValueSpan(result))
        }
        Claim::Overlaps(a, b) => {
            // Overlap is proven by a witness deeper than the accuracy in
            // both: exactly when the shallowness of the intersection is
            // disproven. Disproven only by an intersection proven empty
            // outright: one shown no deeper than the accuracy may still
            // overlap.
            let f = v.doc.g.max(a, b);
            let (p, whole) = v.prove_with(f, &v.proof_options(res, tol));
            let apart = p.exact();
            let mut c = in_window(Status::of(p.status()), whole, Evidence::Overlaps(p, None));
            c.status = match c.status {
                Status::Fails => Status::Holds,
                Status::Holds if apart => Status::Fails,
                _ => Status::Unknown,
            };
            // Depth is distance from outside: every point of a box thinner
            // than twice the accuracy is within the accuracy of leaving it.
            // The intersection lies in its own proven box and in each
            // shape's.
            if c.status == Status::Unknown
                && whole
                && tol > 0.0
                && let (Some(i), Some(ea), Some(eb)) = (v.bounds(f), v.bounds(a), v.bounds(b))
                && let side = (0..3).map(|k| i.1[k].min(ea.1[k]).min(eb.1[k]) - i.0[k].max(ea.0[k]).max(eb.0[k])).fold(f64::INFINITY, f64::min)
                && side < 2.0 * tol
                && let Evidence::Overlaps(_, thin) = &mut c.evidence
            {
                *thin = Some(side.max(0.0));
            }
            // Over ranged params a witness is at one value of them, and
            // says nothing of the rest unless the point itself is in both
            // at every value (interval bounds of the field there over the
            // ranges: a pin's own point, in the frame the pin rests in,
            // while its hole turns about it). Proven apart at every value
            // fails.
            if c.status == Status::Holds && !ranged_reads(v, &[a, b]).is_empty() {
                let stays = match &c.evidence {
                    Evidence::Overlaps(Proof { verdict: query::Verdict::Disproven { witness, .. }, .. }, _) => {
                        stays_inside(&Tape::new(&v.doc.g, f), &v.ctx, *witness, &ranged_reads(v, &[a, b]))
                    }
                    _ => false,
                };
                if !stays {
                    (c.status, c.note) = (Status::Unknown, Some(AT_ONE_VALUE));
                }
            }
            c
        }
    };
    c.ranged = ranged_reads(v, &roots);
    c.accuracy = c.accuracy.map(|_| tol);
    Ok(c)
}

/// A claim reading params set as estimates that gauge their own error
/// (ParamProvider::error). Every band it reads is carried, however small,
/// or the claim is unknown: none is dropped for its size. A value claim at
/// a point holds where it holds for every value in their bands (interval
/// bounds; its value then the bands' worst), fails where it fails for every
/// one, and is unknown between, saying so. A claim about shapes is asked
/// over the bands as over ranged params: it holds where it holds for every
/// value, fails where a point that fails it does so for every value, and is
/// unknown between. A metric of a form whose response is linear in what it
/// reads of the model (a static study's loads) has the bands carried into
/// its own answer (ModelView::measure_over), and is not judged here. A band
/// read through another solve, also or only, has no bound; through an
/// estimate's guard (estimate-when), a formula, it is carried where the
/// guard's condition holds over all of it. None: it reads none.
fn estimated(v: &mut ModelView, claim: &Claim, roots: &[u32], tol: f64, opts: ClaimOptions) -> Option<Check> {
    let mut band: Vec<(usize, f64)> = Vec::new();
    let solves = v.doc.solves_read(roots);
    for &i in &solves {
        let s = &v.doc.solves[i];
        let Some(p) = v.doc.ext[s.ext].form.provider(&v.doc) else { continue };
        if let Some(w) = p.error(&v.ctx.params, v.ctx.t) {
            band.extend(s.owns.iter().copied().zip(w));
        }
    }
    // Only what the claim reads moves it: a param its own expression names, by its band; one it
    // reads through another solve, also or only, by no bound, since a band is not carried through a
    // solve. A solve's other outputs, read by nothing here, do not.
    let read = v.doc.params_read(roots);
    // An estimate's guard (estimate-when) is a formula, not a solve: where its condition keeps its sign
    // over every band, it is 0 at every value there, and what it reads is carried through the formula
    // as any arithmetic is. Where the condition may change sign within a band, the estimate does not
    // apply over all of it: that band has no bound, and the guard's reason says why.
    let finite: Vec<(u32, (f64, f64))> = band
        .iter()
        .filter(|b| b.1.is_finite() && v.ctx.p_spans.get(b.0 as u32).is_none())
        .map(|&(k, w)| {
            let b = add_o(I::point(v.ctx.params[k]), I::new(-w, w));
            (k as u32, (b.lo, b.hi))
        })
        .collect();
    let over_bands = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&finite), ..v.ctx.clone() };
    let unbounded = |k: usize| band.iter().any(|b| b.0 == k && !b.1.is_finite());
    let mut crossed: Vec<String> = Vec::new();
    let through: std::collections::HashSet<usize> = solves
        .iter()
        .filter(|&&i| match v.doc.guard_over(i, &over_bands) {
            Some(Ok(())) => v.doc.solves[i].reads.iter().any(|&k| unbounded(k)),
            Some(Err(why)) => {
                crossed.push(why);
                true
            }
            None => true,
        })
        .flat_map(|&i| v.doc.solves[i].reads.iter().copied())
        .collect();
    // What a form linear in the model's values reads of them: its own to carry.
    let carried: Vec<usize> = match claim {
        Claim::Metric { metric: m, .. } if m.shape.is_none() => match v.doc.ext_form(Some(m.package), &m.form) {
            Some(e) if e.form.linear_in_worst(m.metric).is_some() => v.doc.params_read(&e.form.worst_reads().iter().flat_map(|(_, n)| n.iter().copied()).collect::<Vec<u32>>()),
            _ => Vec::new(),
        },
        _ => Vec::new(),
    };
    let mut band: Vec<(usize, f64)> = band
        .into_iter()
        .filter(|b| b.1 != 0.0 && !carried.contains(&b.0))
        .filter_map(|b| match (read.contains(&b.0), through.contains(&b.0)) {
            (true, false) => Some(b),
            // Read through another solve too, directly or not: that solve's answer moves with it by no bound.
            (_, true) => Some((b.0, f64::INFINITY)),
            (false, false) => None,
        })
        .collect();
    // A shape built within an allowance of the one written (Doc::allowances): a bound, not a
    // measure, so a claim it reads holds or fails proven where it does so over all of it. Read
    // directly only: what reads it through a solve, or an integral, a mean or a metric (each an
    // estimate with no bound of its own), is not moved by it here.
    let bounded = band.is_empty();
    let shape_or_value = !matches!(claim, Claim::Metric { .. } | Claim::Integral { .. } | Claim::Mean { .. });
    let allowed: Vec<(usize, f64)> = v.doc.allowances.iter().copied().filter(|a| a.1 != 0.0 && read.contains(&a.0)).collect();
    if shape_or_value {
        band.extend(allowed.iter().copied());
    }
    if band.is_empty() {
        return None;
    }
    // A proven answer over the bands stays one where they are all allowances.
    let accuracy = |c: Option<f64>| if bounded { c.map(|_| tol) } else { None };
    // A claim about shapes, over the bands as over ranged params (asked again here with them
    // ranged, it is any claim over ranged params).
    if !matches!(claim, Claim::Value(..) | Claim::Metric { .. } | Claim::Integral { .. } | Claim::Mean { .. }) && band.iter().all(|b| b.1.is_finite()) {
        if band.iter().all(|b| v.ctx.p_spans.get(b.0 as u32).is_some()) {
            return None;
        }
        if band.iter().all(|b| v.ctx.p_spans.get(b.0 as u32).is_none()) {
            let spans: Vec<(u32, (f64, f64))> = band
                .iter()
                .map(|&(k, w)| {
                    let b = add_o(I::point(v.ctx.params[k]), I::new(-w, w));
                    (k as u32, (b.lo, b.hi))
                })
                .collect();
            let over = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&spans), ..v.ctx.clone() };
            // What fails it at every value: a point deep in what it says is empty at both ends of the
            // bands (found there), that is in it at every value between (interval bounds over the bands).
            let inside = match *claim {
                Claim::Empty(a) => Some(a),
                Claim::Disjoint(a, b) => Some(v.doc.g.max(a, b)),
                Claim::Subset(a, b) => {
                    let out = v.doc.g.neg(b);
                    Some(v.doc.g.max(a, out))
                }
                _ => None,
            };
            // The claim with each band's param at a value (`ends`: one per band).
            let at_ends = |v: &mut ModelView, n: u32, ends: &[f64]| {
                let at: Vec<(Var, u32)> = band.iter().zip(ends).map(|(&(k, _), &x)| (Var::Param(k as u32), v.doc.g.constant(x))).collect();
                v.doc.g.subst(n, &at)
            };
            // Where what it says is empty moves one way with each band's param at every point
            // (Tape::param_way: through sums, min, max, negations and the like only), it is emptiest
            // with each at the end that raises it and fullest at the other: it holds at every value
            // where it holds at its fullest, and a point in it at its emptiest is in it at every
            // value. Else every value is asked at once, and its ends are the bands' all low and all
            // high.
            let ways: Option<Vec<i8>> = inside.and_then(|f| {
                let tape = Tape::new(&v.doc.g, f);
                band.iter().map(|&(k, _)| tape.param_way(k as u32).map_or(Some(1), |w| (w != 0).then_some(w))).collect()
            });
            let end = |k: usize, w: f64, up: bool| v.ctx.params[k] + if up { w } else { -w };
            let corners: Vec<Vec<f64>> = match &ways {
                Some(ways) => vec![band.iter().zip(ways).map(|(&(k, w), &d)| end(k, w, d > 0)).collect()],
                None => [false, true].iter().map(|&up| band.iter().map(|&(k, w)| end(k, w, up)).collect()).collect(),
            };
            let asked = match (&ways, claim) {
                (Some(ways), c) => {
                    let fullest: Vec<f64> = band.iter().zip(ways).map(|(&(k, w), &d)| end(k, w, d < 0)).collect();
                    let c = match *c {
                        Claim::Empty(a) => Claim::Empty(at_ends(v, a, &fullest)),
                        Claim::Disjoint(a, b) => Claim::Disjoint(at_ends(v, a, &fullest), at_ends(v, b, &fullest)),
                        Claim::Subset(a, b) => Claim::Subset(at_ends(v, a, &fullest), at_ends(v, b, &fullest)),
                        _ => unreachable!("a claim with what it says is empty"),
                    };
                    decide(v, c, opts)
                }
                (None, _) => v.at_point(over.clone(), |v| decide(v, claim.clone(), opts)),
            };
            if let Ok(c) = &asked
                && c.status == Status::Holds
            {
                return Some(Check { accuracy: accuracy(c.accuracy), ranged: Vec::new(), ..c.clone() });
            }
            // A point fails at every value where the field's bounds over the bands at it say so.
            if let Claim::Point { f, at, inside: is_in } = *claim {
                let p = at.map(|n| v.doc.scalar(n, &v.ctx));
                let iv = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, f), &over, p.map(I::point), None);
                if (is_in && iv.lo >= 0.0) || (!is_in && iv.hi < 0.0) {
                    return Some(Check { accuracy: accuracy(Some(tol)), ..checked(Status::Fails, Evidence::Point { at: p, bounds: (iv.lo, iv.hi) }) });
                }
            }
            if let Some(f) = inside {
                // The witness read at its own moment and ranged params' values, over the bands.
                let at = |t: f64, params: &[(u32, f64)]| {
                    let mut at = Ctx { t, t_span: None, ..over.clone() };
                    for &(k, x) in params {
                        at.params[k as usize] = x;
                        at.p_spans.narrow(k, (x, x));
                    }
                    at
                };
                let fails = |v: &ModelView, found: &Verdict| match found {
                    Verdict::Disproven { witness, t, params, .. } => eval_interval(&Tape::new(&v.doc.g, f), &at(*t, params), witness.map(I::point), None).hi < 0.0,
                    _ => false,
                };
                // The point the claim asked over the bands found, if it is in it at every value.
                if let Ok(c) = &asked
                    && let Evidence::Empty(found) = &c.evidence
                    && fails(v, &found.verdict)
                {
                    return Some(Check { accuracy: accuracy(Some(tol)), ..checked(Status::Fails, Evidence::Empty(found.clone())) });
                }
                // Else one searched for at its ends (its emptiest, or both), as the claim's own search
                // is asked (a subset's inside a and not inside b), with what is left of one search's
                // boxes: none where the search at its fullest used them all (or the time limit
                // passed). Its emptiest is inside its fullest, deeper: a point there is one that
                // search was given the boxes to find.
                let spent = ways.is_some()
                    && matches!(&asked, Ok(c) if c.status == Status::Unknown
                        && matches!(&c.evidence, Evidence::Empty(p) if matches!(p.verdict, Verdict::Unknown { ended: query::Ended::Limit | query::Ended::Stopped, .. })));
                if !spent {
                    let found = match *claim {
                        Claim::Subset(a, b) => {
                            let ends: Vec<(u32, u32)> = corners.iter().map(|e| (at_ends(v, a, e), at_ends(v, b, e))).collect();
                            let a = v.doc.reduce(ink_core::graph::BinOp::Max, &ends.iter().map(|e| e.0).collect::<Vec<_>>());
                            let b = v.doc.reduce(ink_core::graph::BinOp::Min, &ends.iter().map(|e| e.1).collect::<Vec<_>>());
                            let both = subset_field(&mut v.doc.g, a, b);
                            v.prove_except(both, Some((a, b)), opts.res, tol).0
                        }
                        _ => {
                            let ends: Vec<u32> = corners.iter().map(|e| at_ends(v, f, e)).collect();
                            let both = v.doc.reduce(ink_core::graph::BinOp::Max, &ends);
                            v.prove(both, opts.res, tol).0
                        }
                    };
                    if fails(v, &found.verdict) {
                        return Some(Check { accuracy: accuracy(Some(tol)), ..checked(Status::Fails, Evidence::Empty(found)) });
                    }
                }
            }
            // A clearance fails at every value where it fails with each shape at its smallest: the
            // distance between two shapes only shrinks as either grows. Each band's param at the end
            // that raises both fields (one reading it alone, its own); none where one moves either
            // way, or the two each their own way.
            if let Claim::Clearance(a, b, d) = *claim {
                let (ta, tb) = (Tape::new(&v.doc.g, a), Tape::new(&v.doc.g, b));
                let smallest: Option<Vec<f64>> = band
                    .iter()
                    .map(|&(k, w)| {
                        let up = match (ta.param_way(k as u32), tb.param_way(k as u32)) {
                            (Some(0), _) | (_, Some(0)) => return None,
                            (Some(x), Some(y)) if x != y => return None,
                            (Some(x), _) | (_, Some(x)) => x > 0,
                            (None, None) => return Some(v.ctx.params[k]),
                        };
                        Some(v.ctx.params[k] + if up { w } else { -w })
                    })
                    .collect();
                if let Some(e) = smallest {
                    let c = Claim::Clearance(at_ends(v, a, &e), at_ends(v, b, &e), d);
                    if let Ok(c) = decide(v, c, opts)
                        && c.status == Status::Fails
                    {
                        return Some(Check { accuracy: accuracy(c.accuracy), ..c });
                    }
                }
            }
            // Not settled at its fullest: that is why, not the bands (its search's own answer, as
            // where it ended at its limit of boxes).
            if ways.is_some()
                && let Ok(c) = asked
                && c.status == Status::Unknown
            {
                return Some(c);
            }
        }
    }
    let mut at = None;
    if let Claim::Value(f, cmp, bound) = *claim
        && ranged_reads(v, &[f]).is_empty()
        && band.iter().all(|b| b.1.is_finite())
    {
        let c = Comparison { cmp, bound, tol };
        let spans: Vec<(u32, (f64, f64))> = band
            .iter()
            .map(|&(k, w)| {
                let b = add_o(I::point(v.ctx.params[k]), I::new(-w, w));
                (k as u32, (b.lo, b.hi))
            })
            .collect();
        let over = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&spans), ..v.ctx.clone() };
        let tape = Tape::new(&v.doc.g, f);
        let iv = ink_core::eval::eval_interval(&tape, &over, [I::point(0.0); 3], None);
        let x = ink_core::eval::eval_point(&tape, &v.ctx, [0.0; 3]);
        if c.holds_on(iv) {
            let worst = if matches!(cmp, Cmp::Lt | Cmp::Le) { iv.hi } else { iv.lo };
            return Some(Check { accuracy: accuracy(Some(tol)), ..checked(Status::Holds, Evidence::Value(worst)) });
        }
        if c.fails_on(iv) || x.is_nan() {
            return Some(Check { accuracy: accuracy(Some(tol)), ..checked(Status::Fails, Evidence::Value(x)) });
        }
        at = Some(x);
    }
    let named = |b: &[(usize, f64)]| b.iter().map(|&(k, w)| format!("{} {:.6} within {w:.6}", v.doc.params[k].name, v.ctx.params[k])).collect::<Vec<_>>().join(", ");
    let (allowances, estimates): (Vec<(usize, f64)>, Vec<(usize, f64)>) = band.iter().partition(|b| allowed.iter().any(|a| a.0 == b.0));
    let what = match (allowances.is_empty(), estimates.is_empty()) {
        (true, _) => format!("it reads an estimate with a band either way ({}), from its own measure, not a bound (its note says how)", named(&estimates)),
        (false, true) => format!("it reads a shape built within an allowance of the one written ({}), a bound on how far the shape built may be from it", named(&allowances)),
        (false, false) => format!(
            "it reads an estimate with a band either way ({}), from its own measure, not a bound (its note says how), and a shape built within an allowance of the one written ({})",
            named(&estimates),
            named(&allowances)
        ),
    };
    let why = format!(
        "{what}. {}",
        match at {
            Some(x) => format!("Its value {x:.6} neither holds nor fails over all of that band"),
            None => "It neither holds at every value in that band nor fails at every one, or it reads the band through a solve, which carries none".into(),
        }
    );
    let why = match crossed.is_empty() {
        true => why,
        false => format!("{why}. An estimate it reads applies only under a condition the band may cross: {}", crossed.join("; ")),
    };
    Some(checked(Status::Unknown, Evidence::Unsolved(InkError::msg(why))))
}
