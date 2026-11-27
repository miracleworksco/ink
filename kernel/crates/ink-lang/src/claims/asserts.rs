use super::*;

/// Check the document's asserts (only the one named `only`, if given; only
/// those written in the file loaded, not in its imports, if `own`), at the
/// view's params and time. Each is reported as it is decided (`progress`).
pub fn run_asserts(v: &mut ModelView, only: Option<&str>, own: bool, opts: ClaimOptions) -> Result<AssertReport, InkError> {
    run_asserts_with(v, only, own, opts, None, |_| None)
}

/// An answer's note on what it asked at one moment though the question
/// reads t (claims, pairs: each named, with the moment it recorded): said
/// once, for the answer. None: nothing was.
pub fn one_moment_note(asked: &[(String, f64)]) -> Option<String> {
    let mut at: Vec<(f64, Vec<&str>)> = Vec::new();
    for (what, t) in asked {
        match at.iter_mut().find(|m| m.0 == *t) {
            Some(m) => m.1.push(what),
            None => at.push((*t, vec![what])),
        }
    }
    let each: Vec<String> = at.iter().map(|(t, ws)| format!("{} at t = {t}", ws.join(", "))).collect();
    (!each.is_empty()).then(|| format!("answered at one moment only, though each reads t and no (t-span) covers it: {}. Add (t-span a b) or use --t-span", each.join("; ")))
}

/// `run_asserts`, each answer reported with its row of the caller's answer,
/// `row` (so a run stopped before its answer gives the rows done as the
/// whole answer would). At a param's `case` (its values, a run of a
/// declared case), each line it reports says so, and it reports no total:
/// its caller's counts every case.
pub fn run_asserts_with(
    v: &mut ModelView,
    only: Option<&str>,
    own: bool,
    opts: ClaimOptions,
    case: Option<&serde_json::Value>,
    row: impl Fn(&AssertResult) -> Option<serde_json::Value>,
) -> Result<AssertReport, InkError> {
    v.scope_question()?;
    let mut not_checked = v.doc.unchecked_claims()?;
    if let Some(n) = only {
        let source = v
            .doc
            .asserts
            .iter()
            .find(|x| x.name == n)
            .map(|x| (&x.file, x.line, v.doc.is_own(x)))
            .or_else(|| not_checked.iter().find(|x| x.name == n).map(|x| (&x.file, x.line, true)));
        match source {
            None => {
                let names = v.doc.asserts.iter().map(|x| x.name.as_str()).chain(not_checked.iter().map(|x| x.name.as_str()));
                let hint = crate::choices(n, "Its requirements", names);
                return Err(InkError::msg(format!("no requirement named '{n}'. {hint}")));
            }
            Some((file, line, false)) if own => {
                return Err(InkError::msg(format!("requirement '{n}' is imported (from {file}:{line}): --imports asks its imports' too")));
            }
            Some(_) => {}
        }
    }
    // An unchecked placement belongs to the run even if its rating came from an import.
    not_checked.retain(|x| only.is_none_or(|n| n == x.name));
    let imported_not_run = own.then(|| v.doc.asserts.iter().filter(|x| !v.doc.is_own(x)).count());
    type Listed = (String, String, String, Claim, Option<u32>, (Option<f64>, Option<f64>, Option<crate::units::Quantity>), Vec<(String, u32)>, Vec<usize>, crate::Scope);
    let list: Vec<Listed> = v
        .doc
        .asserts
        .iter()
        .filter(|x| only.is_none_or(|n| n == x.name) && (!own || v.doc.is_own(x)))
        .map(|x| {
            let source = format!("{}:{}", x.file, x.line);
            (x.name.clone(), x.text.clone(), source, x.claim.clone(), x.bound, (x.res, x.accuracy, x.quantity), x.terms.clone(), x.fixed.clone(), x.scope)
        })
        .collect();
    let mut asserts: Vec<AssertResult> = Vec::new();
    // A template's generated limits and relations remain addressable by
    // both names. Within this one run (one model and set of cases/params),
    // reuse their proof only for exactly the same question and term fields.
    let mut proofs = std::collections::HashMap::new();
    if case.is_none() {
        progress(serde_json::json!({"total": list.len() + not_checked.len(), "of": "requirements", "kept": true}));
        // What each will cost, for a caller that may stop the run (ink
        // mcp): its answer then says which were the costly ones. Said once,
        // by the run at the values: a case's run is its caller's.
        if ink_kernel::progress::reporting() {
            let ks: Vec<usize> = (0..v.doc.asserts.len()).filter(|&k| list.iter().any(|l| l.0 == v.doc.asserts[k].name)).collect();
            let fields: Vec<u32> = ks.iter().flat_map(|&k| crate::cost::sketched(v, k)).collect();
            let sk = crate::cost::Sketches::surveyed(v, &fields);
            let mut each = serde_json::Map::new();
            for k in ks {
                let name = v.doc.asserts[k].name.clone();
                each.insert(name, crate::cost::of_assert(v, &sk, k, opts.res).map_or(serde_json::Value::Null, |c| c.json()));
            }
            progress(serde_json::json!({"estimates": each}));
        }
    }
    let report_done = |a: &AssertResult| {
        let mut done = serde_json::json!({"done": a.name, "status": a.check.status.as_str(), "ms": a.ms});
        if let Some(c) = case {
            done["case"] = c.clone();
        }
        if let Some(r) = row(a) {
            // Its row's verdict, where the row gives one, is what a stopped answer counts.
            if let Some(verdict) = r.get("verdict") {
                done["status"] = verdict.clone();
            }
            done["row"] = r;
        }
        progress(done);
    };
    for (name, claim, source, c, bound, (own_res, tol, quantity), terms, fixed, file) in list {
        // Over the motion of the file it is written in, imported or not
        // (ModelView::motion_of).
        v.motion_of(file).map_err(|e| InkError { message: format!("requirement '{name}' ({source}): {}", e.message), ..e })?;
        // The bound at the view's params (solved ones included, wherever
        // the solve is in the file).
        let c = v.doc.claim_at(&c, bound, &v.ctx);
        let roots: Vec<u32> = v.doc.claim_roots(&c).into_iter().chain(bound).collect();
        // What it sweeps, said as it starts: the ranged params it reads.
        let mut running = serde_json::json!({"running": name});
        if let Some(c) = case {
            running["case"] = c.clone();
        }
        let ranged = ranged_reads(v, &roots);
        if !ranged.is_empty() {
            running["ranged"] = ranged.iter().map(|&(k, s)| (v.doc.params[k as usize].name.clone(), serde_json::json!([s.0, s.1]))).collect::<serde_json::Map<_, _>>().into();
        }
        progress(running);
        let res = claim_res(own_res, opts.res);
        let accuracy_by = accuracy_by(v, &c, tol, opts.accuracy);
        let start = std::time::Instant::now();
        let count = || {
            use std::sync::atomic::Ordering::Relaxed;
            (ink_kernel::cache::MISSES.load(Relaxed), ink_kernel::cache::HITS.load(Relaxed))
        };
        let before = count();
        // Generated checks compare a margin with zero. Restrict sharing to
        // ordinary fields: a provider or solve can change its inputs when
        // asked. Compare node identities, never their values at one point.
        let key = match c {
            Claim::Value(f, Cmp::Ge, 0.0) if !terms.is_empty() && bound.is_none() && !own && only.is_none() && !ink_kernel::stop::reached() => {
                let fields: Vec<_> = std::iter::once(f).chain(terms.iter().map(|(_, n)| *n)).collect();
                v.doc.solves_read(&fields).is_empty().then(|| {
                    (
                        f,
                        terms.clone(),
                        fixed.clone(),
                        tol.map(f64::to_bits),
                        quantity,
                        own_res.map(f64::to_bits),
                        v.ctx.t.to_bits(),
                        v.motion.map(|(a, b)| (a.to_bits(), b.to_bits())),
                    )
                })
            }
            _ => None,
        };
        let previous = key.as_ref().and_then(|k| proofs.get(k)).map(|&i: &usize| &asserts[i]);
        let reused = previous.is_some();
        let operands = match &c {
            Claim::Disjoint(p, q) | Claim::Subset(p, q) | Claim::Touches(p, q) => (vec![*p, *q], false),
            Claim::Clearance(p, q, _) => (vec![*p, *q], true),
            Claim::Empty(f) | Claim::Point { f, .. } => (vec![*f], false),
            _ => (Vec::new(), false),
        };
        let (check, terms) = if let Some(a) = previous {
            (a.check.clone(), a.terms.clone())
        } else {
            let answer = asserted(v, c, bound, res.map(|r| r.0), (tol, quantity), &fixed, &terms, opts)
                .map_err(|e| InkError { message: format!("requirement '{name}' ({source}): {}", e.message), ..e })?;
            if answer.0.accuracy.is_some()
                && let Some(k) = key
            {
                proofs.insert(k, asserts.len());
            }
            answer
        };
        let ms = start.elapsed().as_secs_f64() * 1e3;
        let after = count();
        let cached = reused || (after.0 == before.0 && after.1 > before.1 && !matches!(check.evidence, Evidence::Stopped));
        let declared_at = v.doc.claim_declarations.get(&name).cloned();
        let at_witness = bodies_at(v, &operands.0, operands.1, &check);
        let a = AssertResult { name, claim, source, declared_at, res, check, ms, cached, terms, accuracy_by, at_witness };
        report_done(&a);
        let one_moment = a.check.one_moment.is_some();
        asserts.push(a);
        // The answer's note as it grows, so a run stopped before its
        // answer (ink mcp) gives it with the rows done.
        if one_moment {
            progress(serde_json::json!({"answer": {"note": assert_note(&asserts)}}));
        }
    }
    for limit in not_checked {
        let result = AssertResult {
            name: limit.name,
            claim: limit.claim.into(),
            source: format!("{}:{}", limit.file, limit.line),
            declared_at: limit.declared_at,
            check: checked(Status::Unknown, Evidence::NotChecked(limit.reason)).estimate(),
            ms: 0.0,
            res: None,
            terms: Vec::new(),
            cached: false,
            accuracy_by: None,
            at_witness: Vec::new(),
        };
        report_done(&result);
        asserts.push(result);
    }
    Ok(AssertReport { asserts, imported_not_run })
}

/// One claim asked as `ink require` asks it (`c` at the view's params,
/// Doc::claim_at; `bound` its bound's field): what it reads provided while
/// it is checked, and no longer; unknown where it rests on a solve with no
/// value here, or was read at load with a ranged param as one number
/// (`fixed`); a bound reading a ranged param moved into the claim; at cells
/// of `res` (claim_res) and the claim's accuracy and quantity, if it gives
/// them. With the values of `terms` where its value is. Every claim asked
/// of a design is asked here: an assert's, an optimization's requirement.
#[allow(clippy::too_many_arguments)]
pub fn asserted(
    v: &mut ModelView,
    mut c: Claim,
    bound: Option<u32>,
    res: Option<f64>,
    (tol, quantity): (Option<f64>, Option<crate::units::Quantity>),
    fixed: &[usize],
    terms: &[(String, u32)],
    opts: ClaimOptions,
) -> Result<(Check, Vec<(String, f64)>), InkError> {
    let roots = v.doc.claim_motion_inputs(&c, bound);
    // Resting on a solve with no value here (one depending on t, over a
    // span; one that does not solve at this moment), it is not decided:
    // unknown, with the solve's error; the rest are.
    // What is provided only where read (an effort) is provided while
    // the claim reading it is checked, and no longer.
    let wanted = v.doc.wanted.clone();
    let refused = v.want(&roots).err();
    // A moving claim can retry a point refusal at each step. set_moment
    // cannot re-solve over parameter ranges. Keep a refusal for that scope.
    let unsolved = refused.or_else(|| v.refusal(&roots)).filter(|_| !v.ctx.p_spans.is_empty() || !v.motion.is_some_and(|span| moves(v, &roots, span)));
    // Read at load with a ranged param as one number (AssertDef::fixed),
    // it is not decided over the range.
    let fixed: Vec<usize> = fixed.iter().copied().filter(|&k| v.ctx.p_spans.get(k as u32).is_some()).collect();
    let unsolved = unsolved.or_else(|| {
        (!fixed.is_empty()).then(|| InkError {
            hint: Some(v.doc.at_one_value(&fixed, &v.ctx.p_spans)),
            ..InkError::msg(format!(
                "it reads ranged {} as one number on loading (its relation's kind reads its parts' shapes). No proof covers the range",
                fixed.iter().map(|&k| v.doc.params[k].name.as_str()).collect::<Vec<_>>().join(", ")
            ))
        })
    });
    // A bound that reads a ranged param varies with it: a value claim
    // compares the difference with 0 (one expression, so a param in
    // both keeps its dependency, and a witness fails at its own
    // values); other bounds are numbers by definition. Either way the
    // bound is then in the claim, not beside it.
    let mut note = None;
    let mut beside = bound;
    if let Some(b) = bound.filter(|&b| !ranged_reads(v, &[b]).is_empty()) {
        match c {
            Claim::Value(f, cmp, _) => {
                c = Claim::Value(v.doc.g.sub(f, b), cmp, 0.0);
                beside = None;
                note = Some("its bound reads a ranged param, so the value checked is the expression minus the bound, against 0");
            }
            Claim::Mean { parts, frames, cmp, .. } => {
                let parts = parts.into_iter().map(|(s, w, hs)| (s, w, hs.into_iter().map(|h| v.doc.g.sub(h, b)).collect())).collect();
                c = Claim::Mean { parts, frames, cmp, bound: 0.0 };
                beside = None;
                note = Some("its bound reads a ranged param, so the mean checked is that of the field minus the bound, against 0");
            }
            _ => {
                let ks: Vec<usize> = ranged_reads(v, &[b]).iter().map(|(k, _)| *k as usize).collect();
                let names: Vec<&str> = ks.iter().map(|&k| v.doc.params[k].name.as_str()).collect();
                v.doc.wanted = wanted;
                return Err(InkError {
                    hint: Some(format!("a bound is one number. Write the varying part into the claim, e.g. (value (- e b) > 0), or {}", v.doc.at_one_value(&ks, &v.ctx.p_spans))),
                    ..InkError::msg(format!("its bound reads ranged param(s) {}", names.join(", ")))
                });
            }
        }
    }
    // A run's accuracy is a length, for the claims about shapes.
    let tol = tol.or(opts.accuracy.filter(|_| c.accuracy_is_length() || v.doc.shape_accuracy(&c)));
    let mut check = match unsolved {
        // Past the call's time limit, the rest are not begun.
        _ if ink_kernel::stop::reached() => {
            Check { quantity: v.doc.accuracy_quantity(&c), ranged: ranged_reads(v, &roots), ..checked(Status::Unknown, Evidence::Stopped).estimate() }
        }
        Some(e) => Check { quantity: v.doc.accuracy_quantity(&c), ..checked(Status::Unknown, Evidence::Unsolved(e)).estimate() },
        None => ask_bounded(v, c, beside, bound, ClaimOptions { res, accuracy: tol, ..opts })?.wait(v)?,
    };
    check.note = check.note.or(note);
    // Asked at one moment though it reads t, no moment asked: recorded
    // with the answer, at the moment it was asked.
    if check.motion.is_none() {
        check.one_moment = one_moment_of(v, &roots);
    }
    check.quantity = check.quantity.or(quantity);
    if let Some(b) = bound {
        for r in ranged_reads(v, &[b]) {
            if !check.ranged.contains(&r) {
                check.ranged.push(r);
            }
        }
        check.ranged.sort_by_key(|r| r.0);
    }
    // The claim's terms where its value is: at its witness, else at
    // the given params and t.
    let at = match &check.evidence {
        Evidence::ValueSpan(ValueSpan { witness: Some((at, _)), .. }) => {
            let mut c = v.ctx.clone();
            // Over a motion with t ranged, its moment is its param t's.
            c.t = check.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(at.t);
            for &(k, x) in &at.params {
                c.params[k as usize] = x;
            }
            Some(c.point())
        }
        _ if check.motion.as_ref().is_some_and(|m| m.motion.t.is_none()) => None,
        Evidence::Value(_) => {
            let mut c = v.ctx.clone();
            c.t = check.motion.as_ref().and_then(|m| m.motion.t).unwrap_or(c.t);
            Some(c)
        }
        _ => None,
    };
    let terms = at.map(|c| terms_at(v, terms, c)).unwrap_or_default();
    v.doc.wanted = wanted;
    Ok((check, terms))
}
