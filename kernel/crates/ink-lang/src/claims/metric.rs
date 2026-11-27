use super::*;

/// A claim's status on a metric, the one decision every path of
/// `metric_claim` takes: unknown where a reason not to decide was found
/// (`refused`: no value, a part left out that would change it, a grid that
/// cannot see the part or check it); else on the interval of the value's
/// band where it has one (a proven value takes none), or of the value
/// alone, widened on each side by the grid error measured (`margin`, if
/// any): it holds where both ends hold, fails where both fail, and is
/// otherwise unknown, saying, of a band, what dominates it.
fn metric_status(c: Comparison, measured: &Measured, margin: Option<f64>, refused: bool) -> Decided {
    let (Ok(x), false) = (&measured.value, refused) else { return Decided { status: Status::Unknown, crossing: None } };
    let e = margin.unwrap_or(0.0);
    let band = measured.band.as_ref().filter(|_| !measured.proven);
    let read = measured.read_error;
    let (lo, hi) = band.map_or(spread(*x, e, read), |b| (spread(b.low, e, read).0, spread(b.high, e, read).1));
    let status = match (point_status(c, lo), point_status(c, hi)) {
        (a, z) if a == z => a,
        _ => Status::Unknown,
    };
    let crossing = band.filter(|_| status == Status::Unknown).map(|b| {
        let num = sig;
        let parts = b.parts.iter().cloned().chain(margin.map(|e| ("grid".to_string(), e)));
        let dominant = parts.max_by(|a, b| a.1.total_cmp(&b.1)).map_or("its band".to_string(), |p| p.0);
        let widened = margin.map_or(String::new(), |e| format!(" with its grid error ({} to {} measured, grid error {})", num(b.low), num(b.high), num(e)));
        let narrows = b.narrows.iter().find(|n| n.0 == dominant).map_or("(res r) finer narrows each", |n| n.1.as_str());
        (
            format!("its measured band, {} to {}{widened}, crosses the bound; most of the band's width comes from {dominant}", num(lo), num(hi)),
            format!("to narrow the band: {narrows}; or move the bound past the band"),
        )
    });
    Decided { status, crossing }
}

/// `metric_status`'s answer: the status, and where a band crossing the bound left
/// it unknown, why and what would decide it.
struct Decided {
    status: Status,
    crossing: Option<(String, String)>,
}
/// A claim on a metric: `metric_claim_at`, once; or, of a form loaded from the model over a motion
/// whose metric is linear in those loads (`ModelView::moments_over`), once at every moment whose loads
/// differ, each on its own grids, and the least favourable of those answers: a no before a not settled
/// before a yes, and of two alike the one whose value reads worse. So no part of a motion answers
/// worse than the whole, by whatever its answer is decided on (its value, what it leaves out, what its
/// grid cannot see, a buckle).
pub(super) fn metric_claim(v: &mut ModelView, m: &FormMetric, c: Comparison, res: Option<f64>, accuracy: Option<f64>, excluding: &[&str]) -> Result<Check, InkError> {
    let moments = match v.moments_over(m) {
        Ok(moments) => moments,
        Err(e) if !m.form.is_empty() => return Ok(checked(Status::Unknown, Evidence::Unsolved(e)).estimate()),
        Err(e) => return Err(e),
    };
    if moments == 0 {
        return metric_claim_at(v, m, c, res, accuracy, excluding);
    }
    let rank = |s: Status| match s {
        Status::Fails => 0,
        Status::Unknown => 1,
        Status::Holds => 2,
    };
    let opts = RunOpts { res, tol: LADDER_TOL, accuracy, split_as: None, both_ways: false };
    let mut least: Option<(usize, Check)> = None;
    for k in 0..moments {
        v.pin(Some(k));
        let check = metric_claim_at(v, m, c, res, accuracy, excluding);
        v.pin(None);
        let check = check?;
        let less = match &least {
            None => true,
            Some((at, old)) => rank(check.status) < rank(old.status) || check.status == old.status && v.worse_moment(m, opts, k, *at),
        };
        if less {
            least = Some((k, check));
        }
    }
    Ok(least.expect("a moment").1)
}

/// A package's metric against its comparison: an estimate (its accuracy
/// moves the comparison, and guarantees nothing; exact where the metric is
/// proven, Measured::proven), and unknown where the
/// metric has no value (a failed solve; nothing to measure over). Forms
/// keep their own (res ...): the geometry's res is not theirs (cells a
/// hundredth of a millimetre would not fit a study); `res` overrides (a
/// shape's metric: the claim's).
///
/// A shape's metric on a grid (Package::shape_res) checks its grid error
/// on cells of 2h and 4h. Where the call has a time limit, those are
/// measured first, coarsest first: a limit passing before cells of h
/// finish leaves the finest that did, with their grid error from the next
/// coarser, an estimate that holds or fails only where that error cannot
/// move it past its bound (`Short`). Unstopped, it reads what it would
/// read after, and answers the same.
fn metric_claim_at(v: &mut ModelView, m: &FormMetric, c: Comparison, res: Option<f64>, accuracy: Option<f64>, excluding: &[&str]) -> Result<Check, InkError> {
    let opts = |res| RunOpts { res, tol: LADDER_TOL, accuracy, split_as: None, both_ways: false };
    let limited = if ink_kernel::stop::limited() { v.shape_res(m, res) } else { None };
    let mut coarser: Vec<(f64, Measured)> = Vec::new();
    if let Some(h) = limited {
        for k in [4.0, 2.0] {
            let measured = match v.measure(m, opts(Some(k * h))) {
                Ok(m) => m,
                Err(e) if !m.form.is_empty() => Measured::none(e.message),
                Err(e) => return Err(e),
            };
            coarser.push((k, measured));
        }
    }
    let mut measured = match v.measure(m, opts(res)) {
        Ok(m) => m,
        Err(e) if !m.form.is_empty() => return Ok(checked(Status::Unknown, Evidence::Unsolved(e)).estimate()),
        Err(e) => return Err(e),
    };
    // A value on a grid, holding or failing, is checked on a ladder of
    // grids (`ladder`), and answered on its finest; but for one whose read
    // is proven to bound the true value (`ReadError`) where the claim
    // answers alike at every value that allows, which a finer grid could
    // not change.
    let on_grid = |m: &Measured| !m.proven && m.unseen.is_none() && m.res.is_some() && m.value.is_ok();
    let proven_decides = |m: &Measured| {
        let allows = |x: f64, r: ReadError| I::new(r.span(x).0, r.span(x).1);
        let decides = |x: f64, r: ReadError| c.holds_on(allows(x, r)) || c.fails_on(allows(x, r));
        m.read_error.filter(|r| r.proven).is_some_and(|r| m.value.as_ref().is_ok_and(|&x| decides(x, r)))
    };
    let ladder = (on_grid(&measured) && !proven_decides(&measured)).then(|| ladder(v, m, accuracy, &measured));
    if let Some(Ok(l)) = &ladder {
        measured = l.finest.clone();
    }
    let head = v.doc.ext_form(Some(m.package), &m.form).map_or("form", |e| e.head);
    let num = sig;
    // A grid's value that saw the part.
    let seen = |m: &Measured| m.unseen.is_none().then(|| m.value.clone().ok()).flatten();
    let mut undecided: Option<(Vec<String>, Vec<String>)> = None;
    let mut status = match measured.value {
        Ok(x) => point_status(c, x),
        Err(_) => Status::Unknown,
    };
    // What the answer is decided on (`metric_status`): whether a reason not to
    // decide was found (a refusal), and the grid error a value is widened
    // by, where one was measured.
    let banded = |m: &Measured| !m.proven && m.band.is_some();
    let mut refused = false;
    let mut margin: Option<f64> = None;
    // Stopped before cells of h finished: the finest coarser that did.
    let stopped = match limited {
        Some(h) if measured.value.is_err() && ink_kernel::stop::reached() => {
            let finest = coarser.iter().rposition(|(_, c)| seen(c).is_some());
            let mut short = Short { res: h, error: None };
            if let Some(i) = finest {
                let (k, at) = coarser[i].clone();
                let x = seen(&at).expect("seen");
                let next = coarser[..i].last().filter(|(_, c)| seen(c).is_some());
                let compared = next.map(|(_, c)| grids_error(&[&at, c], v.metric_order(m)));
                short.error = compared.as_ref().and_then(|e| e.as_ref().ok()).map(|e| e.0);
                let said = match (&compared, next) {
                    (Some(Ok((e, _))), Some((k2, c))) => format!("grid error about {} (two grids: {} on cells of {} mm)", num(*e), num(seen(c).expect("seen")), k2 * h),
                    (Some(Err(why)), _) => format!("its grid error not known: {why}"),
                    _ => "the coarsest, so its grid error is not known".into(),
                };
                let why = format!("the call's time limit passed before cells of {h} mm finished: {} is its value on cells of {} mm, {said}", num(x), k * h);
                // (a read proven to bound the true value needs none where
                // the claim answers alike at every value it allows)
                let proven = proven_decides(&at);
                margin = short.error.filter(|_| !proven);
                refused = short.error.is_none() && !proven;
                // (a band, widened by that error, is decided on below)
                let band = banded(&at);
                status = match short.error {
                    _ if band && !refused => point_status(c, x),
                    _ if proven => point_status(c, x),
                    Some(e) if point_status(c, spread(x, e, at.read_error).0) == point_status(c, spread(x, e, at.read_error).1) => point_status(c, spread(x, e, at.read_error).0),
                    _ => Status::Unknown,
                };
                measured = at;
                match (status, short.error) {
                    _ if band && !refused => measured.warnings.push(why),
                    (Status::Unknown, Some(_)) => undecided.get_or_insert_with(Default::default).0.push(format!("{why}, more than its margin to the bound")),
                    (Status::Unknown, None) => undecided.get_or_insert_with(Default::default).0.push(why),
                    _ => measured.warnings.push(why),
                }
            } else {
                refused = true;
            }
            if status == Status::Unknown && (refused || !banded(&measured)) {
                let u = undecided.get_or_insert_with(Default::default);
                if u.0.is_empty() {
                    u.0.push(format!("the call's time limit passed before any cells finished (cells of {h} mm asked)"));
                }
                u.1.push(format!("the same call with a larger timeout_s goes on from the cache to cells of {h} mm"));
            }
            Some(short)
        }
        _ => None,
    };
    // What the metric leaves out, where counting it would answer
    // otherwise: the claim is undecided by what was measured, and says
    // what was left out, where, and how to decide it; unless it leaves
    // that out on purpose, when the answer stands and says so.
    // A proven value's claim is exact, unless it leaves out on purpose
    // what would change it.
    let mut exact = measured.proven;
    if let Ok(x) = measured.value.clone() {
        let measured_status = status;
        for ex in &measured.excluded {
            // Counted with its grid error, as the value is checked: what a
            // grid leaves out is read on that grid, no truer than the value
            // (but for a singular one, which has none to carry).
            let (e, side) = match &ladder {
                Some(Ok(l)) if !ex.singular => l.excluded_error(ex.key),
                _ => (0.0, 0.0),
            };
            // Counting adds what was left out, on one side of the value: its
            // error never carries it past the value. It is read as the value
            // is, so it may miss as its read may (`ReadError`, measured).
            let (ra, rb) = measured.read_error.filter(|r| !r.proven).map_or((ex.value, ex.value), |r| r.span(ex.value));
            let (below, above) = (ex.value - ra, rb - ex.value);
            let (lo, hi) = match side {
                s if s > 0.0 => ((ex.value - e - below).max(x), ex.value + e + above),
                s if s < 0.0 => (ex.value - e - below, (ex.value + e + above).min(x)),
                _ => (ra, rb),
            };
            let counted = [ex.value, lo, hi].map(|v| point_status(c, v)).into_iter().find(|&s| s != measured_status).unwrap_or(measured_status);
            if counted == measured_status && !(ex.possible && ex.value.is_infinite()) {
                continue;
            }
            let num = |x: f64| if e > 0.0 { format!("{}, its grid error about {}", num(x), num(e)) } else { num(x) };
            let at = ex.at.map_or(String::new(), |p| format!(" at [{:.3}, {:.3}, {:.3}]", p[0], p[1], p[2]));
            // The verdict with it counted is its value's own; where only an end of its error's span
            // differs, counting it may change the verdict, and is not said to.
            let effect = if ex.possible {
                "Counting it may change the verdict".to_string()
            } else if ex.singular && (ex.value == 0.0 || ex.value.is_infinite()) {
                // A value with no bound (a stress's, its safety factor 0): no number is counted, the bound is.
                format!("Counted, it is unbounded, and the verdict is {}", counted.as_str())
            } else if point_status(c, ex.value) == counted {
                format!("Including its value ({}{at}), the verdict is {}", num(ex.value), counted.as_str())
            } else {
                format!("Including its value ({}{at}), which may lie from {} to {}, the verdict may be {}", num(ex.value), sig(lo), sig(hi), counted.as_str())
            };
            if excluding.contains(&ex.key) {
                exact = false;
                measured.warnings.push(format!("this claim {}, leaving out on purpose (excluding {}): {}. {effect}", measured_status.as_str(), ex.key, ex.what));
            } else {
                status = Status::Unknown;
                refused = true;
                let u = undecided.get_or_insert_with(Default::default);
                u.0.push(format!("it {} on a value that leaves out {}. {effect}", measured_status.as_str(), ex.what));
                let offered = v.doc.registry.exclusions(m.metric).iter().any(|(key, _)| *key == ex.key);
                // Where it goes: inside the claim, after its bound.
                u.1.push(if offered {
                    format!(
                        "{}; or leave it out on purpose, inside the claim after its bound: (require NAME ({} {} ... (excluding {}))), and the claim is then on what is read without it",
                        ex.decide, m.metric, m.form, ex.key
                    )
                } else {
                    ex.decide.clone()
                });
            }
        }
    }
    // A grid's answer, holding or failing, stands only where its grid saw
    // the part, and its grid error, estimated from coarser grids, cannot
    // move it past its bound (`grid_error`); coarser grids that do not see
    // the part, or solve nothing, cannot check it, and leave it unknown.
    let mut grid: Option<String> = None;
    let mut grid_err: Option<(f64, &'static str)> = None;
    // What its value on the grid says, before a refusal: a banded value
    // (a buckling factor's) may fail there.
    let read_as = match status {
        Status::Holds => "it holds",
        Status::Fails => "it fails",
        _ => "its band crosses the bound",
    };
    // A solve that could not see the part decides nothing either way: a gap
    // glued shut joins what the part keeps apart, and a thin feature held
    // too stiff reads too little motion, which fails a true `>` as surely as
    // it holds a false `<`. A value and a band alike (a band crossing the
    // bound, nothing).
    let reads = match banded(&measured) {
        true if refused => Status::Unknown,
        true => metric_status(c, &measured, None, false).status,
        false => status,
    };
    let where_res = if !m.form.is_empty() { format!(" in ({head} {} ...)", m.form) } else { String::new() };
    if (reads != Status::Unknown || banded(&measured) && !refused)
        && let Some(why) = measured.unseen.clone()
    {
        status = Status::Unknown;
        refused = true;
        let u = undecided.get_or_insert_with(Default::default);
        u.0.push(format!("{read_as} on a grid that cannot see the part: {why}"));
        u.1.push(format!(
            "refine (res r){where_res}: to under four times the width of a gap or wall it cannot see (split to an eighth, two of its cells across it), or as it says for a curve; or close the gap, thicken the wall or round the curve in the model"
        ));
    }
    if (status != Status::Unknown || banded(&measured) && !refused)
        && let Some(l) = ladder
    {
        match l {
            Err((h, why)) => {
                status = Status::Unknown;
                refused = true;
                let u = undecided.get_or_insert_with(Default::default);
                u.0.push(format!(
                    "{read_as} on cells of {h} mm, but no three grids, each twice as coarse as the last, the coarsest no coarser than four times these, all read the part: {why}"
                ));
                // One number to refine to: the curve's own where the reason gives it, else half these cells.
                u.1.push(if why.contains("refine (res r) to") {
                    if m.form.is_empty() {
                        "refine it as that says, or give the part walls and rounds the cells can hold".to_string()
                    } else {
                        format!("refine it as that says, in ({head} {} ...), or give the part walls and rounds the cells can hold", m.form)
                    }
                } else if m.form.is_empty() {
                    format!("refine it, (res r) at {} mm or finer, or give the part walls and rounds the cells can hold", h / 2.0)
                } else {
                    format!("refine it, (res {}) in ({head} {} ...) or finer, or give the part walls and rounds the cells can hold", h / 2.0, m.form)
                });
            }
            Ok(l) => {
                let ([x, _, _], [h, _, _]) = (l.values, l.cells);
                let [x2, x4] = [l.values[1], l.values[2]];
                let [_, h2, h4] = &l.named;
                let (e, how) = l.error();
                grid_err = Some((e, how));
                margin = Some(e);
                let read = format!("{how}: {} on cells of {h2}, {} on cells of {h4}", num(x2), num(x4));
                grid = Some(match e.is_finite() {
                    true => format!("its grid error about {} ({read})", num(e)),
                    false => format!("no grid error read from its grids ({read})"),
                });
                let (lo, hi) = spread(x, e, measured.read_error);
                if !banded(&measured) && (point_status(c, lo) != status || point_status(c, hi) != status) {
                    status = Status::Unknown;
                    let u = undecided.get_or_insert_with(Default::default);
                    // The read's error on the side whose end crosses the bound.
                    let off = measured.read_error.filter(|r| !r.proven).map_or(String::new(), |r| {
                        let (a, b) = r.span(x);
                        let sides: Vec<String> = [(point_status(c, lo) != status, x - a, "under"), (point_status(c, hi) != status, b - x, "over")]
                            .into_iter()
                            .filter(|s| s.0)
                            .map(|(_, d, way)| format!("{} {way} it", num(d)))
                            .collect();
                        format!(", and the value may lie up to {} ({})", sides.join(" or "), r.how)
                    });
                    if e.is_finite() {
                        u.0.push(format!("{read_as} at {} on cells of {h} mm, but its grid error is about {} ({read}){off}, more than its margin to the bound", num(x), num(e)));
                        u.1.push(if m.form.is_empty() {
                            format!("refine it, (res r) at {} mm or finer, until its value moves by less than its margin", h / 2.0)
                        } else {
                            format!("refine it, (res {}) in ({head} {} ...) or finer, until its value moves by less than its margin", h / 2.0, m.form)
                        });
                    } else {
                        // How far from converging, in its own numbers: each change against the last, and the margin.
                        let (d1, d2) = (x - x2, x2 - x4);
                        let swing = if d1 * d2 < 0.0 { ", the other way" } else { "" };
                        let margin = (x - c.bound).abs();
                        u.0.push(format!(
                            "{read_as} at {} on cells of {h} mm, but no grid error can be read from its grids ({read}): its last change, {}, is {} times the one before it, {}{swing}, \
                             and three grids converge only where each change is smaller than the last; its margin to the bound is {}, {} times its last change. {}",
                            num(x),
                            num(d1.abs()),
                            num(d1.abs() / d2.abs().max(f64::MIN_POSITIVE)),
                            num(d2.abs()),
                            num(margin),
                            num(margin / d1.abs().max(f64::MIN_POSITIVE)),
                            match l.unbounded() {
                                Some(why) => format!("It may grow without bound: {why}"),
                                None => "A value that moves away as the cells shrink may go on moving".to_string(),
                            }
                        ));
                        let refine = if m.form.is_empty() { format!("(res r) at {} mm", h / 2.0) } else { format!("(res {}) in ({head} {} ...)", h / 2.0, m.form) };
                        u.1.push(format!(
                            "refine it, {refine} (about 8 times these cells), and it is read on cells of {}, {h} and {} mm: settled there if its next change is smaller than {}",
                            h / 2.0,
                            2.0 * h,
                            num(d1.abs())
                        ));
                    }
                }
            }
        }
    }
    // One decision for every path: a refusal is unknown; else the value's
    // band (or the value), widened on each side by the grid error measured.
    let decided = metric_status(c, &measured, margin, refused);
    if decided.status == Status::Unknown
        && !refused
        && let Some(why) = decided.crossing
    {
        let u = undecided.get_or_insert_with(Default::default);
        u.0.push(why.0);
        u.1.push(why.1);
    }
    debug_assert!(banded(&measured) || decided.status == status, "decide {:?} against {status:?}", decided.status);
    status = decided.status;
    let undecided = undecided.map(|(why, how)| (why.join("; "), how.join("; ")));
    let quantity = v.doc.registry.metric_quantity(m.metric);
    let held = measured.held.clone();
    // What an estimate rests on: its cells, and its grid check where it had one.
    let estimate = (!exact).then(|| {
        let on = measured.res.map_or("a numerical value, not a bound".to_string(), |h| format!("a value on cells of {h} mm, not a bound"));
        let read = measured.read_error.filter(|_| status != Status::Unknown);
        let proven = read.filter(|r| r.proven && proven_decides(&measured));
        let proven = proven.map(|r| format!("{on}; bounding the true value ({}), so it {} at any value that allows, with no grid check", r.how, status.as_str()));
        let measured_read = read.filter(|r| !r.proven).map_or(String::new(), |r| format!("; and its read's error ({})", r.how));
        proven.or_else(|| grid.map(|g| format!("{on}; {g}{measured_read}"))).unwrap_or(on)
    });
    let check = Check { held, ..checked(status, Evidence::Metric { measured, head, form: m.form.clone(), quantity, undecided, stopped, estimate, grid_error: grid_err }) };
    Ok(if exact { check } else { check.estimate() })
}
