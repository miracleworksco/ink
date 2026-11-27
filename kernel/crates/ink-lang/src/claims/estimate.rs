use super::*;

/// A part of an integral claim: its shape, its integrands (None meaning 1),
/// and its body's own frame, where it is a placed body.
pub(super) type Part = (u32, Vec<Option<u32>>, Option<Placed>);

/// Integrals summed over parts, estimated, and the claim's status on them.
pub(super) struct Summed {
    /// One per integrand: the sum over the parts.
    pub(super) sums: Vec<f64>,
    /// The coarsest cells a part stopped at (0: none had material).
    pub(super) res: f64,
    /// Every part is proven bounded (else the sums are of the rest).
    pub(super) whole: bool,
    pub(super) status: Status,
}

/// Sums over `parts` of each integrand's integral over the part's shape,
/// estimated at the view's params and moment; the sums, the coarsest
/// cells a part stopped at, and whether every part is proven bounded (else
/// the sums are of the rest). A placed body is integrated in its own frame
/// (frame::in_frame), times 1 / |det R|: copies of a body are one
/// question there, answered once (the cache).
fn estimate_parts(v: &mut ModelView, parts: &[Part], res: Option<f64>, decided: impl Fn(&[f64], &[f64]) -> bool) -> (Vec<f64>, Vec<f64>, Vec<f64>, f64, bool) {
    let n = parts.first().map_or(0, |p| p.1.len());
    // Each part as integrated (its frame's shape and integrands), alike
    // parts once, with how many and their scale.
    let mut alike: Vec<((u32, Vec<Option<u32>>), f64)> = Vec::new();
    for (shape, gs, frame) in parts {
        let local = frame.and_then(|f| in_frame(v, &f, gs)).filter(|f| v.whole(f.shape));
        let (key, scale) = match local {
            Some(f) => ((f.shape, f.gs), f.inv_det.mid()),
            None => ((*shape, gs.clone()), 1.0),
        };
        match alike.iter_mut().find(|a| a.0 == key) {
            Some(a) => a.1 += scale,
            None => alike.push((key, scale)),
        }
    }
    let whole = alike.iter().all(|((s, _), _)| v.whole(*s));
    let mut each = Vec::new();
    for ((shape, gs), scale) in alike {
        let Some((lo, hi, _)) = v.bounds(shape).filter(|_| v.whole(shape)) else { continue };
        let tapes: Vec<Option<Tape>> = gs.iter().map(|g| g.map(|g| Tape::new(&v.doc.g, g))).collect();
        each.push((v.clipped(shape), tapes, lo, hi, scale));
    }
    // Measured as `ink mass` measures a mass, to 16 and 4 times its
    // accuracy first, stopping where the bounds proven so far decide the
    // claim (at `res` once, where given); one part at a time: each
    // measurement is parallel, and holds its boxes until it ends.
    let steps = measure_steps(res);
    let (mut sums, mut coarsest): (Vec<f64>, f64) = (vec![0.0; n], 0.0);
    let (mut lo, mut hi) = (vec![0.0; n], vec![0.0; n]);
    for (i, &until) in steps.iter().enumerate() {
        (lo, hi) = (vec![0.0; n], vec![0.0; n]);
        (sums, coarsest) = (vec![0.0; n], 0.0);
        for (tape, tapes, blo, bhi, scale) in &each {
            let gs: Vec<Option<&Tape>> = tapes.iter().map(Option::as_ref).collect();
            let m = query::measured(tape, &gs, false, &v.ctx, *blo, *bhi, until);
            for k in 0..n {
                (sums[k], lo[k], hi[k]) = (sums[k] + scale * m.sums[k], lo[k] + scale * m.lo[k], hi[k] + scale * m.hi[k]);
            }
            coarsest = coarsest.max(m.res);
        }
        if i + 1 == steps.len() || (whole && decided(&lo, &hi)) {
            break;
        }
    }
    (sums, lo, hi, coarsest, whole)
}

/// How far an integral claim's integrals are measured, in turn, until their
/// proven bounds decide it: at `res` once, where given; else to 16 and 4
/// times `ink mass`'s accuracy, and on to a sixteenth of it where the
/// bounds have not decided (the estimate no longer decides, and a mean's
/// bounds, a ratio of two integrals', are wider than either's).
pub(crate) fn measure_steps(res: Option<f64>) -> Vec<query::Until> {
    match res {
        Some(r) => vec![query::Until::Res(r)],
        None => [16.0, 4.0, 1.0, 0.25, 0.0625, 0.015625].map(|k| query::Until::Accuracy(k * query::ACCURACY)).to_vec(),
    }
}

/// Values of a span an estimate is taken at, evenly spaced, its ends
/// included: a motion's moments (motion::measured), a ranged param's
/// values where it is the only one read (`estimate_scope`). One count for
/// both: on the catalog's, tests' and examples' 247 motion models, 17
/// moments and 34 gave every claim the same answer, 17 in 85% of the
/// measured claims' time; 33 values of a band doubled a holding band
/// estimate's (examples/desk-lamp.ink's swung claims, 1.7 s to 3.3 s).
pub const SAMPLES: usize = 17;

/// The most combinations of declared cases one `ink require` asks.
pub const MOST_CASES: usize = 256;

/// Why a param's `(cases ...)` are refused where a combination of the
/// cases declares `here` of them and the params' values `first`, and how
/// to make the declarations agree.
pub fn cases_count_changes(name: &str, first: usize, here: usize) -> String {
    let cases = |k: usize| if k == 1 { "1 case".to_string() } else { format!("{k} cases") };
    format!(
        "{name} declares {} at this combination but {first} at the params' values: a param's (cases ...) name as many values at every combination, \
         their values may differ. Declare it with every case at each, a value chosen by the other params where it differs ((cases 1 (if (== a 1) 2 1)), say), \
         or check each as its own file",
        cases(here)
    )
}

/// Why a param's `(cases ...)` are refused where it is declared at some
/// combinations of the cases and not others (`where_`: where it is and
/// is not), and what to do instead.
pub fn cases_not_everywhere(name: &str, where_: &str) -> String {
    format!(
        "{name}'s (cases ...) are refused: {name} {where_}, so not every combination of the cases declares it (a `for` or `if` over a param makes it). \
         Cases go on params every combination declares, and a param declared at some takes its value there. To check other values of it, \
         name them as cases of a param every combination declares (the count the `for` reads, say), or check each as its own file"
    )
}

/// SAMPLES values of `span`, evenly spaced, its ends included (one, where
/// it is a single value).
pub fn samples(span: (f64, f64)) -> Vec<f64> {
    spaced(span, SAMPLES)
}

/// n values of `span`, evenly spaced, its ends included (one, where it is
/// a single value or n is 1).
fn spaced((lo, hi): (f64, f64), n: usize) -> Vec<f64> {
    if lo >= hi || n < 2 {
        return vec![lo];
    }
    (0..n).map(|i| if i + 1 == n { hi } else { lo + (hi - lo) * i as f64 / (n - 1) as f64 }).collect()
}

/// `estimate_parts` over the claim's scope, with the claim's status on the
/// sums (`status`): at the view's params where its fields read no ranged
/// param; else at a grid of their values, each param's span at evenly
/// spaced values, ends included (SAMPLES of one; of n, about SAMPLES in
/// all, at least both ends of each). It fails at the
/// first value where it fails (the sums there), holds where it holds at
/// every one (the sums at the last), else is unknown: an estimate, which a
/// violation between the values escapes. What the fields read is what is
/// sampled: of a part in its frame, its shape and integrands there, and
/// its placement.
pub(super) fn estimate_scope(v: &mut ModelView, parts: &[Part], res: Option<f64>, status: impl Fn(&[f64]) -> Status) -> (Summed, Option<Sampled>) {
    let mut fields: Vec<u32> = Vec::new();
    for (s, gs, f) in parts {
        match f.and_then(|f| Some((in_frame(v, &f, gs)?, f))).filter(|(i, _)| v.whole(i.shape)) {
            // In its frame: its shape and integrands there, and its
            // placement where its determinant varies (a stretch).
            Some((i, f)) => {
                if i.inv_det.hi - i.inv_det.lo > 1e-9 * i.inv_det.hi {
                    fields.extend(f.at.r.iter().flatten());
                }
                fields.extend(std::iter::once(i.shape).chain(i.gs.into_iter().flatten()));
            }
            None => fields.extend(std::iter::once(*s).chain(gs.iter().flatten().copied())),
        }
    }
    let read = ranged_reads(v, &fields);
    let summed = |v: &mut ModelView| {
        // Decided where the status is the same at every corner of the
        // integrals' proven bounds (the first and the rest apart: a mean's
        // weight and moments).
        let decided = |lo: &[f64], hi: &[f64]| {
            let corner = |a: bool, b: bool| status(&(0..lo.len()).map(|k| if (k == 0 && a) || (k > 0 && b) { hi[k] } else { lo[k] }).collect::<Vec<f64>>());
            let st = corner(false, false);
            st != Status::Unknown && [(true, true), (true, false), (false, true)].iter().all(|&(a, b)| corner(a, b) == st)
        };
        // Decided by the proven bounds alone: the estimate is the value shown, and decides nothing (it has
        // no bound of its own).
        let (sums, lo, hi, res, whole) = estimate_parts(v, parts, res, decided);
        let status = if whole && decided(&lo, &hi) { status(&lo) } else { Status::Unknown };
        Summed { sums, res, whole, status }
    };
    if read.is_empty() {
        return (summed(v), None);
    }
    let per = ((SAMPLES as f64).powf(1.0 / read.len() as f64).floor() as usize).max(2);
    let each: Vec<Vec<f64>> = read.iter().map(|&(_, s)| spaced(s, per)).collect();
    let values: usize = each.iter().map(Vec::len).product();
    let mut result: Option<Summed> = None;
    for k in 0..values {
        let mut rest = k;
        let at: Vec<(u32, f64)> = read
            .iter()
            .zip(&each)
            .map(|(&(p, _), xs)| {
                let x = xs[rest % xs.len()];
                rest /= xs.len();
                (p, x)
            })
            .collect();
        let ctx = at.iter().fold(v.ctx.clone(), |c, &(p, x)| c.with_param(p, (x, x)));
        let s = v.at_point(ctx, |v| summed(v));
        if s.status == Status::Fails {
            return (s, Some(Sampled { values: k + 1, failing: Some(at) }));
        }
        // An undefined sample prevents a hold, but a later failure still has a witness.
        if result.as_ref().is_none_or(|r| r.status != Status::Unknown) {
            result = Some(s);
        }
    }
    (result.expect("a value"), Some(Sampled { values, failing: None }))
}

/// `(mean target h (weight w) cmp v)`: the integral of w h over the parts
/// over that of w, for each h (a vec), each compared with v, estimated
/// (`estimate_scope`). Where the weight's integral is not positive there
/// is no mean: unknown.
pub(super) fn mean(v: &mut ModelView, parts: &[(u32, u32, Vec<u32>)], frames: &[Option<Placed>], c: Comparison, res: Option<f64>) -> Check {
    let split: Vec<Part> = parts
        .iter()
        .enumerate()
        .map(|(k, (s, w, hs))| {
            let moments = hs.iter().map(|&h| Some(v.doc.g.mul(*w, h)));
            (*s, std::iter::once(Some(*w)).chain(moments).collect(), frames.get(k).copied().flatten())
        })
        .collect();
    let each = |i: &[f64]| -> Vec<(f64, Status)> { i[1..].iter().map(|a| a / i[0]).map(|m| (m, if i[0] > 0.0 { point_status(c, m) } else { Status::Unknown })).collect() };
    let (s, sampled) = estimate_scope(v, &split, res, |i| all_of(each(i).into_iter().map(|m| m.1)));
    let means = each(&s.sums);
    let sampled_scope = sampled.is_some();
    let mut c = checked(s.status, Evidence::Mean { weight: s.sums[0], means, res: s.res, sampled });
    if sampled_scope {
        c = c.estimate();
    }
    if s.sums[0] <= 0.0 {
        c.note = Some("the weight's integral over the target is not positive (a target with no material has no mean): no mean to compare");
    }
    c.unbounded = !s.whole;
    c
}

/// An answer by a method that guarantees it (`check` writes its
/// accuracy); `Check::estimate` makes one an estimate.
pub(crate) fn checked(status: Status, evidence: Evidence) -> Check {
    Check {
        status,
        accuracy: Some(0.0),
        asked: None,
        quantity: None,
        evidence,
        in_window: None,
        unbounded: false,
        ranged: Vec::new(),
        note: None,
        inexact: Default::default(),
        motion: None,
        one_moment: None,
        held: None,
    }
}
