//! Claims' answers as JSON (`ink assert`'s rows, and wherever else a
//! claim is asked as it asks one: an optimization's requirements): its
//! status, how far off it could be, and what decided it.

use crate::claims::{Check, Evidence, Status};
use crate::command::{bounds, inexact_fields, measured, n, v3};
use ink_core::V3;
use ink_kernel::query::{Ended, Verdict};
use ink_kernel::ranged::Proof;
use serde_json::{Value as J, json};

/// How a claim's hint says to look finer.
pub const ASK_FINER: &str = "A smaller (accuracy x) looks finer";

/// The names of the document's params, and of those the answer is about
/// (the ranged params its fields read: a witness names only these).
pub struct Names<'a> {
    pub all: &'a [String],
    pub read: &'a [u32],
}

/// Ranged params' values (a witness's) that the answer is about, by name.
fn params_json(names: &Names, ps: &[(u32, f64)]) -> J {
    J::Object(ps.iter().filter(|(k, _)| names.read.contains(k)).map(|(k, v)| (names.all[*k as usize].clone(), n(*v))).collect())
}

/// Ranged params' spans that the answer is about, by name: those a model
/// named, not a (ranged) no def names (std's hidden makes one per view),
/// which spans what it was written with whatever the run.
fn spans_json(names: &Names, ps: &[(u32, (f64, f64))]) -> J {
    J::Object(
        ps.iter()
            .filter(|(k, _)| names.read.contains(k) && !crate::unnamed_range(&names.all[*k as usize]))
            .map(|(k, (a, b))| (names.all[*k as usize].clone(), bounds(*a, *b)))
            .collect(),
    )
}

/// An answer's accuracy (`status` its status): a number only where it
/// holds or fails guaranteed (`asked`, None for an estimate); an unknown
/// gives the accuracy it was undecided at as `accuracy_asked`.
fn accuracy_fields(r: &mut J, status: Status, asked: Option<f64>) {
    let o = r.as_object_mut().expect("an answer is an object");
    o.remove("accuracy_asked");
    match (status, asked) {
        (Status::Unknown, Some(t)) => {
            o.insert("accuracy".into(), J::Null);
            o.insert("accuracy_asked".into(), n(t));
        }
        (_, t) => {
            o.insert("accuracy".into(), t.map_or(J::Null, n));
        }
    }
}

/// A proof's fields after its status: its accuracy (`accuracy_fields`:
/// the one it was asked to, 0 exactly), the witness (and its ranged
/// params' values) or undecided cells, with what can decide them and the
/// accuracy it holds to where the answer is the proof's own (`attained`:
/// Check::attained); and the cell size.
fn proof_fields(r: &mut J, p: &Proof, names: &Names, attained: Option<(f64, bool)>) {
    accuracy_fields(r, Status::of(p.status()), Some(p.tol));
    match &p.verdict {
        Verdict::Proven => {}
        Verdict::Disproven { witness, value, params, .. } => {
            r["witness"] = v3(*witness);
            r["value"] = n(*value);
            if !params.is_empty() {
                r["params"] = params_json(names, params);
            }
        }
        Verdict::Unknown { cells, volume, ended, cell, .. } => {
            undecided_fields(r, *cells, *volume, *cell);
            let holds = attained_fields(r, attained, "no point of what violates the claim is deeper than that");
            r["hint"] = json!(undecided_hint(p.tol, *ended, &holds, !names.read.is_empty()));
        }
    }
    r["res"] = n(p.res);
}

/// What an unknown still holds to (Check::attained), as its answer gives
/// it: `accuracy_attained`, or `accuracy_attained_in_window` where only
/// the window was examined; and what its hint says of it, `means` being
/// what holding to it says.
fn attained_fields(r: &mut J, attained: Option<(f64, bool)>, means: &str) -> String {
    match attained {
        Some((a, false)) => {
            r["accuracy_attained"] = n(a);
            format!(". It holds to {} mm: {means}", up(a))
        }
        Some((a, true)) => {
            r["accuracy_attained_in_window"] = n(a);
            format!(". In the 1 m cube at the origin it holds to {} mm ({means} there), not past it", up(a))
        }
        None => String::new(),
    }
}

/// A bound for a hint's prose: x rounded up to two significant figures,
/// so what it bounds stays bounded.
fn up(x: f64) -> String {
    if !(x.is_finite() && x > 0.0) {
        return format!("{x}");
    }
    let e = x.log10().floor() as i32 - 1;
    let digits = (-e).max(0) as usize;
    let mut k = (x / 10f64.powi(e)).ceil();
    loop {
        let s = format!("{:.*}", digits, k * 10f64.powi(e));
        if s.parse::<f64>().is_ok_and(|v| v >= x) {
            return s;
        }
        k += 1.0;
    }
}

/// Undecided cells: how many, their volume, and one of them (the largest:
/// its centre and size).
fn undecided_fields(r: &mut J, cells: usize, volume: f64, cell: Option<[V3; 2]>) {
    r["undecided_cells"] = json!(cells);
    r["undecided_volume"] = n(volume);
    if let Some([lo, hi]) = cell {
        r["undecided_cell"] = json!({"centre": v3(std::array::from_fn(|k| 0.5 * lo[k] + 0.5 * hi[k])), "size": v3(std::array::from_fn(|k| hi[k] - lo[k]))});
    }
}

/// What can decide an emptiness proof left undecided, asked to accuracy
/// `tol` (0: exactly, where the claim takes none or its shape is thin), by
/// why its search ended (`ended`): a box at its finest cells stayed
/// undecided (the first the search, depth first, reached); it tested the
/// boxes it may; or the time limit; `holds`, what it says of the
/// accuracy it holds to (attained_fields); `ranged`, where the claim reads
/// ranged params (each is halved at most PARAM_SPLITS times in a box).
fn undecided_hint(tol: f64, ended: Ended, holds: &str, ranged: bool) -> String {
    let params = narrower(ranged);
    let first = "undecided_cell, the first cell the search refined to its finest size, stays undecided";
    match (tol > 0.0, ended) {
        (true, Ended::Finest) => format!(
            "{first} at half the claim's accuracy. What violates it there is about as deep as the accuracy, or too thin to hold a cube twice that size. The search ends there, the rest unrefined{holds}. {ASK_FINER}{params}"
        ),
        (false, Ended::Finest) => format!(
            "{first} at the finest cells (res), asked exactly. Its bounds there neither settle the claim nor show a point failing it, so the search ends there, the rest unrefined{holds}. \
             A finer (res r) decides boxes whose bounds narrow with their size. No res decides a claim whose value reaches its bound (surfaces that touch), or one reading a slope (diff) \
             where its formula switches branch (min, max, abs, clamp): there the bounds span every rate between the branches', at any size{params}"
        ),
        (_, Ended::Limit) => format!(
            "the search stopped at its limit of {} boxes, or would have passed four times that at its rate (surfaces within the accuracy of each other over a large area){holds}. A larger (accuracy x) decides it sooner{params}",
            ink_kernel::query::BOXES
        ),
        (_, Ended::Stopped) => STOPPED.into(),
    }
}

/// What narrows a search over ranged params (`ranged`: it reads some).
fn narrower(ranged: bool) -> String {
    match ranged {
        true => format!(". Narrower ranged params (each halved at most {} times in a box) decide it sooner", ink_kernel::scope::PARAM_SPLITS),
        false => String::new(),
    }
}

/// Why a search the call's time limit stopped is undecided, and what to do.
const STOPPED: &str = "the call's time limit (timeout_s) stopped the search. undecided_cells are the boxes it left (for a sweep or an integral, its last full pass). A larger timeout_s checks it again, and fewer ranged params, narrower ranges or a smaller region cost less";

/// Where an estimate over ranged params was taken (claims::Sampled): how
/// many values of them (evenly spaced, ends included), and those where it
/// fails.
fn sampled_fields(r: &mut J, sampled: Option<&crate::claims::Sampled>, names: &Names) {
    let Some(s) = sampled else { return };
    r["sampled_values"] = json!(s.values);
    if let Some(at) = &s.failing {
        r["params"] = params_json(names, at);
    }
}

/// A proof's answer (`attained`: Check::attained, where it is the
/// question's own).
pub fn proof_json(p: &Proof, names: &Names, attained: Option<(f64, bool)>) -> J {
    let mut r = json!({"status": Status::of(p.status()).as_str()});
    proof_fields(&mut r, p, names, attained);
    r
}

/// A claim's answer as JSON: its status, how far off it could be (its
/// accuracy: a number, guaranteed; null, an estimate), and what decided
/// it.
pub fn check_json(c: &Check, all: &[String]) -> J {
    let read: Vec<u32> = c.ranged.iter().map(|r| r.0).collect();
    let names = &Names { all, read: &read };
    let mut r = json!({"status": c.status.as_str()});
    match &c.evidence {
        Evidence::Empty(p) => {
            proof_fields(&mut r, p, names, c.attained());
            // How deep the witness lies in what violates the claim (both
            // shapes, for disjoint): the first point found, not the deepest.
            if let Verdict::Disproven { value, .. } = &p.verdict {
                r["depth"] = n(-value);
            }
        }
        Evidence::Vacuous(note) => r["note"] = json!(note),
        Evidence::Unexamined => {}
        Evidence::ClearanceSearch { verdict, pair, res } => {
            match verdict {
                Verdict::Proven => {}
                Verdict::Disproven { witness, value, params, .. } => {
                    r["witness"] = v3(*witness);
                    r["value"] = n(*value);
                    if !params.is_empty() {
                        r["params"] = params_json(names, params);
                    }
                }
                Verdict::Unknown { cells, volume, ended, cell, .. } => {
                    undecided_fields(&mut r, *cells, *volume, *cell);
                    // What it holds to (no pair closer than the clearance
                    // less that), over all of the shapes or in the window
                    // (Check::attained).
                    let holds = attained_fields(&mut r, c.attained(), "no pair of points is closer than the clearance less that");
                    let near = "the gap is too close to the clearance to decide, or the shapes' bounds cannot separate them";
                    let smaller = "If the fit allows, claim a slightly smaller clearance: coarser cells decide it";
                    let finest = up((*cell).map_or(*res, |[lo, hi]| ink_kernel::query::max_extent(lo, hi)));
                    r["hint"] = json!(match ended {
                        Ended::Limit => format!(
                            "The search stopped at its limit of {} boxes, or would have passed four times that at its rate. On the cells it reached, {near}{holds}. {smaller}{}",
                            ink_kernel::query::BOXES,
                            narrower(!names.read.is_empty())
                        ),
                        Ended::Stopped => STOPPED.into(),
                        _ => format!("At the search's finest cells, about {finest} mm across, {near}{holds}. {ASK_FINER}. {smaller}"),
                    });
                }
            }
            if let Some((q, p)) = pair {
                r["witness_pair"] = json!({"in_second": v3(*q), "in_first": v3(*p)});
                // How far apart the pair is (value: rounded up): the shapes
                // come this close, so their smallest gap is no larger (the
                // search stops at the first pair closer than the clearance,
                // not the closest). Only a violation has one.
                if let Verdict::Disproven { value, .. } = verdict {
                    r["value"] = n(*value);
                    r["gap"] = n((0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt());
                }
            }
            r["res"] = n(*res);
        }
        Evidence::Integral { integrals, res, sampled } => {
            let one = |(x, s): &(f64, Status)| json!({"estimate": n(*x), "status": s.as_str()});
            match integrals.as_slice() {
                [i] => {
                    // The claim's number, where the Results table says a
                    // number is: one field's integral (several: each under `integrals`).
                    r["value"] = n(i.0);
                    r["integral"] = one(i);
                }
                all => r["integrals"] = J::Array(all.iter().map(one).collect()),
            }
            sampled_fields(&mut r, sampled.as_ref(), names);
            r["res"] = n(*res);
        }
        Evidence::Mean { weight, means, res, sampled } => {
            let one = |(x, s): &(f64, Status)| json!({"estimate": n(*x), "status": s.as_str()});
            match means.as_slice() {
                [m] => {
                    // The claim's number, where the Results table says a
                    // number is: one field's mean (several: each under `means`).
                    r["value"] = n(m.0);
                    r["mean"] = one(m);
                }
                all => r["means"] = J::Array(all.iter().map(one).collect()),
            }
            r["weight"] = json!({"estimate": n(*weight)});
            sampled_fields(&mut r, sampled.as_ref(), names);
            r["res"] = n(*res);
        }
        Evidence::Metric { measured: v, head, form, undecided, .. } => {
            match &v.value {
                Ok(x) => r["value"] = n(*x),
                Err(why) => {
                    r["value"] = J::Null;
                    r["unavailable"] = json!(why);
                }
            }
            if let Some(at) = v.at {
                r["at"] = v3(at);
            }
            if !form.is_empty() {
                r[*head] = json!(form);
            }
            if let Some(h) = v.res {
                r["res"] = n(h);
            }
            // What it leaves out: each with its value counted (null:
            // unbounded) and where.
            let excludes: Vec<J> = v
                .excluded
                .iter()
                .map(|ex| {
                    let mut e = json!({"key": ex.key, "value": if ex.value.is_finite() { n(ex.value) } else { J::Null }, "what": ex.what});
                    if let Some(at) = ex.at {
                        e["at"] = v3(at);
                    }
                    e
                })
                .collect();
            if !excludes.is_empty() {
                r["excludes"] = J::Array(excludes);
            }
            if let Some((why, how)) = undecided {
                r["undecided"] = json!(why);
                r["hint"] = json!(how);
            }
            if !v.warnings.is_empty() {
                r["warnings"] = json!(v.warnings);
            }
        }
        Evidence::Value(v) => r["value"] = n(*v),
        Evidence::ValueSpan(result) => {
            r["bounds_over_scope"] = bounds(result.range.0, result.range.1);
            r["pieces"] = json!(result.pieces);
            if let Some((at, v)) = &result.witness {
                if !at.params.is_empty() {
                    r["params"] = params_json(names, &at.params);
                }
                r["value"] = n(*v);
            }
        }
        Evidence::NotCovered(why) => r["note"] = json!(why),
        Evidence::Stepped => {}
        Evidence::Outside(why) => r["unavailable"] = json!(why),
        Evidence::Stopped => r["hint"] = json!("not checked: the call's time limit (timeout_s) had passed. A call of its own (name) checks it"),
        Evidence::Unsolved(e) => {
            r["unavailable"] = json!(e.message);
            if let Some(h) = &e.hint {
                r["hint"] = json!(h);
            }
        }
        Evidence::Overlaps(p, thin) => {
            proof_fields(&mut r, p, names, None);
            if let Some(side) = thin {
                r["hint"] = json!(format!(
                    "what both shapes hold lies in a proven box {} thick, under twice the accuracy ({}): no point of it is deeper than the accuracy, so the overlap cannot be shown at it. \
                     Ask an (accuracy x) under {}, or test a point with (probe x y z), decided exactly",
                    measured(*side),
                    measured(p.tol),
                    measured(side / 2.0)
                ));
            }
        }
        Evidence::Touches { overlap, gap } => {
            match overlap.as_ref() {
                Evidence::Empty(p) => proof_fields(&mut r, p, names, None),
                Evidence::Vacuous(note) => r["note"] = json!(note),
                _ => {}
            }
            if let Evidence::ClearanceSearch { pair: Some((q, p)), verdict: Verdict::Disproven { value, .. }, .. } = gap.as_ref() {
                r["contact"] = json!({"in_second": v3(*q), "in_first": v3(*p), "distance": n(*value)});
            }
        }
        Evidence::Point { at, bounds } => {
            r["point"] = v3(*at);
            r["value"] = crate::command::bounds(bounds.0, bounds.1);
        }
    }
    // The answer's accuracy: null, unknown or an estimate.
    accuracy_fields(&mut r, c.status, c.asked);
    // How a guaranteed emptiness held: outright, or only to the accuracy
    // (with the deepest violating point its loose cells' centres showed).
    match c.held.filter(|_| c.holds_guaranteed()) {
        Some(crate::claims::Held::Outright) => r["held"] = json!("outright"),
        Some(crate::claims::Held::ToAccuracy { near }) => {
            r["held"] = json!("to the accuracy");
            if let Some((depth, at)) = near {
                r["near"] = json!({"depth": depth});
                if let Some(p) = at {
                    r["near"]["at"] = v3(p);
                }
            }
        }
        None => {}
    }
    if let Some(m) = &c.motion {
        motion_fields(&mut r, m);
        left_fields(&mut r, m, all);
    }
    if let Some(w) = c.in_window {
        r["in_window"] = json!(w.as_str());
    }
    if c.unbounded {
        r["unbounded"] = json!(true);
    }
    let spans = spans_json(names, &c.ranged);
    if spans.as_object().is_some_and(|o| !o.is_empty()) {
        r["param_spans"] = spans;
    }
    if let Some(note) = c.note {
        r["note"] = json!(note);
    }
    // Asked at that moment only, though it reads t (the answer's note
    // says so once).
    if let Some(t) = c.one_moment {
        r["one_moment"] = n(t);
    }
    inexact_fields(&mut r, &c.inexact, all);
    r
}

/// A claim measured at moments where t ranged did not decide it: why
/// (`motion.not_ranged`), and where the search ran and was left undecided,
/// its own answer (`motion.t_ranged`: the hint, undecided cell and the
/// params it ranged, t among them), and a `hint` saying both: the estimate
/// says what would make it a proof, as an unknown does.
fn left_fields(r: &mut J, s: &crate::motion::Stepped, all: &[String]) {
    let (Some(why), false) = (&s.left, s.bound) else { return };
    let mut hint = format!("an estimate: checked only at the moments in motion.steps, since {why}");
    if let Some(a) = &s.left_answer {
        let inner = check_json(a, all);
        let mut o = serde_json::Map::new();
        for k in ["hint", "note", "undecided_cell", "undecided_cells", "accuracy_attained", "param_spans", "res"] {
            if let Some(x) = inner.get(k).filter(|x| !x.is_null()) {
                o.insert(k.into(), x.clone());
            }
        }
        if a.note == Some(crate::claims::AT_ONE_VALUE) {
            // A contact (or overlap) found at one moment: what would decide
            // it is a claim that can hold at every moment.
            let others = o.get("param_spans").and_then(J::as_object).is_some_and(|p| p.keys().any(|k| k != "t"));
            let at = if others { "one moment and value of the ranged params" } else { "one moment" };
            hint = format!(
                "measured only at the moments in motion.steps. A contact or overlap found at {at} says nothing of the others, so ranging t cannot prove it. \
                 To prove it, claim what can hold at every moment (a clearance, disjoint, a value), or ask it at rest with why the motion can't change it"
            );
        } else if let Some(h) = inner.get("hint").or(inner.get("note")).and_then(J::as_str) {
            hint += &format!(" (motion.t_ranged): {h}");
        }
        r["motion"]["t_ranged"] = J::Object(o);
    }
    if r.get("hint").is_none_or(J::is_null) {
        r["hint"] = json!(hint);
    }
}

/// How a claim (or pair) over a motion was checked (`how`): with t ranged
/// (guaranteed), asked once where it rides one moving instance, or measured
/// at moments (an estimate, a failure aside, with why t ranged did not
/// decide it); a disproof's moment, and what was left undecided.
pub fn motion_fields(r: &mut J, s: &crate::motion::Stepped) {
    // The deciding step's own hint, from its evidence.
    let own = r.get("hint").and_then(J::as_str).map(str::to_string);
    let m = &s.motion;
    let spans = |v: &[(f64, f64)]| J::Array(v.iter().take(8).map(|&(a, b)| bounds(a, b)).collect());
    let mut o = json!({"span": bounds(m.span.0, m.span.1)});
    if s.riding {
        o["how"] = json!("asked once");
    } else if s.ranged {
        o["how"] = json!("t ranged");
    } else {
        o["how"] = json!("measured");
        o["steps"] = json!(m.steps);
        match &s.basis {
            // Checked at the moments a provider measures what it reads at.
            Some(basis) => o["basis"] = json!(basis),
            None => o["step"] = n(s.step.unwrap_or((m.span.1 - m.span.0) / m.steps.max(1) as f64)),
        }
    }
    if let Some(why) = &s.left {
        o["not_ranged"] = json!(why);
    }
    if !m.undecided.is_empty() {
        o["undecided"] = spans(&m.undecided);
    }
    // The moment of a failure.
    if let Some(t) = m.t {
        r["t"] = n(t);
    }
    if !m.loose.is_empty() {
        o["loose"] = spans(&m.loose);
        r["hint"] = json!(match s.motion.steps >= crate::motion::MAX_STEPS || s.basis.is_some() {
            false => "the pieces were halved until their parts moved within the precision asked. A smaller (accuracy x) halves them further, a larger one accepts a deeper overlap",
            true => "the budget ran out before the pieces were decided. Check a shorter --t-span (one part of the motion) at a time",
        });
    }
    if !m.unchecked.is_empty() {
        o["unchecked"] = spans(&m.unchecked);
    }
    // Unknown where the limit left steps or moments undecided: how many,
    // and the rest's own hint (a dead point's).
    if m.stopped > 0 && m.status == ink_kernel::ranged::Status::Unknown && !m.undecided.is_empty() {
        let shown = if m.undecided.len() > 8 { " (the first 8 listed)" } else { "" };
        let mut h = match m.stopped < m.undecided.len() {
            false => format!("the call's time limit (timeout_s) stopped it: those in motion.undecided{shown} were not decided"),
            true => format!(
                "the call's time limit (timeout_s) stopped {} of the {} in motion.undecided{shown} before they were decided. The others are undecided for their own reason",
                m.stopped,
                m.undecided.len()
            ),
        };
        h += ". A larger timeout_s checks those it stopped again, and a shorter --t-span costs less";
        if let Some(own) = own.filter(|_| m.stopped < m.undecided.len()) {
            h += &format!(". The rest: {own}");
        }
        r["hint"] = json!(h);
    }
    r["motion"] = o;
}

#[cfg(test)]
mod tests {
    /// A bound's prose is never below it, at decimal boundaries most of all.
    #[test]
    fn a_bound_rounds_up_in_prose() {
        for base in [0.011, 0.012, 0.1, 1.0, 2.5, 47.0, 0.023, 9.9, 100.0] {
            let mut x = base;
            for _ in 0..50 {
                assert!(super::up(x).parse::<f64>().unwrap() >= x, "{x} printed as {}", super::up(x));
                x = x.next_up();
            }
        }
    }
}
