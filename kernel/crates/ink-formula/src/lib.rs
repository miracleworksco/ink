//! The formula tool: a fact read off a maker's table at another fact. It
//! plugs into the claims core and asks it for the fact the table is read
//! at; it knows nothing else of ink. A unit is a word it only compares:
//! whoever builds a table gives it in the units its facts are in.

use ink_claims::{About, Claims, Fact, Strength, Tool, Unavailable, Verdict, resting, same};

/// A curve as its maker gave it: points, x rising, in these two units.
///
/// At a point it says that point's y. Between two points it says the y on
/// the straight line between them, and stands only behind this: y is
/// somewhere between the two points' ys. That is an assumption, so the
/// answer is an estimate and lists the gap under what it left out. Before
/// the first point and after the last it says nothing.
///
/// `line`: the maker gives the straight lines themselves (a DC motor's
/// torque from stall to free speed), so between two points it says the
/// line's y and stands behind it.
///
/// `below`: what it says before its first point, where whoever built it
/// says the curve is no lower there than at that point (a stepper's
/// pull-out torque toward rest, under its holding torque): the first
/// point's y, an estimate, saying so. `Err`: it says nothing there, and
/// why (the reason its builder gives).
pub struct Table {
    pub x: String,
    pub y: String,
    pub points: Vec<(f64, f64)>,
    pub line: bool,
    /// The points are read off the maker's drawn curve, each to within this much of its y: between two of
    /// them the curve is their line, this much either side, an estimate. None: each point is a measurement
    /// alone (a datasheet's table), and between two the value is anywhere between theirs.
    pub drawn: Option<f64>,
    pub below: Result<(), String>,
}

/// The fact `about` is `table`'s y where its x is the fact `at`, or with
/// `size` that fact's size (its number without its sign).
pub struct Formula {
    pub about: About,
    pub table: Table,
    pub at: About,
    pub size: bool,
    pub from: String,
}

pub struct Formulas(pub Vec<Formula>);

impl Tool for Formulas {
    fn ask(&self, about: &About, claims: &mut Claims) -> Option<Fact> {
        let formula = self.0.iter().find(|f| &f.about == about)?;
        let x = claims.answer(&formula.at);
        let mut fact = resting(about, &formula.from, [&x].into_iter());
        let points = &formula.table.points;
        let outside = Verdict::NotAvailable(Unavailable::Outside);
        let (x_unit, first) = (&formula.table.x, points.first().map_or(f64::NAN, |p| p.0));
        // Why it has no answer, in words, beside the verdict: an author reads which it was.
        let mut why: Option<String> = None;
        // Whether it was read before its first point, as that point's y.
        let mut before = false;
        let mut read = || {
            if let Some(why) = x.held.filter(|v| v.open()) {
                return Err(why); // what stopped it passes up, as its own reason
            }
            // Another unit: said as the unit it is, never read as if it fitted.
            if x.value.is_some() && x.unit != *x_unit {
                let stated = if x.unit.is_empty() { "no unit".to_string() } else { x.unit.clone() };
                why = Some(format!("its table reads '{}' in {x_unit}, and it is stated in {stated}", formula.at.what));
                return Err(Verdict::NotAvailable(Unavailable::Units));
            }
            // A yes/no, or a table out of order.
            let (low, high) = x.range().filter(|_| points.windows(2).all(|w| w[0].0 < w[1].0)).ok_or(outside)?;
            let (low, high) = match formula.size {
                true => (if low <= 0.0 && high >= 0.0 { 0.0 } else { low.abs().min(high.abs()) }, low.abs().max(high.abs())),
                false => (low, high),
            };
            // The last point at or before the low end, and the first at or after the high end (`same`: rounding is no distance).
            let from = match (points.iter().rposition(|p| p.0 <= low || same(p.0, low)), &formula.table.below) {
                (Some(from), _) => from,
                (None, Ok(())) if !points.is_empty() => {
                    before = true;
                    0
                }
                (None, below) => {
                    let said = below.as_ref().err().map(|why| format!(": {why}")).unwrap_or_default();
                    why = Some(format!("its table begins at {first} {x_unit}, and '{}' is as low as {low} {x_unit}{said}", formula.at.what));
                    return Err(outside);
                }
            };
            let Some(to) = points.iter().position(|p| p.0 >= high || same(p.0, high)) else {
                let last = points.last().map_or(f64::NAN, |p| p.0);
                why = Some(format!("its table ends at {last} {x_unit}, and '{}' is as high as {high} {x_unit}", formula.at.what));
                return Err(outside);
            };
            // Points nearer each other than rounding tells apart can each be the one: all of them are taken.
            Ok((low, high, from.min(to), from.max(to)))
        };
        // The line between the two points about x (the end's value past an end), each point weighed, so no difference of two ys overflows.
        let line = |x: f64| {
            let i = points.iter().rposition(|p| p.0 <= x).unwrap_or(0);
            let ((x0, y0), (x1, y1)) = (points[i], points[(i + 1).min(points.len() - 1)]);
            let t = if x1 > x0 { ((x - x0) / (x1 - x0)).clamp(0.0, 1.0) } else { 0.0 };
            y0 * (1.0 - t) + y1 * t
        };
        let found = read();
        match found {
            Ok((low, high, from, to)) => {
                // What it stands behind: on a line, the line over the range; else the points about it.
                let drawn = formula.table.drawn;
                let ys: Vec<f64> = match formula.table.line || drawn.is_some() {
                    true => [line(low), line(high)].into_iter().chain(points.get(from + 1..to).unwrap_or(&[]).iter().map(|p| p.1)).collect(),
                    false => points[from..=to].iter().map(|p| p.1).collect(),
                };
                let (least, most) = (ys.iter().copied().fold(f64::INFINITY, f64::min), ys.iter().copied().fold(f64::NEG_INFINITY, f64::max));
                // A drawn curve is read to its reading error, either side.
                let (least, most) = (least - drawn.unwrap_or(0.0), most + drawn.unwrap_or(0.0));
                fact.unit = formula.table.y.clone();
                fact.value = x.value.map(|x| line(if formula.size { x.abs() } else { x }));
                fact.band = (least < most).then_some((least, most));
                if ![least, most].into_iter().chain(fact.value).all(f64::is_finite) {
                    // A table whose numbers are none says nothing.
                    (fact.value, fact.band, fact.unit, fact.held, fact.strength) = (None, None, String::new(), Some(outside), None);
                    return Some(fact);
                }
                if let Some(read) = drawn {
                    fact.strength = Some(Strength::Estimated);
                    // To six figures: a unit's conversion leaves rounding in it.
                    let shown = if read > 0.0 { format!("{:.5e}", read).parse::<f64>().unwrap_or(read) } else { read };
                    fact.notes.push(format!("the maker's drawn curve, read between its points as their line, to within {shown} {} either side", formula.table.y));
                } else if !formula.table.line && from < to {
                    fact.strength = Some(Strength::Estimated);
                    fact.left_out.push(format!("the table between {} and {} {}: taken as between its values there", points[from].0, points[to].0, formula.table.x));
                }
                if before {
                    fact.strength = Some(Strength::Estimated);
                    fact.left_out.push(format!("the table below its first point, {first} {x_unit}: taken as its value there"));
                }
            }
            Err(open) => {
                (fact.held, fact.strength) = (Some(open), None);
                fact.notes.extend(why);
            }
        }
        Some(fact)
    }
}
