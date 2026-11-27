//! The gearmotor 17HS19-1684S-PG27's pull-out curve, read at a speed.
//! data/ is copied whole from the catalog's record: `motor.curves[0]`,
//! torque at the gearbox output (N m) against its speed (rpm), ten points
//! read off the maker's chart.

use ink_claims::*;
use ink_formula::*;
use serde_json::Value;

const SHAFT: &str = "opener/motor/shaft";

fn pull_out() -> Table {
    let motor: Value = serde_json::from_str(include_str!("data/17HS19-1684S-PG27.json")).unwrap();
    let curve = &motor["motor"]["curves"][0];
    assert_eq!(
        (curve["kind"].as_str(), curve["x"].as_str(), curve["y"].as_str(), curve["at"]["shaft"].as_str()),
        (Some("pull-out"), Some("speed-rpm"), Some("torque-Nm"), Some("output"))
    );
    Table {
        below: Err(String::new()),
        line: false,
        drawn: None,
        x: "rpm".into(),
        y: "N*m".into(),
        points: curve["points"].as_array().unwrap().iter().map(|p| (p[0].as_f64().unwrap(), p[1].as_f64().unwrap())).collect(),
    }
}

fn speed(rpm: f64) -> Fact {
    Fact::stated(about("speed", SHAFT), rpm, "rpm", "stated here")
}

/// The torque the motor can give with `facts` stated, and whether `torque` is under it.
fn can(facts: Vec<Fact>) -> (Fact, Fact) {
    let formulas = Formulas(vec![Formula {
        about: about("pull-out-torque", SHAFT),
        table: pull_out(),
        at: about("speed", SHAFT),
        size: false,
        from: "17HS19-1684S-PG27 motor curves 0".into(),
    }]);
    let under = Rule { about: about("under", ""), says: Says::Ordered(vec![Term::Read(about("torque", SHAFT)), Term::Read(about("pull-out-torque", SHAFT))]), from: String::new() };
    let design = Design { facts, rules: vec![under] };
    let tools: &[&dyn Tool] = &[&formulas];
    let mut claims = Claims::new(&design, tools);
    (claims.answer(&about("pull-out-torque", SHAFT)), claims.answer(&about("under", "")))
}

#[test]
fn the_torque_against_the_real_curve() {
    // 22 rpm is between the curve's points at 19.5522 and 22.3453 rpm, 6.313 and 5.878 N m.
    let (can, under) = can(vec![speed(22.0), Fact::stated(about("torque", SHAFT), 3.4, "N*m", "stated here")]);
    assert!((can.value.unwrap() - 5.932).abs() < 0.001 && can.unit == "N*m", "{can:?}");
    assert_eq!(can.band, Some((5.878, 6.313)), "it stands only behind: between the two points' torques");
    assert_eq!(can.strength, Some(Strength::Estimated));
    assert_eq!(can.left_out, ["the table between 19.5522 and 22.3453 rpm: taken as between its values there"]);
    assert_eq!((&can.used[..], &can.from[..]), (&[about("speed", SHAFT)][..], "17HS19-1684S-PG27 motor curves 0"));
    // A torque under the lower of the two is under the curve, and the answer says how far to trust that.
    assert_eq!((under.held, under.strength, &under.left_out), (Some(Verdict::Yes), Some(Strength::Estimated), &can.left_out));
}

#[test]
fn what_a_table_says_at_between_and_past_its_points() {
    let can = |speed: Fact| can(vec![speed]).0;
    let outside = Some(Verdict::NotAvailable(Unavailable::Outside));

    // At a point: that point's torque, exactly, and nothing left out.
    let at = can(speed(22.3453));
    assert_eq!((at.value, at.band, at.strength, at.left_out.len()), (Some(5.878), None, None, 0), "as stated: the speed is");
    let (first, last) = (can(speed(2.7932)), can(speed(27.9317)));
    assert_eq!((first.value, first.band, last.value, last.band), (Some(6.964), None, Some(4.789), None));

    // Past either end it says nothing, however near; and what reads it is then not available.
    for rpm in [0.0, 2.79, 27.94, 30.0] {
        let past = can(speed(rpm));
        assert_eq!((past.held, past.value, past.strength), (outside, None, None), "{rpm}");
    }
    assert_eq!(self::can(vec![speed(30.0), Fact::stated(about("torque", SHAFT), 1.0, "N*m", "")]).1.held, outside);

    // A speed with a band takes every point from the one at or below its low end to the one at or above its high end:
    // 6 to 12 rpm spans the points at 5.5863 (7.182), 8.3795 (7.182), 11.1727 (6.857) and 13.9658 (6.747).
    let wide = can(Fact { band: Some((6.0, 12.0)), ..speed(9.0) });
    assert_eq!(wide.band, Some((6.747, 7.182)));
    assert_eq!(wide.left_out, ["the table between 5.5863 and 13.9658 rpm: taken as between its values there"]);
    assert_eq!(can(Fact { band: Some((20.0, 28.0)), ..speed(22.0) }).held, outside, "a band that runs past the last point");

    // A speed in another unit is refused, not converted: 2.3 rad/s is 22 rpm, and 22 of anything else is not.
    let units = Some(Verdict::NotAvailable(Unavailable::Units));
    let other = can(Fact { unit: "rad/s".into(), ..speed(2.3) });
    assert_eq!((other.held, other.notes.as_slice()), (units, &["its table reads 'speed' in rpm, and it is stated in rad/s".to_string()][..]));
    let bare = can(Fact { unit: "".into(), ..speed(22.0) });
    assert_eq!((bare.held, bare.notes.as_slice()), (units, &["its table reads 'speed' in rpm, and it is stated in no unit".to_string()][..]));
    // Before its first point (2.7932 rpm) it says nothing, and says where it begins; told its curve is no lower
    // there, it reads its first value, an estimate that says so.
    let slow = can(speed(1.0));
    assert!(slow.held == outside && slow.notes[0].starts_with("its table begins at 2.7932 rpm, and 'speed' is as low as 1 rpm"), "{slow:?}");
    // A yes/no and a missing speed pass on why.
    assert_eq!(can(Fact { value: None, held: Some(Verdict::Yes), ..speed(22.0) }).held, outside);
    assert_eq!(self::can(vec![]).0.held, Some(Verdict::NotAvailable(Unavailable::NoSource)));
    assert_eq!(can(Fact { held: Some(Verdict::NotSettled(Unsettled::Lacks)), ..speed(22.0) }).held, Some(Verdict::NotSettled(Unsettled::Lacks)));

    // A table out of order is refused.
    let mut table = pull_out();
    table.points.swap(3, 4);
    let formulas = Formulas(vec![Formula { about: about("t", ""), table, at: about("speed", SHAFT), size: false, from: String::new() }]);
    let design = Design { facts: vec![speed(22.0)], rules: vec![] };
    assert_eq!(Claims::new(&design, &[&formulas]).answer(&about("t", "")).held, outside);
    // It offers only its own tables.
    assert!(formulas.ask(&about("speed", SHAFT), &mut Claims::new(&design, &[])).is_none());

    // A table that is its maker's line (made up: 8 N*m at rest to none at 100 rpm) says the line's value and stands behind it.
    let dc = |speed: Fact| {
        let table = Table { below: Err(String::new()), line: true, drawn: None, x: "rpm".into(), y: "N*m".into(), points: vec![(0.0, 8.0), (100.0, 0.0)] };
        let formulas = Formulas(vec![Formula { about: about("t", ""), table, at: about("speed", SHAFT), size: false, from: String::new() }]);
        Claims::new(&Design { facts: vec![speed], rules: vec![] }, &[&formulas]).answer(&about("t", ""))
    };
    let on = dc(speed(25.0));
    assert_eq!((on.value, on.band, on.strength, on.left_out.len()), (Some(6.0), None, None, 0));
    assert_eq!(dc(Fact { band: Some((25.0, 50.0)), ..speed(30.0) }).band, Some((4.0, 6.0)));
    assert_eq!(dc(speed(101.0)).held, outside);

    // A line of three points (made up: flat at 8 N*m to 50 rpm, then down to none at 100): a range of speeds that
    // spans the corner takes the line at its two ends and the corner between them.
    let knee = |speed: Fact, size: bool| {
        let table = Table { below: Err(String::new()), line: true, drawn: None, x: "rpm".into(), y: "N*m".into(), points: vec![(0.0, 2.0), (50.0, 8.0), (100.0, 0.0)] };
        let formulas = Formulas(vec![Formula { about: about("t", ""), table, at: about("speed", SHAFT), size, from: String::new() }]);
        Claims::new(&Design { facts: vec![speed], rules: vec![] }, &[&formulas]).answer(&about("t", ""))
    };
    let across = knee(Fact { band: Some((25.0, 75.0)), ..speed(40.0) }, false);
    assert!(same(across.value.unwrap(), 6.8), "{across:?}");
    assert_eq!((across.band, across.left_out.len()), (Some((4.0, 8.0)), 0), "5 at 25, 8 at the corner, 4 at 75");
    assert_eq!(knee(Fact { band: Some((60.0, 75.0)), ..speed(70.0) }, false).band, Some((4.0, 6.4)), "no corner between");

    // Read at a size: a speed the other way round reads the same place; signed, it is before the table's first point.
    assert_eq!((knee(speed(-75.0), true).value, knee(speed(-75.0), false).held), (Some(4.0), outside));
    assert_eq!(knee(Fact { band: Some((-75.0, 25.0)), ..speed(-10.0) }, true).band, Some((2.0, 8.0)), "sizes from 0 to 75");

    // A speed at the table's end but for the rounding of its unit is at its end: 0.1 rpm in deg/s and back.
    let end = 0.1_f64 * 6.0 / 6.0 * 1000.0;
    assert!(end != 100.0 && same(end, 100.0));
    assert_eq!((knee(speed(end), false).value, knee(speed(100.0000001), false).held), (Some(0.0), outside));
    // And at its start: a table from 0.3 rpm, read at 0.1 + 0.2.
    let start = Table { below: Err(String::new()), line: true, drawn: None, x: "rpm".into(), y: "N*m".into(), points: vec![(0.3, 8.0), (100.0, 0.0)] };
    let formulas = Formulas(vec![Formula { about: about("t", ""), table: start, at: about("speed", SHAFT), size: false, from: String::new() }]);
    let at = |rpm: f64| Claims::new(&Design { facts: vec![speed(rpm)], rules: vec![] }, &[&formulas]).answer(&about("t", ""));
    assert!(at(0.3 - 5e-17).value.is_some_and(|y| same(y, 8.0)) && at(0.2999).held == outside, "{:?}", at(0.3 - 5e-17));
    // A table told its curve is no lower before its first point than at it (a stepper's pull-out torque toward
    // rest) reads that value there: an estimate, which says what it took. Across the point, the band of both.
    let toward_rest = Table { below: Ok(()), line: false, drawn: None, x: "rpm".into(), y: "N*m".into(), points: vec![(5.0, 2.0), (10.0, 3.0)] };
    let formulas = Formulas(vec![Formula { about: about("t", ""), table: toward_rest, at: about("speed", SHAFT), size: true, from: String::new() }]);
    let at = |fact: Fact| Claims::new(&Design { facts: vec![fact], rules: vec![] }, &[&formulas]).answer(&about("t", ""));
    let rest = at(speed(0.0));
    assert_eq!((rest.value, rest.band, rest.strength), (Some(2.0), None, Some(Strength::Estimated)), "{rest:?}");
    assert_eq!(rest.left_out, ["the table below its first point, 5 rpm: taken as its value there"]);
    let across = at(Fact { band: Some((0.0, 7.0)), ..speed(3.0) });
    assert_eq!((across.band, across.left_out.len()), (Some((2.0, 3.0)), 2), "{across:?}");
    // One told why it is not read there says so.
    let not_below = Table {
        below: Err("its first value is above its holding torque".into()),
        line: false,
        drawn: None,
        x: "rpm".into(),
        y: "N*m".into(),
        points: vec![(5.0, 2.0), (10.0, 3.0)],
    };
    let formulas = Formulas(vec![Formula { about: about("t", ""), table: not_below, at: about("speed", SHAFT), size: true, from: String::new() }]);
    let rest = Claims::new(&Design { facts: vec![speed(0.0)], rules: vec![] }, &[&formulas]).answer(&about("t", ""));
    assert_eq!(rest.notes, ["its table begins at 5 rpm, and 'speed' is as low as 0 rpm: its first value is above its holding torque"]);
}

/// A datasheet's torque-speed curve is a drawn line, and its points are read off it to a stated error: between
/// two points it is their line, that error either side, an estimate that says so. So it can say no between
/// chart points, where "anywhere between the two" could not. Points that are measurements each alone (a table's
/// rows) keep the wide reading.
#[test]
fn a_drawn_curve_is_its_points_line_to_its_reading_error() {
    let table = |drawn| Table { below: Ok(()), line: false, drawn, x: "rpm".into(), y: "N*m".into(), points: vec![(10.0, 2.0), (20.0, 1.0), (30.0, 0.8)] };
    let at = |drawn, fact: Fact| {
        let formulas = Formulas(vec![Formula { about: about("t", ""), table: table(drawn), at: about("speed", SHAFT), size: true, from: String::new() }]);
        Claims::new(&Design { facts: vec![fact], rules: vec![] }, &[&formulas]).answer(&about("t", ""))
    };
    // Half way between 2 and 1 N m: 1.5, to 0.01 either side. As lone measurements: anywhere from 1 to 2.
    let drawn = at(Some(0.01), speed(15.0));
    assert_eq!((drawn.value, drawn.strength), (Some(1.5), Some(Strength::Estimated)));
    assert!(drawn.band.is_some_and(|(low, high)| same(low, 1.49) && same(high, 1.51)), "{drawn:?}");
    assert_eq!(drawn.notes, ["the maker's drawn curve, read between its points as their line, to within 0.01 N*m either side"]);
    assert!(drawn.left_out.is_empty(), "{:?}", drawn.left_out);
    let alone = at(None, speed(15.0));
    assert_eq!((alone.value, alone.band), (Some(1.5), Some((1.0, 2.0))));
    // On a point: its value, to the same error. Over a range of speeds: the line's least and most over it.
    assert!(at(Some(0.01), speed(20.0)).band.is_some_and(|(low, high)| same(low, 0.99) && same(high, 1.01)));
    let over = at(Some(0.01), Fact { band: Some((15.0, 25.0)), ..speed(20.0) });
    assert!(over.band.is_some_and(|(low, high)| same(low, 0.89) && same(high, 1.51)), "{over:?}");
    // Below its first point, where it may be read: that point's value, to its error, and it says what it took.
    let rest = at(Some(0.01), speed(0.0));
    assert!(rest.value == Some(2.0) && rest.band.is_some_and(|(low, high)| same(low, 1.99) && same(high, 2.01)), "{rest:?}");
    assert_eq!(rest.left_out, ["the table below its first point, 10 rpm: taken as its value there"]);
}
