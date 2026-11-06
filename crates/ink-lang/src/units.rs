//! Units. The language reads each quantity in one unit (lengths mm,
//! pressures MPa, ...: `Quantity::unit`), so a bare number has one meaning
//! everywhere. A number may say its unit instead, `1in`, `50psi`, `68degF`,
//! read as the language's (`UNITS`); a form that takes a known quantity
//! refuses a number written in another's (`misfit`). Answers are in the
//! language's units.

use crate::sexp::{SK, Sexp};

/// A kind of quantity, and the one unit the language reads it in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quantity {
    Length,
    Area,
    Volume,
    Mass,
    Density,
    Force,
    Pressure,
    Torque,
    /// A temperature (a point on the scale: 68degF is 20 degC), not a
    /// change of one.
    Temperature,
    Time,
    Frequency,
    Power,
    Speed,
    AngularSpeed,
    Acceleration,
    FlowRate,
    Viscosity,
    Conductivity,
    Capacity,
    FilmCoefficient,
    HeatFlux,
    Stiffness,
    /// A moment of inertia (mass times area).
    Inertia,
}

use Quantity::*;

impl Quantity {
    /// What it is, for messages: "a length".
    pub fn name(self) -> &'static str {
        match self {
            Length => "a length",
            Area => "an area",
            Volume => "a volume",
            Mass => "a mass",
            Density => "a density",
            Force => "a force",
            Pressure => "a pressure or stress",
            Torque => "a torque",
            Temperature => "a temperature",
            Time => "a time",
            Frequency => "a frequency",
            Power => "a power",
            Speed => "a speed",
            AngularSpeed => "a turning rate",
            Acceleration => "an acceleration",
            FlowRate => "a flow rate",
            Viscosity => "a viscosity",
            Conductivity => "a thermal conductivity",
            Capacity => "a specific heat",
            FilmCoefficient => "a convection coefficient",
            HeatFlux => "a heat flux",
            Stiffness => "a stiffness",
            Inertia => "a moment of inertia",
        }
    }

    /// The unit the language reads and answers it in (metric): a name in
    /// `UNITS` with scale 1.
    pub fn unit(self) -> &'static str {
        match self {
            Length => "mm",
            Area => "mm2",
            Volume => "mm3",
            Mass => "g",
            Density => "g/mm3",
            Force => "N",
            Pressure => "MPa",
            Torque => "N*mm",
            Temperature => "degC",
            Time => "s",
            Frequency => "Hz",
            Power => "W",
            Speed => "mm/s",
            AngularSpeed => "deg/s",
            Acceleration => "m/s2",
            FlowRate => "mm3/s",
            Viscosity => "Pa*s",
            Conductivity => "W/m/K",
            Capacity => "J/kg/K",
            FilmCoefficient => "W/m2/K",
            HeatFlux => "W/mm2",
            Stiffness => "N/mm",
            Inertia => "g*mm2",
        }
    }
}

/// How a unit's number becomes the language's: a power of ten (which
/// moves the decimal point, so `2.7g/cm3` is the same double as `0.0027`),
/// a factor, or for a temperature scale `(x - zero) * factor`.
#[derive(Clone, Copy, Debug)]
pub enum Scale {
    Pow10(i32),
    Factor(f64),
    Offset { zero: f64, factor: f64 },
}

#[derive(Debug)]
pub struct Unit {
    pub name: &'static str,
    pub quantity: Quantity,
    pub scale: Scale,
}

const fn u(name: &'static str, quantity: Quantity, scale: Scale) -> Unit {
    Unit { name, quantity, scale }
}
const fn p(name: &'static str, quantity: Quantity, k: i32) -> Unit {
    u(name, quantity, Scale::Pow10(k))
}
const fn f(name: &'static str, quantity: Quantity, k: f64) -> Unit {
    u(name, quantity, Scale::Factor(k))
}

// The international inch (25.4 mm), foot, avoirdupois pound (453.59237 g),
// pound-force (the pound under standard gravity, 9.80665 m/s^2), and the
// International Table BTU (1055.05585262 J).
const IN: f64 = 25.4;
const FT: f64 = 304.8;
const LB: f64 = 453.59237;
const LBF: f64 = 4.4482216152605;
const BTU: f64 = 1055.05585262;
const F_PER_K: f64 = 1.8;

/// Every unit a number may be written in, read as the language's unit of
/// its quantity (the first of each quantity's here). Torque takes either
/// order (`N*m`, `in*lbf`); temperatures are temperatures, not changes.
pub const UNITS: &[Unit] = &[
    p("mm", Length, 0),
    p("cm", Length, 1),
    p("m", Length, 3),
    p("um", Length, -3),
    f("in", Length, IN),
    f("ft", Length, FT),
    f("thou", Length, IN / 1000.0),
    f("mil", Length, IN / 1000.0),
    p("mm2", Area, 0),
    p("cm2", Area, 2),
    p("m2", Area, 6),
    f("in2", Area, IN * IN),
    f("ft2", Area, FT * FT),
    p("mm3", Volume, 0),
    p("cm3", Volume, 3),
    p("m3", Volume, 9),
    p("mL", Volume, 3),
    p("L", Volume, 6),
    f("in3", Volume, IN * IN * IN),
    f("ft3", Volume, FT * FT * FT),
    p("g", Mass, 0),
    p("mg", Mass, -3),
    p("kg", Mass, 3),
    f("lb", Mass, LB),
    f("oz", Mass, LB / 16.0),
    p("g/mm3", Density, 0),
    p("g/cm3", Density, -3),
    p("kg/m3", Density, -6),
    f("lb/in3", Density, LB / (IN * IN * IN)),
    f("lb/ft3", Density, LB / (FT * FT * FT)),
    p("N", Force, 0),
    p("kN", Force, 3),
    f("lbf", Force, LBF),
    p("MPa", Pressure, 0),
    p("Pa", Pressure, -6),
    p("kPa", Pressure, -3),
    p("GPa", Pressure, 3),
    p("bar", Pressure, -1),
    f("psi", Pressure, LBF / (IN * IN)),
    f("ksi", Pressure, 1e3 * LBF / (IN * IN)),
    f("Msi", Pressure, 1e6 * LBF / (IN * IN)),
    p("N*mm", Torque, 0),
    p("N*m", Torque, 3),
    f("lbf*in", Torque, LBF * IN),
    f("lbf*ft", Torque, LBF * FT),
    f("in*lbf", Torque, LBF * IN),
    f("ft*lbf", Torque, LBF * FT),
    p("degC", Temperature, 0),
    u("degF", Temperature, Scale::Offset { zero: 32.0, factor: 1.0 / F_PER_K }),
    u("K", Temperature, Scale::Offset { zero: 273.15, factor: 1.0 }),
    p("s", Time, 0),
    p("ms", Time, -3),
    f("min", Time, 60.0),
    f("h", Time, 3600.0),
    p("Hz", Frequency, 0),
    p("kHz", Frequency, 3),
    p("W", Power, 0),
    p("mW", Power, -3),
    p("kW", Power, 3),
    f("hp", Power, 745.699_871_582_270_2),
    f("BTU/h", Power, BTU / 3600.0),
    p("mm/s", Speed, 0),
    p("m/s", Speed, 3),
    f("in/s", Speed, IN),
    f("ft/s", Speed, FT),
    f("ft/min", Speed, FT / 60.0),
    p("deg/s", AngularSpeed, 0),
    f("rpm", AngularSpeed, 6.0),
    f("rad/s", AngularSpeed, 180.0 / std::f64::consts::PI),
    p("m/s2", Acceleration, 0),
    p("mm/s2", Acceleration, -3),
    f("ft/s2", Acceleration, FT / 1000.0),
    f("in/s2", Acceleration, IN / 1000.0),
    p("mm3/s", FlowRate, 0),
    p("mL/s", FlowRate, 3),
    f("mL/min", FlowRate, 1e3 / 60.0),
    p("L/s", FlowRate, 6),
    f("L/min", FlowRate, 1e6 / 60.0),
    f("in3/s", FlowRate, IN * IN * IN),
    // US gallons (231 in^3) and cubic feet per minute.
    f("gpm", FlowRate, 231.0 * IN * IN * IN / 60.0),
    f("cfm", FlowRate, FT * FT * FT / 60.0),
    p("Pa*s", Viscosity, 0),
    p("mPa*s", Viscosity, -3),
    p("cP", Viscosity, -3),
    p("W/m/K", Conductivity, 0),
    f("BTU/h/ft/degF", Conductivity, BTU / 3600.0 / (FT / 1000.0) * F_PER_K),
    p("J/kg/K", Capacity, 0),
    f("BTU/lb/degF", Capacity, BTU / (LB / 1000.0) * F_PER_K),
    p("W/m2/K", FilmCoefficient, 0),
    f("BTU/h/ft2/degF", FilmCoefficient, BTU / 3600.0 / (FT / 1000.0 * FT / 1000.0) * F_PER_K),
    p("W/mm2", HeatFlux, 0),
    f("W/in2", HeatFlux, 1.0 / (IN * IN)),
    p("N/mm", Stiffness, 0),
    f("lbf/in", Stiffness, LBF / IN),
    p("g*mm2", Inertia, 0),
    p("kg*m2", Inertia, 9),
    f("lb*in2", Inertia, LB * IN * IN),
];

/// A unit by its name.
pub fn unit(name: &str) -> Option<&'static Unit> {
    UNITS.iter().find(|u| u.name == name)
}

/// The units a model may write, by quantity, for messages: "lengths mm cm
/// ...; ...".
pub fn listing() -> String {
    let mut out: Vec<(Quantity, Vec<&str>)> = Vec::new();
    for u in UNITS {
        match out.iter_mut().find(|(q, _)| *q == u.quantity) {
            Some((_, v)) => v.push(u.name),
            None => out.push((u.quantity, vec![u.name])),
        }
    }
    out.iter().map(|(q, v)| format!("{} {}", q.name().split_once(' ').unwrap().1, v.join(" "))).collect::<Vec<_>>().join("; ")
}

/// The units of `q` a model may write, for messages: "MPa Pa kPa ...".
pub fn units_of(q: Quantity) -> String {
    UNITS.iter().filter(|u| u.quantity == q).map(|u| u.name).collect::<Vec<_>>().join(" ")
}

/// A number written with a unit (`digits` its decimal, e.g. "2.7"), in the
/// language's unit. A power of ten joins the decimal's own exponent: one
/// correctly rounded parse.
pub fn read(digits: &str, u: &Unit) -> Option<f64> {
    Some(match u.scale {
        Scale::Pow10(p) => {
            let (m, e) = match digits.find(['e', 'E']) {
                Some(i) => (&digits[..i], digits[i + 1..].parse::<i32>().ok()?),
                None => (digits, 0),
            };
            format!("{m}e{}", e.saturating_add(p)).parse::<f64>().ok()?
        }
        _ => from(digits.parse::<f64>().ok()?, u),
    })
}

/// A value in unit `u`, in the language's unit of its quantity.
pub fn from(x: f64, u: &Unit) -> f64 {
    match u.scale {
        Scale::Pow10(p) => x * 10f64.powi(p),
        Scale::Factor(k) => x * k,
        Scale::Offset { zero, factor } => (x - zero) * factor,
    }
}

/// The unit expression `e` is written in, and where: of a number in it
/// written with one, looking only through sums (`+ - neg abs min max`,
/// select's branches; a product's units combine, and are not read), and
/// through names, by `named` (a def or param written in a unit).
pub fn written<'a>(e: &'a Sexp, named: &dyn Fn(&str) -> Option<&'static Unit>) -> Vec<(&'a Sexp, &'static Unit)> {
    match &e.k {
        SK::Num(_, Some((u, _))) => vec![(e, u)],
        SK::Sym(s) => named(s).map(|u| (e, u)).into_iter().collect(),
        SK::List(items) => {
            let Some((head, args)) = items.split_first() else { return Vec::new() };
            let args = match head.sym() {
                Some("+" | "-" | "neg" | "abs" | "min" | "max") => args,
                Some("select") if args.len() == 3 => &args[1..],
                _ => return Vec::new(),
            };
            args.iter().flat_map(|a| written(a, named)).collect()
        }
        _ => Vec::new(),
    }
}

/// Why expression `e`, where a form takes `q` (`what`: the form, for the
/// message), cannot be one, and where: a number in it written in another
/// quantity's unit (`written`). None: it fits, or says no unit.
pub fn misfit<'a>(e: &'a Sexp, q: Quantity, what: &str, named: &dyn Fn(&str) -> Option<&'static Unit>) -> Option<(&'a Sexp, String)> {
    let (at, u) = written(e, named).into_iter().find(|(_, u)| u.quantity != q)?;
    let mut m = wrong(u, q, what);
    if let Some(n) = at.sym() {
        m = format!("{m} ({n} is written in {})", u.name);
    }
    Some((at, m))
}

/// "pressure takes a pressure (MPa; or with its unit: ...), but lbf is a force".
pub fn wrong(u: &Unit, q: Quantity, what: &str) -> String {
    format!("{what} takes {} ({}, or a number with its unit: {}), but {} is {}", q.name(), q.unit(), units_of(q), u.name, u.quantity.name())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_quantity_has_its_units() {
        for u in UNITS {
            let q = u.quantity;
            // The language's unit is scale 1.
            assert_eq!(read("1", unit(q.unit()).unwrap()), Some(1.0), "{}", q.unit());
        }
        // One name, one unit.
        for (i, a) in UNITS.iter().enumerate() {
            assert!(UNITS[i + 1..].iter().all(|b| a.name != b.name), "{} twice", a.name);
        }
    }

    #[test]
    fn reads_against_references() {
        let r = |s: &str, q: Quantity| {
            let u = unit(s.trim_start_matches(|c: char| c.is_ascii_digit() || c == '.')).unwrap();
            assert_eq!(u.quantity, q);
            read(s.trim_end_matches(u.name), u).unwrap()
        };
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-6 * b.abs();
        // Exact by definition.
        assert_eq!(r("1in", Length), 25.4);
        assert_eq!(r("1ft", Length), 304.8);
        assert_eq!(r("1lb", Mass), 453.59237);
        assert_eq!(r("2.5m", Length), 2500.0);
        assert_eq!(r("68degF", Temperature), 20.0);
        assert_eq!(r("212degF", Temperature), 100.0);
        assert_eq!(r("273.15K", Temperature), 0.0);
        assert_eq!(r("3000rpm", AngularSpeed), 18000.0);
        // Against published conversions (NIST SP 811).
        assert!(close(r("1psi", Pressure), 6.894757e-3));
        assert!(close(r("1ksi", Pressure), 6.894757));
        assert!(close(r("1lbf", Force), 4.448222));
        assert!(close(r("1lbf*ft", Torque), 1355.818));
        assert!(close(r("1in*lbf", Torque), 112.9848));
        assert!(close(r("1oz", Mass), 28.34952));
        assert!(close(r("1hp", Power), 745.6999));
        assert!(close(r("1BTU/h", Power), 0.2930711));
        assert!(close(r("1BTU/h/ft/degF", Conductivity), 1.730735));
        assert!(close(r("1BTU/h/ft2/degF", FilmCoefficient), 5.678263));
        assert!(close(r("1BTU/lb/degF", Capacity), 4186.8));
        assert!(close(r("1gpm", FlowRate), 63090.20));
        assert!(close(r("1cfm", FlowRate), 471947.4));
        assert!(close(r("1ft/s2", Acceleration), 0.3048));
    }
}
