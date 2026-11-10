//! Units. The language reads each quantity in one unit (lengths mm,
//! pressures MPa, ...: `Quantity::unit`), so a bare number has one meaning
//! everywhere. A number may say its unit instead, `1in`, `50psi`, `68degF`,
//! read as the language's (`UNITS`); a form that takes a known quantity
//! refuses a value said in another's (`Said`). Answers are in the
//! language's units.

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

/// Value `x`, in the language's unit, written in unit `u` so that it
/// reads back exactly (a --param an answer gives back): in `u` where its
/// number does, else in the language's unit of its quantity.
pub fn write(x: f64, u: &'static Unit) -> String {
    let w = format!("{}{}", to(x, u), u.name);
    match crate::sexp::number(&w) {
        Ok((y, _)) if y == x => w,
        _ => format!("{x}{}", u.quantity.unit()),
    }
}

/// Of two units values are said in that are combined (a sum's terms, a
/// range's ends), the one the result is said in: one with a zero of its
/// own (degF, K) over one without, as the result reads a point on its
/// scale. Err: two quantities, "are one quantity, but ...".
pub fn kept(a: &'static Unit, b: &'static Unit) -> Result<&'static Unit, String> {
    if a.quantity != b.quantity {
        return Err(format!("are one quantity, but {} is {} and {} is {}", a.name, a.quantity.name(), b.name, b.quantity.name()));
    }
    let offset = |u: &Unit| matches!(u.scale, Scale::Offset { .. });
    Ok(if offset(b) && !offset(a) { b } else { a })
}

/// A value given a param from outside the model (`--param m=50g`, an
/// optimizer's design): a number in the language's unit, and the unit it
/// was given in, if any. A param given a unit reads as said in it
/// (Doc::param_given), wherever the value travels. `text`: the value as
/// written, which a text param (`(param name "Ada")`) takes; one that is
/// no number has a NaN `value`, refused by a number param.
#[derive(Clone, Debug)]
pub struct Given {
    pub value: f64,
    pub unit: Option<&'static Unit>,
    pub text: Option<std::sync::Arc<str>>,
    /// A param's declared case by its place among them (from 0), as `ink
    /// assert` asks it: the load takes that case's value where it reads
    /// the param (its expression there, every param set before it), and
    /// at a load given one, a text a param does not allow is recorded
    /// (Doc::refused), not an error.
    pub case: Option<usize>,
}

impl PartialEq for Given {
    fn eq(&self, o: &Given) -> bool {
        (self.value == o.value || self.value.is_nan() && o.value.is_nan()) && self.unit.map(|u| u.name) == o.unit.map(|u| u.name) && self.text == o.text && self.case == o.case
    }
}

impl From<f64> for Given {
    /// A number in the language's unit, saying none.
    fn from(value: f64) -> Given {
        Given { value, unit: None, text: None, case: None }
    }
}

impl Given {
    /// Text only: what a text param takes, no number.
    pub fn text(s: &str) -> Given {
        Given { value: f64::NAN, unit: None, text: Some(s.into()), case: None }
    }

    /// A param's k-th declared case (`case`), its value the load's.
    pub fn case(k: usize) -> Given {
        Given { value: f64::NAN, unit: None, text: None, case: Some(k) }
    }
}

/// A value in the language's unit of its quantity, in unit `u` (`from`'s
/// inverse).
pub fn to(x: f64, u: &Unit) -> f64 {
    match u.scale {
        Scale::Pow10(p) => x / 10f64.powi(p),
        Scale::Factor(k) => x / k,
        Scale::Offset { zero, factor } => x / factor + zero,
    }
}

/// The number written with a unit that a value's quantity was said by,
/// carried with the value as it is computed (Value::Field): through names,
/// lets, calls and the branch taken, and through sums (`+ - neg abs min
/// max`, select's branches), which take one quantity; a product's units
/// combine, and it carries none. A slot taking a known quantity checks it.
#[derive(Clone, Debug)]
pub struct Said {
    pub unit: &'static Unit,
    /// The number as written, in `unit`.
    pub written: f64,
    pub file: std::sync::Arc<str>,
    pub line: u32,
}

impl Said {
    /// "50psi (part.ink:4)"; given on the command line, "50psi (--param p)".
    pub fn cite(&self) -> String {
        match self.line {
            0 => format!("{}{} ({})", self.written, self.unit.name, self.file),
            n => format!("{}{} ({}:{n})", self.written, self.unit.name, self.file),
        }
    }

    /// Why a value said by this, where a form takes `q` (`what`: the form,
    /// for the message), cannot be one. None: it is one.
    pub fn misfit(&self, q: Quantity, what: &str) -> Option<String> {
        (self.unit.quantity != q).then(|| format!("{}, written {}", wrong(self.unit, q, what), self.cite()))
    }
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
    fn a_value_is_written_back_in_its_unit() {
        // As a --param an answer gives back: it reads back exactly.
        let back = |x: f64, u: &str| {
            let w = write(x, unit(u).unwrap());
            (crate::sexp::number(&w).unwrap().0 == x).then_some(w)
        };
        assert_eq!(back(read("50", unit("psi").unwrap()).unwrap(), "psi").as_deref(), Some("50psi"));
        assert_eq!(back(20.0, "degF").as_deref(), Some("68degF"));
        // Not exact in its unit: in the language's, of its quantity.
        assert_eq!(back(0.1 + 0.2, "in").map(|w| w.ends_with("mm") || w.ends_with("in")), Some(true));
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
