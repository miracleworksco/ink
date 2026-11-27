//! Units. The language reads each quantity in one unit (lengths mm,
//! pressures MPa, ...: `Quantity::unit`), so a bare number has one meaning
//! everywhere. A number may say its unit instead, `1in`, `50psi`, `68degF`,
//! read as the language's (`UNITS`); a form that takes a known quantity
//! refuses a value said in another's (`Said`). Answers are in the
//! language's units.

/// A kind of quantity, and the one unit the language reads it in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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
    /// An angle, in degrees: what every form that takes one reads a bare
    /// number as.
    Angle,
    AngularSpeed,
    /// How fast a turning rate changes: what a turning profile ramps at.
    AngularAcceleration,
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
    /// What a value is that is none of the above: a product or quotient
    /// of such a dimension (a force times a time, seconds per degree; its
    /// powers are kept, `Said::powers`), or a value whose unit an
    /// operation could not follow (the sine of a length). It has no unit
    /// to be read in, is taken where a bare number is, and is never a pure
    /// number again (`OTHER`).
    Other,
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
            Angle => "an angle",
            AngularSpeed => "a turning rate",
            AngularAcceleration => "a turning acceleration",
            Acceleration => "an acceleration",
            FlowRate => "a flow rate",
            Viscosity => "a viscosity",
            Conductivity => "a thermal conductivity",
            Capacity => "a specific heat",
            FilmCoefficient => "a convection coefficient",
            HeatFlux => "a heat flux",
            Stiffness => "a stiffness",
            Inertia => "a moment of inertia",
            Other => "a quantity the language has no unit for",
        }
    }

    /// `e`, a bare number worked out, said in this quantity's unit, as a
    /// hint writes it: `(* 1N e)`. A temperature is a point on its scale,
    /// which is never multiplied: its scale's zero is added, `(+ 0degC e)`.
    pub fn say(self, e: &str) -> String {
        match self {
            Temperature => format!("(+ 0degC {e})"),
            q => format!("(* 1{} {e})", q.unit()),
        }
    }

    /// `x`, a value of this quantity, as a bare number in its unit, as a
    /// hint writes it: `(/ x 1N)`. A temperature is first made the
    /// difference from its scale's zero, an amount that divides:
    /// `(/ (- x 0degC) (- 1degC 0degC))`.
    pub fn unsay(self, x: &str) -> String {
        match self {
            Temperature => format!("(/ (- {x} 0degC) (- 1degC 0degC))"),
            q => format!("(/ {x} 1{})", q.unit()),
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
            Angle => "deg",
            AngularSpeed => "deg/s",
            AngularAcceleration => "deg/s2",
            Acceleration => "m/s2",
            FlowRate => "mm3/s",
            Viscosity => "Pa*s",
            Conductivity => "W/m/K",
            Capacity => "J/kg/K",
            FilmCoefficient => "W/m2/K",
            HeatFlux => "W/mm2",
            Stiffness => "N/mm",
            Inertia => "g*mm2",
            Other => "?",
        }
    }
}

/// A quantity's dimension: its powers of kilogram, metre, second, kelvin
/// and degree of angle, and how many of the SI unit of those powers (with
/// degrees for the angle) one of the language's unit of it is. What a
/// product, a quotient, a power or a root of quantities is (`quantity_of`).
/// An angle is a dimension of its own, in degrees, as every form that
/// takes one reads it; where a product needs it in radians is one rule
/// (lib.rs product_said).
pub fn dimension(q: Quantity) -> ([i32; 5], f64) {
    match q {
        Length => ([0, 1, 0, 0, 0], 1e-3),
        Area => ([0, 2, 0, 0, 0], 1e-6),
        Volume => ([0, 3, 0, 0, 0], 1e-9),
        Mass => ([1, 0, 0, 0, 0], 1e-3),
        Density => ([1, -3, 0, 0, 0], 1e6),
        Force => ([1, 1, -2, 0, 0], 1.0),
        Pressure => ([1, -1, -2, 0, 0], 1e6),
        Torque => ([1, 2, -2, 0, 0], 1e-3),
        Temperature => ([0, 0, 0, 1, 0], 1.0),
        Time => ([0, 0, 1, 0, 0], 1.0),
        Frequency => ([0, 0, -1, 0, 0], 1.0),
        Power => ([1, 2, -3, 0, 0], 1.0),
        Speed => ([0, 1, -1, 0, 0], 1e-3),
        Angle => ([0, 0, 0, 0, 1], 1.0),
        AngularSpeed => ([0, 0, -1, 0, 1], 1.0),
        AngularAcceleration => ([0, 0, -2, 0, 1], 1.0),
        Acceleration => ([0, 1, -2, 0, 0], 1.0),
        FlowRate => ([0, 3, -1, 0, 0], 1e-9),
        Viscosity => ([1, -1, -1, 0, 0], 1.0),
        Conductivity => ([1, 1, -3, -1, 0], 1.0),
        Capacity => ([0, 2, -2, -1, 0], 1.0),
        FilmCoefficient => ([1, 0, -3, -1, 0], 1.0),
        HeatFlux => ([1, 0, -3, 0, 0], 1e6),
        Stiffness => ([1, 0, -2, 0, 0], 1e3),
        Inertia => ([1, 2, 0, 0, 0], 1e-9),
        // Its own powers are its value's (Said::dims).
        Other => ([0; 5], 1.0),
    }
}

/// The quantity of this dimension, if the language has one.
pub fn quantity_of(dim: [i32; 5]) -> Option<Quantity> {
    const ALL: [Quantity; 25] = [
        Length,
        Area,
        Volume,
        Mass,
        Density,
        Force,
        Pressure,
        Torque,
        Temperature,
        Time,
        Frequency,
        Power,
        Speed,
        Angle,
        AngularSpeed,
        AngularAcceleration,
        Acceleration,
        FlowRate,
        Viscosity,
        Conductivity,
        Capacity,
        FilmCoefficient,
        HeatFlux,
        Stiffness,
        Inertia,
    ];
    ALL.into_iter().find(|q| dimension(*q).0 == dim)
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

impl Unit {
    /// Whether a step of it is a kelvin (K, degC; not degF): what a bare
    /// number added to a zero of it (`(+ 0degC n)`) is read in.
    pub fn steps_as_kelvins(&self) -> bool {
        matches!(self.scale, Scale::Pow10(0) | Scale::Factor(1.0) | Scale::Offset { factor: 1.0, .. })
    }
}

/// The unit of `Quantity::Other`: none. Not in `UNITS`, so never written.
pub static OTHER: Unit = Unit { name: "?", quantity: Quantity::Other, scale: Scale::Factor(1.0) };

/// The unit of a temperature difference (a temperature less a temperature,
/// in any scale): kelvins, an amount, where a temperature said in K, degC
/// or degF is a point on a scale. A difference scales and multiplies; a
/// point does not. Not in `UNITS`: a difference has no literal.
pub static DIFFERENCE: Unit = Unit { name: "K", quantity: Quantity::Temperature, scale: Scale::Factor(1.0) };

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
    p("deg", Angle, 0),
    f("rad", Angle, 180.0 / std::f64::consts::PI),
    p("deg/s", AngularSpeed, 0),
    f("rpm", AngularSpeed, 6.0),
    f("rad/s", AngularSpeed, 180.0 / std::f64::consts::PI),
    p("deg/s2", AngularAcceleration, 0),
    f("rad/s2", AngularAcceleration, 180.0 / std::f64::consts::PI),
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

/// A number written in unit `u` (its digits' exact value `r`), exactly in the language's unit: where the
/// unit is a power of ten, or a factor written as a short decimal (25.4 to the inch, 304.8 to the foot,
/// 453.59237 to the pound); a factor of pi (a degree) or an offset scale has none.
pub fn exact(r: ink_core::exact::Ratio, u: &Unit) -> Option<ink_core::exact::Ratio> {
    use ink_core::exact::Ratio;
    match u.scale {
        Scale::Pow10(p) => r.times(if p >= 0 { Ratio::new(10i128.checked_pow(p as u32)?, 1)? } else { Ratio::new(1, 10i128.checked_pow(p.unsigned_abs())?)? }),
        Scale::Factor(k) => {
            let written = format!("{k}");
            let digits = written.trim_start_matches(['-', '0', '.']).replace('.', "");
            (digits.len() <= 12).then(|| Ratio::of_decimal(&written)).flatten().and_then(|f| r.times(f))
        }
        Scale::Offset { .. } => None,
    }
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

/// What a bare number is one of: a length in mm, so an area and a volume, and, where a pose is taken, an angle in
/// degrees.
pub(crate) fn bare(q: Quantity) -> bool {
    matches!(q, Length | Area | Volume | Angle)
}

/// A declaration can mention units even when its result has no known
/// dimension (a product). That alone requires units on authored inputs;
/// it does not give the result a pressure exemption or guess its quantity.
#[derive(Clone)]
pub(crate) struct UnitRequirement {
    written: String,
    zero_ok: bool,
    /// Where the result's dimension is unknown, explicit vary bounds supply
    /// the unit for optimizer steps and their replay.
    pub(crate) bound_unit: Option<&'static Unit>,
}

impl UnitRequirement {
    pub(crate) fn of(expr: &crate::sexp::Sexp, value: f64, said: Option<&Said>) -> Option<Self> {
        use crate::sexp::{SK, Sexp, Written};
        fn visit(e: &Sexp, units: &mut Vec<&'static Unit>) {
            match &e.k {
                SK::Num(_, Some(Written::Here(u, _)), _) => units.push(u),
                SK::Num(_, Some(Written::From(s)), _) => units.push(s.unit),
                SK::List(xs) => xs.iter().for_each(|x| visit(x, units)),
                _ => {}
            }
        }
        let mut units = Vec::new();
        visit(expr, &mut units);
        units.extend(said.map(|s| s.unit));
        // A length, an area, a volume and an angle in the language's own unit take a bare number, which is that.
        // Every other quantity is given with its unit, as it is written with one everywhere.
        if units.is_empty() || (said.is_some_and(|s| bare(s.unit.quantity)) && units.iter().all(|u| u.name == u.quantity.unit())) {
            return None;
        }
        let written = if expr.list().is_some() { expr.to_string() } else { said.map_or_else(|| expr.to_string(), |s| s.display_value(value)) };
        Some(Self { written, bound_unit: None, zero_ok: units.iter().all(|u| !matches!(u.scale, Scale::Offset { .. })) })
    }

    pub(crate) fn carried(s: &Said, value: f64) -> Option<Self> {
        (s.unit.name != s.unit.quantity.unit() || !bare(s.unit.quantity)).then(|| Self {
            written: s.display_value(value),
            bound_unit: None,
            zero_ok: !matches!(s.unit.scale, Scale::Offset { .. }),
        })
    }

    pub(crate) fn check(&self, name: &str, input: &str, explicit: bool, zero: bool) -> Result<(), String> {
        if explicit || (zero && self.zero_ok) {
            return Ok(());
        }
        Err(format!("param {name} is written {}; give {input} a unit", self.written))
    }
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

/// What to write where a frequency (anything in 1/s: cycles a second, and equally a speed over a radius) meets
/// a turning rate. The two ways a 1/s arises differ by 2 pi, so the fix names both: no hint picks one.
pub const FREQUENCY_AS_RATE: &str = "A value in 1/s is read as a frequency, never as a turning rate: a cycle is not a radian. Say which it is: a speed over \
    a radius is radians a second, (* VALUE 1rad); a count of turns a second, (* VALUE 360deg)";
/// And a turning rate where a frequency is taken.
pub const RATE_AS_FREQUENCY: &str = "A turning rate is never read as a frequency: a cycle is not a radian. Its turns a second are (/ VALUE 360deg); its \
    radians a second, (/ VALUE 1rad)";

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
/// max`, select's branches), which take one quantity. Scaling by bare
/// numeric literals preserves it. A slot taking a known quantity checks it.
#[derive(Clone, Debug)]
pub struct Said {
    /// None when the value always states a unit. Otherwise a field,
    /// negative where a unit was explicit and zero where the input was
    /// bare. A select can have one branch of each, both of this quantity.
    pub(crate) explicit_when: Option<ink_core::graph::NodeId>,
    pub unit: &'static Unit,
    /// The number as written, in `unit`; NaN where it was worked out.
    pub written: f64,
    pub file: std::sync::Arc<str>,
    pub line: u32,
    /// For `Quantity::Other`: its dimension's powers and the SI value of
    /// one of its number (units::dimension), where they are known.
    pub powers: Option<([i32; 5], f64)>,
}

impl Said {
    /// A value worked out at `file:line`, said in `unit` (for `OTHER`, of `powers`).
    pub(crate) fn worked(unit: &'static Unit, powers: Option<([i32; 5], f64)>, explicit_when: Option<ink_core::graph::NodeId>, file: &str, line: u32) -> Said {
        Said { explicit_when, unit, written: f64::NAN, file: std::sync::Arc::from(file), line, powers }
    }

    /// Its dimension's powers and the SI value of one of its number; None where it is not known.
    pub fn dims(&self) -> Option<([i32; 5], f64)> {
        if self.unit.quantity == Quantity::Other { self.powers } else { Some(dimension(self.unit.quantity)) }
    }

    /// A temperature that is a point on a scale, not a difference.
    pub fn is_point(&self) -> bool {
        self.unit.quantity == Quantity::Temperature && !std::ptr::eq(self.unit, &DIFFERENCE)
    }

    /// Show the default in the author's unit, including when it was computed
    /// from several terms. A written decimal can be more exact than undoing
    /// the unit conversion in floating point.
    pub(crate) fn display_value(&self, value: f64) -> String {
        let original = self.written.to_string();
        let number = if read(&original, self.unit) == Some(value) { original } else { to(value, self.unit).to_string() };
        format!("{number}{}", self.unit.name)
    }

    /// "50psi (part.ink:4)"; given on the command line, "50psi (--param p)".
    pub fn cite(&self) -> String {
        // A quantity worked out (a product's), not written: no number to show.
        if self.written.is_nan() {
            return format!("{} worked out at {}:{}", self.unit.quantity.name(), self.file, self.line);
        }
        match self.line {
            0 => format!("{}{} ({})", self.written, self.unit.name, self.file),
            n => format!("{}{} ({}:{n})", self.written, self.unit.name, self.file),
        }
    }

    /// What the number of a value said by this is times, read as quantity
    /// `q` where a form takes one (`what`: the form, for the message); Err,
    /// why it cannot be one.
    ///
    /// Arithmetic carries an angle's power and never drops it, so this is
    /// where a radian appears: a value that is `q` times a power of angle
    /// (a radius times a turning rate, where a speed is taken) is read as
    /// `q` with the angle in radians. Not where `q` is itself an angle or a
    /// turning rate, and not as a frequency: a cycle is not a radian, so a
    /// frequency and a turning rate are never read as each other.
    pub fn read_as(&self, q: Quantity, what: &str) -> Result<f64, String> {
        use Quantity::{AngularSpeed, Frequency, Other};
        let own = self.unit.quantity;
        if own == q {
            return Ok(1.0);
        }
        if own != Other {
            let turn = match (own, q) {
                (Frequency, AngularSpeed) => FREQUENCY_AS_RATE,
                (AngularSpeed, Frequency) => RATE_AS_FREQUENCY,
                _ => "",
            };
            let turn = if turn.is_empty() { String::new() } else { format!(". {turn}") };
            return Err(format!("{}, written {}{turn}", wrong(self.unit, q, what), self.cite()));
        }
        // A value whose unit an operation could not follow (the sine of a length, an exponential of a time) is no
        // number of `q`'s unit: the unit that went in is still in it, unsaid.
        if self.powers.is_none() && q != Other {
            return Err(format!(
                "{what} takes {} ({}), but the value worked out at {}:{} lost its unit in a function that takes a pure number: divide what goes in by its unit first ((/ x 1mm), (/ t 1s))",
                q.name(),
                q.unit(),
                self.file,
                self.line
            ));
        }
        // What the language has no unit for is taken where a bare number is, but for an angle left in `q`.
        let ((of, one), Some((d, si))) = (dimension(q), self.powers) else { return Ok(1.0) };
        if d[..4] != of[..4] {
            // Where a bare number is taken (a length, a pure number), so is this. Elsewhere its dimension is known
            // and is not `q`'s: its number is in no unit of `q` (a product with a force made on the way put it in
            // g, mm and s), so reading it as one would be wrong by that scale with nothing said.
            if matches!(q, Quantity::Length | Quantity::Area | Quantity::Volume | Quantity::Angle | Other) {
                return Ok(1.0);
            }
            return Err(format!(
                "{what} takes {} ({}), but the value worked out at {}:{} is {}, not {}: its formula's units do not come out as one. Say each number's unit where it is written",
                q.name(),
                q.unit(),
                self.file,
                self.line,
                base_units(d),
                q.name()
            ));
        }
        if of[4] != 0 || q == Frequency {
            let turn = if q == Frequency { ": a cycle is not a radian, so say which is meant by dividing by 360deg (turns) or 1rad" } else { "" };
            return Err(format!(
                "{what} takes {} ({}), but the value worked out at {}:{} is {} with another power of angle{turn}",
                q.name(),
                q.unit(),
                self.file,
                self.line,
                q.name()
            ));
        }
        Ok(si * (std::f64::consts::PI / 180.0).powi(d[4]) / one)
    }

    /// The quantity a value said by this is read as where nothing says
    /// which (a fact, a printed unit), and what its number is times: its
    /// own; or, where it is a quantity times a power of angle, that
    /// quantity with the angle in radians (`read_as`). None: no quantity.
    pub fn read(&self) -> Option<(Quantity, f64)> {
        if self.unit.quantity != Quantity::Other {
            return Some((self.unit.quantity, 1.0));
        }
        let mut rest = self.powers?.0;
        rest[4] = 0;
        let q = quantity_of(rest).filter(|q| !matches!(q, Quantity::Frequency | Quantity::Temperature))?;
        Some((q, self.read_as(q, "").ok()?))
    }

    /// Why a value said by this, where a form takes `q` (`what`: the form,
    /// for the message), cannot be one as it is. None: it is one.
    pub fn misfit(&self, q: Quantity, what: &str) -> Option<String> {
        match self.read_as(q, what) {
            Err(why) => Some(why),
            Ok(times) if times != 1.0 => Some(format!(
                "{what} takes {}, but the value worked out at {}:{} is {} times an angle. Write (/ VALUE 1rad) to take the angle in radians",
                q.name(),
                self.file,
                self.line,
                q.name()
            )),
            Ok(_) => None,
        }
    }
}

/// A dimension in the units a value of no named quantity has its number in: "g mm2/s2/deg".
pub(crate) fn base_units(d: [i32; 5]) -> String {
    let term = |name: &str, p: i32| if p.abs() == 1 { name.to_string() } else { format!("{name}{}", p.abs()) };
    let names = ["g", "mm", "s", "K", "deg"];
    let over: Vec<String> = (0..5).filter(|&i| d[i] > 0).map(|i| term(names[i], d[i])).collect();
    let under: Vec<String> = (0..5).filter(|&i| d[i] < 0).map(|i| term(names[i], d[i])).collect();
    let top = if over.is_empty() { "1".to_string() } else { over.join(" ") };
    if under.is_empty() { top } else { format!("{top}/{}", under.join("/")) }
}

/// "pressure takes a pressure (MPa; or with its unit: ...), but lbf is a force".
pub fn wrong(u: &Unit, q: Quantity, what: &str) -> String {
    format!("{what} takes {} ({}, or a number with its unit: {}), but {} is {}", q.name(), q.unit(), units_of(q), u.name, u.quantity.name())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_dimension_names_one_quantity_and_the_units_agree_with_it() {
        // Every unit's quantity is found again from its dimension: no two share one.
        for u in UNITS {
            assert_eq!(quantity_of(dimension(u.quantity).0), Some(u.quantity), "{}", u.name);
        }
        // The SI factors, by units whose SI value is known: 1 N = 1 kg m/s2; 1 MPa mm2 = 1 N; 1 N m = 1 N*m.
        let si = |x: f64, name: &str| from(x, unit(name).unwrap()) * dimension(unit(name).unwrap().quantity).1;
        for (a, b) in [(si(1.0, "kg") * si(1.0, "m/s2"), si(1.0, "N")), (si(1.0, "MPa") * si(1.0, "mm2"), si(1.0, "N")), (si(1.0, "N") * si(1.0, "m"), si(1.0, "N*m"))] {
            assert!((a - b).abs() <= 1e-12 * b, "{a} {b}");
        }
        assert!((si(1.0, "in") - 0.0254).abs() < 1e-15 && (si(1.0, "lbf") - LBF).abs() < 1e-12 && si(1.0, "rpm") == 6.0 && si(2.0, "kW") == 2000.0);
    }

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
