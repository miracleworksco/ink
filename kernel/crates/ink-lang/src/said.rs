//! Units as written: the unit a value is said in, what it must fit (a
//! quantity, a drive's change), and a param's given value.

use super::*;

impl Doc {
    /// The cell size of a `(res h)` clause (an assert's, a sampled field's,
    /// a study's): a number known at load time, above 0 and finite, as for
    /// `--res`.
    pub(super) fn cell_size(&mut self, e: &Sexp, env: &mut Env) -> R<f64> {
        let f = self.field_of(e, Quantity::Length, "res", env)?;
        let h = self.const_of(e, f)?;
        self.check_cell(e, h)
    }

    /// A joint's drive `v`, written `e`: an angle for a turn, a length for a slide. A bare number is one (degrees,
    /// mm: a pose); a value with a unit must be one, or a length times an angle for a slide (a radius turned
    /// through it, the angle in radians). Anything else is refused, saying what to write: a number times t is a
    /// time, not a place.
    pub(super) fn drives(&mut self, v: &mut Value, e: &Sexp, turn: bool) -> R<()> {
        let Value::Field(f, Some(s)) = &*v else { return Ok(()) };
        let (q, what, unit, rate) = if turn { (Quantity::Angle, "a turn's drive", "deg", "deg/s") } else { (Quantity::Length, "a slide's drive", "mm", "mm/s") };
        let own = s.unit.quantity;
        let like = s.dims().is_some_and(|(d, _)| d[..4] == units::dimension(q).0[..4]);
        if own == q || (own == Quantity::Other && like) {
            let times = s.read_as(q, what).map_err(|m| self.err(e, m))?;
            if times != 1.0 {
                let k = self.g.constant(times);
                let u = units::unit(q.unit()).expect("a quantity's unit is a unit");
                *v = Value::Field(self.g.mul(*f, k), Some(Rc::new(units::Said::worked(u, None, s.explicit_when, &s.file, s.line))));
            }
            return Ok(());
        }
        let is = match s.dims() {
            Some(_) if own == Quantity::Time => format!("a time (its t is written at {}:{})", s.file, s.line),
            Some(_) if own != Quantity::Other => format!("{} ({})", own.name(), s.cite()),
            Some((d, _)) => format!("{} (worked out at {}:{})", units::base_units(d), s.file, s.line),
            None => format!("no quantity (its unit was lost in a function, at {}:{})", s.file, s.line),
        };
        let hint = if own == Quantity::Time {
            format!("t is in seconds, so a number times t is a time: say the rate's unit, (* 360{rate} t), or the travel's, (* 90{unit} (/ t 0.5s))")
        } else {
            format!("say each number's unit where it is written: a rate (* 360{rate} t), a travel (* 90{unit} (sin (* 180deg/s t))), a time (- t 0.3s)")
        };
        Err(self.err_hint(e, format!("{what} takes {} ({unit}), and {e} is {is}", q.name()), Some(hint)))
    }

    /// Value `v` of expression `e` read in a slot that takes quantity `q`
    /// (`what`: the slot, for the message): Err where it, or one of a vec
    /// of it, was said in another quantity's unit; and what is `q` times a
    /// power of angle becomes `q`, its angle in radians (units::Said::read_as).
    pub(super) fn fits(&mut self, v: &mut Value, e: &Sexp, q: Quantity, what: &str) -> R<()> {
        // A bare number is a length in mm (so its square an area and its cube a volume) or a pure number (and, where a
        // pose is taken, degrees): where another
        // quantity is read it must say its unit, or one meant in another unit is read wrong with nothing said.
        // Nothing is no quantity's: a bare 0, and no limit (an infinity), stand.
        if let Value::Field(f, None) = &*v
            && !matches!(q, Quantity::Length | Quantity::Area | Quantity::Volume | Quantity::Angle | Quantity::Other)
            && self.g.as_const(*f).is_none_or(|k| k != 0.0 && k.is_finite())
        {
            let unit = q.unit();
            let fix = match &e.k {
                SK::Num(..) => format!("say its unit ({e}{unit}, if it is in {unit})"),
                _ => format!("say it where its number is written (5{unit}), or of the value worked out: {}", q.say(&e.to_string())),
            };
            return Err(self.err_hint(
                e,
                format!("{what} takes {}, and {e} is a bare number: {fix}", q.name()),
                Some(format!("a bare number is a length in mm or a pure number; {} says its unit: {}", q.name(), units::units_of(q))),
            ));
        }
        match v {
            Value::Field(f, Some(s)) => {
                let times = s.read_as(q, what).map_err(|m| self.err(e, m))?;
                if times != 1.0 {
                    let k = self.g.constant(times);
                    let unit = units::unit(q.unit()).expect("a quantity's unit is a unit");
                    *v = Value::Field(self.g.mul(*f, k), Some(Rc::new(units::Said::worked(unit, None, s.explicit_when, &s.file, s.line))));
                }
            }
            Value::Vec(items) if !items.is_empty() => {
                let mut items = (**items).clone();
                for item in &mut items {
                    self.fits(item, e, q, what)?;
                }
                *v = Value::list(items);
            }
            _ => {}
        }
        Ok(())
    }

    /// Value `v` of `e` is a change of a quantity (an accuracy, a
    /// tolerance's deviation): a unit with a zero of its own (degF, K)
    /// reads a point on its scale (1degF is -17.2 degC), not a change, so
    /// it is refused.
    pub(super) fn fits_change(&self, v: &Value, e: &Sexp, what: &str) -> R<()> {
        let zero = said(v).into_iter().find_map(|s| match s.unit.scale {
            units::Scale::Offset { factor, .. } => Some((s, factor)),
            _ => None,
        });
        match zero {
            Some((s, k)) => Err(self.err_hint(
                e,
                format!("{what} is a change, and {} reads a point on its scale (written {})", s.unit.name, s.cite()),
                Some(format!("write the change in {}: one of 1{} is {}{}", s.unit.quantity.unit(), s.unit.name, (k * 1e6).round() / 1e6, s.unit.quantity.unit())),
            )),
            None => Ok(()),
        }
    }

    /// What param `name`'s value is said by: the unit `--param` gave it in
    /// (`overrides`), which must be the quantity its default (`said`, at
    /// `at`) was said in, else the default's.
    pub(super) fn param_given(&mut self, name: &str, default: f64, said: Option<Rc<units::Said>>, at: &Sexp) -> R<Option<Rc<units::Said>>> {
        let Some((_, given)) = self.overrides.iter().find(|(k, _)| k == name) else { return Ok(said) };
        let Some(unit) = given.unit else {
            let inherited = said.as_deref().and_then(|s| units::UnitRequirement::carried(s, default));
            if let Some(rule) = self.param_requirements.get(name).or(inherited.as_ref()) {
                let value = given.text.as_deref().map(str::to_owned).unwrap_or_else(|| given.value.to_string());
                // A symmetric range's midpoint is zero even when neither end is.
                let zero = given.value == 0.0 && given.text.as_deref().is_none_or(|t| t.split(':').all(|e| e.trim().parse::<f64>() == Ok(0.0)));
                rule.check(name, &format!("{name}={value}"), false, zero).map_err(|m| self.err(at, m))?;
            }
            return Ok(said.map(|s| Rc::new(units::Said { explicit_when: Some(self.g.constant(0.0)), ..(*s).clone() })));
        };
        let v = given.value;
        if let Some(s) = &said
            && s.unit.quantity != unit.quantity
        {
            return Err(self.err(at, format!("--param {name}: {name} is {}, written {}, but {} is {}", s.unit.quantity.name(), s.cite(), unit.name, unit.quantity.name())));
        }
        Ok(Some(Rc::new(units::Said { explicit_when: None, unit, written: units::to(v, unit), file: std::sync::Arc::from(format!("--param {name}")), line: 0, powers: None })))
    }

    /// Expression `e` as a field where a slot takes quantity `q` (`what`:
    /// the slot), its value checked as `fits` does.
    pub(super) fn field_of(&mut self, e: &Sexp, q: Quantity, what: &str, env: &mut Env) -> R<NodeId> {
        let mut v = self.eval(e, env)?;
        self.fits(&mut v, e, q, what)?;
        self.as_field(v, e)
    }

    /// What a builtin's value is said by, and what its number is times to
    /// be in that unit: the one place a number gets its unit from an
    /// operation (apply_builtin calls it for every builtin, however it was
    /// reached: written, applied, folded, reduced or mapped).
    ///
    /// - A product, a quotient, a square, a root and a whole-number power:
    ///   the quantity their operands make (`product_said`).
    /// - A sum, a difference, a least, a greatest, `neg`, `abs`, `floor`
    ///   and `select`'s branches: their terms' one dimension (`one_said`).
    /// - A sine, cosine or tangent of an angle, its number of degrees or
    ///   radians (`rad->deg`, `deg->rad`), and `atan2` of two values of one
    ///   dimension: a pure number.
    /// - Anything else given a value with a unit (an exponential, the sine
    ///   of a length): no quantity, marked so (units::OTHER), never a pure
    ///   number that a later product could give a unit.
    pub(super) fn said_of(&mut self, name: &str, args: &[Value], at: &Sexp) -> R<(Option<Rc<units::Said>>, f64)> {
        let dims: Vec<Option<Option<([i32; 5], f64)>>> = args.iter().map(|a| if let Value::Field(_, Some(s)) = a { Some(s.dims()) } else { None }).collect();
        if dims.iter().all(Option::is_none) {
            return Ok((None, 1.0));
        }
        let (file, line) = (self.file.clone(), at.line);
        let lost = || Some(Rc::new(units::Said::worked(&units::OTHER, None, None, &file, line)));
        Ok(match name {
            "*" | "/" | "square" | "sqrt" | "pow" => return self.product_said(name, args, at),
            "+" | "-" | "neg" | "abs" | "min" | "max" | "emax" | "floor" => (self.one_said(name, args, at)?, 1.0),
            "select" if args.len() == 3 => {
                let said = self.one_said(name, &args[1..], at)?;
                if let (Some(s), [Value::Field(condition, _), a, b]) = (&said, args) {
                    let (left, right) = (self.explicit_unit(a), self.explicit_unit(b));
                    // Both written the one way: that way, whatever the condition (which a fact's value keeps unfolded).
                    let which = if left == right { left } else { self.g.select(*condition, left, right) };
                    let explicit_when = (self.g.as_const(which) != Some(-1.0)).then_some(which);
                    (Some(Rc::new(units::Said { explicit_when, ..(**s).clone() })), 1.0)
                } else {
                    (said, 1.0)
                }
            }
            // Of an angle (degrees, which a turning rate times a time is): a pure number.
            "sin" | "cos" | "tan" | "deg->rad" | "rad->deg" if matches!(dims[0], Some(Some((d, _))) if d == [0, 0, 0, 0, 1]) => (None, 1.0),
            // Of another quantity: its number of degrees is no angle's.
            "deg->rad" | "rad->deg" if let Some(Some(Some(_))) = dims.first() => {
                let Value::Field(_, Some(s)) = &args[0] else { unreachable!("it has a dimension") };
                let x = at.list().and_then(|l| l.get(1)).map_or_else(|| "its value".to_string(), |x| x.to_string());
                return Err(self.err_hint(
                    at,
                    format!(
                        "'{name}' takes an angle, and {x} is {}",
                        if s.unit.quantity == Quantity::Other { "no angle".to_string() } else { format!("{} ({})", s.unit.quantity.name(), s.cite()) }
                    ),
                    Some("an angle is a bare number of degrees, or a number with deg or rad; a rate times t is one, (* 90deg/s t)".into()),
                ));
            }
            // Of a time: a number times t, its rate left bare. Radians a second, turns, degrees? Nothing says.
            "sin" | "cos" | "tan" if matches!(dims[0], Some(Some((d, _))) if d == [0, 0, 1, 0, 0]) => {
                let x = at.list().and_then(|l| l.get(1)).map_or_else(|| "X".to_string(), |x| x.to_string());
                return Err(self.err_hint(
                    at,
                    format!("'{name}' takes an angle, and {x} is a time: t is in seconds, so say how fast the angle turns"),
                    Some(format!("write the rate with its unit, ({name} (* 180deg/s t)) or ({name} (* 2rad/s t)); a bare number here was radians a second")),
                ));
            }
            "atan2" if matches!((dims[0], dims[1]), (Some(Some((a, _))), Some(Some((b, _)))) if a == b) => (None, 1.0),
            _ => (lost(), 1.0),
        })
    }

    /// What a product's, a quotient's, a square's, a square root's or a
    /// whole-number power's value is said by, and what its number is times
    /// to be in that unit. The operands that carry a unit give their
    /// dimensions (units::Said::dims), added and subtracted, with the SI
    /// value of one of each; bare operands are pure numbers. The powers are
    /// carried through whatever is between and judged at the end: a
    /// quantity of the language's is said in its unit (kg times m/s2 is
    /// newtons: the number, worked out in grams, is times 0.001); nothing
    /// left (mm over mm, a metre over a millimetre: 1000) is a pure number;
    /// anything else (N s, seconds per degree) is marked no
    /// quantity, its powers kept and its number as worked out. An angle's
    /// power is one of them, never dropped or converted here, so a product
    /// is one value however it is grouped; the radian appears where the
    /// value is read (units::Said::read_as). A
    /// temperature that is a point on a scale (K, degC, degF) is an error
    /// in any of these; a temperature difference is an amount like any.
    pub(super) fn product_said(&mut self, name: &str, args: &[Value], at: &Sexp) -> R<(Option<Rc<units::Said>>, f64)> {
        let operands = if name == "pow" { &args[..1] } else { args };
        let saids: Vec<Option<&Rc<units::Said>>> = operands.iter().map(|a| if let Value::Field(_, s) = a { s.as_ref() } else { None }).collect();
        if let Some(t) = saids.iter().flatten().find(|s| s.is_point()) {
            let why = format!(
                "a temperature is a point on its scale, not an amount: {} cannot be multiplied, divided or raised, in K, degC or degF alike. A difference of two temperatures can",
                t.cite()
            );
            // A bare number added to a zero is read as kelvins (degC's steps): in degF it is scaled by one of its own.
            let (u, step) = (t.unit.name, |n: &str| if t.unit.steps_as_kelvins() { n.to_string() } else { format!("(* {n} (- 1{} 0{}))", t.unit.name, t.unit.name) });
            let how = format!(
                "a temperature worked out is a bare number with its scale's zero added: (+ 0{u} {}); an amount of one (a rise, what a rate is per) is a difference, which scales: (* 80 (- 1{u} 0{u}))",
                step("(* 80 2)")
            );
            return Err(self.err_hint(at, why, Some(how)));
        }
        let explicit_when = saids.iter().flatten().find_map(|s| s.explicit_when);
        let other = |d: &Doc, powers| Rc::new(units::Said::worked(&units::OTHER, powers, explicit_when, &d.file, at.line));
        // Each operand's power in the result: a quotient's after the first are under; a square, a root, a whole power.
        let whole = match (name, args.get(1)) {
            ("pow", Some(Value::Field(n, None))) => self.g.as_const(*n).filter(|x| x.fract() == 0.0 && x.abs() < 64.0).map(|x| x as i32),
            _ => None,
        };
        let power = |k: usize| -> Option<(i32, i32)> {
            Some(match name {
                "*" => (1, 1),
                "/" => (if k == 0 { 1 } else { -1 }, 1),
                "square" => (2, 1),
                "sqrt" => (1, 2),
                _ => (whole?, 1),
            })
        };
        let (mut dim, mut si) = ([0i32; 5], 1.0f64);
        for (k, said) in saids.iter().enumerate() {
            let Some(said) = said else { continue };
            // A dimension that is not known, or a power ink cannot follow (not whole, or worked out): no quantity.
            let (Some((d, f)), Some((up, under))) = (said.dims(), power(k)) else { return Ok((Some(other(self, None)), 1.0)) };
            for i in 0..5 {
                if (d[i] * up) % under != 0 {
                    return Ok((Some(other(self, None)), 1.0));
                }
                dim[i] = dim[i].saturating_add(d[i] * up / under);
            }
            si *= f.powf(f64::from(up) / f64::from(under));
        }
        let carriers: Vec<&Rc<units::Said>> = saids.iter().flatten().copied().collect();
        let Some(first) = carriers.first() else { return Ok((None, 1.0)) };
        // One quantity scaled by pure numbers: it stays as it was said.
        if carriers.len() == 1 && first.unit.quantity != Quantity::Other && Some(dim) == first.dims().map(|d| d.0) && matches!(name, "*" | "/") {
            return Ok((Some((*first).clone()), 1.0));
        }
        if dim == [0; 5] {
            return Ok((None, si));
        }
        // No quantity of the language's: its number is put in the base units (g, mm, s, K, deg), whatever units its
        // operands were in and however they were grouped (a newton made on the way is 10^6 g mm/s2), so it is one
        // number. An angle left in what is otherwise a quantity stays: it is read in radians where the value is
        // used (units::Said::read_as).
        let Some(quantity) = units::quantity_of(dim) else {
            let base = 1e-3f64.powi(dim[0].saturating_add(dim[1]));
            return Ok((Some(other(self, Some((dim, base)))), si / base));
        };
        // A temperature made by arithmetic is an amount of it: a difference.
        let unit = if quantity == Quantity::Temperature { &units::DIFFERENCE } else { units::unit(quantity.unit()).expect("a quantity's unit is a unit") };
        Ok((Some(Rc::new(units::Said::worked(unit, None, explicit_when, &self.file, at.line))), si / units::dimension(quantity).1))
    }

    pub(super) fn explicit_unit(&mut self, v: &Value) -> NodeId {
        match v {
            Value::Field(_, Some(s)) => s.explicit_when.unwrap_or_else(|| self.g.constant(-1.0)),
            _ => self.g.constant(0.0),
        }
    }

    /// What values that are one quantity (a sum's terms, a range's ends)
    /// are said by. Every term whose dimension is known is of one
    /// dimension, an error at `at` where two differ. The result is said as
    /// the first that is a quantity of the language's (a temperature less a
    /// temperature, a difference; a point where any term is one); else as
    /// the first whose powers are known; else it is marked no quantity. A
    /// bare term is taken in the language's unit, as it always was.
    pub(super) fn one_said(&mut self, what: &str, terms: &[Value], at: &Sexp) -> R<Option<Rc<units::Said>>> {
        let said: Vec<&Rc<units::Said>> = terms.iter().filter_map(|v| if let Value::Field(_, Some(s)) = v { Some(s) } else { None }).collect();
        let Some(any) = said.first() else { return Ok(None) };
        let known: Vec<&Rc<units::Said>> = said.iter().copied().filter(|s| s.dims().is_some()).collect();
        if let Some(first) = known.first()
            && let Some(s) = known.iter().find(|s| s.dims().map(|d| d.0) != first.dims().map(|d| d.0))
        {
            let m = format!("'{what}' takes values of one quantity, but {} is {} and {} is {}", first.cite(), first.unit.quantity.name(), s.cite(), s.unit.quantity.name());
            // Arithmetic keeps an angle's power: one that differs only by it says how to take the angle out.
            let angle = first.dims().zip(s.dims()).is_some_and(|(a, b)| a.0[..4] == b.0[..4]);
            // A time beside what bare numbers made: t is seconds, and the numbers worked with it say no unit.
            let time = [first, s].iter().any(|x| x.unit.quantity == Quantity::Time) && [first, s].iter().any(|x| x.unit.quantity == Quantity::Other);
            let hint = if angle {
                Some("they differ by an angle: divide the one that has it by 1rad to take the angle in radians".to_string())
            } else if time {
                Some("t is in seconds, so what is worked out with it says its units too: a time (0.3s), a rate (1800deg/s, 40mm/s), an acceleration (9000deg/s2, 400mm/s2), a travel (90deg, 100mm)".to_string())
            } else {
                None
            };
            return Err(self.err_hint(at, m, hint));
        }
        // A term whose unit was lost (a function of a value with one) has no dimension to compare: the whole has none
        // either, and never takes another term's unit.
        if said.iter().any(|s| s.dims().is_none()) {
            return Ok(Some(Rc::new(units::Said::worked(&units::OTHER, None, None, &self.file, at.line))));
        }
        let named: Vec<&Rc<units::Said>> = known.iter().copied().filter(|s| s.unit.quantity != Quantity::Other).collect();
        let Some(first) = named.first() else { return Ok(Some((*known.first().unwrap_or(any)).clone())) };
        let mut kept = *first;
        for s in &named {
            // Said in the unit units::kept keeps (degF over degC), a point on a scale over a difference.
            if (kept.is_point() == s.is_point() && units::kept(kept.unit, s.unit).is_ok_and(|u| !std::ptr::eq(u, kept.unit))) || (s.is_point() && !kept.is_point()) {
                kept = s;
            }
        }
        // A bare term beside a quantity is a number nothing says the unit of: `(- t 300)` meant in ms is 300 s. A
        // length and an angle take one (mm, degrees: what a bare number is); 0 and no limit are no quantity's.
        // A zero written with its unit says the unit of what is added to it: (+ 0deg/s rate), (+ 0degC (* 2 x)).
        let says = matches!(what, "+" | "-") && terms.iter().any(|v| matches!(v, Value::Field(f, Some(z)) if z.written == 0.0 && self.g.as_const(*f).is_some()));
        if !says
            && !matches!(first.unit.quantity, Quantity::Length | Quantity::Area | Quantity::Volume | Quantity::Angle)
            && let Some(k) = terms.iter().position(|v| matches!(v, Value::Field(f, None) if self.g.as_const(*f).is_none_or(|x| x != 0.0 && x.is_finite())))
        {
            let unit = first.unit.quantity.unit();
            // A temperature difference is an amount: what is added to it is one too, never a point (degC) that
            // would make the sum a temperature.
            let difference = first.unit.quantity == Quantity::Temperature && !first.is_point();
            let say = |e: &str| if difference { format!("(* {e} (- 1degC 0degC))") } else { first.unit.quantity.say(e) };
            let unsay = |x: &str| if difference { format!("(/ {x} (- 1degC 0degC))") } else { first.unit.quantity.unsay(x) };
            let written = at.list().map(|l| &l[1..]).filter(|l| l.len() == terms.len() || (matches!(what, "select" | "value") && l.len() == terms.len() + 1));
            // a select's condition comes first; a value claim's comparison is between its two
            let term = written.and_then(|l| l.get(if what == "value" { 2 * k } else { k + usize::from(what == "select" && l.len() == terms.len() + 1) }));
            let (shown, fix) = match term {
                Some(e) if matches!(e.k, SK::Num(..)) && difference => (e.to_string(), format!("say it is a difference too: {}", say(&e.to_string()))),
                Some(e) if matches!(e.k, SK::Num(..)) => (e.to_string(), format!("say its unit ({e}{unit}, if it is in {unit})")),
                // A count is a pure number, and no unit makes it one's size.
                Some(e) if e.list().and_then(|l| l.first()).and_then(Sexp::sym) == Some("len") => {
                    (e.to_string(), "len counts a list's items, a pure number; a vector's size, which keeps its unit, is std's (norm V), or (len3 X Y Z)".to_string())
                }
                Some(e) => (
                    e.to_string(),
                    format!(
                        "it says no unit, so it is a length in mm or a pure number. If it is {}, say so: {}; if the other is meant as a plain number, take its unit out: {}",
                        first.unit.quantity.name(),
                        say(&e.to_string()),
                        unsay("X")
                    ),
                ),
                None => ("one of its terms".to_string(), format!("say each term's unit (5{unit})")),
            };
            return Err(self.err_hint(
                at,
                format!(
                    "{} takes values of one quantity, but {} is {} and {shown} is a bare number: {fix}",
                    if what == "value" { "a value claim".to_string() } else { format!("'{what}'") },
                    first.cite(),
                    first.unit.quantity.name()
                ),
                Some(format!("a bare number is a length in mm or a pure number; {} says its unit: {}", first.unit.quantity.name(), units::units_of(first.unit.quantity))),
            ));
        }
        // A temperature less a temperature is how far apart they are: an amount, in kelvins whatever their scales.
        let difference = what == "-" && terms.len() == 2 && named.len() == 2 && named.iter().all(|s| s.is_point());
        // A bare addend is in the native unit. It cannot inherit the
        // explicit-unit exemption from another addend.
        let mut mask = self.g.constant(-1.0);
        for term in terms {
            if matches!(what, "+" | "-" | "min" | "max") && matches!(term, Value::Field(f, None) if self.g.as_const(*f) == Some(0.0)) {
                continue;
            }
            let other = self.explicit_unit(term);
            mask = self.g.max(mask, other);
        }
        let explicit_when = (self.g.as_const(mask) != Some(-1.0)).then_some(mask);
        if difference {
            return Ok(Some(Rc::new(units::Said::worked(&units::DIFFERENCE, None, explicit_when, &self.file, at.line))));
        }
        if explicit_when != kept.explicit_when {
            return Ok(Some(Rc::new(units::Said { explicit_when, ..(**kept).clone() })));
        }
        Ok(Some(kept.clone()))
    }

    pub(super) fn check_cell(&self, e: &Sexp, h: f64) -> R<f64> {
        if !(h > 0.0 && h.is_finite()) {
            return Err(self.err(e, format!("(res h) expects a cell size above 0 and finite, got {h}")));
        }
        Ok(h)
    }
}
