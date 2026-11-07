//! S-expression reader. Positions are tracked so every error can point at source.

use crate::units::{self, Quantity, Unit};
use std::fmt;

#[derive(Clone, Debug)]
pub enum SK {
    /// A number in the language's unit of its quantity, and the unit it
    /// says, if any (`Written`). It prints as written.
    Num(f64, Option<Written>),
    Sym(String),
    Str(String),
    List(Vec<Sexp>),
}

/// The unit a number says: written in it (`1in`: in, and 1), or a value
/// a `for` made a number of, said by an earlier one (units::Said).
#[derive(Clone, Debug)]
pub enum Written {
    Here(&'static Unit, f64),
    From(Box<units::Said>),
}

#[derive(Clone, Debug)]
pub struct Sexp {
    pub k: SK,
    pub line: u32,
    pub col: u32,
}

impl Sexp {
    pub fn sym(&self) -> Option<&str> {
        match &self.k {
            SK::Sym(s) => Some(s),
            _ => None,
        }
    }
    pub fn list(&self) -> Option<&[Sexp]> {
        match &self.k {
            SK::List(v) => Some(v),
            _ => None,
        }
    }
    /// A list's head symbol and the rest: a clause `(name args...)`. The
    /// head is "" where it is not a symbol or the list is empty (a clause
    /// nothing matches); None: not a list.
    pub fn head_args(&self) -> Option<(&str, &[Sexp])> {
        Some(match self.list()?.split_first() {
            Some((h, rest)) => (h.sym().unwrap_or(""), rest),
            None => ("", &[]),
        })
    }

    pub fn num(&self) -> Option<f64> {
        match self.k {
            SK::Num(v, _) => Some(v),
            _ => None,
        }
    }
    pub fn to_json(&self) -> serde_json::Value {
        match &self.k {
            SK::Num(v, _) => serde_json::json!(v),
            SK::Sym(s) | SK::Str(s) => serde_json::json!(s),
            SK::List(v) => serde_json::Value::Array(v.iter().map(|e| e.to_json()).collect()),
        }
    }
}

impl fmt::Display for Sexp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.k {
            SK::Num(v, None) => write!(f, "{v}"),
            SK::Num(_, Some(Written::Here(u, w))) => write!(f, "{w}{}", u.name),
            SK::Num(v, Some(Written::From(s))) => write!(f, "{}{}", units::to(*v, s.unit), s.unit.name),
            SK::Sym(s) => write!(f, "{s}"),
            SK::Str(s) => write!(f, "{s:?}"),
            SK::List(v) => {
                write!(f, "(")?;
                for (i, e) in v.iter().enumerate() {
                    if i > 0 {
                        write!(f, " ")?;
                    }
                    write!(f, "{e}")?;
                }
                write!(f, ")")
            }
        }
    }
}

#[derive(Debug)]
pub struct ParseError {
    pub msg: String,
    pub line: u32,
    pub col: u32,
}

/// The deepest forms nest: far past any model, and short of the stack (each
/// later pass over a form, and dropping it, recurses as deep).
pub const MAX_NESTING: usize = 512;

pub fn parse(src: &str) -> Result<Vec<Sexp>, ParseError> {
    let mut p = Reader { c: src.chars().collect(), i: 0, line: 1, col: 1, depth: 0 };
    let mut out = Vec::new();
    loop {
        p.skip_ws();
        if p.i >= p.c.len() {
            return Ok(out);
        }
        out.push(p.expr()?);
    }
}

struct Reader {
    c: Vec<char>,
    i: usize,
    line: u32,
    col: u32,
    /// The forms open around the one being read.
    depth: usize,
}

impl Reader {
    fn peek(&self) -> Option<char> {
        self.c.get(self.i).copied()
    }

    fn bump(&mut self) -> Option<char> {
        let ch = self.peek()?;
        self.i += 1;
        if ch == '\n' {
            self.line += 1;
            self.col = 1;
        } else {
            self.col += 1;
        }
        Some(ch)
    }

    fn skip_ws(&mut self) {
        while let Some(ch) = self.peek() {
            if ch == ';' {
                while let Some(c) = self.peek() {
                    if c == '\n' {
                        break;
                    }
                    self.bump();
                }
            } else if ch.is_whitespace() {
                self.bump();
            } else {
                break;
            }
        }
    }

    fn err(&self, msg: impl Into<String>, line: u32, col: u32) -> ParseError {
        ParseError { msg: msg.into(), line, col }
    }

    fn expr(&mut self) -> Result<Sexp, ParseError> {
        self.skip_ws();
        let (line, col) = (self.line, self.col);
        match self.peek() {
            None => Err(self.err("unexpected end of input", line, col)),
            Some('(') => {
                if self.depth >= MAX_NESTING {
                    return Err(self.err(format!("forms nest more than {MAX_NESTING} deep"), line, col));
                }
                self.bump();
                self.depth += 1;
                let mut items = Vec::new();
                loop {
                    self.skip_ws();
                    match self.peek() {
                        None => {
                            return Err(self.err(format!("unclosed '(' opened at {line}:{col}"), line, col));
                        }
                        Some(')') => {
                            self.bump();
                            self.depth -= 1;
                            return Ok(Sexp { k: SK::List(items), line, col });
                        }
                        _ => items.push(self.expr()?),
                    }
                }
            }
            Some(')') => Err(self.err("unexpected ')'", line, col)),
            Some('"') => {
                self.bump();
                let mut s = String::new();
                loop {
                    match self.bump() {
                        None => return Err(self.err("unterminated string", line, col)),
                        Some('"') => break,
                        Some('\\') => match self.bump() {
                            Some('n') => s.push('\n'),
                            Some('t') => s.push('\t'),
                            Some(c @ ('"' | '\\')) => s.push(c),
                            Some(c) => return Err(self.err(format!("unknown escape '\\{c}' in a string: the escapes are \\n \\t \\\" \\\\ (a path can use /)"), line, col)),
                            None => return Err(self.err("unterminated string", line, col)),
                        },
                        Some(c) => s.push(c),
                    }
                }
                Ok(Sexp { k: SK::Str(s), line, col })
            }
            Some(_) => {
                let mut s = String::new();
                while let Some(c) = self.peek() {
                    if c.is_whitespace() || matches!(c, '(' | ')' | '"' | ';') {
                        break;
                    }
                    s.push(c);
                    self.bump();
                }
                let k = if looks_numeric(&s) {
                    let (v, u) = number(&s).map_err(|m| self.err(m, line, col))?;
                    let digits = &s[..s.len() - u.map_or(0, |u| u.name.len())];
                    SK::Num(v, u.map(|u| Written::Here(u, digits.parse().unwrap_or(v))))
                } else {
                    SK::Sym(s)
                };
                Ok(Sexp { k, line, col })
            }
        }
    }
}

/// A number as the language reads it: a decimal, optionally followed by a
/// unit (`25.4mm`, `1in`, `1270kg/m3`: units::UNITS), read in the
/// language's unit of its quantity, and that unit. Refused: past the f64
/// range (never silently infinite or zero), an unknown unit.
pub fn number(s: &str) -> Result<(f64, Option<&'static Unit>), String> {
    // The decimal: sign, digits and point, and an exponent where digits
    // follow its e. The rest is the unit.
    let b = s.as_bytes();
    let mut i = usize::from(matches!(b.first(), Some(b'+' | b'-')));
    while i < b.len() && (b[i].is_ascii_digit() || b[i] == b'.') {
        i += 1;
    }
    if i < b.len() && matches!(b[i], b'e' | b'E') {
        let j = i + 1 + usize::from(matches!(b.get(i + 1), Some(b'+' | b'-')));
        if b.get(j).is_some_and(u8::is_ascii_digit) {
            i = j;
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
        }
    }
    let (digits, rest) = s.split_at(i);
    let malformed = || format!("malformed number '{s}'");
    let u = match rest {
        "" => None,
        r if r.starts_with(['e', 'E', '.', '+', '-']) || r.starts_with(|c: char| c.is_ascii_digit()) => return Err(malformed()),
        r => Some(
            units::unit(r)
                .ok_or_else(|| format!("unknown unit '{r}' in '{s}'. Write the number bare (lengths in mm, pressures in MPa, ...) or with one of: {}", units::listing()))?,
        ),
    };
    let v = match u {
        Some(u) => units::read(digits, u),
        None => digits.parse::<f64>().ok(),
    }
    .ok_or_else(malformed)?;
    // A literal past the f64 range is not silently infinite, or zero.
    let mantissa = digits.split(['e', 'E']).next().unwrap_or("");
    if !v.is_finite() || (v == 0.0 && mantissa.bytes().any(|c| matches!(c, b'1'..=b'9')) && !matches!(u.map(|u| u.scale), Some(units::Scale::Offset { .. }))) {
        return Err(format!("number '{s}' is outside the range of a double"));
    }
    Ok((v, u))
}

/// A number written for quantity `q` (a command's option, say): its value
/// in the language's unit, refused if it says another quantity's unit.
pub fn number_of(s: &str, q: Quantity, what: &str) -> Result<f64, String> {
    let (v, u) = number(s)?;
    match u {
        Some(u) if u.quantity != q => Err(units::wrong(u, q, what)),
        _ => Ok(v),
    }
}

fn looks_numeric(s: &str) -> bool {
    let b = s.as_bytes();
    match b[0] {
        b'0'..=b'9' => true,
        b'-' | b'+' | b'.' => b.len() > 1 && (b[1].is_ascii_digit() || (b[1] == b'.' && b.len() > 2 && b[2].is_ascii_digit())),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_nested_forms_and_reports_unclosed() {
        let v = parse("(def a (+ 1 -2.5e1)) ; comment\n(sym \"s\\\"q\")").unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].to_string(), "(def a (+ 1 -25))");
        assert_eq!(v[1].list().unwrap()[1].to_string(), "\"s\\\"q\"");
        assert_eq!(v[1].line, 2);
        assert!(parse("(a (b c)").unwrap_err().msg.contains("unclosed"));
    }

    /// Forms nested past the stack are refused, not a crash (a million `(`
    /// would overflow it).
    #[test]
    fn deep_nesting_is_an_error() {
        let deep = |n: usize| format!("{}{}", "(".repeat(n), ")".repeat(n));
        assert!(parse(&deep(MAX_NESTING)).is_ok());
        let e = parse(&deep(1_000_000)).unwrap_err();
        assert!(e.msg.contains("nest more than 512 deep"), "{}", e.msg);
    }

    #[test]
    fn reads_units_as_the_languages() {
        let v = |s: &str| number(s).unwrap().0;
        // A power of ten moves the decimal point: the same double as the
        // number written in g/mm^3.
        assert_eq!(v("2.7g/cm3").to_bits(), 0.0027f64.to_bits());
        assert_eq!(v("1270kg/m3").to_bits(), 0.00127f64.to_bits());
        assert_eq!(v("7.85e3kg/m3").to_bits(), 0.00785f64.to_bits());
        assert_eq!(v("2.5e-3m").to_bits(), 2.5f64.to_bits());
        assert!((v("0.284lb/in3") - 0.007861).abs() < 1e-6);
        // A unit's quantity is kept; a bare number has none.
        assert_eq!(v("-0.25in"), -6.35);
        assert_eq!(number("5N*m").unwrap().1.unwrap().quantity, Quantity::Torque);
        assert!(number("5").unwrap().1.is_none());
        assert_eq!(v("1e3"), 1000.0);
        assert!(number("1x2").unwrap_err().contains("unknown unit 'x2'"));
        assert!(number("1e").unwrap_err().contains("malformed"));
        assert!(number("1.2.3").unwrap_err().contains("malformed"));
        assert!(number("1e999in").unwrap_err().contains("outside the range"));
        // -459.67degF is 0 K: zero, not out of range.
        assert!((v("-459.67degF") + 273.15).abs() < 1e-9);
        assert!(number_of("50lbf", Quantity::Pressure, "pressure").unwrap_err().contains("lbf is a force"));
        // A number in a unit prints as written, and holds the language's.
        let e = &parse("(force 0 0 -50lbf)").unwrap()[0];
        assert_eq!(e.to_string(), "(force 0 0 -50lbf)");
        assert!((e.list().unwrap()[3].num().unwrap() + 222.411).abs() < 1e-3);
        assert_eq!(number_of("2", Quantity::Pressure, "pressure").unwrap(), 2.0);
    }
}
