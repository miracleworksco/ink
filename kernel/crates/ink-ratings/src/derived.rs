//! A part's formulas as a tool of the claims core: a fact a part works out
//! from other facts about it, worked out when it is asked for.
//!
//! A body says them in `(attr derives (vec (vec NAME FN UNIT ...) ...))`:
//! the fact NAME, on the part where it is placed, is FN of the facts its
//! parameters name there, each read in the quantity of the UNIT written for
//! it. FN is an ink function, called once when the tool is made
//! (`ink_lang::Doc::formula`), so the fact's units, its band and its
//! numbers that are none are the language's own.
//!
//! Each fact FN reads is asked of the core, so it is what the core says of
//! it from every source: stated on the part, stated from an assembly above
//! it, or offered by another tool. The fact is open with the reason of the
//! first of those that is open (not stated, sources disagree), refused for
//! its units where one is not in the quantity its parameter reads (never
//! read as if it were), and has the band the bands of those it reads give
//! it. Where FN says where it holds (`(vec VALUE OK WHY)`), it is refused
//! elsewhere, with WHY. A design that asks nothing makes no tool and calls
//! no FN.

use ink_claims::{About, Claims, Fact, Tool, Unavailable, Verdict, about, resting};
use ink_lang::view::ModelView;
use ink_lang::{Formula, Stated, Value};
use std::cell::RefCell;

pub struct Derived<'a> {
    view: &'a RefCell<ModelView>,
    formulas: Vec<(About, Formula)>,
}

impl<'a> Derived<'a> {
    /// The tool for the formulas of every part at or under what the design's requirements name, on a design
    /// loaded by `crate::load`, which it shares with the tools asked beside it. Make it before them: each fact a
    /// formula reads gets a variable of its own in the design.
    pub fn new(view: &'a RefCell<ModelView>) -> Result<Derived<'a>, String> {
        let view_now = &mut *view.borrow_mut();
        let mut formulas: Vec<(About, Formula)> = Vec::new();
        for within in view_now.doc.requires.iter().filter_map(|r| r.within.clone()).collect::<Vec<_>>() {
            for item in view_now.doc.items(&within)? {
                let Some((_, derives)) = item.attrs.iter().find(|(name, _)| name == "derives") else { continue };
                let Value::Vec(derives) = derives.clone() else { return Err(format!("{}: (attr derives ...) is a vec of (vec NAME FN UNIT ...)", item.label)) };
                for entry in derives.iter() {
                    let formula = view_now.doc.formula(entry, &item.label)?;
                    let fact = about(&formula.what, &item.label);
                    if formulas.iter().any(|(a, _)| *a == fact) {
                        // Under two requirements: made once. (Its variables are made twice; the second set is unused.)
                        continue;
                    }
                    formulas.push((fact, formula));
                }
            }
        }
        view_now.ctx.params.resize(view_now.doc.params.len(), 0.0);
        Ok(Derived { view, formulas })
    }
}

impl Tool for Derived<'_> {
    fn ask(&self, asked: &About, claims: &mut Claims) -> Option<Fact> {
        let (_, formula) = self.formulas.iter().find(|(a, _)| a == asked)?;
        let reads: Vec<Fact> = formula.reads.iter().map(|(what, ..)| claims.answer(&about(what, &asked.on))).collect();
        let mut fact = resting(asked, &format!("{}:{}", formula.file, formula.line), reads.iter());
        let refuse = |mut fact: Fact, why: Verdict, note: Option<String>| {
            (fact.held, fact.strength) = (Some(why), None);
            fact.notes.extend(note);
            Some(fact)
        };
        // One that reads nothing is a number its part states (a record's limit): as stated, not proven.
        if reads.is_empty() {
            fact.strength = None;
        }
        // What it reads that is open leaves it open, for that reason.
        if let Some(open) = reads.iter().find_map(|f| f.held.filter(|v| v.open())) {
            return refuse(fact, open, None);
        }
        for (read, (what, _, unit)) in reads.iter().zip(&formula.reads) {
            if read.value.is_none() {
                return refuse(fact, Verdict::NotAvailable(Unavailable::Outside), Some(format!("'{what}' is no number")));
            }
            if read.unit != *unit {
                let said = |unit: &str| if unit.is_empty() { "no unit".to_string() } else { unit.to_string() };
                let turn = match [read.unit.as_str(), unit] {
                    ["Hz", "deg/s"] => format!(". {}", ink_lang::units::FREQUENCY_AS_RATE),
                    ["deg/s", "Hz"] => format!(". {}", ink_lang::units::RATE_AS_FREQUENCY),
                    _ => String::new(),
                };
                return refuse(
                    fact,
                    Verdict::NotAvailable(Unavailable::Units),
                    Some(format!("its formula reads '{what}' in {}, and it is stated in {}{turn}", said(unit), said(&read.unit))),
                );
            }
        }
        // The formula at the values read, each over its band.
        let view = &mut *self.view.borrow_mut();
        let saved = view.ctx.clone();
        let mut spans: Vec<(usize, (f64, f64))> = view.ctx.p_spans.iter().filter(|s| !s.derived).map(|s| (s.k as usize, s.span)).collect();
        for (read, (_, k, _)) in reads.iter().zip(&formula.reads) {
            view.ctx.params[*k] = read.value.expect("checked");
            spans.extend(read.band.map(|band| (*k, band)));
        }
        let worked = (|| {
            if reads.iter().any(|f| f.band.is_some()) {
                view.ctx.p_spans = view.scope_spans(&spans).map_err(|e| (false, e.message))?;
            }
            let of = |value| Stated { what: String::new(), under: Vec::new(), frame: None, value, unit: None, file: formula.file.clone(), line: formula.line };
            if let Some((ok, why)) = &formula.only {
                // It holds where OK is 1 for every value read.
                let (at, over) = view.stated(&of(*ok)).map_err(|why| (false, why))?;
                if at != 1.0 || over.is_some_and(|band| band != (1.0, 1.0)) {
                    return Err((true, why.clone()));
                }
            }
            view.stated(&of(formula.value)).map_err(|why| (false, why))
        })();
        view.ctx = saved;
        match worked {
            Ok((value, band)) => {
                (fact.value, fact.band, fact.unit) = (Some(value), band, formula.unit.unwrap_or("").to_string());
                Some(fact)
            }
            Err((said, why)) => refuse(fact, Verdict::NotAvailable(Unavailable::Outside), Some(if said { why } else { format!("its value is no number: {why}") })),
        }
    }
}
