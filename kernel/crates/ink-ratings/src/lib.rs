//! Where the claims core meets `.ink`: what a loaded design asks
//! (`require`), what it states (`fact`) and the ratings its parts carry,
//! as the core's `Design` and the formula tool's tables. Nothing is
//! evaluated here but the numbers the design wrote, and nothing at all for
//! a design that asks nothing.
//! A claim is the prover's to answer; its answer is made a fact of the
//! core's kind here (asserts.rs).
//!
//! A part carries ratings as data in its body's `ratings` attribute, an
//! object by name:
//!
//! ```json
//! { "side-load": { "ordered": [{ "read": "radial-load", "on": "shaft" }, { "number": 100, "unit": "N" }],
//!                  "from": "gearbox table" } }
//! ```
//!
//! `ordered` is the core's `Says::Ordered`: each term at most the next. A
//! term is a number with its unit (`"number"` the number, or `{"record":
//! [KEY, ...]}` where the part's record holds it: said once, there, and
//! read from the body's `spec` attribute; a table's `points` alike); a fact to read, on the part or (`on`) on
//! one of its frames; or a table read at a fact:
//! `{ "table": { "x": "rpm", "y": "N*m", "points": [[x, y], ...] }, "at": { "read": "speed", "on": "shaft" } }`
//! (with `"line": true` where the maker gives the lines between the points; with `"drawn": E` where the
//! points are read off the maker's drawn curve to within E of its y: the curve is then its points' line, E
//! either side, an estimate; with neither, the points are measurements each alone, and between two of them
//! the value is anywhere between theirs; with
//! `"below": "first"` where the curve is no lower before its first point than at it, so that value is read
//! there, an estimate; or `"below": "none: WHY"`, the reason it is not).
//! For a part placed at `opener/motor` the rating becomes the rule
//! `side-load on opener/motor`, reading `radial-load on opener/motor shaft`:
//! a part's place is its path, a frame's the part's path, a space and the
//! frame's name, so no frame shares a place with an instance. A table's
//! value is a fact at the part's path, a space and `/N`, which neither is.
//!
//! A rating with `"sizes": true` limits sizes: each fact it reads is read
//! without its sign (a load either way, a speed either way round). One
//! without it compares the numbers as they are.
//!
//! Units are ink's own (`ink_lang::units`), and every number is handed on
//! in the language's unit of its quantity (N, N*mm, deg/s): a rating's by
//! the unit it names, a fact's by the unit ink kept for its value by its
//! one arithmetic, an angle left in a quantity taken in radians (torque
//! times a turning rate is watts). A fact with no unit (a bare number, no
//! quantity of ink's), held
//! against a rating that has one, is refused by the core, not read as if
//! it had it. A number that is none in ink's units is an error.

mod asserts;
mod derived;
pub mod duty;
pub use asserts::fact;
pub use derived::Derived;

use ink_claims::{About, Design, Fact, Rule, Says, Term, Unavailable, Verdict, about};
use ink_formula::{Formula, Formulas, Table};
use ink_lang::view::ModelView;
use ink_lang::{Doc, InkError, Value, command, units};
use serde_json::Value as Json;

/// The model `args` name, as a command's are (FILE, `--param NAME=V` and
/// `NAME=LO:HI`), ready to ask: each param at its value, each tolerance,
/// `(ranged lo hi)` and range given as a range. A param the file does not
/// have is an error (the load's, and the scope's), not the default.
pub fn load(registry: &'static ink_lang::Registry, args: &command::Args) -> Result<ModelView, InkError> {
    let file = args.file.as_ref().ok_or_else(|| InkError::msg("no file"))?;
    let given = command::param_args(args)?.values;
    let mut view = ModelView::load_for_claims(registry, std::path::Path::new(file), &given, 0.0)?;
    let scope = command::claim_scope(&view, args)?;
    view.set_claim_scope(&scope)?;
    Ok(view)
}

/// The design the file's requirements ask about: a rule for each `require`
/// of ratings (named as written, on nothing), the ratings of the parts under them,
/// and every `fact` of the assemblies they name, each at the view's params
/// with the band it has over the ranged ones. Data that is not a rating is
/// an error, never a rating skipped.
pub fn design(view: &mut ModelView) -> Result<(Design, Formulas), String> {
    let (mut facts, mut rules, mut formulas, mut worked, mut least) = (vec![], vec![], vec![], vec![], vec![]);
    let mut tops: Vec<String> = vec![];
    for require in view.doc.requires.clone() {
        // A requirement that states a claim has no rule: the prover offers its fact (asserts.rs).
        let Some(within) = &require.within else { continue };
        let says = if require.excluding.is_empty() { Says::Within(within.clone()) } else { Says::WithinBut(within.clone(), require.excluding.clone()) };
        rules.push(Rule { about: about(&require.name, ""), says, from: format!("{}:{}", require.file, require.line) });
        ratings(&mut view.doc, within, &mut rules, &mut formulas, &mut worked, &mut least)?;
        let top = within.split('/').next().unwrap().to_string();
        if !tops.contains(&top) {
            stated(view, &top, &top, &mut facts);
            tops.push(top);
        }
    }
    // A fact a part works out, stated by the design as well, would be read as stated and the part's formula never asked.
    if let Some(fact) = facts.iter().find(|f| worked.contains(&f.about)) {
        return Err(format!("{} ({}): the part works this out from other facts (its derives); a design does not state it", fact.about, fact.from));
    }
    // The same of what is of a whole duty (`duty::READS`: a mean speed, a duty's load): the duty tool works each
    // out from a fact at every moment, and a typed one would stand in for the motion's (a bearing turning at
    // 1667 deg/s said to turn at 10 rpm passed a life it fails).
    if let Some((fact, of)) = facts.iter().find_map(|f| Some((f, duty::READS.iter().find(|(what, _)| *what == f.about.what)?.1))) {
        return Err(format!(
            "{} ({}): '{}' is worked out over the motion from '{of}' at each moment; a design does not state it. State what it is the duty of, (fact {of} ...), which is held beside the model's",
            fact.about, fact.from, fact.about.what
        ));
    }
    // A fact a rating reads as no less than some least (a required life is not negative): one stated under it is
    // no number of that kind, and the rating reading it has nothing to compare (it would pass any limit).
    for fact in facts.iter_mut() {
        let Some((_, low, unit)) = least.iter().find(|(about, ..)| *about == fact.about) else { continue };
        let said = fact.band.map_or(fact.value, |(low, _)| Some(low));
        if fact.unit == *unit && said.is_some_and(|x| x < *low) {
            fact.notes.push(format!(
                "it is stated as {} {unit}, and what reads it takes it to be at least {low} {unit}: a requirement of this kind is a size (a slip of sign?)",
                said.unwrap_or_default()
            ));
            (fact.held, fact.value, fact.band) = (Some(Verdict::NotAvailable(Unavailable::Outside)), None, None);
        }
    }
    // A rating and a fact of one name at one place would be one thing to the core, and the fact would answer for the rating.
    if let Some(rule) = rules.iter().find(|r| facts.iter().any(|f| f.about == r.about)) {
        return Err(format!("{}: a rating and a fact stated there have this one name; the rating would not be checked. Rename one", rule.about));
    }
    // A frequency held against a turning rate is refused for its unit: the fact says why, and what to write.
    for rule in &rules {
        let Says::Ordered(terms) = &rule.says else { continue };
        let read: Vec<&About> = terms.iter().filter_map(|t| if let Term::Read(a) | Term::Size(a) = t { Some(a) } else { None }).collect();
        let numbers = terms.iter().filter_map(|t| if let Term::Number(_, unit) = t { Some(unit.clone()) } else { None });
        let units: Vec<String> = numbers.chain(facts.iter().filter(|f| read.contains(&&f.about)).map(|f| f.unit.clone())).collect();
        if ["Hz", "deg/s"].iter().all(|u| units.iter().any(|unit| unit == u)) {
            for fact in facts.iter_mut().filter(|f| read.contains(&&f.about) && ["Hz", "deg/s"].contains(&f.unit.as_str())) {
                fact.notes.push(if fact.unit == "Hz" { ink_lang::units::FREQUENCY_AS_RATE } else { ink_lang::units::RATE_AS_FREQUENCY }.to_string());
            }
        }
    }
    // The same of a table read at a turning rate (a motor's curve at its shaft's speed), or at a frequency.
    for formula in &formulas {
        for fact in facts.iter_mut().filter(|f| f.about == formula.at) {
            let hint = match [fact.unit.as_str(), formula.table.x.as_str()] {
                ["Hz", "deg/s"] => ink_lang::units::FREQUENCY_AS_RATE,
                ["deg/s", "Hz"] => ink_lang::units::RATE_AS_FREQUENCY,
                _ => continue,
            };
            if !fact.notes.iter().any(|n| n == hint) {
                fact.notes.push(hint.to_string());
            }
        }
    }
    Ok((Design { facts, rules }, Formulas(formulas)))
}

/// The facts a design itself must state for the ratings of the parts at or under `target`: those a rating
/// reads as what is required of the part (`least`: a life), each at its place, with the unit its rating
/// writes it in. No tool offers one, so an answer that lacks one says the line to write.
pub fn design_states(doc: &mut Doc, target: &str) -> Vec<(About, String)> {
    let mut out = vec![];
    for item in doc.items(target).unwrap_or_default() {
        let Some((_, Value::Data(data))) = item.attrs.iter().find(|(name, _)| name == "ratings") else { continue };
        for rating in data.value().as_object().into_iter().flat_map(|by_name| by_name.values()) {
            for term in rating["ordered"].as_array().into_iter().flatten() {
                if let (Some(what), Some(_), Some(unit)) = (term["read"].as_str(), term.get("least"), term["unit"].as_str()) {
                    let on = term["on"].as_str().map_or(item.label.clone(), |frame| format!("{} {frame}", item.label));
                    out.push((about(what, &on), unit.to_string()));
                }
            }
        }
    }
    out
}

/// The facts the parts at or under `target` take: what their ratings and their formulas read, each at its
/// place (a part's path, or a frame of it), less what a part works out itself. A fact stated at a place that
/// takes none of its name is read by nothing.
pub fn takes(doc: &mut Doc, target: &str) -> Result<Vec<About>, String> {
    let (mut rules, mut formulas, mut worked) = (vec![], vec![], vec![]);
    ratings(doc, target, &mut rules, &mut formulas, &mut worked, &mut vec![])?;
    let mut read: Vec<About> = vec![];
    for rule in &rules {
        let Says::Ordered(terms) = &rule.says else { continue };
        read.extend(terms.iter().filter_map(|t| if let Term::Read(a) | Term::Size(a) = t { Some(a.clone()) } else { None }));
    }
    read.extend(formulas.iter().map(|f| f.at.clone()));
    for item in doc.items(target)? {
        let Some((_, Value::Vec(derives))) = item.attrs.iter().find(|(name, _)| name == "derives") else { continue };
        for (_, reads) in derives.iter().filter_map(Doc::derived_reads) {
            // A duty fact is the duty of a fact the place then takes too (a mean speed, of its speed).
            let duty = reads.iter().filter_map(|name| duty::READS.iter().find(|(what, _)| what == name).map(|(_, of)| of.to_string()));
            read.extend(duty.collect::<Vec<_>>().iter().map(|name| about(name, &item.label)));
            read.extend(reads.iter().map(|name| about(name, &item.label)));
        }
    }
    // Not what a part works out (its formulas' own facts) or a table's limit, which nobody states.
    let made = |a: &About| worked.contains(a) || formulas.iter().any(|f| f.about == *a);
    let mut takes: Vec<About> = vec![];
    for a in read.into_iter().filter(|a| !made(a)) {
        if !takes.contains(&a) {
            takes.push(a);
        }
    }
    Ok(takes)
}

/// What a design states and what its parts take, whether or not anything is required of it (`ink check`):
/// under the file's root (with none, each of its assemblies that none of them places), every stated fact, and the facts each place takes (`takes`).
pub fn stated_and_taken(view: &mut ModelView) -> Result<(Design, Vec<About>), String> {
    // The file's own root, where it has one; else each of its own assemblies that none of them places.
    let doc = &view.doc;
    let own: Vec<&ink_lang::Assembly> = doc.assemblies.iter().filter(|a| doc.is_own_file(&a.file) && !a.instances.is_empty()).collect();
    let placed = |name: &str| own.iter().any(|a| a.instances.iter().any(|i| i.child == name || doc.assembly(&i.child).and_then(|c| c.of.as_deref()) == Some(name)));
    let tops: Vec<String> = match own.iter().find(|a| a.name == "root") {
        Some(root) => vec![root.name.clone()],
        None => own.iter().filter(|a| !placed(&a.name)).map(|a| a.name.clone()).collect(),
    };
    let (mut facts, mut taken) = (vec![], vec![]);
    for top in &tops {
        stated(view, top, top, &mut facts);
        taken.extend(takes(&mut view.doc, top)?);
    }
    Ok((Design { facts, rules: vec![] }, taken))
}

/// What is said of a stated fact nothing reads, `names` what its place takes (`unread`).
pub fn unread_said(fact: &Fact, names: &[String]) -> String {
    let taken = if names.is_empty() { "that place takes no fact".to_string() } else { format!("that place takes: {}", names.join(", ")) };
    format!("the fact '{}' stated on {} ({}) is read by nothing: no rating or formula there takes that name; {taken}", fact.about.what, fact.about.on, fact.from)
}

/// The stated facts of `design` that nothing reads, each with the names its place does take (`takes`): where
/// its own place takes none, those its part takes at its frames, as `NAME on FRAME`. A misspelled fact is one.
pub fn unread<'d>(design: &'d Design, takes: &[About]) -> Vec<(&'d Fact, Vec<String>)> {
    let part = |on: &str| on.split(' ').next().unwrap_or_default().to_string();
    let mut out = vec![];
    for fact in design.facts.iter().filter(|f| !takes.contains(&f.about)) {
        let mut names: Vec<String> = takes.iter().filter(|t| t.on == fact.about.on).map(|t| t.what.clone()).collect();
        if names.is_empty() {
            let here = part(&fact.about.on);
            names = takes.iter().filter(|t| part(&t.on) == here).map(|t| t.on.split_once(' ').map_or(t.what.clone(), |(_, frame)| format!("{} on {frame}", t.what))).collect();
        }
        out.push((fact, names));
    }
    out
}

/// The facts stated on the instances of `assembly`, placed at `path`, and under them.
fn stated(view: &ModelView, assembly: &str, path: &str, facts: &mut Vec<Fact>) {
    let doc = &view.doc;
    let Some(placed) = doc.assembly(assembly) else { return };
    for instance in &placed.instances {
        let path = format!("{path}/{}", instance.name);
        for fact in &instance.facts {
            // Down to the part it is about, then to the frame's innermost name.
            let (mut part, mut path) = (&instance.child, path.clone());
            for under in &fact.under {
                part = &doc.assembly(part).and_then(|a| a.instances.iter().find(|i| i.name == *under)).expect("checked at load").child;
                path = format!("{path}/{under}");
            }
            let place = match &fact.frame {
                Some(frame) => frame_place(doc, part, &path, frame),
                None => path,
            };
            let about = about(&fact.what, &place);
            let no_unit = fact.unit.is_none().then(|| {
                "ink has no unit for this value: it is a bare number, or no quantity ink names (also with its angles taken in radians), or a root or function dropped its unit"
                    .to_string()
            });
            let from = format!("{}:{}", fact.file, fact.line);
            facts.push(match view.stated(fact) {
                Ok((value, band)) => Fact { band, notes: no_unit.into_iter().collect(), ..Fact::stated(about, value, fact.unit.unwrap_or(""), &from) },
                // A value that is no number is no fact.
                Err(why) => Fact { about, from, held: Some(Verdict::NotAvailable(Unavailable::Outside)), notes: vec![format!("its value is no number: {why}")], ..Fact::default() },
            });
        }
        stated(view, &instance.child, &path, facts);
    }
}

/// The place `frame` of the part `child` placed at `path` is, by its
/// innermost name: an assembly's frame written as exactly a frame of one
/// of its instances is that instance's (a motor's `shaft` is its rotor's),
/// which is where a rating on that part reads it.
fn frame_place(doc: &Doc, child: &str, path: &str, frame: &str) -> String {
    let (mut child, mut path, mut frame) = (child, path.to_string(), frame);
    while let Some((_, of, Some(its))) = doc.assembly(child).and_then(|a| a.frame_on.iter().find(|(name, ..)| name == frame)) {
        for instance in of[0].split('/') {
            let Some(next) = doc.assembly(child).and_then(|a| a.instances.iter().find(|i| i.name == instance)) else { return format!("{path} {frame}") };
            (child, path) = (&next.child, format!("{path}/{instance}"));
        }
        let frames = doc.body(child).map(|b| &b.frames).or(doc.assembly(child).map(|a| &a.frames)).expect("an instance's child");
        frame = &frames[*its as usize].0;
    }
    format!("{path} {frame}")
}

/// `x` said in `unit`, in the language's unit of its quantity, and that unit's word.
fn base(x: f64, unit: &str) -> Result<(f64, &'static str), String> {
    let known = units::unit(unit).ok_or_else(|| format!("'{unit}' is no unit of ink's"))?;
    let value = units::from(x, known);
    if !value.is_finite() {
        return Err(format!("{x} {unit} is no number in ink's unit ({})", known.quantity.unit()));
    }
    Ok((value, known.quantity.unit()))
}

/// The ratings of every part under `target`, placed: the rules, and the tables they read.
/// `least`: each fact a rating reads as no less than a number (`{"read": NAME, "least": N, "unit": U}`), with
/// that number in the language's unit.
fn ratings(doc: &mut Doc, target: &str, rules: &mut Vec<Rule>, formulas: &mut Vec<Formula>, works: &mut Vec<About>, least: &mut Vec<(About, f64, String)>) -> Result<(), String> {
    for item in doc.items(target)? {
        // The facts the part works out (its `derives`): theirs to say, not a rating's name nor the design's to state.
        let worked = item.attrs.iter().find(|(name, _)| name == "derives").map(|(_, derives)| ink_lang::Value::strings_first(derives)).unwrap_or_default();
        works.extend(worked.iter().map(|name| about(name, &item.label)));
        let Some((_, carried)) = item.attrs.iter().find(|(name, _)| name == "ratings") else { continue };
        let Value::Data(data) = carried else { return Err(format!("{}: its ratings are not data", item.label)) };
        let Some(by_name) = data.value().as_object() else { return Err(format!("{}: its ratings are not an object by name", item.label)) };
        if let Some(name) = by_name.keys().find(|name| worked.contains(name)) {
            return Err(format!("{}: '{name}' is a rating and a fact the part works out (its derives). Rename one", item.label));
        }
        for (name, rating) in by_name {
            let rule = about(name, &item.label);
            if rules.iter().any(|r| r.about == rule) {
                continue; // under two requirements
            }
            let bad = |what: &str| format!("{}: rating '{name}': {what}", item.label);
            let from = rating["from"].as_str().ok_or_else(|| bad("no `from`"))?;
            let size = rating["sizes"] == true;
            let place = |term: &Json| -> Result<About, String> {
                let what = term["read"].as_str().ok_or_else(|| bad("a term is a number with its unit, a read, or a table at a read"))?;
                Ok(about(what, &term["on"].as_str().map_or(item.label.clone(), |frame| format!("{} {frame}", item.label))))
            };
            // A number or a table's points the part's record already holds is named, not copied: `{"record":
            // [KEY, ...]}` is what the record (the body's `spec` attribute) has there.
            let record = item.attrs.iter().find(|(name, _)| name == "spec").and_then(|(_, spec)| if let Value::Data(spec) = spec { Some(spec.value()) } else { None });
            let held = |named: &Json| -> Result<Json, String> {
                let Some(path) = named.get("record").and_then(Json::as_array) else { return Ok(named.clone()) };
                let record = record.ok_or_else(|| bad("it reads its record, and its body carries none: (attr spec RECORD)"))?;
                let said: Vec<String> = path.iter().map(|k| k.as_str().map_or_else(|| k.to_string(), String::from)).collect();
                let found = path.iter().try_fold(record, |at, key| match (key.as_str(), key.as_u64()) {
                    (Some(key), _) => at.get(key),
                    (_, Some(index)) => at.get(index as usize),
                    _ => None,
                });
                found.cloned().ok_or_else(|| bad(&format!("its record has no {}", said.join(" "))))
            };
            let (mut terms, mut tables) = (vec![], 0);
            for term in rating["ordered"].as_array().ok_or_else(|| bad("no `ordered`"))? {
                // A fact read as no less than a number: what a design requires of the part (a life), a size.
                if let (Some(_), Some(low)) = (term.get("read"), term.get("least")) {
                    let (Some(low), Some(unit)) = (low.as_f64(), term["unit"].as_str()) else { return Err(bad("`least` is a number, with its `unit`")) };
                    let (low, unit) = base(low, unit).map_err(|why| bad(&why))?;
                    least.push((place(term)?, low, unit.into()));
                }
                terms.push(if let (Some(n), Some(unit)) = (held(&term["number"])?.as_f64(), term["unit"].as_str()) {
                    let (n, unit) = base(n, unit).map_err(|why| bad(&why))?;
                    Term::Number(n, unit.into())
                } else if let Some(table) = term.get("table") {
                    let points = held(&table["points"])?;
                    let points = points.as_array().ok_or_else(|| bad("a table has points"))?;
                    let mut read = vec![];
                    for p in points {
                        // Either number of a point may be the record's (a DC motor's line: its stall torque at rest, none at its free speed).
                        let (Some(x), Some(y)) = (held(&p[0])?.as_f64(), held(&p[1])?.as_f64()) else { return Err(bad("a point is [x, y]")) };
                        read.push((x, y));
                    }
                    let points: Vec<(f64, f64)> = read;
                    let (Some(x), Some(y)) = (table["x"].as_str(), table["y"].as_str()) else { return Err(bad("a table names its x and y units")) };
                    let (xu, yu) = (base(0.0, x).map_err(|why| bad(&why))?.1, base(0.0, y).map_err(|why| bad(&why))?.1);
                    let points = points.iter().map(|p| Ok((base(p.0, x)?.0, base(p.1, y)?.0))).collect::<Result<_, String>>().map_err(|why| bad(&why))?;
                    // The table's value is a fact of its own, which the rule reads, at a place no part or frame has.
                    tables += 1;
                    let limit = about(&format!("{name}'s limit"), &format!("{} /{tables}", item.label));
                    // Before its first point: its first value where the part says so (`"below": "first"`), else
                    // nothing, with the part's reason where it gives one (`"below": "none: WHY"`).
                    let below = match table["below"].as_str() {
                        Some("first") => Ok(()),
                        Some(why) => Err(why.strip_prefix("none: ").unwrap_or(why).to_string()),
                        None => Err(String::new()),
                    };
                    // A drawn curve's reading error, in the table's y unit: `"drawn": E`.
                    let drawn = match held(&table["drawn"])?.as_f64() {
                        Some(e) if e >= 0.0 => Some(base(e, y).map_err(|why| bad(&why))?.0),
                        Some(_) => return Err(bad("a drawn curve's reading error is not negative")),
                        None if table.get("drawn").is_some() => return Err(bad("`drawn` is the reading error of a drawn curve's points, a number in the table's y unit")),
                        None => None,
                    };
                    let table = Table { x: xu.into(), y: yu.into(), points, line: table["line"] == true, drawn, below };
                    formulas.push(Formula { about: limit.clone(), table, at: place(&term["at"])?, size, from: from.into() });
                    Term::Read(limit)
                } else if size {
                    Term::Size(place(term)?)
                } else {
                    Term::Read(place(term)?)
                });
            }
            rules.push(Rule { about: rule, says: Says::Ordered(terms), from: from.into() });
        }
    }
    Ok(())
}
