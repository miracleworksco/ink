//! The claims core. A fact is a value with its story. A rule is a
//! requirement: how one fact is worked out from others. One loop answers
//! what is asked: is it stated, can a rule work it out, else ask the tools.
//! Nothing here knows geometry or physics, and nothing runs until asked.
//! Where a tool reads its facts at each of several values of a variable (a
//! motion's moments, a design's load cases), `over` answers at each with
//! the same loop, and gives the worst.

use std::collections::HashMap;
use std::fmt;

/// What a fact is about: a statement or a quantity (`what`), and the part
/// and place it is on (`on`, empty when it is on none).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct About {
    pub what: String,
    pub on: String,
}

pub fn about(what: &str, on: &str) -> About {
    About { what: what.into(), on: on.into() }
}

impl fmt::Display for About {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.on.is_empty() { write!(f, "{}", self.what) } else { write!(f, "{} on {}", self.what, self.on) }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    Yes,
    No,
    /// Not decided, and something could decide it.
    NotSettled(Unsettled),
    /// Nothing here can answer it as asked; asking again will not help.
    NotAvailable(Unavailable),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unsettled {
    /// The design lacks something its author can add.
    Lacks,
    /// A source tried and stopped without deciding.
    Undecided,
    /// The band reaches across the limit.
    Straddles,
    /// Its sources differ. `Claims::sources` shows each.
    Disagree,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unavailable {
    /// No fact states it, no rule works it out and no tool offers it.
    NoSource,
    /// A source offers it, but not for this design or not as asked.
    Outside,
    /// It rests on itself.
    Circular,
    /// It checks nothing: every rating where nothing carries one, or a
    /// comparison of fewer than two terms.
    Empty,
    /// It compares numbers in different units.
    Units,
}

impl Verdict {
    pub fn open(self) -> bool {
        matches!(self, Verdict::NotSettled(_) | Verdict::NotAvailable(_))
    }
}

impl fmt::Display for Verdict {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str(match self {
            Verdict::Yes => "yes",
            Verdict::No => "no",
            Verdict::NotSettled(Unsettled::Lacks) => "not settled: the design lacks something",
            Verdict::NotSettled(Unsettled::Undecided) => "not settled: a source could not decide",
            Verdict::NotSettled(Unsettled::Straddles) => "not settled: the band reaches across the limit",
            Verdict::NotSettled(Unsettled::Disagree) => "not settled: sources disagree",
            Verdict::NotAvailable(Unavailable::NoSource) => "not available: nothing states it or offers it",
            Verdict::NotAvailable(Unavailable::Outside) => "not available: outside what its source covers",
            Verdict::NotAvailable(Unavailable::Circular) => "not available: it rests on itself",
            Verdict::NotAvailable(Unavailable::Empty) => "not available: it checks nothing",
            Verdict::NotAvailable(Unavailable::Units) => "not available: it compares different units",
        })
    }
}

/// How strong a fact is, weakest first. A stated fact has a source and no
/// strength, and neither has what rests only on stated facts: proven is
/// for a tool that proves. In that order a stated fact comes after an
/// estimate (`rank`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strength {
    /// From a method that can be wrong.
    Estimated,
    /// A statement proven only by the allowance of the accuracy it was
    /// asked to.
    ToAccuracy,
    /// Shown to hold, under the model `from` names.
    Proven,
}

/// Over how much of one variable a fact was shown, least first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Covers {
    One,
    /// At this many values; nothing is claimed between them.
    These(usize),
    Every,
}

/// A value and its story. An answer to a requirement is one of these.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Fact {
    pub about: About,
    /// The number, when it has one, and the unit it is in. The core only
    /// asks whether two units are the same word.
    pub value: Option<f64>,
    pub unit: String,
    /// Whether it held, when it is a statement; or why it has no answer.
    pub held: Option<Verdict>,
    /// The numbers inside its statement, by name.
    pub terms: Vec<(String, f64)>,
    /// The low and high its source stands behind. Accuracy is this and
    /// nothing else: a number with no band is exact, every comparison is of
    /// bands, and the accuracy a statement was asked to is one of its terms.
    pub band: Option<(f64, f64)>,
    pub strength: Option<Strength>,
    /// For each variable it reads, how much of it was covered.
    pub covers: Vec<(String, Covers)>,
    /// What it did not look at that is not a variable.
    pub left_out: Vec<String>,
    /// Where it came from: a line, a datasheet, a tool and its model.
    pub from: String,
    /// A comparison's room: the least distance from a term's high end up
    /// to the next term's low end, in its unit (under 0 where they cross).
    pub margin: Option<f64>,
    /// The facts it rests on.
    pub used: Vec<About>,
    pub notes: Vec<String>,
}

impl Fact {
    pub fn stated(about: About, value: f64, unit: &str, from: &str) -> Fact {
        Fact { about, value: Some(value), unit: unit.into(), from: from.into(), ..Fact::default() }
    }

    pub fn open(&self) -> bool {
        self.held.is_some_and(Verdict::open)
    }

    /// Where its number lies: its band, or exactly its value.
    pub fn range(&self) -> Option<(f64, f64)> {
        if self.open() { None } else { self.band.or(self.value.map(|v| (v, v))) }
    }
}

/// A requirement: how the fact `about` is worked out from other facts. One
/// on a part is a rating: data the part carries until something asks.
#[derive(Clone, Debug, PartialEq)]
pub struct Rule {
    pub about: About,
    pub says: Says,
    pub from: String,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Says {
    /// Each term is at most the next: `x <= limit`, `lo <= x <= hi`.
    Ordered(Vec<Term>),
    /// Every rule on the part at this path, or on a part under it, holds.
    Within(String),
    /// The same, but for the parts at or under these paths, left out on purpose: no rule on them is asked.
    WithinBut(String, Vec<String>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum Term {
    /// A number and its unit.
    Number(f64, String),
    Read(About),
    /// A fact's size: its number without its sign.
    Size(About),
}

/// What the language hands up: the stated facts and the rules.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Design {
    pub facts: Vec<Fact>,
    pub rules: Vec<Rule>,
}

/// A plug. It returns `None` for what it does not offer, else a fact, which
/// may be open with its reason. It may ask for facts while it works.
pub trait Tool {
    fn ask(&self, about: &About, claims: &mut Claims) -> Option<Fact>;
    /// Whether what it gives is the design's own model worked out (a solver's reading of its motion and
    /// loads), and then which way such a fact is the worse: what a stated fact about the same thing is held
    /// beside (`Claims::beside`), so that no typed number makes an answer better than the model's.
    fn measures(&self) -> Option<Worse> {
        None
    }
}

/// Which way a fact the model works out is the worse for what reads it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Worse {
    /// A demand: a load, a torque, a speed, a stress. More of it is worse.
    More,
    /// What there is to meet a demand: a life, a margin. Less of it is worse.
    Less,
}

/// One round of asking. `tools` are tried in order, so the caller puts the
/// cheapest first. A round ends: a tool is asked for a thing once, and
/// again only if it worked its answer out from other facts and one of the
/// round's facts has since been given by a later tool; that happens at
/// most once for each thing and each tool after the first that gave it.
pub struct Claims<'a> {
    design: &'a Design,
    tools: &'a [&'a dyn Tool],
    /// What has been worked out, as the facts stand now: forgotten whole
    /// when a later tool gives one of them.
    found: HashMap<About, Vec<Fact>>,
    /// What is being found now, outermost first, and the outermost of them
    /// that something inside asked for again.
    finding: Vec<About>,
    circle: usize,
    /// Each tool's reply for each thing it was asked, and how many facts
    /// had been given by a later tool when it was asked.
    replies: HashMap<(usize, About), (Option<Fact>, usize)>,
    /// For each fact a tool gave, which: the tools before it are passed
    /// over for it from then on.
    gave: HashMap<About, usize>,
    /// How many times a later tool has given a fact in place of an earlier one's.
    moved: usize,
    /// Whether a tool is being asked only to be held beside a stated fact (`beside`), and whether what one
    /// worked out has joined a stated fact as a source: a case's tool has then given something the answer rests on.
    weighing: bool,
    contended: bool,
}

impl<'a> Claims<'a> {
    pub fn new(design: &'a Design, tools: &'a [&'a dyn Tool]) -> Self {
        Claims {
            design,
            tools,
            found: HashMap::new(),
            finding: vec![],
            circle: usize::MAX,
            replies: HashMap::new(),
            gave: HashMap::new(),
            moved: 0,
            weighing: false,
            contended: false,
        }
    }

    /// The fact about `about`. Sources that differ in what they say (yes or
    /// no, value, unit or band) leave it not settled; sources that say the
    /// same give one answer, as strong as the weakest of them.
    pub fn answer(&mut self, about: &About) -> Fact {
        let sources = self.sources(about);
        let says = |f: &Fact| (f.held, f.value, f.unit.clone(), f.band);
        if sources[1..].iter().any(|s| says(s) != says(&sources[0])) {
            // Each source, with what it says and where it is from.
            let mut notes: Vec<String> = sources.iter().map(|s| format!("{}: {}", if s.from.is_empty() { "a source" } else { &s.from }, said(s))).collect();
            notes.extend(sources.iter().flat_map(|s| s.notes.iter().filter(|n| n.starts_with("the model works out ")).cloned()));
            return Fact { about: about.clone(), held: Some(Verdict::NotSettled(Unsettled::Disagree)), notes, ..Fact::default() };
        }
        let mut fact = sources[0].clone();
        sources[1..].iter().for_each(|s| weaken(&mut fact, s));
        fact
    }

    /// Every source of a fact, never empty: each stated fact; else what
    /// each rule about it works out; else the tools in the order given,
    /// from the one that gave it if one has, until one settles it (if none
    /// does, the first reply); else "no source".
    pub fn sources(&mut self, about: &About) -> Vec<Fact> {
        if let Some(found) = self.found.get(about) {
            return found.clone();
        }
        let open = |why| vec![Fact { about: about.clone(), held: Some(Verdict::NotAvailable(why)), ..Fact::default() }];
        if let Some(outer) = self.finding.iter().position(|a| a == about) {
            self.circle = self.circle.min(outer);
            return open(Unavailable::Circular);
        }
        self.finding.push(about.clone());
        let (design, moved) = (self.design, self.moved);
        let mut found: Vec<Fact> = design.facts.iter().filter(|f| &f.about == about).cloned().collect();
        if !found.is_empty() {
            self.beside(about, &mut found);
        }
        if found.is_empty() {
            found = design.rules.iter().filter(|r| &r.about == about).map(|r| self.work(r)).collect();
        }
        if found.is_empty() {
            for i in self.gave.get(about).copied().unwrap_or(0)..self.tools.len() {
                match self.reply(i, about) {
                    Some(fact) if !fact.open() => {
                        found = vec![fact];
                        self.gave.insert(about.clone(), i);
                        break;
                    }
                    Some(fact) if found.is_empty() => found.push(fact),
                    _ => {}
                }
            }
        }
        if found.is_empty() {
            found = open(Unavailable::NoSource);
        }
        // A yes or a number that rests on an open fact is open, for its reason: no line says yes beside one
        // that is not settled or not available. (A no stands on the one thing that fails; and what nothing
        // states or offers is no fact here, only a name its source cites: an old answer's inputs.)
        let unknown = Verdict::NotAvailable(Unavailable::NoSource);
        for fact in found.iter_mut().filter(|f| f.held == Some(Verdict::Yes) || (f.held.is_none() && f.value.is_some())) {
            if let Some(why) = fact.used.clone().iter().find_map(|used| self.answer(used).held.filter(|v| v.open() && *v != unknown)) {
                (fact.held, fact.value, fact.band, fact.strength) = (Some(why), None, None, None);
                fact.unit.clear();
            }
        }
        self.finding.pop();
        // What was found while something outside it was still being found
        // may rest on that thing's stand-in, and what was found while a
        // fact moved to a later tool may rest on the earlier one's: not kept.
        match self.circle.cmp(&self.finding.len()) {
            std::cmp::Ordering::Less => return found,
            std::cmp::Ordering::Equal => self.circle = usize::MAX,
            std::cmp::Ordering::Greater => {}
        }
        if self.moved == moved {
            self.found.insert(about.clone(), found.clone());
        }
        found
    }

    /// The stated facts `found` about `about`, held beside what the first tool that measures the model
    /// (`Tool::measures`) works out about it in the same unit, so that a typed number never makes an answer
    /// better than the model's. Sizes are compared, toward the worse way the tool says (`Worse`; in the words
    /// below, of a demand, where more is worse):
    /// - stated at or past the model's most, it stands: an allowance is the designer's to state;
    /// - stated short of the model's least, the model's fact joins it as a source, and they disagree;
    /// - stated anywhere inside the model's band, it does not narrow it toward the better: the fact is the
    ///   model's, its band from the worse of the two least to the worse of the two most (a designer may be
    ///   more wary than the model, never less), and says what was stated and where.
    ///
    /// Stated facts that differ among themselves disagree as they are: nothing is weighed.
    fn beside(&mut self, about: &About, found: &mut Vec<Fact>) {
        let says = |f: &Fact| (f.held, f.value, f.unit.clone(), f.band);
        if found[1..].iter().any(|f| says(f) != says(&found[0])) {
            return;
        }
        self.weighing = true;
        let measuring: Vec<(usize, Worse)> = (0..self.tools.len()).filter_map(|i| Some((i, self.tools[i].measures()?))).collect();
        let model = measuring.into_iter().find_map(|(i, worse)| Some((self.reply(i, about).filter(|f| !f.open() && f.value.is_some())?, worse)));
        self.weighing = false;
        let Some((mut model, worse)) = model else { return };
        // A range's sizes: the least and the most its number can be, without its sign.
        let sizes = |f: &Fact| f.range().map(|(low, high)| (if low <= 0.0 && high >= 0.0 { 0.0 } else { low.abs().min(high.abs()) }, low.abs().max(high.abs())));
        let Some((least, most)) = sizes(&model) else { return };
        // Each stated fact in the model's unit: whether it leaves room, falls short, or lies inside; and the
        // worst end any inside one reaches to.
        let (mut short, mut inside, mut reach, mut wary) = (false, vec![], if worse == Worse::More { most } else { least }, if worse == Worse::More { least } else { most });
        for stated in found.iter().filter(|f| f.unit == model.unit) {
            let Some((low, high)) = sizes(stated) else { continue };
            let (room, lacks) = match worse {
                Worse::More => (low >= most * (1.0 - 1e-9), high < least * (1.0 - 1e-9)),
                Worse::Less => (high <= least * (1.0 + 1e-9), low > most * (1.0 + 1e-9)),
            };
            short |= lacks;
            if !room && !lacks {
                (reach, wary) = if worse == Worse::More { (reach.max(high), wary.max(low)) } else { (reach.min(low), wary.min(high)) };
                inside.push(format!("stated as {} ({}): inside what the model works out, which a typed number does not narrow", said(stated), stated.from));
            }
        }
        if short {
            let more = if worse == Worse::More { "more" } else { "less" };
            model.notes.push(format!(
                "the model works out {more} than is stated: state no {} than it, or leave the fact out to be answered from the model",
                if worse == Worse::More { "less" } else { "more" }
            ));
            found.push(model);
            self.contended = true;
        } else if !inside.is_empty() {
            // The model's band, each end the worse of its own and the stated one's: as sizes, then with its sign.
            let (low, high) = model.range().expect("it has sizes");
            let (near, far) = if worse == Worse::More { (wary, reach) } else { (reach, wary) };
            let band = if least == 0.0 && near == 0.0 {
                (if low.abs() > high.abs() { -far } else { low }, if low.abs() > high.abs() { high } else { far })
            } else if high.abs() >= low.abs() {
                (high.signum() * near, high.signum() * far)
            } else {
                (low.signum() * far, low.signum() * near)
            };
            if band != (low, high) || model.band.is_some() {
                model.band = Some(band);
                // Its number is no better than its band's better end.
                model.value = model.value.map(|v| v.clamp(band.0, band.1));
            }
            model.notes.extend(inside);
            *found = vec![model];
            self.contended = true;
        }
    }

    /// The facts the design states about `about`, as stated: for a tool whose reading rests on what the
    /// design says of something else (how a load statics leaves open is shared).
    pub fn stated(&self, about: &About) -> Vec<&Fact> {
        self.design.facts.iter().filter(|f| &f.about == about).collect()
    }

    /// Whether a tool is being asked only to be held beside a stated fact: what it says then is no part of
    /// the answer unless it contends.
    pub fn weighing(&self) -> bool {
        self.weighing
    }

    /// Tool `i`'s reply about `about`: asked once, and again only where it
    /// worked its reply out from other facts and a fact has moved since.
    fn reply(&mut self, i: usize, about: &About) -> Option<Fact> {
        let key = (i, about.clone());
        if let Some((reply, at)) = self.replies.get(&key)
            && (*at == self.moved || reply.as_ref().is_none_or(|f| f.used.is_empty()))
        {
            return reply.clone();
        }
        let at = self.moved;
        let reply = self.tools[i].ask(about, self);
        // One worked out while something outside it was being found may rest on that thing's stand-in: not kept.
        if self.circle >= self.finding.len() {
            self.replies.insert(key, (reply.clone(), at));
        }
        reply
    }

    fn work(&mut self, rule: &Rule) -> Fact {
        // Worked again while a fact it rests on moved to a later tool as it was worked.
        loop {
            let moved = self.moved;
            let fact = match &rule.says {
                Says::Ordered(terms) => {
                    fn reads(t: &Term) -> Option<&About> {
                        if let Term::Read(a) | Term::Size(a) = t { Some(a) } else { None }
                    }
                    let read: Vec<Option<Fact>> = terms.iter().map(|t| reads(t).map(|a| self.answer(a))).collect();
                    let fact = ordered(rule, terms, &read);
                    // Too close to call on what a tool gave: the next tool that settles one of the facts gives it.
                    if moved == self.moved && fact.held == Some(Verdict::NotSettled(Unsettled::Straddles)) {
                        terms.iter().filter_map(reads).any(|a| self.next_tool(a, &mut vec![]));
                    }
                    fact
                }
                Says::Within(_) | Says::WithinBut(..) => {
                    let (path, but): (&String, &[String]) = match &rule.says {
                        Says::WithinBut(path, but) => (path, but),
                        Says::Within(path) => (path, &[]),
                        Says::Ordered(_) => unreachable!("a rule of ratings"),
                    };
                    let at = |on: &str, path: &str| on == path || on.strip_prefix(path).is_some_and(|rest| rest.starts_with('/'));
                    let under = |on: &str| at(on, path) && !but.iter().any(|left| at(on, left));
                    let carried: Vec<&About> = self.design.rules.iter().map(|r| &r.about).filter(|a| under(&a.on) && **a != rule.about).collect();
                    let answers: Vec<Fact> = carried.into_iter().map(|a| self.answer(a)).collect();
                    // An answer that is neither yes nor no (a number under a rating's name) is no answer.
                    let held: Vec<Verdict> = answers.iter().map(|a| a.held.unwrap_or(Verdict::NotAvailable(Unavailable::Outside))).collect();
                    // A no is decided by one rating that fails, the strongest of them; a yes by all of them.
                    let fails = answers.iter().filter(|a| a.held == Some(Verdict::No)).max_by_key(|a| rank(a.strength));
                    let mut fact = resting(&rule.about, &rule.from, fails.map_or(&answers[..], std::slice::from_ref).iter());
                    fact.used = answers.iter().map(|a| a.about.clone()).collect();
                    let first = |want: fn(&Verdict) -> bool| held.iter().copied().find(want);
                    let none = held.is_empty().then_some(Verdict::NotAvailable(Unavailable::Empty));
                    settle(&mut fact, none.or(first(|v| *v == Verdict::No)).or(first(|v| v.open())).unwrap_or(Verdict::Yes));
                    fact
                }
            };
            if moved == self.moved {
                return fact;
            }
        }
    }

    /// Whether the fact about `a`, as found, is or rests on the one about `on`.
    fn rests_on(&self, a: &About, on: &About, seen: &mut Vec<About>) -> bool {
        a == on
            || (!seen.contains(a) && {
                seen.push(a.clone());
                self.found.get(a).into_iter().flatten().flat_map(|f| &f.used).any(|u| self.rests_on(u, on, seen))
            })
    }

    /// Has the first tool after the one that gave the fact about `about`
    /// that settles it give it, or failing that does so for a fact it
    /// used. True where one did: everything worked out is then forgotten,
    /// to be worked out again from the facts as they now stand.
    fn next_tool(&mut self, about: &About, seen: &mut Vec<About>) -> bool {
        if seen.contains(about) {
            return false;
        }
        seen.push(about.clone());
        let used: Vec<About> = self.found.get(about).into_iter().flatten().flat_map(|f| f.used.clone()).collect();
        for i in self.gave.get(about).map_or(self.tools.len(), |i| i + 1)..self.tools.len() {
            // Not to a tool that works it out from a fact that rests on it: that would be a circle.
            if self.reply(i, about).is_some_and(|f| !f.open() && !f.used.iter().any(|u| self.rests_on(u, about, &mut vec![]))) {
                self.gave.insert(about.clone(), i);
                self.moved += 1;
                self.found.clear();
                // Its reply was worked out from the facts as they stand after the move.
                if let Some(kept) = self.replies.get_mut(&(i, about.clone())) {
                    kept.1 = self.moved;
                }
                return true;
            }
        }
        used.iter().any(|u| self.next_tool(u, seen))
    }
}

/// A case's tool as `over` asks it: whether it gave a fact, and its facts
/// without their coverage of the variable (inside its case each is the
/// whole of it; what an answer then still says of the variable is what the
/// other facts it rests on cover).
struct Case<'a>(&'a dyn Tool, &'a str, std::cell::Cell<bool>);

impl Tool for Case<'_> {
    fn measures(&self) -> Option<Worse> {
        self.0.measures()
    }

    fn ask(&self, about: &About, claims: &mut Claims) -> Option<Fact> {
        let mut fact = self.0.ask(about, claims)?;
        // What it could not say, and what it is asked only to be weighed beside a stated fact, is no part of
        // an answer over the cases.
        if !claims.weighing() && !fact.open() {
            self.2.set(true);
        }
        fact.covers.retain(|(variable, _)| variable != self.1);
        Some(fact)
    }
}

/// The fact about `asked` where a tool reads its facts at each of several
/// values of `variable` (a motion's moments, a design's load cases): `at`
/// gives that tool at case 0, 1, ... and None past the last, and it is
/// asked in `tools`' order before the tool at `slot`. The fact is answered
/// at every case by the one loop, and the answer over them is the worst
/// case's: a no if any is no, else open if any is open, else yes; among
/// those the one with the least margin (one with no margin first: how bad
/// it is is not known), its terms and notes that case's. With which case
/// it was. It is as strong as the weakest case and leaves out what any
/// did. It covers `variable` at as many values as there were cases, and
/// no further than any other fact it rests on covers it: a limit known at
/// one moment is known at one. Where no case's tool gave anything the
/// answer reads, it does not vary: it is the first's, as the tools alone
/// give it.
pub fn over<'t>(design: &Design, asked: &About, variable: &str, tools: &[&dyn Tool], slot: usize, at: &dyn Fn(usize) -> Option<Box<dyn Tool + 't>>) -> (Fact, usize) {
    let mut answers: Vec<Fact> = vec![];
    let mut gave = false;
    for case in 0.. {
        let Some(tool) = at(case) else { break };
        let case = Case(tool.as_ref(), variable, std::cell::Cell::new(false));
        let with: Vec<&dyn Tool> = tools[..slot].iter().copied().chain([&case as &dyn Tool]).chain(tools[slot..].iter().copied()).collect();
        let mut claims = Claims::new(design, &with);
        answers.push(claims.answer(asked));
        gave |= case.2.get() || claims.contended;
    }
    if !gave {
        return (answers.into_iter().next().unwrap_or_else(|| Claims::new(design, tools).answer(asked)), 0);
    }
    // Worst first: a no, then what is open, then the rest; in each, no margin before the least.
    let badness = |f: &Fact| {
        let rank = match f.held {
            Some(Verdict::No) => 0,
            Some(v) if v.open() => 1,
            _ => 2,
        };
        (rank, f.margin.is_some(), f.margin.unwrap_or(0.0))
    };
    let order = |a: &Fact, b: &Fact| {
        let (a, b) = (badness(a), badness(b));
        (a.0, a.1).cmp(&(b.0, b.1)).then(a.2.total_cmp(&b.2))
    };
    let worst = (0..answers.len()).min_by(|&a, &b| order(&answers[a], &answers[b])).expect("a case gave a fact");
    let mut fact = answers[worst].clone();
    // Which of the values this is, said: not the largest demand nor the lowest limit, but where they come
    // nearest, each at the end of its band that is the worse.
    if let Some(room) = fact.margin.filter(|_| answers.len() > 1 && fact.held == Some(Verdict::Yes)) {
        fact.notes.push(format!(
            "shown at the one of the {} values of {variable} that leaves the least room, {}, from the top of what is asked to the bottom of its limit, their bands counted: \
             a larger demand elsewhere has more room under its limit there",
            answers.len(),
            (room * 1e3).round() / 1e3
        ));
        // And the largest of what is asked, where that is another value's: the peak a reader looks for.
        let first = |f: &Fact| f.terms.first().map(|(name, value)| (name.clone(), value.abs()));
        let most = answers.iter().filter_map(first).max_by(|a, b| a.1.total_cmp(&b.1));
        if let (Some((name, most)), Some((_, here))) = (most, first(&fact))
            && most > here
        {
            fact.notes.push(format!("the largest {name} at any of them is {}, in the unit shown below", (most * 1e3).round() / 1e3));
        }
    }
    answers.iter().for_each(|other| weaken(&mut fact, other));
    if fact.open() {
        fact.strength = None;
    }
    // As many values as cases, and no more than what else it rests on covers.
    let cases = if answers.len() == 1 { Covers::One } else { Covers::These(answers.len()) };
    let others = fact.covers.iter().find(|(v, _)| v == variable).map(|(_, covers)| *covers);
    fact.covers.retain(|(v, _)| v != variable);
    fact.covers.push((variable.into(), others.map_or(cases, |others| others.min(cases))));
    (fact, worst)
}

/// Where a strength stands, weakest first: an estimate, a stated number
/// (none), proven to an accuracy, proven.
fn rank(strength: Option<Strength>) -> u8 {
    strength.map_or(1, |s| [0, 2, 3][s as usize])
}

/// Numbers no more than four floating-point steps apart are equal: what a
/// unit's conversion or a table's arithmetic rounds away. Nothing looser.
pub fn same(a: f64, b: f64) -> bool {
    // Finite numbers of one sign are in the order of their bits, a step apart each.
    a == b || (a.is_finite() && b.is_finite() && (a > 0.0) == (b > 0.0) && a != 0.0 && b != 0.0 && a.to_bits().abs_diff(b.to_bits()) <= 4)
}

/// A fact worked out from others, before its value or verdict: it rests
/// on them, is as strong as the weakest (proven where all are, or there
/// are none), covers each variable as the least of them does, and leaves
/// out what any of them did.
pub fn resting<'f>(about: &About, from: &str, facts: impl Iterator<Item = &'f Fact>) -> Fact {
    let mut fact = Fact { about: about.clone(), from: from.into(), strength: Some(Strength::Proven), ..Fact::default() };
    for f in facts {
        fact.used.push(f.about.clone());
        weaken(&mut fact, f);
    }
    fact
}

/// Leaves `fact` no stronger than `other`, covering no more, leaving out no less.
/// What a fact says, in a few words: its number and unit with its band, or why it has none.
fn said(fact: &Fact) -> String {
    let unit = if fact.unit.is_empty() { String::new() } else { format!(" {}", fact.unit) };
    let number = |x: f64| ((x * 1e3).round() / 1e3).to_string();
    match (fact.held.filter(|v| v.open()), fact.value, fact.band) {
        (Some(why), ..) => why.to_string(),
        (None, Some(value), Some((low, high))) => format!("{}{unit}, between {} and {}{unit}", number(value), number(low), number(high)),
        (None, Some(value), None) => format!("{}{unit}", number(value)),
        (None, None, _) => fact.held.map_or("nothing".to_string(), |v| v.to_string()),
    }
}

fn weaken(fact: &mut Fact, other: &Fact) {
    if rank(other.strength) < rank(fact.strength) {
        fact.strength = other.strength;
    }
    for (var, covers) in &other.covers {
        match fact.covers.iter_mut().find(|(v, _)| v == var) {
            Some((_, ours)) => *ours = (*ours).min(*covers),
            None => fact.covers.push((var.clone(), *covers)),
        }
    }
    // (A motion read in thousands of steps leaves out thousands of stretches, the same at every moment:
    // looked up, not searched.)
    if other.left_out.is_empty() || other.left_out == fact.left_out {
        return;
    }
    let ours: std::collections::HashSet<&str> = fact.left_out.iter().map(String::as_str).collect();
    let new: Vec<String> = other.left_out.iter().filter(|l| !ours.contains(l.as_str())).cloned().collect();
    fact.left_out.extend(new);
}

fn settle(fact: &mut Fact, held: Verdict) {
    fact.held = Some(held);
    if held.open() {
        fact.strength = None;
    }
}

/// `terms` in order, each read term as the fact found for it.
fn ordered(rule: &Rule, terms: &[Term], read: &[Option<Fact>]) -> Fact {
    let mut fact = resting(&rule.about, &rule.from, read.iter().flatten());
    let (mut ranges, mut units) = (vec![], vec![]);
    for (term, read) in terms.iter().zip(read) {
        let range = match (term, read) {
            (Term::Number(n, unit), _) => {
                units.push(unit);
                Some((*n, *n))
            }
            (Term::Read(a) | Term::Size(a), Some(f)) => {
                fact.terms.extend(f.value.map(|v| (a.to_string(), v)));
                units.push(&f.unit);
                // A size: the least and greatest distance from nothing in its range.
                f.range().map(|(low, high)| match term {
                    Term::Size(_) => (if low <= 0.0 && high >= 0.0 { 0.0 } else { low.abs().min(high.abs()) }, low.abs().max(high.abs())),
                    _ => (low, high),
                })
            }
            _ => None,
        };
        let Some(range) = range else {
            // An open fact passes its reason up; a statement is no number.
            settle(&mut fact, read.as_ref().and_then(|f| f.held).filter(|v| v.open()).unwrap_or(Verdict::NotAvailable(Unavailable::Outside)));
            return fact;
        };
        ranges.push(range);
    }
    if ranges.len() >= 2 && units.windows(2).all(|w| w[0] == w[1]) {
        // Ends that are the same number, both without end too, are no distance apart; a distance that is no
        // number is no margin.
        let rooms: Vec<f64> = ranges.windows(2).map(|w| if w[1].0 == w[0].1 { 0.0 } else { w[1].0 - w[0].1 }).collect();
        fact.margin = rooms.iter().copied().min_by(f64::total_cmp).filter(|_| rooms.iter().all(|room| !room.is_nan()));
    }
    let held = if ranges.len() < 2 {
        Verdict::NotAvailable(Unavailable::Empty)
    } else if units.windows(2).any(|w| w[0] != w[1]) {
        Verdict::NotAvailable(Unavailable::Units)
    } else if ranges.windows(2).any(|w| w[0].0 > w[1].1 && !same(w[0].0, w[1].1)) {
        Verdict::No
    } else if ranges.windows(2).all(|w| w[0].1 <= w[1].0 || same(w[0].1, w[1].0)) {
        Verdict::Yes
    } else {
        Verdict::NotSettled(Unsettled::Straddles)
    };
    settle(&mut fact, held);
    fact
}
