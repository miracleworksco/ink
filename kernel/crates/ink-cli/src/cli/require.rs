//! `ink require FILE`: what the file's requirements answer. Each one: yes,
//! no, not settled or not available, with why, how strong, over how much
//! of what it reads, and what it rests on. JSON, as every command answers;
//! `--plain` prints the same answer in words. A file that requires nothing
//! answers nothing and runs nothing.
//!
//! A claim (`(require NAME CLAIM)`, and those a limit clause or a relation
//! makes) is the prover's: every claim of the file is checked in its one
//! run (`claims`), at the params' values and at each declared case, and its
//! row carries the prover's evidence (a witness, an undecided cell, hints)
//! beside the verdict its answer is as a fact of the core
//! (`ink_ratings::fact`).
//!
//! A requirement of ratings (`within-ratings`) is the core's: every tool
//! this distribution has is asked, the cheapest first: what the design
//! states (the core's own first source), the tables its parts carry
//! (ink-formula), the formulas its parts carry (ink-ratings' Derived),
//! then the distribution's solvers (`Solvers`: the kernel has none;
//! simulation hands its motion's efforts and its studies), each started
//! only by a fact nothing before it settled.
//!
//! A requirement that reads a motion's torque, force or speed is answered
//! at every moment the motion's run read (the core's `over`): it is yes
//! where it is yes at each, and what is under it is the moment with the
//! least room, each rating its own.

use crate::*;
use ink_claims::{About, Claims, Covers, Design, Fact, Says, Strength, Term, Tool, Verdict, about, over};
use std::cell::RefCell;

/// A distribution's solver tools for the requirements of one assembly.
pub trait Solvers {
    /// Its tools, asked after what is stated, tabled, worked out and proven.
    fn tools(&self) -> Vec<&dyn Tool>;
    /// The tool for its motion's moment `moment` (0, 1, ...), None past
    /// the last (and at every one where it has no motion).
    fn at(&self, moment: usize) -> Option<Box<dyn Tool + '_>>;
    /// The time each of its motion's moments weighs in a mean over it (the step it is the middle of; a
    /// span's end, none), in their order; empty where it has no motion. Asking starts its run.
    fn steps(&self) -> Vec<f64> {
        vec![]
    }
    /// How far its answers reach, asked once the requirements are answered: it starts no run.
    fn coverage(&self) -> Coverage {
        Coverage::default()
    }
}

/// How far a solver's model says its answers reach (`Solvers::coverage`): what a yes over a path does not
/// say unless it is said.
#[derive(Default)]
pub struct Coverage {
    /// The span of t its run was over, where it ran over one.
    pub span: Option<(f64, f64)>,
    /// Why the written motion is not over at the span's end, where the model can tell.
    pub goes_on: Option<String>,
    /// The parts its model puts a demand on that a maker rates (what drives, what a coupling passes force
    /// through, what is named as holding something), by instance path, each with why.
    pub loaded: Vec<(String, String)>,
    /// The joints that hold a coupling's pull, by instance path, each with why: a load across or along the
    /// joint, which a rating of its part checks only if it reads one.
    pub holds: Vec<(String, String)>,
    /// The loads stated that its model does not apply, each said.
    pub loads_left: Vec<String>,
    /// What the motion needs that nothing in the design supplies or holds (an effort no actuator gives, a
    /// body moved through no joint, a load on a part on no joint), each said: the answers leave it out.
    pub unheld: Vec<String>,
    /// What its answers assume because the design does not say (a coupling with no efficiency: lossless),
    /// each said once. An answer resting on one is an estimate.
    pub assumed: Vec<String>,
    /// Where the solver's own books do not balance (what it handed out against what is needed, power in
    /// against work done), each said: its answers there are not to be relied on.
    pub unbalanced: Vec<String>,
    /// How many joints at how many moments had their own balance (weight, loads, pulls, mass times
    /// acceleration against the reaction) added up apart from the solve.
    pub balanced: usize,
    /// And how many of those had their moments (about the world's origin, with the rate their angular
    /// momentum changes at) added up too.
    pub balanced_moments: usize,
}

/// What a requirement of ratings covered and did not: of the parts under its path, of its motion, of the
/// loads stated.
#[derive(Default)]
struct Covered {
    span: Option<(f64, f64)>,
    goes_on: Option<String>,
    /// The parts under the path that carry a rating.
    rated: Vec<String>,
    /// Those the model loads that carry none, each with why it counts; and those of them the design leaves out on purpose.
    unrated: Vec<(String, String)>,
    excluded: Vec<String>,
    /// Those of the loaded parts with no rating whose record says its maker states no limit (a body's
    /// `no-rating` attribute, in its words): nothing can check them, which is said apart from an oversight.
    states_none: Vec<(String, String)>,
    /// The parts that carry ratings and hold a coupling's pull that none of their ratings reads (a motor rated
    /// for its torque, under its belt's pull across its shaft), each with why.
    pull_unread: Vec<(String, String)>,
    /// How many other bodies under the path carry no rating: the model puts no demand on them a maker rates.
    others: usize,
    loads_left: Vec<String>,
    unheld: Vec<String>,
    assumed: Vec<String>,
    unbalanced: Vec<String>,
    balanced: usize,
    balanced_moments: usize,
}

impl Covered {
    /// Whether a yes over the path is a yes about all of it: no loaded part unrated, the motion over, every load applied.
    fn whole(&self) -> bool {
        self.unrated.is_empty()
            && self.states_none.is_empty()
            && self.pull_unread.is_empty()
            && self.goes_on.is_none()
            && self.loads_left.is_empty()
            && self.unheld.is_empty()
            && self.unbalanced.is_empty()
    }
}

/// What makes a distribution's solvers for the assembly named, on the
/// design every tool shares (`ink_cli::main`'s `solvers`).
pub type Plug = for<'a> fn(&'a RefCell<ModelView>, &str) -> Result<Box<dyn Solvers + 'a>, String>;

/// A fact asked, as answered: its line, and what it rests on.
struct Node {
    asked: About,
    /// "the size of " where a rating reads it without its sign.
    lead: &'static str,
    fact: Fact,
    /// The numbers inside a statement no rule here made.
    terms: Vec<(String, f64)>,
    kind: Kind,
    under: Vec<Under>,
    /// What its reader should look at though the answer stands: a fact stated under it that nothing reads.
    warnings: Vec<String>,
    /// Where the prover answered it: its evidence as it gives it (a witness, an undecided cell, hints, its
    /// value and accuracy, how it held over a motion), by key.
    proof: Option<J>,
    /// Of a requirement of ratings: what it covered.
    covered: Option<Covered>,
}

enum Kind {
    /// A comparison: each of `under` at most the next.
    Ordered,
    /// Every rating under a path: `under` is each.
    Ratings,
    /// Anything else: `under` is what it used.
    Used,
}

enum Under {
    Number(f64, String),
    Fact(Node),
}

/// How a fact is asked for: of the tools, and of one motion's moments (if
/// any) in their place among them.
struct Asking<'a> {
    design: &'a Design,
    tools: &'a [&'a dyn Tool],
    slot: usize,
    motion: Option<&'a dyn Solvers>,
}

pub(crate) fn cmd_require(view: Session, a: &Args) -> Result<J, InkError> {
    let flag = |k: &str| match a.get(k) {
        None => Ok(false),
        Some("true") => Ok(true),
        Some(v) => Err(InkError::msg(format!("--{k} takes no value, got '{v}' (write the file first: ink require FILE --{k})"))),
    };
    let (imports, plain) = (flag("imports")?, flag("plain")?);
    let mut view = view;
    // The file's own requirements are what is asked unless --imports: of ratings as of claims. (An imported
    // assembly's requirement of its own motor is that file's to ask, with its own motion and gravity.)
    let all = view.doc.requires.len();
    if !imports {
        let own: Vec<bool> = view.doc.requires.iter().map(|r| r.within.is_none() || view.doc.is_own_file(&r.file)).collect();
        let mut own = own.into_iter();
        view.doc.requires.retain(|_| own.next().unwrap());
    }
    let left = all - view.doc.requires.len();
    // --name: one requirement, of ratings or a claim.
    let only = a.get("name");
    let rated = |name: &str| view.doc.requires.iter().any(|r| r.within.is_some() && r.name == name);
    let (ratings, proved) = (only.is_none_or(rated), !only.is_some_and(rated));
    // The claims, by the prover: at the params' values and at each combination of the cases they declare,
    // each claim's answer the worst of those.
    let (mut nodes, mut said) = (vec![], json!({}));
    // The parts no requirement of ratings asks over (`answer`).
    let mut unasked: Vec<(String, String)> = vec![];
    // Whether its parts are clear of each other at the start of its motion (`placement`), where ratings are asked.
    let mut placed: Option<Result<Placed, String>> = None;
    if proved {
        let found = claims(&mut view, a, only, !imports, &|answer, params| row(&proven(answer, params)))?;
        let (answers, params) = &found.at_values;
        nodes = answers.iter().map(|answer| proven(answer, params)).collect();
        for (at, answers, params) in &found.at_cases {
            worst(&mut nodes, answers.iter().map(|answer| proven(answer, params)).collect(), at.clone());
        }
        if let Some(cases) = &found.cases {
            said_of_cases(&mut nodes, cases);
        }
        // Each row's evidence is the claim's row of the run: with how it came out at each case, where there are any.
        said = found.out;
        let rows = said["requirements"].take();
        for node in &mut nodes {
            node.proof = rows.as_array().and_then(|rows| rows.iter().find(|r| r["name"] == node.asked.what.as_str())).cloned().or(node.proof.take());
        }
        // The file's own claims are what is asked unless --imports: how many of its imports' were left, where any.
        let fields = said.as_object_mut().expect("an answer");
        fields.remove("requirements");
        fields.remove("scope");
        if fields.get("imported_not_run") == Some(&json!(0)) {
            fields.remove("imported_not_run");
        }
    }
    if left > 0 {
        said["imported_not_run"] = json!(said["imported_not_run"].as_u64().unwrap_or(0) + left as u64);
    }
    // The requirements of ratings, by the core and its tools: at the params' values and at each case alike.
    if ratings && view.doc.requires.iter().any(|r| r.within.is_some() && only.is_none_or(|n| n == r.name)) {
        // A requirement of ratings is answered over a span, the file's or --t-span's: a moment asked is not read, and is refused.
        if a.get("t").is_some() {
            return Err(InkError::msg("ink require answers a requirement of ratings over the file's (t-span) or --t-span A:B, not at one moment: --t is not read"));
        }
        let mut at_cases: Vec<(serde_json::Map<String, J>, Session)> = vec![];
        let cases = each_case(
            &view,
            a,
            |_| {},
            |_| {},
            |at, v, _| {
                at_cases.push((at.clone(), v));
                Ok(())
            },
        )?;
        placed = Some(placement(&mut view));
        let (mut rated, not_asked) = answer(view, only)?;
        unasked = not_asked;
        if let Some(cases) = &cases {
            for (at, v) in at_cases {
                worst(&mut rated, answer(v, only)?.0, at);
            }
            said_of_cases(&mut rated, cases);
            said["cases"] = json!(cases.asked);
            for (key, those) in [("cases_not_allowed", &cases.not_allowed), ("cases_not_loaded", &cases.not_loaded)] {
                if !those.is_empty() {
                    said[key] = json!(those);
                }
            }
        }
        rated.append(&mut nodes);
        nodes = rated;
    } else if !proved {
        let names = view.doc.requires.iter().map(|r| r.name.as_str());
        return Err(InkError::msg(format!("no requirement named '{}'. {}", only.unwrap_or_default(), ink::lang::choices(only.unwrap_or_default(), "Its requirements", names))));
    }
    let mut counts = serde_json::Map::new();
    for node in &nodes {
        let said = node.fact.held.map_or("no answer", verdict);
        counts.insert(said.into(), json!(counts.get(said).and_then(J::as_u64).unwrap_or(0) + 1));
    }
    // A yes shown at one value of something it reads (one moment of a motion no span covers) is a yes there,
    // and the file is not all yes.
    let one: Vec<&str> =
        nodes.iter().filter(|n| n.fact.held == Some(Verdict::Yes) && n.fact.covers.iter().any(|(_, c)| *c == Covers::One)).map(|n| n.asked.what.as_str()).collect();
    // Every answer yes, and at every value asked of what it reads; and each of them proven, none an estimate
    // or only as stated. With nothing asked there is nothing to say of all: null, and why.
    let complete = said.get("cases_not_loaded").is_none();
    // A yes over a path that left a part the model loads unrated, a stated load unapplied, or the rest of a
    // motion unread is a yes about less than was asked: said by name, and the file is not all yes.
    let part: Vec<String> = nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.unrated.iter().map(|(part, _)| part.clone())).collect();
    let on: Vec<&str> = nodes.iter().filter(|n| n.covered.as_ref().is_some_and(|c| c.goes_on.is_some())).map(|n| n.asked.what.as_str()).collect();
    let loads: Vec<&String> = nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| &c.loads_left).collect();
    let unheld: Vec<String> = nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.unheld.clone()).collect();
    let clashes = placed.as_ref().and_then(|p| p.as_ref().ok()).map_or(0, |p| p.clashes.len());
    PLACEMENT_SAID.with(|said| {
        *said.borrow_mut() = match &placed {
            Some(Ok(_)) if clashes == 0 => "no clash at the start pose".to_string(),
            Some(Ok(_)) => format!("{clashes} clash{} in how the parts are placed", if clashes == 1 { "" } else { "es" }),
            _ => "placement not checked".to_string(),
        }
    });
    let whole = clashes == 0 && nodes.iter().all(|n| n.covered.as_ref().is_none_or(Covered::whole));
    let all_yes = complete && whole && unasked.is_empty() && one.is_empty() && nodes.iter().all(|n| n.fact.held == Some(Verdict::Yes));
    let proven = all_yes && nodes.iter().all(|n| matches!(n.fact.strength, Some(Strength::Proven | Strength::ToAccuracy)));
    let mut out = if nodes.is_empty() && complete {
        json!({"all_yes": null, "all_proven": null, "counts": counts, "nothing_required": "nothing was asked: the file states no requirement (require forms, limit clauses, relations), or none was asked for"})
    } else {
        json!({"all_yes": all_yes, "all_proven": proven, "counts": counts})
    };
    if !one.is_empty() {
        out["yes_at_one_value"] = json!(one);
    }
    let none: Vec<String> = nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.states_none.iter().map(|(part, _)| part.clone())).collect();
    let pulls: Vec<String> = nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.pull_unread.iter().map(|(part, _)| part.clone())).collect();
    match &placed {
        Some(Ok(p)) => {
            out["placement"] = json!({"checked": "at the start of the motion", "clashes": p.clashes, "loose": p.loose});
            out["not_checked_here"] = json!([PLACED_AFTER]);
        }
        Some(Err(why)) => out["not_checked_here"] = json!([format!("{PLACEMENT} (not checked here: {why})")]),
        None if nodes.iter().any(|n| n.asked.what.starts_with("relate:")) => out["not_checked_here"] = json!([PLACEMENT]),
        None => {}
    }
    for (key, those) in [
        ("parts_under_no_requirement", json!(unasked.iter().map(|(part, why)| json!({"part": part, "why": why})).collect::<Vec<_>>())),
        ("parts_not_rated", json!(dedup(part))),
        ("parts_with_no_limit_stated", json!(dedup(none))),
        ("loads_no_rating_reads", json!(dedup(pulls))),
        ("books_do_not_balance", json!(dedup(nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.unbalanced.clone()).collect()))),
        ("not_driven_or_held", json!(dedup(unheld))),
        ("assuming", json!(dedup(nodes.iter().filter_map(|n| n.covered.as_ref()).flat_map(|c| c.assumed.clone()).collect()))),
        ("motion_goes_on", json!(on)),
        ("loads_not_applied", json!(dedup(loads.into_iter().cloned().collect()))),
    ] {
        if those.as_array().is_some_and(|a| !a.is_empty()) {
            out[key] = those;
        }
    }
    // The limits a package declares that this run cannot check, by name, said once for the answer.
    let unchecked: Vec<(&str, &str)> = nodes.iter().filter_map(|n| Some((n.asked.what.as_str(), n.proof.as_ref()?["not_checked"].as_str()?))).collect();
    unchecked_fields(&mut out, unchecked.into_iter());
    extend(&mut out, said);
    if plain {
        let mut text = String::new();
        // Whether the parts are where they would be built is no requirement's to say: said first, with what checks it.
        match &placed {
            Some(Ok(p)) if p.clashes.is_empty() && p.loose.is_empty() => text += PLACED_CLEAR,
            Some(Ok(p)) if p.clashes.is_empty() => text += &format!("PLACEMENT: no two parts overlap at the start of the motion. {} {PLACED_AFTER}", p.loose()),
            Some(Ok(p)) => {
                text +=
                    &format!("NOT COVERED: the parts are not placed as built: {}. {}{PLACED_AFTER}", p.clashes.join("; "), p.loose() + if p.loose.is_empty() { "" } else { " " })
            }
            Some(Err(why)) => text += &format!("{PLACEMENT} (not checked here: {why})"),
            None if nodes.iter().any(|n| n.asked.what.starts_with("relate:")) => text += PLACEMENT,
            None => {}
        }
        if placed.is_some() || nodes.iter().any(|n| n.asked.what.starts_with("relate:")) {
            text += "\n";
        }
        if !unasked.is_empty() {
            let each: Vec<String> = unasked.iter().map(|(part, why)| format!("{part} ({why})")).collect();
            text += &format!(
                "NOT COVERED by any requirement: {} part{} of the design, so a yes to each requirement is no yes about the design: {}. Ask over the whole of it, \
                 (require NAME (within-ratings ROOT)), or leave a part out on purpose, (excluding PART ...)\n",
                unasked.len(),
                if unasked.len() == 1 { "" } else { "s" },
                each.join("; ")
            );
        }
        let mut said = std::collections::HashSet::new();
        nodes.iter().for_each(|n| words(&mut text, n, 0, &mut said));
        out["plain"] = json!(text);
    }
    // --keys: of each requirement, its name and what is picked.
    out["requirements"] = J::Array(Records::new(a).picked(nodes.iter().map(row).collect(), &["name"], &[])?);
    Ok(out)
}

/// The most parts a requirement's placement is checked for unasked (`placement`): a small machine's clash
/// check at one pose takes a second or two; a robot of 800 parts takes minutes, and is left to `ink interfere`.
const PLACEMENT_PARTS: usize = 60;

/// Whether the design's parts are clear of each other where the file places them, at the start of its
/// motion: `ink interfere`'s own check (boxes first, then the prover on the pairs that may meet), at that one
/// pose. Ok: the pairs found overlapping, each said. Err: why it was not checked (more than
/// `PLACEMENT_PARTS` parts, a file with no root). The view is left as it was.
fn placement(view: &mut Session) -> Result<Placed, String> {
    // A machine that moves: a file with a span of t. (A rig of parts asked as stated, side by side or one over
    // another to compare them, is no assembly to check.)
    if view.motion.or(view.doc.t_span).is_none() {
        return Err("the file states no motion, (t-span A B): it is read as parts asked as stated, not as a machine assembled".into());
    }
    let items = ink::lang::interfere::default_items(view).map_err(|e| e.message)?;
    if items.len() > PLACEMENT_PARTS {
        return Err(format!("it places {} parts, over the {PLACEMENT_PARTS} whose placement is checked unasked", items.len()));
    }
    let (motion, t) = (view.motion, view.ctx.t);
    let start = motion.or(view.doc.t_span).map_or(t, |(a, _)| a);
    let check = |view: &mut Session| -> Result<Placed, InkError> {
        view.set_motion(None)?;
        view.set_moment(start)?;
        let target = view.doc.target("root");
        let relations: Vec<(String, String)> = target.as_ref().map_or_else(Vec::new, |t| view.doc.relations_of(t)).into_iter().map(|r| (r.a, r.b)).collect();
        let loose = ink::lang::loose::loose(view, &items, &relations, LOOSE_GAP)?;
        let plan = ink::lang::interfere::plan(view, items, target.as_ref(), None, None)?;
        let pairs = plan.run(view, None, &|_| None)?;
        let clashes = pairs.iter().filter(|p| p.overlaps()).map(|p| format!("{} and {} overlap in solid", p.items[0], p.items[1])).collect();
        Ok(Placed { clashes, loose })
    };
    let found = check(view).map_err(|e| e.message);
    // (As it was: its moment, then its motion.)
    let back = view.set_moment(t).and_then(|()| view.set_motion(motion)).map_err(|e| e.message);
    found.and_then(|found| back.map(|()| found))
}

/// What the placement check found at the start of the motion: the pairs overlapping in solid, each said, and
/// the groups of parts that touch nothing of the rest (ink::lang::loose). A clash is not covered and the file
/// is not all yes; a loose part is named, and changes no verdict: a rig may hang a load where its arm would be.
struct Placed {
    clashes: Vec<String>,
    loose: Vec<Vec<String>>,
}

impl Placed {
    /// The loose groups in words, or nothing.
    fn loose(&self) -> String {
        if self.loose.is_empty() {
            return String::new();
        }
        let each: Vec<String> = self.loose.iter().map(|g| if g.len() == 1 { g[0].clone() } else { format!("{} (touching only each other)", g.join(", ")) }).collect();
        let n: usize = self.loose.iter().map(Vec::len).sum();
        let (s, es, it, they) = if n == 1 { ("", "es", "it", "it was") } else { ("s", "", "them", "they were") };
        format!(
            "LOOSE: {n} part{s} touch{es} no other part, none within {LOOSE_GAP} mm and no relation joining {it}: {}. Nothing holds {it}, or what does is no part in \
             the file (a pin, a standoff, a belt): a requirement answers for {it} as if {they} held.",
            each.join("; "),
        )
    }
}

/// How near two parts come to touch, mm: a running fit's clearance is less, and a part left where a standoff
/// or pin held it, more.
const LOOSE_GAP: f64 = 1.0;

/// The first line where the parts were found clear of each other at the start of the motion, and what is left.
const PLACED_CLEAR: &str = "PLACEMENT: no two parts overlap at the start of the motion. None is proven more than 1 mm from all the rest. Not checked here: the rest \
    of the motion (`ink interfere FILE` proves it), and what holds each part where it touches (`ink render FILE --out look.png` shows it)";
const PLACED_AFTER: &str = "Not checked here: the rest of the motion (`ink interfere FILE`), and what holds each part where it touches (`ink render FILE --out look.png`)";

/// What no requirement here checks, said first in every answer that asks ratings or relates parts.
const PLACEMENT: &str = "PLACEMENT NOT CHECKED: a requirement answers for the parts as they are joined in the file, not for where they sit. Parts may pass through \
    each other, hang in mid-air or miss the face they mount on (gear teeth at their phase too) and every number below stand. `ink interfere FILE` proves the parts clear of \
    each other over the motion; `ink render FILE --out look.png` shows them";

thread_local! {
    /// What a headline says of placement in this answer: set once its check has run or been skipped.
    static PLACEMENT_SAID: RefCell<String> = RefCell::new("placement not checked".into());
}

/// A requirement's first line after its verdict: what the verdict is of. For a requirement of ratings, how
/// many ratings it checked and how much it did not cover, assumed or left unchecked; for an open one, the
/// first thing under it with no answer, and why. Empty where a yes or a no is the whole of it.
fn headline(node: &Node) -> String {
    let Some(covered) = &node.covered else { return String::new() };
    let fact = &node.fact;
    if fact.held.is_some_and(|v| v.open()) {
        // The first open thing, followed down to what has no answer.
        fn first_open(node: &Node) -> Option<&Node> {
            let open = node.under.iter().filter_map(|u| if let Under::Fact(n) = u { Some(n) } else { None }).find(|n| n.fact.held.is_some_and(|v| v.open()));
            match open {
                Some(n) => first_open(n).or(Some(n)),
                None => None,
            }
        }
        let open = node.under.iter().filter(|u| matches!(u, Under::Fact(n) if n.fact.held.is_some_and(|v| v.open()))).count();
        return match first_open(node) {
            Some(leaf) => {
                let why = leaf.fact.held.map(|v| v.to_string()).unwrap_or_default();
                let why = why.split_once(": ").map_or(why.as_str(), |(_, why)| why);
                format!(
                    ": {open} of its {} ratings {} no answer, first for want of `{}` ({why}); the line to write is under it",
                    node.under.len(),
                    if open == 1 { "has" } else { "have" },
                    leaf.asked
                )
            }
            None if !covered.unbalanced.is_empty() => ": its ratings answered, and the solver's books do not balance (below)".into(),
            None => String::new(),
        };
    }
    let not_covered =
        covered.unrated.len() + covered.states_none.len() + covered.pull_unread.len() + covered.unheld.len() + covered.loads_left.len() + usize::from(covered.goes_on.is_some());
    let mut parts = vec![format!("for the {} rating{} checked", node.under.len(), if node.under.len() == 1 { "" } else { "s" })];
    if not_covered > 0 {
        parts.push(format!("{not_covered} thing{} not covered", if not_covered == 1 { "" } else { "s" }));
    }
    if !covered.assumed.is_empty() {
        parts.push(format!("{} assumed", covered.assumed.len()));
    }
    parts.push(PLACEMENT_SAID.with(|said| said.borrow().clone()));
    if !covered.whole() {
        parts.push("not all yes".into());
    }
    format!(": {}", parts.join("; "))
}

fn dedup(mut lines: Vec<String>) -> Vec<String> {
    let mut seen = std::collections::HashSet::new();
    lines.retain(|l| seen.insert(l.clone()));
    lines
}

/// A requirement's answer over its cases is the worst: where it is worse at the case `at` than so far (a no
/// over an open answer over a yes), that case's answer, saying which; one made only there, its own.
fn worst(nodes: &mut Vec<Node>, there: Vec<Node>, at: serde_json::Map<String, J>) {
    let rank = |n: &Node| match n.fact.held {
        Some(Verdict::Yes) => 0,
        Some(Verdict::No) => 2,
        _ => 1,
    };
    let case = J::Object(at).to_string();
    for mut node in there {
        match nodes.iter_mut().find(|n| n.asked == node.asked) {
            Some(here) if rank(&node) > rank(here) => {
                node.fact.notes.push(format!("the worst of its cases: at the case {case}"));
                // No stronger than it was shown at any case.
                node.fact.strength = weaker(node.fact.strength, here.fact.strength);
                *here = node;
            }
            Some(here) => here.fact.strength = weaker(here.fact.strength, node.fact.strength),
            None => {
                node.fact.notes.push(format!("made only at the case {case}, not at the params' values"));
                nodes.push(node);
            }
        }
    }
}

/// The weaker of two strengths: an estimate, then what is only stated (none), then proven to an accuracy, then proven.
fn weaker(a: Option<Strength>, b: Option<Strength>) -> Option<Strength> {
    let rank = |s: Option<Strength>| match s {
        Some(Strength::Estimated) => 0,
        None => 1,
        Some(Strength::ToAccuracy) => 2,
        Some(Strength::Proven) => 3,
    };
    if rank(b) < rank(a) { b } else { a }
}

/// What each answer says of the cases asked: at how many values of each param with cases it was shown, and
/// which combinations it was not asked at. One not loaded leaves a yes not settled.
fn said_of_cases(nodes: &mut [Node], cases: &CasesAsked) {
    let mut values: Vec<(String, Vec<&J>)> = vec![];
    for (param, value) in cases.asked.iter().filter_map(J::as_object).flatten() {
        match values.iter_mut().find(|(p, _)| p == param) {
            Some((_, seen)) => seen.extend((!seen.contains(&value)).then_some(value)),
            None => values.push((param.clone(), vec![value])),
        }
    }
    for node in nodes {
        for (param, seen) in &values {
            node.fact.covers.retain(|(variable, _)| variable != param);
            node.fact.covers.push((param.clone(), Covers::These(seen.len() + 1)));
        }
        node.fact.left_out.extend(cases.not_allowed.iter().map(|c| format!("the case {}: not asked there, {}", c["case"], c["why"].as_str().unwrap_or_default())));
        if !cases.not_loaded.is_empty() {
            node.fact.left_out.extend(cases.not_loaded.iter().map(|c| format!("the case {c}: not loaded before the call's time limit passed")));
            if node.fact.held == Some(Verdict::Yes) {
                node.fact.held = Some(Verdict::NotSettled(ink_claims::Unsettled::Undecided));
            }
        }
    }
}

/// Every requirement of ratings of the design `view` (only the one named `only`), answered by the core.
fn answer(view: Session, only: Option<&str>) -> Result<(Vec<Node>, Vec<(String, String)>), InkError> {
    let mut view = view;
    let (design, formulas) = ink_ratings::design(&mut view).map_err(InkError::msg)?;
    // A requirement of ratings names an assembly (the first of its path), whose motion and studies answer it.
    let within: Vec<(About, String)> =
        view.doc.requires.iter().filter(|r| only.is_none_or(|n| n == r.name)).filter_map(|r| Some((about(&r.name, ""), r.within.clone()?))).collect();
    // What each leaves out on purpose, and the bodies placed under its path.
    let excluding: Vec<Vec<String>> = view.doc.requires.iter().filter(|r| only.is_none_or(|n| n == r.name) && r.within.is_some()).map(|r| r.excluding.clone()).collect();
    let mut placed: Vec<Vec<String>> = vec![];
    // And those whose record states no limit, each with its words: `(attr no-rating "...")` on the body.
    let mut no_rating: Vec<(String, String)> = vec![];
    for (_, path) in &within {
        let items = view.doc.items(path).map_err(InkError::msg)?;
        for item in &items {
            if let Some((_, Value::Str(why))) = item.attrs.iter().find(|(name, _)| name == "no-rating") {
                no_rating.push((item.label.clone(), why.clone()));
            }
        }
        placed.push(items.into_iter().map(|item| item.label).collect());
    }
    // An exclusion does what it says or is an error: one that names nothing placed leaves nothing out.
    for (k, left) in excluding.iter().enumerate().flat_map(|(k, those)| those.iter().map(move |left| (k, left))) {
        if !placed[k].iter().any(|body| body == left || body.strip_prefix(left.as_str()).is_some_and(|rest| rest.starts_with('/'))) {
            return Err(InkError::msg(format!("(excluding {left}) of '{}' names nothing placed: no part is at or under {left}, so it leaves nothing out", within[k].0.what)));
        }
    }
    // The whole design, whatever the requirements name: every body the assemblies they name place that carries a
    // rating, and the paths every requirement of the file asks over (not only the one asked for here).
    let asked_over: Vec<String> = view.doc.requires.iter().filter_map(|r| r.within.clone()).collect();
    let left_out: Vec<String> = view.doc.requires.iter().filter(|r| r.within.is_some()).flat_map(|r| r.excluding.clone()).collect();
    let mut carrying: Vec<String> = vec![];
    {
        let mut tops: Vec<&str> = vec![];
        for path in &asked_over {
            let top = path.split('/').next().unwrap_or_default();
            if !tops.contains(&top) {
                tops.push(top);
                for item in view.doc.items(top).map_err(InkError::msg)? {
                    let rated = item.attrs.iter().any(|(name, value)| name == "ratings" && matches!(value, Value::Data(d) if d.value().as_object().is_some_and(|o| !o.is_empty())));
                    if rated && !carrying.contains(&item.label) {
                        carrying.push(item.label);
                    }
                }
            }
        }
    }
    let asked: Vec<(About, String)> = within.iter().map(|(a, path)| (a.clone(), path.split('/').next().unwrap_or_default().to_string())).collect();
    // A fact stated where nothing takes its name (no rating or formula of the part there reads it) is said under
    // each requirement over its place: a misspelled fact changes no answer, and must not pass unseen.
    let mut unread: Vec<(String, String)> = vec![];
    let mut takes = vec![];
    {
        for (_, path) in &within {
            takes.extend(ink_ratings::takes(&mut view.doc, path).map_err(InkError::msg)?);
        }
        for (fact, names) in ink_ratings::unread(&design, &takes) {
            unread.push((fact.about.on.clone(), ink_ratings::unread_said(fact, &names)));
        }
    }
    // The facts only the design can state (a required life), with their units: an answer that lacks one says the line.
    let mut must: Vec<(About, String)> = vec![];
    for (_, path) in &within {
        must.extend(ink_ratings::design_states(&mut view.doc, path));
    }
    let view = RefCell::new(view);
    // What a part works out from other facts about it (its `derives`), made before the tools that read the design.
    let derived = ink_ratings::Derived::new(&view).map_err(InkError::msg)?;
    // The solvers share the one design, each on an assembly the requirements name.
    let mut tops: Vec<&String> = vec![];
    asked.iter().for_each(|(_, top)| tops.extend((!tops.contains(&top)).then_some(top)));
    let solvers = match crate::solvers() {
        Some(plug) => tops.iter().map(|top| plug(&view, top)).collect::<Result<Vec<_>, _>>().map_err(InkError::msg)?,
        None => vec![],
    };
    let mut base: Vec<&dyn Tool> = vec![&formulas, &derived];
    // A motion's moments are asked here: after what is stated and tabled, before its peaks.
    let slot = base.len();
    solvers.iter().for_each(|s| base.extend(s.tools()));
    // What is of a whole duty and of no moment (a life's load and speed), over each assembly's motion; with
    // no solver, of what is stated alone.
    let duties: Vec<ink_ratings::duty::Duty> = match solvers.is_empty() {
        true => vec![ink_ratings::duty::Duty::new(&design, &base, slot, None, Box::new(|_| None), Box::new(Vec::new))],
        false => (tops.iter().zip(&solvers))
            .map(|(top, solver)| ink_ratings::duty::Duty::new(&design, &base, slot, Some(top.as_str()), Box::new(|moment| solver.at(moment)), Box::new(|| solver.steps())))
            .collect(),
    };
    let mut tools = base.clone();
    tools.extend(duties.iter().map(|d| d as &dyn Tool));
    let mut nodes = vec![];
    for (asked, top) in &asked {
        // A requirement's moments are those of the motion of the assembly it names.
        let motion = tops.iter().position(|t| *t == top).and_then(|i| solvers.get(i)).map(|s| &**s);
        let mut node = show(&Asking { design: &design, tools: &tools, slot, motion }, asked, 0);
        say_what_to_state(&mut node, &must);
        let path = within.iter().find(|(a, _)| a == asked).map_or("", |(_, path)| path.as_str());
        let under = |on: &str| on == path || on.strip_prefix(path).is_some_and(|rest| rest.starts_with(['/', ' ']));
        node.warnings = unread.iter().filter(|(on, _)| under(on)).map(|(_, said)| said.clone()).collect();
        // What it covered: the parts under its path that carry a rating, those the model loads that carry
        // none, its span, and the loads its model left out. Asked now, the run made (or never needed).
        let k = within.iter().position(|(a, _)| a == asked).expect("a requirement asked");
        let at_or_under = |on: &str, part: &str| on == part || on.strip_prefix(part).is_some_and(|rest| rest.starts_with(['/', ' ']));
        let left = |part: &str| excluding[k].iter().any(|left| at_or_under(part, left));
        let rated = |part: &str| !left(part) && design.rules.iter().any(|r| matches!(r.says, Says::Ordered(_)) && at_or_under(&r.about.on, part) && !left(&r.about.on));
        let coverage = motion.map(|m| m.coverage()).unwrap_or_default();
        let mut covered = Covered {
            span: coverage.span,
            goes_on: coverage.goes_on,
            loads_left: coverage.loads_left,
            unheld: coverage.unheld,
            assumed: coverage.assumed,
            unbalanced: coverage.unbalanced,
            balanced: coverage.balanced,
            balanced_moments: coverage.balanced_moments,
            ..Covered::default()
        };
        covered.rated = placed[k].iter().filter(|body| rated(body)).cloned().collect();
        for (part, why) in coverage.loaded.into_iter().filter(|(part, _)| under(part) && !rated(part)) {
            // Every body under it that says anything says its maker states no limit: nothing to hold it against.
            let none: Vec<&(String, String)> = no_rating.iter().filter(|(body, _)| at_or_under(body, &part)).collect();
            match left(&part) {
                true => {}
                false if !none.is_empty() => {
                    let said: Vec<&str> = none.iter().map(|(_, why)| why.as_str()).collect();
                    covered.states_none.push((part, format!("{why}; {}", dedup(said.iter().map(|s| s.to_string()).collect()).join("; "))));
                }
                false => covered.unrated.push((part, why)),
            }
        }
        // A rated part holding a coupling's pull, with no rating that reads a load across or along it.
        let reads_load = |part: &str| takes.iter().any(|a: &About| at_or_under(&a.on, part) && (a.what.ends_with("-load") || a.what.starts_with("moment-") || a.what == "tension"));
        for (part, why) in coverage.holds.into_iter().filter(|(part, _)| under(part) && rated(part) && !reads_load(part)) {
            covered.pull_unread.push((part, why));
        }
        // Left out on purpose: said as written, each path once.
        covered.excluded = excluding[k].clone();
        let named = |body: &String| covered.unrated.iter().chain(&covered.states_none).map(|(part, _)| part).chain(&covered.excluded).any(|part| at_or_under(body, part));
        covered.others = placed[k].iter().filter(|body| !rated(body) && !named(body)).count();
        // Books that do not balance: the loads are not all accounted for, and a yes resting on them is none.
        if !covered.unbalanced.is_empty() && node.fact.held == Some(Verdict::Yes) {
            (node.fact.held, node.fact.strength) = (Some(Verdict::NotSettled(ink_claims::Unsettled::Undecided)), None);
        }
        node.covered = Some(covered);
        nodes.push(node);
    }
    // What no requirement of the file asks over and none leaves out on purpose: a part carrying a rating, or
    // one the model loads. A yes to each requirement is then no yes about the design.
    let at_or_under = |on: &str, part: &str| on == part || on.strip_prefix(part).is_some_and(|rest| rest.starts_with(['/', ' ']));
    let accounted = |part: &str| asked_over.iter().chain(&left_out).any(|path| at_or_under(part, path));
    let mut unasked: Vec<(String, String)> = carrying.into_iter().filter(|part| !accounted(part)).map(|part| (part, "it carries ratings".to_string())).collect();
    for solver in &solvers {
        for (part, why) in solver.coverage().loaded {
            // (A loaded part that holds a rated body is that body's line.)
            if !accounted(&part) && !unasked.iter().any(|(p, _)| at_or_under(p, &part) || at_or_under(&part, p)) {
                unasked.push((part, why));
            }
        }
    }
    Ok((nodes, unasked))
}

/// Where a fact the design must state (`must`: each with its rating's unit) is stated by nothing, its node
/// under `node` says the line to write and where it goes, once for each place.
fn say_what_to_state(node: &mut Node, must: &[(About, String)]) {
    for under in &mut node.under {
        let Under::Fact(read) = under else { continue };
        say_what_to_state(read, must);
        if read.fact.held != Some(Verdict::NotAvailable(ink_claims::Unavailable::NoSource)) || !read.fact.notes.is_empty() {
            continue;
        }
        let Some((_, unit)) = must.iter().find(|(about, _)| *about == read.asked) else { continue };
        let (path, frame) = read.asked.on.split_once(' ').map_or((read.asked.on.as_str(), None), |(path, frame)| (path, Some(frame)));
        let name = path.rsplit('/').next().unwrap_or(path);
        let what = &read.asked.what;
        let fact = frame.map_or(format!("(fact {what} 10000{unit})"), |frame| format!("(fact {what} {frame} 10000{unit})"));
        read.fact.notes.push(format!(
            "only the design can say it: write it in the instance that places {path}, (instance {name} ... {fact}), with your own number in place of 10000 ({unit})"
        ));
    }
}

/// The prover's answer to a claim as a requirement's: the fact it is of the core (`ink_ratings::fact`), and
/// its evidence as the prover gives it.
fn proven(answer: &ink::claims::AssertResult, params: &[String]) -> Node {
    let fact = ink_ratings::fact(params, answer);
    Node {
        asked: fact.about.clone(),
        lead: "",
        terms: fact.terms.clone(),
        fact,
        kind: Kind::Used,
        under: vec![],
        proof: Some(assert_json(answer, params)),
        warnings: vec![],
        covered: None,
    }
}

/// The fact about `asked` over every moment of the motion (`ink_claims::
/// over`); then what it rests on. A requirement's ratings are each
/// answered so, with its own worst moment. A comparison's terms are those
/// of the moment it is the answer of, read from the one round of asking
/// that answered it there: what is under an answer is what decided it, a
/// fact a later tool settled as that tool gave it.
fn show(asking: &Asking, asked: &About, depth: usize) -> Node {
    let at = |moment: usize| asking.motion.and_then(|m| m.at(moment));
    let (fact, moment) = over(asking.design, asked, "t", asking.tools, asking.slot, &at);
    if let Some(Says::Within(_) | Says::WithinBut(..)) = asking.design.rules.iter().find(|r| r.about == *asked).map(|r| &r.says) {
        let under = fact.used.iter().map(|used| Under::Fact(show(asking, used, depth + 1))).collect();
        return Node { asked: asked.clone(), lead: "", fact, terms: vec![], kind: Kind::Ratings, under, proof: None, warnings: vec![], covered: None };
    }
    let tool = at(moment);
    let with: Vec<&dyn Tool> = asking.tools[..asking.slot].iter().copied().chain(tool.as_deref()).chain(asking.tools[asking.slot..].iter().copied()).collect();
    let mut claims = Claims::new(asking.design, &with);
    claims.answer(asked);
    let (kind, terms, under) = under(&mut claims, asking.design, asked, depth);
    Node { asked: asked.clone(), lead: "", fact, terms, kind, under, proof: None, warnings: vec![], covered: None }
}

/// What the fact about `asked` rests on, as `claims` found it: a
/// comparison's terms in order, anything else what it used, each with what
/// it rests on in turn.
fn under(claims: &mut Claims, design: &Design, asked: &About, depth: usize) -> (Kind, Vec<(String, f64)>, Vec<Under>) {
    if depth > 8 {
        return (Kind::Used, vec![], vec![]);
    }
    let fact = claims.answer(asked);
    let shown = |claims: &mut Claims, read: &About, lead: &'static str| {
        let fact = claims.answer(read);
        let (kind, terms, under) = under(claims, design, read, depth + 1);
        Under::Fact(Node { asked: read.clone(), lead, fact, terms, kind, under, proof: None, warnings: vec![], covered: None })
    };
    match design.rules.iter().find(|r| r.about == *asked).map(|r| &r.says) {
        Some(Says::Ordered(terms)) => {
            let each = terms.iter().map(|term| match term {
                Term::Number(n, unit) => Under::Number(*n, unit.clone()),
                Term::Read(read) => shown(claims, read, ""),
                Term::Size(read) => shown(claims, read, "the size of "),
            });
            (Kind::Ordered, vec![], each.collect())
        }
        Some(Says::Within(_) | Says::WithinBut(..)) => (Kind::Ratings, vec![], fact.used.iter().map(|used| shown(claims, used, "")).collect()),
        _ => (Kind::Used, fact.terms.clone(), fact.used.iter().map(|used| shown(claims, used, "")).collect()),
    }
}

fn number(x: f64, unit: &str) -> String {
    format!("{} {unit}", (x * 1000.0).round() / 1000.0).trim_end().to_string()
}

/// A verdict's word, without its reason.
fn verdict(v: Verdict) -> &'static str {
    match v {
        Verdict::Yes => "yes",
        Verdict::No => "no",
        Verdict::NotSettled(_) => "not settled",
        Verdict::NotAvailable(_) => "not available",
    }
}

fn strength(fact: &Fact) -> Option<&'static str> {
    match (fact.strength, fact.held) {
        (Some(Strength::Proven), _) => Some("proven"),
        (Some(Strength::ToAccuracy), _) => Some("proven to its accuracy"),
        (Some(Strength::Estimated), _) => Some("estimated"),
        // A yes or no that rests on what the design states, which nothing proved.
        (None, Some(Verdict::Yes | Verdict::No)) => Some("as stated"),
        (None, _) => None,
    }
}

/// Over how much of each variable it reads a yes or a no was shown.
fn covers(fact: &Fact) -> impl Iterator<Item = &(String, Covers)> {
    fact.covers.iter().filter(|_| fact.held.is_some_and(|v| !v.open()))
}

/// A node as its row of the answer.
fn row(node: &Node) -> J {
    let fact = &node.fact;
    let mut r = json!({"name": node.asked.what});
    if node.covered.is_some() {
        let said = match fact.held {
            Some(Verdict::Yes) => "YES".to_string(),
            Some(Verdict::No) => "NO".to_string(),
            Some(open) => open.to_string(),
            None => String::new(),
        };
        r["headline"] = json!(format!("{said}{}", headline(node)));
    }
    if !node.asked.on.is_empty() {
        r["on"] = json!(node.asked.on);
    }
    if !node.lead.is_empty() {
        r["size"] = json!(true);
    }
    match (fact.held, fact.value) {
        (Some(v), _) => {
            r["verdict"] = json!(verdict(v));
            if let Some((_, why)) = v.to_string().split_once(": ") {
                r["why"] = json!(why);
            }
        }
        (None, Some(value)) => {
            r["value"] = n(value);
            r["unit"] = json!(fact.unit);
        }
        (None, None) => {}
    }
    if let Some((low, high)) = fact.band {
        r["band"] = json!([n(low), n(high)]);
    }
    if let Some(s) = strength(fact) {
        r["strength"] = json!(s);
    }
    let covered: serde_json::Map<String, J> = covers(fact)
        .map(|(variable, covers)| {
            let how = match covers {
                Covers::One => json!("one value"),
                Covers::These(values) => json!(values),
                Covers::Every => json!("every value"),
            };
            (variable.clone(), how)
        })
        .collect();
    if !covered.is_empty() {
        r["covers"] = J::Object(covered);
    }
    if !fact.from.is_empty() {
        r["from"] = json!(fact.from);
    }
    for (key, lines) in [("warnings", &node.warnings), ("notes", &fact.notes), ("left_out", &fact.left_out)] {
        if !lines.is_empty() {
            r[key] = json!(lines);
        }
    }
    match &node.proof {
        // The prover's evidence, by its own keys (its `status` its own word for what it found, as `ink interfere`'s pairs have); its name and its line are the row's already.
        Some(J::Object(proof)) => {
            for (key, value) in proof.iter().filter(|(key, _)| !["name", "source"].contains(&key.as_str())) {
                r[key.as_str()] = value.clone();
            }
        }
        _ if !node.terms.is_empty() => r["terms"] = J::Object(node.terms.iter().map(|(k, x)| (k.clone(), n(*x))).collect()),
        _ => {}
    }
    if let Some(covered) = &node.covered {
        if let Some((a, b)) = covered.span {
            r["span"] = json!({"t": [n(a), n(b)]});
        }
        let unrated: Vec<J> = covered.unrated.iter().map(|(part, why)| json!({"part": part, "why": why})).collect();
        let none: Vec<J> = covered.states_none.iter().map(|(part, why)| json!({"part": part, "why": why})).collect();
        if !covered.pull_unread.is_empty() {
            r["loads_no_rating_reads"] = json!(covered.pull_unread.iter().map(|(part, why)| json!({"part": part, "why": why})).collect::<Vec<_>>());
        }
        r["parts"] = json!({"rated": covered.rated, "not_rated": unrated, "no_limit_stated": none, "left_out_on_purpose": covered.excluded, "others_not_rated": covered.others});
        if let Some(why) = &covered.goes_on {
            r["motion_goes_on"] = json!(why);
        }
        if !covered.unheld.is_empty() {
            r["not_driven_or_held"] = json!(covered.unheld);
        }
        if !covered.assumed.is_empty() {
            r["assuming"] = json!(covered.assumed);
        }
        if !covered.unbalanced.is_empty() {
            r["books_do_not_balance"] = json!(covered.unbalanced);
        }
        r["books_balanced"] = json!({"forces": covered.balanced, "moments": covered.balanced_moments});
        if !covered.loads_left.is_empty() {
            r["loads_not_applied"] = json!(covered.loads_left);
        }
    }
    match node.kind {
        Kind::Ordered => r["at_most"] = json!(true),
        Kind::Ratings => r["ratings"] = json!(node.under.len()),
        Kind::Used => {}
    }
    if !node.under.is_empty() {
        let each = node.under.iter().map(|u| match u {
            Under::Number(x, unit) => json!({"number": n(*x), "unit": unit}),
            Under::Fact(node) => row(node),
        });
        r["rests_on"] = J::Array(each.collect());
    }
    r
}

/// A node in words: its line and notes, then what it rests on, indented.
/// A node's lines in words, then what is under it. Each note and each thing left out is said once in an
/// answer, where it is first met (top down: the requirement's own line first); lines alike but for their
/// numbers (the same thing left out at one moment after another) are said by their first, with how many more.
/// A fact's own moment (`at t = ...`) is always said. The JSON rows carry every line.
fn words(text: &mut String, node: &Node, depth: usize, said: &mut std::collections::HashSet<String>) {
    let (pad, fact) = ("  ".repeat(depth), &node.fact);
    let mut line = format!("{pad}{}{}: ", node.lead, node.asked);
    line += &match (fact.held, fact.value) {
        (Some(Verdict::Yes), _) => "YES".to_string(),
        (Some(Verdict::No), _) => "NO".to_string(),
        (Some(open), _) => open.to_string(),
        (None, Some(value)) => number(value, &fact.unit),
        (None, None) => "no number".to_string(),
    };
    if let Some((low, high)) = fact.band {
        line += &format!(", between {} and {}", number(low, ""), number(high, &fact.unit));
    }
    line += &strength(fact).map(|s| format!(", {s}")).unwrap_or_default();
    for (variable, covers) in covers(fact) {
        line += &match covers {
            Covers::One => format!(", at one value of {variable}"),
            Covers::These(values) => format!(", at {values} values of {variable}"),
            Covers::Every => format!(", at every value of {variable}"),
        };
        // Over what span of it, where its motion's run says.
        if let Some((a, b)) = node.covered.as_ref().and_then(|c| c.span).filter(|_| variable == "t") {
            line += &format!(" from {} to {}", number(a, ""), number(b, ""));
        }
    }
    line += &headline(node);
    if !fact.from.is_empty() {
        line += &format!("  ({})", fact.from);
    }
    *text += &format!("{line}\n");
    // A fact shown in full once is, where it is read again, its line and its moment: what it rests on is above.
    let again = !matches!(node.kind, Kind::Ratings) && !node.under.is_empty() && !said.insert(format!("fact: {}", node.asked));
    let leaf_again = node.under.is_empty() && !said.insert(format!("fact: {}", node.asked));
    if again || leaf_again {
        fact.notes.iter().filter(|l| l.starts_with("at t = ")).for_each(|l| *text += &format!("{pad}  note: {l}\n"));
        if again {
            *text += &format!("{pad}  (what it rests on: as above)\n");
        }
        return;
    }
    // What it did not cover, first: a yes is of the rest.
    if let Some(covered) = &node.covered {
        let loaded = covered.unrated.len() + covered.excluded.len() + covered.states_none.len();
        if !covered.unrated.is_empty() {
            let each: Vec<String> = covered.unrated.iter().map(|(part, why)| format!("{part} ({why})")).collect();
            *text += &format!(
                "{pad}  NOT COVERED: {} of the parts the model loads under it carr{} no rating, so nothing here checks {}: {}\n",
                covered.unrated.len(),
                if covered.unrated.len() == 1 { "ies" } else { "y" },
                if covered.unrated.len() == 1 { "it" } else { "them" },
                each.join("; ")
            );
        }
        if !covered.states_none.is_empty() {
            let each: Vec<String> = covered.states_none.iter().map(|(part, why)| format!("{part} ({why})")).collect();
            *text += &format!(
                "{pad}  NOT COVERED, and nothing could: {} of the parts the model loads under it {} no limit stated by {} maker, so there is nothing to hold {} against: {}\n",
                covered.states_none.len(),
                if covered.states_none.len() == 1 { "has" } else { "have" },
                if covered.states_none.len() == 1 { "its" } else { "their" },
                if covered.states_none.len() == 1 { "it" } else { "them" },
                each.join("; ")
            );
        }
        if !covered.pull_unread.is_empty() {
            let each: Vec<String> = covered.pull_unread.iter().map(|(part, why)| format!("{part} ({why})")).collect();
            *text += &format!(
                "{pad}  NOT COVERED: a load on {} that none of its ratings reads, so nothing here checks it (its record states no limit for a load across or along its \
                 shaft or travel): {}\n",
                if covered.pull_unread.len() == 1 { "a rated part" } else { "rated parts" },
                each.join("; ")
            );
        }
        covered.unbalanced.iter().for_each(|l| *text += &format!("{pad}  THE BOOKS DO NOT BALANCE: {l}. The loads above are not all accounted for: this is no yes\n"));
        covered.assumed.iter().for_each(|l| *text += &format!("{pad}  ASSUMING {l}\n"));
        covered.unheld.iter().for_each(|l| *text += &format!("{pad}  NOT COVERED: {l}\n"));
        covered.loads_left.iter().for_each(|l| *text += &format!("{pad}  NOT COVERED: {l}\n"));
        covered.goes_on.iter().for_each(|l| *text += &format!("{pad}  NOT COVERED: {l}\n"));
        if !covered.excluded.is_empty() {
            *text += &format!("{pad}  left out on purpose (excluding): {}\n", covered.excluded.join(", "));
        }
        *text += &format!(
            "{pad}  covers: {} part{} carrying ratings ({}){}{}\n",
            covered.rated.len(),
            if covered.rated.len() == 1 { "" } else { "s" },
            if covered.rated.is_empty() { "none".to_string() } else { covered.rated.join(", ") },
            if loaded > 0 { format!("; of the {} it loads without one, {} left out on purpose", loaded, covered.excluded.len()) } else { String::new() },
            if covered.others > 0 {
                format!("; {} other bod{} under it carry none, and the model puts no demand a maker rates on them", covered.others, if covered.others == 1 { "y" } else { "ies" })
            } else {
                String::new()
            },
        );
    }
    node.warnings.iter().for_each(|l| *text += &format!("{pad}  warning: {l}\n"));
    // What the prover says beside its answer (what a study was run with, what a claim left out on purpose).
    for l in node.proof.as_ref().and_then(|p| p["warnings"].as_array()).into_iter().flatten().filter_map(J::as_str) {
        *text += &format!("{pad}  note: {l}\n");
    }
    // A claim with no value: the band it was decided on, and why, said on its own row whatever rows above said.
    let proof = node.proof.as_ref().filter(|p| p["value"].is_null());
    let why = proof.and_then(|p| p["unavailable"].as_str());
    if let Some(p) = proof
        && let (Some(lo), Some(hi)) = (p["band"]["low"].as_f64(), p["band"]["high"].as_f64())
    {
        let (lo, hi) = (number(lo, ""), number(hi, ""));
        *text += &if lo == hi { format!("{pad}  with: value {lo}\n") } else { format!("{pad}  with: value between {lo} and {hi}\n") };
    }
    if let Some(why) = why.filter(|w| !fact.notes.iter().any(|l| l == w)) {
        *text += &format!("{pad}  note: {why}\n");
    }
    for (what, lines) in [("note", &fact.notes), ("left out", &fact.left_out)] {
        // Alike but for their numbers: one line, and a count.
        let alike = |l: &str| l.chars().map(|c| if c.is_ascii_digit() { '#' } else { c }).collect::<String>();
        let mut groups: Vec<(String, &String, usize)> = vec![];
        for l in lines.iter().filter(|l| l.starts_with("at t = ") || Some(l.as_str()) == why || said.insert(format!("{what}: {l}"))) {
            match groups.iter_mut().find(|(key, ..)| *key == alike(l)) {
                Some((_, _, more)) => *more += 1,
                None => groups.push((alike(l), l, 0)),
            }
        }
        for (_, l, more) in groups {
            let more = if more > 0 { format!("  (and {more} more like it, with other numbers)") } else { String::new() };
            *text += &format!("{pad}  {what}: {l}{more}\n");
        }
    }
    match node.kind {
        // A requirement says how many ratings its parts carry, then each: what no line names was not checked.
        Kind::Ratings => *text += &format!("{pad}  the {} rating{} its parts carry, each checked:\n", node.under.len(), if node.under.len() == 1 { "" } else { "s" }),
        Kind::Used if !node.terms.is_empty() => {
            *text += &format!("{pad}  with: {}\n", node.terms.iter().map(|(name, x)| format!("{name} {}", number(*x, ""))).collect::<Vec<_>>().join(", "))
        }
        _ => {}
    }
    for (i, under) in node.under.iter().enumerate() {
        if i > 0 && matches!(node.kind, Kind::Ordered) {
            *text += &format!("{pad}  at most\n");
        }
        match under {
            Under::Number(x, unit) => *text += &format!("{pad}  {}\n", number(*x, unit)),
            Under::Fact(node) => words(text, node, depth + 1, said),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink_claims::{Rule, Unsettled};

    /// A tool that gives `load` as a number, or between two.
    struct Gives(&'static str, f64, Option<(f64, f64)>);

    impl Tool for Gives {
        fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
            (asked.what == "load").then(|| Fact { band: self.2, strength: Some(Strength::Estimated), ..Fact::stated(asked.clone(), self.1, "N", self.0) })
        }
    }

    /// What is printed under an answer is the source that decided it: a first tool's band across the limit
    /// goes on to a second, whose number settles it, and it is the second's line that is shown (an outside
    /// read found the first's, a band across the limit, printed under the yes). The row says what the words say.
    #[test]
    fn the_source_printed_is_the_one_that_decided() {
        let ok = about("ok", "");
        let rule = Rule { about: ok.clone(), says: Says::Ordered(vec![Term::Read(about("load", "")), Term::Number(5.0, "N".into())]), from: "a limit".into() };
        let design = Design { facts: vec![], rules: vec![rule] };
        let (coarse, fine) = (Gives("a coarse tool", 5.0, Some((4.0, 6.0))), Gives("a finer tool", 4.5, None));
        let said = |tools: &[&dyn Tool]| {
            let node = show(&Asking { design: &design, tools, slot: 0, motion: None }, &ok, 0);
            let mut text = String::new();
            words(&mut text, &node, 0, &mut Default::default());
            (node.fact.held, text, row(&node))
        };
        let (held, text, row) = said(&[&coarse, &fine]);
        assert_eq!((held, text.as_str()), (Some(Verdict::Yes), "ok: YES, estimated  (a limit)\n  load: 4.5 N, estimated  (a finer tool)\n  at most\n  5 N\n"));
        let rests = json!([{"name": "load", "value": 4.5, "unit": "N", "strength": "estimated", "from": "a finer tool"}, {"number": 5.0, "unit": "N"}]);
        assert_eq!(row, json!({"name": "ok", "verdict": "yes", "strength": "estimated", "from": "a limit", "at_most": true, "rests_on": rests}));
        // Alone, the coarse tool's band is what there is, and it is not settled.
        let (held, text, row) = said(&[&coarse]);
        assert!(held == Some(Verdict::NotSettled(Unsettled::Straddles)) && text.contains("load: 5 N, between 4 and 6 N, estimated  (a coarse tool)"), "{text}");
        assert!(row["verdict"] == "not settled" && row["why"] == "the band reaches across the limit" && row["rests_on"][0]["band"] == json!([4.0, 6.0]), "{row}");
    }
}
