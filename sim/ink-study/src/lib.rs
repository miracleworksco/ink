//! The study plug: what a part bears under a `(study ...)`'s loads, as
//! facts for the claims core.
//!
//! The solver is ink-fem, unchanged, entered as `ink solve` enters it
//! (ink-sim's `SolvedStudy::run`, which keeps a solve by the content of its
//! question, in memory and on disk). On each instance under the target
//! whose body's shape is the study's, at its path P, it offers `max-stress
//! on P` (von Mises, MPa) and `max-displacement on P` (mm): the language's
//! units, which are the solver's own, so a part's rating (ink-ratings) can
//! be held against them.
//!
//! One grid shows no error of its own, so no number is given off one: a
//! value is the one on the study's own cells checked on cells twice and
//! four times as coarse (or finer, where those cannot see the part), by the
//! comparison `ink require` makes (`ink_lang::claims::ladder`). Its band is
//! the grid error that shows, and beyond it what the solver has measured
//! its read may miss (a peak stress's, fitted on reads with known answers).
//! The first solve is started the first time a fact is asked for, and each
//! grid's is kept.
//!
//! Every fact is an estimate. Where the solver says its number is not the
//! part's (a solve that failed, a feature its cells could not see, loads
//! that may buckle it, grids that do not show how it moves) the fact is
//! open and says why, and carries no number. A solve is at one value of
//! each param: where the study reads a ranged one no fact is given and
//! nothing is solved, as the old claim on a study says.
//!
//! Several studies of one part are its load cases, and a part's stress
//! or deflection is its worst under any of them: `Studies` is every static
//! elastic study of a target's parts as one tool, which gives the largest,
//! and nothing where one of them is open. A study alone (`Studied`) says
//! which others it leaves out.

use ink_claims::{About, Claims, Fact, Strength, Tool, Unavailable, Unsettled, Verdict};
use ink_lang::view::ModelView;
use ink_lang::{FormMetric, Outcome, RunOpts, claims, command};
use ink_sim::{SolvedStudy, Study, StudyPhysics};
use std::cell::{Cell, OnceCell, RefCell};
use std::collections::HashMap;

/// What it offers, each the solver's metric of that name, and its unit.
const OFFERS: [(&str, &str); 2] = [("max-stress", "MPa"), ("max-displacement", "mm")];

pub struct Studied<'a> {
    view: &'a RefCell<ModelView>,
    study: String,
    /// The instances under the target whose shape the study solves, by path.
    places: Vec<String>,
    /// The other studies of the same shape: load cases its facts leave out.
    left_out: Vec<String>,
    /// The ranged params the study reads, each with its range: a solve is at one value of each.
    ranged: Vec<String>,
    run: OnceCell<Result<SolvedStudy, String>>,
    runs: Cell<usize>,
    /// The fact found for each thing it offers, less its place.
    found: RefCell<HashMap<&'static str, Fact>>,
}

impl<'a> Studied<'a> {
    /// The plug for study `study` (the file's only one where None) and the
    /// parts under the assembly or instance path `target` of a loaded
    /// design, which it shares with the other tools asked beside it. It
    /// reads which parts the study is of, and solves nothing.
    pub fn new(shared: &'a RefCell<ModelView>, study: Option<&str>, target: &str) -> Result<Studied<'a>, String> {
        let view = &mut *shared.borrow_mut();
        let name = view.doc.picked::<Study>("study", study).map_err(|e| e.message)?.name.clone();
        let (shape, left_out) = {
            let st = view.doc.form::<Study>(&name).expect("picked");
            if !matches!(st.physics, StudyPhysics::Elastic { .. }) || st.modes.is_some() || st.transient.is_some() {
                return Err(format!("study '{name}' is not a static elastic study: it has no peak stress or deflection under loads"));
            }
            let others: Vec<String> = view.doc.forms::<Study>().filter(|o| o.shape == st.shape && o.name != st.name).map(|o| format!("'{}'", o.name)).collect();
            (st.shape, if others.is_empty() { vec![] } else { vec![format!("the other studies of this shape: {}", others.join(", "))] })
        };
        let places: Vec<String> = view.doc.items(target)?.into_iter().filter(|item| item.shape == shape).map(|item| item.label).collect();
        if places.is_empty() {
            return Err(format!("study '{name}' is of no part placed under {target}"));
        }
        // A ranged param the study reads is solved at its one value.
        let ranged = view.doc.params_read(&ink_sim::study_inputs(&view.doc, &name)).into_iter().filter(|&k| view.ctx.p_spans.get(k as u32).is_some());
        let ranged = ranged.map(|k| format!("'{}'", view.doc.params[k].name)).collect();
        Ok(Studied { view: shared, study: name, places, left_out, ranged, run: OnceCell::new(), runs: Cell::new(0), found: RefCell::default() })
    }

    /// How many times a solve on the study's own cells has been started.
    pub fn runs(&self) -> usize {
        self.runs.get()
    }

    fn run(&self) -> &Result<SolvedStudy, String> {
        self.run.get_or_init(|| {
            self.runs.set(self.runs.get() + 1);
            let view = &mut *self.view.borrow_mut();
            // As `ink solve` does: what the study rests on is provided first.
            command::solved(view).map_err(|e| e.message)?;
            let inputs = view.doc.form_inputs(view.doc.ext_form(Some("sim"), &self.study).expect("picked"));
            view.want(&inputs).map_err(|e| e.message)?;
            // What it reads of the model, at the motion's worst moment.
            let ctx = match view.at_worst(Some("sim"), &self.study).map_err(|e| e.message)? {
                // Over a motion its loads vary, and one moment's solve is no answer for the part: the study's
                // own claim bounds every moment (`ModelView::measure`), and this fact gives none.
                Some(_) if view.motion.or(view.doc.t_span).is_some() => {
                    return Err(format!(
                        "study '{0}' is loaded from the model over a motion, and this fact would be of one moment of it. Ask the study itself, which answers for every \
                         moment: (require NAME (max-stress {0} < LIMIT)), (require NAME (max-displacement {0} < LIMIT))",
                        self.study
                    ));
                }
                Some((ctx, _)) => ctx,
                None => view.ctx.clone(),
            };
            SolvedStudy::run(&mut view.doc, &ctx, &self.study, &RunOpts::default()).map_err(|e| e.message)
        })
    }

    /// The fact about `what` (one of OFFERS, in `unit`), less its place.
    fn fact(&self, what: &'static str, unit: &str) -> Fact {
        let mut fact = Fact { unit: unit.into(), from: format!("ink solve --study {}", self.study), ..Fact::default() };
        let open = |mut fact: Fact, why, note: String| {
            (fact.held, fact.notes) = (Some(why), vec![note]);
            fact
        };
        let (outside, undecided) = (Verdict::NotAvailable(Unavailable::Outside), Verdict::NotSettled(Unsettled::Undecided));
        let run = match self.run() {
            Ok(run) => run,
            Err(why) => return open(fact, outside, why.clone()),
        };
        let Some(mut read) = run.metric(what) else { return open(fact, outside, format!("the solve has no {what}")) };
        // The value on the finest of three grids that read the part, and the grid error they show.
        let mut grid = (0.0, String::new());
        if read.value.is_ok() && read.unseen.is_none() {
            let metric = FormMetric { package: "sim", form: self.study.clone(), metric: what, shape: None };
            match claims::ladder(&mut self.view.borrow_mut(), &metric, None, &read) {
                Ok(ladder) => {
                    let ((error, how), [_, x2, x4], [_, h2, h4]) = (ladder.error(), ladder.values, &ladder.named);
                    grid =
                        (error, format!("its grid error about {} {unit} ({how}: {} on cells of {h2}, {} on cells of {h4})", claims::sig(error), claims::sig(x2), claims::sig(x4)));
                    read = ladder.finest;
                }
                Err((h, why)) => {
                    return open(
                        fact,
                        undecided,
                        format!("no three grids, each twice as coarse as the last, the coarsest no coarser than four times its {h} mm, all read the part: {why}"),
                    );
                }
            }
        }
        let value = match &read.value {
            Ok(value) => *value,
            Err(why) => return open(fact, undecided, why.clone()),
        };
        // The number is another part's where the cells could not see this one.
        if let Some(why) = &read.unseen {
            return open(fact, undecided, format!("its cells could not see the part: {why}"));
        }
        // And no number at all where the linear solve may not be the part's response.
        if let Some(e) = read.excluded.iter().find(|e| e.possible) {
            return open(fact, undecided, format!("{}; {}", e.what, e.decide));
        }
        let cells = match (read.res, read.finest) {
            (Some(h), Some(f)) => format!("cells of {h} mm, split to {f} mm where a feature is thin"),
            (h, _) => format!("cells of {} mm", h.unwrap_or(run.sol.h)),
        };
        // What its read may miss, as the solver has measured it: beyond the grid error, never in its place.
        let fitted = read.read_error.filter(|r| !r.proven);
        let (low, high) = claims::widened(value, grid.0);
        let (under, over) = fitted.map_or((value, value), |r| r.span(value));
        (fact.value, fact.band, fact.strength) = (Some(value), Some((low - (value - under), high + (over - value))), Some(Strength::Estimated));
        fact.left_out.clone_from(&self.left_out);
        fact.left_out.extend(read.excluded.iter().map(|e| format!("{}: {} {unit} with it", e.what, claims::sig(e.value))));
        fact.notes.push(format!("a linear static finite-cell solve on {cells}"));
        fact.notes.push(grid.1);
        fact.notes.extend(fitted.map(|r| r.how.to_string()));
        fact.notes.extend(read.warnings);
        fact
    }
}

impl Tool for Studied<'_> {
    /// A stress or a deflection is a demand: a typed one is held beside the study's.
    fn measures(&self) -> Option<ink_claims::Worse> {
        Some(ink_claims::Worse::More)
    }

    fn ask(&self, about: &About, _: &mut Claims) -> Option<Fact> {
        let &(what, unit) = OFFERS.iter().find(|(what, _)| *what == about.what)?;
        if !self.places.contains(&about.on) {
            return None;
        }
        if !self.ranged.is_empty() {
            let why = format!(
                "a study is solved at one value of each param, and this one reads ranged {}: that value says nothing of the rest of the range. Ask at one value (--param NAME=V)",
                self.ranged.join(", ")
            );
            return Some(Fact {
                about: about.clone(),
                held: Some(Verdict::NotAvailable(Unavailable::Outside)),
                notes: vec![why],
                from: format!("ink solve --study {}", self.study),
                ..Fact::default()
            });
        }
        let kept = self.found.borrow().get(what).cloned();
        let fact = kept.unwrap_or_else(|| {
            let fact = self.fact(what, unit);
            self.found.borrow_mut().insert(what, fact.clone());
            fact
        });
        Some(Fact { about: about.clone(), ..fact })
    }
}

/// Every static elastic study of the parts under a target, as one tool: a
/// part's load cases.
pub struct Studies<'a>(Vec<Studied<'a>>);

impl<'a> Studies<'a> {
    /// The tool for `target` of a loaded design: each of its studies that
    /// is a static elastic one of a part placed there. It solves nothing.
    pub fn new(view: &'a RefCell<ModelView>, target: &str) -> Studies<'a> {
        let names: Vec<String> = view.borrow().doc.forms::<Study>().map(|s| s.name.clone()).collect();
        let mut studies: Vec<Studied<'a>> = names.iter().filter_map(|name| Studied::new(view, Some(name), target).ok()).collect();
        // Together they leave no study of a shape out.
        studies.iter_mut().for_each(|s| s.left_out.clear());
        Studies(studies)
    }

    /// How many solves on a study's own cells have been started.
    pub fn runs(&self) -> usize {
        self.0.iter().map(Studied::runs).sum()
    }
}

impl Tool for Studies<'_> {
    /// A stress or a deflection is a demand: a typed one is held beside the study's.
    fn measures(&self) -> Option<ink_claims::Worse> {
        Some(ink_claims::Worse::More)
    }

    /// The largest of the studies' facts about it, by the high end of its
    /// band; the first open one where one is: a load case not known is a
    /// part not known.
    fn ask(&self, about: &About, claims: &mut Claims) -> Option<Fact> {
        let mut facts: Vec<Fact> = self.0.iter().filter_map(|study| study.ask(about, claims)).collect();
        if let Some(open) = facts.iter().position(Fact::open) {
            return Some(facts.swap_remove(open));
        }
        let high = |f: &Fact| f.range().map_or(f64::NEG_INFINITY, |r| r.1);
        let worst = (0..facts.len()).max_by(|&a, &b| high(&facts[a]).total_cmp(&high(&facts[b])))?;
        let mut fact = facts.swap_remove(worst);
        if !facts.is_empty() {
            fact.notes.push(format!("the worst of its {} studies; the others: {}", facts.len() + 1, facts.iter().map(|f| f.from.as_str()).collect::<Vec<_>>().join(", ")));
            for other in facts {
                fact.left_out.extend(other.left_out.into_iter().filter(|l| !fact.left_out.contains(l)).collect::<Vec<_>>());
            }
        }
        Some(fact)
    }
}
