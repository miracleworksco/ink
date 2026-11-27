//! A duty: what a part bears over a whole motion, where the fact asked is
//! of the duty and of no one moment (a bearing's life: ISO 281's
//! equivalent load and mean speed over the cycle).
//!
//! A rating read over a motion is answered at every moment and is the
//! worst one's (`ink_claims::over`). That is right for a limit (a load, a
//! speed, a torque under a curve) and wrong for a life, which no moment
//! has: read per moment it is the life as if the part ran for ever at its
//! worst instant, and at a moment of rest it is no life at all. So a part's
//! formula for a life reads two facts this tool offers at the part's place:
//!
//! - `mean-speed`: the mean over the motion's time of the size of `speed`
//!   there;
//! - `duty-load`: the mean of `equivalent-load` there to the power p,
//!   weighed by the turns made at each moment (the size of `speed` times
//!   its time), and its p-th root: P = (sum P_i^p n_i dt_i / sum n_i
//!   dt_i)^(1/p), ISO 281's equivalent load of a varying duty; p the fact
//!   `life-exponent` there (3 for balls, 10/3 for rollers: the part's);
//! - `rms-torque`: the root of the mean over the motion's time of `torque`
//!   there squared: what a maker's continuous rating is of (a gearbox's
//!   recommended continuous torque, a motor's heating), where its
//!   momentary one is of each moment.
//!
//! Each is read at every moment the motion's run read, each moment weighed
//! by the time of the step it is the middle of (the midpoint rule: the
//! span's two ends weigh nothing). They are estimates, with the band the
//! moments' own bands give, and say how many moments they are of.
//!
//! With no motion there is no duty to work out: `mean-speed` is the speed
//! the design states there, and where it states none it is not available,
//! saying so (a still pose has no life, and the speed a model happens to
//! have at one moment is not a duty); `duty-load` is the load there as it
//! is.

use ink_claims::{About, Claims, Covers, Design, Fact, Strength, Tool, Unavailable, Verdict, about};
use std::cell::RefCell;
use std::collections::HashMap;

/// What a duty fact is the duty of: the fact read at each moment.
pub const READS: [(&str, &str); 3] = [("mean-speed", "speed"), ("duty-load", "equivalent-load"), ("rms-torque", "torque")];

type At<'a> = Box<dyn Fn(usize) -> Option<Box<dyn Tool + 'a>> + 'a>;

pub struct Duty<'a> {
    design: &'a Design,
    /// The tools a moment's facts are asked of, a moment's own asked before the one at `slot` (as `over` asks).
    tools: &'a [&'a dyn Tool],
    slot: usize,
    /// The places it answers at: those at or under this path (None: any).
    top: Option<String>,
    /// The tool of the motion's moment k, and the time each moment weighs (none: no motion).
    at: At<'a>,
    steps: Box<dyn Fn() -> Vec<f64> + 'a>,
    kept: RefCell<HashMap<About, Fact>>,
}

impl<'a> Duty<'a> {
    pub fn new(design: &'a Design, tools: &'a [&'a dyn Tool], slot: usize, top: Option<&str>, at: At<'a>, steps: Box<dyn Fn() -> Vec<f64> + 'a>) -> Duty<'a> {
        Duty { design, tools, slot, top: top.map(String::from), at, steps, kept: RefCell::default() }
    }

    /// With no motion: the speed stated, or none; the load as it is.
    fn still(&self, asked: &About, of: &About) -> Fact {
        // The speed: what is stated, or what the part works out from what is stated. The load: as the design has
        // it at its one pose, the solvers' reading of that pose among its sources.
        let pose = (self.at)(0);
        let with: Vec<&dyn Tool> = self.tools[..self.slot].iter().copied().chain(pose.as_deref()).chain(self.tools[self.slot..].iter().copied()).collect();
        let mut read = Claims::new(self.design, if asked.what == "mean-speed" { &self.tools[..self.slot] } else { &with }).answer(of);
        if asked.what == "mean-speed" && read.held == Some(Verdict::NotAvailable(Unavailable::NoSource)) {
            read.held = Some(Verdict::NotAvailable(Unavailable::Outside));
            read.notes = vec![format!(
                "a life is of a duty: it needs a motion (a (t-span) the design turns over) or a stated speed ((fact speed ...) on {}). At one still pose there is none",
                of.on
            )];
        }
        Fact { about: asked.clone(), used: vec![of.clone()], ..read }
    }

    /// Over the motion's moments, each weighing `steps[k]` of time.
    fn moving(&self, asked: &About, of: &About, steps: &[f64]) -> Fact {
        let speed = about("speed", &asked.on);
        let open = |fact: &Fact, k: usize| Fact {
            about: asked.clone(),
            used: vec![of.clone()],
            notes: [vec![format!("at moment {k} of the motion")], fact.notes.clone()].concat(),
            ..fact.clone()
        };
        // Each moment: what is read there (its number and band), the turns it weighs, and what it left out.
        let (mut read, mut left_out, mut unit) = (vec![], vec![], String::new());
        for (k, step) in steps.iter().enumerate() {
            let Some(tool) = (self.at)(k) else { break };
            let with: Vec<&dyn Tool> = self.tools[..self.slot].iter().copied().chain([tool.as_ref()]).chain(self.tools[self.slot..].iter().copied()).collect();
            let mut claims = Claims::new(self.design, &with);
            // (A torque's duty reads no speed: its own stands in, weighing nothing.)
            let x = claims.answer(of);
            let n = if asked.what == "rms-torque" { x.clone() } else { claims.answer(&speed) };
            // Where neither is read off the motion (both stated, or worked out from what is stated), the duty
            // is what is stated: nothing is a mean, and nothing is estimated.
            if k == 0 && x.strength.is_none() && n.strength.is_none() {
                return self.still(asked, of);
            }
            for fact in [&x, &n] {
                if fact.open() || fact.value.is_none() {
                    return open(fact, k);
                }
                left_out.extend(fact.left_out.iter().filter(|l| !left_out.contains(*l)).cloned().collect::<Vec<_>>());
            }
            let sizes = |f: &Fact| {
                f.range().map(|(low, high)| (if low <= 0.0 && high >= 0.0 { 0.0 } else { low.abs().min(high.abs()) }, low.abs().max(high.abs()))).expect("it has a value")
            };
            unit = x.unit.clone();
            read.push((x.value.expect("checked").abs(), sizes(&x), n.value.expect("checked").abs(), *step));
        }
        let time: f64 = read.iter().map(|r| r.3).sum();
        let mut fact = Fact {
            about: asked.clone(),
            unit,
            strength: Some(Strength::Estimated),
            covers: vec![("t".into(), Covers::These(read.len()))],
            left_out,
            // It rests on every moment's facts, each weighed here (an open one is returned above): none is named
            // as what it rests on, since one moment's would be shown as if it were the duty's.
            used: vec![],
            ..Fact::default()
        };
        fact.from = "its duty over the motion".into();
        if asked.what == "mean-speed" {
            // The size of the speed, over time.
            let mean = |of: &dyn Fn(&(f64, (f64, f64), f64, f64)) -> f64| read.iter().map(|r| of(r) * r.3).sum::<f64>() / time;
            let (value, low, high) = (mean(&|r| r.0), mean(&|r| r.1.0), mean(&|r| r.1.1));
            (fact.value, fact.band) = (Some(value), (low < high).then_some((low, high)));
            fact.notes.push(format!("the mean of the size of its speed over the motion's time, read at {} moments", read.len()));
            return fact;
        }
        if asked.what == "rms-torque" {
            // The root of the mean square over time, of its size at each moment and of each end of its band.
            let rms = |of: &dyn Fn(&(f64, (f64, f64), f64, f64)) -> f64| (read.iter().map(|r| of(r).powi(2) * r.3).sum::<f64>() / time).sqrt();
            let (value, low, high) = (rms(&|r| r.0), rms(&|r| r.1.0), rms(&|r| r.1.1));
            (fact.value, fact.band) = (Some(value), (low < high).then_some((low, high)));
            fact.notes.push(format!("the root of the mean of its torque squared over the motion's time, read at {} moments", read.len()));
            return fact;
        }
        let exponent = Claims::new(self.design, &self.tools[..self.slot]).answer(&about("life-exponent", &asked.on));
        let Some(p) = exponent.value.filter(|p| !exponent.open() && *p > 0.0) else {
            return Fact {
                about: asked.clone(),
                used: vec![about("life-exponent", &asked.on)],
                notes: vec!["a duty's load is its loads' mean to its part's life exponent, which the part does not give".into()],
                ..exponent
            };
        };
        let turns: f64 = read.iter().map(|r| r.2 * r.3).sum();
        if turns <= 0.0 {
            // It does not turn: no load wears it. The greatest is given, which a life at rest does not read.
            fact.value = Some(read.iter().map(|r| r.1.1).fold(0.0, f64::max));
            fact.notes.push(format!("it does not turn over the motion: the greatest load of {} moments", read.len()));
            return fact;
        }
        let mean = |of: &dyn Fn(&(f64, (f64, f64), f64, f64)) -> f64| (read.iter().map(|r| of(r).powf(p) * r.2 * r.3).sum::<f64>() / turns).powf(1.0 / p);
        let (value, low, high) = (mean(&|r| r.0), mean(&|r| r.1.0), mean(&|r| r.1.1));
        (fact.value, fact.band) = (Some(value), (low < high).then_some((low, high)));
        fact.notes.push(format!("the mean of its load to the power {p} by the turns made at each of {} moments, and its root (ISO 281, a varying duty)", read.len()));
        fact
    }
}

impl Tool for Duty<'_> {
    fn ask(&self, asked: &About, _: &mut Claims) -> Option<Fact> {
        let (_, of) = READS.iter().find(|(what, _)| *what == asked.what)?;
        if self.top.as_ref().is_some_and(|top| asked.on != *top && !asked.on.strip_prefix(top.as_str()).is_some_and(|rest| rest.starts_with(['/', ' ']))) {
            return None;
        }
        if let Some(kept) = self.kept.borrow().get(asked) {
            return Some(kept.clone());
        }
        let of = about(of, &asked.on);
        let steps = (self.steps)();
        let fact = if steps.is_empty() { self.still(asked, &of) } else { self.moving(asked, &of, &steps) };
        self.kept.borrow_mut().insert(asked.clone(), fact.clone());
        Some(fact)
    }
}
