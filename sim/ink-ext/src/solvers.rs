//! The solver tools of this distribution for `ink require` (ink-cli's
//! `Solvers`): for the assembly a requirement names, its motion's efforts
//! (ink-effort: the peaks, and each moment its run read) and its parts'
//! stress and deflection under their studies, the worst of each part's
//! (ink-study). Each is started only by a fact nothing before it settled.

use ink::lang::view::ModelView;
use ink_claims::Tool;
use std::cell::RefCell;

struct Sim<'a> {
    effort: ink_effort::Effort<'a>,
    studies: ink_study::Studies<'a>,
}

impl ink_cli::Solvers for Sim<'_> {
    fn tools(&self) -> Vec<&dyn Tool> {
        vec![&self.effort, &self.studies]
    }
    fn at(&self, moment: usize) -> Option<Box<dyn Tool + '_>> {
        self.effort.at(moment).map(|tool| Box::new(tool) as Box<dyn Tool + '_>)
    }
    fn steps(&self) -> Vec<f64> {
        self.effort.steps()
    }
    fn coverage(&self) -> ink_cli::Coverage {
        let (span, goes_on) = self.effort.span();
        ink_cli::Coverage {
            span,
            goes_on,
            loaded: self.effort.loaded(),
            holds: self.effort.holds(),
            unheld: self.effort.unheld(),
            assumed: self.effort.assumed(),
            unbalanced: self.effort.books(),
            balanced: self.effort.balanced().0,
            balanced_moments: self.effort.balanced().1,
            loads_left: self.effort.loads_left(),
        }
    }
}

/// The solvers for the assembly `top` of the design `view`, which every tool asked shares.
pub fn solvers<'a>(view: &'a RefCell<ModelView>, top: &str) -> Result<Box<dyn ink_cli::Solvers + 'a>, String> {
    Ok(Box::new(Sim { effort: ink_effort::Effort::new(view, top)?, studies: ink_study::Studies::new(view, top) }))
}
