//! What this release holds back, refused by name: a model or a command
//! that uses it reads that it is not in this release, never that it is
//! unknown or misspelt. Each is a package not yet ready to be relied on,
//! kept and tested apart; it joins a release as a package again once it
//! is.

use ink::lang::command::{Args, Command, Run};
use ink::lang::sexp::Sexp;
use ink::lang::{Form, InkError, Package, Reader, Registry, Value};
use std::sync::Arc;

/// A package held back: the heads it owned, each refused.
pub struct Withheld {
    /// The package's name (`flow`).
    pub name: &'static str,
    /// What it did, for the refusal (`flow (lattice Boltzmann on the GPU, with conjugate heat)`).
    pub what: &'static str,
    /// Its top-level forms, values and commands.
    pub forms: &'static [&'static str],
    pub values: &'static [&'static str],
    pub commands: &'static [Command],
}

/// Every package this release holds back.
pub static WITHHELD: [&Withheld; 2] = [&FLOW, &LEAST_TIME];

static FLOW: Withheld = Withheld { name: "flow", what: "flow (lattice Boltzmann on the GPU, with conjugate heat)", forms: &["flow"], values: &[], commands: &[refused("flow")] };

static LEAST_TIME: Withheld =
    Withheld { name: "least-time", what: "the least time to run a path within its actuators' limits", forms: &[], values: &["least-time"], commands: &[] };

impl Withheld {
    fn says(&self, used: &str) -> String {
        format!("{used} is not in this release of ink: {} is held back until it is ready to rely on", self.what)
    }
}

/// A command of a withheld package: listed, so `ink help` and MCP say it
/// is held back, and refused whatever it is given.
const fn refused(name: &'static str) -> Command {
    Command {
        name,
        options: &[],
        help: "Not in this release of ink: held back until it is ready to rely on.",
        tool: "Not in this release of ink: held back until it is ready to rely on. Every call is refused, saying so.",
        run: Run::Words { run: refuse, here: |_| true },
    }
}

fn refuse(_: &'static Registry, a: &Args, _: &(dyn Fn() -> bool + Sync)) -> Result<serde_json::Value, InkError> {
    let w = WITHHELD.iter().find(|w| w.commands.iter().any(|c| c.name == a.cmd)).expect("a withheld package's command");
    Err(InkError::msg(w.says(&format!("`ink {}`", a.cmd))))
}

impl Package for Withheld {
    fn name(&self) -> &'static str {
        self.name
    }

    fn forms(&self) -> &'static [&'static str] {
        self.forms
    }

    fn values(&self) -> &'static [&'static str] {
        self.values
    }

    fn commands(&self) -> &'static [Command] {
        self.commands
    }

    fn read(&self, r: &mut Reader, form: &Sexp, head: &str, _: &[Sexp]) -> Result<Arc<dyn Form>, InkError> {
        Err(r.err(form, self.says(&format!("({head} ...)"))))
    }

    fn value(&self, r: &mut Reader, e: &Sexp, head: &str, _: &[Sexp]) -> Result<Value, InkError> {
        Err(r.err(e, self.says(&format!("({head} ...)"))))
    }
}
