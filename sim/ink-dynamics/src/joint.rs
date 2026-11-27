//! `(effort-limits LO HI [(accuracy x)])`, `(actuator [LO HI])` and
//! `(stiffness K)` on an instance with a motion: the joint clauses this
//! package owns (Package::joint_clauses), kept with the joint
//! (Joint::clauses).
//!
//! A joint is an actuator where something drives it: the solver works out
//! the effort it supplies. `(effort-limits LO HI)` says that and what it
//! can give, and claims it (`effort-limits:PATH`). `(actuator LO HI)` says
//! the same and claims nothing: the solver still knows what the actuator
//! can give (to share a coordinate's effort among its actuators, to clamp a
//! command, and as `ink effort`'s margin), and what
//! is asked of it is its part's rating (a stepper's pull-out curve, read
//! from the same record), with `within-ratings`. `(actuator)` states none.

use ink_core::graph::{NodeId, Var};
use ink_lang::sexp::Sexp;
use ink_lang::{InkError, Joint, Quantity, Reader};
use std::any::Any;
use std::rc::Rc;
use std::sync::Arc;

/// An actuator's `(effort-limits LO HI [(accuracy x)])`: the effort it can
/// supply along its motion, N mm for a turn and N for a slide, fields of
/// params and t (a torque-speed line reads the drive's rate). This package
/// claims it at each placement (`effort-limits:PATH`), stated where the
/// motion is set, as limits are (`file`, `line`).
#[derive(Clone)]
pub struct EffortLimits {
    pub lo: NodeId,
    pub hi: NodeId,
    pub accuracy: Option<f64>,
    pub file: Rc<str>,
    pub scope: ink_lang::Scope,
    pub line: u32,
    /// Where the rating itself was declared.
    pub declared_at: String,
    /// Its limits are claimed (`effort-limits:PATH`). Not: an `(actuator
    /// ...)`'s, which the solver uses and no claim holds (with none stated,
    /// `lo` and `hi` are no limit: minus and plus infinity).
    pub claimed: bool,
}

/// A joint's `(stiffness K)`: what holds it at its drive's value, N mm per
/// radian a turn, N per mm a slide (a stepper's holding stiffness, a
/// spring's rate), read by its modes and simulations; and where it is
/// written (file:line).
#[derive(Clone)]
pub struct Stiffness {
    pub k: NodeId,
    pub site: String,
}

/// What `effort-limits` keeps (Rc is not Send: the file is a String there).
struct KeptEffort {
    declared_at: String,
    lo: NodeId,
    hi: NodeId,
    accuracy: Option<f64>,
    claimed: bool,
}

/// A joint's actuator, if it is one: what it can give, and whether that is claimed.
pub fn effort(j: &Joint) -> Option<EffortLimits> {
    let c = j.clauses.iter().find(|c| c.package == "dynamics" && matches!(c.head, "effort-limits" | "actuator"))?;
    let e = c.kept.downcast_ref::<KeptEffort>()?;
    Some(EffortLimits { lo: e.lo, hi: e.hi, accuracy: e.accuracy, file: c.file.clone(), scope: c.scope, line: c.line, declared_at: e.declared_at.clone(), claimed: e.claimed })
}

/// A joint's stiffness, if it states one.
pub fn stiffness(j: &Joint) -> Option<&Stiffness> {
    j.clause::<Stiffness>("stiffness")
}

/// Read `effort-limits` or `stiffness` (Package::clause): a turn's effort
/// limits are torques, a slide's forces; a turn's stiffness is a torque per
/// radian (N*mm), a slide's a force per length (N/mm).
pub fn read(r: &mut Reader, c: &Sexp, head: &str, args: &[Sexp]) -> Result<Arc<dyn Any + Send + Sync>, InkError> {
    let turn = r.motion().expect("a joint clause is read with its motion");
    if head == "stiffness" {
        let [k] = args else { return Err(r.err(c, "usage: (stiffness K), K N*mm per radian a turn, N per mm a slide")) };
        let (q, what) = if turn { (Quantity::Torque, "a turn's stiffness, a torque per radian") } else { (Quantity::Stiffness, "a slide's stiffness") };
        let n = r.field_of(k, q, what)?;
        if r.graph().depends_on(n, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
            return Err(r.err(c, "a joint's stiffness is a number or a field of params, not of x y z or t"));
        }
        return Ok(Arc::new(Stiffness { k: n, site: r.site(c) }));
    }
    let claimed = head == "effort-limits";
    let usage = if claimed {
        "usage: (effort-limits LO HI [(accuracy x)])"
    } else {
        "usage: (actuator [LO HI]): the joint is driven, and the solver works out the effort it supplies; LO HI, what it can give, which no claim holds \
         ((effort-limits LO HI) claims it)"
    };
    if !claimed && args.is_empty() {
        let (lo, hi) = (r.graph().constant(f64::NEG_INFINITY), r.graph().constant(f64::INFINITY));
        return Ok(Arc::new(KeptEffort { lo, hi, accuracy: None, declared_at: r.site(c), claimed }));
    }
    let [lo, hi, tol @ ..] = args else { return Err(r.err(c, usage)) };
    if tol.len() > usize::from(claimed) {
        return Err(r.err(c, usage));
    }
    let q = if turn { Quantity::Torque } else { Quantity::Force };
    let (l, h) = (r.field_of(lo, q, head)?, r.field_of(hi, q, head)?);
    // An actuator's envelope may follow its speed: a field of t.
    let space = |v: Var| matches!(v, Var::X | Var::Y | Var::Z);
    if r.graph().depends_on(l, space) || r.graph().depends_on(h, space) {
        return Err(r.err(c, format!("an actuator's limits ({head}) are fields of params and t (the drive's rate), not of x y z")));
    }
    if let (Some(a), Some(b)) = (r.graph().as_const(l), r.graph().as_const(h))
        && a > b
    {
        return Err(r.err(c, format!("{head} {a} > {b}: write ({head} LO HI) with LO <= HI")));
    }
    let accuracy = tol.first().map(|t| r.accuracy(t, Some(q))).transpose()?;
    Ok(Arc::new(KeptEffort { lo: l, hi: h, accuracy, declared_at: r.site(c), claimed }))
}
