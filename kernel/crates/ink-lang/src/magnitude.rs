//! Physical ranges on values consumed by a form. Keep the field so a new
//! parameter scope or moment checks the same quantity the form will use.
use super::*;
use ink_core::interval::I;
use ink_kernel::ranged::{Cmp, Comparison, ProofOptions, Status, value_span};

/// A physical quantity's allowed interval, optionally with an isolated
/// massless zero. `why` explains an out-of-range value in the model's units.
#[derive(Clone, Copy, Debug)]
pub struct Magnitude {
    pub name: &'static str,
    pub range: (f64, f64),
    pub strict: bool,
    pub zero: bool,
    pub why: fn(f64) -> String,
}

impl Magnitude {
    pub(super) fn covers(self, v: I) -> bool {
        let (lo, hi) = self.range;
        (self.zero && v == I::point(0.0)) || if self.strict { lo < v.lo && v.hi < hi } else { lo <= v.lo && v.hi <= hi }
    }

    pub(super) fn expected(self) -> String {
        let (a, b) = self.range;
        if self.strict {
            format!("({a}, {b})")
        } else if self.zero {
            format!("{a}..{b} or zero")
        } else {
            format!("{a}..{b}")
        }
    }

    pub fn error(self, v: f64) -> Option<String> {
        (!self.covers(I::point(v))).then(|| format!("{} ({} = {v}, expected {})", (self.why)(v), self.name, self.expected()))
    }

    // For fading density, the lower material bound applies to the peak.
    // Nonnegativity and the upper bound still apply everywhere.
    fn margin(self, g: &mut Graph, f: NodeId) -> NodeId {
        let (lo, hi) = (g.constant(if self.zero { 0.0 } else { self.range.0 }), g.constant(self.range.1));
        let (above, below) = (g.sub(f, lo), g.sub(hi, f));
        g.min(above, below)
    }
}

pub(super) struct Guard {
    rule: Magnitude,
    field: NodeId,
    margin: NodeId,
    /// The part and the subset of its material where a spatial property fails.
    material: Option<(NodeId, NodeId)>,
    // Material where density reaches the floor, and where it is positive.
    density_regions: Option<(NodeId, NodeId)>,
    at: InkError,
}

impl Doc {
    pub(super) fn guard_magnitude(&mut self, at: &Sexp, field: NodeId, shape: Option<NodeId>, rule: Magnitude) -> R<()> {
        if let Some(guard) = self.make_guard(self.err(at, ""), field, shape, rule)? {
            self.magnitudes.push(guard);
        }
        Ok(())
    }

    /// Check the attribute selected by `ink mass --density NAME` in this
    /// command's context, without registering a rule on other uses of it.
    pub fn check_density_attribute(&mut self, key: &str, who: &str, attrs: &[(String, Value)], shape: NodeId, field: NodeId, ctx: &Ctx) -> R<()> {
        let site = format!("--density {key} on {who}");
        if let Some((_, Value::Field(_, Some(s)))) = attrs.iter().find(|(name, _)| name == key)
            && let Some(why) = s.misfit(Quantity::Density, &site)
        {
            return Err(InkError::msg(why));
        }
        let at = InkError::msg(&site);
        let guard = self.make_guard(at, field, Some(shape), materials::DENSITY).map_err(|e| InkError::msg(format!("{site}: {}", e.message)))?;
        if let Some(why) = guard.and_then(|g| g.refusal(self, ctx)) {
            return Err(InkError::msg(format!("{site}: {why}")));
        }
        Ok(())
    }

    fn make_guard(&mut self, at: InkError, field: NodeId, shape: Option<NodeId>, rule: Magnitude) -> R<Option<Guard>> {
        if let Some(v) = self.g.as_const(field) {
            return rule.error(v).map_or(Ok(None), |message| Err(InkError { message, ..at }));
        }
        let spatial = self.g.value_depends_on(field, |v| matches!(v, Var::X | Var::Y | Var::Z));
        let margin = rule.margin(&mut self.g, field);
        let material = if spatial {
            let shape = shape.ok_or_else(|| InkError { message: format!("{} must be a scalar", rule.name), ..at.clone() })?;
            // Material is the open set shape < 0. For an open property
            // interval a margin of exactly zero also violates the bound.
            let violation = if rule.strict {
                let inside = self.g.neg(margin);
                let (valid, invalid) = (self.g.constant(1.0), self.g.constant(-1.0));
                self.g.select(inside, valid, invalid)
            } else {
                margin
            };
            Some((shape, self.g.max(shape, violation)))
        } else {
            None
        };
        let density_regions = material.filter(|_| rule.zero).map(|(shape, _)| {
            let floor = self.g.constant(rule.range.0);
            let above = self.g.sub(field, floor);
            let (outside, inside) = (self.g.constant(1.0), self.g.constant(-1.0));
            let reaches = self.g.select(above, outside, inside);
            let neg = self.g.neg(field);
            (self.g.max(shape, reaches), self.g.max(shape, neg))
        });
        Ok(Some(Guard { rule, field, margin, material, density_regions, at }))
    }

    pub(super) fn check_magnitudes(&self, ctx: &Ctx, motion: Option<(f64, f64)>, failed: &[usize], at_load: bool, roots: Option<&[NodeId]>) -> R<()> {
        if self.magnitudes.is_empty() {
            return Ok(());
        }
        let read = roots.map(|roots| {
            // A solve hides its inputs behind the parameter it sets. Include the
            // forms read by each solve, all the way back to the original inputs.
            let mut fields = roots.to_vec();
            for i in self.solves_read(roots) {
                fields.extend(self.form_inputs(&self.ext[self.solves[i].ext]));
            }
            let mut read: HashSet<_> = self.g.reachable(&fields).into_iter().collect();
            // A placed attribute is expressed in new coordinates. Check its
            // original field on its original material, also for nested placements.
            loop {
                let placed = self.placed_fields.iter().map(|(&(f, _), &w)| (f, w));
                // And one a claim moves itself, (move (attr-of b density) ...): a warp's copy.
                let originals: Vec<_> = placed.chain(self.warped.iter().copied()).filter_map(|(f, w)| (read.contains(&w) && !read.contains(&f)).then_some(f)).collect();
                if originals.is_empty() {
                    break;
                }
                read.extend(self.g.reachable(&originals));
            }
            read
        });
        for guard in &self.magnitudes {
            if read.as_ref().is_some_and(|read| !read.contains(&guard.field)) {
                continue;
            }
            let roots: Vec<_> = std::iter::once(guard.field).chain(guard.material.map(|(shape, _)| shape)).collect();
            let solves = self.solves_read(&roots);
            if solves.iter().any(|i| failed.contains(i) || (!self.provided_always(*i) && self.solves[*i].owns.iter().any(|&k| !ctx.params[k].is_finite()))) {
                continue;
            }
            let on_t = solves.iter().any(|&i| self.solves[i].on_t) || roots.iter().any(|&f| self.g.depends_on(f, |v| v == Var::T));
            if at_load && on_t {
                continue;
            }
            let over = motion.filter(|_| on_t).map(|span| Ctx { t_span: Some(span), t: eval::span_mid(span), ..ctx.clone() });
            let ctx = over.as_ref().unwrap_or(ctx);
            if motion.is_some() && solves.iter().any(|&i| self.solves[i].on_t && self.solves[i].owns.iter().any(|&k| !ctx.p_spans.is_derived(k as u32))) {
                // These solves are checked at the moments they are run.
                continue;
            }
            if let Some(message) = guard.refusal(self, ctx) {
                return Err(InkError { message, ..guard.at.clone() });
            }
        }
        Ok(())
    }
}

impl Guard {
    fn refusal(&self, doc: &Doc, ctx: &Ctx) -> Option<String> {
        let why = self.outside(doc, ctx)?;
        let inputs: Vec<_> = doc
            .params_read(&[self.field])
            .into_iter()
            .map(|k| {
                let name = &doc.params[k].name;
                match ctx.p_spans.get(k as u32) {
                    Some((lo, hi)) => format!("{name}={lo}..{hi}"),
                    None => format!("{name}={}", ctx.params[k]),
                }
            })
            .collect();
        let message = if inputs.is_empty() { why } else { format!("{why} (from {})", inputs.join(", ")) };
        Some(message)
    }

    fn outside(&self, doc: &Doc, ctx: &Ctx) -> Option<String> {
        let field = Tape::new(&doc.g, self.field);
        let b = [if self.material.is_some() { I::ENTIRE } else { I::point(0.0) }; 3];
        let mut range = eval::eval_interval_scope(&field, ctx, b, None);
        if self.rule.covers(range) {
            return None;
        }
        if range.lo == range.hi {
            return self.rule.error(range.lo);
        }
        let margin = Tape::new(&doc.g, self.margin);
        let cmp = Comparison::exact(if self.rule.strict { Cmp::Gt } else { Cmp::Ge }, 0.0);
        let opts = ProofOptions { res: None, finest: 0.0, depth: ink_kernel::scope::PARAM_SPLITS, tol: 0.0, boxes: ink_kernel::query::BOXES };
        let envelope = cmp.holds_on(eval::eval_interval_cases(&margin, ctx, b));
        if let Some((shape, violation)) = self.material {
            if envelope && !self.rule.zero {
                return None;
            }
            let reach = doc.reach(shape, ctx);
            let tape = |f| std::sync::Arc::new(Tape::new(&doc.g, f));
            // With the outer bounds proved everywhere, any material point
            // in the search window can witness the positive density floor.
            if !reach.whole
                && envelope
                && let Some((reaches, _)) = self.density_regions
                && query::prove_empty(&tape(reaches), ctx, reach.lo, reach.hi, None, doc.finest(), 0.0).disproven()
            {
                return None;
            }
            if reach.whole {
                let (lo, hi, _) = query::bounds_quick(&tape(shape), ctx, reach.lo, reach.hi)?;
                let b = std::array::from_fn(|k| I::new(lo[k], hi[k]));
                range = eval::eval_interval_cases(&field, ctx, b);
                if self.rule.covers(range) {
                    return None;
                }
                let empty = |f| query::prove_empty(&tape(f), ctx, lo, hi, None, doc.finest(), 0.0);
                if envelope || cmp.holds_on(eval::eval_interval_cases(&margin, ctx, b)) || empty(violation).exact() {
                    match self.density_regions {
                        None => return None,
                        Some((reaches, positive)) if (range.hi >= self.rule.range.0 && empty(reaches).disproven()) || empty(positive).exact() => return None,
                        _ => {}
                    }
                }
            }
        } else {
            let status = |field: &Tape, cmp| value_span(field, ctx, [0.0; 3], cmp, &opts).status;
            if (envelope || status(&margin, cmp) == Status::Proven)
                && (!self.rule.zero
                    || status(&field, Comparison::exact(Cmp::Lt, self.rule.range.0)) == Status::Disproven
                    || status(&field, Comparison::exact(Cmp::Le, 0.0)) == Status::Proven)
            {
                return None;
            }
        }
        let peak = if self.rule.zero { " (the positive lower limit applies to the peak)" } else { "" };
        Some(format!(
            "{} has field bounds {}..{}; expected {}{peak}. Check the expression and write its inputs in the documented units",
            self.rule.name,
            range.lo,
            range.hi,
            self.rule.expected()
        ))
    }
}
