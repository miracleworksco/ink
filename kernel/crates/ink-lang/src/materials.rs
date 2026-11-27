//! Engineering materials and fluids: the library std/materials.ink, plain IR
//! data with explicit grades, source comments and compatibility aliases.
//! Another spelling (another case, spaces or underscores for hyphens, or a
//! name the library lists in `(def other-names (data (aluminium
//! aluminum-6061) ...))`) is refused with the name to write. The kernel
//! holds no table of its own; it reads the library's defs by name. Evidence
//! strings describe the input; they do not change a claim's verdict kind.

use super::*;
use ink_core::eval::eval_interval;
use ink_core::interval::{I, add_o, div_o, mul_o, sub_o};
use std::sync::Arc;
use std::sync::OnceLock;

/// A published mass defines a uniform effective density over what has
/// none of its own: one body's (`(attr mass M)` on it), or the bodies an
/// assembly places (`(attr mass M)` on it), each placed as often as the
/// assembly places it. Bodies there with a density of their own (a steel
/// rotor) weigh as it says, and the rest share what remains. It is
/// measured at the configuration the mass is published for: the params'
/// defaults, or those `(attr mass-params ...)` gives, but the params the
/// mass reads, which keep their values (it is published at each); so a
/// size param scales the mass, as it scales the volume. A body's motion over t is a pose:
/// measured at each moment, the mass stays. A derived param, so every
/// consumer sees one density field; its volumes are measured, so claims
/// reading it are estimates, never proofs; one that does not settle still
/// sets it, and the answers say so.
#[derive(Clone)]
struct PublishedMass {
    name: String,
    /// The shapes sharing the density, each in its own coordinates, and
    /// how many times each is placed.
    shared: Vec<(Arc<Tape>, f64)>,
    shared_shapes: Vec<(NodeId, f64)>,
    /// Bodies with a density of their own: shape, density, times placed.
    own: Vec<KnownBody>,
    mass: Tape,
    /// The published configuration: params and their values.
    at: Vec<(usize, f64)>,
    /// The moment it is published at, where a body sharing it deforms (KeptShare): its density is
    /// found there, and each body keeps the mass it has there. None: at every moment.
    moment: Option<f64>,
    owns: [usize; 1],
    uses_t: bool,
}

#[derive(Clone)]
struct KnownBody {
    shape_id: NodeId,
    density_param: Option<usize>,
    shape: Arc<Tape>,
    density: Arc<Tape>,
    copies: f64,
    /// For a solid and density invariant under every turn about z, their
    /// section at y=0 and the density there times the cylindrical Jacobian.
    radial: Option<(Arc<Tape>, Arc<Tape>)>,
}

impl KnownBody {
    fn new(g: &mut Graph, s: NodeId, rho: NodeId, copies: f64) -> Self {
        use ink_core::symmetry::{Order, turn_symmetry};
        let invariant = |n| {
            let s = turn_symmetry(g, n);
            s.order == Order::Every && s.slack == 0.0
        };
        let radial = (invariant(s) && invariant(rho)).then(|| {
            let zero = g.constant(0.0);
            let section = g.subst(s, &[(Var::Y, zero)]);
            let density = g.subst(rho, &[(Var::Y, zero)]);
            let x = g.var(Var::X);
            let full_turn = g.constant(360.0);
            let tau = g.un_real(UnOp::Rad, full_turn);
            let circumference = g.mul(tau, x);
            let weighted = g.mul(density, circumference);
            (Arc::new(Tape::new(g, section)), Arc::new(Tape::new(g, weighted)))
        });
        let density_param = match g.nodes[rho as usize] {
            graph::Node::Var(Var::Param(k)) => Some(k as usize),
            _ => None,
        };
        Self { shape_id: s, density_param, shape: Arc::new(Tape::new(g, s)), density: Arc::new(Tape::new(g, rho)), copies, radial }
    }

    /// Even unweighed, the material occupies at most its box.
    fn box_mass(&self, ctx: &Ctx, (lo, hi): (V3, V3)) -> I {
        let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
        let volume = b.iter().fold(I::point(self.copies), |v, side| mul_o(v, sub_o(I::point(side.hi), I::point(side.lo))));
        mul_o(I::new(0.0, volume.hi), eval_interval(&self.density, ctx, b, None))
    }
}

struct KnownMeasurement<'a> {
    body: &'a KnownBody,
    region: (V3, V3),
    estimate: f64,
    range: I,
    /// None until the body has been weighed.
    accuracy: Option<f64>,
    limit: Option<String>,
}

impl KnownMeasurement<'_> {
    fn error(&self) -> f64 {
        (self.estimate - self.range.lo).max(self.range.hi - self.estimate)
    }

    fn priority(&self) -> f64 {
        // Each integration starts over. Prefer coarser measurements when
        // further narrowing a body would cost more for the remaining error.
        self.error() * self.accuracy.unwrap_or(query::ACCURACY)
    }

    fn can_refine(&self) -> bool {
        self.range.lo < self.range.hi && self.accuracy.is_none_or(|a| a == query::ACCURACY || (a > f64::EPSILON && self.limit.is_none()))
    }
}

/// The same measured volume occurs in both the denominator of the
/// density and the mass of these bodies. For a complete set of copies
/// at the same configuration, that volume cancels exactly.
#[derive(Debug, PartialEq)]
pub(crate) struct Allocation {
    shapes: Vec<(NodeId, f64)>,
    configuration: Vec<(u32, f64, Option<(f64, f64)>)>,
    estimate: f64,
    mass: I,
}

impl Allocation {
    fn copies(&self, bodies: &[KnownBody], param: usize, ctx: &Ctx) -> Option<f64> {
        if !self.configuration.iter().all(|&(k, value, span)| ctx.params[k as usize] == value && ctx.p_spans.get(k) == span) {
            return None;
        }
        let mut counts = std::collections::BTreeMap::<NodeId, f64>::new();
        for b in bodies.iter().filter(|b| b.density_param == Some(param)) {
            let n = counts.entry(b.shape_id).or_default();
            *n = graph::exact_bin(BinOp::Add, *n, b.copies)?;
        }
        let mut copies = None;
        for &(shape, n) in &self.shapes {
            let ratio = graph::exact_bin(BinOp::Div, counts.remove(&shape)?, n)?;
            if copies.is_some_and(|k| k != ratio) {
                return None;
            }
            copies = Some(ratio);
        }
        copies.filter(|_| counts.is_empty())
    }
}

/// A deforming body's density where it shares a published mass found at `moment`: the mass it has
/// there (the shared density `of` times its volume there) over its volume at each moment, so it keeps
/// its mass as it deforms, as the bodies beside it keep theirs.
#[derive(Clone)]
struct KeptShare {
    name: String,
    shape: Arc<Tape>,
    of: usize,
    moment: f64,
    owns: [usize; 1],
}

impl KeptShare {
    fn volume(&self, ctx: &Ctx) -> Result<(f64, I), String> {
        let reach = crate::reach_of(&self.shape, DEFAULT_DOMAIN, ctx);
        if !reach.whole {
            return Err(format!("{}: a published mass needs bounded bodies", self.name));
        }
        let Some((lo, hi, _)) = query::bounds(&self.shape, ctx, reach.lo, reach.hi) else { return Ok((0.0, I::point(0.0))) };
        let v = query::measured(&self.shape, &[None], true, ctx, lo, hi, query::Until::Accuracy(query::ACCURACY));
        Ok((v.sums[0], I::new(v.lo[0], v.hi[0])))
    }
}

impl Form for KeptShare {
    fn name(&self) -> &str {
        &self.name
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for KeptShare {
    fn name(&self) -> &str {
        &self.name
    }
    fn owns(&self) -> &[usize] {
        &self.owns
    }
    fn uses_t(&self) -> bool {
        true
    }
    fn on_demand(&self) -> bool {
        true
    }
    fn warm_starts(&self) -> bool {
        false
    }
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        self.provide_with_bounds(params, t, &[], &mut Vec::new())
    }
    fn provide_with_bounds(&self, params: &mut [f64], t: f64, inputs: &[ProvidedParam], bounds: &mut Vec<ProvidedParam>) -> Result<Option<String>, Refused> {
        let ((at, at_range), (now, now_range)) = (self.volume(&Ctx::new(self.moment, params.to_vec()))?, self.volume(&Ctx::new(t, params.to_vec()))?);
        if !(now > 0.0 && now_range.lo > 0.0) {
            return Err(Refused::from(format!("{}: the body has no volume at t = {t} to keep its mass in", self.name)));
        }
        let rho = params[self.of];
        let rho_range = inputs.iter().find(|b| b.param == self.of).map_or(I::point(rho), |b| b.range);
        params[self.owns[0]] = rho * at / now;
        bounds.push(ProvidedParam { param: self.owns[0], range: div_o(mul_o(rho_range, at_range), now_range), allocation: None });
        Ok(None)
    }
}

impl Form for PublishedMass {
    fn name(&self) -> &str {
        &self.name
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for PublishedMass {
    fn name(&self) -> &str {
        &self.name
    }
    fn owns(&self) -> &[usize] {
        &self.owns
    }
    fn uses_t(&self) -> bool {
        self.uses_t
    }
    /// Measured where a mass is read, not by every command on the model.
    fn on_demand(&self) -> bool {
        true
    }
    fn warm_starts(&self) -> bool {
        false
    }

    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        self.provide_with_bounds(params, t, &[], &mut Vec::new())
    }

    /// A mass written with a tolerance (a ranged param only its mass reads, no shape and not its
    /// published configuration): its density at the mass's value, enclosed for every value of the
    /// mass by its bounds over them (`provide_with_bounds` from the ranged params' spans). None where
    /// a shape reads one: its volume is measured at single values only.
    fn over_ranged(&self, params: &[f64], t: f64, ranged: &[(usize, (f64, f64))]) -> Option<Result<(f64, (f64, f64)), Refused>> {
        let mass_reads = self.mass.params_read();
        let shapes: Vec<u32> =
            self.shared.iter().flat_map(|(s, _)| s.params_read()).chain(self.own.iter().flat_map(|b| b.shape.params_read().into_iter().chain(b.density.params_read()))).collect();
        let only_mass = ranged.iter().all(|&(k, _)| mass_reads.contains(&(k as u32)) && !shapes.contains(&(k as u32)) && !self.at.iter().any(|a| a.0 == k));
        if !only_mass {
            return None;
        }
        let inputs: Vec<ProvidedParam> = ranged.iter().map(|&(k, (lo, hi))| ProvidedParam { param: k, range: I::new(lo, hi), allocation: None }).collect();
        let (mut values, mut bounds) = (params.to_vec(), Vec::new());
        Some(self.provide_with_bounds(&mut values, t, &inputs, &mut bounds).and_then(|_| {
            let k = self.owns[0];
            match bounds.iter().find(|b| b.param == k) {
                Some(b) if b.range.lo.is_finite() && b.range.hi.is_finite() => Ok((values[k], (b.range.lo, b.range.hi))),
                _ => Err(Refused::from(format!("{}: its density is not enclosed over the mass's tolerance", self.name))),
            }
        }))
    }

    fn provide_with_bounds(&self, params: &mut [f64], t: f64, inputs: &[ProvidedParam], bounds: &mut Vec<ProvidedParam>) -> Result<Option<String>, Refused> {
        let (ctx, mass, mass_bounds) = self.at(params, t, inputs)?;
        let (own, own_bounds, mut unsettled) = self.known_mass(&ctx, mass, mass_bounds, inputs, true)?;
        let (mut volume, mut volume_bounds) = (0.0, I::point(0.0));
        for (s, n) in &self.shared {
            let (v, range, why) = self.measure(s, &ctx)?;
            volume += n * v;
            volume_bounds = add_o(volume_bounds, mul_o(I::point(*n), range));
            unsettled = unsettled.or(why);
        }
        if !(volume.is_finite() && volume > 0.0) {
            return Err(Refused::from(self.empty()));
        }
        if own_bounds.lo > mass_bounds.hi {
            return Err(Refused::from(self.outweighed(own_bounds.lo, mass)));
        }
        let remainder = sub_o(mass_bounds, own_bounds);
        let assigned = I::new(remainder.lo.max(0.0), remainder.hi);
        let range = div_o(assigned, volume_bounds);
        let density = (mass - own).max(0.0) / volume;
        if !DENSITY.covers(I::point(density)) {
            return Err(Refused::from(format!(
                "{}: {mass} g published mass implies density {density} g/mm³ after accounting for {} g over {volume} mm³. Its bounds {}..{} do not establish the expected {}. Check the volume and write the published mass with its unit, e.g. (attr mass 2.7kg)",
                self.name,
                own.abs(),
                range.lo,
                range.hi,
                DENSITY.expected()
            )));
        }
        if volume_bounds.lo > 0.0 {
            let mut ks: Vec<u32> = self.shared.iter().flat_map(|(s, _)| s.params_read()).collect();
            ks.sort_unstable();
            ks.dedup();
            let allocation = Allocation {
                shapes: self.shared_shapes.clone(),
                configuration: ks.into_iter().map(|k| (k, ctx.params[k as usize], ctx.p_spans.get(k))).collect(),
                estimate: (mass - own).max(0.0),
                mass: assigned,
            };
            bounds.push(ProvidedParam { param: self.owns[0], range, allocation: Some(Arc::new(allocation)) });
        }
        params[self.owns[0]] = density;
        if remainder.lo < 0.0 {
            return Ok(Some("the bounds on the known bodies straddle the published mass; feasibility is undecided".into()));
        }
        Ok(unsettled.map(|why| format!("the known mass or shared volume is not measured to the requested accuracy: {why}. The density uses the finest measurement reached")))
    }

    /// Without measuring the volume it shares where it can: its mass, its
    /// bodies' bounds, a volume above 0 (`holds_material`, else measured),
    /// and what its bodies with a density of their own weigh, measured only
    /// where their density's bound over their box could pass the mass.
    fn refusal(&self, params: &[f64], t: f64) -> Option<Refused> {
        self.feasibility(params, t, Some(&[])).err()
    }

    fn feasibility(&self, params: &[f64], t: f64, inputs: Option<&[ProvidedParam]>) -> Result<Option<String>, Refused> {
        let Some(inputs) = inputs else {
            return Ok(Some("published mass feasibility depends on measured inputs without bounds".into()));
        };
        let (ctx, mass, mass_bounds) = self.at(params, t, inputs)?;
        // A volume above 0: a cell of one body proven inside, else
        // the bodies measured.
        let mut solid = false;
        for (s, _) in &self.shared {
            if let Some((lo, hi)) = self.bounds(s, &ctx)? {
                solid = solid || holds_material(s, &ctx, lo, hi);
            }
        }
        if !solid {
            let mut volume = I::point(0.0);
            for (s, n) in &self.shared {
                volume = add_o(volume, mul_o(I::point(*n), self.measure(s, &ctx)?.1));
            }
            if volume.hi <= 0.0 {
                return Err(Refused::from(self.empty()));
            }
            if !(volume.lo > 0.0 && volume.hi.is_finite()) {
                return Ok(Some(format!("{}: the measurement cannot establish positive finite volume", self.name)));
            }
        }
        let (mut most, mut unboxed) = (I::point(0.0), false);
        for body in &self.own {
            match self.bounds(&body.shape, &ctx)? {
                Some(region) => most = add_o(most, I::point(body.box_mass(&ctx, region).hi)),
                // Unbounded, it may weigh anything: the boxes' density below says nothing.
                None => unboxed = true,
            }
        }
        // Its bodies' boxes hold at most so much: a mass that implies a
        // density past any material's over them is refused here, as
        // weighing it would (`ink check` says so unweighed).
        let mut room = I::point(0.0);
        for (s, n) in &self.shared {
            match self.bounds(s, &ctx)? {
                Some((lo, hi)) => room = add_o(room, mul_o(I::point(*n), I::point((0..3).map(|k| hi[k] - lo[k]).product()))),
                None => room = I::new(0.0, f64::INFINITY),
            }
        }
        let least = (mass_bounds.lo - most.hi).max(0.0) / room.hi;
        if !unboxed && room.hi > 0.0 && least > DENSITY_RANGE.1 {
            return Err(Refused::from(format!(
                "{}: {mass} g published mass implies density at least {least} g/mm³: its bodies' boxes hold at most {} mm³{}. That is past the expected {}. Check the volume and write the published mass with its unit, e.g. (attr mass 2.7kg)",
                self.name,
                room.hi,
                if most.hi > 0.0 { format!(", less at most {} g its other bodies weigh", most.hi) } else { String::new() },
                DENSITY.expected()
            )));
        }
        if most.hi <= mass_bounds.lo {
            return Ok(None);
        }
        let (_, own, _) = self.known_mass(&ctx, mass, mass_bounds, inputs, false)?;
        if own.lo > mass_bounds.hi {
            return Err(Refused::from(self.outweighed(own.lo, mass)));
        }
        Ok((own.hi > mass_bounds.lo).then(|| format!("{}: the known bodies' mass bounds [{}, {}] g straddle the published mass", self.name, own.lo, own.hi)))
    }
}

impl PublishedMass {
    /// The context it is measured in (its configuration) and its mass.
    fn at(&self, params: &[f64], t: f64, inputs: &[ProvidedParam]) -> Result<(Ctx, f64, I), String> {
        let mut values = params.to_vec();
        for &(k, v) in &self.at {
            values[k] = v;
        }
        let mut ctx = Ctx::new(self.moment.unwrap_or(t), values);
        let bands: Vec<_> = inputs.iter().map(|b| (b.param as u32, (b.range.lo, b.range.hi))).collect();
        ctx.p_spans = ctx.p_spans.with_derived(&bands);
        let mass = eval_point(&self.mass, &ctx, [0.0; 3]);
        if !mass.is_finite() || mass < 0.0 {
            return Err(format!("{}: mass must be finite and nonnegative (grams)", self.name));
        }
        let bounds = eval_interval(&self.mass, &ctx, [I::point(0.0); 3], None);
        Ok((ctx, mass, bounds))
    }

    /// Enclose the sum, refining first where the most uncertainty comes from.
    /// Checking feasibility only needs to separate the bounds from the total.
    /// Deriving a density also needs accuracy relative to the mass left over.
    fn known_mass(&self, ctx: &Ctx, mass: f64, mass_bounds: I, inputs: &[ProvidedParam], need_remainder: bool) -> Result<(f64, I, Option<String>), String> {
        let (mut assigned, mut assigned_bounds) = (0.0, I::point(0.0));
        let mut accounted = Vec::new();
        for input in inputs {
            if let Some(a) = &input.allocation
                && let Some(n) = a.copies(&self.own, input.param, ctx)
            {
                assigned += n * a.estimate;
                assigned_bounds = add_o(assigned_bounds, mul_o(I::point(n), a.mass));
                accounted.push(input.param);
            }
        }
        let mut measurements = Vec::new();
        for body in self.own.iter().filter(|b| !b.density_param.is_some_and(|p| accounted.contains(&p))) {
            let Some(region) = self.bounds(&body.shape, ctx)? else { continue };
            let (estimate, range, accuracy, limit) = if need_remainder {
                let (estimate, range, limit) = Self::measure_in(&body.shape, Some(&body.density), ctx, region);
                (body.copies * estimate, mul_o(I::point(body.copies), range), Some(query::ACCURACY), limit)
            } else {
                let range = body.box_mass(ctx, region);
                (range.mid(), range, None, None)
            };
            measurements.push(KnownMeasurement { body, region, estimate, range, accuracy, limit });
        }
        loop {
            let (mut own, mut range) = (assigned, assigned_bounds);
            let (mut fixed, mut active) = (assigned_bounds, I::point(0.0));
            for m in &measurements {
                own += m.estimate;
                range = add_o(range, m.range);
                if m.can_refine() {
                    active = add_o(active, m.range);
                } else {
                    fixed = add_o(fixed, m.range);
                }
            }
            let goal = query::ACCURACY * (mass - own).abs();
            let error = (own - range.lo).max(range.hi - own);
            let fits = range.hi <= mass_bounds.lo;
            let why = measurements.iter().find_map(|m| m.limit.clone());
            if range.lo > mass_bounds.hi || (need_remainder && error <= goal) || (!need_remainder && fits) {
                return Ok((own, range, why));
            }

            // Even with the other bodies measured exactly, the fixed bounds
            // can prevent a verdict or swamp the accuracy of any remainder.
            let least_upper = add_o(I::point(fixed.hi), I::point(active.lo)).lo;
            let greatest_lower = add_o(I::point(fixed.lo), I::point(active.hi)).hi;
            let can_decide = !fits && (least_upper <= mass_bounds.lo || greatest_lower > mass_bounds.hi);
            let biggest_remainder = (mass - range.lo).abs().max((mass - range.hi).abs());
            let can_settle = need_remainder && fixed.hi - fixed.lo <= 2.0 * query::ACCURACY * biggest_remainder;
            let stopped = ink_kernel::stop::reached();
            let next =
                if (can_decide || can_settle) && !stopped { measurements.iter_mut().filter(|m| m.can_refine()).max_by(|a, b| a.priority().total_cmp(&b.priority())) } else { None };
            let Some(m) = next else {
                let why = why.or_else(|| {
                    (need_remainder && !measurements.is_empty())
                        .then(|| if stopped { "the call's time limit passed" } else { "the remainder is below measurement precision" }.into())
                });
                return Ok((own, range, why));
            };

            let target = if need_remainder { goal } else { (mass_bounds.lo - own).max(own - mass_bounds.hi).max(0.0) };
            // Clear the error target in one measurement where possible.
            // Otherwise ask for the fractional reduction the sum needs.
            let accuracy = m.accuracy.map_or(query::ACCURACY, |a| {
                let other_error = ((own - range.lo) - (m.estimate - m.range.lo)).max((range.hi - own) - (m.range.hi - m.estimate));
                let desired = if other_error < target { target - other_error } else { m.error() * (target / error) };
                (desired.min(m.error() * 0.5) / m.estimate.abs()).min(a * 0.5).max(f64::EPSILON)
            });
            let body = m.body;
            let (lo, hi) = m.region;
            let (shape, density, lo, hi) = match &body.radial {
                Some((section, weight)) => (section, weight, [0.0, 0.0, lo[2]], [hi[0], 1.0, hi[2]]),
                None => (&body.shape, &body.density, lo, hi),
            };
            let v = query::measured(shape, &[Some(density)], false, ctx, lo, hi, query::Until::Accuracy(accuracy));
            let measured = mul_o(I::point(body.copies), I::new(v.lo[0], v.hi[0]));
            let (lo, hi) = (m.range.lo.max(measured.lo), m.range.hi.min(measured.hi));
            m.range = if lo <= hi { I::new(lo, hi) } else { measured };
            m.estimate = (body.copies * v.sums[0]).clamp(m.range.lo, m.range.hi);
            m.accuracy = Some(accuracy);
            m.limit = refinement_limit(&v);
        }
    }

    /// A body's box (None: empty), or why a mass over it cannot be measured.
    fn bounds(&self, shape: &Arc<Tape>, ctx: &Ctx) -> Result<Option<(V3, V3)>, String> {
        let reach = crate::reach_of(shape, DEFAULT_DOMAIN, ctx);
        if !reach.whole {
            return Err(format!("{}: a published mass needs bounded bodies", self.name));
        }
        Ok(query::bounds(shape, ctx, reach.lo, reach.hi).map(|(lo, hi, _)| (lo, hi)))
    }

    /// A shared body's volume, as `ink mass` measures it.
    fn measure(&self, shape: &Arc<Tape>, ctx: &Ctx) -> Result<(f64, I, Option<String>), String> {
        let Some(region) = self.bounds(shape, ctx)? else { return Ok((0.0, I::point(0.0), None)) };
        Ok(Self::measure_in(shape, None, ctx, region))
    }

    fn measure_in(shape: &Arc<Tape>, rho: Option<&Tape>, ctx: &Ctx, (lo, hi): (V3, V3)) -> (f64, I, Option<String>) {
        let v = weighed(shape, rho, ctx, lo, hi, query::Until::Accuracy(query::ACCURACY));
        (v.sums[0], I::new(v.lo[0], v.hi[0]), v.unsettled())
    }

    fn empty(&self) -> String {
        format!("{}: a published mass needs a body with positive finite volume", self.name)
    }

    fn outweighed(&self, own: f64, mass: f64) -> String {
        format!("{}: its bodies with a density of their own weigh {own:.4} g, more than its published {mass} g", self.name)
    }
}

/// A bounded body's moments under `rho` (None: its volume) in its box
/// `lo..hi`, as `ink mass` (ModelView::moments) and a published mass
/// (PublishedMass::measure) both weigh it: on the field itself, so the
/// volume a published density divides by is the one each body's share is
/// weighed over, and the mass of a body beside them with its own density
/// is the one taken from it; the shares sum to the published mass. Clipped
/// to its reach a bounded field changes no sign, only its tape, and with
/// it the search.
pub fn weighed(shape: &Arc<Tape>, rho: Option<&Tape>, ctx: &Ctx, lo: V3, hi: V3, until: query::Until) -> query::Measured {
    query::measured(shape, &[rho, None], true, ctx, lo, hi, until)
}

/// A stable sampled estimate can still have wide proven bounds.
fn refinement_limit(m: &query::Measured) -> Option<String> {
    m.unsettled().or_else(|| m.sampled.map(|_| "the remaining uncertainty exceeds the integration budget".into()))
}

/// Whether a shape in `lo..hi` holds a cell of material (so a volume
/// above 0): a point the emptiness search finds inside, confirmed, and a
/// box round it a billionth of the shape's size bounded below 0. False
/// where that is not shown.
fn holds_material(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> bool {
    let scale = query::max_extent(lo, hi);
    let floor = 1e-6 * scale;
    if !(floor > 0.0 && floor.is_finite()) {
        return false;
    }
    let query::Verdict::Disproven { witness, .. } = query::prove_empty(shape, ctx, lo, hi, None, floor, 0.0).verdict else { return false };
    let Some(p) = query::confirm_inside(shape, ctx, witness, scale) else { return false };
    let e = 1e-9 * scale;
    eval_interval(shape, ctx, p.map(|x| I::new(x - e, x + e)), None).all_inside()
}

/// Why an assembly's published mass is not used, where every body it
/// places has a density of its own.
pub const UNUSED_MASS: &str = "gives no body a density: every body it places has its own (a material, a density, or an inner assembly's mass), \
    so they weigh what those say and this mass is not used";

pub const LIB: &str = include_str!("../../../std/materials.ink");

/// The library's def of the names it refuses: data, not a material.
const OTHER_NAMES: &str = "other-names";

/// Canonical entries, compatibility aliases (direct targets only), and
/// refused spellings (name, the one to write).
struct Library {
    entries: Vec<(String, serde_json::Value)>,
    aliases: Vec<(String, String)>,
    instead: Vec<(String, String)>,
}

fn library() -> &'static Library {
    static T: OnceLock<Library> = OnceLock::new();
    T.get_or_init(|| {
        let forms = parse(LIB).expect("std/materials.ink parses");
        let mut lib = Library { entries: Vec::new(), aliases: Vec::new(), instead: Vec::new() };
        for f in &forms {
            // (def name (data (key value) ...)): pairs of a symbol and a value.
            let Some([def, name, value]) = f.list() else { panic!("std/materials.ink: {f}: every form is (def name (data ...))") };
            if let (Some("def"), Some(name), Some(to)) = (def.sym(), name.sym(), value.sym()) {
                lib.aliases.push((name.to_string(), to.to_string()));
                continue;
            }
            let (Some("def"), Some(name), Some([data, pairs @ ..])) = (def.sym(), name.sym(), value.list()) else {
                panic!("std/materials.ink: {f}: every form is (def name (data ...))")
            };
            assert!(data.sym() == Some("data"), "std/materials.ink: {name}: every form is (def name (data ...))");
            let pairs = pairs.iter().filter_map(|p| match p.list() {
                Some([k, v]) => Some((k.sym()?, v)),
                _ => None,
            });
            if name == OTHER_NAMES {
                lib.instead.extend(pairs.filter_map(|(k, to)| Some((k.to_string(), to.sym()?.to_string()))));
            } else {
                // Plain data includes the condition and evidence labels. Keep
                // strings as well as scalars, just as the IR's data reader does.
                let m = pairs.map(|(k, v)| (k.to_string(), v.to_json())).collect();
                lib.entries.push((name.to_string(), serde_json::Value::Object(m)));
            }
        }
        for (name, to) in &lib.aliases {
            assert!(!lib.entries.iter().any(|(n, _)| n == name), "std/materials.ink: duplicate alias {name}");
            assert!(lib.entries.iter().any(|(n, _)| n == to), "std/materials.ink: alias {name}: no material {to}");
        }
        for (name, to) in &lib.instead {
            assert!(lib.entries.iter().any(|(n, _)| n == to) || lib.aliases.iter().any(|(n, _)| n == to), "std/materials.ink: {OTHER_NAMES}: ({name} {to}): no material {to}");
        }
        lib
    })
}

/// Solid materials' names (those with a Young's modulus), for hints.
pub fn names() -> Vec<&'static str> {
    library().entries.iter().filter(|(_, v)| v.get("E").is_some()).map(|(n, _)| n.as_str()).collect()
}

/// Fluids' names (those with a viscosity), for hints.
pub fn fluids() -> Vec<&'static str> {
    library().entries.iter().filter(|(_, v)| v.get("viscosity").is_some()).map(|(n, _)| n.as_str()).collect()
}

/// Properties, condition and evidence of a material or fluid, with the
/// requested `name` and resolved `canonical-name`.
pub fn lookup(name: &str) -> Option<serde_json::Value> {
    let lib = library();
    let canonical = lib.aliases.iter().find(|(n, _)| n == name).map_or(name, |(_, to)| to);
    lib.entries.iter().find(|(n, _)| n == canonical).map(|(n, v)| {
        let mut v = v.clone();
        v["name"] = serde_json::json!(name);
        v["canonical-name"] = serde_json::json!(n);
        v
    })
}

/// A body's property `key` from its `(attr material ...)` (a library name,
/// or data with the keys), if it has one with that key.
pub fn material_property(attrs: &[(String, Value)], key: &str) -> Option<f64> {
    let m = match &attrs.iter().find(|(k, _)| k == "material")?.1 {
        Value::Str(n) => lookup(n)?,
        Value::Data(d) => d.value().clone(),
        _ => return None,
    };
    m.get(key)?.as_f64()
}

/// A body's attribute `key`, as every reader takes it (a claim's `(attr
/// key)`, `(by-body ...)`, `ink mass`): its own
/// `(attr key ...)`, else its material's `key`. None where it has neither;
/// an error, naming the target as `who` (Doc::named), where its own is not
/// a number or a field.
pub fn attr(g: &mut Graph, who: &str, attrs: &[(String, Value)], key: &str) -> Result<Option<NodeId>, InkError> {
    Ok(attr_value(g, who, attrs, key)?.map(|(f, _)| f))
}

type ScalarAttr = (NodeId, Option<Rc<units::Said>>);

fn attr_value(g: &mut Graph, who: &str, attrs: &[(String, Value)], key: &str) -> Result<Option<ScalarAttr>, InkError> {
    match attrs.iter().find(|(k, _)| k == key) {
        Some((_, Value::Field(f, said))) => Ok(Some((*f, said.clone()))),
        Some((_, other)) => Err(InkError::msg(format!("attr '{key}' of {who} is a {}, not a number or a field", other.kind()))),
        None => Ok(material_property(attrs, key).map(|v| (g.constant(v), None))),
    }
}

/// `attr`, where the body must have it: else `missing`.
pub fn required(g: &mut Graph, who: &str, attrs: &[(String, Value)], key: &str) -> Result<NodeId, InkError> {
    attr(g, who, attrs, key)?.ok_or_else(|| missing(who, attrs, key))
}

/// The error every reader gives for a target `who` (Doc::named) whose
/// body lacks attribute `key`: the attributes it has, and for a density,
/// what gives one.
pub fn missing(who: &str, attrs: &[(String, Value)], key: &str) -> InkError {
    let names: Vec<&str> = attrs.iter().map(|(a, _)| a.as_str()).collect();
    let has = if names.is_empty() { "none".to_string() } else { names.join(", ") };
    let give = if key == "density" {
        ". Give its body (attr material NAME) or (attr density 7.85g/cm3) in the file that defines it, or give an assembly placing it a published (attr mass M)"
    } else {
        ""
    };
    InkError { hint: Some(format!("its attrs: {has}{give}")), ..InkError::msg(format!("{who} has no attr '{key}'")) }
}

/// A density's unit is g/mm^3 everywhere (a body's attribute, the
/// library, studies, fluids), so a body's mass over its mm^3 is in grams;
/// a number may say its unit, `7.85g/cm3`, `7850kg/m3` (sexp::UNITS). The bounds
/// are a plausibility check. Unit slips inside them are indistinguishable
/// from other materials. Zero is massless.
pub const DENSITY_RANGE: (f64, f64) = (5e-8, 0.023);

/// Supported thermal conductivities, W/(m K), from insulating materials
/// through high-purity conductors. A plausibility check in the native unit.
pub const CONDUCTIVITY: crate::Magnitude = crate::Magnitude {
    name: "thermal conductivity (W/m/K)",
    range: (1e-5, 1e5),
    strict: false,
    zero: false,
    why: |v| format!("conductivity {v} is outside the supported material range; write its unit, e.g. 0.2W/m/K"),
};

/// Supported specific heats in J/(kg K). The low end includes the
/// normalized heat equation; room-temperature hydrogen is about 14,300.
pub const CAPACITY: crate::Magnitude = crate::Magnitude {
    name: "specific heat capacity (J/kg/K)",
    range: (1e-3, 2e4),
    strict: false,
    zero: false,
    why: |v| format!("capacity {v} is outside the supported specific heat range; write its unit, e.g. 4184J/kg/K"),
};

pub const DENSITY: crate::Magnitude = crate::Magnitude { name: "density (g/mm^3)", range: DENSITY_RANGE, strict: false, zero: true, why: density_message };

/// Why `d` cannot be a density in g/mm^3, if it cannot.
pub fn density_error(d: f64) -> Option<String> {
    DENSITY.error(d)
}

fn density_message(d: f64) -> String {
    let (_, hi) = DENSITY_RANGE;
    let why = if d < 0.0 {
        "is negative".to_string()
    } else if d > hi {
        format!("is denser than any material in g/mm^3, ink's density unit (osmium is 0.0226). If it is kg/m^3, write {d}kg/m3; if g/cm^3, {d}g/cm3")
    } else {
        "is lighter than any material in g/mm^3, ink's density unit (hydrogen is 9e-8). Converted twice? Write it with its unit, e.g. 7850kg/m3".to_string()
    };
    format!("density {d} {why}")
}

/// A material's density (g/mm^3), refused if it cannot be one: the
/// library's are, but `(data (density 7850))` is not.
pub fn material_density_error(m: &serde_json::Value) -> Option<String> {
    m.get("density")?.as_f64().and_then(density_error).map(|e| format!("material's {e}"))
}

/// A command's `--density`: a number (which may say its unit, `7.85g/cm3`),
/// refused as `density_error` refuses it, or None for an attribute's name.
pub fn density_option(s: &str) -> Option<Result<f64, String>> {
    if !s.starts_with(|c: char| c.is_ascii_digit() || matches!(c, '.' | '-' | '+')) {
        return None;
    }
    Some(crate::sexp::number_of(s, crate::units::Quantity::Density, "this").and_then(|d| density_error(d).map_or(Ok(d), Err)))
}

/// The density `ink mass` weighs target `who`
/// (Doc::named) by, and what it is (for the answer): the attribute asked
/// for by name (`--density NAME`), which it must have (`required`); else
/// its `density`, if it has one (None: it has no mass).
pub fn mass_density(d: &mut Doc, who: &str, attrs: &[(String, Value)], asked: Option<&str>) -> Result<Option<(NodeId, String)>, InkError> {
    let key = asked.unwrap_or("density");
    let f = match asked {
        Some(_) => required(&mut d.g, who, attrs, key)?,
        None => match attr(&mut d.g, who, attrs, key)? {
            Some(f) => f,
            None => return Ok(None),
        },
    };
    let published = d.published.iter().find(|(k, _)| d.g.var(Var::Param(*k as u32)) == f).map(|(_, what)| what.clone());
    let has = |k: &str| attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v);
    let what = match (published, has(key), has("material"), d.g.as_const(f)) {
        (Some(what), ..) if key == "density" => format!(
            "{what}. Volume is measured to {}% where it settles; known masses are refined to that fraction of the remainder, at the configuration it is published for (the params' defaults, or its mass-params, but those the mass reads). Uniform density: \
             centre of mass is the centroid, inertia assumes that uniform distribution",
            query::ACCURACY * 100.0
        ),
        (_, None, Some(m), Some(d)) => {
            let name = if let Value::Str(n) = m { n.as_str() } else { "its material" };
            if key == "density" { format!("{name}: {d} g/mm^3") } else { format!("{name}'s {key}: {d}") }
        }
        _ => format!("attr '{key}'"),
    };
    Ok(Some((f, what)))
}

/// The name to write for one the library does not take: the same name in
/// another case or with spaces or underscores for hyphens, or a name it
/// lists in `other-names`.
pub fn instead(name: &str) -> Option<&'static str> {
    let key = name.to_ascii_lowercase().replace(['_', ' '], "-");
    let lib = library();
    let own = lib.entries.iter().map(|(n, _)| (n, n));
    own.chain(lib.aliases.iter().map(|(n, _)| (n, n))).chain(lib.instead.iter().map(|(n, to)| (n, to))).find(|(n, _)| **n == key).map(|(_, to)| to.as_str())
}

/// A name that is not a material: what is wrong, and what to write.
pub fn refused(name: &str) -> (String, String) {
    match instead(name) {
        Some(to) => (format!("material '{name}' is spelled '{to}'"), format!("write {to}: see canonical grades and compatibility aliases in std/materials.ink")),
        None => (format!("unknown material '{name}'"), format!("materials: {}", names().join(" "))),
    }
}

impl Doc {
    /// A target as a reader names it: its path, and the body it places,
    /// with the file that defines it where that is another (an import).
    pub fn named(&self, it: &Item) -> String {
        match it.instance.as_ref().and_then(|(b, _)| self.body(b)) {
            Some(b) if *b.file != *self.file => format!("'{}' (body '{}' in {})", it.label, b.name, b.file),
            Some(b) if b.name != it.label => format!("'{}' (body '{}')", it.label, b.name),
            _ => format!("'{}'", it.label),
        }
    }

    /// `(attr mass M)` on a body: its published mass in grams, spread
    /// uniformly over it (PublishedMass).
    /// Built with arguments (`of` its definition) and the same mass, it
    /// takes the definition's density, measured at its defaults: a size
    /// passed as an argument scales the mass as a param would.
    pub(super) fn published_mass(&mut self, at: &Sexp, name: &str, of: Option<&str>, shape: NodeId, attrs: &mut Vec<(String, Value)>) -> R<()> {
        let Some((_, mass)) = attrs.iter().find(|(k, _)| k == "mass") else { return Ok(()) };
        let mass = self.total_mass(at, mass)?;
        if attrs.iter().any(|(k, _)| k == "density") || material_property(attrs, "density").is_some() {
            return Err(self.err_hint(
                at,
                "attr mass defines density, so it cannot accompany a density or a material with density",
                Some("for a solid with a known material, keep material density and compare its measured mass with the published mass".into()),
            ));
        }
        let defined = of.and_then(|d| self.body(d)).map(|b| {
            let get = |k: &str| b.attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v.clone());
            (get("mass"), get("density"))
        });
        if let Some((Some(Value::Field(m, _)), Some(density))) = defined
            && m == mass
        {
            attrs.push(("density".into(), density));
            return Ok(());
        }
        let config = self.published_config(at, attrs.iter().find(|(k, _)| k == "mass-params").map(|(_, v)| v))?;
        let k = self.publish(at, name, format!("published mass of '{name}' over its volume"), mass, vec![(shape, 1.0)], Vec::new(), config, None)?;
        // Worked out of a mass in g over a volume in mm3: it says so, as a density written does.
        let said = Rc::new(units::Said::worked(units::unit(Quantity::Density.unit()).expect("a quantity's unit is a unit"), None, None, &self.file, at.line));
        attrs.push(("density".into(), Value::Field(self.g.var(Var::Param(k as u32)), Some(said))));
        Ok(())
    }

    /// A body built with arguments, `of` a definition an assembly's
    /// published mass gave a density (a catalog part's root weighing its
    /// bodies), with no density, material or mass of its own: it takes that
    /// density, measured where the mass is published, so a size passed as
    /// an argument scales the mass as a param would.
    pub(super) fn weighed_as_defined(&self, mut body: Body) -> Body {
        let own = |a: &[(String, Value)]| a.iter().any(|(k, _)| k == "density" || k == "mass") || material_property(a, "density").is_some();
        let Some(def) = body.of.as_deref().and_then(|d| self.body(d)) else { return body };
        if own(&body.attrs) {
            return body;
        }
        if let Some((_, rho)) = def.attrs.iter().find(|(k, _)| k == "density") {
            let mut attrs = (*body.attrs).clone();
            attrs.push(("density".into(), rho.clone()));
            body.attrs = Attrs::new(attrs);
        }
        body
    }

    /// `(attr mass M)` on assembly `asm`, its instances read: the bodies it
    /// places, each as often as it does, weigh M together. Those with a
    /// density of their own keep it; the rest share one (PublishedMass).
    pub(super) fn assembly_mass(&mut self, at: &Sexp, asm: &str, mass: &Value, params: Option<&Value>) -> R<()> {
        let mass = self.total_mass(at, mass)?;
        let mut placed: Vec<(String, f64)> = Vec::new();
        self.count_placed(at, asm, &mut placed)?;
        let (mut shared, mut own, mut names) = (Vec::new(), Vec::new(), Vec::new());
        // A body with a density when this is read keeps it: its own, its
        // material's, or a share of an assembly inside this one (a
        // catalog part's root): this one's mass weighs the rest.
        for (b, n) in placed {
            let (shape, attrs) = self.body(&b).map(|x| (x.shape, x.attrs.clone())).unwrap();
            match attr(&mut self.g, &format!("'{b}'"), &attrs, "density").map_err(|e| self.err(at, e.message))? {
                Some(rho) => own.push((shape, rho, n)),
                None => {
                    shared.push((shape, n));
                    names.push(b);
                }
            }
        }
        // A mass that gives no body a density is not used: said, not
        // ignored (one lighter than what has its own is an error where
        // weighed: no density can make up the difference).
        if shared.is_empty() {
            let w = format!("{}:{}: the (attr mass ...) of assembly '{asm}' {UNUSED_MASS}. Remove it, or the density of the bodies it should weigh", self.file, at.line);
            if !self.load_warnings.contains(&w) {
                self.load_warnings.push(w);
            }
            self.keep_mass(asm, AssemblyMass { mass, density: None });
            return Ok(());
        }
        let config = self.published_config(at, params)?;
        let what = match own.is_empty() {
            true => format!("published mass of '{asm}' over its bodies' volume"),
            false => format!("published mass of '{asm}', less what its bodies with a density of their own weigh, over the others' volume"),
        };
        // A body sharing it that deforms (its shape reads t) keeps the mass it has where the mass is
        // published, t = 0: the density is found there, and the deforming body's is that mass over
        // its volume at each moment (KeptShare). Else the one density is found at every moment: what
        // moves rigidly keeps its volume.
        let deforms: Vec<bool> = shared.iter().map(|&(s, _)| self.g.depends_on(s, |v| v == Var::T)).collect();
        let moment = (deforms.contains(&true) && !self.g.depends_on(mass, |v| v == Var::T)).then_some(0.0);
        let k = self.publish(at, asm, what, mass, shared.clone(), own, config, moment)?;
        let f = self.g.var(Var::Param(k as u32));
        let said = Rc::new(units::Said::worked(units::unit(Quantity::Density.unit()).expect("a quantity's unit is a unit"), None, None, &self.file, at.line));
        for (b, (&(shape, _), deforming)) in names.into_iter().zip(shared.iter().zip(deforms)) {
            let rho = match (deforming, moment) {
                (true, Some(t0)) => {
                    let kept = self.params.len();
                    let name = format!("{asm}.{b}.kept-density");
                    self.add_param(Param { name: name.clone(), default: f64::NAN, tol: None });
                    let share = KeptShare { name, shape: Arc::new(Tape::new(&self.g, shape)), of: k, moment: t0, owns: [kept] };
                    let fields = vec![shape, f];
                    self.keep_form("mass", "mass", at, Arc::new(share), Reads { fields, ..Reads::default() })?;
                    self.published.push((kept, format!("the mass '{b}' has of '{asm}'s published mass at t = {t0}, over its volume at each moment")));
                    self.g.var(Var::Param(kept as u32))
                }
                _ => f,
            };
            let body = self.bodies.iter_mut().find(|x| x.name == b).unwrap();
            let mut attrs = (*body.attrs).clone();
            attrs.push(("density".into(), Value::Field(rho, Some(said.clone()))));
            body.attrs = Attrs::new(attrs);
        }
        self.keep_mass(asm, AssemblyMass { mass, density: Some(f) });
        Ok(())
    }

    fn keep_mass(&mut self, asm: &str, mass: AssemblyMass) {
        self.assemblies.iter_mut().find(|a| a.name == asm).unwrap().mass = Some(mass);
    }

    /// The bodies assembly `asm` places, through its sub-assemblies, and
    /// how many times each.
    fn count_placed(&self, at: &Sexp, asm: &str, out: &mut Vec<(String, f64)>) -> R<()> {
        for i in &self.assembly(asm).unwrap().instances {
            if self.body(&i.child).is_some() {
                match out.iter_mut().find(|(b, _)| *b == i.child) {
                    Some((_, n)) => *n = graph::exact_bin(BinOp::Add, *n, 1.0).ok_or_else(|| self.err(at, "body count exceeds exact measurement range"))?,
                    None => out.push((i.child.clone(), 1.0)),
                }
            } else {
                self.count_placed(at, &i.child, out)?;
            }
        }
        Ok(())
    }

    /// A published mass's value: grams, the same everywhere.
    fn total_mass(&self, at: &Sexp, mass: &Value) -> R<NodeId> {
        let Value::Field(mass, _) = *mass else { return Err(self.err(at, "attr mass must be a number in grams")) };
        if self.g.depends_on(mass, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
            return Err(self.err(at, "attr mass is a total mass, so it cannot vary with x, y or z"));
        }
        // A number at load is checked there; what the mass is spread over
        // is measured only where a mass is read (PublishedMass::on_demand).
        if !self.g.depends_on(mass, |v| v == Var::T) {
            let m = self.load_eval(mass);
            if !m.is_finite() || m < 0.0 {
                return Err(self.err(at, format!("attr mass must be finite and nonnegative (grams), not {m}")));
            }
        }
        Ok(mass)
    }

    /// The configuration a mass is published for: every param read so far
    /// (but those a solve sets) at its default, or at the value
    /// `(attr mass-params (read-json FILE))` gives it, FILE a JSON object
    /// of param values.
    fn published_config(&self, at: &Sexp, given: Option<&Value>) -> R<Vec<(usize, f64)>> {
        let solved: HashSet<usize> = self.solves.iter().flat_map(|s| s.owns.iter().copied()).collect();
        let mut config: Vec<(usize, f64)> = (0..self.params.len()).filter(|k| !solved.contains(k)).map(|k| (k, self.params[k].default)).collect();
        let usage = "attr mass-params is the configuration the mass is published for: (read-json FILE), FILE a JSON object of param values, e.g. {\"length\": 1193.8}";
        let Some(given) = given else { return Ok(config) };
        let Some(given) = (match given {
            Value::Data(d) => d.value().as_object().cloned(),
            _ => None,
        }) else {
            return Err(self.err_hint(at, "bad attr mass-params", Some(usage.into())));
        };
        for (name, v) in given {
            let (Some(k), Some(v)) = (self.param_index(&name), v.as_f64()) else {
                return Err(self.err_hint(at, format!("attr mass-params: '{name}' is not a param with a number"), Some(usage.into())));
            };
            config.retain(|(j, _)| *j != k);
            config.push((k, v));
        }
        Ok(config)
    }

    /// The derived density `NAME.mass-density` of a published mass over
    /// `shared`, less what `own` weighs (PublishedMass): its param.
    fn publish(
        &mut self,
        at: &Sexp,
        name: &str,
        what: String,
        mass: NodeId,
        shared: Vec<(NodeId, f64)>,
        own: Vec<(NodeId, NodeId, f64)>,
        config: Vec<(usize, f64)>,
        moment: Option<f64>,
    ) -> R<usize> {
        let k = self.params.len();
        let derived = format!("{name}.mass-density");
        // Not a number until measured: a reader that misses Doc::want
        // weighs nothing it could take for a mass (0 is a density).
        self.add_param(Param { name: derived.clone(), default: f64::NAN, tol: None });
        let reads: Vec<NodeId> = std::iter::once(mass).chain(shared.iter().map(|s| s.0)).chain(own.iter().flat_map(|o| [o.0, o.1])).collect();
        let uses_t = moment.is_none() && reads.iter().any(|&r| self.g.depends_on(r, |v| v == Var::T));
        let tape = |g: &Graph, n: NodeId| Arc::new(Tape::new(g, n));
        // The params the mass reads keep their values: it is published at
        // each (a toleranced mass, a mass per length).
        let follows = self.params_read(&[mass]);
        let mut shared_counts = std::collections::BTreeMap::<NodeId, f64>::new();
        for &(shape, n) in &shared {
            let count = shared_counts.entry(shape).or_default();
            *count = graph::exact_bin(BinOp::Add, *count, n).ok_or_else(|| self.err(at, "body count exceeds exact measurement range"))?;
        }
        let p = PublishedMass {
            name: derived,
            shared_shapes: shared_counts.into_iter().collect(),
            shared: shared.iter().map(|&(s, n)| (tape(&self.g, s), n)).collect(),
            own: own.iter().map(|&(s, rho, n)| KnownBody::new(&mut self.g, s, rho, n)).collect(),
            mass: Tape::new(&self.g, mass),
            at: config.into_iter().filter(|(k, _)| !follows.contains(k)).collect(),
            moment,
            owns: [k],
            uses_t,
        };
        self.keep_form("mass", "mass", at, Arc::new(p), Reads { fields: reads, ..Reads::default() })?;
        self.published.push((k, what));
        Ok(k)
    }

    /// A material: a name from the library (symbol or string), or data with
    /// the property keys.
    pub(super) fn material_of(&mut self, m: &Sexp, env: &mut Env) -> R<serde_json::Value> {
        let name = match &m.k {
            SK::Sym(s) | SK::Str(s) => Some(s.clone()),
            _ => None,
        };
        if let Some(n) = &name
            && let Some(v) = lookup(n)
        {
            return Ok(v);
        }
        let refuse = |d: &Self, n: &str| {
            let (what, hint) = refused(n);
            d.err_hint(m, what, Some(hint))
        };
        if let Some(n) = name.as_deref().filter(|n| self.def(n).is_none() && env.iter().all(|(k, _)| k != n)) {
            return Err(refuse(self, n));
        }
        match self.eval(m, env)? {
            Value::Data(d) if d.value().is_object() => match material_density_error(d.value()) {
                Some(e) => Err(self.err(m, e)),
                None => Ok(d.value().clone()),
            },
            Value::Str(s) => lookup(&s).ok_or_else(|| refuse(self, &s)),
            _ => Err(self.err_hint(m, "a material is a library name or (material name)", Some(format!("materials: {}", names().join(" "))))),
        }
    }

    /// `(material name)` and `(by-body target key [default])`.
    pub(super) fn eval_materials(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        match name {
            "material" => {
                let [m] = args else { return Err(self.err_hint(e, "bad material", Some(format!("usage: (material name); names: {}", names().join(" "))))) };
                return Ok(Value::Data(crate::Data::new(self.material_of(m, env)?)));
            }
            "by-body" => {
                // A field taking each body's property where that body is:
                // (by-body target key [default]), the body's attribute `key`.
                let usage = "usage: (by-body target key [default]): a field with each body's property where it is (attr key, else its material's key)";
                let (target, key, default) = match args {
                    [t, k] => (t, k, None),
                    [t, k, d] => {
                        let v = self.eval(d, env)?;
                        let Value::Field(f, said) = v else { return Err(self.err(d, "by-body's default must be a number or a field")) };
                        (t, k, Some((f, said)))
                    }
                    _ => return Err(self.err_hint(e, "bad by-body", Some(usage.into()))),
                };
                let target = self.name_of(target)?;
                let key = self.name_of(key)?;
                let (_, items) = self.scoped_items(&target).map_err(|m| self.err(e, m))?;
                // Each body's own value, else the default; the nearest body
                // (smallest field value) wins, inside and out: solvers
                // sample properties slightly past surfaces, where "inside
                // this body, else the rest" would give the wrong
                // neighbour's value.
                let mut acc: Option<(Value, NodeId)> = None;
                for it in items {
                    let who = self.named(&it);
                    let (shape, attrs) = (self.placed_field(&it), self.placed_attrs(&it));
                    let attr = attr_value(&mut self.g, &who, &attrs, &key).map_err(|m| self.err_hint(e, m.message, m.hint.or(Some(usage.into()))))?;
                    // A body's own number is bare (a material's is in ink's unit): a default said in a
                    // unit cannot sit beside it, and the unit is given to the whole field instead.
                    if let (Some((_, None)), Some((d, Some(said)))) = (&attr, &default)
                        && said.unit.quantity != Quantity::Other
                    {
                        let n = self.g.as_const(*d).map_or("DEFAULT".to_string(), |v| v.to_string());
                        let whole = said.unit.quantity.say(&format!("(by-body {target} {key} {n})"));
                        return Err(self.err_hint(
                            e,
                            format!("by-body gives {who}'s '{key}' as a bare number, its own attr's or its material's in ink's unit, so its default is bare too, and the field takes its unit whole: {whole}"),
                            Some(usage.into()),
                        ));
                    }
                    let (value, said) = attr.or_else(|| default.clone()).ok_or_else(|| {
                        let m = missing(&who, &attrs, &key);
                        let message = format!("{}, nor a material that gives one: give it one, or by-body a default for such bodies, (by-body {target} {key} DEFAULT)", m.message);
                        self.err_hint(e, message, m.hint.or(Some(usage.into())))
                    })?;
                    let v = Value::Field(value, said);
                    acc = Some(match acc {
                        None => (v, shape),
                        Some((a, near)) => {
                            let c = self.g.sub(shape, near);
                            (self.apply_builtin("select", vec![Value::Field(c, None), v, a], e)?, self.g.min(near, shape))
                        }
                    });
                }
                let Some((v, _)) = acc else { return Err(self.err(e, format!("'{target}' has no bodies"))) };
                if let Value::Field(f, _) = v {
                    self.want(&[f]);
                }
                return Ok(v);
            }
            _ => {}
        }
        unreachable!()
    }
}

#[cfg(test)]
mod evidence_tests {
    use super::*;

    #[test]
    fn every_scalar_has_a_source_or_an_explicit_gap() {
        let manifest: serde_json::Value = serde_json::from_str(include_str!("../../../std/material-sources.json")).unwrap();
        let materials = manifest["materials"].as_object().unwrap();
        assert_eq!(materials.len(), library().entries.len());
        for (name, value) in &library().entries {
            let record = &materials[name];
            assert_eq!(record["condition"], value["condition"], "{name}");
            let properties = record["properties"].as_object().unwrap();
            let scalars: Vec<_> = value.as_object().unwrap().iter().filter(|(_, v)| v.is_number()).collect();
            assert_eq!(scalars.len(), properties.len(), "{name}");
            for (key, _) in scalars {
                let evidence = &properties[key];
                let kind = evidence["kind"].as_str().unwrap();
                if kind == "unsourced" {
                    assert!(evidence["source"].is_null(), "{name}.{key}");
                    assert!(!["yield", "tensile", "compressive", "bending"].contains(&key.as_str()), "unsupported strength: {name}.{key}");
                } else {
                    let source = &manifest["sources"][evidence["source"].as_str().unwrap()];
                    for field in ["title", "maker", "revision", "url", "retrieved", "sha256"] {
                        assert!(!source[field].as_str().unwrap().is_empty(), "{name}.{key}: {field}");
                    }
                }
                if ["yield", "tensile", "compressive", "bending"].contains(&key.as_str()) {
                    assert!(value["strength-kind"].as_str().unwrap().contains(&format!("{key}: {kind}")), "{name}.{key}");
                }
            }
        }
        for (alias, canonical) in &library().aliases {
            let mut a = lookup(alias).unwrap();
            let c = lookup(canonical).unwrap();
            a["name"] = c["name"].clone();
            assert_eq!(a, c, "{alias}");
            assert!(materials[canonical]["aliases"].as_array().unwrap().contains(&serde_json::json!(alias)));
        }
    }

    #[test]
    fn unsupported_strengths_are_not_silently_used_as_yield() {
        for (material, property) in [
            ("ptfe", "yield"),
            ("copper", "yield"),
            ("polystyrene", "tensile"),
            ("tungsten-carbide", "bending"),
            ("tungsten-carbide", "compressive"),
            ("glass", "tensile"),
            ("borosilicate-glass", "tensile"),
        ] {
            assert!(lookup(material).unwrap().get(property).is_none(), "{material}.{property}");
        }
        assert_eq!(lookup("glass").unwrap()["bending"].as_f64(), Some(45.0));
        assert_eq!(lookup("borosilicate-glass").unwrap()["bending"].as_f64(), Some(25.0));
    }
}
