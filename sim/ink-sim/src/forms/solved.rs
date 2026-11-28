//! A solved study: its results, its error in time, and the solve.
use super::*;

/// A material property of a study: a field, or one value.
#[derive(Clone)]
pub(super) enum Stiffness {
    Field(std::sync::Arc<crate::tape::Tape>, Ctx),
    Value(f64),
}

/// A solved study: its results as fields, its metrics, and the shape
/// derivative of its compliance.
pub struct SolvedStudy {
    pub sol: crate::Solution,
    pub physics: StudyPhysics,
    /// MPa, of an elastic study's material where it is one number.
    pub yield_strength: Option<f64>,
    /// The yield strength as a field (MPa), read at the stress's points.
    pub(super) yield_field: Option<(std::sync::Arc<crate::tape::Tape>, Ctx)>,
    /// E (or k) and nu, each a field or one value (heat's nu 0): by which a
    /// cell is of one material or not, and its stiffest (`mixed`).
    pub(super) stiffness: [Stiffness; 2],
    /// A fix holds a plane of symmetry: a modal study's modes are those
    /// symmetric about it only.
    pub(super) symmetry: Vec<String>,
    pub(super) alternate_low: Option<f64>,
    pub(super) mode_field: bool,
    /// The study's shape, and every other field it depends on (its
    /// recorded inputs: materials, densities, regions, loads).
    pub(super) shape: NodeId,
    /// The shape as its solve read it, where a peak's surface is measured
    /// (`peak_surface`).
    pub(crate) shape_field: (std::sync::Arc<crate::tape::Tape>, Ctx),
    /// Each peak's surface as measured (`peak_surface`): every claim on the
    /// solve asks for the same peak.
    pub(super) surfaces: std::sync::Mutex<Vec<(V3, Option<(f64, f64)>)>>,
    /// A sharp inside corner on the surface its stress is read on, where
    /// one is (`Solution::crease`), found when first asked.
    pub(super) crease: std::sync::OnceLock<Option<V3>>,
    /// The cells of the claim's own grid it was solved for (`split_as`),
    /// else its own.
    pub(super) judge: f64,
    pub(super) inputs: Vec<NodeId>,
    /// A fix moves its region: the compliance, the loads' work, leaves out
    /// the move's (a pressed part reads 0 however stiff), so it has no
    /// value, and its shape derivative is not the boundary formula
    /// (`adjoint`).
    pub moved: bool,
    /// Its fixes' names, in the order the solution indexes them.
    pub fix_names: Vec<String>,
    /// The support margin it was solved with (mm, at least 2 cells):
    /// stress this close to a fix is reported apart.
    pub margin: f64,
}

/// Why a study a fix moves has no compliance.
/// A transient's error in time: its peak against the same run on steps
/// twice and four times as long (`Report::time_levels`), as the claims'
/// grids are compared in space (`claims::grid_error`). Backward Euler is
/// first order in its step, so each difference is first scaled to what
/// it would be on twice the step (a count that halves unevenly, 25 to
/// 13, steps less than twice as long). An estimate, not a bound; none where
/// three counts do not converge (each difference not smaller than the
/// last, the same way): steps that all smear a response faster
/// than they are read (a pulse into a small region of a conductive part,
/// shorter than every step) read it low alike, their peaks growing as the
/// steps shrink, and no error can be read from them (unbounded: a claim on
/// it is unknown, saying to take more steps).
pub struct TimeError {
    pub error: f64,
    /// How it was found (`grid_error`'s words, of step counts).
    pub how: String,
    /// The run's own steps.
    pub steps: usize,
    pub basis: String,
}

/// `TimeError` of a transient's levels, None for no transient. One step
/// (none to compare with) has an unbounded error.
pub fn time_error(levels: &[(usize, f64, f64)], residual: f64) -> Option<TimeError> {
    let sig = crate::lang::claims::sig;
    let &(n, x, _) = levels.last()?;
    // Each run's longest step: a peak at a jump is read by a step that long.
    let h = |k: usize| levels[k].2;
    let m = levels.len();
    let (error, how) = match m {
        1 => (f64::INFINITY, "one step, nothing to compare it with"),
        _ => {
            // Scaled as first order: x2 the value on twice the step,
            // x4 on four times.
            let x2 = x + (levels[m - 2].1 - x) * h(m - 1) / (h(m - 2) - h(m - 1));
            let x4 = (m > 2).then(|| x2 + (levels[0].1 - levels[1].1) * 2.0 * h(m - 1) / (h(0) - h(1)));
            // Peaks apart by less than the solves resolve (their relative residual) are not apart.
            crate::lang::claims::grid_error(x, x2, x4, residual.max(f64::EPSILON) * x.abs())
        }
    };
    let how = how.replace("grids", "step counts");
    let others: Vec<String> = levels[..m - 1].iter().rev().map(|(k, v, _)| format!("{} on {k} steps", sig(*v))).collect();
    let error_is = match error.is_finite() {
        true => format!("its error in time about {} ({how}), an estimate of what its steps miss, not a bound", sig(error)),
        false => format!("no error in time can be read from them ({how}): its steps may all miss what they are too long to follow"),
    };
    let basis = format!("its peak over time, {} on {n} steps, against {}: {error_is}", sig(x), if others.is_empty() { "no other".to_string() } else { others.join(" and ") });
    Some(TimeError { error, how, steps: n, basis })
}

pub const MOVED_COMPLIANCE: &str = "a fix moves its region: compliance, the loads' work, leaves out the move's (a part pressed by a move reads 0 however stiff it is); read the fix's reaction for the force the move takes";

/// Above this rotation (5 degrees) a linear study's kinematics are off by
/// more than the usual discretisation error: linear kinematics drop the
/// motion that goes with a rotation theta (a bent beam's shortening, about
/// theta / 3 of its deflection) and the stiffening that comes with it.
pub(super) const LINEAR_ROTATION_LIMIT: f64 = 5.0 * std::f64::consts::PI / 180.0;

/// The fewest of its cells a concave surface's radius must span for a
/// peak stress on it to be read: on fewer, its cut cells read the peak
/// only in part, and three grids agree on the miss (a bore of radius 1 on
/// 1 mm cells read 0.59 of its peak, its grids within 2%; past their
/// error, 69% over). The loosest at which every reference with a known
/// answer stays within the read's measured error (`stress_read`): its
/// 11 reads at 1.0 cells missed by up to 69%, and every one at 1.5 and
/// more (a 1.5 mm fillet on 1 mm cells first) by at most 3.4%; a hair
/// under 1.5, so a radius of exactly 1.5 cells (a 3 mm fillet on 2 mm
/// cells, read to 0.1%) passes as those references did.
pub const RESOLVED_RADIUS_CELLS: f64 = 1.45;

impl SolvedStudy {
    /// A static elastic study's metric that is a size of its fields, which are linear in its loads
    /// (`Study::linear_in_worst`): its value and where. None: `metric` is not one.
    pub(super) fn linear_value(&self, metric: &str) -> Option<Result<(f64, Option<V3>), String>> {
        let peak = |e: Result<crate::Extreme, String>| e.map(|e| (e.value, Some(e.at)));
        Some(match metric {
            "max-displacement" => peak(self.sol.max_nodal()),
            "max-stress" => peak(self.sol.max_cell_scalar()),
            // With no cell of one material away from the supports, no value:
            // the bound alone (`bound_band`).
            "min-safety-factor" => peak(self.min_safety_factor()).or_else(|why| self.safety_factor_between_reads(false).map(|b| (b.value, Some(b.at))).ok_or(why)),
            _ => return None,
        })
    }

    /// Itself again, to set its solution's fields to other loads' (`Study::of_sums`): what it
    /// found of a peak's surface is found again.
    pub(super) fn copy(&self) -> SolvedStudy {
        SolvedStudy {
            sol: self.sol.clone(),
            physics: self.physics,
            yield_strength: self.yield_strength,
            yield_field: self.yield_field.clone(),
            stiffness: self.stiffness.clone(),
            symmetry: self.symmetry.clone(),
            alternate_low: self.alternate_low,
            mode_field: self.mode_field,
            shape: self.shape,
            shape_field: self.shape_field.clone(),
            surfaces: Default::default(),
            crease: self.crease.clone(),
            judge: self.judge,
            inputs: self.inputs.clone(),
            moved: self.moved,
            fix_names: self.fix_names.clone(),
            margin: self.margin,
        }
    }

    /// A partial solve can be inspected, but cannot answer for the whole part.
    pub(crate) fn dropped_material(&self) -> Option<String> {
        let r = &self.sol.report;
        (r.dropped_cells > 0).then(|| {
            let at = r.dropped_at.map_or(String::new(), |p| format!(", including a cell near {p:?} mm"));
            let action = if matches!(self.physics, StudyPhysics::Heat { .. }) {
                "join it to the cooled part or give it a temperature fix or convection surface"
            } else {
                "hold it or join it to the supported part"
            };
            format!("{} mm^3 of material is not connected to any fix and was dropped{at}; {action}", ink_lang::command::measured(r.dropped_volume))
        })
    }

    pub(super) fn mode_hint(&self) -> &'static str {
        if self.sol.report.possible_cuts.is_empty() && self.sol.report.unstated_holds.is_empty() {
            "solve the whole part to include its antisymmetric vibration modes"
        } else {
            "for a symmetry cut declare (mirror name x coordinate) with its axis and position; for a physical restraint add (support) to the fix"
        }
    }

    pub(crate) fn hidden_modes_warning(&self) -> Option<String> {
        let r = &self.sol.report;
        if !r.possible_cuts.is_empty() {
            let names = r.possible_cuts.iter().map(|i| self.fix_names[i.0].as_str()).collect::<Vec<_>>().join(", ");
            Some(format!("a possible lower mode if the held face ({names}) is a symmetry cut. {}", self.mode_hint()))
        } else if !r.unstated_holds.is_empty() {
            let names = r.unstated_holds.iter().map(|&i| self.fix_names[i].as_str()).collect::<Vec<_>>().join(", ");
            Some(format!(
                "a possible lower mode if the hold ({names}), which holds one direction only and does not say what it is, is a cut of symmetry and no physical support. {}",
                self.mode_hint()
            ))
        } else if !self.symmetry.is_empty() && !r.modes.is_empty() {
            Some(format!("the vibration modes antisymmetric about the plane of symmetry '{}'. {}", self.symmetry.join(", "), self.mode_hint()))
        } else {
            None
        }
    }

    /// Where a peak stress at `at` lies on a surface: that surface's
    /// tightest concave radius within a cell and a half (0 at a crease), and
    /// the size of the leaf cell the peak is read in, `(radius, cell)`
    /// (`surface_radius`, its steps a quarter to a sixteenth of the cell).
    /// None where it lies inside the material, beyond a cell and a half of
    /// any surface (an interface's peak, a section's). A claim's other grid
    /// (`split_as`) surveys what the claim's own would, cells there `judge`
    /// over its own as big: a coarser grid looking further found a sharp
    /// corner its claim's grid did not, and sent the claim to grids 64
    /// times finer. Its radius is judged against its own cells.
    pub fn peak_surface(&self, at: V3) -> Option<(f64, f64)> {
        if let Some(&(_, s)) = self.surfaces.lock().expect("surfaces").iter().find(|(p, _)| *p == at) {
            return s;
        }
        let h = self.sol.cell_size_at(at).unwrap_or(self.sol.h);
        let seen = h * (self.judge / self.sol.h).min(1.0);
        let (t, ctx) = &self.shape_field;
        let s = crate::surface_radius(t, ctx, at, 1.5 * seen, seen / 8.0).map(|(_, r)| (r, h));
        self.surfaces.lock().expect("surfaces").push((at, s));
        s
    }

    /// Solve study `name` (cached by content, see `sim::run_study`).
    pub fn run(doc: &mut Doc, ctx: &Ctx, name: &str, opts: &RunOpts) -> R<SolvedStudy> {
        let (mut setup, physics, ys) = crate::study_setup(doc, ctx, name, opts.res, opts.tol)?;
        (setup.split_as, setup.both_ways) = (opts.split_as, opts.both_ways);
        let mut sol = crate::solve_setup(doc, ctx, name, &setup)?;
        // There is no geometric test for the intent of a zero-displacement
        // roller on a flat cut. Check all admissible interpretations of the
        // legacy syntax. A physical layer has thickness, so each change to
        // a plane trace needs its own prestress. Declared mirrors need only
        // the one static solution. Parallel cuts do not define a finite part.
        let candidates = &sol.report.possible_cuts;
        let mut alternate_low: Option<f64> = None;
        if !candidates.is_empty() && candidates.len() <= 3 && sol.report.buckle.is_some() {
            let mut axes = [false; 3];
            for f in &setup.fixes {
                if let crate::Boundary::Mirror { axis, .. } = f.boundary {
                    axes[axis] = true;
                }
            }
            let mut low = f64::INFINITY;
            let mut admissible = false;
            for mask in 1..(1 << candidates.len()) {
                let mut selected_axes = axes;
                let mut parallel = false;
                for (bit, &(i, axis, at)) in candidates.iter().enumerate() {
                    setup.fixes[i].boundary = if mask & (1 << bit) != 0 {
                        parallel |= std::mem::replace(&mut selected_axes[axis], true);
                        crate::Boundary::Mirror { axis, at }
                    } else {
                        crate::Boundary::Support
                    };
                }
                // Two parallel cuts cannot both generate a finite part.
                // Either one alone can still be a cut with a physical roller
                // at the other end. Check both of those interpretations.
                if parallel {
                    continue;
                }
                let factor = match crate::solve_cached_checked(&setup, ctx) {
                    Err(crate::SolveError::InvalidMirror(_)) => continue,
                    Ok(s) if s.valid().is_ok() && s.report.incomplete.is_empty() && s.report.dropped_cells == 0 => s.report.buckle.and_then(|b| b.answer()),
                    _ => None,
                };
                admissible = true;
                let Some((sector_low, _)) = factor else {
                    low = 0.0;
                    break;
                };
                low = low.min(sector_low);
            }
            if admissible {
                alternate_low = Some(low);
            } else {
                sol.report.possible_cuts.clear();
            }
        }
        let yield_strength = ys.and_then(|f| doc.g.as_const(f));
        let yield_field = ys.map(|f| {
            let f = doc.g.untagged(f);
            (std::sync::Arc::new(crate::tape::Tape::new(&doc.g, f)), ctx.clone())
        });
        let st = doc.form::<Study>(name).expect("a solved study exists");
        let shape = st.shape;
        let inputs: Vec<NodeId> = crate::study_inputs(doc, name).into_iter().filter(|&f| f != shape).collect();
        let fix_names: Vec<String> = st.fixes.iter().map(|f| f.name.clone()).collect();
        let margin = st.support_margin.unwrap_or(crate::SUPPORT_MARGIN_MM).max(2.0 * sol.h);
        let symmetry = st.fixes.iter().filter(|f| matches!(f.boundary, crate::Boundary::Mirror { .. })).map(|f| f.name.clone()).collect();
        let stiffness = [
            setup.coef_const.map_or_else(|| Stiffness::Field(setup.coef.clone(), ctx.clone()), Stiffness::Value),
            match (&setup.nu, setup.physics) {
                (Some(t), _) => Stiffness::Field(t.clone(), ctx.clone()),
                (None, crate::Physics::Elastic { nu }) => Stiffness::Value(nu),
                (None, crate::Physics::Heat) => Stiffness::Value(0.0),
            },
        ];
        let shape_field = (setup.shape.clone(), ctx.clone());
        let judge = opts.split_as.unwrap_or(sol.h);
        Ok(SolvedStudy {
            sol,
            physics,
            yield_strength,
            yield_field,
            stiffness,
            symmetry,
            alternate_low,
            mode_field: false,
            shape,
            shape_field,
            surfaces: Default::default(),
            crease: Default::default(),
            judge,
            inputs,
            moved: st.moved,
            fix_names,
            margin,
        })
    }
}
