//! What a solved study reads of itself: why a metric is unbounded, its
//! bands, its materials' mix, safety factors and warnings.
use super::*;

impl SolvedStudy {
    /// The solution its linear buckling was solved on (its modes the
    /// buckling modes, `report.buckle` the factors): its own, which runs
    /// it on every static elastic solve that converged.
    /// Why `metric`'s value may grow without bound as the cells shrink, or None where nothing in this solve
    /// could make it (`Measured::unbounded`). A load or hold less than a cell across, its claim's cells: a
    /// point's or a line's displacement, temperature and stress have no bound (100 N on a 0.1 mm patch read
    /// 1.022e-4, 3.262e-4, 7.736e-4, then 6.336e-3 on cells of 0.125 mm). A peak stress at a hold's edge or a
    /// sharp corner (its exclusions, singular), where materials meet, or on a curve finer than its cells.
    /// A buckling factor or a frequency: its mode's reach (`Measured::extent`). Any other: not classified.
    pub(super) fn unbounded(&self, metric: &str, excluded: &[crate::lang::Excluded], at: Option<V3>) -> Option<String> {
        let sol = &self.sol;
        // Under a cell across, as the cells of the claim's own grid (`judge`) see it: a load landing in one
        // cell, or a face of the part less than one wide (a point's spread is 0, a cell's width over the root
        // of 12).
        let across = self.judge / 12f64.sqrt();
        let narrow = |s: &V3| 3 - s.iter().filter(|v| **v >= across).count();
        let concentrated = || -> Option<String> {
            let mm = crate::lang::claims::sig(self.judge);
            sol.report
                .load_spread
                .iter()
                .find(|(_, s)| narrow(s) >= 2)
                .map(|(name, s)| format!("the load '{name}' is applied over less than a cell across in {} directions (cells of {mm} mm), where nothing bounds it", narrow(s)))
                .or_else(|| {
                    sol.report.surface_measures.iter().filter(|m| m.kind == "fix" || m.kind == "load").find(|m| m.covered > 0.0 && narrow(&m.spread) >= 2).map(|m| {
                        format!("the {} '{}' acts over less than a cell across in {} directions (cells of {mm} mm), where nothing bounds it", m.kind, m.name, narrow(&m.spread))
                    })
                })
        };
        match metric {
            "max-displacement" | "compliance" | "max-temperature" => concentrated(),
            "min-buckling-factor" | "min-frequency" => None,
            "max-stress" | "min-safety-factor" => excluded
                .iter()
                .find(|e| e.singular)
                .map(|e| e.what.clone())
                .or_else(concentrated)
                .or_else(|| {
                    matches!(self.stiffness[0], Stiffness::Field(..))
                        .then(|| "its material varies, and its peak may sit where materials meet a free surface, where the stress has no bound".to_string())
                })
                .or_else(|| {
                    let cells = at.and_then(|p| self.peak_surface(p)).map_or(f64::INFINITY, |(r, h)| r / h);
                    (cells < RESOLVED_RADIUS_CELLS).then(|| format!("its peak sits on a curve under {RESOLVED_RADIUS_CELLS} of its cells across, read as a corner"))
                }),
            _ => Some(crate::lang::NOT_CLASSIFIED.to_string()),
        }
    }

    pub fn buckling(&self) -> Result<&crate::Solution, String> {
        if self.sol.report.buckle.is_some() {
            return Ok(&self.sol);
        }
        Err("no linear buckling: it is run on a static elastic study's solve, and this one is not that, or did not converge".into())
    }

    /// `min-buckling-factor`'s band: the factor as found (the Rayleigh
    /// quotient of its mode, at least the lowest factor) down to an
    /// estimate of what the eigensolve's residual allows, taking the mode
    /// found as the lowest. A claim decides on it, widened by its grid
    /// error.
    pub(super) fn buckling_band(&self, metric: &str) -> Option<crate::lang::Band> {
        if metric != "min-buckling-factor" {
            return None;
        }
        let (low, high) = self.buckling().ok()?.report.buckle.as_ref()?.answer()?;
        // No mode found: what its search excluded, and up (none at all,
        // where it excluded every factor: its value, infinite, alone).
        if high.is_infinite() {
            return low.is_finite().then(|| crate::lang::Band {
                low,
                high,
                parts: vec![("eigensolve search".into(), f64::INFINITY)],
                basis: format!(
                    "its linear buckling eigensolve found no direction its loads compress: no factor below {}, an estimate of what its residual allows (measured through a multigrid stand-in for the stiffness's inverse, the mode found taken as the lowest), and none searched above it",
                    crate::lang::claims::sig(low)
                ),
                narrows: Vec::new(),
            });
        }
        Some(crate::lang::Band {
            low,
            high,
            parts: vec![("eigensolve".into(), 0.5 * (high - low))],
            basis: format!(
                "its linear buckling eigensolve: its mode's Rayleigh quotient, at least the lowest factor, and {}% below it, an estimate of what its residual allows (measured through a multigrid stand-in for the stiffness's inverse, the mode found taken as the lowest)",
                100.0 * (1.0 - low / high).max(0.0)
            ),
            narrows: Vec::new(),
        })
    }

    /// `min-safety-factor`'s band where no cell away from the supports is
    /// of one material (`mixed`): no value, so from its bound up. A bound is
    /// the least the factor can be, so a claim it holds holds, value or not;
    /// none is failed but by a value, or by a bound past a claim's upper
    /// limit.
    pub(super) fn bound_band(&self, metric: &str) -> Option<crate::lang::Band> {
        if metric != "min-safety-factor" || self.min_safety_factor().is_ok() {
            return None;
        }
        let b = self.safety_factor_between_reads(false)?;
        Some(crate::lang::Band {
            low: b.value,
            high: f64::INFINITY,
            parts: vec![("cells not of one material".into(), f64::INFINITY)],
            basis: "no cell away from the supports is of one material (E, nu or yield one value across it), so none gives a value: the least the factor can be, each such cell's least yield over the stress's bound over it, and up".into(),
            narrows: Vec::new(),
        })
    }

    /// A transient's peak temperature's band: its error in time
    /// (`time_error`), on each side.
    pub(super) fn time_band(&self, metric: &str) -> Option<crate::lang::Band> {
        if metric != "max-temperature" {
            return None;
        }
        let (_, x, _) = self.sol.report.peak?;
        let t = time_error(&self.sol.report.time_levels, self.sol.report.residual)?;
        Some(crate::lang::Band {
            low: x - t.error,
            high: x + t.error,
            parts: vec![("time steps".into(), t.error)],
            basis: t.basis,
            narrows: vec![("time steps".into(), format!("more (steps n) in its (transient ...): {} or more", 2 * t.steps))],
        })
    }

    /// The yield strength at points (MPa): the study's field, NaN without one.
    pub(super) fn yield_at(&self, pts: &[V3]) -> Vec<f64> {
        match &self.yield_field {
            Some((t, ctx)) => {
                let mut v = Vec::new();
                crate::eval::eval_points(t, ctx, pts, &mut v);
                v
            }
            None => vec![f64::NAN; pts.len()],
        }
    }

    /// The cells, as boxes, that are not of one material: their E, nu or
    /// yield not one value across the box, by the fields' bounds (a strip, a
    /// gradient, an interface between bodies' materials). For each, None
    /// where it is of one material, else the least yield in it and the
    /// stiffest it can be, E / (1 + nu) (NaN where unbounded). No point
    /// read in such a cell is its material's: it is bounded only. These are
    /// the cells a finer grid must split to read.
    pub fn mixed(&self, boxes: &[(V3, V3)]) -> Vec<Option<crate::MixedCell>> {
        let range = |(t, ctx): &(std::sync::Arc<crate::tape::Tape>, Ctx), (lo, hi): &(V3, V3)| {
            crate::eval::eval_interval(t, ctx, std::array::from_fn(|a| ink_core::interval::I::new(lo[a], hi[a])), None)
        };
        let of = |p: &Stiffness, b: &(V3, V3)| match p {
            Stiffness::Field(t, ctx) => range(&(t.clone(), ctx.clone()), b),
            Stiffness::Value(v) => ink_core::interval::I::point(*v),
        };
        let varies = |r: ink_core::interval::I| r.lo.is_nan() || r.hi.is_nan() || r.hi > r.lo;
        let heat = matches!(self.physics, StudyPhysics::Heat { .. });
        boxes
            .iter()
            .map(|b| {
                let y = self.yield_field.as_ref().map(|f| range(f, b));
                let [e, nu] = [of(&self.stiffness[0], b), of(&self.stiffness[1], b)];
                let mixed = varies(e) || varies(nu) || y.is_some_and(varies);
                // Its stiffest over the box shrunk by a millionth of its
                // side: a face on its wall (a cell beside a band's face, its
                // reads reaching it) is the material beyond's, not its own.
                let shrunk: (V3, V3) = (std::array::from_fn(|a| b.0[a] + 1e-6 * (b.1[a] - b.0[a])), std::array::from_fn(|a| b.1[a] - 1e-6 * (b.1[a] - b.0[a])));
                let [e, nu] = [of(&self.stiffness[0], &shrunk), of(&self.stiffness[1], &shrunk)];
                mixed.then(|| crate::MixedCell {
                    least: self.yield_strength.or(y.map(|r| r.lo)).unwrap_or(f64::NAN),
                    stiffest: if heat {
                        e.hi
                    } else if nu.lo > -1.0 {
                        e.hi / (1.0 + nu.lo)
                    } else {
                        f64::NAN
                    },
                })
            })
            .collect()
    }

    /// Where a cell is not of one material (`mixed`; `Solution::
    /// ratio_between_reads`), away from the supports (`near` false) or
    /// within their margins: the least safety factor such cells can have,
    /// their least yield over their stress's bound, and where.
    pub fn safety_factor_between_reads(&self, near: bool) -> Option<crate::Extreme> {
        if self.yield_strength.is_none() && self.yield_field.is_none() {
            return None;
        }
        self.sol.ratio_between_reads(near, &|b| self.mixed(b)).map(|e| crate::Extreme { value: 1.0 / e.value, at: e.at })
    }

    /// The least safety factor away from the supports, the yield strength
    /// over the stress point by point (each material's own where it is),
    /// and where.
    pub fn min_safety_factor(&self) -> Result<crate::Extreme, String> {
        self.sol.max_cell_ratio(&|p| self.yield_at(p), &|b| self.mixed(b)).map(|e| crate::Extreme { value: 1.0 / e.value, at: e.at })
    }

    /// How far the solve's own field rises over a peak stress's read (or
    /// one over a safety factor's) between the surface's reads
    /// (`climbed_scalar`), as a ratio: 1 where it rises no further.
    pub fn climbed(&self, metric: &str) -> f64 {
        let ratio = match metric {
            "max-stress" => self.sol.max_cell_scalar().ok().zip(self.sol.climbed_scalar()).map(|(r, c)| c.value / r.value),
            "min-safety-factor" => self
                .sol
                .max_cell_ratio(&|p| self.yield_at(p), &|b| self.mixed(b))
                .ok()
                .zip(self.sol.climbed_ratio(&|p| self.yield_at(p), &|b| self.mixed(b)))
                .map(|(r, c)| c.value / r.value),
            _ => None,
        };
        ratio.filter(|r| r.is_finite() && *r > 1.0).unwrap_or(1.0)
    }

    /// `min_safety_factor` within the support margins.
    pub fn min_safety_factor_near_supports(&self) -> Option<crate::Extreme> {
        self.sol.max_near_supports_ratio(&|p| self.yield_at(p), &|b| self.mixed(b)).map(|e| crate::Extreme { value: 1.0 / e.value, at: e.at })
    }

    /// The largest rotation of a static elastic solution (mode shapes have
    /// no scale).
    pub fn max_rotation(&self) -> Option<crate::Extreme> {
        let linear = matches!(self.physics, StudyPhysics::Elastic { .. }) && self.sol.report.modes.is_empty();
        linear.then(|| self.sol.max_rotation()).flatten()
    }

    /// Where an elastic study's peak stress sits at the edge of a fix's
    /// support margin: what it is beside, how far (the nearest left-out
    /// cell's centre), and why that matters.
    pub fn peak_beside_support(&self, peak: &crate::Extreme) -> Option<(String, f64)> {
        let s = self.sol.nearest_support(peak.at).filter(|s| s.beside)?;
        Some((self.fix_names.get(s.fix).cloned().unwrap_or_default(), s.distance))
    }

    /// The warning a stress metric carries where its peak is at a support
    /// margin's edge.
    pub fn beside_support_warning(&self, peak: &crate::Extreme) -> Option<String> {
        let (fix, d) = self.peak_beside_support(peak)?;
        let n2 = |v: f64| (v * 1000.0).round() / 1000.0;
        Some(format!(
            "the peak stress, {} MPa at {:?}, is at the edge of the excluded zone around fix '{fix}' ({} mm from its nearest left-out cell's centre): \
             the fix holds its region rigidly, and the stress where that region ends on the free surface is singular: it may reach past the support margin, and grows as res shrinks. \
             Check it on finer cells, widen (support-margin mm), or model the support's contact",
            n2(peak.value),
            peak.at.map(n2),
            n2(d)
        ))
    }

    /// A warning when an elastic study's own solution leaves the
    /// small-rotation range its linear kinematics assume.
    pub fn kinematics_warning(&self) -> Option<String> {
        let w = self.max_rotation().filter(|w| w.value > LINEAR_ROTATION_LIMIT)?;
        let n2 = |v: f64| (v * 100.0).round() / 100.0;
        Some(format!(
            "large rotation: {} deg at {:?}, beyond the ~5 deg that linear (small-displacement) kinematics allow. \
             The study is linear: displacements and stresses scale with the load, and the motion that goes with a \
             rotation (a bent beam's shortening) and stress stiffening are missing, so read the results as qualitative",
            n2(w.value.to_degrees()),
            w.at.map(n2)
        ))
    }
}
