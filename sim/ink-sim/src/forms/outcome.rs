//! A solved study's metrics (`Outcome`) and the read error of a peak
//! stress.
use super::*;

/// How far a peak stress may lie over its read beyond its grids' error, as
/// a fraction of the read: `STRESS_OVER_READ_CELLS` over the cells its
/// peak's concave radius spans (`peak_surface`), and at least
/// `STRESS_OVER_READ`. A curve's cut cells read its peak low by a fraction
/// that shrinks as its radius spans more cells, and three grids share the
/// miss, so their differences do not show it all: a flexure's leaf root
/// (a 0.75 mm fillet on 0.25 mm cells, 3 cells per radius) read 1.394
/// against CalculiX's 1.647, 2.5% past its grids' error. Fitted on 150
/// reads (81 inputs) with known answers (closed forms, converged CalculiX)
/// whose grids would decide: past their error the worst lay 2.2% over at
/// 5 cells per radius (11% times the cells), twice that allowed; the floor
/// keeps a bore whose field is not like a distance (8.5% over at 12 cells)
/// unknown, as it was.
pub const STRESS_OVER_READ_CELLS: f64 = 0.22;
pub const STRESS_OVER_READ: f64 = 0.04;
/// How far a peak stress may lie under its read beyond its grids' error:
/// at bonded interfaces, whatever the cells per radius, up to 8.2% (a ring
/// 0.189 mm thick bonded across materials, its least safety factor), 1.2
/// times that allowed.
pub const STRESS_UNDER_READ: f64 = 0.10;

/// The read error of a peak stress or a safety factor (yield over
/// stress, so its error inverted), its peak's concave radius spanning
/// `cells` of the cells it is read in (infinite where none is in reach),
/// and the solve's own field reaching `climbed` times the read between
/// the surface's reads (`climbed_scalar`; 1 where it reaches no more):
/// the fitted share over is of that. A convex edge oblique to the cells is
/// read inside it, where stress in bending is lower (a square bar turned
/// 37 degrees read 4% low on every grid, with its grids' error 0.3%,
/// against a share over of 4%): what its field reaches on the edge counts.
/// A peak on a curve under `RESOLVED_RADIUS_CELLS` of them is not read
/// (its grid cannot see the part), and its allowance is that of one on
/// that many: at a crease (none) an unbounded allowance (a safety factor's
/// not a number) would make what a claim leaves out, the support margin's
/// stress, answer it in place of the crease, and keep it unknown only
/// while a support is near. So both factors are finite: the cells at
/// least that many, the climb's at least 1 and, not finite, none.
fn stress_read(metric: &str, cells: f64, climbed: f64) -> Option<crate::lang::ReadError> {
    let cells = cells.max(RESOLVED_RADIUS_CELLS);
    let climbed = if climbed.is_finite() { climbed.max(1.0) } else { 1.0 };
    let (over, under) = ((1.0 + STRESS_OVER_READ.max(STRESS_OVER_READ_CELLS / cells)) * climbed - 1.0, STRESS_UNDER_READ);
    let how = "a peak stress may lie over its read by 22% over the cells its peak's radius spans (at least 4%), and under it by 10%, beyond its grid error: twice and 1.2 times the worst of 150 reads with known answers; over, of the most its solve's field reaches on the surface between its reads where that is more";
    match metric {
        "max-stress" => Some(crate::lang::ReadError { below: Some(under), above: Some(over), proven: false, how }),
        "min-safety-factor" => Some(crate::lang::ReadError { below: Some(over / (1.0 + over)), above: Some(under / (1.0 - under)), proven: false, how }),
        _ => None,
    }
}

impl Outcome for SolvedStudy {
    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn failure(&self) -> Option<String> {
        self.sol.valid().err()
    }

    fn metric(&self, metric: &str) -> Option<Measured> {
        let sol = &self.sol;
        let peak = |e: Result<crate::Extreme, String>, f: &dyn Fn(f64) -> f64| e.map(|e| (f(e.value), Some(e.at)));
        let measured = match metric {
            // Transient heat: the peak over time.
            // Transient heat: the peak over time, where and when it is.
            "max-temperature" if sol.report.peak.is_some() => sol.valid().map(|_| sol.report.peak.map_or((f64::NAN, None), |p| (p.1, Some(p.2)))),
            "max-temperature" => peak(sol.max_nodal(), &|v| v),
            "max-displacement" | "max-stress" | "min-safety-factor" => self.linear_value(metric).expect("one of its three"),
            "compliance" if self.moved => Err(MOVED_COMPLIANCE.to_string()),
            "compliance" => sol.valid().map(|_| (sol.report.energy, None)),
            "min-frequency" => sol.valid().and_then(|_| sol.report.modes.first().map(|m| (m.0, None)).ok_or_else(|| "the study found no modes".to_string())),
            "min-buckling-factor" => sol.valid().and_then(|_| self.buckling()).and_then(|b| {
                b.valid()?;
                let bk = b.report.buckle.as_ref().ok_or("the solve ran no buckling analysis")?;
                // Infinite where nothing is compressed: no buckling mode.
                bk.answer().map(|(_, f)| (f, f.is_finite().then_some(bk.at))).ok_or_else(|| bk.failure.clone().unwrap_or_else(|| "the buckling analysis found no factor".into()))
            }),
            _ => return None,
        };
        let measured = self.dropped_material().or_else(|| sol.report.incomplete.first().map(|e| e.reason(&sol.report))).map_or(measured, Err);
        let at = measured.as_ref().ok().and_then(|m| m.1);
        let beside = match metric {
            "max-stress" | "min-safety-factor" => sol.max_cell_scalar().ok().and_then(|p| self.beside_support_warning(&p)),
            _ => None,
        };
        // What the metric leaves out, each with the metric's value counting
        // it: the stress peaks, the support margin of each fix that holds
        // rigidly (a plane of symmetry has none); a buckle (of every static
        // metric: a buckled member's stresses and displacements are
        // unbounded); a gap the cells glue
        // shut (`unseen`); the modes antisymmetric about a plane of
        // symmetry. A claim that any would answer otherwise is undecided
        // (`Excluded`). The last cannot be left out on purpose: no package
        // exclusion names it (the whole part decides it).
        let mut excluded = Vec::new();
        if matches!(metric, "max-stress" | "min-safety-factor")
            && let Some(e) = sol.max_near_supports()
        {
            let held: Vec<&str> = self.fix_names.iter().zip(sol.plane_holds.iter().chain(std::iter::repeat(&false))).filter(|(_, p)| !**p).map(|(n, _)| n.as_str()).collect();
            excluded.push(crate::lang::Excluded {
                possible: false,
                key: "support-margin",
                // The metric counting it: the larger peak, or the least
                // safety factor.
                value: if metric == "max-stress" {
                    sol.max_cell_scalar().map_or(e.value, |p| p.value.max(e.value))
                } else {
                    // Each point's own yield, and the bounds of its cells not
                    // of one material: the stress's bound over a cell by a
                    // rigid hold is its singular stress too (a graded bar's
                    // fell from 17.8 to 13.9 as the cells halved).
                    let bound = self.safety_factor_between_reads(true).map_or(f64::INFINITY, |e| e.value);
                    let near = self.min_safety_factor_near_supports().map_or(f64::INFINITY, |e| e.value).min(bound);
                    self.min_safety_factor().map_or(near, |p| p.value.min(near))
                },
                at: Some(e.at),
                what: format!(
                    "the stress within {} mm of fix {} (the support margin: where a rigid hold ends on the free surface its stress is singular, and grows as the cells shrink)",
                    (self.margin * 1e4).round() / 1e4,
                    held.iter().map(|n| format!("'{n}'")).collect::<Vec<_>>().join(", ")
                ),
                decide: "model the support's contact (the washer, pad or mating part in the study, or a softer layer between) and read the stress there".into(),
                singular: true,
            });
        }
        // A sharp inside corner anywhere on the surface read: the elastic
        // stress there has no bound, so no `max-stress <` or
        // `min-safety-factor >` is true of the part, wherever its cells
        // read their peak (a shaft's keyway, its corners sharp, held `< 1e6`
        // by the 160.7 MPa read at its shoulder's fillet).
        if matches!(metric, "max-stress" | "min-safety-factor")
            && measured.is_ok()
            && let Some(at) = *self.crease.get_or_init(|| sol.crease(&self.shape_field.0, &self.shape_field.1, self.judge, RESOLVED_RADIUS_CELLS))
        {
            excluded.push(crate::lang::Excluded {
                possible: false,
                key: "sharp-corners",
                value: if metric == "max-stress" { f64::INFINITY } else { 0.0 },
                at: Some(at),
                what: format!(
                    "a sharp inside corner, or a round under about 0.4 of its cells (a crease to them), at [{:.3}, {:.3}, {:.3}]: wherever load passes it the elastic stress there is unbounded, whatever the cells read",
                    at[0], at[1], at[2]
                ),
                decide: format!(
                    "round the corner and size the cells to the round: a claim compares three grids, the coarsest four times its cells, and each reads a round of at least \
                     {RESOLVED_RADIUS_CELLS} of its own cells, so a round of radius R wants (res r) at R / {} or finer (R 8 mm: 1.3 mm cells; at these {} mm cells, a round of \
                     {} mm). A filleted part solves about as fast as a sharp one",
                    4.0 * RESOLVED_RADIUS_CELLS,
                    self.sol.h,
                    4.0 * RESOLVED_RADIUS_CELLS * self.sol.h
                ),
                singular: true,
            });
        }
        // A buckle: its linear buckling's band reaching 1 (the loads may
        // buckle it), or no factor resolved (whether they do is not known).
        // A safety factor counts it as the least factor that may buckle
        // the part (its band's low end), where that is under its own.
        let sf = measured.as_ref().ok().filter(|_| metric == "min-safety-factor").map(|m| m.0);
        // A yield strength that is not one value across a cell: each such
        // cell bounded (its least yield over its stress's bound over the
        // cell), and the claim that bound would change unknown.
        // Away from the supports (the margins' cells count in theirs).
        if metric == "min-safety-factor"
            && let Some(e) = self.safety_factor_between_reads(false)
        {
            // The metric counting it: the least of the value and the bound.
            let value = measured.as_ref().map_or(e.value, |m| m.0.min(e.value));
            excluded.push(crate::lang::Excluded { possible: false,
                key: "yield-between-reads",
                value,
                at: Some(e.at),
                what: format!(
                    "the cells of {} mm not of one material (E, nu or yield not one value across them: a strip, a gradient, an interface), no point of which is read as its own: counted as each one's least yield against the stress's bound over it (at a cell near {:?})",
                    sol.h,
                    e.at.map(|v| (v * 1e3).round() / 1e3)
                ),
                decide: "refine the cells (res): the bound closes on the factor read as they shrink, and a weak strip is read once they are under half its width".into(),
                singular: false,
            });
        }
        let static_result = matches!(metric, "max-stress" | "min-safety-factor" | "max-displacement" | "compliance");
        if static_result {
            let mut count_buckle = |key, b: &crate::Buckle| {
                let low = b.answer().map_or(0.0, |(lo, _)| lo);
                if low > 1.0 && !sf.is_some_and(|x| low < x) {
                    return;
                }
                let what = match b.answer() {
                    Some((lo, hi)) => {
                        let band =
                            if hi.is_finite() { format!("{} to {}", claims::sig(lo), claims::sig(hi)) } else { format!("at least {} with no upper estimate", claims::sig(lo)) };
                        format!("a buckle: the estimated linear buckling factor is {band} on cells of {} mm; the linear study cannot predict the response past buckling", sol.h)
                    }
                    None => format!("a buckle the eigensolve could not rule out ({})", b.failure.as_deref().unwrap_or("no factor resolved")),
                };
                excluded.push(crate::lang::Excluded {
                    key,
                    possible: true,
                    value: sf.map_or(f64::INFINITY, |x| x.min(low)),
                    at: None,
                    what,
                    decide: "lower the loads, stiffen the part, or refine (res r) to resolve the factor".into(),
                    singular: false,
                });
            };
            if let Some(b) = sol.report.buckle.as_ref() {
                count_buckle("buckling", b);
            }
            if let Some(b) = sol.report.modeled_buckle.as_ref().or_else(|| (!sol.report.possible_cuts.is_empty()).then_some(sol.report.buckle.as_ref()).flatten()) {
                count_buckle("modeled-buckling", b);
            }
        }
        if let Some(what) = self.hidden_modes_warning()
            && (static_result || matches!(metric, "min-frequency" | "min-buckling-factor"))
            && !(static_result && sf.is_none() && self.alternate_low.is_some_and(|low| low > 1.0))
        {
            excluded.push(crate::lang::Excluded {
                key: if static_result { "buckling" } else { "symmetry-modes" },
                possible: true,
                value: if static_result && sf.is_none() { f64::INFINITY } else { self.alternate_low.unwrap_or(0.0).min(measured.as_ref().map_or(f64::INFINITY, |v| v.0)) },
                at: None,
                what,
                decide: self.mode_hint().into(),
                singular: false,
            });
        }
        // A supported part its supports leave free to move as a whole (a
        // forgotten roller, a gap cutting it from them): each such motion a
        // mode at 0 Hz, taken out of the modes. Its lowest frequency counts
        // them unless the claim says it leaves them out; a study with no fix
        // is of a free part, whose motions are none of its vibration.
        let rigid = sol.report.rigid_modes;
        if metric == "min-frequency" && rigid > 0 && !self.fix_names.is_empty() {
            excluded.push(crate::lang::Excluded {
                possible: false,
                key: "rigid-motions",
                value: 0.0,
                at: None,
                what: format!("{rigid} rigid motion{} the supports leave free (a mode at 0 Hz each), taken out of the modes found", if rigid == 1 { "" } else { "s" }),
                decide: "hold the motion (another fix, or one without the free motion)".into(),
                singular: false,
            });
        }
        // A peak stress on a surface its cells do not resolve: read only in
        // part, alike on every grid (`RESOLVED_RADIUS_CELLS`).
        let k = RESOLVED_RADIUS_CELLS;
        let curved = at.filter(|_| matches!(metric, "max-stress" | "min-safety-factor")).and_then(|p| self.peak_surface(p).map(|s| (p, s))).filter(|(_, (r, h))| r / h < k).map(|(p, (r, h))| {
            let sig = crate::lang::claims::sig;
            if r == 0.0 {
                return format!(
                    "its peak stress, at [{:.3}, {:.3}, {:.3}], lies at a sharp inside corner or a round under about 0.4 of its {h} mm cells (a crease to them), where stress is singular or unread: round the corner, or refine (res r) until a round spans 1.45 of them",
                    p[0], p[1], p[2]
                );
            }
            format!(
                "its peak stress, at [{:.3}, {:.3}, {:.3}], lies on a surface curved to a radius of {} mm, {} of its {h} mm cells: cells read a peak on a radius of at least {k} of them, so these read it only in part; refine (res r) to {} mm or finer (a claim's coarsest grid, four times its cells, resolving it), or round the surface",
                p[0],
                p[1],
                p[2],
                sig(r),
                sig(r / h),
                // The claim's res at which its grid four times as coarse
                // resolves it: its leaf cells here are h / sol.h of a grid's.
                sig(sol.h * r / (4.0 * k * h))
            )
        });
        // An unresolved curve at the peak is no reason for finer grids.
        // (An interface not flat along the cells is not on finer ones
        // either.)
        let stepped = sol.report.interface_stepped;
        let finer = curved.is_none() && stepped.is_none();
        let read_x = measured.as_ref().ok().map(|m| m.0);
        let band = measured.as_ref().ok().and_then(|_| self.buckling_band(metric).or_else(|| self.bound_band(metric)).or_else(|| self.time_band(metric)));
        let unbounded = self.unbounded(metric, &excluded, at);
        let extent = matches!(metric, "min-buckling-factor" | "min-frequency").then(|| sol.mode_extent(1)).flatten();
        Some(Measured {
            value: measured.map(|m| m.0),
            at,
            res: Some(sol.h),
            finest: (sol.report.refine_depth > 0).then(|| sol.h / f64::from(1u32 << sol.report.refine_depth)),
            warnings: self.kinematics_warning().into_iter().chain(beside).chain(sol.report.warnings.iter().cloned()).chain(sol.report.incomplete.iter().map(|e| e.reason(&sol.report))).collect(),
            excluded,
            // What the grid cannot see: a gap it glues shut, or a part
            // with no cell at least half material (walls thinner than half
            // a cell: their stiffness is averaged in, their stress unread);
            // a peak stress on a surface curved tighter than its cells read.
            unseen: curved
                .or(sol.report.glued
                .map(|(cells, at, width)| {
                    let what = if width > 0.0 { format!("a gap about {width} mm wide at {at:?}") } else { format!("a gap or hairline whose sides bounds do not show joined, at {at:?},") };
                    format!(
                        "{what} that its {} glue shut ({cells} cells hold material on both sides of air too thin for their samples to part): it solved the part with the gap filled",
                        ink_fem::cells_of(sol.h, sol.report.refine_depth)
                    )
                })
                .or_else(|| {
                    sol.report.thin_left.map(|(cells, at)| {
                        format!(
                            "a feature thinner than two of its {} left: {cells} cells, one centred at {at:?}, hold material or air too thin for them beyond a taper's edge (a fin, a web or a ligament, held too stiff and read only in part; a slit or a clearance, its sides moved as one)",
                            ink_fem::cells_of(sol.h, sol.report.refine_depth)
                        )
                    })
                })
                .or_else(|| sol.report.material_left.map(|(cells, at)| ink_fem::material_left_what(cells, at, &ink_fem::cells_of(sol.h, sol.report.refine_depth))))
                .or_else(|| sol.report.interface_shared.map(|(cells, at, plane)| ink_fem::interface_shared_what(cells, at, plane, &ink_fem::cells_of(sol.h, sol.report.refine_depth))))
                .or_else(|| {
                    stepped.map(|(cells, at)| ink_fem::interface_stepped_what(cells, at, &ink_fem::cells_of(sol.h, sol.report.refine_depth)))
                })
                .or_else(|| sol.report.surface_measures.iter().find_map(|m| m.unseen()))
                .or_else(|| {
                    (sol.report.unresolved_volume >= sol.report.volume).then(|| format!("no cell of {} mm is at least half material: its walls are thinner than half a cell", sol.h))
                })),
            unseen_finer: finer,
            // What the stiffness kept in void moves a reading by, at most (ink_fem::VOID_SHARE): of what a
            // displacement or temperature varies by over the part, an offset carried whole being exact (a
            // block moved 1000 mm, its deflection still growing, read "alike" against 3.2e-5 of 1000 mm; a
            // block held at 1000 degC, its rise still growing, against 3.2e-5 of 1273 K); of the value else.
            resolves: read_x.filter(|x| x.is_finite()).map_or(0.0, |x| {
                ink_fem::VOID_SHARE * if matches!(metric, "max-displacement" | "max-temperature") { sol.span(x) } else { x.abs() }
            }),
            unbounded,
            extent,
            // A reading that states nothing else of its precision (a
            // displacement, a temperature) is known to what the solver can
            // give it (`Solution::precision`): its grids differing by less
            // agree.
            read_error: stress_read(metric, at.filter(|_| matches!(metric, "max-stress" | "min-safety-factor")).and_then(|p| self.peak_surface(p)).map_or(f64::INFINITY, |(r, h)| r / h), self.climbed(metric))
                .or_else(|| band.is_none().then(|| {
                    // (derived for a node's displacement or temperature: the
                    // field's span is theirs. A frequency or a compliance has
                    // its solve's residual alone until one is derived in its
                    // units: a mode's span is its amplitude's, and a
                    // compliance's floor came out 1e-6 over the load in N.)
                    if matches!(metric, "max-displacement" | "max-temperature") {
                        let p = read_x.map_or(f64::INFINITY, |x| sol.precision_at(x));
                        crate::lang::ReadError { below: Some(p), above: Some(p), proven: false, how: "the solver's precision: the stiffness it keeps in void and held material, 1e-6 of the material's, over its cells' share of the part, and its solve's residual, of the field's span over the part (not of an offset it carries whole)" }
                    } else {
                        let p = sol.report.residual.max(0.0);
                        crate::lang::ReadError { below: Some(p), above: Some(p), proven: false, how: "the solve's residual" }
                    }
                })),
            proven: false,
            held: None,
            band,
        })
    }

    /// Compliance: minus the strain energy density on the boundary
    /// (Hadamard's formula, one solve), for boundaries free of loads and
    /// supports. It covers the shape only: a param any other input reads
    /// (a material, a load, a support) is differenced. A study a fix moves
    /// has no compliance.
    fn adjoint(&self, metric: &str) -> Option<Adjoint<'_>> {
        (metric == "compliance" && matches!(self.physics, StudyPhysics::Elastic { .. }) && self.sol.valid().is_ok() && !self.moved).then(|| Adjoint {
            shape: self.shape,
            weight: Box::new(|p| -self.sol.energy_density(p).unwrap_or(0.0)),
            uncovered: self.inputs.clone(),
        })
    }

    fn undecided(&self, result: &str) -> Option<String> {
        if let Some(why) = self.dropped_material().or_else(|| self.sol.report.incomplete.first().map(|e| e.reason(&self.sol.report))) {
            return Some(why);
        }
        if let Some(why) = self.sol.report.surface_measures.iter().find_map(crate::SurfaceMeasure::unseen) {
            return Some(format!("{why}. Refine (res r) or select a surface the cells resolve"));
        }
        if self.mode_field {
            // A normalized eigenvector is meaningful even if the applied
            // loads exceed the buckling load. Its rank still needs all modes.
            return self.hidden_modes_warning();
        }
        let b = self.sol.report.buckle.as_ref()?;
        if self.alternate_low.is_none_or(|low| low <= 1.0)
            && let Some(why) = self.hidden_modes_warning()
        {
            return Some(why);
        }
        if b.answer().is_none_or(|(low, _)| low <= 1.0) {
            return Some("the loads may buckle the part; a linear static field cannot determine the response. Lower the load, stiffen the part or refine (res r) to resolve its buckling factor".into());
        }
        // Its stress as a field is each cell's at its centre: no peak on the surface, none by a support or
        // at a sharp corner, no read's error, which `max-stress` and `min-safety-factor` measure and say.
        matches!(result, "von-mises").then(|| {
            "a study's von-mises field is each cell's stress at its centre, for showing and for shaping: it holds no peak on the surface, by a support or at a sharp corner, and no error \
             of its read. Claim the strength with (max-stress study < bound) or (min-safety-factor study > bound), which measure those"
                .into()
        })
    }

    fn samples(&self, q: &str) -> Option<Samples> {
        self.sol.valid().ok()?;
        let gf = self.sol.grid_fields();
        let (cells, nodes): ([usize; 3], [usize; 3]) = (gf.dims, gf.dims.map(|d| d + 1));
        let nd = gf.node_values.len() / (nodes[0] * nodes[1] * nodes[2]);
        let centres = gf.origin.map(|o| o + 0.5 * gf.h);
        let (values, dims, origin): (Vec<f32>, [usize; 3], V3) = match q {
            "von-mises" | "heat-flux" => (gf.cell_scalar.clone(), cells, centres),
            "displacement" if nd == 3 => (gf.node_values.chunks(3).map(|c| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()).collect(), nodes, gf.origin),
            "displacement-x" | "displacement-y" | "displacement-z" if nd == 3 => {
                let c = (q.as_bytes()[q.len() - 1] - b'x') as usize;
                (gf.node_values.chunks(3).map(|v| v[c]).collect(), nodes, gf.origin)
            }
            "temperature" if nd == 1 => (gf.node_values.clone(), nodes, gf.origin),
            _ => return None,
        };
        // A cell's result covers its material; a node's, all of it.
        let cover = if dims == cells {
            values.iter().zip(&gf.cell_frac).map(|(v, f)| if v.is_nan() { 0.0 } else { *f }).collect()
        } else {
            values.iter().map(|v| if v.is_nan() { 0.0 } else { 1.0 }).collect()
        };
        Some(Samples { origin, h: gf.h, dims, values, cover })
    }

    /// Vibration mode k of a modal study: the solution's shape k (1 mm
    /// peak) and its stresses, which are
    /// per that 1 mm of peak displacement (a mode's amplitude is arbitrary):
    /// its description says so ("mode 2: 821.6 Hz, 1 mm peak").
    fn mode(&self, k: usize) -> Result<(Box<dyn Outcome>, String), String> {
        let sig = |v: f64| {
            if !v.is_finite() {
                return "unbounded".to_string();
            }
            let d = (3 - v.abs().max(1e-300).log10().floor() as i32).clamp(0, 12) as usize;
            format!("{v:.d$}")
        };
        let r = &self.sol.report;
        let (sol, what) = if !r.modes.is_empty() {
            let n = r.modes.len();
            let sol = self.sol.mode(k).ok_or_else(|| format!("it gives {n} vibration mode{} (mode 1 to {n}: (modes n) asks for more)", if n == 1 { "" } else { "s" }))?;
            (sol, format!("mode {k}: {} Hz, 1 mm peak", sig(r.modes[k - 1].0)))
        } else if matches!(self.physics, StudyPhysics::Elastic { .. }) {
            // A static study's modes are its buckling modes.
            let b = self.buckling()?;
            b.valid()?;
            let bk = b.report.buckle.as_ref().ok_or("the solve ran no buckling analysis")?;
            if let Some(f) = &bk.failure {
                return Err(b.report.incomplete.first().map_or_else(|| f.clone(), |e| e.reason(&b.report)));
            }
            let sol = b.mode(k).ok_or("its buckling analysis gives its first buckling mode only, (mode 1)")?;
            (sol, format!("buckling mode 1: factor {}, 1 mm peak", sig(bk.factor)))
        } else {
            return Err("it is no modal study (give it (modes n))".into());
        };
        // The study's mode: every other field is the study's own.
        Ok((
            Box::new(SolvedStudy {
                sol,
                physics: self.physics,
                yield_strength: self.yield_strength,
                inputs: self.inputs.clone(),
                fix_names: self.fix_names.clone(),
                yield_field: self.yield_field.clone(),
                stiffness: self.stiffness.clone(),
                symmetry: self.symmetry.clone(),
                alternate_low: self.alternate_low,
                mode_field: true,
                shape: self.shape,
                shape_field: self.shape_field.clone(),
                surfaces: Default::default(),
                crease: Default::default(),
                judge: self.judge,
                moved: self.moved,
                margin: self.margin,
            }),
            what,
        ))
    }
}
