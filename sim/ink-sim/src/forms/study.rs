//! What a study is: its regions, physics, fixes and loads, the sim
//! package, and the `(study ...)` form.
use super::*;

/// Where a support or load acts: material where a field is negative, or the
/// material just under a named face (found by provenance, optionally filtered
/// by surface normal).
#[derive(Clone, Debug)]
pub enum Region {
    Field(NodeId),
    Face { tag: String, normal: Option<V3>, depth: Option<f64> },
}

impl Region {
    /// A region: a field (material where it is negative), or
    /// `(face name [(normal nx ny nz)] [(depth d)])` for a named face.
    pub(super) fn read(r: &mut Reader, e: &Sexp) -> R<Region> {
        let Some(parts) = e.list().filter(|l| l.first().and_then(|h| h.sym()) == Some("face")) else {
            return Ok(Region::Field(r.field(e)?));
        };
        let usage = "usage: (face name [(normal nx ny nz)] [(depth mm)]); names come from defs and (tag ...)";
        let Some(n) = parts.get(1) else { return Err(r.err_hint(e, "face needs a name", Some(usage.into()))) };
        let tag = match &n.k {
            SK::Sym(s) | SK::Str(s) => s.clone(),
            _ => return Err(r.err_hint(n, "face name must be a name", Some(usage.into()))),
        };
        let Some(tag) = r.tag(&tag) else {
            let hint = closest(&tag, r.graph().tag_names.iter().map(|s| s.as_str())).unwrap_or_else(|| usage.into());
            return Err(r.err_hint(n, format!("no piece built as (tag {tag} ...)"), Some(hint)));
        };
        let (mut normal, mut depth) = (None, None);
        for opt in &parts[2..] {
            match opt.list() {
                Some([h, a, b, c]) if h.sym() == Some("normal") => {
                    let v = [r.number(a)?, r.number(b)?, r.number(c)?];
                    normal = Some(ink_core::v3::direction(v).ok_or_else(|| r.err(opt, "normal must be non-zero (and finite)"))?);
                }
                Some([h, d]) if h.sym() == Some("depth") => depth = Some(r.number(d)?),
                _ => return Err(r.err_hint(opt, "bad face option", Some(usage.into()))),
            }
        }
        Ok(Region::Face { tag, normal, depth })
    }
}

#[derive(Clone, Copy, Debug)]
pub enum StudyPhysics {
    /// E: Young's modulus field (MPa). nu: Poisson's ratio field. density:
    /// kg/m^3 field (the language's g/mm^3 times 1e6). yield: MPa
    /// field (each material's own where it is, as E's: `(by-body root
    /// yield)`).
    Elastic { e: NodeId, nu: NodeId, density: Option<NodeId>, yield_strength: Option<NodeId> },
    /// k: thermal conductivity field, W/(m K).
    /// density (kg/m^3, as Elastic's) and capacity (J/(kg K)) for transient heat.
    Heat { k: NodeId, density: Option<NodeId>, capacity: Option<f64> },
}

/// Dirichlet condition: `dofs` of material points where `region < 0` held at `value`.
pub struct StudyFix {
    pub boundary: crate::Boundary,
    pub name: String,
    pub region: Region,
    pub dofs: [bool; 3],
    /// Elastic: the displacement it moves its region by, per axis (mm,
    /// fields of x y z: `(move ...)`; zero, a plain support). Heat: the
    /// temperature (the first; convection: the ambient's).
    pub value: [NodeId; 3],
    /// Convection: film coefficient W/(m^2 K) to an ambient temperature `value`.
    pub film: Option<f64>,
    /// Elastic: hold only the displacement along the surface normal (a
    /// frictionless support: a shaft on a bearing seat, a part on a table).
    pub normal: bool,
    /// Elastic: rigid motions of the region left free (a pin turns about
    /// its axis, a ball joint about its centre, a guide slides).
    pub free: Vec<crate::Motion>,
    /// Elastic, with free motions: hold the region rigidly but for them
    /// (`(rigid)`), not its mean motion.
    pub rigid: bool,
}

/// Total `vector` (force in N, or [power W, 0, 0]) spread over material in `region`.
pub struct StudyLoad {
    pub name: String,
    pub region: Region,
    pub vector: [f64; 3],
    /// Components as expressions of t (transient studies).
    pub vector_t: Option<[NodeId; 3]>,
    /// A surface pressure (MPa, any field) instead of a total force: the
    /// traction -p n on the part's surface inside the region.
    pub pressure: Option<NodeId>,
    /// A torque (N mm) instead of a force: about the region's centroid.
    pub torque: Option<[f64; 3]>,
    /// Its force's or torque's components where they read what the model works out (a reaction, a
    /// tension): numbers only when the study is asked, at the motion's worst moment (`worst_reads`).
    pub from_model: Option<[NodeId; 3]>,
}

pub struct Study {
    pub name: String,
    pub shape: NodeId,
    pub physics: StudyPhysics,
    pub fixes: Vec<StudyFix>,
    pub loads: Vec<StudyLoad>,
    /// m/s^2, applied through the density field.
    pub gravity: Option<[f64; 3]>,
    pub res: Option<f64>,
    /// Stresses within this distance (mm) of fix regions are reported apart.
    pub support_margin: Option<f64>,
    /// Modal analysis: number of vibration modes.
    pub modes: Option<usize>,
    pub transient: Option<crate::Transient>,
    /// A fix moves its region: the study has no compliance
    /// (MOVED_COMPLIANCE).
    pub moved: bool,
    /// Each material property given, and where it came from: `(E 2100MPa)`
    /// as written, or the `(material name)` that gave it.
    pub material_from: Vec<(String, String)>,
}

/// The sim package: `(study ...)`.
pub struct Sim;

impl Package for Sim {
    fn name(&self) -> &'static str {
        "sim"
    }

    fn forms(&self) -> &'static [&'static str] {
        &["study"]
    }

    fn commands(&self) -> &'static [ink_lang::command::Command] {
        &crate::command::COMMANDS
    }

    fn metrics(&self) -> &'static [(&'static str, &'static str)] {
        &[
            ("max-stress", "the peak von Mises stress (MPa) of an elastic study"),
            ("max-displacement", "the peak displacement (mm) of an elastic study"),
            ("compliance", "the compliance (N mm, the loads' work: lower is stiffer) of an elastic study"),
            (
                "min-safety-factor",
                "the least of the yield strength over the stress, point by point (each material's own yield where it is), of an elastic study whose material has (yield ...)",
            ),
            ("min-frequency", "the lowest natural frequency (Hz) of an elastic study with (modes n)"),
            (
                "min-buckling-factor",
                "the least factor on a static elastic study's loads that buckles it, by linear buckling analysis (an estimate: imperfections lower it, a thin shell's most)",
            ),
            ("max-temperature", "the peak temperature of a heat study (over time, if transient)"),
        ]
    }

    fn exclusions(&self) -> &'static [(&'static str, &'static str)] {
        &[
            ("support-margin", "the stress within the support margin of a fix that holds rigidly, where a clamp's or a pad's edge makes it singular (a symmetry plane has none)"),
            (
                "yield-between-reads",
                "a yield strength that falls between the points the stress is read at (a strip of weaker material narrower than the cells), counted at its least over each such cell",
            ),
            (
                "sharp-corners",
                "a sharp inside corner (or a round under about 0.4 of a cell) on the surface the stress is read on, where the elastic stress is unbounded whatever the cells read: left out, the claim is on the stress the cells read",
            ),
            (
                "rigid-motions",
                "the rigid motions a supported modal study's supports leave free, each a mode at 0 Hz: left out, the claim is on its elastic modes alone (a study with no fix is of a free part, and leaves them out without saying so)",
            ),
            (
                "buckling",
                "a buckle, which a linear static study does not model: where its loads may buckle it (a linear buckling factor whose band reaches 1; for a safety factor, reaches below it), or its linear buckling analysis gave no factor",
            ),
        ]
    }

    fn metric_order(&self, metric: &str) -> f64 {
        // A stress, a safety factor over it, is the displacement's
        // derivative: trilinear cells' converges at order 1.
        match metric {
            "max-stress" | "min-safety-factor" => 1.0,
            _ => 2.0,
        }
    }

    fn metric_quantity(&self, metric: &str) -> Option<Quantity> {
        match metric {
            "max-stress" => Some(Quantity::Pressure),
            "max-displacement" => Some(Quantity::Length),
            // Work, N mm: the unit of a torque.
            "compliance" => Some(Quantity::Torque),
            "min-frequency" => Some(Quantity::Frequency),
            "max-temperature" => Some(Quantity::Temperature),
            _ => None,
        }
    }

    fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
        Ok(std::sync::Arc::new(read_study(r, form, args)?))
    }
}

impl Form for Study {
    fn name(&self) -> &str {
        &self.name
    }

    /// A transient's loads on t are read over its own time, 0 to its duration.
    fn clock(&self) -> Option<(f64, f64)> {
        self.transient.map(|tr| (0.0, tr.duration))
    }

    fn results(&self) -> Vec<ResultSpec> {
        match self.physics {
            StudyPhysics::Heat { .. } => vec![ResultSpec::scalar("temperature", "°C"), ResultSpec::scalar("heat-flux", "W/mm²")],
            StudyPhysics::Elastic { .. } => vec![
                ResultSpec::scalar("von-mises", "MPa"),
                ResultSpec { vector: Some(["displacement-x", "displacement-y", "displacement-z"]), ..ResultSpec::scalar("displacement", "mm") },
                ResultSpec::scalar("displacement-x", "mm"),
                ResultSpec::scalar("displacement-y", "mm"),
                ResultSpec::scalar("displacement-z", "mm"),
            ],
        }
    }

    fn metrics(&self) -> Vec<&'static str> {
        match self.physics {
            StudyPhysics::Heat { .. } => vec!["max-temperature"],
            StudyPhysics::Elastic { yield_strength, .. } => {
                let mut m = vec!["max-stress", "max-displacement", "compliance"];
                m.extend(yield_strength.map(|_| "min-safety-factor"));
                m.extend(self.modes.map(|_| "min-frequency"));
                m.extend(self.modes.is_none().then_some("min-buckling-factor"));
                m
            }
        }
    }

    fn run(&self, doc: &mut Doc, ctx: &Ctx, opts: &RunOpts) -> R<Box<dyn Outcome>> {
        Ok(Box::new(SolvedStudy::run(doc, ctx, &self.name, opts)?))
    }

    fn linear_in_worst(&self, metric: &str) -> Option<bool> {
        // A static elastic solve is linear in its loads: its displacement and stress fields are, and the
        // largest of either is a seminorm of them. (A buckling factor, a frequency, a compliance are not.)
        let stat = matches!(self.physics, StudyPhysics::Elastic { .. }) && self.modes.is_none() && self.transient.is_none();
        match metric {
            "max-displacement" | "max-stress" if stat => Some(false),
            "min-safety-factor" if stat => Some(true),
            _ => None,
        }
    }

    fn besides_worst(&self) -> bool {
        self.gravity.is_some() || self.moved || self.loads.iter().any(|l| l.from_model.is_none())
    }

    fn run_each(&self, doc: &mut Doc, at: &[(Ctx, RunOpts)]) -> R<Vec<Box<dyn Outcome>>> {
        crate::solve_each(doc, at, &self.name)?;
        at.iter().map(|(ctx, opts)| self.run(doc, ctx, opts)).collect()
    }

    fn of_sums(&self, outcomes: &[&dyn Outcome], weights: &[Vec<f64>], metric: &str) -> Option<Vec<Measured>> {
        let solved: Vec<&SolvedStudy> = outcomes.iter().map(|o| o.as_any()?.downcast_ref()).collect::<Option<_>>()?;
        // One copy of a run, its fields set to each sum in turn.
        let mut sum = solved.first()?.copy();
        // What a run says of its own loads' buckling is not the sum's.
        sum.sol.report.warnings.retain(|w| !w.starts_with("buckles") && !w.starts_with("buckling not resolved"));
        (sum.sol.report.modeled_buckle, sum.alternate_low) = (None, None);
        weights
            .iter()
            .map(|row| {
                // A run's own moment: its own answer.
                let own = |i: usize| row.iter().enumerate().all(|(k, w)| (w - f64::from(k == i)).abs() <= 1e-9);
                if let Some(i) = (0..row.len()).find(|&i| own(i)) {
                    return solved[i].metric(metric);
                }
                let terms: Vec<(f64, &crate::Solution)> = row.iter().copied().zip(solved.iter().map(|s| &s.sol)).collect();
                if !sum.sol.set_sum(&terms) {
                    return None;
                }
                // Its buckling factor is no sum: bounded from the runs' own, each way (`buckle_of_sum`).
                sum.sol.report.buckle = Some(match crate::buckle_of_sum(&terms) {
                    Ok(low) => crate::Buckle { low, factor: f64::INFINITY, ..Default::default() },
                    Err(why) => crate::Buckle { failure: Some(format!("not bounded at this moment of the motion: {why}")), ..Default::default() },
                });
                sum.metric(metric)
            })
            .collect()
    }

    fn steady_besides(&self, doc: &mut Doc, at: &[Ctx], opts: &RunOpts) -> Result<(), String> {
        // The solver's question at each moment (`Setup::key`), the model's loads set aside.
        let modelled: Vec<bool> = self.loads.iter().map(|l| l.from_model.is_some()).collect();
        let mut first: Option<(u128, f64)> = None;
        for ctx in at {
            let (mut setup, _, _) = crate::study_setup(doc, ctx, &self.name, opts.res, opts.tol).map_err(|e| e.message)?;
            for (load, _) in setup.loads.iter_mut().zip(&modelled).filter(|(_, m)| **m) {
                (load.total, load.torque) = ([0.0; 3], load.torque.map(|_| [0.0; 3]));
            }
            // The key holds the moment only as far as the study reads it. With no key, nothing says what it reads.
            let Some((key, _)) = setup.key(ctx) else { return Err(format!("study '{}' is not one the solver can tell from another by what it reads", self.name)) };
            match first {
                None => first = Some((key, ctx.t)),
                Some((k, _)) if k == key => {}
                Some((_, t)) => {
                    return Err(format!(
                        "study '{}' is another study at t = {} than at t = {t}: something in it besides the loads it reads of the model follows t or the motion (a pressure, a \
                         material, a support's move, a region or its shape)",
                        self.name, ctx.t
                    ));
                }
            }
        }
        Ok(())
    }

    fn worst_reads(&self) -> Vec<(String, Vec<NodeId>)> {
        self.loads.iter().filter_map(|l| Some((format!("load '{}' ({})", l.name, if l.torque.is_some() { "N*mm" } else { "N" }), l.from_model?.to_vec()))).collect()
    }

    fn unavailable_in(&self, doc: &Doc, ctx: &Ctx, metric: &str) -> Option<String> {
        self.unavailable(metric).or_else(|| self.check_load_clock(doc, ctx).err().map(|e| e.message))
    }

    fn unavailable(&self, metric: &str) -> Option<String> {
        match metric {
            "compliance" if self.moved => Some(MOVED_COMPLIANCE.to_string()),
            _ => None,
        }
    }
}

impl Study {
    pub(crate) fn check_load_clock(&self, doc: &Doc, ctx: &Ctx) -> R<()> {
        if let Some(tr) = self.transient {
            for &f in self.loads.iter().filter_map(|l| l.vector_t.as_ref()).flatten() {
                crate::lang::motion::check_local_clock(doc, ctx, f, (0.0, tr.duration))?;
            }
        }
        Ok(())
    }
}
