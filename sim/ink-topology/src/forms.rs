//! `(topology ...)`: a study's shape as a design domain, its results (the
//! density, and the optimized part as a shape) as fields.

use crate::eval::Ctx;
use crate::graph::NodeId;
use crate::lang::*;
use crate::sexp::Sexp;
use crate::{InkError, sim};
use ink_core::grid::{self, GridData};

type R<T> = Result<T, InkError>;

/// `(topology name (study s) ...)`: the stiffest layout of study `s`'s
/// shape within a volume budget.
pub struct TopologyDef {
    pub name: String,
    /// The study it lays out: a static elastic one whose fixes hold still.
    pub study: String,
    pub volume_fraction: f64,
    pub iterations: usize,
    /// Filter radius, mm (default 1.5 cells).
    pub filter: Option<f64>,
    /// Material here stays solid.
    pub keep: Option<NodeId>,
}

/// The topology package: `(topology ...)`.
pub struct Topology;

impl Package for Topology {
    fn name(&self) -> &'static str {
        "topology"
    }

    fn forms(&self) -> &'static [&'static str] {
        &["topology"]
    }

    fn commands(&self) -> &'static [ink_lang::command::Command] {
        &crate::command::COMMANDS
    }

    fn read(&self, r: &mut Reader, form: &Sexp, _head: &str, args: &[Sexp]) -> R<std::sync::Arc<dyn Form>> {
        Ok(std::sync::Arc::new(read_topology(r, form, args)?))
    }
}

/// (topology name (study s) [(volume-fraction v)] [(iterations n)] [(filter mm)] [(keep field)])
fn read_topology(r: &mut Reader, form: &Sexp, args: &[Sexp]) -> R<TopologyDef> {
    let usage = "usage: (topology name (study s) [(volume-fraction 0.3)] [(iterations 60)] [(filter mm)] [(keep field)]): \
                 a stiff layout, by local search, of static elastic study s's shape using that fraction of it";
    let Some((name, clauses)) = args.split_first() else { return Err(r.err_hint(form, "bad topology", Some(usage.into()))) };
    let n = r.name(name)?;
    r.new_result_name(name, &n)?;
    let mut t = TopologyDef { name: n, study: String::new(), volume_fraction: 0.3, iterations: 60, filter: None, keep: None };
    for c in &r.clauses(clauses)? {
        let (head, v) = match c.list() {
            Some([k, v]) => (r.name(k)?, v),
            _ => return Err(r.err_hint(c, "expected (key value)", Some(usage.into()))),
        };
        match head.as_str() {
            "study" => {
                let s = r.name(v)?;
                // Kept by its public name: the study it read, found again.
                let Some((public, st)) = r.form::<sim::Study>(&s) else { return Err(r.err(v, format!("no study named '{s}' before this topology"))) };
                if !matches!(st.physics, sim::StudyPhysics::Elastic { .. }) || st.modes.is_some() || st.transient.is_some() {
                    return Err(r.err(v, format!("study '{s}': a topology lays out a static elastic study")));
                }
                if st.moved {
                    return Err(r.err(v, format!("study '{s}': a topology minimizes the compliance under loads: its fixes hold still (no (move ...))")));
                }
                t.study = public;
            }
            "volume-fraction" => t.volume_fraction = r.number(v)?,
            "iterations" => {
                let n = r.number(v)?;
                if !(n >= 1.0 && n.fract() == 0.0 && n <= 1e6) {
                    return Err(r.err(v, format!("iterations must be a whole number, at least 1 (got {n})")));
                }
                t.iterations = n as usize;
            }
            "filter" => {
                let f = r.number_of(v, Quantity::Length, "filter")?;
                if !(f > 0.0) {
                    return Err(r.err(v, format!("filter must be a positive length, the density filter's radius (got {f})")));
                }
                t.filter = Some(f);
            }
            "keep" => t.keep = Some(r.field(v)?),
            other => return Err(r.err_hint(c, format!("unknown topology option '{other}' (study volume-fraction iterations filter keep)"), Some(usage.into()))),
        }
    }
    if t.study.is_empty() {
        return Err(r.err_hint(form, "topology needs (study s)", Some(usage.into())));
    }
    if !(t.volume_fraction > 0.0 && t.volume_fraction < 1.0) {
        return Err(r.err(form, "volume-fraction must be between 0 and 1"));
    }
    Ok(t)
}

impl Form for TopologyDef {
    fn name(&self) -> &str {
        &self.name
    }

    /// The density per cell (0..1), and the optimized part as a shape's
    /// field: 1/2 - density, negative where material was kept.
    fn results(&self) -> Vec<ResultSpec> {
        vec![ResultSpec::scalar("density", ""), ResultSpec { shape: true, ..ResultSpec::scalar("solid", "") }]
    }

    fn run(&self, doc: &mut Doc, ctx: &Ctx, opts: &RunOpts) -> R<Box<dyn Outcome>> {
        Ok(Box::new(crate::run(doc, ctx, &self.name, opts.res, opts.tol)?))
    }
}

impl Outcome for crate::Layout {
    fn failure(&self) -> Option<String> {
        self.sol.valid().err()
    }

    fn samples(&self, q: &str) -> Option<Samples> {
        self.sol.valid().ok()?;
        if q != "density" && q != "solid" {
            return None;
        }
        let gf = self.sol.grid_fields();
        let dens = self.sol.cells_on_grid(&self.record.density);
        let values: Vec<f32> = if q == "solid" { dens.iter().map(|d| 0.5 - d).collect() } else { dens };
        // A cell's result covers its material.
        let cover = values.iter().zip(&gf.cell_frac).map(|(v, f)| if v.is_nan() { 0.0 } else { *f }).collect();
        Some(Samples { origin: gf.origin.map(|o| o + 0.5 * gf.h), h: gf.h, dims: gf.dims, values, cover })
    }

    fn field(&self, q: &str) -> Option<GridData> {
        let s = self.samples(q)?;
        // Past the material: 0 (void density for "solid" is +1/2).
        let data = sampled_result(s.values, s.dims, s.origin, s.h, if q == "solid" { 0.5 } else { 0.0 });
        if q != "solid" {
            return Some(data);
        }
        Some(solid_spline(&data))
    }
}

/// Pad the controls with one layer of void. On each outer face at least
/// 5/6 of the spline's weight is on void, so the zero surface closes
/// smoothly inside the grid box.
fn solid_spline(g: &GridData) -> GridData {
    let inner = g.dims();
    let dims = inner.map(|n| n + 2);
    let values = (0..dims.iter().product())
        .map(|i| {
            let q = grid::ijk(i, dims);
            if (0..3).any(|a| q[a] == 0 || q[a] == dims[a] - 1) { 0.5 } else { g.values()[grid::index(q.map(|i| i - 1), inner)] }
        })
        .collect();
    GridData::new(g.origin().map(|o| o - g.h()), g.h(), dims, values).cubic_bspline()
}

#[cfg(test)]
mod tests {
    use super::solid_spline;
    use ink_core::grid::GridData;

    #[test]
    fn solid_closes_inside_its_spline_grid() {
        let g = GridData::new([0.0; 3], 1.0, [3; 3], vec![-0.5; 27]);
        let s = solid_spline(&g);
        assert!((s.component([1.0; 3], 0) + 0.5).abs() < 1e-15);
        for p in [[3.0, 1.0, 1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, 3.0]] {
            assert!(s.component(p, 0) >= 1.0 / 3.0 - 1e-15);
        }
    }
}
