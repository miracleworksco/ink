//! Topology optimization: the `(topology ...)` form, the stiffest layout of
//! a static elastic study's shape within a volume budget (minimum
//! compliance, SIMP). It names a study (ink-sim) and runs its loop through
//! the finite element solver's redesign interface (`ink_fem::Assembly`):
//! each iteration scales the cells' stiffness by their densities and
//! solves, then moves the densities by the method of moving asymptotes
//! (ink-mma). It depends on core, kernel, language, sim and mma, never on
//! another package (../kernel/docs/CORE.md).

// The crates it builds on, at the paths the package's code uses.
use ink_core::{V3, eval, graph, tape};
use ink_kernel::cache;
use ink_lang as lang;
use ink_lang::InkError;
use ink_lang::sexp;
use ink_sim as sim;

mod command;
mod forms;
pub use self::forms::*;

use crate::eval::Ctx;
use crate::tape::Tape;
use rayon::prelude::*;
use std::sync::Arc;

/// A run settles only within its volume budget to this fraction of the
/// domain.
const VOLUME_TOL: f64 = 1e-3;

#[cfg(test)]
thread_local! {
    /// Per iteration: compliance, the cells' energies (times their
    /// stiffness) summed, and the slope's sum over the design (its change
    /// with a uniform change of every design variable).
    static PROBE: std::cell::RefCell<Vec<(f64, f64, f64)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Minimum-compliance topology optimization (SIMP): a density per cell,
/// stiffness density^penalty, a volume budget, a density filter and
/// projection, updates by the method of moving asymptotes. Material in
/// `keep` stays solid.
#[derive(Clone)]
pub struct Settings {
    pub volume_fraction: f64,
    pub iterations: usize,
    /// Filter radius, mm (default 1.5 cells).
    pub filter: Option<f64>,
    pub keep: Option<Arc<Tape>>,
}

/// What a run found, besides the solution at its design.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Record {
    /// Density per cell (0..1), as `Assembly` numbers them.
    pub density: Vec<f64>,
    /// Per design solved (the first, then after each update), its number,
    /// objective (the compliance), volume fraction, the largest change of a
    /// cell's density from the one before, and the projection's steepness
    /// it was solved at. The last row is the design returned.
    pub history: Vec<(usize, f64, f64, f64, f64)>,
    /// Whether the design settled (no density moving by 0.01 at the
    /// steepest projection, within the volume budget) before the iterations
    /// ran out.
    pub settled: bool,
    /// The fraction of the domain held solid (by `keep`, or wholly held by
    /// a fix), and the target
    /// volume fraction. When the first reaches the second nothing is left
    /// to optimize: free material vanishes.
    pub fixed: (f64, f64),
}

/// Runs the optimization on an assembled study, its iterations' solves to
/// `tol` (looser, 1e-4, until the design returned). `ctx` is the one the
/// solve runs in (`keep` is read in it).
pub fn optimize(a: &mut dyn sim::Assembly, top: &Settings, ctx: &Ctx, tol: f64) -> Result<Record, String> {
    let nc = a.cells();
    let h = a.h();
    let cells: Vec<([i64; 3], V3, f64)> = (0..nc).map(|i| a.cell(i)).collect();
    // Solid for good: `keep`, and each cell a fix holds all the material
    // of (passive elements, as Sigmund's and Andreassen et al.'s codes take
    // them). A hold is part of the matrix of the cell its free material is
    // in, scaled with that cell's stiffness: it holds what material the
    // optimizer leaves there. A cell wholly held has none, and its hold is
    // its neighbours' (a fix's edge on the grid's plane): its material
    // held nothing the optimizer saw, so a column fixed below x = 1 on 1
    // mm cells was laid out with none there, and the part it drew was held
    // by nothing (its buckling factor 4e-6, against 0.006 so held). A cell
    // the region cuts keeps its hold and is designed (a plane held through
    // every cell of a plate leaves it all to optimize). A penalty at the
    // solid's over void cells, 1e12 of their stiffness, stalled the solve
    // (a plate held along a plane through every cell: residual 6e-4 after
    // 17,493 iterations).
    let centers: Vec<V3> = cells.iter().map(|c| c.1).collect();
    let mut fixed = a.held();
    if let Some(k) = &top.keep {
        for (i, v) in eval::eval_points_par(k, ctx, &centers).iter().enumerate() {
            if *v < 0.0 {
                fixed[i] = true;
            }
        }
    }
    let vol: Vec<f64> = cells.iter().map(|c| c.2).collect();
    let total: f64 = vol.iter().sum();
    let held = ((0..nc).filter(|&i| fixed[i]).map(|i| vol[i]).sum::<f64>() / total, top.volume_fraction);
    // Filter neighbourhoods: cells (each piece of a cut cell) within
    // radius r, centre to centre.
    let r = top.filter.unwrap_or(1.5 * h).max(h);
    let reach = (r / h).floor() as i64;
    let neighbours: Vec<Vec<(u32, f64)>> = cells
        .iter()
        .map(|&(idx, _, _)| {
            let mut out = Vec::new();
            for dz in -reach..=reach {
                for dy in -reach..=reach {
                    for dx in -reach..=reach {
                        let w = r - h * ((dx * dx + dy * dy + dz * dz) as f64).sqrt();
                        if w <= 0.0 {
                            continue;
                        }
                        out.extend(a.cells_at([idx[0] + dx, idx[1] + dy, idx[2] + dz]).map(|p| (p as u32, w)));
                    }
                }
            }
            out
        })
        .collect();
    // A design variable per cell, a kept cell's held at 1; MMA moves the
    // free ones.
    let free: Vec<usize> = (0..nc).filter(|&i| !fixed[i]).collect();
    let pick = |v: &[f64]| -> Vec<f64> { free.iter().map(|&i| v[i]).collect() };
    let expand = |y: &[f64]| -> Vec<f64> {
        let mut v = vec![1.0; nc];
        free.iter().zip(y).for_each(|(&i, y)| v[i] = *y);
        v
    };
    let mut xs: Vec<f64> = fixed.iter().map(|&f| if f { 1.0 } else { top.volume_fraction }).collect();
    // Three-field scheme: design x, filtered x~ (density filter, radius
    // r), projected x^ (smoothed Heaviside at 1/2, sharpened by doubling
    // beta), stiffness x^^3. The projection drives the result to solid
    // and void, so the extracted part has crisp walls and no grey specks.
    let weight: Vec<f64> = neighbours.iter().map(|n| n.iter().map(|&(_, w)| w).sum()).collect();
    let filter = |v: &[f64]| -> Vec<f64> { (0..nc).into_par_iter().map(|i| neighbours[i].iter().map(|&(j, w)| w * v[j as usize]).sum::<f64>() / weight[i]).collect() };
    // Its transpose, from a cell's slope to the design variables'.
    let filter_t = |g: &[f64]| -> Vec<f64> {
        let gw: Vec<f64> = (0..nc).map(|i| g[i] / weight[i]).collect();
        (0..nc).into_par_iter().map(|j| if fixed[j] { 0.0 } else { neighbours[j].iter().map(|&(i, w)| w * gw[i as usize]).sum() }).collect()
    };
    let project = |beta: f64, t: f64| -> (f64, f64) {
        let (a, d) = (libm::tanh(beta * 0.5), libm::tanh(beta * 0.5) * 2.0);
        let th = libm::tanh(beta * (t - 0.5));
        ((a + th) / d, beta * (1.0 - th * th) / d)
    };
    // Per cell: its density and the projection's slope.
    let physical = |beta: f64, xs: &[f64]| -> (Vec<f64>, Vec<f64>) {
        let xt = filter(xs);
        (0..nc).map(|i| if fixed[i] { (1.0, 0.0) } else { project(beta, xt[i]) }).unzip()
    };
    // Each cell's stiffness multiplier from its density.
    let (penalty, emin) = (3.0, 1e-9);
    let stiffness = |xb: &[f64]| -> Vec<f64> { xb.iter().map(|v| emin + (1.0 - emin) * libm::pow(*v, penalty)).collect() };
    // Updates by the method of moving asymptotes (ink-mma). Optimality
    // criteria (top88's update, move 0.2, damping 1/2) settle while the
    // projection is soft, then at beta 16 flip the cells at its
    // threshold by the whole move limit every iteration for good (an
    // MBB beam's densities moved 4.6e-2 rms per iteration after 100
    // iterations at beta 16; smaller moves or damping shrink the cycle,
    // not end it): MMA's asymptotes close in on a variable that
    // oscillates. Objective and constraint are scaled so their gradients
    // are of order one per variable (the objective by its first value),
    // as its regularisation assumes, whatever the number of variables.
    // The volume is met exactly at every iterate, not only as
    // approximated.
    let mut mma = ink_mma::Mma::new(0.2);
    let mut c0: Option<f64> = None;
    // Beta doubles from 1 to 16 once the design has settled at its value
    // (no physical density moving by 0.01, as Andreassen et al. continue
    // theirs) or after a tenth of the budget there, so the last value
    // has most of it; the run stops once the design settles at 16.
    // Settled is judged on the physical densities, not the design
    // variables: where the projection saturates a design variable may
    // creep (0.03 an iteration) with no density moving.
    let (mut beta, every, mut since) = (1.0, (top.iterations / 10).max(1), 0usize);
    let mut prev: Option<Vec<f64>> = None;
    let (mut history, mut settled) = (Vec::new(), false);
    for it in 1..=top.iterations.max(1) {
        let (xb, dxu) = physical(beta, &xs);
        // The density update needs only a rough displacement field.
        let mult = stiffness(&xb);
        let objective = a.solve(&mult, tol.max(1e-4))?;
        // Element strain energy at unit stiffness (a support's penalty
        // with it).
        let energy = a.energies();
        let volume: f64 = xb.iter().zip(&vol).map(|(a, b)| a * b).sum();
        let change = prev.as_ref().map_or(1.0, |p| p.iter().zip(&xb).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max));
        history.push((it, objective, volume / total, change, beta));
        // Settled at the steepest projection, within its budget (the
        // update meets the volume at the steepness it ran at; one at a
        // steepness just doubled may not, and a cube stopped at 0.99989
        // of 0.99).
        if change < 0.01 && beta >= 16.0 && volume / total <= top.volume_fraction + VOLUME_TOL {
            settled = true;
            break;
        }
        // Sensitivities with respect to the design variables: per cell,
        // through the projection, then the filter.
        let n = nc as f64;
        let c0 = *c0.get_or_insert(objective);
        let slope: Vec<f64> = (0..nc).map(|i| -penalty * (1.0 - emin) * libm::pow(xb[i], penalty - 1.0) * energy[i]).collect();
        let design = |g: &[f64], scale: f64| -> Vec<f64> { filter_t(&g.iter().zip(&dxu).map(|(g, d)| g * d * scale).collect::<Vec<_>>()) };
        #[cfg(test)]
        PROBE.with(|p| {
            let work = mult.iter().zip(&energy).map(|(m, e)| m * e).sum::<f64>();
            p.borrow_mut().push((objective, work, design(&slope, 1.0).iter().sum::<f64>()))
        });
        let (df, dv) = (design(&slope, n / c0), design(&vol, n / total));
        let over = |y: &[f64]| physical(beta, &expand(y)).0.iter().zip(&vol).map(|(a, b)| a * b).sum::<f64>() > top.volume_fraction * total;
        let next = mma.update(&pick(&xs), &pick(&df), &[(n * (volume / total - top.volume_fraction), pick(&dv))], Some(&over));
        xs = expand(&next);
        since += 1;
        // (Not after the last: the design returned is the one this
        // update made, which meets the volume at this steepness.)
        if beta < 16.0 && (change < 0.01 || since >= every) && it < top.iterations {
            (beta, since) = (beta * 2.0, 0);
        }
        prev = Some(xb);
    }
    // Final state at the final (projected) densities: the design
    // returned, which the history's last row describes. Out of budget,
    // that is the design after the last update, a row of its own (its
    // change from the last solved), not the one before it.
    let xb = physical(beta, &xs).0;
    let objective = a.solve(&stiffness(&xb), tol)?;
    let volume = xb.iter().zip(&vol).map(|(a, b)| a * b).sum::<f64>() / total;
    // A run that settled returns the design its last row describes.
    if settled {
        let last = history.last_mut().expect("a settled run solved a design");
        (last.1, last.2) = (objective, volume);
    } else {
        let change = prev.as_ref().map_or(1.0, |p| p.iter().zip(&xb).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max));
        history.push((history.len() + 1, objective, volume, change, beta));
    }
    Ok(Record { density: xb, history, settled, fixed: held })
}

/// A solved topology: the solution at its design, what the run found, and
/// its part.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Layout {
    pub sol: sim::Solution,
    pub record: Record,
    pub part: Part,
}

/// The layout's part, as `(study-field TOPOLOGY solid)` draws it in the
/// study's shape (material where its density is over 1/2), solved as a
/// study of its own: its linear buckling analysis and the stability check
/// on it, which a redesign's density model runs none of (grey cells read at
/// density^3 buckle where no part does). Its band, where the loads may reach
/// it, is the study's own `buckles` warning; why it has none, if it has
/// none. The cells it was solved on.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct Part {
    pub buckle: Option<sim::Buckle>,
    pub failure: Option<String>,
    pub h: f64,
}

/// `Part` of `layout`, a topology of `study` set up as `setup`, solved in
/// `ctx`.
fn part(doc: &mut lang::Doc, ctx: &Ctx, study: &str, mut setup: sim::Setup, layout: &Layout) -> Part {
    use lang::Outcome;
    let failed = |why: String| Part { buckle: None, failure: Some(why), h: setup.h };
    let Some(solid) = layout.field("solid") else { return failed("the layout gave no part".into()) };
    let Some(shape) = doc.form::<sim::Study>(study).map(|st| st.shape) else { return failed(format!("no study '{study}'")) };
    let at = [doc.g.var(graph::Var::X), doc.g.var(graph::Var::Y), doc.g.var(graph::Var::Z)];
    let solid = doc.g.grid(Arc::new(solid), at);
    let drawn = doc.g.max(shape, solid);
    let drawn = doc.g.untagged(drawn);
    setup.shape = Arc::new(Tape::new(&doc.g, drawn));
    match sim::solve_cached(&setup, ctx) {
        Ok(sol) => match (sol.report.failure, sol.report.buckle) {
            (Some(why), _) => failed(why),
            (None, Some(b)) => Part { failure: b.failure.clone(), buckle: Some(b), h: sol.h },
            (None, None) => failed("its solve gave no buckling analysis".into()),
        },
        Err(why) => failed(why),
    }
}

/// A cached layout: the tolerance it was solved to, and the layout.
#[derive(serde::Serialize, serde::Deserialize)]
struct Cached {
    tol: f64,
    layout: Layout,
}

/// Runs topology `name`: its study set up (at `res`, if given) and
/// optimized. Cached by the question as the solver reads it (the study's
/// `Setup::key`, with the settings and `keep`), in memory and on disk; a
/// tighter cached tolerance serves looser requests.
pub fn run(doc: &mut lang::Doc, ctx: &Ctx, name: &str, res: Option<f64>, tol: f64) -> Result<Layout, InkError> {
    let def = doc.pick::<TopologyDef>("topology", Some(name))?;
    let (study, volume_fraction, iterations, filter, keep) = (def.study.clone(), def.volume_fraction, def.iterations, def.filter, def.keep);
    let (setup, _, _) = sim::study_setup(doc, ctx, &study, res, tol)?;
    let keep = keep.map(|k| {
        let k = doc.g.untagged(k);
        Arc::new(Tape::new(&doc.g, k))
    });
    let settings = Settings { volume_fraction, iterations, filter, keep };
    let words = format!("topology {volume_fraction:?} {iterations} {filter:?} {}", settings.keep.is_some());
    let key = setup.key_with(ctx, &words, &settings.keep.iter().map(|k| &**k).collect::<Vec<_>>());
    if let Some((k, _)) = key {
        if let Some((ctol, l)) = doc.solved.get(&k).and_then(|v| v.downcast_ref::<(f64, Layout)>())
            && *ctol <= tol
        {
            return Ok(l.clone());
        }
        if let Some(Cached { tol: ctol, layout: mut l }) = cache::load::<Cached>("topology", k)
            && ctol <= tol
        {
            l.sol.report.cached = true;
            doc.solved.insert(k, Arc::new((ctol, l.clone())));
            return Ok(l);
        }
    }
    // In the context the key projects: what it leaves out is not read.
    let at = key.as_ref().map_or(ctx, |(_, at)| at);
    let mut record = None;
    let sol = sim::solve_redesigned(&setup, at, &mut |a| {
        record = Some(optimize(a, &settings, at, setup.tol)?);
        Ok(())
    })
    .map_err(|m| InkError::msg(format!("topology '{name}': {m}")))?;
    let mut layout = Layout { sol, record: record.expect("a redesigned solve ran its redesign"), part: Part { buckle: None, failure: None, h: setup.h } };
    layout.part = part(doc, at, &study, setup, &layout);
    // A failed solve is not a result to serve again.
    if let Some((k, _)) = key
        && layout.sol.report.failure.is_none()
    {
        cache::store("topology", k, &Cached { tol, layout: layout.clone() });
        doc.solved.insert(k, Arc::new((tol, layout.clone())));
    }
    Ok(layout)
}

#[cfg(test)]
mod tests {
    use super::*;
    use sim::{Fix, Held, Load, Motion, Physics, RegionSel, Setup};

    /// A 12 x 2 x 2 beam on two pins about y (at x = -5.5, and a roller in z
    /// at x = 5.5), 1 N down over its middle, E 1000, nu 0, res 1.
    fn pinned_beam() -> (Setup, Ctx) {
        use crate::graph::{Graph, NodeId, UnOp, Var};
        let mut g = Graph::new();
        let (x, y, z) = (g.var(Var::X), g.var(Var::Y), g.var(Var::Z));
        let c = |g: &mut Graph, v: f64| g.constant(v);
        let less = |g: &mut Graph, a: NodeId, k: f64| {
            let a = g.un(UnOp::Abs, a);
            let k = g.constant(k);
            g.sub(a, k)
        };
        let (bx, by, bz) = (less(&mut g, x, 6.0), less(&mut g, y, 1.0), less(&mut g, z, 1.0));
        let shape = {
            let m = g.max(bx, by);
            g.max(m, bz)
        };
        let five = c(&mut g, 5.0);
        let left = g.add(x, five);
        let right = g.sub(five, x);
        let middle = less(&mut g, x, 0.5);
        let one = c(&mut g, 1000.0);
        let t = |n: NodeId| Arc::new(Tape::new(&g, n));
        let fix = |name: &str, region: NodeId, dofs: [bool; 3], at: f64| Fix {
            boundary: sim::Boundary::Unspecified,
            name: name.into(),
            region: RegionSel::Field(t(region)),
            dofs,
            value: [Held::Number(0.0), Held::Number(0.0), Held::Number(0.0)],
            film: None,
            normal: false,
            free: vec![Motion::turn([at, 0.0, 0.0], [0.0, 1.0, 0.0])],
            rigid: false,
        };
        let setup = Setup {
            physics: Physics::Elastic { nu: 0.0 },
            shape: t(shape),
            coef: t(one),
            coef_const: Some(1000.0),
            nu: None,
            strength: None,
            body: None,
            fixes: vec![fix("pin", left, [true; 3], -5.5), fix("roller", right, [false, false, true], 5.5)],
            loads: vec![Load { name: "down".into(), region: RegionSel::Field(t(middle)), total: [0.0, 0.0, -1.0], total_t: None, pressure: None, torque: None }],
            lo: [-6.0, -1.0, -1.0],
            hi: [6.0, 1.0, 1.0],
            h: 1.0,
            pressure_grid: 1.0,
            tol: 1e-10,
            support_margin: 2.0,
            modes: 0,
            density: None,
            transient: None,
            split_as: None,
            both_ways: false,
        };
        (setup, Ctx::new(0.0, vec![]))
    }

    #[test]
    fn slopes_take_a_pins_turn_as_no_stiffness() {
        // At every design the cells' energies (each times its stiffness)
        // sum to the compliance, and the slope summed over the design is
        // the compliance's change with every design variable at once (the
        // filter one cell: each cell's own). Summing u^T A u over a pin's
        // cells counted its free turn: 61.04 against 0.0213 N mm.
        // The probe is read on the thread the redesign ran on: the solve's.
        let run = |vf: f64| {
            let (s, ctx) = pinned_beam();
            let top = Settings { volume_fraction: vf, iterations: 1, filter: Some(1.0), keep: None };
            let mut probe = None;
            sim::solve_redesigned(&s, &ctx, &mut |a| {
                PROBE.with(|p| p.borrow_mut().clear());
                optimize(a, &top, &ctx, s.tol)?;
                probe = PROBE.with(|p| p.borrow().first().copied());
                Ok(())
            })
            .unwrap();
            probe.unwrap()
        };
        let (c, work, slope) = run(0.5);
        assert!((work / c - 1.0).abs() < 1e-6, "energies {work}, compliance {c}");
        let d = 1e-4;
        let fd = (run(0.5 + d).0 - run(0.5 - d).0) / (2.0 * d);
        assert!((slope / fd - 1.0).abs() < 1e-3, "slope {slope}, differences {fd}");
    }

    #[test]
    fn pin_energies_follow_individual_changes_to_stiffness() {
        for rigid in [false, true] {
            let (mut s, ctx) = pinned_beam();
            for f in &mut s.fixes {
                f.rigid = rigid;
            }
            sim::solve_redesigned(&s, &ctx, &mut |a| {
                let mult: Vec<f64> = (0..a.cells()).map(|i| 0.5 + 0.125 * (i % 5) as f64).collect();
                let c = a.solve(&mult, s.tol)?;
                let e = a.energies();
                let work = mult.iter().zip(&e).map(|(m, e)| m * e).sum::<f64>();
                assert!((work / c - 1.0).abs() < 1e-6, "rigid {rigid}: energies {work}, compliance {c}");
                // A nonuniform perturbation changes the distribution of a
                // hold's weight and the ratio of a lent traction's two scales.
                let delta: Vec<f64> = (0..a.cells()).map(|i| (i % 7) as f64 - 3.0).collect();
                let slope = |t: f64, a: &mut dyn sim::Assembly| a.solve(&mult.iter().zip(&delta).map(|(m, d)| m + t * d).collect::<Vec<_>>(), s.tol);
                let step = 1e-4;
                let fd = (slope(step, a)? - slope(-step, a)?) / (2.0 * step);
                let expected = -e.iter().zip(&delta).map(|(e, d)| e * d).sum::<f64>();
                assert!((fd - expected).abs() < 1e-6 * c, "rigid {rigid}: slope {expected}, differences {fd}");
                Ok(())
            })
            .unwrap();
        }
    }

    /// A small grid's solves run on the threads it keeps busy, not the
    /// caller's every thread (48 cells: one): the study's solve sizes its
    /// threads to its grid.
    #[test]
    fn a_small_grid_runs_on_few_threads() {
        let (s, ctx) = pinned_beam();
        let top = Settings { volume_fraction: 0.5, iterations: 3, filter: Some(1.0), keep: None };
        let wide = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let mut threads = 0;
        wide.install(|| {
            sim::solve_redesigned(&s, &ctx, &mut |a| {
                threads = rayon::current_num_threads();
                optimize(a, &top, &ctx, s.tol).map(|_| ())
            })
        })
        .unwrap();
        assert_eq!(threads, 1);
    }

    /// A kept region's cells come back at density 1, solid, whatever the
    /// budget left the rest (here a beam kept below x = 4, at a 0.3 budget).
    #[test]
    fn kept_cells_come_back_solid() {
        use crate::graph::{Graph, Var};
        let mut g = Graph::new();
        let x = g.var(Var::X);
        let four = g.constant(4.0);
        let end = g.sub(x, four);
        let keep = Arc::new(Tape::new(&g, end));
        let (s, ctx) = pinned_beam();
        let top = Settings { volume_fraction: 0.3, iterations: 2, filter: Some(1.0), keep: Some(keep.clone()) };
        let mut kept = vec![];
        sim::solve_redesigned(&s, &ctx, &mut |a| {
            let r = optimize(a, &top, &ctx, s.tol)?;
            let centers: Vec<V3> = (0..a.cells()).map(|i| a.cell(i).1).collect();
            let inside = eval::eval_points_par(&keep, &ctx, &centers);
            kept = (0..a.cells()).filter(|&i| inside[i] < 0.0).map(|i| r.density[i]).collect();
            Ok(())
        })
        .unwrap();
        assert!(!kept.is_empty());
        assert!(kept.iter().all(|&d| d == 1.0), "{kept:?}");
    }
}
