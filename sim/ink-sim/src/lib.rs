//! The simulation package: the `(study ...)` form, its metrics and results
//! (forms.rs), and its solves: a study read from the document into a setup
//! of the field-based finite element solver (ink-fem), solved and
//! cached by the question as the solver reads it (`Setup::key`).
//! It depends on core, kernel, language and its solver, never on another
//! package (../kernel/docs/CORE.md).

// The crates it builds on, at the paths the package's code uses.
use ink_core::{V3, eval, graph, tape};
use ink_kernel::{cache, query};
use ink_lang as lang;
use ink_lang::InkError;
use ink_lang::sexp;

mod command;
mod forms;
pub use self::forms::*;
// The solver, at the package's paths.
pub use ink_fem::*;

use crate::eval::Ctx;
use crate::tape::Tape;
use std::sync::Arc;

/// Set up and solve a study.
pub fn run_study(doc: &mut lang::Doc, ctx: &Ctx, name: &str, res: Option<f64>, tol: f64) -> Result<(Solution, StudyPhysics, Option<u32>), InkError> {
    let (setup, physics, ys) = study_setup(doc, ctx, name, res, tol)?;
    Ok((solve_setup(doc, ctx, name, &setup)?, physics, ys))
}

/// Solve study `name`'s setup, cached by the solver's question
/// (`Setup::key`: the study's inputs as it reads them, not names), in
/// memory and on disk; a tighter cached tolerance serves looser requests.
pub fn solve_setup(doc: &mut lang::Doc, ctx: &Ctx, name: &str, setup: &Setup) -> Result<Solution, InkError> {
    let tol = setup.tol;
    let key = setup.key(ctx);
    let named = |mut sol: Solution, doc: &lang::Doc| {
        name_report(&mut sol.report, doc.form::<Study>(name).unwrap());
        Ok(sol)
    };
    if let Some(sol) = solved(doc, setup, ctx) {
        // Read from disk earlier in this run, or solved in it.
        if sol.report.cached {
            cache::served()
        } else {
            cache::computed()
        }
        return named(sol.clone(), doc);
    }
    let sol = solve_cached(setup, ctx).map_err(|m| InkError::msg(format!("study '{name}': {m}")))?;
    if let Some((k, _)) = key
        && sol.report.failure.is_none()
    {
        doc.solved.insert(k, Arc::new((tol, sol.clone())));
    }
    named(sol, doc)
}

/// Solves study `name` at each of `at` (moments of one motion: one setup, its loads apart) ahead of
/// `solve_setup`, which then finds each solved: at once, as many as memory leaves room for
/// (`solves_at_once`, by what a solve on these cells this run already holds says one takes; with none,
/// the first alone says), each one solve on the threads all share. One that fails is left for
/// `solve_setup` to say.
pub fn solve_each(doc: &mut lang::Doc, at: &[(Ctx, lang::RunOpts)], name: &str) -> Result<(), InkError> {
    let mut todo: Vec<(u128, Setup, &Ctx)> = Vec::new();
    for (ctx, opts) in at {
        let (mut setup, _, _) = study_setup(doc, ctx, name, opts.res, opts.tol)?;
        (setup.split_as, setup.both_ways) = (opts.split_as, opts.both_ways);
        let Some((k, _)) = setup.key(ctx) else { continue };
        let held = solved(doc, &setup, ctx).is_some();
        if !held && todo.iter().all(|t| t.0 != k) {
            todo.push((k, setup, ctx));
        }
    }
    if todo.len() < 2 {
        return Ok(());
    }
    let keep = |doc: &mut lang::Doc, k: u128, tol: f64, sol: Result<Solution, SolveError>| -> Option<f64> {
        let sol = sol.ok().filter(|s| s.report.failure.is_none())?;
        let need = sol.report.memory_estimate;
        doc.solved.insert(k, Arc::new((tol, sol)));
        Some(need)
    };
    let (k, setup, ctx) = &todo[0];
    let like = doc
        .solved
        .values()
        .filter_map(|v| v.downcast_ref::<(f64, Solution)>())
        .find(|(_, s)| s.h == setup.h && s.bounds == (setup.lo, setup.hi))
        .map(|(_, s)| s.report.memory_estimate);
    let (need, from) = match like {
        Some(need) => (need, 0),
        None => match keep(doc, *k, setup.tol, solve_cached_checked(setup, ctx)) {
            Some(need) => (need, 1),
            None => return Ok(()),
        },
    };
    let rest: Vec<(u128, &Setup, &Ctx)> = todo[from..].iter().map(|(k, s, c)| (*k, s, *c)).collect();
    for (k, tol, sol) in at_once(solves_at_once(need), &rest, |(k, setup, ctx)| (*k, setup.tol, solve_cached_checked(setup, ctx))) {
        keep(doc, k, tol, sol);
    }
    Ok(())
}

/// A setup's solution where this run holds one, solved to its tolerance or tighter: its own, or the
/// one with its buckling analysis both ways (`Setup::both_ways`), which holds all its own does.
fn solved<'a>(doc: &'a lang::Doc, setup: &Setup, ctx: &Ctx) -> Option<&'a Solution> {
    let held = |key: Option<(u128, Ctx)>| doc.solved.get(&key?.0).and_then(|v| v.downcast_ref::<(f64, Solution)>()).filter(|(tol, _)| *tol <= setup.tol).map(|(_, sol)| sol);
    held(setup.key(ctx)).or_else(|| held(setup.key_both_ways(ctx)))
}

/// Solve a setup, cached on disk by its key (`Setup::key`) where it has
/// one. A failed solve is not a result to serve again (a looser tolerance
/// might reach it).
pub fn solve_cached(setup: &Setup, ctx: &Ctx) -> Result<Solution, String> {
    solve_cached_checked(setup, ctx).map_err(|e| e.to_string())
}

fn solve_cached_checked(setup: &Setup, ctx: &Ctx) -> Result<Solution, SolveError> {
    let tol = setup.tol;
    let key = setup.key(ctx);
    if let Some((k, _)) = key
        && let Some(Cached { tol: ctol, solution: mut sol }) = cache::load::<Cached>("studies", k)
        && ctol <= tol
    {
        sol.report.cached = true;
        cache::served();
        return Ok(sol);
    }
    // In the context the key projects: what it leaves out is not read.
    let at = key.as_ref().map_or(ctx, |(_, at)| at);
    let sol = solve_checked(setup, at)?;
    cache::computed();
    if let Some((k, _)) = key
        && sol.report.failure.is_none()
    {
        cache::store("studies", k, &Cached { tol, solution: sol.clone() });
    }
    Ok(sol)
}

/// Puts the study's support and load names on a solve's report, by
/// ordinal: a solve is keyed by content, not names, so a cached one
/// carries those of the study that first solved it.
fn name_report(r: &mut Report, st: &Study) {
    let body = r.applied.len().saturating_sub(st.loads.len()); // gravity's total comes first
    r.applied[body..].iter_mut().zip(&st.loads).for_each(|(e, l)| e.0.clone_from(&l.name));
    r.reactions.iter_mut().zip(&st.fixes).for_each(|(e, f)| e.0.clone_from(&f.name));
    let freed = st.fixes.iter().filter(|f| !f.free.is_empty());
    r.free_motions.iter_mut().zip(freed).for_each(|(e, f)| e.0.clone_from(&f.name));
}

/// The fields a study depends on: shape, material, regions, time-varying
/// loads.
pub fn study_inputs(doc: &lang::Doc, name: &str) -> Vec<u32> {
    use forms::Region;
    let Some(st) = doc.form::<Study>(name) else { return Vec::new() };
    let mut nodes = vec![st.shape];
    match st.physics {
        StudyPhysics::Elastic { e, nu, density, .. } => nodes.extend([e, nu].into_iter().chain(density)),
        StudyPhysics::Heat { k, density, .. } => nodes.extend(std::iter::once(k).chain(density)),
    }
    for r in st.fixes.iter().map(|f| &f.region).chain(st.loads.iter().map(|l| &l.region)) {
        if let Region::Field(f) = r {
            nodes.push(*f);
        }
    }
    nodes.extend(st.fixes.iter().flat_map(|f| f.value).filter(|&v| doc.g.as_const(v).is_none()));
    nodes.extend(st.loads.iter().flat_map(|l| l.vector_t.iter().flatten()));
    nodes.extend(st.loads.iter().flat_map(|l| l.from_model.iter().flatten()));
    nodes.extend(st.loads.iter().filter_map(|l| l.pressure));
    nodes
}

/// A cached study result: the tolerance it was solved to, and the solution.
#[derive(serde::Serialize, serde::Deserialize)]
struct Cached {
    tol: f64,
    solution: Solution,
}

/// Build a solver setup from a `(study ...)` in a document.
/// Units in the IR: mm, N, MPa; density kg/m^3 with gravity m/s^2;
/// conductivity W/(m K), power W. Returns the setup, the study's physics,
/// and its yield strength (if any).
pub fn study_setup(doc: &mut lang::Doc, ctx: &Ctx, name: &str, res: Option<f64>, tol: f64) -> Result<(Setup, StudyPhysics, Option<u32>), InkError> {
    let st = doc.pick::<Study>("study", Some(name))?;
    let (shape, physics, res0, gravity, margin, modes, transient) = (st.shape, st.physics, st.res, st.gravity, st.support_margin, st.modes, st.transient);
    let fixes: Vec<(String, forms::Region, [bool; 3], [u32; 3], Option<f64>, bool, Vec<Motion>, bool, Boundary)> =
        st.fixes.iter().map(|f| (f.name.clone(), f.region.clone(), f.dofs, f.value, f.film, f.normal, f.free.clone(), f.rigid, f.boundary)).collect();
    let loads: Vec<(String, forms::Region, [f64; 3], Option<[u32; 3]>, Option<u32>, Option<[f64; 3]>)> =
        st.loads.iter().map(|l| (l.name.clone(), l.region.clone(), l.vector, l.vector_t, l.pressure, l.torque)).collect();
    // What a load reads of the model, at the moment it is run in (`ModelView::at_worst`): its force, or its torque.
    let modelled: Vec<Option<[f64; 3]>> = st.loads.iter().map(|l| l.from_model.map(|nodes| nodes.map(|n| doc.scalar(n, ctx)))).collect();
    if let Some(l) = st.loads.iter().zip(&modelled).find(|(_, m)| m.is_some_and(|m| m.iter().any(|x| !x.is_finite()))).map(|(l, _)| &l.name) {
        return Err(InkError::msg(format!(
            "study '{name}': load '{l}' reads what the model works out, and has no value here: a study that reads it is asked through `ink require` or `ink solve`, which run it at the motion's worst moment"
        )));
    }
    let loads: Vec<_> = loads
        .into_iter()
        .zip(modelled)
        .map(|((n, r, v, vt, pr, tq), m)| match (m, tq) {
            (Some(m), Some(_)) => (n, r, v, vt, pr, Some(m)),
            (Some(m), None) => (n, r, m, vt, pr, tq),
            (None, _) => (n, r, v, vt, pr, tq),
        })
        .collect();
    st.check_load_clock(doc, ctx)?;
    // All of the shape is solved, so all of it must be in reach: one not proven bounded has material the
    // window may leave out (a solid no fix holds, off beyond it), and its study would answer for the rest.
    let (bounded, reach) = doc.clip_reach(shape, ctx);
    if !reach.whole {
        return Err(InkError::msg(format!(
            "study '{name}': its shape is unbounded or not proven bounded (a half-space, an unclipped lattice, or a field whose bounds do not close far away), so the part it would solve is \
             only what lies in the window: intersect it with a containing shape first, e.g. (max (box 40 40 40) f)"
        )));
    }
    let (lo, hi, _) = query::bounds(&Arc::new(Tape::new(&doc.g, bounded)), ctx, reach.lo, reach.hi).ok_or_else(|| InkError::msg(format!("study '{name}': shape is empty")))?;
    let pressure_grid = res0.unwrap_or(query::max_extent(lo, hi) / 100.0);
    let h = res.unwrap_or(pressure_grid);
    let pad = 2.0 * if loads.iter().any(|l| l.4.is_some()) { pressure_grid } else { h };
    let clipped = doc.clip(bounded, [lo[0] - pad, lo[1] - pad, lo[2] - pad], [hi[0] + pad, hi[1] + pad, hi[2] + pad]);
    // The solver reads names only through a face region's provenance on the
    // shape: every other field it is given without them, so a field a def
    // names is the one written in place (`Setup::key` keys them alike).
    let faces = fixes.iter().map(|f| &f.1).chain(loads.iter().map(|l| &l.1)).any(|r| matches!(r, forms::Region::Face { .. }));
    let tape = |d: &mut lang::Doc, f: u32| {
        let f = d.g.untagged(f);
        Arc::new(Tape::new(&d.g, f))
    };
    let shape_tape = if faces { Arc::new(Tape::new(&doc.g, clipped)) } else { tape(doc, clipped) };
    let region_sel = |d: &mut lang::Doc, r: &forms::Region| match r {
        forms::Region::Field(f) => {
            let f = d.g.untagged(*f);
            let f = pieces_as_distances(d, ctx, f, (lo, hi));
            RegionSel::Field(tape(d, f))
        }
        forms::Region::Face { tag, normal, depth } => RegionSel::Face {
            shape: shape_tape.clone(),
            tag: d.g.tag_names.iter().position(|n| n == tag).expect("tag checked at parse") as u32,
            normal: *normal,
            depth: depth.unwrap_or(1.5 * h),
        },
    };
    let (sim_physics, coef, ys) = match physics {
        // A constant nu is the physics' own; a field is sampled per cell.
        StudyPhysics::Elastic { e, nu, yield_strength, .. } => (Physics::Elastic { nu: doc.g.as_const(nu).unwrap_or(0.3) }, e, yield_strength),
        StudyPhysics::Heat { k, .. } => {
            let scale = doc.g.constant(1e-3); // W/(m K) -> W/(mm K)
            (Physics::Heat, doc.g.mul(k, scale), None)
        }
    };
    let body = match (physics, gravity) {
        // kg/m^3 * m/s^2 * 1e-9 m^3/mm^3 -> N/mm^3
        (StudyPhysics::Elastic { density: Some(rho), .. }, Some(g)) => Some((tape(doc, rho), g.map(|c| c * 1e-9))),
        _ => None,
    };
    let setup = Setup {
        physics: sim_physics,
        shape: shape_tape.clone(),
        coef: tape(doc, coef),
        coef_const: doc.g.as_const(coef),
        nu: match physics {
            StudyPhysics::Elastic { nu, .. } if doc.g.as_const(nu).is_none() => Some(tape(doc, nu)),
            _ => None,
        },
        strength: ys.filter(|&y| doc.g.as_const(y).is_none()).map(|y| tape(doc, y)),
        body,
        fixes: fixes
            .iter()
            .map(|(n, r, d, v, film, normal, free, rigid, boundary)| Fix {
                name: n.clone(),
                region: region_sel(doc, r),
                dofs: *d,
                value: v.map(|f| match doc.g.as_const(f) {
                    Some(x) => Held::Number(x),
                    None => Held::Field(tape(doc, f)),
                }),
                film: film.map(|h| h * 1e-6),
                boundary: *boundary,
                normal: *normal,
                free: free.clone(),
                rigid: *rigid,
            })
            .collect(),
        loads: loads
            .iter()
            .map(|(n, r, v, vt, pr, tq)| Load {
                name: n.clone(),
                region: region_sel(doc, r),
                total: *v,
                total_t: vt.map(|c| c.map(|f| tape(doc, f))),
                pressure: pr.map(|f| tape(doc, f)),
                torque: *tq,
            })
            .collect(),
        lo,
        hi,
        h,
        pressure_grid,
        tol,
        support_margin: margin.unwrap_or(SUPPORT_MARGIN_MM),
        modes: modes.unwrap_or(0),
        density: match physics {
            StudyPhysics::Elastic { density: Some(rho), .. } | StudyPhysics::Heat { density: Some(rho), .. } => Some(tape(doc, rho)),
            _ => None,
        },
        transient,
        split_as: None,
        both_ways: false,
    };
    Ok((setup, physics, ys))
}

/// A region's field with each piece of its maxes and mins (a corner's
/// faces, a box's sides, a bore) at the scale of a distance, so its faces
/// are traced by which is nearest, not by which field is largest (a face
/// written a thousand times itself went untraced where it was the nearer,
/// and a hold's corner moved 0.137 for 0.103; a bore's field a thousandth
/// of itself, 0.0419 for 0.0311). Each piece is divided by a positive
/// constant, so the region and what its bounds prove are the same: a flat
/// piece by its gradient's size, a curved one by its gradient's size where
/// it reads nearest its zero among points across the study's box `span`
/// (a distance times any constant, that constant; an ellipsoid times its
/// least axis, near one). Its gradient's size as a field would loosen its
/// bounds; a piece whose gradient is nowhere read is left as it is.
fn pieces_as_distances(d: &mut lang::Doc, ctx: &Ctx, f: u32, span: (V3, V3)) -> u32 {
    use graph::{BinOp, Node, Var};
    let piece = |d: &mut lang::Doc, f: u32| -> u32 {
        if matches!(d.g.nodes[f as usize], Node::Bin(BinOp::Max | BinOp::Min | BinOp::Emax, ..)) {
            return pieces_as_distances(d, ctx, f, span);
        }
        let grad: Option<Vec<f64>> = [Var::X, Var::Y, Var::Z]
            .into_iter()
            .map(|v| {
                let g = d.derivative(f, v);
                d.g.as_const(g)
            })
            .collect();
        let n = match grad {
            Some(g) => size([g[0], g[1], g[2]]),
            None => curved_scale(&Tape::new(&d.g, f), ctx, span),
        };
        // (Not the model's number, so not noted as one out of range; and
        // over the least normal, so its quotient's gradient stays finite.)
        if n >= f64::MIN_POSITIVE && n.is_finite() && n != 1.0 {
            let n = d.g.own(n);
            d.g.bin(BinOp::Div, f, n)
        } else {
            f
        }
    };
    match d.g.nodes[f as usize] {
        Node::Bin(op @ (BinOp::Max | BinOp::Min | BinOp::Emax), a, b) => {
            let (na, nb) = (piece(d, a), piece(d, b));
            if (na, nb) == (a, b) { f } else { d.g.bin(op, na, nb) }
        }
        _ => f,
    }
}

/// A curved field's gradient's size where it reads nearest its zero, by
/// its value over that size, among 9^3 points across `span` (0 where no
/// point has a gradient: a constant, or one read only where it has none).
/// Within a millionth of one it is one: a distance written plainly is left
/// as it is to the last digit.
fn curved_scale(tape: &Tape, ctx: &Ctx, (lo, hi): (V3, V3)) -> f64 {
    const N: usize = 9;
    let pts: Vec<V3> = (0..N * N * N).map(|i| std::array::from_fn(|a| lo[a] + (hi[a] - lo[a]) * (([i % N, i / N % N, i / (N * N)][a]) as f64 + 0.5) / N as f64)).collect();
    let near = eval::eval_grads(tape, ctx, &pts)
        .into_iter()
        .map(|g| (g[0], size([g[1], g[2], g[3]])))
        .filter(|&(v, n)| v.is_finite() && n >= f64::MIN_POSITIVE && n.is_finite())
        .min_by(|a, b| (a.0 / a.1).abs().total_cmp(&(b.0 / b.1).abs()));
    match near {
        Some((_, n)) if (n - 1.0).abs() > 1e-6 => n,
        Some(_) => 1.0,
        None => 0.0,
    }
}

/// A gradient's size, read over its largest part so that neither its
/// squares' overflow (past 1e154) nor their underflow (under 1e-154) loses it.
fn size(g: V3) -> f64 {
    let m = g.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    if m == 0.0 || !m.is_finite() { m } else { m * g.iter().map(|x| (x / m).powi(2)).sum::<f64>().sqrt() }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A curved piece's scale is its gradient's size near its zero: a
    /// distance times a constant reads that constant, a distance written
    /// plainly one exactly, and a piece with a gradient nowhere it is read
    /// none (left as it is).
    #[test]
    fn a_curved_pieces_scale_is_its_gradients_size_near_its_zero() {
        static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&crate::Sim]);
        let src = "(import \"std\")
(def thin (* 0.001 (- (len2 y z) 1.5)))
(def plain (- (len2 y z) 1.5))
(def wide (* 1000 (- (len2 (len2 x y) z) 1)))
(def level (max 0 (- (len2 y z) 100)))
(def tiny (* 1e-100 (* 1e-100 (- (len2 y z) 1.5))))
(def huge (* 1e100 (* 1e100 (- (len2 y z) 1.5))))
(def subnormal (* 1e-100 (* 1e-100 (* 1e-100 (* 1e-10 (- (len2 y z) 1.5))))))";
        let doc = ink_lang::Doc::load_str_in(&REGISTRY, src, "scales", std::path::Path::new("."), &[]).unwrap();
        let ctx = doc.ctx(&[], 0.0).unwrap();
        let span = ([-2.0, -2.0, -2.0], [2.0, 2.0, 2.0]);
        for (name, want) in [("thin", 0.001), ("plain", 1.0), ("wide", 1000.0), ("level", 0.0), ("tiny", 1e-200), ("huge", 1e200), ("subnormal", 0.0)] {
            let ink_lang::Value::Field(f, _) = doc.def(name).unwrap().clone() else { panic!("{name}: a field") };
            let got = curved_scale(&Tape::new(&doc.g, f), &ctx, span);
            assert!((got - want).abs() <= 1e-12 * want, "{name}: {got} for {want}");
        }
        for (g, want) in [([3e200, 4e200, 0.0], 5e200), ([3e-200, 0.0, 4e-200], 5e-200)] {
            assert!((size(g) / want - 1.0).abs() < 1e-15, "{g:?}: {}", size(g));
        }
        assert_eq!(size([0.0; 3]), 0.0);
    }
}
