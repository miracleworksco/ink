//! A transient heat study (`transient`): the temperatures stepped over
//! its duration, on three time grids, and their peak.

use crate::Held;
use crate::cells::place;
use crate::element::shape_n;
use crate::eval::Ctx;
use crate::multigrid::{Multigrid, pcg};
use crate::occupancy::at_xi;
use crate::system::System;
use crate::tape::Tape;
use crate::transient::{INTEGRAL_TOL, MAX_PIECES, integral_in_t, jumps_in_t};
use crate::{GRAIN, MAX_ITER, Physics};
use ink_core::{V3, interval};
use ink_kernel::mesh::find;
use rayon::prelude::*;

use super::SolveError;
use super::assemble::{Mesh, Model, assemble};
use super::connect::Pieces;
use super::held::held_cells;
use super::holds::{Surfaces, hold_fixes};
use super::loads::Loads;
use super::readout::SurfaceReads;
use super::run::{Rhs, Run, lift as lift_of};
use super::study::Study;
use crate::Transient;

#[allow(clippy::too_many_arguments)]
pub(super) fn transient(
    model: &Model,
    tr: Transient,
    sys: &mut System,
    rhs: &mut Rhs,
    loads: &mut Loads,
    joined: &mut Pieces,
    warnings: &mut Vec<String>,
    node_reliable: &[bool],
    lift: &[Option<f64>],
    reads: &SurfaceReads,
    surf: &Surfaces,
    run: &mut Run,
) -> Result<(), SolveError> {
    let st = model.st;
    let Study { s, ctx, ndof, ref grid, .. } = *st;
    let Mesh { nn, ref remap, ref final_idx, ref final_pos, ref held_nodes, .. } = *model.mesh;
    let Rhs { ref mut f, ref f_held, ref mut f_loads } = *rhs;
    let Loads { ref mut applied, timed: ref timed_loads, .. } = *loads;
    let Pieces { ref mut insulated, ref mut parent, ref keep, .. } = *joined;
    let SurfaceReads { ref own, ref held, ref held_from, .. } = *reads;
    let surface_points = own;
    let Run { ref mut x, ref mut iters, ref mut res, ref mut history, ref mut peak_time, ref mut peak, ref mut time_levels } = *run;
    if s.physics != Physics::Heat {
        return Err("a transient study is a heat study".into());
    }
    let rho = s.density.as_ref().ok_or("a transient study needs a density")?;
    // Heat capacity (J/K) per node.
    let parts: Vec<([f64; 8], f64)> = sys.cells.par_iter().map(|c| model.sample_cell(c, None, Some(rho))).collect();
    let share = model.held_share(&sys.cells);
    let mut m = vec![0.0; nn];
    for (ci, (c, (nodal, _))) in sys.cells.iter().zip(&parts).enumerate() {
        let held = share.get(&ci);
        for k in 0..8 {
            m[c.nodes[k] as usize] += nodal[k] * 1e-9 * tr.capacity * (1.0 - held.map_or(0.0, |h| h[k][0]));
        }
    }
    // The numerical mass floor below must not give an insulated,
    // massless piece an invented heat capacity. Its constant temperature
    // mode needs some real capacity on that same connected piece.
    if !insulated.is_empty() {
        for (old, &new) in remap.iter().enumerate() {
            if new != u32::MAX && m[new as usize] > 0.0 {
                insulated.remove(&find(parent, old as u32));
            }
        }
        if !insulated.is_empty() {
            return Err(
                "an insulated piece has no positive heat capacity: give it a density and specific heat capacity, or connect it to a temperature fix or convection surface".into()
            );
        }
    }
    sys.fold(&mut m, 1);
    let floor = 1e-8 * m.iter().sum::<f64>() / nn as f64;
    let md: Vec<f64> = (0..sys.ndofs).map(|i| m[i / ndof].max(floor)).collect();
    // Each step starts close to its answer; this is ample for time stepping.
    let tol = s.tol.max(1e-6);
    // The loads over a step from a to b: each varying total's mean
    // over it (`integral_in_t`), so a step takes all the heat its
    // sources give in it, a pulse shorter than the step too (read at
    // its end, a 1 s pulse inside a 6 s step put no heat in).
    let mut spread = 0.0f64;
    let mut force_over = |a: f64, b: f64| -> Result<Vec<f64>, String> {
        let mut out = f.clone();
        for (li, unit) in timed_loads.iter() {
            let tt = s.loads[*li].total_t.as_ref().unwrap();
            let mut tot = [0.0; 3];
            for d in 0..ndof.min(3) {
                let (v, sp) = integral_in_t(&tt[d], ctx, a, b).map_err(|why| format!("{} '{}': {why}", if ndof == 1 { "source" } else { "load" }, s.loads[*li].name))?;
                tot[d] = v / (b - a);
                spread = spread.max(sp);
            }
            for (n, u) in unit.iter().enumerate() {
                for d in 0..ndof {
                    out[n * ndof + d] += tot[d] * u;
                }
            }
        }
        Ok(out)
    };
    let force_at = |t: f64| -> Vec<f64> {
        let mut out = f.clone();
        let c = Ctx { t, ..ctx.clone() };
        for (li, unit) in timed_loads.iter() {
            let tt = s.loads[*li].total_t.as_ref().unwrap();
            let tot: [f64; 3] = std::array::from_fn(|d| crate::eval::eval_point(&tt[d], &c, [0.0; 3]));
            for (n, u) in unit.iter().enumerate() {
                for d in 0..ndof {
                    out[n * ndof + d] += tot[d] * u;
                }
            }
        }
        out
    };
    // A held temperature that varies in time (`Study::timed`) is held as
    // it is at each step's end, as backward Euler reads all else there:
    // its holds' terms made again at that moment (`held_at`), their
    // right-hand side in place of the setup's, which is the last
    // moment's, the state the study reports. Read once, at the start, a
    // bar's end held at 100 sin(4.5 deg/s t) degC held 0 throughout, and
    // its middle read 0 at 32 s for 36.6 (NAFEMS T3).
    let timed = st.timed.contains(&true);
    let held_at = |t: f64| -> Result<(Vec<f64>, Vec<Option<f64>>), SolveError> {
        st.hold_at(t);
        // (the surface's meshes as they are, its measures the setup's)
        let mut again = Surfaces { pressure: None, measures: Vec::new(), split: surf.split.clone(), facets: surf.facets.clone(), samples: surf.samples.clone() };
        let mut holds = hold_fixes(st, model.cells, &mut again)?;
        let mut held = held_cells(st, model.cells, &mut holds)?;
        let (sys_t, mut f_t, mesh_t) = assemble(st, model.cells, &holds, &mut held, keep);
        sys_t.fold(&mut f_t, ndof);
        if f_t.len() != f_held.len() {
            return Err("a held temperature made again at a step numbered the nodes otherwise".into());
        }
        f_t.iter_mut().zip(f_held).for_each(|(a, b)| *a -= b);
        let lifted = lift_of(&Model { st, cells: model.cells, holds: &holds, lent: &held.lent, mesh: &mesh_t }, &sys_t);
        Ok((f_t, lifted))
    };
    // The peak temperature over reliable nodes (all, if none are), as
    // `Solution::max_nodal`, and its node.
    let any_reliable = node_reliable.iter().any(|&r| r);
    // (Inside held material, the held value: `held_nodes`; at the part's
    // surface points in a fix's region, its value there: each at the
    // moment `held_t`.)
    let held_reads = || -> (Vec<Option<f64>>, Vec<Option<f64>>) {
        let mut read = vec![None; nn];
        let mut surface = vec![None; held.len()];
        for (fi, fx) in s.fixes.iter().enumerate() {
            let nodes: Vec<u32> = held_nodes.iter().filter(|h| h.1 == fi).map(|h| h.0).collect();
            let at: Vec<V3> = nodes.iter().map(|&n| final_pos[n as usize]).collect();
            nodes.iter().zip(st.values_at(fx, &at)).for_each(|(&n, v)| read[n as usize] = Some(v[0]));
            let on: Vec<(usize, V3)> = held_from.iter().enumerate().filter_map(|(i, h)| h.filter(|h| h.0 == fi).map(|h| (i, h.1))).collect();
            let at: Vec<V3> = on.iter().map(|h| h.1).collect();
            on.iter().zip(st.values_at(fx, &at)).for_each(|(&(i, _), v)| surface[i] = Some(v[0]));
        }
        (read, surface)
    };
    let fixed = (!timed).then(|| {
        let mut read = vec![None; nn];
        held_nodes.iter().for_each(|&(n, _, g)| read[n as usize] = Some(g[0]));
        (read, held.iter().map(|h| h.map(|h| h.1[0])).collect::<Vec<_>>())
    });
    // And at the part's surface points through each one's cell's field,
    // as `max_nodal`: a slab heated at its faces read its last node
    // inside, 66.48 C (65.75 to 67.20) where the faces reach 68.04.
    let hottest = |u: &[f64], (read, surface_held): &(Vec<Option<f64>>, Vec<Option<f64>>)| -> (f64, V3) {
        let at_nodes = (0..nn)
            .filter(|&n| !any_reliable || node_reliable[n])
            .map(|n| (read[n].unwrap_or(u[n]), final_pos[n]))
            .fold((f64::NEG_INFINITY, [0.0; 3]), |a, b| if b.0 > a.0 { b } else { a });
        surface_points
            .iter()
            .zip(surface_held)
            .map(|(&(ci, xi), held)| {
                let c = &sys.cells[ci as usize];
                // (in the cell it lies in, as `max_nodal`; in a fix's
                // region, its held value)
                let xi = xi.map(|v| v.clamp(-2.0, 2.0));
                let (ctr, size) = place(grid, c);
                let v = held.unwrap_or_else(|| (0..8).map(|k| shape_n(k, xi) * u[c.nodes[k] as usize]).sum::<f64>());
                (v, at_xi(ctr, size, xi))
            })
            .filter(|x| x.0.is_finite())
            .fold(at_nodes, |a, b| if b.0 > a.0 { b } else { a })
    };
    // (K + C/dt) T+ = f + C/dt T, f the loads' mean over the step, from
    // the initial temperature with the held values lifted onto their
    // nodes (they hold from t = 0).
    // Each step solves for its change, T+ - T, from zero, so its
    // residual is relative to what drives that step (f - K T), not
    // to the right-hand side, which the held values' penalty dominates
    // (beta A v against the capacity's C/dt T grows as 1/h^2): measured
    // that way, a slab with a held face at 0.125 mm cells took no
    // iterations on most steps and read 74% low.
    // The time check: the same on steps twice and four times as long
    // first (`time_levels`), as many as halving the count allows, then
    // the study's own, whose state is reported. Backward Euler is first
    // order in its step, so their peaks estimate its error in time as
    // the claims' grids do in space (ink-sim's `time_error`).
    // Every run's steps end at each moment a source may jump: a pulse
    // ending inside a step is read at the step's ends, short of its
    // peak, on all three runs alike, and their differences miss it (a
    // copper cube heated for 3.3 s, on 23, 12 and 6 steps of 10 s:
    // 24.52, its band to 24.65, against 24.72).
    let sources: Vec<&Tape> = timed_loads.iter().flat_map(|(li, _)| s.loads[*li].total_t.as_ref().unwrap().iter().map(|t| &**t)).collect();
    let most = tr.steps.max(1);
    // (a jump at the end is no break: the last step ends there anyway)
    let inside = |sp: &(f64, f64)| sp.1 < tr.duration;
    let Some(source_jumps) = jumps_in_t(&sources, ctx, [interval::I::point(0.0); 3], tr.duration, most) else {
        // (One with no value somewhere is refused for that, by its integral.)
        force_over(0.0, tr.duration)?;
        return Err(format!("a source may switch more often than its {most} steps: take at least as many (steps n) as the times its power can jump").into());
    };
    let mut breaks: Vec<f64> = source_jumps.iter().filter(|sp| inside(sp)).map(|&sp| ink_core::eval::span_mid(sp)).collect();
    // So too where a held temperature may jump, over the cells its fix
    // holds (over the part, a hold at its far end switching at t = x mm/s
    // was refused, its switch at every moment the part spans): a step
    // ending at the jump holds the value before it, the next the one after
    // (`held_before`); the last too where it jumps at the end.
    let mut held_jumps: Vec<(f64, f64)> = Vec::new();
    for (fi, fx) in s.fixes.iter().enumerate().filter(|x| st.timed[x.0]) {
        let cells = model.holds.clamps.iter().filter(|c| c.fix == fi).map(|c| place(grid, &model.cells.classified[c.cell]));
        let Some(over) = cells
            .map(|(ctr, size)| std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * size, ctr[a] + 0.5 * size)))
            .reduce(|a: [interval::I; 3], b| std::array::from_fn(|k| interval::I::new(a[k].lo.min(b[k].lo), a[k].hi.max(b[k].hi))))
        else {
            continue;
        };
        let tapes: Vec<&Tape> = fx.value.iter().filter_map(|v| if let Held::Field(t) = v { Some(&**t) } else { None }).collect();
        let Some(found) = jumps_in_t(&tapes, ctx, over, tr.duration, most) else {
            return Err(format!("held temperature '{}' may switch more often than its {most} steps: take at least as many (steps n) as the times it can jump", fx.name).into());
        };
        held_jumps.extend(found);
    }
    breaks.extend(held_jumps.iter().filter(|sp| inside(sp)).map(|&sp| ink_core::eval::span_mid(sp)));
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    if breaks.len() > most {
        return Err(format!("its sources and held temperatures may switch more often than its {most} steps: take at least as many (steps n) as the times they can jump").into());
    }
    // A step from a ending at a held temperature's jump (its break, or the
    // end) reads it before the jump: the float under the jump's span, no
    // earlier than a. (Stepped back by a share of the duration, a pulse of
    // 1e-10 s at 2 s in a run of 1e6 s was read before its start, and
    // held nothing.)
    let held_before = |a: f64, t: f64| held_jumps.iter().find(|sp| sp.0 < t && t <= sp.1).map_or(t, |sp| sp.0.next_down().max(a));
    let mut levels = vec![tr.steps.max(1)];
    while levels.len() < 3 && levels[levels.len() - 1] > 1 {
        levels.push(levels[levels.len() - 1].div_ceil(2));
    }
    levels.reverse();
    let mut ax = vec![0.0; sys.ndofs];
    for (li, &steps) in levels.iter().enumerate() {
        let own = li + 1 == levels.len();
        let start = if timed { held_at(0.0)?.1 } else { lift.to_vec() };
        *x = start.iter().map(|v| v.unwrap_or(tr.initial)).collect();
        sys.expand(x, ndof);
        let now_reads;
        let (v0, n0) = hottest(
            x,
            match &fixed {
                Some(r) => r,
                None => {
                    now_reads = held_reads();
                    &now_reads
                }
            },
        );
        let mut best = (0.0, v0, n0);
        if own {
            history.push((0.0, v0));
        }
        // Steps of duration / steps laid back from each jump and from the
        // end, a shorter one left over after each jump: the step a peak
        // at a jump is read by halves exactly from run to run.
        let dt = tr.duration / steps as f64;
        let mut ends = vec![0.0];
        for w in std::iter::once(0.0).chain(breaks.iter().copied()).chain([tr.duration]).collect::<Vec<_>>().windows(2) {
            // (A whole number of steps, to its quotient's rounding, leaves no sliver over.)
            let k = ((w[1] - w[0]) / dt * (1.0 - 4.0 * f64::EPSILON)).ceil().max(1.0) as usize;
            ends.extend((0..k).rev().map(|j| w[1] - j as f64 * dt).filter(|&e| e > w[0]));
        }
        // A multigrid for the steps of dt, and one for the last other.
        sys.shift = md.iter().map(|c| c / dt).collect();
        let even = Multigrid::build(sys, final_idx);
        let (mut longest, mut other): (f64, Option<(f64, Multigrid)>) = (0.0, None);
        for w in ends.windows(2) {
            let (a, t) = (w[0], w[1]);
            longest = longest.max(t - a);
            let shift: Vec<f64> = md.iter().map(|c| c / (t - a)).collect();
            sys.shift = shift.clone();
            if t - a != dt && other.as_ref().is_none_or(|m| m.0 != t - a) {
                other = Some((t - a, Multigrid::build(sys, final_idx)));
            }
            let mg = if t - a == dt { &even } else { &other.as_ref().unwrap().1 };
            let precond = |r: &[f64]| mg.vcycle(0, sys, r);
            // (K + C/dt) dT = f + C/dt T - (K + C/dt) T
            let mut b = force_over(a, t)?;
            let now = held_before(a, t);
            if timed {
                b.iter_mut().zip(held_at(now)?.0).for_each(|(b, d)| *b += d);
            }
            sys.apply(x, &mut ax);
            b.par_iter_mut().with_min_len(GRAIN).zip(&shift).zip(&*x).zip(&ax).for_each(|(((b, c), x), a)| *b += c * x - a);
            let mut dx = vec![0.0; sys.ndofs];
            let (it, r) = pcg(sys, &[], &precond, &b, &mut dx, tol, MAX_ITER);
            (*iters, *res) = (*iters + it, res.max(r));
            x.par_iter_mut().with_min_len(GRAIN).zip(&dx).for_each(|(x, d)| *x += d);
            sys.expand(x, ndof);
            let now_reads;
            let (v, n) = hottest(
                x,
                match &fixed {
                    Some(r) => r,
                    None => {
                        now_reads = held_reads();
                        &now_reads
                    }
                },
            );
            if v > best.1 {
                best = (t, v, n);
            }
            if own {
                history.push((t, v));
            }
        }
        time_levels.push((ends.len() - 1, best.1, longest));
        if own {
            *peak = Some((best.0, best.1, best.2));
        }
    }
    if spread > INTEGRAL_TOL {
        warnings.push(format!(
            "a source's heat over a time step is known only to {:.2}% (its power bounded over pieces of the step, {MAX_PIECES} at most): it varies too sharply to bound closer",
            100.0 * spread
        ));
    }
    sys.shift = Vec::new();
    st.hold_at(tr.duration);
    // The final state is reported (the peak temperature is in the
    // history).
    let tp = tr.duration;
    *f = force_at(tp);
    *f_loads = f.iter().zip(f_held).map(|(a, b)| a - b).collect();
    for (li, _) in timed_loads.iter() {
        let tt = s.loads[*li].total_t.as_ref().unwrap();
        let c = Ctx { t: tp, ..ctx.clone() };
        let tot: [f64; 3] = std::array::from_fn(|d| crate::eval::eval_point(&tt[d], &c, [0.0; 3]));
        if let Some(e) = applied.iter_mut().find(|e| e.0 == s.loads[*li].name) {
            e.1 = tot;
        }
    }
    *peak_time = Some(tp);
    Ok(())
}
