//! Step 5: the pieces the cells join into (`connect`): an anchor on each, the
//! cells kept, and what leaves the study incomplete.

use crate::cells::{K, MIN_FRACTION, face_neighbours, place};
use crate::occupancy::at_xi;
use crate::report::Incomplete;
use crate::surface::{KeptPiece, enclosed_apart};
use crate::{Boundary, RegionSel};
use ink_core::{V3, interval};
use ink_kernel::mesh::find;

use super::held::HeldCells;
use super::holds::Holds;
use super::settle::Cells;
use super::study::Study;

/// The cells' pieces: which cells are kept (on a piece something holds),
/// what was dropped, and what the study cannot settle or must say.
pub(super) struct Pieces {
    /// The nodes joined into pieces (`find`).
    pub(super) parent: Vec<u32>,
    /// Pieces nothing holds, solved as insulated (transient heat).
    pub(super) insulated: std::collections::HashSet<u32>,
    pub(super) keep: Vec<bool>,
    pub(super) dropped_cells: usize,
    pub(super) dropped_volume: f64,
    pub(super) dropped_at: Option<V3>,
    pub(super) incomplete: Vec<Incomplete>,
    /// Each fix that duplicates a mirror's hold on its cut.
    pub(super) duplicate_cut: Vec<bool>,
    /// What the study says of lips nothing acts on (`thin_lip`).
    pub(super) idle: Vec<String>,
    pub(super) insulated_notes: Vec<String>,
}

#[allow(clippy::type_complexity)]
pub(super) fn connect(st: &Study, cells: &Cells, holds: &Holds, held: &HeldCells, missed_at: Option<([i64; 3], V3, String)>) -> Pieces {
    let Study { s, ctx, ndof, ref mirrors, ref grid, ref xis, .. } = *st;
    let Cells { ref classified, ref cell_nodes, ref node_idx, ref hang, .. } = *cells;
    let slivers_left = &cells.slivers.left;
    let Holds { ref clamps, ref film_points, .. } = *holds;
    let HeldCells { ref free_material, shared_hold, .. } = *held;
    let (thin_lip, absorbed_note) = (held.thin_lip.clone(), held.absorbed_note.clone());
    let classified_at = |idx: [i64; 3]| cells.classified_at(st, idx);

    // Test the actual intersection with each mirror's material trace,
    // including a fix covering only part of it. A volume hold extruded
    // across the cut (e.g. an end clamp) stays a physical end restraint.
    // A normal contact, in contrast, holds the surface itself.
    let mut duplicate_cut = vec![false; s.fixes.len()];
    for &(axis, at) in mirrors {
        let cut_cells: std::collections::HashSet<_> = clamps.iter().filter(|cl| s.fixes[cl.fix].boundary == Boundary::Mirror { axis, at }).map(|cl| cl.cell).collect();
        for cl in clamps.iter().filter(|cl| cut_cells.contains(&cl.cell)) {
            let fx = &s.fixes[cl.fix];
            if duplicate_cut[cl.fix] || fx.boundary != Boundary::Unspecified || !(fx.normal || fx.dofs[axis]) {
                continue;
            }
            let (ctr, size) = place(grid, &classified[cl.cell]);
            let mut face = std::array::from_fn(|a| interval::I::new(ctr[a] - 0.5 * size, ctr[a] + 0.5 * size));
            face[axis] = interval::I::point(at);
            duplicate_cut[cl.fix] = match &fx.region {
                RegionSel::Field(t) => (fx.normal || t.space_read()[axis]) && fx.region.touches_plane(ctx, face),
                RegionSel::Face { normal, .. } => normal.is_none_or(|n| n[axis].abs() >= libm::cos(30f64.to_radians())) && fx.region.touches_plane(ctx, face),
            };
        }
    }

    // 5. A steady solution needs an anchor on each piece. In transient heat,
    // initial temperature and heat capacity determine an insulated piece's
    // temperature too. Keep it before distributing any source over the part.
    // A modal study with no fix keeps every piece: each is free, its rigid
    // motions projected out of the modes (`rigid`). One with fixes drops a
    // piece none holds, as a static study does, and its claims are not
    // settled: a part cut from its supports by a gap would otherwise ring
    // at its free modes, far over its supported ones.
    let keep_free = s.transient.is_some() || s.modes > 0 && s.fixes.is_empty();
    let mut parent: Vec<u32> = (0..node_idx.len() as u32).collect();
    // A hanging node joins what it hangs on.
    let hung = hang.iter().map(|(n, m)| std::iter::once(*n).chain(m.iter().map(|x| x.0)).collect::<Vec<u32>>());
    for nodes in cell_nodes.iter().map(|n| n.to_vec()).chain(hung) {
        let r0 = find(&mut parent, nodes[0]);
        for &n in &nodes[1..] {
            let r = find(&mut parent, n);
            parent[r as usize] = r0;
        }
    }
    let mut anchored = vec![false; node_idx.len()];
    for cl in clamps {
        let r = find(&mut parent, cell_nodes[cl.cell][0]);
        anchored[r as usize] = true;
    }
    let mut insulated = std::collections::HashSet::new();
    let (mut dropped_cells, mut dropped_volume) = (0, 0.0);
    let mut dropped_at = None;
    let multiple = if keep_free { parent.iter().enumerate().filter(|(i, p)| *i == **p as usize).nth(1).is_some() } else { anchored.iter().filter(|v| **v).nth(1).is_some() };
    let mut kept = std::collections::BTreeMap::<u32, KeptPiece>::new();
    let mut incomplete: Vec<Incomplete> = shared_hold.map(|(fixes, at)| Incomplete::SharedHold { fixes, at }).into_iter().collect();
    // (A study with a lip the cells cannot carry, or free slivers held with
    // a fix, where something acts on the lip's face, reads what acts there
    // wrong: its claims are not settled, not only told. A bar held to a
    // millionth short of a cell's plane and pulled on its end, its last
    // 0.004 mm beyond the plane, read 0.00067 N mm for 0.0064 with a band of
    // 0.0019, and its claim held. Where nothing acts on it, the study says
    // the lip is there: `idle`.)
    let idle: Vec<String> = match thin_lip {
        Some((lip, Some(_))) => {
            incomplete.extend(std::iter::once(lip).chain(absorbed_note).map(Incomplete::Held));
            Vec::new()
        }
        Some((lip, None)) => std::iter::once(lip).chain(absorbed_note).collect(),
        None => absorbed_note.into_iter().collect(),
    };
    if let Some((film_fi, ci_cls, q, fi)) = film_points.iter().find_map(|&(film_fi, ci_cls, q)| {
        let g = free_material.get(&ci_cls)?;
        (g.2[0] && g.1[0][0] == 0.0).then(|| holds.held_free_at(st, ci_cls, 0, q).map(|fi| (film_fi, ci_cls, q, fi))).flatten()
    }) {
        incomplete.push(Incomplete::Held(format!(
            "fix '{}' drops or absorbs free material under {}% of a cell near {:?} mm, and film '{}' acts on it: the cells cannot carry what it puts there, so the study is not settled. To settle it, write the hold's region out to the part's face, or refine (res r) to at least {:.4} mm",
            s.fixes[fi].name,
            100.0 * MIN_FRACTION,
            q.map(|v| (v * 1000.0).round() / 1000.0),
            s.fixes[film_fi].name,
            place(grid, &classified[ci_cls]).1 * 0.5
        )));
    }
    // A dropped sliver cell (`slivers_left`) whose kept face-neighbours are
    // all held has no trace piece for `thin_lip` to see: if a non-pressure
    // load or source acts on its free material along a held axis, the study
    // cannot settle.
    if !slivers_left.is_empty() && s.loads.iter().any(|l| l.pressure.is_none()) {
        'dropped_lip: for sl in slivers_left {
            if !sl.moments.as_ref().is_some_and(|m| m[0] > 1e-12 * 8.0 * sl.coef) {
                continue;
            }
            let neighbours: Vec<_> = face_neighbours(classified, sl, classified_at).collect();
            if neighbours.is_empty() || neighbours.iter().any(|ni| free_material.get(ni).is_none_or(|g| g.0[0] != 0.0)) {
                continue;
            }
            let (ctr, size) = place(grid, sl);
            let points: Vec<V3> = xis.iter().enumerate().filter(|&(k, _)| sl.mask >> k & 1 == 1).map(|(_, &xi)| at_xi(ctr, size, xi)).collect();
            for l in s.loads.iter().filter(|l| l.pressure.is_none()) {
                let weights = l.region.occupied(ctx, &points, size / K as f64);
                let Some(fi) = (0..ndof).find_map(|d| {
                    let pushes = l.total_t.is_some() || l.torque.is_some() || l.total[d] != 0.0;
                    if !pushes {
                        return None;
                    }
                    points.iter().zip(&weights).filter(|(_, w)| w.0 > 0.0).find_map(|(&p, _)| {
                        neighbours.iter().find_map(|&ni| {
                            let g = &free_material[&ni];
                            (g.2[d] && g.1[d][0] == 0.0).then(|| holds.held_free_at(st, ni, d, p)).flatten()
                        })
                    })
                }) else {
                    continue;
                };
                incomplete.push(Incomplete::Held(format!(
                        "fix '{}' drops free material under {}% of a cell near {:?} mm, and {} '{}' acts on it: the cells cannot carry what it puts there, so the study is not settled. To settle it, write the hold's region out to the part's face, or refine (res r)",
                        s.fixes[fi].name,
                        100.0 * MIN_FRACTION,
                        ctr.map(|v| (v * 1000.0).round() / 1000.0),
                        if ndof == 1 { "source" } else { "load" },
                        l.name
                    )));
                break 'dropped_lip;
            }
        }
    }
    if let Some((_, at, what)) = missed_at {
        incomplete.push(Incomplete::Separation(format!(
            "material the cells' samples missed lies between them near {at:?} mm, and {what}: a web, foil or ligament thinner than a quarter of a cell the mesh does not hold. Make it at least two cells thick, or refine (res r) to a quarter of its thickness"
        )));
    }
    let mut insulated_notes = Vec::new();
    let keep: Vec<bool> = classified
        .iter()
        .enumerate()
        .map(|(ci, c)| {
            let r = find(&mut parent, cell_nodes[ci][0]);
            let (ctr, side) = place(grid, c);
            if anchored[r as usize] || keep_free {
                let piece = kept.entry(r).or_default();
                if piece.cells.is_empty() || multiple && piece.cells.last() != Some(&c.idx) {
                    piece.cells.push(c.idx);
                }
                piece.volume += c.frac * side.powi(3);
                if !anchored[r as usize] && s.transient.is_some() {
                    insulated.insert(r);
                }
                return true;
            }
            dropped_at.get_or_insert(ctr);
            dropped_cells += 1;
            dropped_volume += c.frac * place(grid, c).1.powi(3);
            false
        })
        .collect();

    for (&root, piece) in &kept {
        let at = grid.pos(piece.cells[0]);
        if insulated.contains(&root) {
            insulated_notes.push(format!("{:.4} mm^3 in cells near {at:?} mm has no temperature fix or convection surface and was solved as insulated", piece.volume));
        }
        if incomplete.is_empty() {
            for (_, other) in kept.range(..root) {
                if !enclosed_apart(&s.shape, ctx, grid, piece, other) {
                    incomplete.push(Incomplete::Separation(format!("the mesh separated material near {at:?} mm from material near {:?} mm, but the shape may join them through a web the cells missed. Refine (res r) across the possible connection", grid.pos(other.cells[0]))));
                    break;
                }
            }
        }
    }
    Pieces { parent, insulated, keep, dropped_cells, dropped_volume, dropped_at, incomplete, duplicate_cut, idle, insulated_notes }
}
