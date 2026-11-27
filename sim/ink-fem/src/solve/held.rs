//! The held cells' terms (`held_cells`): what each hold's edge reads in
//! the cells beside it (`Lent`), the free slivers behind it, the lips it
//! leaves, and each cell's free stiffness with its holds' terms (`build`).

use crate::buckle::HeldGeometry;
use crate::cells::{ALPHA, K, LEND, MIN_FRACTION, locate, place};
use crate::element::{element_basis, full_moments, laminate_extended, laminate_matrix, shape_n, signs, traction, uncondensed};
use crate::eval::eval_points;
use crate::hold;
use crate::laminate::Laminate;
use crate::occupancy::{Classified, at_xi, spread_moments};
use crate::{Load, Physics, RegionSel};
use ink_core::{V3, eval, interval, v3};
use rayon::prelude::*;

use super::holds::{CellHolds, Clamp, HeldPart, Holds, clamp, material_of};
use super::settle::Cells;
use super::study::Study;
use super::{PENALTY, SolveError};

// Each cell holding held material, or the boundary of the free material
// with it: its free material's stiffness (the cell's own where it holds
// none), and its holds' terms (`hold::nitsche`) with that. A held
// region's material keeps its free axes' stiffness, and on its held
// ones, ALPHA's spring to the held value, as void keeps ALPHA of the
// material: no node is left without stiffness.

/// What the held cells' terms read: the cells, the tables of their
/// regions' material (`CellHolds`), and the free slivers behind pieces of
/// a hold's edge that are not traced (`slivers_behind`), set once those
/// pieces are known.
#[allow(clippy::type_complexity)]
struct Reads<'a> {
    st: &'a Study<'a>,
    cells: &'a Cells,
    tables: &'a std::collections::HashMap<usize, CellHolds>,
    slivers: std::sync::OnceLock<std::collections::HashMap<usize, Vec<(usize, [f64; 27])>>>,
}

impl Reads<'_> {
    /// Cell ci's material (moments, coefficient-weighted) in the
    /// combinations of its regions `keep` takes, given their fixes; where no
    /// region reaches it, all of it.
    /// (and a cell's free slivers behind pieces of a hold's edge that are
    /// not traced, each with its fix's material: `slivers`, set once the
    /// pieces are known)
    fn material_where(&self, ci: usize, keep: &dyn Fn(&[usize], u32) -> bool) -> [f64; 27] {
        let Reads { tables, ref slivers, .. } = *self;
        let classified = &self.cells.classified;
        let Some(t) = tables.get(&ci) else { return material_of(&classified[ci]) };
        let mut f = [0.0; 27];
        for (_, m) in t.pieces.iter().filter(|p| keep(&t.fixes, p.0)) {
            f.iter_mut().zip(m).for_each(|(f, v)| *f += v);
        }
        for (fi, m) in slivers.get().and_then(|s| s.get(&ci)).into_iter().flatten() {
            let held = t.fixes.iter().position(|f| f == fi).is_some_and(|b| keep(&t.fixes, 1 << b));
            let sign = f64::from(u8::from(held)) - f64::from(u8::from(keep(&t.fixes, 0)));
            f.iter_mut().zip(m).for_each(|(f, v)| *f += sign * v);
        }
        f
    }

    fn free_of(&self, ci: usize) -> [f64; 27] {
        self.material_where(ci, &|_, mask| mask == 0)
    }

    /// Who holds each axis of a combination of regions' material: a fix that
    /// moves it with free motions, where one holds that axis, else the first.
    fn owners(&self, t: &CellHolds, mask: u32) -> [Option<usize>; 3] {
        let Study { s, ndof, .. } = *self.st;
        let holding = |d: usize| (0..t.fixes.len()).filter(move |&b| d < ndof && mask >> b & 1 == 1 && s.fixes[t.fixes[b]].dofs[d]).map(|b| t.fixes[b]);
        std::array::from_fn(|d| holding(d).find(|&f| !s.fixes[f].free.is_empty()).or_else(|| holding(d).next()))
    }

    /// Material free of a set of axes: no region holding one of them holds it.
    /// A hold's edge sees that of the axes it holds there (a clamp's edge in
    /// a guide of y alone holds x and z, and sees the guide's material).
    fn free_on(&self, axes: [bool; 3]) -> impl Fn(&[usize], u32) -> bool + '_ {
        let Study { s, ndof, .. } = *self.st;
        move |fixes: &[usize], mask: u32| (0..fixes.len()).all(|b| mask >> b & 1 == 0 || !(0..ndof).any(|d| axes[d] && s.fixes[fixes[b]].dofs[d]))
    }

    /// The axes a hold's edges in a cell hold (its fix's where it has none).
    fn traced(&self, p: &HeldPart) -> [bool; 3] {
        let s = self.st.s;
        if p.traces.is_empty() {
            return s.fixes[p.fix].dofs;
        }
        let mut a = [false; 3];
        p.traces.iter().for_each(|t| (0..3).for_each(|d| a[d] |= t.p[d][d] != 0.0));
        a
    }

    fn phys_of(&self, c: &Classified) -> Physics {
        let s = self.st.s;
        match s.physics {
            Physics::Elastic { .. } => Physics::Elastic { nu: c.nu },
            p => p,
        }
    }

    /// A cell's free material's share of it, of the axes given (`free_on`).
    fn share(&self, k: usize, axes: [bool; 3]) -> f64 {
        let classified = &self.cells.classified;
        self.material_where(k, &self.free_on(axes))[0] / (8.0 * classified[k].coef)
    }

    /// The cell a piece's free normal leaves its own by (from a hair inside
    /// the cell across its other faces: a piece on the part's face found the
    /// air beyond).
    fn over(&self, ci: usize, t: &hold::Trace) -> Option<usize> {
        let (st, cells) = (self.st, self.cells);
        let (grid, classified) = (&st.grid, &cells.classified);
        let classified_at = |idx: [i64; 3]| cells.classified_at(st, idx);
        let (ctr, size) = place(grid, &classified[ci]);
        let q = at_xi(ctr, size, t.xi);
        let exit = (0..3).filter(|&a| t.n[a] != 0.0).map(|a| 0.5 * size * (1.0 + t.n[a].signum() * t.xi[a]) / t.n[a].abs()).fold(f64::INFINITY, f64::min);
        let probe: V3 = std::array::from_fn(|a| q[a] + 1e-9 * (ctr[a] - q[a]) - t.n[a] * (exit.max(0.0) + 1e-6 * size));
        let idx: [i64; 3] = std::array::from_fn(|a| ((probe[a] - grid.origin[a]) / grid.h).floor() as i64);
        locate(classified, classified_at(idx), grid, probe).map(|x| x.0)
    }
}

// A traced piece's traction is read partly in its own cell and partly in
// a cell beside it (`Lent`): all its own while the piece lies at the
// middle of its cell or deeper along its free normal, and the more the
// neighbour's (its field carried back to the piece) as the free material
// beyond the piece, or its cell's free share, goes to nothing; and its
// constant is the two cells' in the same shares. The piece's held value
// and test function stay its own cell's, so the hold is consistent; each
// cell's constant is twice what the tractions read in it (its own
// pieces', by their shares, and those lent to it) take of its free
// material's energy, so the system is positive definite; and as a
// piece's own share falls faster than its cell's free material's hold on
// its traction weakens, that constant stays bounded however little is
// free. A piece crossing a cell's plane is one piece on both sides of
// it. (Read in its own cell alone, a cell's constant grew as 1 / its
// free share: a corner block's reading jumped 2% as its face crossed a
// plane; handed whole to the neighbour, the test function the
// neighbour's too, a bar read 0.4% long.)
pub(super) struct Lent {
    pub(super) from: usize,
    pub(super) to: usize,
    pub(super) fix: usize,
    /// The piece, in its own cell.
    pub(super) t: hold::Trace,
    /// The neighbour's share of its boundary.
    pub(super) w: f64,
    /// The neighbour's traction on the held axes there, per its dof (`ndof` x `ne`).
    pub(super) pf: Vec<f64>,
}

/// The held cells' terms, and what the study is told of the holds.
pub(super) struct HeldCells {
    pub(super) lent: Vec<Lent>,
    /// A held cell's free material's stiffness, and its geometry (`HeldGeometry`).
    pub(super) free_matrix: std::collections::HashMap<usize, Vec<f64>>,
    pub(super) free_material: std::collections::HashMap<usize, HeldGeometry>,
    /// The holds' edges' area (mm^2).
    pub(super) held_edge: f64,
    /// Material two fixes hold on different axes, one moving it (`Incomplete::SharedHold`).
    pub(super) shared_hold: Option<([usize; 2], V3)>,
    /// A lip a hold leaves (`thin_lip`), and what acts on it.
    pub(super) thin_lip: Option<(String, Option<String>)>,
    /// Free slivers held with a fix (`absorbed`, `slivers_behind`), said.
    pub(super) absorbed_note: Option<String>,
}

pub(super) fn held_cells(st: &Study, cells: &Cells, holds: &mut Holds) -> Result<HeldCells, SolveError> {
    let Study { s, ref grid, .. } = *st;
    let classified = &cells.classified;
    let Holds { ref mut clamps, ref mut held_parts, ref tables, .. } = *holds;
    let reads = Reads { st, cells, tables, slivers: std::sync::OnceLock::new() };
    // Material two fixes hold on different axes, one moving it with free
    // motions: their stiffness between those axes is in each one's terms at
    // the other's given value alone (moves are), not where free motions
    // move it: two fixes of x and of y on one pivot read 0.0302 mm against
    // 0.0585 as one fix of both. Refused (`Incomplete::SharedHold`).
    let mut held_cells: Vec<usize> = tables.keys().copied().collect();
    held_cells.sort_unstable();
    let shared_hold: Option<([usize; 2], V3)> = held_cells.iter().find_map(|&ci| {
        let (c, t) = (&classified[ci], &tables[&ci]);
        t.pieces.iter().filter(|p| p.0 != 0 && p.1[0] > MIN_FRACTION * 8.0 * c.coef).find_map(|p| {
            let mut held: Vec<usize> = reads.owners(t, p.0).iter().flatten().copied().collect();
            held.sort_unstable();
            held.dedup();
            (held.len() > 1 && held.iter().any(|&f| !s.fixes[f].free.is_empty())).then(|| ([held[0], held[1]], place(grid, c).0.map(|v| (v * 1000.0).round() / 1000.0)))
        })
    });
    let held_edge: f64 = held_parts.iter().flat_map(|p| p.traces.iter().map(|t| t.w)).sum();
    let (lent, untraced) = lend(&reads, held_parts);
    // Free material under the floor of a cut cell's material (MIN_FRACTION:
    // a cell's material under it is dropped) that no cell beside lends a
    // field to carries none of its own: its piece of the hold's edge is not
    // traced (above), and a cell whose free material is all of it, no piece
    // left, is held whole with the region, and the study says how much.
    // (Kept, a bar held to 1e-6 short of a cell's plane, its last 0.004 mm
    // beyond the plane dropped as under the floor, hung the end's load on
    // the sliver: its nodes ran to 0.18 mm and its compliance read 26.28 N
    // mm for 0.0064. A lip of 0.001 mm between a held slab and the bar's end
    // face, its face traced and read in the bar's cell, held that cell's
    // field at both the slab's faces: 0.76% short.)
    let absorbed: std::collections::HashMap<usize, f64> = held_parts
        .chunk_by(|a, b| a.cell == b.cell)
        .filter(|on| on.iter().any(|p| p.vol > 0.0) && on.iter().all(|p| p.traces.is_empty()))
        .filter_map(|on| {
            let (ci, c) = (on[0].cell, &classified[on[0].cell]);
            let free = reads.free_of(ci)[0];
            (free > 1e-9 * material_of(c)[0] && free < MIN_FRACTION * 8.0 * c.coef).then(|| (ci, free / (8.0 * c.coef) * place(grid, c).1.powi(3)))
        })
        .collect();
    // In a cell that keeps other pieces, the sliver behind each piece not
    // traced is held with its fix: the material from the piece along its
    // free normal to the part's face or the cell's, taken from the cell's
    // free material. (Left free, a lip 0.001 mm thick between a held slab
    // and a bar's end face, in the cell of the slab's other face, was
    // strained as the bar's field carried on through it: the bar read
    // 0.005% short, its claim's band 0.0004% wide.)
    let (behind, behind_vol) = slivers_behind(&reads, &untraced, &absorbed);
    reads.slivers.set(behind).expect("set once");
    // A lip: free material between a piece of a hold's edge and the part's
    // face, thinner along the piece's free normal than a sample of its cell,
    // where the face there faces the piece more than it runs along it (a
    // hold's face that leaves the part at an angle bounds a wedge joined to
    // the rest beside it, no lip). The cell's one field serves the lip and
    // whatever else it holds, so what acts on the lip's face is read wrong
    // with grids that agree: a lip 0.01 mm thick pulled on its face moved 16
    // times too far, +-70%, on a bar turned 3.7 degrees; one 0.09 thick on
    // 0.5 mm cells read its work 0.33% under with a band of 0.22%. Where a
    // load or a film acts on the lip's face the study is not settled
    // (`acted`: its name); where none does, the lip idles, and the study
    // says it is there.
    let thin_lip = thin_lip(&reads, held_parts, &untraced, &absorbed);
    // (Said where more than a cell's volume's rounding: a rod held under its
    // own round face said it held 1.4e-16 mm^3.)
    let held_free: f64 = absorbed.values().chain(behind_vol.values()).sum();
    let absorbed_note = (held_free > 64.0 * f64::EPSILON * grid.h.powi(3)).then(|| {
        let ci = *absorbed.keys().chain(behind_vol.keys()).min().expect("a cell");
        let fix = held_parts.iter().find(|p| p.cell == ci).map_or(0, |p| p.fix);
        format!(
            "fix '{}' and the fixes beside it also hold {:.3e} mm^3 of free material, in {} cells near {:?} mm: free slivers under {}% of a cell with no free material beside them, which the cells cannot carry (as material that thin is dropped), so the study reads as if they were held. Write the hold's region out to the part's face, or refine (res r)",
            s.fixes[fix].name,
            held_free,
            absorbed.len() + behind_vol.len(),
            place(grid, &classified[ci]).0.map(|v| (v * 1000.0).round() / 1000.0),
            100.0 * MIN_FRACTION
        )
    });
    let mut lent_to: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for (i, l) in lent.iter().enumerate() {
        lent_to.entry(l.to).or_default().push(i);
    }
    // The nodes only held material reaches, on the axes it holds (every cell
    // round them all held there): no free material's field needs them, so
    // they are held hard, PENALTY's spring at each (lumped, so it ties no
    // other node). ALPHA's alone gave way to whatever else acts on them (a
    // pin held at 100 C all through, cooled by its air, convected 1% of its
    // heat at 100 C). A cell is all held on an axis where none of its
    // material is free on it: a wall's edge cell whose rest a guide of
    // another axis holds is free along the wall's (taken as held on every
    // axis either holds, its far face's nodes held a bar to that face).
    let inert = inert(&reads, held_parts, &absorbed);
    let mut with_holds: Vec<usize> = held_parts.iter().map(|p| p.cell).chain(lent.iter().map(|l| l.to)).collect();
    with_holds.sort_unstable();
    with_holds.dedup();
    let cell_terms = CellTerms { reads: &reads, held_parts, lent: &lent, lent_to: &lent_to, absorbed: &absorbed, inert: &inert };
    // Each cell's constant: twice its bound. One its free material cannot
    // bound at all (its eigenproblem singular) is no answer: refused by
    // name, not solved to NaN.
    #[allow(clippy::type_complexity)]
    let first: Vec<(usize, Option<(Vec<f64>, HeldGeometry)>, Vec<Clamp>, f64, bool)> = with_holds.par_iter().map(|&ci| cell_terms.build(ci, None)).collect();
    let constants: std::collections::HashMap<usize, f64> = first.iter().map(|x| (x.0, 2.0 * x.3)).collect();
    if let Some(&ci) = with_holds.iter().find(|ci| !constants[ci].is_finite()) {
        let fix = held_parts.iter().find(|p| p.cell == ci).map(|p| p.fix).or_else(|| lent_to.get(&ci).map(|l| lent[l[0]].fix)).unwrap_or(0);
        return Err(format!(
            "fix '{}' holds the part to an edge near {:?} mm in a cell whose free material cannot bound its hold's strength: give the hold a region whose edge lies in the part's free material, or refine (res r)",
            s.fixes[fix].name,
            place(grid, &classified[ci]).0.map(|v| (v * 1000.0).round() / 1000.0)
        )
        .into());
    }
    let built: Vec<(usize, Option<(Vec<f64>, HeldGeometry)>, Vec<Clamp>, f64, bool)> =
        first.into_par_iter().map(|x| if x.4 { x } else { cell_terms.build(x.0, Some(&constants)) }).collect();
    let mut free_matrix: std::collections::HashMap<usize, Vec<f64>> = std::collections::HashMap::new();
    let mut free_material: std::collections::HashMap<usize, HeldGeometry> = std::collections::HashMap::new();
    for (ci, kf, cls, _, _) in built {
        if let Some((kf, geometry)) = kf {
            free_matrix.insert(ci, kf);
            free_material.insert(ci, geometry);
        }
        clamps.extend(cls);
    }
    Ok(HeldCells { lent, free_matrix, free_material, held_edge, shared_hold, thin_lip, absorbed_note })
}

/// Each traced piece's lender (`Lent`), and the pieces left untraced: (cell, fix, piece).
fn lend(reads: &Reads, held_parts: &mut [HeldPart]) -> (Vec<Lent>, Vec<(usize, usize, hold::Trace)>) {
    let Reads { st, cells, tables, .. } = *reads;
    let Study { s, ndof, ne, ref grid, ref distributing, .. } = *st;
    let classified = &cells.classified;
    let classified_at = |idx: [i64; 3]| cells.classified_at(st, idx);
    let basis = element_basis(s.physics);
    let mut lent: Vec<Lent> = Vec::new();
    let mut untraced: Vec<(usize, usize, hold::Trace)> = Vec::new();
    for p in held_parts.iter_mut().filter(|p| !p.traces.is_empty()) {
        p.lenders = vec![usize::MAX; p.traces.len()];
        // (a distributing hold's pieces carry its weight alone, no traction:
        // nothing of theirs to read in another cell)
        if distributing[p.fix] {
            continue;
        }
        let axes = reads.traced(p);
        let (ctr, size) = place(grid, &classified[p.cell]);
        let own = reads.material_where(p.cell, &reads.free_on(axes))[0] / material_of(&classified[p.cell])[0];
        // The cells round it, of any size (a grid cell's children beside a
        // whole one), and the share of each free of every hold. (Free of this
        // hold's axes alone, a cell a guide on another axis holds all through
        // lent: its energy has none of what the lent traction draws of that
        // axis, its constant came out 850,000 of a cell's, and carried into
        // the pieces' own it took a solve 1,532 iterations for 36.)
        let home = classified[p.cell].idx;
        let free = |k: usize| reads.share(k, [true; 3]);
        let lenders: Vec<(usize, V3, f64, f64)> = (0..27)
            .flat_map(|o| classified_at(std::array::from_fn(|a| home[a] + [o % 3, (o / 3) % 3, o / 9][a] as i64 - 1)))
            .filter(|&k| k != p.cell && (classified[k].lam.is_none() || !tables.contains_key(&k)))
            .map(|k| (k, place(grid, &classified[k]).0, place(grid, &classified[k]).1, free(k)))
            .collect();
        let mut none: Vec<bool> = vec![false; p.traces.len()];
        for ((t, lender), none) in p.traces.iter_mut().zip(p.lenders.iter_mut()).zip(none.iter_mut()) {
            // Its own share: (twice its place along its free normal, of the
            // cell's reach along it; no more than twice the free share of
            // the cell's material, a corner's wedge being thinner than its
            // place says)^2. (Of the cell, not of its material, a cell the
            // part's face cut lent what it could carry: a cantilever's two
            // coarser grids moved 0.04%, and its claim's three grids no
            // longer read as converging.) Any power from 1 keeps the cell's
            // constant bounded as its free material goes to nothing (what a
            // traction takes of the energy behind it goes as 1 / its depth);
            // squared, the own share's part of the constant goes to nothing
            // too. Cubed reads within 0.003% of it; to the first power a
            // bar held 0.2 of a cell from a plane read 0.03% off its
            // 0.125 mm cells' reading for 0.007%.
            // (`depth`: 1 where the piece lies at its cell's middle. A whole
            // cell cut flat at the piece, a share f of it free, has the piece
            // at 2 f: so twice the free share of the cell's material is the
            // same measure by volume, and under twice the material floor is
            // free material under the floor.)
            let reach: f64 = t.n.iter().map(|v| v.abs()).sum();
            let depth = (1.0 + v3::dot(t.n, t.xi) / reach).min(2.0 * own).clamp(0.0, 1.0);
            let keep = depth.powi(2);
            if keep >= 1.0 {
                continue;
            }
            // The lender: of the cells at least `LEND` free, the one its free
            // normal leaves this one by, else the nearest (a lender's constant
            // goes as the inverse of its free share: at a half, about twice a
            // whole cell's); no further from the lender than a cell of the
            // lender's size beside it reaches (a small cell's field carried
            // further is no bound). With none, the cell its normal leaves by
            // where that holds more free material than its own: a lip of the
            // part beyond a hold's face, thinner than half a cell, is what
            // the piece bounds (with no lender, a face on the lip's cell's
            // plane was held by nothing, and a millionth short of it by its
            // own sliver at 330,000 times a cell's constant). None: the piece
            // stays its own.
            let q = at_xi(ctr, size, t.xi);
            let over = reads.over(p.cell, t);
            // (and on the piece's free side: the free material it bounds lies
            // there. A cell beside it along the hold's face lent its own:
            // the face of a lip 0.001 mm thick, between a held slab and a
            // bar's end, read the bar's side of the slab, 0.5% short.)
            // (3: in the lender's own coordinates, 1 is its face and 3 a
            // whole cell of its size beyond it, which is as far as the cell
            // beside a piece lies; 1e-9 a rounding)
            let near = |l: &&(usize, V3, f64, f64)| v3::dot(v3::sub(l.1, q), t.n) < 0.0 && (0..3).all(|a| ((q[a] - l.1[a]) * 2.0 / l.2).abs() <= 3.0 + 1e-9);
            let Some(&(k, kc, sk, _)) = lenders
                .iter()
                .filter(|l| l.3 >= LEND)
                .filter(near)
                .min_by(|a, b| (Some(b.0) == over).cmp(&(Some(a.0) == over)).then(v3::norm(v3::sub(a.1, q)).total_cmp(&v3::norm(v3::sub(b.1, q)))))
                .or_else(|| lenders.iter().filter(near).find(|l| Some(l.0) == over && l.3 > free(p.cell) && l.3 >= MIN_FRACTION))
            else {
                // (no cell lends it one, the free material behind it is
                // under the floor of a cut cell's material, half its depth
                // in cells, and no cell on its free side holds material free
                // of any axis it holds: it bounds nothing that carries a
                // field. Where one does, held on another axis by another fix
                // and so no lender, the piece is that material's hold and
                // stays its own: dropped, a hold of x, or a wall of all
                // three axes, a millionth short of a plane a guide of y
                // began at held x nowhere.)
                let free_beyond = (0..ndof).any(|d| {
                    let mut only = [false; 3];
                    only[d] = true;
                    t.p[d][d] != 0.0 && lenders.iter().filter(near).any(|l| reads.share(l.0, only) >= MIN_FRACTION)
                });
                *none = depth < 2.0 * MIN_FRACTION && !free_beyond;
                continue;
            };
            let b = &classified[k];
            let f = traction(basis, reads.phys_of(b), &reads.material_where(k, &reads.free_on(axes)), b.coef, sk, b.lam.as_deref())(
                std::array::from_fn(|a| (q[a] - kc[a]) * 2.0 / sk),
                t.n,
            );
            let pf = (0..ndof * ne).map(|r| (0..ndof).map(|d| t.p[r / ne][d] * f[d * ne + r % ne]).sum()).collect();
            lent.push(Lent { from: p.cell, to: k, fix: p.fix, t: *t, w: t.w * (1.0 - keep), pf });
            t.wt = t.w * keep;
            *lender = k;
        }
        untraced.extend(p.traces.iter().zip(&none).filter(|x| *x.1).map(|x| (p.cell, p.fix, *x.0)));
        let mut at = none.iter();
        p.traces.retain(|_| !*at.next().expect("a trace's"));
        let mut at = none.iter();
        p.lenders.retain(|_| !*at.next().expect("a trace's"));
    }
    (lent, untraced)
}

#[allow(clippy::type_complexity)]
fn slivers_behind(
    reads: &Reads,
    untraced: &[(usize, usize, hold::Trace)],
    absorbed: &std::collections::HashMap<usize, f64>,
) -> (std::collections::HashMap<usize, Vec<(usize, [f64; 27])>>, std::collections::HashMap<usize, f64>) {
    let Study { s, ctx, ref grid, .. } = *reads.st;
    let classified = &reads.cells.classified;
    let mut behind: std::collections::HashMap<usize, Vec<(usize, [f64; 27])>> = std::collections::HashMap::new();
    let mut behind_vol: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    for (ci, fi, t) in untraced.iter().filter(|x| !absorbed.contains_key(&x.0)) {
        let c = &classified[*ci];
        let (ctr, size) = place(grid, c);
        let q = at_xi(ctr, size, t.xi);
        let exit = (0..3).filter(|&a| t.n[a] != 0.0).map(|a| 0.5 * size * (1.0 + t.n[a].signum() * t.xi[a]) / t.n[a].abs()).fold(f64::INFINITY, f64::min).max(0.0);
        let inside = |d: f64| eval::eval_point(&s.shape, ctx, std::array::from_fn(|a| q[a] - t.n[a] * d)) <= 0.0;
        let (mut lo, mut hi) = (0.0, exit);
        if !inside(exit) {
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                if inside(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }
        let (tau, to_xi) = (hi, 2.0 / size);
        if !(tau > 0.0 && tau.is_finite()) {
            continue;
        }
        let at: V3 = std::array::from_fn(|a| t.xi[a] - 0.5 * tau * to_xi * t.n[a]);
        let k = (tau * to_xi).powi(2) / 12.0;
        let cov = [k * t.n[0] * t.n[0], k * t.n[1] * t.n[1], k * t.n[2] * t.n[2], k * t.n[0] * t.n[1], k * t.n[1] * t.n[2], k * t.n[0] * t.n[2]];
        let m = spread_moments(at, &cov).map(|v| v * t.w * tau * to_xi.powi(3) * c.coef);
        behind.entry(*ci).or_default().push((*fi, m));
        *behind_vol.entry(*ci).or_insert(0.0) += t.w * tau;
    }
    // (no more than the cell's free material: a sliver wider than its piece
    // says, its face tilted, is all of it)
    for (ci, v) in behind.iter_mut() {
        let (free, all) = (reads.free_of(*ci)[0], v.iter().map(|x| x.1[0]).sum::<f64>());
        if all > free {
            v.iter_mut().for_each(|x| x.1.iter_mut().for_each(|m| *m *= free / all));
        }
    }
    (behind, behind_vol)
}

fn thin_lip(
    reads: &Reads,
    held_parts: &[HeldPart],
    untraced: &[(usize, usize, hold::Trace)],
    absorbed: &std::collections::HashMap<usize, f64>,
) -> Option<(String, Option<String>)> {
    let st = reads.st;
    let Study { s, ctx, ndof, h, ref grid, ref distributing, .. } = *st;
    let classified = &reads.cells.classified;
    {
        let pieces: Vec<(usize, usize, hold::Trace)> =
            held_parts.iter().filter(|p| !distributing[p.fix]).flat_map(|p| p.traces.iter().map(move |t| (p.cell, p.fix, *t))).chain(untraced.iter().copied()).collect();
        // (those after are not traced: their lips are held with their fix)
        let traced_pieces = pieces.len() - untraced.len();
        let sample = |ci: usize| place(grid, &classified[ci]).1 / K as f64;
        let from = |x: &(usize, usize, hold::Trace), d: f64| -> V3 {
            let (ctr, size) = place(grid, &classified[x.0]);
            let q = at_xi(ctr, size, x.2.xi);
            std::array::from_fn(|a| q[a] - x.2.n[a] * d)
        };
        // (Where the part ends within a cell of the piece: the face's place,
        // by halving, to a point of the part there. Within a sample for what
        // acts on the lip's face; within a cell for what acts through it.)
        let mut beyond = Vec::new();
        eval_points(&s.shape, ctx, &pieces.iter().map(|x| from(x, K as f64 * sample(x.0))).collect::<Vec<V3>>(), &mut beyond);
        let ends: Vec<(usize, f64)> = (0..pieces.len())
            .filter(|&i| beyond[i] > 0.0)
            .map(|i| {
                let (mut lo, mut hi) = (0.0, K as f64 * sample(pieces[i].0));
                for _ in 0..40 {
                    let mid = 0.5 * (lo + hi);
                    if eval::eval_point(&s.shape, ctx, from(&pieces[i], mid)) <= 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                (i, lo)
            })
            // (a piece on the part's face itself bounds no material: a hold
            // written to end at a face, to a rounding)
            .filter(|&(i, d)| d > 1e-9 * sample(pieces[i].0))
            .collect();
        let faces = crate::eval::eval_grads(&s.shape, ctx, &ends.iter().map(|&(i, d)| from(&pieces[i], d)).collect::<Vec<V3>>());
        let facing: Vec<(usize, f64)> = ends
            .iter()
            .zip(&faces)
            .filter(|((i, _), g)| v3::direction([g[1], g[2], g[3]]).is_some_and(|out| -v3::dot(out, pieces[*i].2.n) > std::f64::consts::FRAC_1_SQRT_2))
            .map(|(e, _)| *e)
            .collect();
        // A lip is of one thickness: half of it or more, by its edge's
        // area, within a tenth of its thickest. Where a rounded member
        // leaves a held face (a fillet's root, tangent to the face) or a
        // face leaves it at a shallow angle, the free material between
        // runs from nothing up to the thickness asked about (a twentieth
        // of a fillet's root within a tenth of its thickest, a tenth of a
        // flat wedge): joined to the thick of the part, which the cells
        // carry. (Taken for lips, every clamp with a filleted member
        // leaving it was told of one: a flexure's, 0.024 mm at its roots.)
        let of_one_thickness = |within: &dyn Fn(usize, f64) -> bool| -> Vec<(usize, f64)> {
            let mut by_fix: std::collections::BTreeMap<usize, Vec<(usize, f64)>> = std::collections::BTreeMap::new();
            facing.iter().filter(|x| within(x.0, x.1)).for_each(|&(i, d)| by_fix.entry(pieces[i].1).or_default().push((i, d)));
            by_fix
                .into_values()
                .filter_map(|mut of| {
                    of.sort_by(|a, b| a.1.total_cmp(&b.1));
                    let (all, thickest) = (of.iter().map(|x| pieces[x.0].2.w).sum::<f64>(), of.last()?.1);
                    let mut below = 0.0;
                    let middle = of.iter().find(|x| {
                        below += pieces[x.0].2.w;
                        below >= 0.5 * all
                    })?;
                    (middle.1 > 0.9 * thickest).then_some(of)
                })
                .flatten()
                .collect()
        };
        // Thinner than a sample: it may not carry what a pressure or a film
        // puts on its face (and it is told of where it idles). Thinner than
        // a cell, and traced: no node lies in it, so a force or a source
        // spread through its thickness, whose field bends across it, it
        // cannot carry either, though a pressure on its face, whose field
        // is straight across it, it can (a lip 0.3 mm thick on 0.5 mm cells:
        // 0.12 +- 1.5e-6 for 0.12 pressed on its face; 0.0452 +- 0.0023 for
        // 0.0533 with 40 N over a band inside it). What a load spread over a
        // region puts on a lip held with its fix, material under the floor,
        // is its share of the load by volume, nothing: there only a
        // pressure or a film acts (a held temperature written a millionth
        // inside a part's sides, the part heated all through, was left
        // unsettled over lips 1e-6 mm thick).
        let lips = of_one_thickness(&|i, d| d < sample(pieces[i].0));
        let wide: Vec<(usize, f64)> = of_one_thickness(&|i, _| i < traced_pieces).into_iter().chain(lips.iter().copied().filter(|x| x.0 < traced_pieces)).collect();
        // A traced, grid-aligned lip between parallel faces in one cell
        // can carry a surface pressure with a field linear across it.
        // Require the free material to survive the floor in a cell cut
        // by the hold. A trace on a cell plane says nothing about how
        // much material was dropped from the cells beside that plane.
        let carries_pressure = |(i, d): (usize, f64)| {
            let (ci, _, trace) = pieces[i];
            let c = &classified[ci];
            if i >= traced_pieces
                || trace.wt != trace.w
                || !trace.n.iter().zip(trace.xi).any(|(n, xi)| n.abs() >= 1.0 - 1e-12 && xi.abs() < 1.0 - 1e-9)
                || absorbed.contains_key(&ci)
                || reads.free_of(ci)[0] < MIN_FRACTION * 8.0 * c.coef
            {
                return false;
            }
            let (ctr, size) = place(grid, c);
            let face = from(&pieces[i], d);
            let g = crate::eval::eval_grads(&s.shape, ctx, &[face])[0];
            d >= MIN_FRACTION * size
                && !st.stiffness_varies(ctr, size)
                && v3::direction([g[1], g[2], g[3]]).is_some_and(|out| -v3::dot(out, trace.n) >= 1.0 - 1e-12)
                && (0..3).all(|a| (face[a] - ctr[a]).abs() <= (0.5 + 1e-9) * size)
        };
        let face_of = |set: &[(usize, f64)]| -> Vec<V3> { set.iter().map(|&(j, d)| from(&pieces[j], d)).collect() };
        let holds = |r: &RegionSel, at: &[V3]| r.occupied(ctx, at, 1e-6 * h).iter().any(|o| o.0 > 0.0);
        // (A region's field over the line from a piece to the lip's face;
        // and the load on an axis the lip's fix holds: what a force along a
        // plate puts on it does not pass through its thickness to a roller
        // under it. A plate one cell thick on a symmetry layer, pulled in
        // its plane, was left unsettled as a lip.)
        let through = |l: &Load| -> Option<(usize, f64)> {
            let pushes = |fi: usize| l.total_t.is_some() || l.torque.is_some() || (0..ndof).any(|d| l.total[d] != 0.0 && s.fixes[fi].dofs[d]);
            let r = &l.region;
            wide.iter().copied().filter(|x| pushes(pieces[x.0].1)).find(|&(j, d)| match r {
                RegionSel::Field(_) => {
                    let (a, b) = (from(&pieces[j], 0.0), from(&pieces[j], d));
                    r.touches_plane(ctx, std::array::from_fn(|k| interval::I::new(a[k].min(b[k]), a[k].max(b[k]))))
                }
                RegionSel::Face { .. } => holds(r, &[from(&pieces[j], d)]),
            })
        };
        let name = |l: &Load| format!("{} '{}'", if ndof == 1 { "source" } else { "load" }, l.name);
        let write = |ci: usize| {
            format!(
                "write the hold's region out to the part's face, or leave at least two cells ({:.4} mm) of free material between its edge and the face, or refine (res r) until the lip is two cells thick",
                2.0 * place(grid, &classified[ci]).1
            )
        };
        let told = |set: &[(usize, f64)], i: usize, thin: f64, than: &str, size: f64| {
            format!(
                "fix '{}' leaves free material {:.3e} mm thick between its edge and the part's face near {:?} mm (over {:.3e} mm^2 of its edge), thinner than a {than} of the cells there ({:.4} mm)",
                s.fixes[pieces[i].1].name,
                thin,
                from(&pieces[i], 0.0).map(|v| (v * 1000.0).round() / 1000.0),
                set.iter().map(|&(j, _)| pieces[j].2.w).sum::<f64>(),
                size
            )
        };
        // What acts on a thin lip's face: a pressure's region or a film's,
        // holding a point of it. (Another fix holding that face makes it no
        // free face: the two fixes are a sliver apart, as a held side
        // written to stop 2e-6 mm under a held top is, and the study says
        // so.)
        let on_face = face_of(&lips);
        let fixes: std::collections::HashSet<usize> = lips.iter().map(|&(j, _)| pieces[j].1).collect();
        let other =
            |film: bool| s.fixes.iter().enumerate().find(|(j, f)| !fixes.contains(j) && f.film.is_some() == film && holds(&f.region, &on_face)).map(|(_, f)| f.name.clone());
        let pressure_lips: Vec<_> = lips.iter().copied().filter(|&lip| !carries_pressure(lip)).collect();
        let pressure_faces = face_of(&pressure_lips);
        let pressed = s
            .loads
            .iter()
            .find(|l| l.pressure.is_some() && holds(&l.region, &pressure_faces))
            .map(|l| (&pressure_lips, name(l)))
            .or_else(|| other(true).map(|f| (&lips, format!("film '{f}'"))));
        let spread = s.loads.iter().filter(|l| l.pressure.is_none()).find_map(|l| through(l).map(|at| (name(l), at)));
        let thinnest = pressure_lips.iter().copied().min_by(|a, b| a.1.total_cmp(&b.1));
        match (thinnest, pressed, spread) {
            (_, Some((set, by)), _) => {
                let (i, thin) = set.iter().copied().min_by(|a, b| a.1.total_cmp(&b.1)).expect("an acted-on lip");
                let ci = pieces[i].0;
                Some((
                    format!(
                        "{}, and {by} acts on it: the cells cannot carry what it puts there, so the study is not settled. To settle it, {}",
                        told(set, i, thin, "sample", sample(ci)),
                        write(ci)
                    ),
                    Some(by),
                ))
            }
            (_, _, Some((by, (i, thick)))) => {
                let ci = pieces[i].0;
                Some((
                    format!(
                        "{}, and {by} acts within it: no node of the cells lies in it to carry what is spread through it, so the study is not settled. To settle it, {}",
                        told(&wide, i, thick, "cell", K as f64 * sample(ci)),
                        write(ci)
                    ),
                    Some(by),
                ))
            }
            (Some((i, thin)), None, None) => {
                let ci = pieces[i].0;
                let lip = told(&pressure_lips, i, thin, "sample", sample(ci));
                Some(match other(false) {
                    Some(fix) => (
                        format!(
                            "{lip}, which fix '{fix}' holds: the two fixes are that far apart, and what passes between them through it is not resolved. To be rid of it, {}",
                            write(ci)
                        ),
                        None,
                    ),
                    None => (format!("{lip}: no load or film acts on it, so it idles, read as the cells can. To be rid of it, {}", write(ci)), None),
                })
            }
            (None, _, None) => None,
        }
    }
}

fn inert(reads: &Reads, held_parts: &[HeldPart], absorbed: &std::collections::HashMap<usize, f64>) -> Vec<[bool; 3]> {
    let Study { s, ndof, .. } = *reads.st;
    let Cells { ref classified, ref cell_nodes, ref node_idx, .. } = *reads.cells;
    {
        let mut wholly: std::collections::HashMap<usize, [bool; 3]> = std::collections::HashMap::new();
        for on in held_parts.chunk_by(|a, b| a.cell == b.cell) {
            let ci = on[0].cell;
            let whole = material_of(&classified[ci])[0];
            let absorbed = absorbed.contains_key(&ci);
            if absorbed || reads.free_of(ci)[0] <= 1e-9 * whole {
                let mut axes = [false; 3];
                for p in on.iter().filter(|p| p.vol > 0.0) {
                    (0..ndof).for_each(|d| axes[d] |= s.fixes[p.fix].dofs[d]);
                }
                if !absorbed {
                    for (d, a) in axes.iter_mut().enumerate().take(ndof) {
                        let mut only = [false; 3];
                        only[d] = true;
                        *a &= reads.material_where(ci, &reads.free_on(only))[0] <= 1e-9 * whole;
                    }
                }
                wholly.insert(ci, axes);
            }
        }
        let mut inert = vec![[true; 3]; node_idx.len()];
        for (ci, nodes) in cell_nodes.iter().enumerate() {
            let axes = wholly.get(&ci).copied().unwrap_or([false; 3]);
            for &n in nodes {
                (0..3).for_each(|d| inert[n as usize][d] &= axes[d]);
            }
        }
        inert
    }
}

/// What a held cell's terms read (`CellTerms::build`).
struct CellTerms<'a> {
    reads: &'a Reads<'a>,
    held_parts: &'a [HeldPart],
    lent: &'a [Lent],
    lent_to: &'a std::collections::HashMap<usize, Vec<usize>>,
    absorbed: &'a std::collections::HashMap<usize, f64>,
    inert: &'a [[bool; 3]],
}

impl CellTerms<'_> {
    /// (A cell's terms, its pieces at the cells' constants; with none given,
    /// for its bound alone.)
    #[allow(clippy::type_complexity)]
    fn build(&self, ci: usize, constants: Option<&std::collections::HashMap<usize, f64>>) -> (usize, Option<(Vec<f64>, HeldGeometry)>, Vec<Clamp>, f64, bool) {
        let CellTerms { reads, held_parts, lent, lent_to, absorbed, inert } = *self;
        let Reads { st, cells, tables, ref slivers } = *reads;
        let Study { s, ndof, ne, ref grid, ref distributing, .. } = *st;
        let Cells { ref classified, ref cell_nodes, .. } = *cells;
        let basis = element_basis(s.physics);
        let fm = full_moments();
        let on = &held_parts[held_parts.partition_point(|p| p.cell < ci)..held_parts.partition_point(|p| p.cell <= ci)];
        let (c, table) = (&classified[ci], tables.get(&ci));
        let (ctr, size) = place(grid, c);
        let half = 0.5 * size;
        let phys = match s.physics {
            Physics::Elastic { .. } => Physics::Elastic { nu: c.nu },
            p => p,
        };
        // The cell's stiffness for the material of moments m
        // (coefficient-weighted): its laminate's where it is one
        // (`laminate`, of that material's share of each: `lam_where`),
        // else isotropic; and its whole cell's, times its coefficient.
        // (A sliver of material, under MIN_FRACTION of a cell, without
        // its incompatible modes: condensed on a piece that thin, they
        // ran to 1e15, and the solve broke down.)
        let stiff = |m: &[f64; 27], lam: Option<&Laminate>| -> Vec<f64> {
            let sliver = m[0] < MIN_FRACTION * 8.0 * c.coef;
            match (lam, phys) {
                (Some(lam), Physics::Elastic { .. }) if sliver => uncondensed(&laminate_extended(lam, &m.map(|v| v / c.coef)).expect("an elastic basis").0, half),
                (None, Physics::Elastic { nu }) if sliver => uncondensed(&basis.extended(nu, m).expect("an elastic basis").0, half),
                (Some(lam), _) => laminate_matrix(lam, s.physics, &m.map(|v| v / c.coef), half),
                (None, _) => basis.matrix(phys, m, half),
            }
        };
        let lam_where = |keep: &dyn Fn(&[usize], u32) -> bool| -> Option<Laminate> {
            let l = c.lam.as_deref()?;
            let (v, a) = table.map_or((0.0, 0.0), |t| t.shares.iter().filter(|x| keep(&t.fixes, x.0)).fold((0.0, 0.0), |(v, a), x| (v + x.1, a + x.2)));
            Some(Laminate { fa: if v > 0.0 { a / v } else { l.fa }, ..l.clone() })
        };
        let full: Vec<f64> = match c.lam.as_deref() {
            Some(lam) => laminate_matrix(lam, s.physics, &fm, half),
            None => basis.matrix(phys, &fm, half).iter().map(|v| v * c.coef).collect(),
        };
        let mat = material_of(c);
        let is_absorbed = absorbed.contains_key(&ci);
        let free_where = |keep: &dyn Fn(&[usize], u32) -> bool| if is_absorbed { [0.0; 27] } else { reads.material_where(ci, keep) };
        let (free, lam_free) = (free_where(&|_, mask| mask == 0), lam_where(&|_, mask| mask == 0));
        let holds = on.iter().any(|p| p.vol > 0.0);
        let any_free = free[0] > 1e-9 * mat[0];
        let kf: Vec<f64> = match (holds, &c.moments) {
            (false, None) => full.clone(),
            (false, Some(m)) => stiff(m, c.lam.as_deref()).iter().zip(&full).map(|(k, f)| k + ALPHA * f).collect(),
            (true, _) if any_free => stiff(&free, lam_free.as_ref()).iter().zip(&full).map(|(k, f)| k + ALPHA * f).collect(),
            (true, _) => full.iter().map(|f| ALPHA * f).collect(),
        };
        let corners: Vec<V3> = (0..8).map(|i| at_xi(ctr, size, signs(i))).collect();
        let mut terms: Vec<hold::Terms> = on.iter().map(|_| hold::Terms::zero(ndof)).collect();
        // The cell's stiffness on its free axes: its free material's, and
        // held material's on the axes left free where what holds it does
        // not move (moving with a pivot, it bears less: its block's
        // complement over the held axes; counted whole, a beam on a
        // roller's pivot did not converge).
        let mut k_act = kf.clone();
        // A held region's free axes stay material: their stiffness, and
        // its coupling to the held ones. Material in a combination of
        // regions is held on the axes one of them holds, by the first of
        // them that does on each, and free on the rest (including slivers
        // under MIN_FRACTION, which `stiff` builds without incompatible
        // modes, and any sliver moved into a fix by `behind`).
        if let Some(tb) = table {
            let index = |fi: usize| on.iter().position(|p| p.fix == fi);
            let mut at: std::collections::HashMap<usize, Vec<f64>> = std::collections::HashMap::new();
            let mut held_pieces: Vec<(u32, [f64; 27])> = tb.pieces.iter().filter(|p| p.0 != 0).copied().collect();
            if !is_absorbed {
                for &(fi, m) in slivers.get().and_then(|s| s.get(&ci)).into_iter().flatten() {
                    if let Some(b) = tb.fixes.iter().position(|&f| f == fi) {
                        let mask = 1u32 << b;
                        match held_pieces.iter_mut().find(|p| p.0 == mask) {
                            Some(p) => p.1.iter_mut().zip(m).for_each(|(a, b)| *a += b),
                            None => held_pieces.push((mask, m)),
                        }
                    }
                }
            }
            for (mask, m) in held_pieces.iter().filter(|p| p.1[0] > 1e-12 * 8.0 * c.coef) {
                let owner = reads.owners(tb, *mask);
                if (0..ndof).all(|d| owner[d].is_some()) {
                    continue;
                }
                for &fi in owner.iter().flatten() {
                    at.entry(fi).or_insert_with(|| st.values_at(&s.fixes[fi], &corners).iter().flat_map(|v| v[..ndof].to_vec()).collect());
                }
                let first = index(tb.fixes[mask.trailing_zeros() as usize]);
                let kb = stiff(m, lam_where(&|_, b| b == *mask).as_ref());
                let still = (0..tb.fixes.len()).all(|b| mask >> b & 1 == 0 || s.fixes[tb.fixes[b]].free.is_empty());
                for i in 0..ne {
                    for j in 0..ne {
                        let v = kb[i * ne + j];
                        match (owner[i % ndof], owner[j % ndof]) {
                            (None, None) => {
                                if let Some(k) = first {
                                    terms[k].kh[i * ne + j] += v;
                                }
                                if still {
                                    k_act[i * ne + j] += v;
                                }
                            }
                            (None, Some(o)) => {
                                if let Some(k) = index(o) {
                                    terms[k].kuw_mut()[i * ne + j] += v;
                                    terms[k].gu[i] -= v * at[&o][j];
                                }
                            }
                            (Some(o), Some(o2)) => {
                                if let Some(k) = index(o) {
                                    if o == o2 {
                                        terms[k].kww_mut()[i * ne + j] += v;
                                    }
                                    terms[k].gw[i] += v * at[&o2][j];
                                }
                            }
                            (Some(_), None) => {}
                        }
                    }
                }
            }
        }
        // Each hold's boundary is the free material's on its axes: its
        // traction is the material's no region holding one of the axes it
        // holds there holds (`free_on`)
        // (a guide holding a part along one axis all through left a
        // support beside it, along another, no free material: its
        // constant was the void's, and a beam it held buckled 13% over
        // Euler).
        let fluxes: Vec<Option<Box<dyn Fn(V3, V3) -> Vec<f64> + Sync + '_>>> = on
            .iter()
            .map(|p| {
                let keep = reads.free_on(reads.traced(p));
                let only_free = table.is_none_or(|t| t.pieces.iter().all(|q| keep(&t.fixes, q.0) == (q.0 == 0)));
                let (gf, lam) = if only_free { (free, lam_free.clone()) } else { (free_where(&keep), lam_where(&keep)) };
                (!p.traces.is_empty()).then(|| traction(basis, phys, &gf, c.coef, size, lam.as_ref()))
            })
            .collect();
        let traces: Vec<hold::Trace> = on.iter().enumerate().flat_map(|(k, p)| p.traces.iter().map(move |t| hold::Trace { hold: k, ..*t })).collect();
        // Each piece at its constant: its cell's and its lender's, in its
        // shares (`Lent`). The cell's bound: the largest its own pieces'
        // tractions and those lent to it take of its free axes' energy.
        // (A cell none of whose pieces lends takes its own constant alone,
        // which its own bound gives: its terms are made once, at 1, and set
        // to it (`Terms::at`). One that lends waits for the others'
        // constants and is made again with them: the tractions, which are
        // the cost, were read twice in every held cell, and a bracket's
        // set-up took 0.36 s for 0.20.)
        let alone = constants.is_none() && on.iter().all(|p| p.lenders.iter().all(|&l| l == usize::MAX));
        let constant = |k: usize| constants.and_then(|g| g.get(&k)).copied().unwrap_or(0.0);
        let gammas: Vec<f64> = on
            .iter()
            .flat_map(|p| {
                p.traces.iter().zip(&p.lenders).map(move |(t, &l)| {
                    let keep = if l == usize::MAX { 1.0 } else { t.wt / t.w };
                    if alone { 1.0 } else { keep * constant(ci) + (1.0 - keep) * constant(l) }
                })
            })
            .collect();
        let mut gram = vec![0.0; ne * ne];
        let mut held: Vec<Option<hold::Terms>> = Vec::new();
        if !traces.is_empty() {
            let flux = |k: usize, xi: V3, n: V3| fluxes[k].as_ref().map_or_else(|| vec![0.0; ndof * ne], |f| f(xi, n));
            (held, gram) = hold::nitsche(ndof, &traces, &flux, on.len(), &gammas);
        }
        for l in lent_to.get(&ci).into_iter().flatten().map(|&i| &lent[i]) {
            for ij in 0..ne * ne {
                gram[ij] += l.w * (0..ndof).map(|a| l.pf[a * ne + ij / ne] * l.pf[a * ne + ij % ne]).sum::<f64>();
            }
        }
        // (the bound is the first pass's)
        let bound = if constants.is_none() && gram.iter().any(|v| *v != 0.0) { hold::largest(ndof, &k_act, &gram) } else { 0.0 };
        if constants.is_none() && !alone {
            return (ci, None, Vec::new(), bound, false);
        }
        for (t, n) in terms.iter_mut().zip(held) {
            if let Some(n) = n {
                t.add(&if alone { n.at(2.0 * bound) } else { n });
            }
        }
        for (p, t) in on.iter().zip(terms.iter_mut()) {
            let fx = &s.fixes[p.fix];
            t.add(&hold::penalty(ndof, &fx.dofs, if distributing[p.fix] { &[0.0; 64] } else { &p.m }, &p.g, ALPHA * c.coef / (size * size)));
            // Its nodes only held material reaches (`inert`), held hard.
            if !any_free && p.vol > 0.0 {
                let held_at: Vec<usize> = (0..8).filter(|&k| (0..ndof).any(|d| fx.dofs[d] && inert[cell_nodes[ci][k] as usize][d])).collect();
                let values = st.values_at(fx, &held_at.iter().map(|&k| corners[k]).collect::<Vec<V3>>());
                for (&k, v) in held_at.iter().zip(&values) {
                    let (mut m, mut g) = ([0.0; 64], vec![0.0; ne]);
                    m[k * 8 + k] = 1.0;
                    let axes: [bool; 3] = std::array::from_fn(|d| d < ndof && fx.dofs[d] && inert[cell_nodes[ci][k] as usize][d]);
                    (0..ndof).filter(|&d| axes[d]).for_each(|d| g[k * ndof + d] = v[d]);
                    t.add(&hold::penalty(ndof, &axes, &m, &g, PENALTY * c.coef * size / 8.0));
                }
            }
        }
        let spring = ALPHA * c.coef / (size * size);
        let mut clamps: Vec<Clamp> =
            on.iter().zip(terms).map(|(p, t)| Clamp { vol: p.vol, hn: p.hn, whole: p.whole || is_absorbed, spring: Some((spring, p.m)), ..clamp(ci, p.fix, t) }).collect();
        // Material two fixes hold, where one leaves motions free (a
        // pivot, a slide) and the other none: the other holds that one's
        // free motions there, by a penalty on the held value alone, over
        // what they share: a hinge's slab, its bottom on a floor that
        // holds it, turned freely (its pull's moment took a quarter of
        // the floor's reaction). The free material is not touched: held
        // over the shared material, its field there (the free field's
        // extension in a cut cell) clamped a plate's simply supported
        // edges where guides crossed them (buckling 2.92 times the load
        // for 2.26). Two fixes that both leave motions free hold the
        // shared material each, by its penalty, as one solve of both
        // their motions would.
        if let Some(tb) = table {
            for (a, pa) in on.iter().enumerate() {
                for pb in on[a + 1..].iter().filter(|pb| pb.fix != pa.fix) {
                    let (fa, fb) = (&s.fixes[pa.fix], &s.fixes[pb.fix]);
                    let (Some(ba), Some(bb)) = (tb.fixes.iter().position(|&f| f == pa.fix), tb.fixes.iter().position(|&f| f == pb.fix)) else { continue };
                    if fa.free.is_empty() && fb.free.is_empty() {
                        continue;
                    }
                    let shared: Vec<&(u32, V3, f64, V3)> = tb.held.iter().filter(|x| x.0 >> ba & 1 == 1 && x.0 >> bb & 1 == 1).collect();
                    let vol: f64 = shared.iter().map(|x| x.2).sum();
                    if vol <= 0.0 {
                        continue;
                    }
                    let mut m = [0.0; 64];
                    for x in &shared {
                        (0..64).for_each(|ij| m[ij] += x.2 * shape_n(ij / 8, x.1) * shape_n(ij % 8, x.1));
                    }
                    let beta = PENALTY * c.coef * size / vol;
                    let at: Vec<V3> = shared.iter().map(|x| x.3).collect();
                    if !fa.free.is_empty() && !fb.free.is_empty() {
                        for p in [pa, pb] {
                            let fx = &s.fixes[p.fix];
                            let mut g = vec![0.0; ne];
                            for (x, v) in shared.iter().zip(st.values_at(fx, &at)) {
                                for i in 0..8 {
                                    (0..ndof).filter(|&d| fx.dofs[d]).for_each(|d| g[i * ndof + d] += x.2 * shape_n(i, x.1) * v[d]);
                                }
                            }
                            clamps.push(clamp(ci, p.fix, hold::penalty(ndof, &fx.dofs, &m, &g, beta)));
                        }
                        continue;
                    }
                    // `moving` leaves motions free; `by` holds them on the
                    // axes both hold: beta |P (w - g_by)|^2 / 2 over the shared.
                    let (moving, by) = if fa.free.is_empty() { (pb.fix, pa.fix) } else { (pa.fix, pb.fix) };
                    let (fm, fby) = (&s.fixes[moving], &s.fixes[by]);
                    let axes: [bool; 3] = std::array::from_fn(|d| d < ndof && fm.dofs[d] && fby.dofs[d]);
                    if !axes.iter().any(|a| *a) {
                        continue;
                    }
                    let (vm, vb) = (st.values_at(fm, &at), st.values_at(fby, &at));
                    let mut g = vec![0.0; ne];
                    for ((x, gm), gb) in shared.iter().zip(&vm).zip(&vb) {
                        for i in 0..8 {
                            (0..ndof).filter(|&d| axes[d]).for_each(|d| g[i * ndof + d] += x.2 * shape_n(i, x.1) * (gm[d] - gb[d]));
                        }
                    }
                    // The held value's alone: no stiffness on u.
                    let t = hold::penalty(ndof, &axes, &m, &g, beta).held_value_alone();
                    clamps.push(Clamp { other: Some(by), ..clamp(ci, moving, t) });
                }
            }
        }
        // What its held material leaves free (its moments, of unit
        // coefficient), on each axis held: the cell's geometric
        // stiffness is its free material's on those, the material outside
        // the regions that hold that axis (a guide holding a beam along y
        // all through left its supports' cells no material free along z
        // either: it buckled 13% over Euler).
        let mut axes = [false; 3];
        for p in on.iter().filter(|p| p.vol > 0.0) {
            (0..ndof).for_each(|d| axes[d] |= s.fixes[p.fix].dofs[d] || is_absorbed);
        }
        let per_axis: [[f64; 27]; 3] = std::array::from_fn(|d| {
            if !axes[d] {
                return [0.0; 27];
            }
            let f = free_where(&|fixes, mask| (0..fixes.len()).all(|b| mask >> b & 1 == 0 || !s.fixes[fixes[b]].dofs[d]));
            if f[0] > 1e-9 * mat[0] { f.map(|v| v / c.coef) } else { [0.0; 27] }
        });
        (ci, holds.then(|| (kf, (free.map(|v| if any_free { v / c.coef } else { 0.0 }), per_axis, axes))), clamps, bound, true)
    }
}
