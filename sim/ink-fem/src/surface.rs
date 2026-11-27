//! The part's surface as facets: where pressures and other surface loads
//! act, and how much of a region's surface the cells carry.

use crate::cells::Placed;
use crate::eval::{Ctx, eval_interval, eval_points};
use crate::hold;
use crate::octree::{Cell, Grid, root_for};
use crate::tape::Tape;
use crate::{RegionSel, Setup};
use ink_core::{V3, eval, interval, tape, v3};
use rayon::prelude::*;
use std::sync::Arc;

/// The part's surface as the facets of a closed mesh half a cell fine:
/// (centroid, outward normal, area). Facets tile the surface, so an
/// integral over any patch of it sums to the patch's area whatever the
/// grid: samples from separate tangent planes do not tile, and their ripple
/// (1e-3 of a pressure load, with a cos 2 theta part) visibly bent a thin
/// ring.
fn surface_facets(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3, h: f64) -> Vec<Facet> {
    let m = crate::mesh::mesh(shape, ctx, lo, hi, 0.5 * h);
    m.tris.iter().filter_map(|t| Facet::new(t.map(|i| m.verts[i as usize]), None, 0.5 * h)).collect()
}

/// Cells belonging to one connected component of the numerical model.
#[derive(Default)]
pub(crate) struct KeptPiece {
    pub(crate) cells: Vec<[i64; 3]>,
    pub(crate) volume: f64,
}

/// Try to enclose one component in air. Work only in the pair's combined
/// box, padded by one cell. Cross a grid face unless the original shape
/// is strictly positive over its entire area. Escaping the box is
/// inconclusive: a missed web could run outside and come back. Trying
/// the inner component also handles a pin surrounded by its housing.
pub(crate) fn enclosed_apart(shape: &Tape, ctx: &Ctx, grid: &Grid, a: &KeptPiece, b: &KeptPiece) -> bool {
    use std::collections::HashSet;
    let lo: [i64; 3] = std::array::from_fn(|d| a.cells.iter().chain(&b.cells).map(|i| i[d]).min().unwrap() - 1);
    let hi: [i64; 3] = std::array::from_fn(|d| a.cells.iter().chain(&b.cells).map(|i| i[d]).max().unwrap() + 1);
    let enclosed = |inside: &KeptPiece, outside: &KeptPiece| {
        let target: HashSet<_> = outside.cells.iter().copied().collect();
        let mut seen: HashSet<_> = inside.cells.iter().copied().collect();
        if !seen.is_disjoint(&target) {
            return false;
        }
        let mut todo = inside.cells.clone();
        while let Some(i) = todo.pop() {
            let bounds = grid.cell_box(Cell { i, size: 1 });
            for axis in 0..3 {
                for sign in [-1, 1] {
                    let mut j = i;
                    j[axis] += sign;
                    if seen.contains(&j) {
                        continue;
                    }
                    let mut face = bounds;
                    face[axis] = interval::I::point(if sign < 0 { bounds[axis].lo } else { bounds[axis].hi });
                    if eval_interval(shape, ctx, face, None).lo > 0.0 {
                        continue;
                    }
                    if target.contains(&j) || j[axis] < lo[axis] || j[axis] > hi[axis] {
                        return false;
                    }
                    seen.insert(j);
                    todo.push(j);
                }
            }
        }
        true
    };
    let (small, large) = if a.cells.len() <= b.cells.len() { (a, b) } else { (b, a) };
    enclosed(small, large) || enclosed(large, small)
}

/// When the shape's interval on every lattice cell lies on one side of
/// zero, its boundary is the exposed faces of the occupied cells. Tile
/// those faces directly: contouring them rounds corners and changes the
/// pressure applied at the edge of a symmetry cut.
pub(crate) fn lattice_surface(s: &Setup, ctx: &Ctx) -> Option<Vec<Facet>> {
    // Decide on the claim's reference grid even when its comparison grid
    // is too coarse to line up with the faces (or happens to line up better).
    let step = s.pressure_grid;
    if (0..3).map(|a| ((s.hi[a] - s.lo[a]) / step).ceil() + 2.0).product::<f64>() > 400_000_000.0 {
        return None;
    }
    let lo = s.lo.map(|x| ((x / step).floor() - 1.0) * step);
    let hi = s.hi.map(|x| ((x / step).ceil() + 1.0) * step);
    let (grid, root) = root_for(lo, hi, step, 1);
    // A box touching both a face and the centre of an absolute value
    // can have a flat derivative bound. A shorter interval separates the
    // zero at the face from the interior kink. An unresolved zero plateau
    // falls back to the contour on every comparison grid.
    fn filled(t: &Tape, ctx: &Ctx, bounds: [interval::I; 3], depth: u8) -> bool {
        if t.zeros_thin(&eval::eval_interval_grad_every(t, ctx, bounds)) {
            return true;
        }
        depth > 0
            && (0..8).all(|corner| {
                let child = std::array::from_fn(|a| {
                    let r = bounds[a];
                    let mid = (r.lo + r.hi) * 0.5;
                    if corner & (1 << a) == 0 { interval::I::new(r.lo, mid) } else { interval::I::new(mid, r.hi) }
                });
                filled(t, ctx, child, depth - 1)
            })
    }
    fn boundary(t: &Tape, ctx: &Ctx, grid: &Grid, cell: Cell, cells: &mut Vec<[i64; 3]>) -> Option<()> {
        let range = eval_interval(t, ctx, grid.cell_box(cell), None);
        if range.lo >= 0.0 || range.hi < 0.0 {
            return Some(());
        }
        if cell.size == 1 {
            // A nonpositive interval alone also admits a plateau at zero
            // that occupies part of the cell. Only tile proved material.
            if !(range.hi <= 0.0 && filled(t, ctx, grid.cell_box(cell), 2)) {
                return None;
            }
            cells.push(cell.i);
        } else {
            for child in cell.children() {
                boundary(t, ctx, grid, child, cells)?;
            }
        }
        Some(())
    }
    let mut cells = Vec::new();
    boundary(&s.shape, ctx, &grid, root, &mut cells)?;
    // Split each exact face at the solve's own cell lines, even for a
    // non-nested --res override. Its nodal weights then lie in one cell.
    let spans = |low: f64, high: f64| {
        let first = (low / s.h).floor() as i64 + 1;
        let end = (high / s.h).ceil() as i64;
        let mut cuts = vec![low];
        cuts.extend((first..end).map(|i| i as f64 * s.h).filter(|x| *x > low && *x < high));
        cuts.push(high);
        cuts
    };
    let mut facets = Vec::new();
    for idx in cells {
        for axis in 0..3 {
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            for sign in [-1, 1] {
                let mut neighbor = idx;
                neighbor[axis] += sign;
                if eval_interval(&s.shape, ctx, grid.cell_box(Cell { i: neighbor, size: 1 }), None).lo < 0.0 {
                    continue;
                }
                let mut corner = idx;
                corner[axis] += i64::from(sign > 0);
                let origin = grid.pos(corner);
                let mut opposite = corner;
                opposite[u] += 1;
                opposite[v] += 1;
                let end = grid.pos(opposite);
                let rows = spans(origin[u], end[u]);
                let cols = spans(origin[v], end[v]);
                let mut normal = [0.0; 3];
                normal[axis] = sign as f64;
                for row in rows.windows(2) {
                    for col in cols.windows(2) {
                        let corners = [(0, 0), (1, 0), (1, 1), (0, 1)].map(|(a, b)| {
                            let mut p = origin;
                            p[u] = row[a];
                            p[v] = col[b];
                            p
                        });
                        let mid = std::array::from_fn(|k| (corners[0][k] + corners[2][k]) * 0.5);
                        for j in 0..4 {
                            facets.extend(Facet::new([mid, corners[j], corners[(j + 1) % 4]], Some(normal), s.h.min(step)));
                        }
                    }
                }
            }
        }
    }
    Some(facets)
}

/// The half-spaces of a field whose negative region is an intersection
/// of planes over `bounds`, with each plane relative to `mid`.
pub(crate) fn half_spaces(tape: &Tape, ctx: &Ctx, mid: V3, bounds: [interval::I; 3]) -> Option<Vec<[f64; 4]>> {
    use crate::graph::{BinOp, UnOp};
    use crate::tape::Op;

    let ssa = &tape.ssa;
    let over = eval::eval_interval_grad_every(tape, ctx, bounds);
    let constant = |i: usize| {
        let g = &over[i];
        (g.l == 0.0 && g.v.lo == g.v.hi && g.v.lo.is_finite()).then_some(g.v.lo)
    };
    // Each term is sign * instruction + offset <= 0. Only split
    // operations for which this is equivalent to intersecting the terms.
    let mut pending = vec![(ssa.len() - 1, 1.0, 0.0)];
    let mut seen = Vec::new();
    let mut terms = Vec::new();
    while let Some(term @ (i, sign, offset)) = pending.pop() {
        if seen.contains(&term) {
            continue;
        }
        // A shared subtree can be reached with many offsets. Leave
        // those shapes to the contour before the walk grows too large.
        if seen.len() >= 2 * ssa.len() {
            return None;
        }
        seen.push(term);
        let ins = &ssa[i];
        let (a, b) = (ins.a as usize, ins.b as usize);
        match ins.op {
            Op::Tag(_) | Op::Lipschitz => pending.push((a, sign, offset)),
            Op::Bin(BinOp::Max) if sign > 0.0 => pending.extend([(a, sign, offset), (b, sign, offset)]),
            Op::Bin(BinOp::Emax) if sign > 0.0 && offset >= 0.0 => pending.extend([(a, sign, offset), (b, sign, offset)]),
            Op::Un(UnOp::Abs) if sign > 0.0 && offset < 0.0 => pending.extend([(a, sign, offset), (a, -sign, offset)]),
            _ if constant(i).is_some_and(|v| sign * v + offset < 0.0) => {}
            _ if over[i].l.is_finite() && over[i].d.iter().all(|d| d.lo == d.hi && d.lo.is_finite()) => terms.push(term),
            Op::Bin(BinOp::Sub) if constant(b).is_some() => pending.push((a, sign, offset - sign * constant(b)?)),
            Op::Bin(BinOp::Add) if constant(b).is_some() => pending.push((a, sign, offset + sign * constant(b)?)),
            Op::Bin(BinOp::Add) if constant(a).is_some() => pending.push((b, sign, offset + sign * constant(a)?)),
            _ => return None,
        }
    }
    terms.sort_unstable_by_key(|t| t.0);
    let mut ids: Vec<usize> = terms.iter().map(|t| t.0).collect();
    ids.dedup();
    let grads = eval::eval_grads_at(tape, ctx, &[mid], &ids);
    let mut planes: Vec<[f64; 4]> = Vec::new();
    for (i, sign, offset) in terms {
        let [v, x, y, z] = grads[ids.binary_search(&i).ok()?][0];
        let len = v3::norm([x, y, z]);
        let plane = [(sign * v + offset) / len, sign * x / len, sign * y / len, sign * z / len];
        if !plane.iter().all(|v| v.is_finite()) {
            return None;
        }
        // The tightest constraint in each direction owns the face.
        match planes.iter_mut().find(|p| p[1..] == plane[1..]) {
            Some(p) => p[0] = p[0].max(plane[0]),
            None => planes.push(plane),
        }
    }
    Some(planes)
}

/// A convex polyhedron's pressure faces. The contour can join two faces
/// across a taper before their true intersection, or miss the tip when
/// neither face crosses a mesh edge there. Intersect the planes themselves
/// and tile those polygons for the same pressure quadrature as the contour.
pub(crate) fn polyhedral_surface(s: &Setup, ctx: &Ctx) -> Option<Vec<Facet>> {
    let mid = std::array::from_fn(|a| s.lo[a] + 0.5 * (s.hi[a] - s.lo[a]));
    let side = (0..3).map(|a| s.hi[a] - s.lo[a]).fold(0.0, f64::max) + 2.0 * s.h;
    if !side.is_finite() || (0..3).any(|a| mid[a] + s.h == mid[a]) {
        return None;
    }
    let bounds = mid.map(|x| interval::I::new(x - 0.5 * side, x + 0.5 * side));
    let planes = half_spaces(&s.shape, ctx, mid, bounds)?;
    let mut facets = Vec::new();
    for (i, plane) in planes.iter().enumerate() {
        let normal = [plane[1], plane[2], plane[3]];
        let mut face = hold::plane_piece(plane, side, &[], &|_, _| true);
        for clip in planes[..i].iter().chain(&planes[i + 1..]) {
            face = hold::clip(&face, clip, 0.0);
        }
        let mut todo = vec![face];
        while let Some(poly) = todo.pop() {
            if poly.len() < 3 {
                continue;
            }
            let limits: [(f64, f64); 3] = std::array::from_fn(|a| poly.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p[a]), hi.max(p[a]))));
            let (axis, &(lo, hi)) = limits.iter().enumerate().max_by(|(_, a), (_, b)| (a.1 - a.0).total_cmp(&(b.1 - b.0)))?;
            // Keep the facet and region quadrature local to the grid,
            // including where the face cuts diagonally across many cells.
            if hi - lo > 0.5 * s.h {
                let at = lo + 0.5 * (hi - lo);
                todo.extend([-1.0, 1.0].map(|way| clip_poly(&poly, axis, at, way)));
            } else {
                facets.extend((2..poly.len()).filter_map(|k| Facet::new([poly[0], poly[k - 1], poly[k]].map(|p| v3::add(mid, p)), Some(normal), 0.5 * s.h)));
            }
        }
    }
    (!facets.is_empty()).then_some(facets)
}

/// The part of each planar facet in a region of planar half-spaces.
pub(crate) fn clip_facets(facets: &[Facet], mid: V3, planes: &[[f64; 4]], tie: f64) -> Vec<Facet> {
    facets
        .par_iter()
        .flat_map_iter(|f| {
            let corners = f.corners.map(|p| v3::sub(p, mid));
            if planes.iter().all(|p| corners.iter().all(|c| p[0] + v3::dot([p[1], p[2], p[3]], *c) <= tie)) {
                return vec![*f];
            }
            let mut poly = corners.to_vec();
            for p in planes {
                poly = hold::clip(&poly, p, tie);
                if poly.len() < 3 {
                    break;
                }
            }
            (2..poly.len()).filter_map(|k| Facet::new([poly[0], poly[k - 1], poly[k]].map(|p| v3::add(mid, p)), Some(f.normal), f.step)).collect()
        })
        .collect()
}

/// A facet of the part's surface (`surface_facets`): its centroid,
/// outward normal, area and corners, and the step of the mesh it is of.
#[derive(Clone, Copy)]
pub(crate) struct Facet {
    pub(crate) at: V3,
    pub(crate) normal: V3,
    pub(crate) area: f64,
    pub(crate) corners: [V3; 3],
    pub(crate) step: f64,
}

impl Facet {
    /// The triangle's facet, its normal its own or the one it was cut from.
    fn new([a, b, c]: [V3; 3], normal: Option<V3>, step: f64) -> Option<Facet> {
        let n = v3::normal([a, b, c]);
        let len = v3::norm(n);
        (len > 0.0).then(|| Facet {
            at: std::array::from_fn(|k| (a[k] + b[k] + c[k]) / 3.0),
            normal: normal.unwrap_or(n.map(|x| x / len)),
            area: 0.5 * len,
            corners: [a, b, c],
            step,
        })
    }

    /// Its part in the box `lo..hi`, as facets of the same plane.
    pub(crate) fn clipped(&self, lo: V3, hi: V3) -> Vec<Facet> {
        let mut poly = self.corners.to_vec();
        for a in 0..3 {
            poly = clip_poly(&poly, a, lo[a], 1.0);
            poly = clip_poly(&poly, a, hi[a], -1.0);
        }
        (2..poly.len()).filter_map(|k| Facet::new([poly[0], poly[k - 1], poly[k]], Some(self.normal), self.step)).collect()
    }

    /// Its surface measure along the axes: the area it shows looking along
    /// each, summed (sum |n_a| dA), as lines along them count it
    /// (`surface_measure`).
    pub(crate) fn axial(&self) -> f64 {
        self.area * self.normal.iter().map(|n| n.abs()).sum::<f64>()
    }
}

/// The cosine under which two of a facet's samples' surface normals make it a rim's (`rim_pieces`): 30
/// degrees apart, past a smooth face's turn across a facet.
pub(crate) const RIM: f64 = 0.866;

/// Each of `pts` carried along `n` onto the surface of `shape`: the field's
/// zero, bracketed within `reach` either way (out from inside, in from
/// outside) and halved to rounding. A facet's points lie a little off the
/// surface where its corners were drawn in beside an edge.
pub(crate) fn onto_surface(shape: &Tape, ctx: &Ctx, pts: &[V3], n: V3, reach: f64) -> Vec<V3> {
    let f0 = crate::eval::eval_points_par(shape, ctx, pts);
    let (mut lo, mut hi): (Vec<f64>, Vec<f64>) = f0.iter().map(|v| if *v < 0.0 { (0.0, reach) } else { (-reach, 0.0) }).unzip();
    let at = |c: &V3, t: f64| -> V3 { v3::add(*c, v3::scale(n, t)) };
    for _ in 0..SURFACE_HALVINGS {
        if lo.iter().zip(&hi).all(|(a, b)| b - a <= 4.0 * f64::EPSILON * (1.0 + reach)) {
            break;
        }
        let mid: Vec<V3> = pts.iter().zip(lo.iter().zip(&hi)).map(|(c, (a, b))| at(c, 0.5 * (a + b))).collect();
        let fm = crate::eval::eval_points_par(shape, ctx, &mid);
        for (k, v) in fm.iter().enumerate() {
            let m = 0.5 * (lo[k] + hi[k]);
            // Material below the zero: inside moves the bracket out along n, outside in.
            if (*v < 0.0) == (f0[k] < 0.0) {
                if f0[k] < 0.0 { lo[k] = m } else { hi[k] = m }
            } else if f0[k] < 0.0 {
                hi[k] = m
            } else {
                lo[k] = m
            }
        }
    }
    pts.iter().zip(lo.iter().zip(&hi)).map(|(c, (a, b))| at(c, 0.5 * (a + b))).collect()
}

/// Of facet `f` across a rim, the pieces whose points on the surface region `r` holds: the facet halved
/// `RIM_DEPTH` times, each piece's centre projected along the facet's normal onto the surface of `shape`
/// (the field's zero, bracketed within the facet's step either way and halved to rounding), each with its
/// point there and its area on the facet (a pressure takes the surface's area from it).
pub(crate) fn rim_pieces(f: &Facet, shape: &Tape, ctx: &Ctx, r: &RegionSel, h: f64) -> Vec<(V3, f64)> {
    let mut tris = vec![f.corners];
    for _ in 0..RIM_DEPTH {
        tris = tris
            .into_iter()
            .flat_map(|c| {
                let mid = |p: V3, q: V3| -> V3 { std::array::from_fn(|k| 0.5 * (p[k] + q[k])) };
                let (ab, bc, ca) = (mid(c[0], c[1]), mid(c[1], c[2]), mid(c[2], c[0]));
                [[c[0], ab, ca], [ab, c[1], bc], [ca, bc, c[2]], [ab, bc, ca]]
            })
            .collect();
    }
    let area = f.area / tris.len() as f64;
    let centres: Vec<V3> = tris.iter().map(|c| std::array::from_fn(|k| (c[0][k] + c[1][k] + c[2][k]) / 3.0)).collect();
    let on = onto_surface(shape, ctx, &centres, f.normal, f.step.max(1e-9 * h));
    let held = r.holds_on_surface(ctx, &on, h);
    on.into_iter().zip(held).filter(|(_, held)| *held).map(|(p, _)| (p, area)).collect()
}

/// Halvings of a rim facet `rim_pieces` takes its pieces from: 256 pieces.
const RIM_DEPTH: u32 = 4;

/// Of each facet the region holds whole (`whole`, not a rim's), the strips of its face that lie past
/// its edges where the mesh ends short of the face or past it: each strip's point and its area,
/// negative where the edge lies past the face. The mesh's rim is a polygon on a curved edge (a chord
/// whose lens the face holds and no facet), and its vertices on a sharp edge are where its cells put
/// them, off the edge (a quarter tube's end ended 6 to 17 microns short of its cuts): a tube's end took
/// its pressure 0.14% short, its quarter 0.26%. An edge is the mesh's end on the face where a point just
/// past it lies on no other facet of that face, pressed or not (where a region ends within a face, the
/// facets past it are the face's, their share the region's: no end of the mesh); there the face's reach past the edge is found at
/// its ends and its middle (`face_reach`), and the strip's area is Simpson's rule of the three: exact
/// for a parabola, so for an arc's lens but for the fifth power of the chord.
pub(crate) fn face_strips(pts: &[Facet], whole: &[bool], shape: &Tape, ctx: &Ctx, r: &RegionSel, h: f64) -> Vec<Vec<(V3, f64)>> {
    use std::collections::HashMap;
    let step = pts.iter().map(|f| f.step).fold(0.0, f64::max).max(1e-9 * h);
    let key = |p: V3| -> [i64; 3] { p.map(|x| (x / step).floor() as i64) };
    let mut near: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (i, f) in pts.iter().enumerate() {
        let (a, b) = (
            key(std::array::from_fn(|k| f.corners.iter().map(|c| c[k]).fold(f64::INFINITY, f64::min))),
            key(std::array::from_fn(|k| f.corners.iter().map(|c| c[k]).fold(f64::NEG_INFINITY, f64::max))),
        );
        for x in a[0]..=b[0] {
            for y in a[1]..=b[1] {
                for z in a[2]..=b[2] {
                    near.entry([x, y, z]).or_default().push(i);
                }
            }
        }
    }
    // Whether facet j, of the same face (turned less than RIM), holds point q in its triangle (within
    // its plane's rounding). A rim's sliver along a sharp edge (a tube's hole) is of neither face.
    let holds = |j: usize, n: V3, q: V3| -> bool {
        let f = &pts[j];
        if v3::dot(f.normal, n) < RIM || v3::dot(v3::sub(q, f.corners[0]), f.normal).abs() > 0.5 * f.step {
            return false;
        }
        let [a, b, c] = f.corners;
        // (about the triangle's own turn: a lattice face's corners run either way round its normal)
        let turn = v3::cross(v3::sub(b, a), v3::sub(c, a));
        let side = |p: V3, s: V3| v3::dot(v3::cross(v3::sub(s, p), v3::sub(q, p)), turn) / (v3::norm(v3::sub(s, p)) * v3::norm(turn)).max(f64::MIN_POSITIVE);
        let e = 1e-9 * f.step;
        side(a, b) >= -e && side(b, c) >= -e && side(c, a) >= -e
    };
    // Each end of the mesh on a face: (facet, edge's ends, its outward direction in the facet's plane).
    let (near, holds) = (&near, &holds);
    let ends: Vec<(usize, V3, V3, V3)> = (0..pts.len())
        .into_par_iter()
        .filter(|&i| whole[i])
        .flat_map_iter(|i| {
            let f = &pts[i];
            (0..3).filter_map(move |k| {
                let (a, b, c) = (f.corners[k], f.corners[(k + 1) % 3], f.corners[(k + 2) % 3]);
                let mut out = v3::direction(v3::cross(v3::sub(b, a), f.normal))?;
                if v3::dot(out, v3::sub(c, a)) > 0.0 {
                    out = v3::scale(out, -1.0);
                }
                let mid = v3::scale(v3::add(a, b), 0.5);
                let q = v3::add(mid, v3::scale(out, 1e-6 * f.step));
                let held = near.get(&key(q)).is_some_and(|js| js.iter().any(|&j| j != i && holds(j, f.normal, q)));
                (!held).then_some((i, a, b, out))
            })
        })
        .collect();
    let mut strips = vec![Vec::new(); pts.len()];
    if ends.is_empty() {
        return strips;
    }
    // The face's reach past each edge at its two ends and its middle.
    let at: Vec<(V3, V3, V3, f64)> = ends.iter().flat_map(|&(i, a, b, out)| [a, v3::scale(v3::add(a, b), 0.5), b].map(|p| (p, out, pts[i].normal, pts[i].step))).collect();
    let reach = face_reach(&at, shape, ctx, r, h);
    for (e, &(i, a, b, out)) in ends.iter().enumerate() {
        let [s0, sm, s1] = [reach[3 * e], reach[3 * e + 1], reach[3 * e + 2]];
        let area = v3::norm(v3::sub(b, a)) * (s0 + 4.0 * sm + s1) / 6.0;
        if area != 0.0 {
            strips[i].push((v3::add(v3::scale(v3::add(a, b), 0.5), v3::scale(out, 0.5 * sm)), area));
        }
    }
    strips
}

/// How far the face of each point's facet (`at`: the point, the outward direction in the facet's
/// plane, its normal and its step) reaches past the point along that direction: negative where the
/// point lies past the face's end. A point is on the face where the surface along the facet's normal
/// within its step is found, turned less than RIM from it, and the region holds it there. The reach is
/// halved to rounding within a step either way; 0 where the face does not end within it (a crack in the
/// mesh, not its end).
fn face_reach(at: &[(V3, V3, V3, f64)], shape: &Tape, ctx: &Ctx, r: &RegionSel, h: f64) -> Vec<f64> {
    let on_face = |ts: &[f64]| -> Vec<bool> {
        let q: Vec<V3> = at.iter().zip(ts).map(|((p, out, _, _), t)| v3::add(*p, v3::scale(*out, *t))).collect();
        // The surface along the normal, bracketed within the step either way and halved to rounding.
        let (lo, hi): (Vec<V3>, Vec<V3>) = q.iter().zip(at).map(|(q, (_, _, n, s))| (v3::sub(*q, v3::scale(*n, *s)), v3::add(*q, v3::scale(*n, *s)))).unzip();
        let (flo, fhi) = (crate::eval::eval_points_par(shape, ctx, &lo), crate::eval::eval_points_par(shape, ctx, &hi));
        let found: Vec<bool> = flo.iter().zip(&fhi).map(|(a, b)| *a <= 0.0 && *b > 0.0).collect();
        let (mut a, mut b): (Vec<f64>, Vec<f64>) = at.iter().map(|(_, _, _, s)| (-s, *s)).unzip();
        for _ in 0..SURFACE_HALVINGS {
            let mid: Vec<V3> = q.iter().zip(at).zip(a.iter().zip(&b)).map(|((q, (_, _, n, _)), (a, b))| v3::add(*q, v3::scale(*n, 0.5 * (a + b)))).collect();
            for (k, v) in crate::eval::eval_points_par(shape, ctx, &mid).into_iter().enumerate() {
                let m = 0.5 * (a[k] + b[k]);
                if v <= 0.0 { a[k] = m } else { b[k] = m }
            }
        }
        let on: Vec<V3> = q.iter().zip(at).zip(a.iter().zip(&b)).map(|((q, (_, _, n, _)), (a, b))| v3::add(*q, v3::scale(*n, 0.5 * (a + b)))).collect();
        let grads = crate::eval::eval_grads(shape, ctx, &on);
        let held = r.holds_on_surface(ctx, &on, h);
        (0..at.len()).map(|k| found[k] && held[k] && v3::direction([grads[k][1], grads[k][2], grads[k][3]]).is_some_and(|g| v3::dot(g, at[k].2) >= RIM)).collect()
    };
    let start = on_face(&vec![0.0; at.len()]);
    // From the point toward the face's end: outward from on it, inward from past it.
    let far: Vec<f64> = at.iter().zip(&start).map(|((_, _, _, s), on)| if *on { *s } else { -s }).collect();
    let ends = on_face(&far);
    let (mut a, mut b): (Vec<f64>, Vec<f64>) = (vec![0.0; at.len()], far.clone());
    for _ in 0..SURFACE_HALVINGS {
        let mid: Vec<f64> = a.iter().zip(&b).map(|(a, b)| 0.5 * (a + b)).collect();
        for (k, on) in on_face(&mid).into_iter().enumerate() {
            if on == start[k] { a[k] = mid[k] } else { b[k] = mid[k] }
        }
    }
    (0..at.len()).map(|k| if ends[k] == start[k] { 0.0 } else { 0.5 * (a[k] + b[k]) }).collect()
}

/// The least cosine between a facet's normal and the surface's at a sample for a pressure there to
/// take the surface's normal and area (`patch / cos`): past it (a facet across a sharp edge, its sample
/// on the other face) the facet's own, so no sample's area grows more than threefold.
pub(crate) const ON_SURFACE: f64 = 1.0 / 3.0;

/// A facet's six points and each one's area of it: the triangle rule of
/// degree four (Dunavant), which integrates a trilinear cell shape function
/// on a planar facet within one cell, and that times an affine pressure.
pub(crate) fn six_points(f: &Facet) -> impl Iterator<Item = (V3, f64)> + '_ {
    [(0.108103018168070, 0.223381589678011), (0.816847572980459, 0.109951743655322)].into_iter().flat_map(move |(corner, weight)| {
        (0..3).map(move |j| (std::array::from_fn(|a| corner * f.corners[j][a] + 0.5 * (1.0 - corner) * (f.corners[(j + 1) % 3][a] + f.corners[(j + 2) % 3][a])), weight * f.area))
    })
}

/// A facet's pieces, each within one cell: cut at the planes of the grid's
/// cells (`origin`, side `h`), and a piece whose cell is smaller (`leaf`:
/// the side of the cell taking a point of the surface, None where none
/// does) cut again at the planes of that size. A cell's shape functions
/// are trilinear within it and bend at its planes, so what is integrated
/// against them over the surface is integrated piece by piece.
pub(crate) fn facet_in_cells(f: &Facet, origin: V3, h: f64, leaf: &dyn Fn(V3) -> Option<f64>) -> Vec<Facet> {
    let (mut out, mut todo) = (Vec::new(), vec![(*f, h)]);
    while let Some((piece, size)) = todo.pop() {
        let at = |way: fn(f64, f64) -> f64, start: f64| -> [i64; 3] {
            std::array::from_fn(|a| ((piece.corners.iter().map(|p| p[a]).fold(start, way) - origin[a]) / size).floor() as i64)
        };
        let (i0, i1) = (at(f64::min, f64::INFINITY), at(f64::max, f64::NEG_INFINITY));
        for i in (i0[0]..=i1[0]).flat_map(|x| (i0[1]..=i1[1]).flat_map(move |y| (i0[2]..=i1[2]).map(move |z| [x, y, z]))) {
            let lo: V3 = std::array::from_fn(|a| origin[a] + i[a] as f64 * size);
            for part in piece.clipped(lo, lo.map(|v| v + size)) {
                match leaf(part.at) {
                    Some(side) if side < 0.75 * size => todo.push((part, side)),
                    _ => out.push(part),
                }
            }
        }
    }
    out
}

/// A planar facet's parts where `field` is not positive (`within`) and
/// where it is (`without`), cut along the field's zero: straight between
/// the zeros on its sides, exact where the field is affine over the facet
/// (its value at the facet's centre the mean of its corners', to a
/// billionth of their size). Where it is not (a curved region's edge, the
/// crease of a max of two) the facet is taken in four and each again,
/// `depth` times, and the last by its sides' zeros found by halving: exact
/// still where the zero runs straight through it (a crease beside it).
pub(crate) fn facet_cut(f: &Facet, field: &dyn Fn(V3) -> f64, depth: u8, within: &mut Vec<Facet>, without: &mut Vec<Facet>) {
    let v = f.corners.map(field);
    let affine = (field(f.at) - (v[0] + v[1] + v[2]) / 3.0).abs() <= 1e-9 * v.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    if !affine && depth > 0 {
        let mid: [V3; 3] = std::array::from_fn(|k| std::array::from_fn(|a| 0.5 * (f.corners[k][a] + f.corners[(k + 1) % 3][a])));
        for corners in [[f.corners[0], mid[0], mid[2]], [mid[0], f.corners[1], mid[1]], [mid[2], mid[1], f.corners[2]], [mid[0], mid[1], mid[2]]] {
            if let Some(part) = Facet::new(corners, Some(f.normal), f.step) {
                facet_cut(&part, field, depth - 1, within, without);
            }
        }
        return;
    }
    // Where the field crosses zero from a (not positive at one end) to b.
    let zero = |a: V3, b: V3, va: f64, vb: f64| -> V3 {
        let along = |s: f64| -> V3 { std::array::from_fn(|i| a[i] + s * (b[i] - a[i])) };
        if affine {
            return along(va / (va - vb));
        }
        let (mut lo, mut hi) = if va <= 0.0 { (0.0, 1.0) } else { (1.0, 0.0) };
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if field(along(mid)) <= 0.0 { lo = mid } else { hi = mid }
        }
        along(0.5 * (lo + hi))
    };
    let crossings: [Option<V3>; 3] = std::array::from_fn(|k| ((v[k] <= 0.0) != (v[(k + 1) % 3] <= 0.0)).then(|| zero(f.corners[k], f.corners[(k + 1) % 3], v[k], v[(k + 1) % 3])));
    for (keep, out) in [(true, within), (false, without)] {
        let mut poly: Vec<V3> = Vec::new();
        for k in 0..3 {
            if (v[k] <= 0.0) == keep {
                poly.push(f.corners[k]);
            }
            poly.extend(crossings[k]);
        }
        out.extend((2..poly.len()).filter_map(|k| Facet::new([poly[0], poly[k - 1], poly[k]], Some(f.normal), f.step)));
    }
}

/// How many times `facet_cut` quarters a facet a field is not affine over.
pub(crate) const CUT_DEPTH: u8 = 3;

/// Integrate a varying pressure on a complete planar facet. The six-point
/// triangle rule (Dunavant, degree four) integrates an affine pressure
/// times a trilinear cell shape function. Bound the departure from an
/// affine field by the range of each derivative times the facet span.
/// Both the integral and the positive-weight rule can differ from the
/// affine integral by this remainder. A jump has only the value bound.
pub(crate) fn integrate_pressure(pt: &Tape, ctx: &Ctx, facet: &Facet, tolerance: f64, subdivide: bool, samples: &mut Vec<(V3, f64, f64)>) -> bool {
    // Keep each clipped convex polygon intact until quadrature. Triangulating
    // after every cut repeatedly splits the same narrow strip at its kinks.
    const MAX_VISITS: usize = 256;
    let mut todo = vec![facet.corners.to_vec()];
    let mut choices = vec![tape::BOTH; pt.len()];
    let mut visits = 0;
    while let Some(poly) = todo.pop() {
        visits += 1;
        if visits > MAX_VISITS {
            return false;
        }
        let triangles: Vec<_> = (2..poly.len()).filter_map(|k| Facet::new([poly[0], poly[k - 1], poly[k]], Some(facet.normal), facet.step)).collect();
        let area: f64 = triangles.iter().map(|p| p.area).sum();
        if area == 0.0 {
            continue;
        }
        // Half the error allowance is proportional to area. The other
        // half is divided among at most MAX_VISITS leaves, so a vanishing
        // sliver at a kink need not resolve its pointwise derivative.
        let allowed = 0.5 * tolerance * (1.0 + facet.area / (area * MAX_VISITS as f64));
        let bounds = std::array::from_fn(|a| {
            let x = poly.iter().map(|p| p[a]);
            interval::I::new(x.clone().fold(f64::INFINITY, f64::min), x.fold(f64::NEG_INFINITY, f64::max))
        });
        let p = eval_interval(pt, ctx, bounds, Some(&mut choices));
        let mut error = p.hi - p.lo;
        let mut blur = bounds.map(|r| r.hi - r.lo);
        if error > allowed && pt.continuous(pt.len() - 1, &choices) {
            let d = eval::eval_interval_grad(pt, ctx, bounds).d;
            blur = std::array::from_fn(|a| if blur[a] == 0.0 { 0.0 } else { blur[a] * (d[a].hi - d[a].lo) });
            error = error.min(2.0 * blur.iter().sum::<f64>());
        }
        if error <= allowed {
            samples.extend(triangles.iter().flat_map(six_points).map(|(point, area)| (point, eval::eval_point(pt, ctx, point), area)));
        } else {
            if !subdivide {
                return false;
            }
            let a = (0..3).max_by(|&a, &b| blur[a].total_cmp(&blur[b])).unwrap();
            let mid = (bounds[a].lo + bounds[a].hi) * 0.5;
            if !(mid > bounds[a].lo && mid < bounds[a].hi) {
                return false;
            }
            for side in [-1.0, 1.0] {
                let part = clip_poly(&poly, a, mid, side);
                if part.len() >= 3 {
                    todo.push(part);
                }
            }
        }
    }
    true
}

/// The part of a polygon on the side `way` of the plane `x[axis] = at`
/// (Sutherland-Hodgman).
fn clip_poly(poly: &[V3], axis: usize, at: f64, way: f64) -> Vec<V3> {
    let mut out = Vec::with_capacity(poly.len() + 1);
    for (k, p) in poly.iter().enumerate() {
        let q = poly[(k + 1) % poly.len()];
        let (dp, dq) = (way * (p[axis] - at), way * (q[axis] - at));
        if dp >= 0.0 {
            out.push(*p);
        }
        if (dp >= 0.0) != (dq >= 0.0) {
            let t = dp / (dp - dq);
            out.push(std::array::from_fn(|a| p[a] + t * (q[a] - p[a])));
        }
    }
    out
}

/// Split grid cells and those beside them, at the finest level among the
/// split ones beside each (a face can lie in a cell holding only a sliver,
/// dropped, beside a split one), sorted by index (`part_surface`).
pub(crate) fn split_cells<T: Placed>(cells: &[T], dims: [i64; 3]) -> Vec<([i64; 3], u8)> {
    let mut split: Vec<([i64; 3], u8)> = cells
        .iter()
        .filter(|c| c.level() > 0)
        .flat_map(|c| {
            let (i, l) = (c.idx(), c.level());
            (0..27).map(move |d| ([i[0] + d % 3 - 1, i[1] + d / 3 % 3 - 1, i[2] + d / 9 - 1], l))
        })
        .filter(|c| (0..3).all(|a| (0..dims[a]).contains(&c.0[a])))
        .collect();
    split.sort_unstable();
    split.reverse();
    split.dedup_by_key(|c| c.0);
    split.reverse();
    split
}

/// The part's surface as facets meshed at half the cells they lie in
/// (`surface_facets`): half the grid's cell, and in and beside split grid
/// cells (`split_cells`) half the finest child of those near them (in
/// blocks, `SPLIT_BLOCK`), each mesh cut exactly at the grid cells' faces,
/// so they tile the surface between them. Every integral over the surface
/// (pressures, films, frictionless holds), the
/// test of a symmetry plane and the stress reads take it: meshed at half
/// the grid's cells, a fin 0.189 mm thick lying between the mesh's planes
/// (0.5 mm cells, split to 0.0625) had no facet, and a pressure on it was
/// lost.
pub(crate) fn part_surface(shape: &Arc<Tape>, ctx: &Ctx, grid: &Grid, lo: V3, hi: V3, split: &[([i64; 3], u8)]) -> Vec<Facet> {
    let h = grid.h;
    let index = |p: V3| -> [i64; 3] { std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64) };
    // The grid cell a piece of surface is of: the one holding the material
    // just under its centroid, so a piece on a cell's face (a box's face
    // on the grid's plane) is one cell's, not both sides'.
    let owner = |f: &Facet| index(std::array::from_fn(|a| f.at[a] - 1e-6 * h * f.normal[a]));
    let is_split = |i: [i64; 3]| split.binary_search_by(|c| c.0.cmp(&i)).is_ok();
    // A facet's pieces in the grid cells `keep` takes, each its owner's.
    let pieces = |f: &Facet, keep: &dyn Fn([i64; 3]) -> bool| -> Vec<Facet> {
        let a = index(std::array::from_fn(|k| f.corners.iter().fold(f64::INFINITY, |m, p| m.min(p[k])) - 1e-6 * h));
        let b = index(std::array::from_fn(|k| f.corners.iter().fold(f64::NEG_INFINITY, |m, p| m.max(p[k])) + 1e-6 * h));
        let cells: Vec<[i64; 3]> = (a[0]..=b[0]).flat_map(|x| (a[1]..=b[1]).flat_map(move |y| (a[2]..=b[2]).map(move |z| [x, y, z]))).collect();
        if cells.iter().all(|&i| keep(i)) {
            return vec![*f];
        }
        cells
            .into_iter()
            .filter(|&i| keep(i))
            .flat_map(|i| {
                let (lo, hi) = grid.cell_bounds(Cell { i, size: 1 });
                f.clipped(lo, hi).into_iter().filter(move |p| owner(p) == i)
            })
            .collect()
    };
    let base = surface_facets(shape, ctx, lo, hi, h);
    let mut out: Vec<Facet> = base.par_iter().flat_map_iter(|f| pieces(f, &|i| !is_split(i))).collect();
    // The split cells in blocks of `SPLIT_BLOCK` grid cells to a side, each
    // block's meshed at once at half its finest child, past its split
    // cells' box by a child and cut to them: a box's own mesh closes the
    // surface across its faces, and a mesh of each cell alone, so padded,
    // meshed 8 times its cell's volume where cells split once.
    let mut blocks: Vec<([i64; 3], [i64; 3], u8)> = split.iter().map(|&(i, l)| (i.map(|x| x.div_euclid(SPLIT_BLOCK)), i, l)).collect();
    blocks.sort_unstable();
    let blocks: Vec<&[([i64; 3], [i64; 3], u8)]> = blocks.chunk_by(|a, b| a.0 == b.0).collect();
    out.par_extend(blocks.par_iter().flat_map_iter(|cells| {
        let level = cells.iter().map(|c| c.2).max().unwrap_or(0);
        let step = h / f64::from(1u32 << level);
        let lo: V3 = std::array::from_fn(|a| cells.iter().map(|c| grid.pos(c.1)[a]).fold(f64::INFINITY, f64::min) - step);
        let hi: V3 = std::array::from_fn(|a| cells.iter().map(|c| grid.pos(c.1)[a] + h).fold(f64::NEG_INFINITY, f64::max) + step);
        let ours = |i: [i64; 3]| cells.binary_search_by(|c| c.1.cmp(&i)).is_ok();
        surface_facets(shape, ctx, lo, hi, step).iter().flat_map(|f| pieces(f, &ours)).collect::<Vec<Facet>>()
    }));
    out
}

/// Grid cells to a side of the blocks split cells are meshed in
/// (`part_surface`).
pub(crate) const SPLIT_BLOCK: i64 = 4;

/// Sub-lines across a region's rim, per line, that `surface_measure`
/// counts a line near it by.
const RIM_LINES: usize = 16;

/// What lines through the cells count of a region's surface, per patch of
/// `SPLIT_BLOCK` grid cells to a side (`surface_measure`): the surface
/// measure along the axes (sum |n_a| dA, as `Facet::axial`), the variance
/// of its count (a crossing more or less per line), and the most its
/// rounding at the region's rims may be off by.
#[derive(Clone, Copy, Default)]
pub(crate) struct Patch {
    pub(crate) lines: f64,
    pub(crate) var: f64,
    pub(crate) rim: f64,
}

/// The region's surface measure along the axes by lines, apart from the
/// mesh, per patch (`Patch`): in each grid cell holding or beside the
/// cells' material (`cells`) whose box may hold both the surface and the
/// region, lines along each axis half the cell apart, or half the finest
/// child of the split cells of its block (`part_surface`), sampled halfway
/// between the mesh's planes (the grid's mesh's from `mesh_lo`, the split
/// cells' from the grid's origin), so a sheet between two planes thicker
/// than half their spacing crosses the lines or the mesh, if not both. Each
/// crossing of the surface counts as the line's cross-section where its
/// zero (found by false position) is in the cell's box and the region
/// (Cauchy-Crofton). The grid cells' boxes tile space round the material, so the
/// lines see the surface there, each piece once, a face in a sliver its
/// cells dropped too. A crossing whose cross-section the region's rim may
/// cut (a field region, its value at the crossing within its rise across
/// the cross-section) is counted again on `RIM_LINES` sub-lines across the
/// rim, each its share: a strip a line or two wide is measured, not read
/// as no line or a whole one, to half a sub-line at each rim it crosses;
/// so is one a line runs down the middle of (its rims tie there, and its
/// test at the cross-section's edges differs from the middle's).
#[allow(clippy::too_many_arguments)]
pub(crate) fn surface_measure<T: Placed>(
    shape: &Tape,
    ctx: &Ctx,
    grid: &Grid,
    dims: [i64; 3],
    cells: &[T],
    split: &[([i64; 3], u8)],
    region: &RegionSel,
    mesh_lo: V3,
) -> Vec<([i64; 3], Patch)> {
    let h = grid.h;
    let boxed = |i: [i64; 3], n: i64| -> [crate::interval::I; 3] {
        let (lo, hi) = (grid.pos(i.map(|x| x - n)), grid.pos(i.map(|x| x + 1 + n)));
        std::array::from_fn(|a| crate::interval::I::new(lo[a], hi[a]))
    };
    let meets = |b: [crate::interval::I; 3]| {
        let f = eval_interval(shape, ctx, b, None);
        f.lo <= 0.0 && f.hi >= 0.0 && region.may_touch(ctx, b)
    };
    // The level each grid cell's surface is meshed at: its block's finest
    // split cell's, where it is split or beside one (`part_surface`).
    let mut block_level: Vec<([i64; 3], u8)> = split.iter().map(|&(i, l)| (i.map(|x| x.div_euclid(SPLIT_BLOCK)), l)).collect();
    block_level.sort_unstable();
    block_level.reverse();
    block_level.dedup_by_key(|b| b.0);
    block_level.reverse();
    let level_of = |i: [i64; 3]| -> u8 {
        if split.binary_search_by(|c| c.0.cmp(&i)).is_err() {
            return 0;
        }
        let b = i.map(|x| x.div_euclid(SPLIT_BLOCK));
        block_level.binary_search_by(|c| c.0.cmp(&b)).map_or(0, |k| block_level[k].1)
    };
    let field = match region {
        RegionSel::Field(t) => Some(t),
        RegionSel::Face { .. } => None,
    };
    let mut held: Vec<[i64; 3]> = cells.iter().map(|c| c.idx()).collect();
    held.sort_unstable();
    held.dedup();
    let mut near: Vec<[i64; 3]> = held
        .into_par_iter()
        .filter(|&i| meets(boxed(i, 1)))
        .flat_map_iter(|i| (0..27).map(move |d| [i[0] + d % 3 - 1, i[1] + d / 3 % 3 - 1, i[2] + d / 9 - 1]))
        .filter(|i| (0..3).all(|a| (0..dims[a]).contains(&i[a])))
        .collect();
    near.par_sort_unstable();
    near.dedup();
    let mut per: Vec<([i64; 3], Patch)> = near
        .into_par_iter()
        .filter(|&i| meets(boxed(i, 0)))
        .map(|i| {
            let level = level_of(i);
            let n = 2usize << level;
            let s = h / n as f64;
            let (lo, hi) = (grid.pos(i), grid.pos(i.map(|x| x + 1)));
            // The mesh's planes along each axis lie at multiples of s from
            // `o`; the lines' points halfway between them, from before the
            // box to past it.
            let o = if level == 0 { mesh_lo } else { grid.origin };
            let first: [i64; 3] = std::array::from_fn(|a| ((lo[a] - o[a]) / s - 1.5).floor() as i64);
            let count: [usize; 3] = std::array::from_fn(|a| (((hi[a] - o[a]) / s + 0.5).ceil() as i64 - first[a] + 1) as usize);
            let mut pts: Vec<(V3, usize)> = Vec::new();
            let mut runs: Vec<std::ops::Range<usize>> = Vec::new();
            for axis in 0..3 {
                let (b, d) = ((axis + 1) % 3, (axis + 2) % 3);
                for j in 0..n {
                    for k in 0..n {
                        let start = pts.len();
                        for t in 0..count[axis] {
                            let mut p = lo;
                            p[axis] = o[axis] + (first[axis] as f64 + t as f64 + 0.5) * s;
                            p[b] += (j as f64 + 0.5) * s;
                            p[d] += (k as f64 + 0.5) * s;
                            pts.push((p, axis));
                        }
                        runs.push(start..pts.len());
                    }
                }
            }
            let at: Vec<V3> = pts.iter().map(|p| p.0).collect();
            let mut f = Vec::new();
            eval_points(shape, ctx, &at, &mut f);
            let mut segs: Vec<(V3, f64, V3, f64, usize)> = Vec::new();
            for r in &runs {
                for t in r.start..r.end - 1 {
                    let (fa, fb) = (f[t], f[t + 1]);
                    if fa.is_finite() && fb.is_finite() && (fa < 0.0) != (fb < 0.0) {
                        segs.push(if fa < 0.0 { (at[t], fa, at[t + 1], fb, pts[t].1) } else { (at[t + 1], fb, at[t], fa, pts[t].1) });
                    }
                }
            }
            // Each crossing's zero, and those of this box: the box holding
            // the material just under it, as a facet's is its cell's
            // (`part_surface`), so a face on the box's face is one box's.
            let inside_box = |p: V3| (0..3).all(|a| p[a] >= lo[a] && p[a] < hi[a]);
            let zeros = crossings(shape, ctx, &segs);
            let under = |k: usize, z: V3| -> V3 {
                let (a, axis) = (segs[k].0, segs[k].4);
                let mut q = z;
                q[axis] += 1e-6 * h * (a[axis] - z[axis]).signum();
                q
            };
            let mut patch = Patch::default();
            let mine: Vec<(usize, V3)> = zeros.iter().enumerate().filter(|(k, z)| inside_box(under(*k, **z))).map(|(k, z)| (k, *z)).collect();
            if mine.is_empty() {
                return (i.map(|x| x.div_euclid(SPLIT_BLOCK)), patch);
            }
            let w = s * s;
            let at: Vec<V3> = mine.iter().map(|m| m.1).collect();
            let inside = region.holds_on_surface(ctx, &at, h);
            // Crossings whose cross-section the rim may cut: its sub-lines
            // across it, along the axis across the line it rises most along.
            let rims: Vec<(usize, usize)> = match field {
                Some(t) => {
                    // The region's rise along the surface (its gradient less
                    // its part along the surface's normal), over the patch of
                    // surface the line's cross-section meets (wider the more
                    // the surface slants to the line); or the region's own
                    // test at the cross-section's edges, half a line either
                    // way across it, differing from its centre's: where two
                    // rims tie (a strip's middle down the line) its gradient
                    // along the surface is none.
                    let g = crate::eval::eval_grads(t, ctx, &at);
                    let sn = crate::eval::eval_grads(shape, ctx, &at);
                    let edges: Vec<V3> = mine
                        .iter()
                        .flat_map(|&(k, z)| {
                            let axis = segs[k].4;
                            [(axis + 1) % 3, (axis + 2) % 3].into_iter().flat_map(move |c| {
                                [-0.5, 0.5].map(|u| {
                                    let mut p = z;
                                    p[c] += u * s;
                                    p
                                })
                            })
                        })
                        .collect();
                    let at_edges = region.holds_on_surface(ctx, &edges, h);
                    mine.iter()
                        .zip(g.iter().zip(&sn))
                        .zip(at_edges.chunks(4).zip(&inside))
                        .filter_map(|((&(k, _), (g, sn)), (e, &centre))| {
                            let axis = segs[k].4;
                            let (b, d) = ((axis + 1) % 3, (axis + 2) % 3);
                            let n = v3::direction([sn[1], sn[2], sn[3]]).unwrap_or([0.0; 3]);
                            let gn = g[1] * n[0] + g[2] * n[1] + g[3] * n[2];
                            let gt: V3 = std::array::from_fn(|a| g[1 + a] - gn * n[a]);
                            let reach = 0.5 * s * (gt[b].abs() + gt[d].abs() + gt[axis].abs()) / n[axis].abs().max(0.2);
                            let differ = [e[0..2].iter().filter(|x| **x != centre).count(), e[2..4].iter().filter(|x| **x != centre).count()];
                            let across = match differ {
                                [x, y] if x > y => b,
                                [x, y] if y > x => d,
                                _ if gt[b].abs() >= gt[d].abs() => b,
                                _ => d,
                            };
                            (g[0].abs() <= reach || differ != [0, 0]).then_some((k, across))
                        })
                        .collect()
                }
                None => Vec::new(),
            };
            let near_rim: std::collections::HashSet<usize> = rims.iter().map(|r| r.0).collect();
            for ((k, _), inn) in mine.iter().zip(&inside) {
                if !near_rim.contains(k) && *inn {
                    patch.lines += w;
                    patch.var += w * w;
                }
            }
            if !rims.is_empty() {
                // Each rim crossing's line, moved across to each sub-line's
                // place, sampled a step either side of its segment: the
                // surface's crossing moves along a slanting line as it moves
                // across. Its crossing nearest the line's.
                let sub: Vec<[V3; 4]> = rims
                    .iter()
                    .flat_map(|&(k, across)| {
                        let (a, _, b, _, _) = segs[k];
                        (0..RIM_LINES).map(move |m| {
                            let off = ((m as f64 + 0.5) / RIM_LINES as f64 - 0.5) * s;
                            std::array::from_fn(|t| {
                                let u = t as f64 - 1.0;
                                let mut p: V3 = std::array::from_fn(|c| a[c] + u * (b[c] - a[c]));
                                p[across] += off;
                                p
                            })
                        })
                    })
                    .collect();
                let ends: Vec<V3> = sub.iter().flatten().copied().collect();
                eval_points(shape, ctx, &ends, &mut f);
                let picked: Vec<Option<(V3, f64, V3, f64, usize)>> = sub
                    .iter()
                    .zip(f.chunks(4))
                    .map(|(p, v)| {
                        [1usize, 0, 2].into_iter().find_map(|t| {
                            let (fa, fb) = (v[t], v[t + 1]);
                            (fa.is_finite() && fb.is_finite() && (fa < 0.0) != (fb < 0.0)).then(|| if fa < 0.0 { (p[t], fa, p[t + 1], fb, 0) } else { (p[t + 1], fb, p[t], fa, 0) })
                        })
                    })
                    .collect();
                let crossed: Vec<bool> = picked.iter().map(Option::is_some).collect();
                let sub: Vec<(V3, f64, V3, f64, usize)> = picked.iter().map(|p| p.unwrap_or(([0.0; 3], -1.0, [0.0; 3], 1.0, 0))).collect();
                let z = crossings(shape, ctx, &sub);
                let held = region.holds_on_surface(ctx, &z, h);
                let ws = w / RIM_LINES as f64;
                for line in held.chunks(RIM_LINES).zip(crossed.chunks(RIM_LINES)) {
                    let ins: Vec<bool> = line.0.iter().zip(line.1).map(|(h, c)| *h && *c).collect();
                    let k = ins.iter().filter(|x| **x).count() as f64;
                    patch.lines += k * ws;
                    // Half a sub-line at each rim it crosses.
                    let rims = ins.windows(2).filter(|p| p[0] != p[1]).count() as f64;
                    patch.rim += rims * 0.5 * ws;
                }
            }
            (i.map(|x| x.div_euclid(SPLIT_BLOCK)), patch)
        })
        .collect();
    per.sort_unstable_by_key(|p| p.0);
    let mut out: Vec<([i64; 3], Patch)> = Vec::new();
    for (k, p) in per {
        match out.last_mut() {
            Some(last) if last.0 == k => {
                last.1.lines += p.lines;
                last.1.var += p.var;
                last.1.rim += p.rim;
            }
            _ => out.push((k, p)),
        }
    }
    out.retain(|p| p.1.lines > 0.0 || p.1.rim > 0.0);
    out
}

/// The zeros of the field on segments (inside, its value, outside, its
/// value), by false position (Illinois): exact where the field is linear
/// along the segment, to far below its length in a few steps where it is
/// smooth.
pub(crate) fn crossings(shape: &Tape, ctx: &Ctx, segs: &[(V3, f64, V3, f64, usize)]) -> Vec<V3> {
    let mut c: Vec<(V3, f64, V3, f64, i8)> = segs.iter().map(|s| (s.0, s.1, s.2, s.3, 0)).collect();
    let guess = |c: &[(V3, f64, V3, f64, i8)]| -> Vec<V3> {
        c.iter()
            .map(|&(a, fa, b, fb, _)| {
                let t = if fa == fb { 0.5 } else { (fa / (fa - fb)).clamp(0.0, 1.0) };
                std::array::from_fn(|k| a[k] + t * (b[k] - a[k]))
            })
            .collect()
    };
    let mut f = Vec::new();
    for _ in 0..6 {
        let m = guess(&c);
        eval_points(shape, ctx, &m, &mut f);
        for ((c, m), &v) in c.iter_mut().zip(m).zip(&f) {
            if !v.is_finite() {
                continue;
            }
            if v < 0.0 {
                (c.0, c.1) = (m, v);
                if c.4 < 0 {
                    c.3 *= 0.5;
                }
                c.4 = -1;
            } else {
                (c.2, c.3) = (m, v);
                if c.4 > 0 {
                    c.1 *= 0.5;
                }
                c.4 = 1;
            }
        }
    }
    guess(&c)
}

/// Four equal-area samples on each triangle. On an exact lattice face
/// these already lie on the surface, even for a zero-thickness selector.
pub(crate) fn facet_samples(facets: &[Facet]) -> Vec<V3> {
    const W: [[f64; 3]; 4] = [[4.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0], [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0], [1.0 / 6.0, 1.0 / 6.0, 4.0 / 6.0], [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]];
    facets.iter().flat_map(|f| W.map(|w| std::array::from_fn(|k| w[0] * f.corners[0][k] + w[1] * f.corners[1][k] + w[2] * f.corners[2][k]))).collect()
}

/// The share of the surface patch `tri` (a facet's corners on the surface) that region `t` holds, by
/// area, and that share's centroid (the patch's own where it holds none): halved where the region's
/// bounds over a part's box (thickened by `slack`, how far the surface may lie off the patch) hold its
/// zero, to `SHARE_DEPTH` halvings, a part still undecided there taken by its centre's sign. A region
/// narrower than the facet (a strip of a face) takes its area, not the share of four points it holds.
pub(crate) fn region_share(t: &Tape, ctx: &Ctx, tri: [V3; 3], slack: f64, tie: f64) -> (f64, V3) {
    let area = |c: &[V3; 3]| 0.5 * v3::norm(v3::normal(*c));
    let centre = |c: &[V3; 3]| -> V3 { std::array::from_fn(|k| (c[0][k] + c[1][k] + c[2][k]) / 3.0) };
    let whole = area(&tri);
    if !(whole > 0.0) {
        return (0.0, centre(&tri));
    }
    let (mut held, mut moment) = (0.0, [0.0; 3]);
    let mut todo = vec![(tri, 0u32)];
    while let Some((c, depth)) = todo.pop() {
        let b: [interval::I; 3] = std::array::from_fn(|k| {
            let (lo, hi) = (c[0][k].min(c[1][k]).min(c[2][k]), c[0][k].max(c[1][k]).max(c[2][k]));
            interval::I::new(lo - slack, hi + slack)
        });
        let r = eval_interval(t, ctx, b, None);
        let inside = if r.hi <= tie {
            Some(true)
        } else if r.lo > tie {
            Some(false)
        } else if depth >= SHARE_DEPTH {
            Some(crate::eval::eval_point(t, ctx, centre(&c)) <= tie)
        } else {
            None
        };
        match inside {
            Some(true) => {
                let (a, m) = (area(&c), centre(&c));
                held += a;
                (0..3).for_each(|k| moment[k] += a * m[k]);
            }
            Some(false) => {}
            None => {
                let mid = |p: V3, q: V3| -> V3 { std::array::from_fn(|k| 0.5 * (p[k] + q[k])) };
                let (ab, bc, ca) = (mid(c[0], c[1]), mid(c[1], c[2]), mid(c[2], c[0]));
                todo.extend([([c[0], ab, ca], depth + 1), ([ab, c[1], bc], depth + 1), ([ca, bc, c[2]], depth + 1), ([ab, bc, ca], depth + 1)]);
            }
        }
    }
    if held > 0.0 { ((held / whole).min(1.0), moment.map(|m| m / held)) } else { (0.0, centre(&tri)) }
}

/// Halvings of a facet `region_share` reads a region's edge across it to: its edge then placed to
/// 1/256 of the facet, an area within about 1/128 of the facet's along each edge, taken by sign.
const SHARE_DEPTH: u32 = 8;

/// The most halvings `surface_samples` brackets a facet's surface point with: past a double's precision.
const SURFACE_HALVINGS: usize = 64;

/// Points of the part's surface for each facet: the centroids of its four
/// halved triangles, each moved along the facet's normal onto the surface
/// (the field's zero, bracketed within half the facet's mesh step, a
/// quarter of its cell, and halved to the rounding of its point). A mesh vertex by an edge sits off both faces
/// (0.005 to 0.012 mm at half a cell of 0.5), farther than a region drawn
/// under a face may be deep: tested on the mesh, a facet of that face read
/// as partly outside a region 0.001 deep, and one beyond the part as
/// partly in it.
pub(crate) fn surface_samples(shape: &Tape, ctx: &Ctx, facets: &[Facet]) -> Vec<V3> {
    let pts = facet_samples(facets);
    let dirs: Vec<V3> = facets.iter().flat_map(|f| [f.normal; 4]).collect();
    let at = |s: &[f64]| -> Vec<V3> { pts.iter().zip(&dirs).zip(s).map(|((p, n), s)| v3::add(*p, v3::scale(*n, *s))).collect() };
    let f0 = crate::eval::eval_points_par(shape, ctx, &pts);
    // Toward the surface: out along the normal from inside, in from outside.
    let cap = facets.iter().flat_map(|f| [0.5 * f.step; 4]);
    let far: Vec<f64> = f0.iter().zip(cap).map(|(f, cap)| if *f < 0.0 { cap } else { -cap }).collect();
    let f1 = crate::eval::eval_points_par(shape, ctx, &at(&far));
    // A bracket [lo, hi] of the zero, where the far point is across it.
    let (mut lo, mut hi): (Vec<f64>, Vec<f64>) = far.iter().map(|&s| if s > 0.0 { (0.0, s) } else { (s, 0.0) }).unzip();
    let crossed: Vec<bool> = f0.iter().zip(&f1).map(|(a, b)| a.is_finite() && b.is_finite() && (*a < 0.0) != (*b < 0.0)).collect();
    // Halved until each bracket is a rounding of its point: a region drawn as a slab a micron under a
    // face (NAFEMS's, by benchmarks/standard) holds the face's points then, where twelve halvings left
    // them 6e-5 of a cell off it, on either side, and the region took none of the face.
    for _ in 0..SURFACE_HALVINGS {
        if lo.iter().zip(&hi).zip(&pts).all(|((a, b), p)| b - a <= 4.0 * f64::EPSILON * (1.0 + v3::norm(*p))) {
            break;
        }
        let mid: Vec<f64> = lo.iter().zip(&hi).map(|(a, b)| 0.5 * (a + b)).collect();
        let fm = crate::eval::eval_points_par(shape, ctx, &at(&mid));
        for i in 0..mid.len() {
            if fm[i] < 0.0 {
                lo[i] = mid[i];
            } else {
                hi[i] = mid[i];
            }
        }
    }
    let s: Vec<f64> = (0..pts.len()).map(|i| if crossed[i] { 0.5 * (lo[i] + hi[i]) } else { 0.0 }).collect();
    at(&s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cells::ALL;
    use crate::report::SurfaceMeasure;
    use ink_core::graph;

    /// A grid cell of a test's surface (`Placed`), whole.
    struct Whole([i64; 3]);
    impl Placed for Whole {
        fn idx(&self) -> [i64; 3] {
            self.0
        }
        fn level(&self) -> u8 {
            0
        }
        fn sub(&self) -> [u16; 3] {
            [0; 3]
        }
        fn mask(&self) -> u64 {
            ALL
        }
    }

    #[test]
    fn a_region_s_surface_is_meshed_at_its_cells_and_checked_by_lines() {
        // A sheet 0.189 thick lying on z = 0 off a block, on 0.5 mm cells
        // split to 0.0625 round it, pressed over every face above z 0.0945:
        // 149.1 mm2 of it (the sheet's top 64, the block's top 24, its upper
        // sides 57.7, the sheet's upper sides 3.4). Meshed at half the grid's
        // cells, the mesh's planes at z 0 and 0.25 miss the sheet: its
        // facets cover 85 mm2, and the lines say they fall short.
        let mut g = graph::Graph::new();
        let xyz = [graph::Var::X, graph::Var::Y, graph::Var::Z].map(|v| g.var(v));
        let mut cuboid = |c: V3, half: V3| {
            let sides: Vec<graph::NodeId> = (0..3)
                .flat_map(|a| {
                    let (lo, hi) = (g.constant(c[a] - half[a]), g.constant(c[a] + half[a]));
                    [g.sub(xyz[a], hi), g.sub(lo, xyz[a])]
                })
                .collect();
            sides[1..].iter().fold(sides[0], |m, s| g.max(m, *s))
        };
        let (block, sheet) = (cuboid([2.0, 0.0, 0.0], [2.0, 3.0, 3.0]), cuboid([10.0, 0.0, 0.0945], [10.0, 2.0, 0.0945]));
        let part = g.min(block, sheet);
        let shape = Arc::new(Tape::new(&g, part));
        let (z, at) = (g.var(graph::Var::Z), g.constant(0.0945));
        let above = g.sub(at, z);
        let region = RegionSel::Field(Arc::new(Tape::new(&g, above)));
        let ctx = Ctx::new(0.0, Vec::new());
        let h: f64 = 0.5;
        let (lo, hi): (V3, V3) = ([-0.5, -3.5, -3.5], [20.5, 3.5, 3.5]);
        let first: [i64; 3] = std::array::from_fn(|a| (lo[a] / h).floor() as i64 - 1);
        let dims: [i64; 3] = std::array::from_fn(|a| (hi[a] / h).ceil() as i64 + 1 - first[a]);
        let origin: V3 = std::array::from_fn(|a| first[a] as f64 * h);
        let (grid, _) = root_for(origin, std::array::from_fn(|a| origin[a] + dims[a] as f64 * h), h, 1);
        let index = |p: V3| -> [i64; 3] { std::array::from_fn(|a| ((p[a] - grid.origin[a]) / h).floor() as i64) };
        let cells: Vec<Whole> = (0..dims[0] * dims[1] * dims[2])
            .map(|k| [k % dims[0], k / dims[0] % dims[1], k / (dims[0] * dims[1])])
            .filter(|&i| eval_interval(&shape, &ctx, grid.cell_box(Cell { i, size: 1 }), None).lo < 0.0)
            .map(Whole)
            .collect();
        // The sheet's cells split three times, as a solve splits them.
        let sheet_cells: Vec<SheetCell> = cells.iter().filter(|c| grid.pos(c.0)[0] >= 4.0 && grid.pos(c.0)[2].abs() < 0.25).map(|c| SheetCell(c.0)).collect();
        let split = split_cells(&sheet_cells, dims);
        assert!(!split.is_empty() && index([10.0, 0.0, 0.1]) == split[split.binary_search_by(|c| c.0.cmp(&index([10.0, 0.0, 0.1]))).unwrap()].0);
        let measure = |split: &[([i64; 3], u8)]| {
            let facets = part_surface(&shape, &ctx, &grid, lo, hi, split);
            let samples = surface_samples(&shape, &ctx, &facets);
            let took: Vec<(Facet, f64)> = facets.iter().zip(region.on_surface(&ctx, &facets, &samples, h)).filter(|(_, w)| w.0 > 0.0).map(|(f, w)| (*f, w.0)).collect();
            SurfaceMeasure::of("top", "load", &took, surface_measure(&shape, &ctx, &grid, dims, &cells, split, &region, lo), &grid)
        };
        let area = 64.0 + 24.0 + 2.0 * (6.0 + 4.0) * (3.0 - 0.0945) - 4.0 * 0.0945 + (2.0 * 16.0 + 4.0) * 0.0945;
        let fine = measure(&split);
        assert!((fine.lines - area).abs() < 0.01 * area && (fine.covered - area).abs() < 0.01 * area && fine.shortfall == 0.0, "the split cells' mesh: {fine:?}, the faces {area}");
        // The grid's mesh alone, the lines as they lie for the split
        // cells: the sheet's top (64) is missed, and the lines say so.
        let facets = part_surface(&shape, &ctx, &grid, lo, hi, &[]);
        let samples = surface_samples(&shape, &ctx, &facets);
        let took: Vec<(Facet, f64)> = facets.iter().zip(region.on_surface(&ctx, &facets, &samples, h)).filter(|(_, w)| w.0 > 0.0).map(|(f, w)| (*f, w.0)).collect();
        let coarse = SurfaceMeasure::of("top", "load", &took, surface_measure(&shape, &ctx, &grid, dims, &cells, &split, &region, lo), &grid);
        assert!(coarse.covered < 0.6 * area && coarse.shortfall > 0.9 * 64.0 && coarse.shortfall < 1.1 * 64.0, "the grid's mesh: {coarse:?}");
    }

    #[test]
    fn a_feature_s_surface_missed_shows_in_its_patch_however_large_the_region() {
        // A region over 100 patches, 10 mm2 in each by the lines; the
        // facets cover all but one. Over the region that is 1% short, under
        // the lines' noise there; in its patch it is all of it.
        let grid = Grid::new([0.0; 3], 1.0);
        let side = SPLIT_BLOCK as f64;
        let patch = |k: i64| [k, 0, 0];
        let lines: Vec<([i64; 3], Patch)> = (0..100).map(|k| (patch(k), Patch { lines: 10.0, var: 10.0 * 0.0625, rim: 0.0 })).collect();
        // A square of 10 mm2 facing z in each patch but the 50th.
        let e = 10f64.sqrt();
        let square = |k: i64| {
            let x = k as f64 * side + 0.4;
            let c = [[x, 0.4, 0.5], [x + e, 0.4, 0.5], [x + e, 0.4 + e, 0.5], [x, 0.4 + e, 0.5]];
            [Facet::new([c[0], c[1], c[2]], None, 0.5).unwrap(), Facet::new([c[0], c[2], c[3]], None, 0.5).unwrap()]
        };
        let took: Vec<(Facet, f64)> = (0..100).filter(|&k| k != 50).flat_map(square).map(|f| (f, 1.0)).collect();
        let m = SurfaceMeasure::of("r", "load", &took, lines, &grid);
        assert!((m.covered - 990.0).abs() < 1e-9 && m.lines == 1000.0, "{m:?}");
        // Over the region: 10 short, within 3 x 0.79 x 10 = 23.7 of noise.
        assert!(1000.0 - 990.0 < m.uncertainty, "{m:?}");
        // Round its patch (it and the two beside it on the line): 10 short,
        // past 3 x 1.37 of noise.
        assert!((m.shortfall - (10.0 - 3.0 * (3.0 * 0.625f64).sqrt())).abs() < 1e-9, "{m:?}");
        // A seam: what the 50th patch's facets lack, the 51st's hold (the
        // mesh crossing their face off the surface). Nothing is missed.
        let mut seam: Vec<(Facet, f64)> = (0..100).filter(|&k| k != 50).flat_map(square).map(|f| (f, 1.0)).collect();
        seam.extend(square(51).map(|f| (f, 1.0)));
        let lines: Vec<([i64; 3], Patch)> = (0..100).map(|k| (patch(k), Patch { lines: 10.0, var: 10.0 * 0.0625, rim: 0.0 })).collect();
        let m = SurfaceMeasure::of("r", "load", &seam, lines, &grid);
        assert!(m.shortfall == 0.0, "{m:?}");
    }

    /// A grid cell split three times (`split_cells`).
    struct SheetCell([i64; 3]);
    impl Placed for SheetCell {
        fn idx(&self) -> [i64; 3] {
            self.0
        }
        fn level(&self) -> u8 {
            3
        }
        fn sub(&self) -> [u16; 3] {
            [0; 3]
        }
        fn mask(&self) -> u64 {
            ALL
        }
    }
}
