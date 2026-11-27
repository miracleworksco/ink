//! Work and resolution budgets for view meshes.

use super::mesh::MeshHow;
use crate::{HashMap, V3, max_extent, measured};

/// Meshes for seeing as planned and made (make_meshes), for
/// `cost_warnings`: each target's mesh, each mesh's first target, its res
/// and the triangles it sampled (all its meshes: a client's frames counted
/// as its first until it makes them and says what they sampled),
/// `view_floor`'s parts and floor, each mesh's box and planned res, and
/// what the meshes made cost the budget (`work`).
pub struct Budget {
    pub jobs: HashMap<String, usize>,
    pub names: Vec<String>,
    pub made: Vec<(f64, f64)>,
    pub parts: Vec<(f64, f64, f64)>,
    pub floor: f64,
    pub boxes: Vec<(V3, V3, f64, Option<f64>)>,
    pub spent: f64,
}

impl Budget {
    /// What the meshes cost (`cost_warnings`).
    pub fn warnings(&self) -> Vec<String> {
        cost_warnings(&self.parts, self.floor, &self.boxes, &self.made, &self.names.iter().map(String::as_str).collect::<Vec<_>>())
    }
}

/// GLB's meshes, for seeing at a screen's resolution: coarser than for
/// printing, and simplified to a quarter of res, within VIEW_WORK (in
/// triangles, as `work` counts them).
pub const SEEN: MeshHow = MeshHow { div: 128.0, tol_per_res: 0.25, measure: false, budget: VIEW_WORK / TRIANGLE_WORK, by_area: true, finest: 0.0, seen: true };

/// The machine instructions a model's view meshes may retire (`work`):
/// 30 to 70 CPU s of a Zen 5 core, which retires 0.7 to 1.5e10 of them a
/// second meshing (more where its field is most of the work).
const VIEW_WORK: f64 = 4.4e11;

/// The triangles a model's GLB meshes are planned to sample when first
/// made, as their surfaces' areas predict them (`view_floor`): the rest of
/// SEEN's budget is for meshes made again.
pub(super) const VIEW_TRIANGLES: f64 = 2.5e6;

/// The machine instructions a view mesh retires a triangle it samples,
/// besides reading its field and paint (`work`): marching its cells,
/// simplifying and shading it.
pub(super) const TRIANGLE_WORK: f64 = 8.0e4;

/// The machine instructions reading a field or paint retires for each
/// instruction Tape::work counts of it: the reading's own bookkeeping, and
/// the descent's passes and tapes it does not count.
pub(super) const READ_WORK: f64 = 1.8;

/// The triangles a surface samples per cell of it, as predicted: a flat
/// face's two (a folded surface samples more).
pub(super) const SAMPLED_PER_CELL: f64 = 2.0;

pub(super) fn box_area(lo: V3, hi: V3) -> f64 {
    let d = ink::v3::sub(hi, lo);
    2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0])
}

/// The model's floor for its view meshes: 0 where each at its own res
/// samples at most VIEW_TRIANGLES, else the least `f` at which they do,
/// each mesh at f times the cube root of the times it is made (a
/// flipbook's frames) where that is coarser than its own res. That is
/// the least error over what is seen for the triangles: res h over area A
/// is off by about h over A, each mesh is seen once at a time, and it
/// samples 2 A / h^2 each time it is made, so the sum of A h within the
/// budget is least where every h^3 is in proportion to the times made.
/// So small detailed parts are coarsened first, large ones are not made
/// finer, and a mesh made 40 times counts as about 3.4 of its size, not
/// 40. `parts`: each mesh's triangles predicted at res 1 (all the times it
/// is made), its own res, and the times it is made.
pub(super) fn view_floor(parts: &[(f64, f64, f64)]) -> f64 {
    let sampled = |f: f64| parts.iter().map(|&(c, h, t)| c / h.max(f * t.cbrt()).powi(2)).sum::<f64>();
    if sampled(0.0) <= VIEW_TRIANGLES {
        return 0.0;
    }
    // At hi every mesh's res is hi or coarser: within VIEW_TRIANGLES.
    let mut lo = parts.iter().map(|p| p.1 / p.2.cbrt()).fold(f64::INFINITY, f64::min);
    let mut hi = (parts.iter().map(|p| p.0).sum::<f64>() / VIEW_TRIANGLES).sqrt();
    for _ in 0..60 {
        let mid = (lo * hi).sqrt();
        if sampled(mid) <= VIEW_TRIANGLES { hi = mid } else { lo = mid }
    }
    hi
}

/// Warnings for each mesh that sampled more than half of the triangles
/// the meshes sampled, where they sampled more than VIEW_TRIANGLES, or
/// whose cost alone brought another mesh's res up a step of the ladder
/// (only a mesh predicted at over 1 - 2^(-1/2) of them can: without it the
/// floor falls at most as the square root of what is left): it, its
/// frames, res and size, and what it sampled. `parts` as `view_floor`
/// takes them, `boxes` each mesh's box and planned res, `made` each mesh's
/// res and the triangles it sampled (its frames' as many as its first's),
/// `names` each mesh's first target. And, where they sampled more than
/// VIEW_TRIANGLES, how many were made finer than planned.
pub(super) fn cost_warnings(parts: &[(f64, f64, f64)], floor: f64, boxes: &[(V3, V3, f64, Option<f64>)], made: &[(f64, f64)], names: &[&str]) -> Vec<String> {
    let planned: Vec<f64> = parts.iter().zip(boxes).map(|(p, b)| p.0 / (b.2 * b.2)).collect();
    let (planned_all, all) = (planned.iter().sum::<f64>(), made.iter().map(|m| m.1).sum::<f64>());
    let mega = |t: f64| format!("{:.2} M", t / 1e6);
    // A mesh's res at floor f.
    let at = |k: usize, f: f64| ladder(parts[k].1.max(f * parts[k].2.cbrt()));
    (0..parts.len())
        .filter_map(|j| {
            let most = all > VIEW_TRIANGLES && made[j].1 > 0.5 * all;
            let without = (floor > 0.0 && planned[j] > (1.0 - 0.5f64.sqrt()) * planned_all)
                .then(|| view_floor(&parts.iter().enumerate().filter(|&(k, _)| k != j).map(|(_, p)| *p).collect::<Vec<_>>()))
                .filter(|&w| (0..parts.len()).any(|k| k != j && at(k, w) < at(k, floor)));
            (most || without.is_some()).then(|| {
                let what = match parts[j].2 {
                    1.0 => "its mesh".to_string(),
                    t => format!("its {} frames", t - 1.0),
                };
                let mut w = format!(
                    "'{}': {what} at res {} over {} mm sampled {} of the {} triangles the GLB's meshes sampled ({:.0}%)",
                    names[j],
                    measured(made[j].0),
                    measured(max_extent(boxes[j].0, boxes[j].1)),
                    mega(made[j].1),
                    mega(all),
                    100.0 * made[j].1 / all
                );
                if let Some(w0) = without {
                    let others = if w0 == 0.0 { "each at its own res".to_string() } else { format!("res {}", measured(ladder(w0))) };
                    w += &format!("; it brings the other meshes up to res {} (without it, {others})", measured(ladder(floor)));
                }
                w
            })
        })
        .chain(
            // Past the budget: what drawing them took.
            {
                let finer: Vec<f64> = made.iter().zip(boxes).filter(|(m, b)| m.0 < b.2).map(|(m, _)| m.0).collect();
                (all > VIEW_TRIANGLES && !finer.is_empty()).then(|| {
                    format!(
                        "the GLB's meshes sampled {} triangles, {:.1} times the {} they are planned within: made finer than planned (to res {} \
                         at the finest), {} of its {} meshes drew what the floor (res {}) brought past their features",
                        mega(all),
                        all / VIEW_TRIANGLES,
                        mega(VIEW_TRIANGLES),
                        measured(finer.iter().fold(f64::INFINITY, |a, &r| a.min(r))),
                        finer.len(),
                        made.len(),
                        measured(ladder(floor))
                    )
                })
            },
        )
        .collect()
}

/// `h` rounded up to the ladder of view resolutions, 2^(k/4) mm: a part
/// gets the same res, so the same cached mesh, in models whose floors
/// differ a little.
pub(super) fn ladder(h: f64) -> f64 {
    const STEP: [f64; 4] = [1.0, 1.189207115002721, std::f64::consts::SQRT_2, 1.681792830507429];
    let k = (4.0 * h.log2() - 1e-9).ceil() as i32;
    2f64.powi(k.div_euclid(4)) * STEP[k.rem_euclid(4) as usize]
}

/// What making a view mesh costs, in triangles sampled: TRIANGLE_WORK
/// machine instructions each, its triangles' own and its reading its field
/// and paint (READ_WORK an instruction Tape::work counts): its grid's,
/// crossings' and search's samples and its proofs over boxes
/// (mesh::Unresolved::sampled and read), splitting, simplifying and
/// shading it (`field_work`), and its paint (`paint_work`). Counted, each
/// ratio holds to within WORK_TOLERANCE on one core
/// (the_view_budget_charges_what_meshes_count).
pub(super) fn work(v: &ink::seen::Viewed) -> f64 {
    let u = &v.unresolved;
    let read = (u.sampled + u.read + v.field_work + v.paint_work) as f64;
    v.uniform_triangles as f64 + READ_WORK * read / TRIANGLE_WORK
}
