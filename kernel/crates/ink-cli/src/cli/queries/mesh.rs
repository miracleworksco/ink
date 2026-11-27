//! Planning, building, and placing meshes.

use super::budget::{Budget, SAMPLED_PER_CELL, box_area, ladder, view_floor, work};
use super::target::bounded_local;
use crate::{Arc, Args, CliView, HashMap, InkError, Instant, Item, Mesh, Rigid, Session, Tape, V3, json, max_extent, measured, mesh, n, paint, union_bounds};
use ink::lang::look::Look;
use ink::view::Local;

pub struct Meshed {
    pub label: String,
    /// The mesh in the body's own frame, shared by all instances of that body.
    pub mesh: Arc<Mesh>,
    pub colors: Arc<Vec<[f32; 3]>>,
    /// For seeing: the field's normal at each triangle's corners
    /// (`corner_normals`); else empty.
    pub normals: Arc<Vec<[f32; 3]>>,
    /// Where each vertex's paint (colour, metallic, roughness) is read: for
    /// seeing, under the surface (`ink::seen::Viewed::paint_at`); else the
    /// vertices.
    pub paint_at: Arc<Vec<V3>>,
    /// How it was made and how far it lies from the surface.
    pub quality: Arc<Quality>,
    /// Placement of this instance (None: the mesh is already in world space).
    pub place: Option<Rigid>,
    /// How much farther than in the body's frame the placed mesh may be
    /// from the instance's surface (ink::view::Local::off; 0 unplaced).
    pub off: f64,
}

pub struct Quality {
    /// The sampling, and the accuracy it was simplified to (0: uniform).
    pub res: f64,
    pub tol: f64,
    /// Triangles at res, before simplification.
    pub uniform_triangles: usize,
    /// The farthest a vertex at res moved from the result (<= tol).
    pub moved: f64,
    /// Measured for `ink mesh`: the sampled distance to the surface.
    pub deviation: Option<ink::simplify::Deviation>,
    /// Sign changes between samples, which the mesh leaves out.
    pub unresolved: mesh::Unresolved,
}

impl Meshed {
    /// The mesh in world coordinates.
    pub(super) fn world(&self) -> Mesh {
        match &self.place {
            None => (*self.mesh).clone(),
            Some(rg) => Mesh { verts: self.mesh.verts.iter().map(|&v| rg.apply(v)).collect(), tris: self.mesh.tris.clone() },
        }
    }
}

/// How `mesh_items` meshes: cells across a target without --res; the
/// accuracy as a fraction of res without --accuracy (0: uniform); whether to
/// measure the distance to the surface (`measure`).
pub struct MeshHow {
    pub div: f64,
    pub tol_per_res: f64,
    pub measure: bool,
    /// What the meshes may cost in all, in triangles sampled (`work`), each
    /// as often as it is made (every one `drawn` makes of it), whether made
    /// or found in the cache, and a flipbook's frames as predicted: meshes
    /// are made again only within it, so they are the model's alone, not
    /// its host's, threads', load's or cache's.
    pub budget: f64,
    /// For seeing, each body's own res by its surface's area (else by its
    /// box, Plan::own).
    pub by_area: bool,
    /// The finest res a mesh is made again at (0: any).
    pub finest: f64,
    /// For seeing, not measuring: a body is meshed once and placed by its
    /// placement's numbers even where they are not proven exact (a turn by
    /// a cosine of t), so a moving part keeps one mesh.
    pub seen: bool,
}

/// Mesh every target (`plan_meshes`, then `make_meshes`).
pub(super) fn mesh_items(s: &mut Session, a: &Args, items: &[Item], how: MeshHow) -> Result<(Vec<Meshed>, Option<Budget>), InkError> {
    let plan = plan_meshes(s, a, items, &how)?;
    make_meshes(s, a, plan, &how)
}

/// What `make_meshes` makes: each distinct mesh once (`jobs`: its item,
/// holding the field meshed, and its first target's place in the items),
/// its box, res and bend (`boxes`: triangles kept short along x, for a
/// strip a client bends round arcs as tight as it), the times it is made
/// (`times`: 1, or a client's frames of it, which the budget predicts), and
/// each target's mesh and placement (`targets`: its label, mesh, rigid
/// placement, and its instance's frame with its mesh's box). A client may
/// change what a mesh is (its field, box and res) before it is made.
pub struct Plan {
    pub jobs: Vec<(Item, usize)>,
    pub boxes: Vec<(V3, V3, f64, Option<f64>)>,
    pub times: Vec<f64>,
    pub targets: Vec<(String, usize, Option<Rigid>, Option<(Local, V3, V3)>)>,
    /// How many frames a client makes of the meshes it makes more than
    /// once (`times` over 1), from the triangles they sample and keep at
    /// one moment: what the budget predicts for them.
    pub frames: Option<Box<dyn Fn(usize, usize) -> usize>>,
    /// --res, and --max-res (infinite without).
    pub res: Option<f64>,
    pub coarsest: f64,
    div: f64,
    floor: f64,
    seen: bool,
}

impl Plan {
    /// The res a box of its own gives a mesh without --res.
    pub fn own(&self, lo: V3, hi: V3) -> f64 {
        let h = match self.seen {
            true => (box_area(lo, hi) / (6.0 * self.div * self.div)).sqrt(),
            false => (max_extent(lo, hi) / self.div).max(self.floor),
        };
        h.min(self.coarsest)
    }
}

/// What to mesh for the targets. Without --res, each object gets its own
/// resolution (its size / `how.div`, but no finer than the whole scene's
/// size / (8 * how.div), so small parts in huge assemblies stay
/// affordable; and no coarser than --max-res). Instances of one body share
/// a single mesh made in the body's frame (ModelView::local), at the res of
/// its size there, however each is turned. For seeing (`how.seen`) a
/// body's own res gives its box as many cells over its surface as a cube
/// `how.div` cells across (so a long part is not coarse across), and a
/// placement that reads t places the body at each motion sample too (the
/// GLB), so it must be rigid at each (ModelView::rigid).
pub fn plan_meshes(s: &mut Session, a: &Args, items: &[Item], how: &MeshHow) -> Result<Plan, InkError> {
    let res = a.cell("res")?;
    let coarsest = a.cell("max-res")?.unwrap_or(f64::INFINITY);
    let tol = a.f64("accuracy")?;
    if tol.is_some_and(|t| !(t >= 0.0 && t.is_finite())) {
        return Err(InkError::msg("--accuracy must be 0 (uniform) or positive and finite"));
    }
    let moves = |s: &Session, it: &Item| it.instance.as_ref().is_some_and(|(_, xf)| ink::claims::reads_t(s, &xf.r.iter().flatten().chain(&xf.t).copied().collect::<Vec<_>>()));
    let moments = match how.seen && items.iter().any(|it| moves(s, it)) {
        true => motion(s)?.1,
        false => vec![],
    };
    let locals: Vec<Option<(u32, Rigid, Option<Local>)>> = items
        .iter()
        .map(|it| {
            let local = bounded_local(s, it)?;
            Ok(match (s.rigid(it, &moments), local) {
                (Some(_), Some(l)) => Some((l.shape, l.at, Some(l))),
                (Some((body, rg)), None) if how.seen => {
                    let shape = s.doc.body(&body).unwrap().shape;
                    // If only the placed field has a bound, mesh it in world coordinates.
                    s.whole(shape).then_some((shape, rg, None))
                }
                _ => None,
            })
        })
        .collect::<Result<_, InkError>>()?;
    let bounds: Vec<Option<(V3, V3, f64)>> = items
        .iter()
        .zip(&locals)
        .map(|(it, l)| match l {
            Some((shape, rg, _)) if how.seen => s.bounds(*shape).map(|(lo, hi, h)| {
                let (a, b) = rg.aabb(lo, hi);
                (a, b, h)
            }),
            _ => s.item_bounds(it),
        })
        .collect();
    let all: Vec<(V3, V3, f64)> = bounds.iter().flatten().copied().collect();
    // Empty: proven so (a max of shapes apart, a subtraction of all of it).
    let empty = || {
        let names: Vec<&str> = items.iter().map(|it| it.label.as_str()).collect();
        InkError {
            hint: Some("a shape is empty where nothing is inside it: a max (intersection) of shapes that do not overlap, say; a union is their min".into()),
            ..InkError::msg(format!("nothing to mesh: every target is empty ({})", names.join(", ")))
        }
    };
    let Some((ulo, uhi)) = union_bounds(&all) else { return Err(empty()) };
    let floor = if how.seen { 0.0 } else { max_extent(ulo, uhi) / (8.0 * how.div) };
    let mut plan = Plan { jobs: Vec::new(), boxes: Vec::new(), times: Vec::new(), targets: Vec::new(), frames: None, res, coarsest, div: how.div, floor, seen: how.seen };
    let mut shared: HashMap<(String, u64), usize> = HashMap::new();
    for (k, ((it, b), l)) in items.iter().zip(&bounds).zip(&locals).enumerate() {
        let Some((lo, hi, _)) = *b else { continue };
        let h = res.unwrap_or(plan.own(lo, hi));
        let cells = max_extent(lo, hi) / h;
        if cells > 20000.0 {
            return Err(InkError::msg(format!("--res {h} gives {cells:.0} cells across '{}'; the limit is 20000", it.label)));
        }
        // Instanced: mesh the body once in its own frame.
        if let Some((shape, rg, local)) = *l {
            let body = it.instance.as_ref().unwrap().0.clone();
            let Some((llo, lhi, _)) = s.bounds(shape) else { continue };
            let h = res.unwrap_or(plan.own(llo, lhi));
            let j = *shared.entry((body.clone(), h.to_bits())).or_insert_with(|| {
                plan.jobs.push((Item { label: body.clone(), shape, attrs: it.attrs.clone(), instance: None, links: Vec::new() }, k));
                plan.boxes.push((llo, lhi, h, None));
                plan.jobs.len() - 1
            });
            plan.targets.push((it.label.clone(), j, Some(rg), local.map(|l| (l, llo, lhi))));
            continue;
        }
        plan.jobs.push((Item { label: it.label.clone(), shape: s.doc.placed_field(it), attrs: s.doc.placed_attrs(it), instance: None, links: it.links.clone() }, k));
        plan.boxes.push((lo, hi, h, None));
        plan.targets.push((it.label.clone(), plan.jobs.len() - 1, None, None));
    }
    if plan.targets.is_empty() {
        return Err(empty());
    }
    plan.times = vec![1.0; plan.jobs.len()];
    Ok(plan)
}

/// Make the meshes planned (`plan_meshes`). For seeing (`how.seen`, no
/// --res) a mesh is no finer than the model's floor times the cube root of
/// the times it is made (`view_floor`: its meshes within VIEW_TRIANGLES,
/// as their surfaces' areas predict them, ink::seen::surface), on the
/// ladder (`ladder`), and what it leaves out is searched for within its
/// share of the budget (`search`). A mesh that leaves material out is made
/// again under the thinnest of it, within `how.budget` (not one at
/// --max-res); what is left is reported (`mesh_warnings`). Each
/// mesh is reported as it is made (ink::claims::progress); the budget, as
/// planned and made, is returned for `cost_warnings`.
pub fn make_meshes(s: &mut Session, a: &Args, plan: Plan, how: &MeshHow) -> Result<(Vec<Meshed>, Option<Budget>), InkError> {
    let (res, coarsest, tol) = (plan.res, plan.coarsest, a.f64("accuracy")?);
    let Plan { jobs, mut boxes, times, targets, frames, div, .. } = plan;
    // For seeing, each mesh's surface (ink::seen::surface): its area
    // predicts what it samples (a body of thin strips over a large box
    // samples its strips', not its box's). Where none is found, its box's.
    use rayon::prelude::*;
    let sized = how.seen && res.is_none();
    let areas: Vec<f64> = match sized {
        true => {
            let fields: Vec<Arc<Tape>> = jobs.iter().map(|(it, _)| Arc::new(Tape::new(&s.doc.g, it.shape))).collect();
            fields.par_iter().zip(&boxes).map(|(t, b)| ink::seen::surface(t, &s.ctx, b.0, b.1)).collect()
        }
        false => vec![0.0; jobs.len()],
    };
    // For seeing, a mesh's own res gives its surface as many cells as a
    // cube `div` cells across has (Plan::own gives its box's, before its
    // surface is found): a belt looped round a large box is not coarse.
    if sized && how.by_area {
        for (b, &area) in boxes.iter_mut().zip(&areas).filter(|(_, a)| **a > 0.0) {
            b.2 = b.2.min((area / (6.0 * div * div)).sqrt());
        }
    }
    // Each mesh's triangles sampled at res 1 (all the times it is made), as
    // its area predicts them, and as its box does.
    let at_1: Vec<f64> = areas.iter().zip(&boxes).zip(&times).map(|((&area, b), t)| t * SAMPLED_PER_CELL * if area > 0.0 { area } else { box_area(b.0, b.1) }).collect();
    // For seeing, no mesh finer than the model's floor times the cube root
    // of the times it is made: its meshes' triangles sampled, as their
    // areas predict them, within VIEW_TRIANGLES.
    let parts: Vec<(f64, f64, f64)> = (0..jobs.len()).map(|j| (at_1[j], boxes[j].2, times[j])).collect();
    let floor = if sized { view_floor(&parts) } else { 0.0 };
    if sized {
        for (b, p) in boxes.iter_mut().zip(&parts) {
            b.2 = ladder(p.1.max(floor * p.2.cbrt())).min(coarsest);
        }
    }
    // Its deviation is sampled within 2 res of the mesh (mesh_tape).
    let plan: Vec<(String, usize, Option<Rigid>, f64)> =
        targets.into_iter().map(|(label, j, place, local)| (label, j, place, local.map_or(0.0, |(l, lo, hi)| l.off(lo, hi, 2.0 * boxes[j].2)))).collect();
    // Every mesh at once, each on its own thread (simplifying one is
    // serial, and most of the time).
    let tapes: Vec<Arc<Tape>> = jobs.iter().zip(&boxes).map(|((it, _), &(lo, hi, h, _))| tape_of(s, it, lo, hi, h)).collect();
    let paints: Vec<Vec<Arc<Tape>>> = jobs.iter().map(|(it, _)| paint_tapes(s, &it.attrs)).collect();
    let ctx = &s.ctx;
    let labels: Vec<&str> = jobs.iter().map(|(it, _)| it.label.as_str()).collect();
    // Whether a mesh's res is the plan's (no --res or --accuracy, and under
    // --max-res). A mesh's search for what it leaves out samples as much
    // again as its grid is planned to (mesh::mesh_checked; its own alone,
    // so it is the same mesh in any model); for seeing, grid and search no
    // more than the whole budget pays for (past that its grid is made
    // coarser, and says so).
    let free = |j: usize| sized && tol.is_none() && boxes[j].2 < coarsest;
    let limit = |j: usize| {
        let planned = at_1[j] / (boxes[j].2 * boxes[j].2);
        match how.seen {
            true => mesh::Limit { bricks: mesh::Limit::bricks_for(how.budget), search: mesh::Limit::samples_for(planned.min(how.budget)) },
            false => mesh::Limit { bricks: usize::MAX, search: mesh::Limit::samples_for(planned) },
        }
    };
    let seen_how = |j: usize, h: f64| ink::seen::How { res: h, limit: limit(j), tol, tol_per_res: how.tol_per_res, bend: boxes[j].3, export: !how.seen };
    // Each mesh (ink::seen::mesh, reported as made: `drawn`), what making
    // it cost (`work`, the same whether made or found in the cache), with
    // what `how` measures of it.
    let mesh = |j: usize, h: f64| {
        let (lo, hi, _, _) = boxes[j];
        let v = drawn(&tapes[j], &paints[j], ctx, lo, hi, seen_how(j, h), labels[j]);
        let cost = work(&v);
        let mut m = viewed(v);
        if how.measure {
            measure(&m.0, &mut m.2, &tapes[j], ctx, lo);
        }
        (m, cost)
    };
    // For seeing, each mesh is reported as it is made: it is cached then.
    let mut total = REPORTED.load(std::sync::atomic::Ordering::Relaxed) + jobs.len();
    if how.seen {
        ink::claims::progress(json!({"total": total, "of": "meshes", "kept": true}));
    }
    let first: Vec<_> = (0..jobs.len()).into_par_iter().map(|j| mesh(j, boxes[j].2)).collect();
    // What the meshes cost (`work`), made or found in the cache alike
    // (`how.budget`), and each kept's own. A client's frames of a mesh
    // (Plan::frames) are made after: they are predicted, not counted.
    let mut spent: f64 = first.iter().map(|m| m.1).sum();
    let mut own: Vec<f64> = first.iter().map(|m| m.1).collect();
    let mut made: Vec<(Mesh, Vec<[f32; 3]>, Quality, Vec<V3>)> = first.into_iter().map(|m| m.0).collect();
    // A mesh that leaves material out (mesh::Unresolved: a part, or a
    // stretch of one, thinner than res) is made again at the ladder's res
    // under the thinnest of it (no finer than `how.finest`, nor than the
    // cells' limit across its box allows, as --res). Those
    // drawing the most of what is left out for what they cost first (each
    // cell left out holds a sheet across it, two faces of res squared),
    // while what the meshes cost, a client's frames as predicted and what
    // the next is predicted to cost stay within `how.budget`; again while
    // one still leaves material out. A gap
    // thinner than res (two surfaces closed into one, a ball resting in its
    // race) is not: it looks as it would apart, and where it thins to
    // nothing no res resolves it. What is left is reported.
    if how.seen && res.is_none() && tol.is_none() {
        loop {
            // A mesh a client makes again per frame (Plan::frames), after
            // this: each frame predicted at its mesh's.
            let live: Vec<usize> = (0..jobs.len()).filter(|&j| times[j] > 1.0).collect();
            let (sampled, kept) = live.iter().fold((0, 0), |(a, b), &j| (a + made[j].2.uniform_triangles, b + made[j].0.tris.len()));
            let frames = match &frames {
                Some(f) if !live.is_empty() => f(sampled, kept) as f64 * live.iter().map(|&j| own[j]).sum::<f64>(),
                _ => 0.0,
            };
            let mut room = how.budget - spent - frames;
            // Each to make again, its res, what it costs (its cost now per
            // cell, over the cells at the finer res, and at least what its
            // surface's area predicts there: a part left out whole drew
            // none) and the area it would draw that it leaves out.
            let mut again: Vec<(usize, f64, f64, f64)> = (0..jobs.len())
                .filter_map(|j| {
                    let q = &made[j].2;
                    let u = &q.unresolved;
                    let h = ladder(u.thinnest.min(q.res) / 2f64.powf(0.25)).max(how.finest).max(ladder(max_extent(boxes[j].0, boxes[j].1) / 20000.0));
                    let drawn = 2.0 * u.material as f64 * q.res * q.res;
                    (free(j) && u.thinnest > 0.0 && u.refused == 0.0 && h < q.res).then(|| (j, h, (times[j] * own[j] * (q.res / h).powi(2)).max(at_1[j] / (h * h)), drawn))
                })
                .collect();
            again.sort_by(|a, b| (a.2 / a.3).total_cmp(&(b.2 / b.3)));
            again.retain(|&(_, _, more, _)| {
                more <= room && {
                    room -= more;
                    true
                }
            });
            if again.is_empty() {
                break;
            }
            total += again.len();
            ink::claims::progress(json!({"total": total, "of": "meshes", "kept": true}));
            let finer: Vec<_> = again.par_iter().map(|&(j, h, _, _)| (j, mesh(j, h))).collect();
            for (j, (m, c)) in finer {
                spent += c;
                (made[j], own[j]) = (m, c);
            }
        }
    }
    // For seeing, what each mesh was planned at and made at (a client's
    // frames counted as its first until it makes them).
    let budget = sized.then(|| Budget {
        jobs: plan.iter().map(|p| (p.0.clone(), p.1)).collect(),
        names: (0..jobs.len()).map(|j| plan.iter().find(|p| p.1 == j).map_or(jobs[j].0.label.clone(), |p| p.0.clone())).collect(),
        made: made.iter().zip(&times).map(|(m, t)| (m.2.res, t * m.2.uniform_triangles as f64)).collect(),
        parts,
        floor,
        boxes: boxes.clone(),
        spent,
    });
    let made: Vec<(Arc<Mesh>, Arc<Vec<[f32; 3]>>, Arc<Vec<[f32; 3]>>, Arc<Vec<V3>>, Arc<Quality>)> = made
        .into_iter()
        .zip(&jobs)
        .map(|((m, normals, q, at), (it, k))| {
            let colors = s.colors(it, &at, *k);
            (Arc::new(m), Arc::new(colors), Arc::new(normals), Arc::new(at), Arc::new(q))
        })
        .collect();
    // A target whose body gives no colour takes the palette's by its place
    // among the targets: two sharing a mesh (two instances of a screw) are
    // told apart as two meshes are.
    let palette = |k: usize| ink::lang::look::PALETTE[k % ink::lang::look::PALETTE.len()];
    let unpainted: Vec<bool> = jobs.iter().map(|(it, _)| ink::lang::look::Look::of(&it.attrs).color.is_none()).collect();
    Ok((
        plan.into_iter()
            .enumerate()
            .map(|(k, (label, j, place, off))| {
                let (mesh, mut colors, normals, paint_at, quality) = made[j].clone();
                if unpainted[j] {
                    colors = Arc::new(vec![palette(k); colors.len()]);
                }
                Meshed { label, mesh, colors, normals, paint_at, quality, place, off }
            })
            .collect(),
        budget,
    ))
}

/// A view mesh of a body (ink::seen::mesh), reported as made or found in
/// the cache (`label`).
pub fn drawn(tape: &Arc<Tape>, paint: &[Arc<Tape>], ctx: &ink::eval::Ctx, lo: V3, hi: V3, how: ink::seen::How, label: &str) -> ink::seen::Viewed {
    let t = Instant::now();
    let v = ink::seen::mesh(tape, paint, ctx, lo, hi, how);
    reported(&format!("{label} at res {}", measured(how.res)), v.cached, t);
    v
}

/// A viewing mesh as `mesh_tape` gives one.
fn viewed(v: ink::seen::Viewed) -> (Mesh, Vec<[f32; 3]>, Quality, Vec<V3>) {
    let q = Quality { res: v.res, tol: v.tol, uniform_triangles: v.uniform_triangles, moved: v.moved, deviation: None, unresolved: v.unresolved };
    (Mesh { verts: v.verts, tris: v.tris }, v.normals, q, v.paint_at)
}

/// The paint values of an item's look (its colour's channels, metallic,
/// roughness, transmission) that vary over space.
pub fn paint_tapes(s: &Session, attrs: &ink::lang::Attrs) -> Vec<Arc<Tape>> {
    let fs = Look::of(attrs).channels().iter().flatten().copied().collect::<Vec<_>>();
    let spatial = |f: u32| s.doc.g.depends_on(f, |v| matches!(v, ink::graph::Var::X | ink::graph::Var::Y | ink::graph::Var::Z));
    fs.into_iter().filter(|&f| spatial(f)).map(|f| Arc::new(Tape::new(&s.doc.g, f))).collect()
}

/// Metallic, roughness and transmission at each vertex (0..1; the look's
/// defaults where it gives none).
pub fn looks(s: &Session, attrs: &ink::lang::Attrs, verts: &[V3]) -> Vec<[f32; 3]> {
    use ink::lang::look::{METALLIC, ROUGHNESS, TRANSMISSION};
    let look = Look::of(attrs);
    paint(s, [look.metallic, look.roughness, look.transmission], [METALLIC, ROUGHNESS, TRANSMISSION], verts)
}

/// An item's field clipped to its box and a cell beyond.
pub fn tape_of(s: &mut Session, it: &Item, lo: V3, hi: V3, h: f64) -> Arc<Tape> {
    let f = s.doc.clip(it.shape, [lo[0] - h, lo[1] - h, lo[2] - h], [hi[0] + h, hi[1] + h, hi[2] + h]);
    Arc::new(Tape::new(&s.doc.g, f))
}

/// Its sampled distance from the surface (ink::simplify::deviation,
/// sampled within 2 res of it, every quarter of its tolerance), for `ink
/// mesh`.
fn measure(m: &Mesh, q: &mut Quality, tape: &Arc<Tape>, ctx: &ink::eval::Ctx, lo: V3) {
    let h = q.res;
    // The field near the surface, a tape per block of 8 cells.
    let local = ink::octree::Local::new(tape, ctx, lo, 8.0 * h);
    let step = if q.tol > 0.0 { q.tol } else { h / 30.0 } / 4.0;
    q.deviation = Some(ink::simplify::deviation(&|p: &[V3]| local.eval(p), &|p: &[V3]| local.grads(p), m, h, step, 2.0 * h));
}

/// Meshes reported done in this process (`reported`), which a later total
/// counts on from (a GLB's flipbook frames after its meshes).
pub static REPORTED: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Report a mesh done (ink::claims::progress): made, or found in the
/// cache, since `t`.
pub fn reported(what: &str, cached: bool, t: Instant) {
    REPORTED.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let status = if cached { "cached" } else { "made" };
    ink::claims::progress(json!({"done": what, "status": status, "ms": n(t.elapsed().as_secs_f64() * 1e3)}));
}

/// Motion samples per second of the span, with their limits.
const MOTION_RATE: f64 = 30.0;
const MOTION_SAMPLES: (usize, usize) = (60, 1200);

/// Moments over the span: `rate` a second, between lo and hi of them but
/// at most n (and at least 1), from its start to its end.
pub fn span_times(span: (f64, f64), n: usize, (lo, hi): (usize, usize), rate: f64) -> Vec<f64> {
    let len = span.1 - span.0;
    let n = ((rate * len).ceil() as usize).clamp(lo, hi).min(n).max(1);
    (0..=n).map(|i| span.0 + len * i as f64 / n as f64).collect()
}

/// The moments a GLB samples motion at, and the context at each.
pub fn motion(s: &Session) -> Result<(Vec<f64>, Vec<ink::eval::Ctx>), InkError> {
    let times = span_times(s.doc.t_span.unwrap_or((0.0, 1.0)), usize::MAX, MOTION_SAMPLES, MOTION_RATE);
    let ctxs = times.iter().map(|&t| s.doc.ctx(&s.overrides, t).map_err(|e| InkError::msg(format!("at t = {t}: {e}")))).collect::<Result<_, _>>()?;
    Ok((times, ctxs))
}
