//! Load-time calls into kernel services: bounds, grid baking,
//! redistancing, STL import.

use super::*;
use ink_core::interval::{self, I};

/// What a closed mesh's samples are (import-stl).
const SIGNED: &str = "the distance to the mesh's triangles, negative inside the solid they enclose (a point is inside when a line from it \
    to far away crosses the mesh an odd number of times: whichever way its facets face; where it crosses itself, overlaps cancel): the solid's \
    signed distance unless triangles overlap in one plane (as where touching bodies share a face: parity cancels them, the distance still counts \
    them); sampled in floating point, with no error bound: measured_error is the deviation measured";

/// A sampled grid as cached: its points, values and what is known of its
/// source.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
struct Snap {
    origin: V3,
    h: f64,
    dims: [usize; 3],
    values: Vec<f32>,
    band: Option<(f64, Option<f64>)>,
    measured: Option<(f64, usize)>,
    note: Option<String>,
    /// The box a region import is cut to (sampled_grid).
    cut: Option<(V3, V3)>,
    /// A mesh that bounds no solid (said at load).
    open: bool,
}

/// The grids this process has sampled or read, by key (Doc::grid_key), up
/// to KEPT_BYTES of values: a load again (at a case) or a form asking for
/// one again (a letter met twice) takes it here, the disk cache on or off.
static KEPT: std::sync::Mutex<(usize, Vec<(u128, std::sync::Arc<Snap>)>)> = std::sync::Mutex::new((0, Vec::new()));
const KEPT_BYTES: usize = 1 << 29;

impl Snap {
    fn kept(key: u128) -> Option<Snap> {
        KEPT.lock().unwrap().1.iter().find(|k| k.0 == key).map(|k| (*k.1).clone())
    }

    fn keep(key: u128, s: &Snap) {
        let mut kept = KEPT.lock().unwrap();
        let bytes = 4 * s.values.len();
        if kept.0 + bytes <= KEPT_BYTES && !kept.1.iter().any(|k| k.0 == key) {
            kept.0 += bytes;
            kept.1.push((key, std::sync::Arc::new(s.clone())));
        }
    }
}

impl Doc {
    /// The names (tag ids) body `body` (its shape `shape`) is built from,
    /// each of which should form some of its surface (a name
    /// feature::surface_names proves forms none is hidden inside another
    /// shape or cut away). Another body used in this one (its
    /// `shape-of`) counts as one construction, under its own name: the
    /// names inside it are its own, checked on it, and a cut by a whole
    /// part need not touch every feature of that part.
    pub fn named_parts(&self, body: &str, shape: NodeId) -> Vec<u32> {
        use crate::graph::Node;
        let g = &self.g;
        let ids: HashMap<&str, u32> = g.tag_names.iter().enumerate().map(|(t, n)| (n.as_str(), t as u32)).collect();
        let own = ids.get(body).copied();
        let bodies: HashSet<u32> = self.bodies.iter().filter_map(|b| ids.get(b.name.as_str()).copied()).collect();
        let mut named = Vec::new();
        let mut seen = HashSet::new();
        let mut stack = vec![shape];
        while let Some(n) = stack.pop() {
            if !seen.insert(n) {
                continue;
            }
            let node = &g.nodes[n as usize];
            if let Node::Tag(t, _) = *node
                && Some(t) != own
            {
                named.push(t);
                if bodies.contains(&t) {
                    continue;
                }
            }
            let (ch, k) = node.children();
            stack.extend_from_slice(&ch[..k]);
        }
        named.sort_unstable();
        named.dedup();
        named
    }

    /// `f` clipped to its reach, which must hold all of it: a service
    /// that needs the whole shape (its box, a grid over it) errors on one
    /// that is not proven bounded.
    fn clip_bounded(&mut self, e: &Sexp, name: &str, f: NodeId, ctx: &Ctx) -> R<(NodeId, Reach)> {
        let (c, reach) = self.clip_reach(f, ctx);
        if !reach.whole {
            return Err(self.err_hint(
                e,
                format!("{name}: the shape is unbounded or not proven bounded (a half-space, an unclipped lattice, or a field interval bounds cannot bound at infinity)"),
                Some("intersect it with a containing shape first, e.g. (max (box 40 40 40) f)".into()),
            ));
        }
        Ok((c, reach))
    }

    /// The cache key of a sampled grid: its content, the source (a field as a function, or an STL file's bytes)
    /// and the resolution. A field's fingerprint folds in every value it
    /// reads, a package's result (`study-field`) included, so a changed
    /// solver changes the key through the values it gives: the key need
    /// not name the package.
    /// `inner`: the sampled data the field rests on (its own bound does not
    /// reach that), which the field's fingerprint need not show (a param a
    /// solve over samples sets).
    pub(super) fn grid_key(&self, name: &str, res: Option<f64>, region: Option<(V3, V3)>, field: Option<NodeId>, stl: &[u8], params: &[f64], inner: &[String]) -> u128 {
        let mut fp = crate::graph::Fingerprint::default();
        // Each part counted or tagged, so no two differ only in where one
        // ends and the next begins.
        fp.str(name).f64(res.unwrap_or(-1.0)).word(inner.len() as u64);
        inner.iter().for_each(|w| _ = fp.str(w));
        fp.word(region.is_some() as u64);
        if let Some((lo, hi)) = region {
            lo.iter().chain(&hi).for_each(|c| _ = fp.f64(*c));
        }
        match field {
            Some(f) => self.g.fingerprint(&[f], params, false, fp.word(1)),
            None => {
                fp.word(0).word(stl.len() as u64).raw(stl);
            }
        }
        crate::cache::key(fp.value())
    }

    /// `(grid-bake f)`, `(redistance f)`, `(import-stl "file.stl")`: a
    /// sampled field, cached on disk, saying what it was sampled from and,
    /// where that can be bounded soundly, how far it may be from it
    /// (grid::Source; docs/IR.md, Load-time programming, Sampled fields).
    fn sampled_grid(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        use rayon::prelude::*;
        let usage = "usage: (grid-bake field [(res mm)]), (redistance field [(res mm)]) or (import-stl \"file.stl\" [(res mm)] [(region (x0 y0 z0) (x1 y1 z1))])";
        let Some((src, opts)) = args.split_first() else { return Err(self.err_hint(e, "missing argument", Some(usage.into()))) };
        let (mut res, mut region) = (None, None);
        for o in opts {
            match o.list() {
                Some([h, v]) if h.sym() == Some("res") => res = Some(self.cell_size(v, env)?),
                Some([h, a, b]) if h.sym() == Some("region") && name == "import-stl" => {
                    let (lo, hi) = (self.const_v3(a)?, self.const_v3(b)?);
                    if (0..3).any(|k| !(lo[k] < hi[k])) {
                        return Err(self.err_hint(o, "region: its first corner must be below its second on every axis", Some(usage.into())));
                    }
                    region = Some((lo, hi));
                }
                _ => return Err(self.err_hint(o, "bad option", Some(usage.into()))),
            }
        }
        // The grid's step over a box of largest extent ext: (res ...), else
        // ext / 128, and never under ext 1e-4.
        let step = |ext: f64| res.unwrap_or(ext / 128.0).max(ext * 1e-4);
        let at_form = format!("{}:{}", self.file, e.line);
        // The source: a field, or the bytes of an STL file.
        // The sampled data a field reads (directly, or through a solve
        // setting a param it reads), which its own bound does not reach.
        let mut inner: Vec<String> = Vec::new();
        let mut from: Vec<std::sync::Arc<crate::grid::Source>> = Vec::new();
        let (field, stl, what) = if name == "import-stl" {
            let SK::Str(rel) = &src.k else { return Err(self.err(src, "import-stl path must be a string")) };
            let path = self.dir().join(rel);
            let bytes =
                std::fs::read(&path).inspect_err(|_| self.files_read.push((path.clone(), None))).map_err(|er| self.err(src, format!("cannot read {}: {er}", shown(&path))))?;
            self.read(&shown(&path), Some(&path), &bytes);
            let what = match region {
                Some((lo, hi)) => format!("{rel} in the region ({} {} {}) ({} {} {})", lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]),
                None => rel.clone(),
            };
            (None, bytes, what)
        } else {
            let f = self.field(src, env)?;
            self.bake(&[f]);
            from = self.inexact(&[f], &[]).samples;
            inner = from.iter().map(|s| s.what.clone()).collect();
            let mut text = src.to_string();
            if text.chars().count() > 80 {
                text = text.chars().take(77).chain("...".chars()).collect();
            }
            (Some(f), Vec::new(), text)
        };
        let ctx = self.load_ctx_for(&field.into_iter().collect::<Vec<_>>());
        let kind = match name {
            "import-stl" => "mesh",
            "grid-bake" => "bake",
            _ => "redistance",
        };
        // Sampled grids are cached on disk, so every command that loads the
        // file does not redo them; with what is known of their source.
        let key = self.grid_key(name, res, region, field, &stl, &ctx.params, &inner);
        let finish = |d: &mut Doc, c: Snap| -> Value {
            // The res asked, under the grid's floor (step), was coarsened.
            if let Some(r) = res.filter(|&r| c.h > r) {
                let smaller = if name == "import-stl" { "sample a smaller (region ...)" } else { "sample a smaller field" };
                let w = format!("{name} '{what}' ({at_form}): (res {r}) is under 1e-4 of the grid's largest extent, so it is sampled at {}; {smaller} to sample finer", c.h);
                if !d.load_warnings.contains(&w) {
                    d.load_warnings.push(w);
                }
            }
            if c.open {
                let w = format!("import-stl '{what}' ({at_form}): {}", c.note.as_deref().unwrap_or_default());
                if !d.load_warnings.contains(&w) {
                    d.load_warnings.push(w);
                }
            }
            let source = crate::grid::Source {
                what: what.clone(),
                at: at_form.clone(),
                kind,
                res: c.h,
                band: c.band.map(|(err, outside)| crate::grid::Band { err, outside }),
                measured: c.measured,
                note: c.note,
                unit: "",
                span: None,
                component: None,
                // (a refusal of what it was sampled from is its own)
                undecided: crate::inexact::undecided_from(&from),
            };
            let data = std::sync::Arc::new(crate::grid::GridData::new(c.origin, c.h, c.dims, c.values).with_source(source));
            // what it was sampled from, as `(bounds f)` keeps it
            if let Some(s) = data.source().filter(|_| !from.is_empty()) {
                d.found_from.push((s.clone(), from.clone()));
            }
            let at = [d.g.var(Var::X), d.g.var(Var::Y), d.g.var(Var::Z)];
            let mut f = d.g.grid(data, at);
            // A region import is the mesh's solid inside the box it was cut
            // to: positive beyond it (the box's own field, 1-Lipschitz).
            if let Some((lo, hi)) = c.cut {
                for k in 0..3 {
                    let (l, h) = (d.g.constant(lo[k]), d.g.constant(hi[k]));
                    let below = d.g.sub(l, at[k]);
                    let above = d.g.sub(at[k], h);
                    f = d.g.max(f, below);
                    f = d.g.max(f, above);
                }
            }
            Value::Field(f, None)
        };
        if let Some(c) = Snap::kept(key).or_else(|| crate::cache::load::<Snap>("grids", key)) {
            crate::cache::served();
            Snap::keep(key, &c);
            return Ok(finish(self, c));
        }
        // Sampling a field needs it bounded (clip_bounded), and a grid of
        // sane size: checked before any meshing or sampling.
        let sized = |d: &Doc, lo: V3, hi: V3| -> R<()> {
            let ext = crate::query::max_extent(lo, hi);
            let h = step(ext);
            let pts: f64 = (0..3).map(|a| (hi[a] - lo[a]) / h + 7.0).product();
            if pts > 64e6 {
                return Err(d.err(
                    e,
                    format!("{name}: a grid of {pts:.2e} points is too large. Increase (res ...){}", if name == "import-stl" { ", or sample a (region ...)" } else { "" }),
                ));
            }
            Ok(())
        };
        // What to sample, and where.
        enum Src {
            Field(std::sync::Arc<Tape>),
            Mesh(crate::grid::MeshDistance),
            /// Distance to a mesh of the field's surface, signed by the field.
            Surface(crate::grid::MeshDistance, std::sync::Arc<Tape>),
        }
        let (source, lo, hi) = if name == "redistance" {
            let f = field.unwrap();
            let clipped = self.clip_bounded(e, name, f, &ctx)?;
            let tape = std::sync::Arc::new(Tape::new(&self.g, clipped.0));
            let (lo, hi, _) = crate::query::bounds(&tape, &ctx, clipped.1.lo, clipped.1.hi).ok_or_else(|| self.err(e, "redistance of an empty field"))?;
            sized(self, lo, hi)?;
            let ext = crate::query::max_extent(lo, hi);
            // Mesh twice as fine as the grid, so the surface is
            // captured well below the grid's own interpolation error.
            let hm = 0.5 * step(ext);
            let m = crate::mesh::mesh_in_time(&tape, &ctx, lo, hi, hm).ok_or_else(|| self.stopped(e))?;
            if m.tris.is_empty() {
                return Err(self.err(e, "redistance: the field has no surface"));
            }
            let tris: Vec<[V3; 3]> = m.tris.iter().map(|t| t.map(|i| m.verts[i as usize])).collect();
            (Src::Surface(crate::grid::MeshDistance::in_time(&tris).ok_or_else(|| self.stopped(e))?, tape), lo, hi)
        } else if name == "grid-bake" {
            let f = field.unwrap();
            let clipped = self.clip_bounded(e, name, f, &ctx)?;
            let tape = std::sync::Arc::new(Tape::new(&self.g, clipped.0));
            let (lo, hi, _) = crate::query::bounds(&tape, &ctx, clipped.1.lo, clipped.1.hi).ok_or_else(|| self.err(e, "grid-bake of an empty field"))?;
            sized(self, lo, hi)?;
            (Src::Field(std::sync::Arc::new(Tape::new(&self.g, f))), lo, hi)
        } else {
            let tris = crate::grid::read_stl(&stl).map_err(|m| self.err(src, m))?;
            if tris.is_empty() {
                return Err(self.err(src, "STL has no triangles"));
            }
            let md = crate::grid::MeshDistance::in_time(&tris).ok_or_else(|| self.stopped(e))?;
            if md.triangles() == 0 {
                return Err(self.err(src, format!("STL has no triangles of nonzero area ({} degenerate facet{})", tris.len(), if tris.len() == 1 { "" } else { "s" })));
            }
            let (mut lo, mut hi) = md.bounds();
            // A region samples only the part of the mesh's box inside it.
            if let Some((rlo, rhi)) = region {
                lo = std::array::from_fn(|k| lo[k].max(rlo[k]));
                hi = std::array::from_fn(|k| hi[k].min(rhi[k]));
                if (0..3).any(|k| lo[k] > hi[k]) {
                    return Err(self.err(e, "import-stl: the region holds none of the mesh's box"));
                }
            }
            sized(self, lo, hi)?;
            (Src::Mesh(md), lo, hi)
        };
        let ext = crate::query::max_extent(lo, hi);
        let h = step(ext);
        let pad = 3.0 * h;
        // A mesh's grid is offset by an irrational-looking fraction of a
        // step, so its lines miss the mesh's edges and vertices at round
        // coordinates (where its side is unclear: MeshDistance::sample).
        let shift = if matches!(source, Src::Mesh(_)) { 0.381_966_011_250_105_1 * h } else { 0.0 };
        let origin: V3 = std::array::from_fn(|a| lo[a] - pad - shift);
        let dims: [usize; 3] = std::array::from_fn(|a| (((hi[a] + pad - origin[a]) / h).ceil() as usize + 1).max(2));
        let total = dims[0] * dims[1] * dims[2];
        if total > 64_000_000 {
            return Err(self.err(e, format!("grid of {dims:?} points is too large. Increase (res ...)")));
        }
        let point = |n: usize| -> V3 {
            let [i, j, k] = crate::grid::ijk(n, dims);
            [origin[0] + h * i as f64, origin[1] + h * j as f64, origin[2] + h * k as f64]
        };
        let pts: Vec<V3> = (0..total).map(point).collect();
        // A mesh's lines signed by its facets' normals, not by parity
        // (MeshDistance::sample).
        let mut unclear = 0;
        // Each run of points taken only while the call's time limit has not
        // passed (`stop`): past it, the load stops here.
        let in_time = |each: &(dyn Fn(usize, &[V3]) -> Vec<f64> + Sync)| -> Option<Vec<f64>> {
            let runs: Vec<Option<Vec<f64>>> = pts.par_chunks(4096).enumerate().map(|(k, c)| (!ink_kernel::stop::reached()).then(|| each(4096 * k, c))).collect();
            runs.into_iter().collect::<Option<Vec<_>>>().map(|r| r.concat())
        };
        let samples: Vec<f64> = match &source {
            Src::Field(t) => in_time(&|_, c| {
                let mut out = Vec::with_capacity(c.len());
                crate::eval::eval_points(t, &ctx, c, &mut out);
                out
            }),
            Src::Mesh(md) => {
                let (v, u) = md.sample(origin, h, dims);
                unclear = u;
                Some(v)
            }
            Src::Surface(md, t) => ink_kernel::octree::grid_inside(t, &ctx, origin, h, dims).and_then(|(inside, _)| {
                in_time(&|at, c| {
                    // The mesh lies within a point's neighbour's distance
                    // plus the step between them (and a hair, for rounding).
                    let mut last: Option<(V3, f64)> = None;
                    c.iter()
                        .zip(&inside[at..])
                        .map(|(p, i)| {
                            let reach = last.map_or(f64::INFINITY, |(q, d)| (d + ink_core::v3::norm(ink_core::v3::sub(*p, q))) * (1.0 + 1e-6) + 1e-9);
                            let d = md.distance_within(*p, reach).abs();
                            last = Some((*p, d));
                            if *i { -d } else { d }
                        })
                        .collect()
                })
            }),
        }
        .filter(|_| !ink_kernel::stop::reached())
        .ok_or_else(|| self.stopped(e))?;
        // A grid keeps f32s: a finite sample past their range would be
        // infinite, and the shape lost.
        if let Some(v) = samples.iter().find(|v| v.is_finite() && !(**v as f32).is_finite()) {
            return Err(self.err_hint(e, format!("{name}: a sample is {v:e}, past the range of the grid's numbers (about 3.4e38)"), Some("scale the field down".into())));
        }
        let values: Vec<f32> = samples.into_iter().map(|v| v as f32).collect();
        // How far the source may be from the samples: of a field, each
        // sample's distance from the field's real value (an interval at
        // the point) plus L times how far a grid point computed in floats
        // may be from where the interpolant puts it (generously: 1e-9
        // (1 + the largest coordinate)).
        let scale = (0..3).map(|a| origin[a].abs().max((origin[a] + h * (dims[a] - 1) as f64).abs())).fold(0.0, f64::max);
        let geometry_round = 1e-9 * (1.0 + scale);
        // Probes for the deviation actually measured: the centres of the
        // cells the surface may cross (a corner within a cell diagonal of
        // 0), where interpolation errs most, at most 250 000 of them, evenly.
        let cells: [usize; 3] = std::array::from_fn(|a| dims[a] - 1);
        let at = |i: usize, j: usize, k: usize| values[crate::grid::index([i, j, k], dims)] as f64;
        let near: Vec<([usize; 3], V3, f64)> = (0..cells[0] * cells[1] * cells[2])
            .into_par_iter()
            .filter_map(|n| {
                let [i, j, k] = crate::grid::ijk(n, cells);
                let c: Vec<f64> = (0..8).map(|m| at(i + (m & 1), j + ((m >> 1) & 1), k + (m >> 2))).collect();
                let reach = c.iter().fold(f64::INFINITY, |m, v| m.min(v.abs()));
                (reach < 1.8 * h)
                    .then(|| ([i, j, k], [origin[0] + h * (i as f64 + 0.5), origin[1] + h * (j as f64 + 0.5), origin[2] + h * (k as f64 + 0.5)], c.iter().sum::<f64>() / 8.0))
            })
            .collect();
        let stride = near.len().div_ceil(250_000).max(1);
        let probes: Vec<([usize; 3], V3, f64)> = near.into_iter().step_by(stride).collect();
        // The largest |samples - truth| over the probes the truth is known
        // at (given each probe's cell and point), and how many those are.
        let deviation = |truth: &(dyn Fn([usize; 3], V3) -> Option<f64> + Sync)| -> Option<(f64, usize)> {
            let (m, k) = probes
                .par_iter()
                .filter_map(|(c, p, i)| if ink_kernel::stop::reached() { None } else { truth(*c, *p).map(|t| (i - t).abs()) })
                .fold(|| (0.0, 0), |(m, k), d| (f64::max(m, d), k + 1))
                .reduce(|| (0.0, 0), |a, b| (a.0.max(b.0), a.1 + b.1));
            (k > 0).then_some((m, k))
        };
        // The samples read at points of a mesh (where its distance is 0):
        // its vertices, where sharp features are, and centroids, at most
        // 250 000, evenly.
        let interpolant = |p: V3| -> f64 {
            let mut idx = [0usize; 3];
            let mut t = [0.0; 3];
            for a in 0..3 {
                let u = ((p[a] - origin[a]) / h).clamp(0.0, (dims[a] - 1) as f64);
                idx[a] = (u.floor() as usize).min(dims[a] - 2);
                t[a] = u - idx[a] as f64;
            }
            (0..8)
                .map(|m: usize| {
                    let o = [m & 1, (m >> 1) & 1, m >> 2];
                    let w: f64 = (0..3).map(|a| if o[a] == 1 { t[a] } else { 1.0 - t[a] }).product();
                    w * at(idx[0] + o[0], idx[1] + o[1], idx[2] + o[2])
                })
                .sum()
        };
        let on_mesh = |md: &crate::grid::MeshDistance| -> (f64, usize) {
            // Those in the part of the mesh's box sampled.
            let pts: Vec<V3> = md.surface_points().into_iter().filter(|p| (0..3).all(|a| lo[a] <= p[a] && p[a] <= hi[a])).collect();
            let stride = pts.len().div_ceil(250_000).max(1);
            let pts: Vec<V3> = pts.into_iter().step_by(stride).collect();
            (pts.par_iter().map(|p| interpolant(*p).abs()).reduce(|| 0.0, f64::max), pts.len())
        };
        // A field that reads sampled data stands for what that was sampled
        // from, which its own bound does not reach (`inner`).
        let (band, measured, note, cut) = match &source {
            Src::Field(t) if !inner.is_empty() => (
                None,
                deviation(&|_, p| Some(crate::eval::eval_point(t, &ctx, p))),
                Some(format!("the field baked reads sampled data ({}): no error bound to what that was sampled from", inner.join(", "))),
                None,
            ),
            Src::Mesh(md) => {
                let whole = region.is_none_or(|(rlo, rhi)| {
                    let (mlo, mhi) = md.bounds();
                    (0..3).all(|k| rlo[k] <= mlo[k] && mhi[k] <= rhi[k])
                });
                let odd = md.odd_edges();
                let note = match (odd, unclear) {
                    (0, 0) => SIGNED.to_string(),
                    (0, k) => format!("{SIGNED}; {k} lines of the grid meet an edge or vertex of the mesh: their points are signed by its facets' normals"),
                    _ => format!("the mesh is not closed ({odd} edges with an odd number of triangles), so it bounds no solid: sampled as its facets face"),
                };
                // A region import is cut to the region, kept inside the
                // grid by half its margin so the cut never reaches the
                // grid's ends (where a source not whole is unknown).
                let cut = (!whole).then(|| {
                    let (rlo, rhi) = region.unwrap();
                    (std::array::from_fn(|k| rlo[k].max(lo[k] - 0.5 * pad)), std::array::from_fn(|k| rhi[k].min(hi[k] + 0.5 * pad)))
                });
                // The deviation measured: at cell centres near the surface,
                // the samples against the source (MeshDistance::signed); at
                // points of the mesh, against 0.
                let measured = deviation(&|_, p| Some(md.signed(p))).map(|(m, k)| {
                    let (s, j) = on_mesh(md);
                    (m.max(s), k + j)
                });
                (None, measured, Some(note), cut)
            }
            Src::Field(t) => {
                // A field with a slope bound L over the grid box: the
                // interpolation bound, plus how far each sample is from the
                // real-number value of the field (its interval at the point),
                // plus L times how far a grid point, computed in floats, may
                // be from where the interpolant puts it (generously).
                let top: V3 = std::array::from_fn(|a| origin[a] + h * (dims[a] - 1) as f64);
                let gb: [I; 3] = std::array::from_fn(|a| I::new(origin[a], top[a]));
                let g = crate::eval::eval_interval_grad(t, &ctx, gb);
                let l = interval::lipschitz(&g);
                let band = l.is_finite().then(|| {
                    let round = pts
                        .par_iter()
                        .zip(&values)
                        .map(|(p, v)| {
                            if ink_kernel::stop::reached() {
                                return 0.0;
                            }
                            let r = crate::eval::eval_interval(t, &ctx, p.map(I::point), None);
                            (*v as f64 - r.lo).abs().max((r.hi - *v as f64).abs())
                        })
                        .reduce(|| 0.0, f64::max);
                    let place = interval::mul_up(l, geometry_round);
                    (interval::add_up(interval::add_up(crate::grid::Source::interpolation_bound(h, l), round.next_up()), place), None)
                });
                let note = match band {
                    Some(_) => format!("the field's slope is at most {l:.6} over the grid (interval bounds)"),
                    None => "the field's slope is not bounded over the grid: no error bound".into(),
                };
                (band, deviation(&|_, p| Some(crate::eval::eval_point(t, &ctx, p))), Some(note), None)
            }
            Src::Surface(..) => (
                None,
                None,
                Some("the distance to a mesh of the field's surface (meshed at half the grid spacing), whose distance from that surface is not bounded: no error bound".into()),
                None,
            ),
        };
        // Past the time limit a pass above stopped short (what it measured
        // is not the grid's): the load stops here, the grid not kept.
        if ink_kernel::stop::reached() {
            return Err(self.stopped(e));
        }
        let open = matches!(&source, Src::Mesh(md) if md.odd_edges() > 0);
        let snap = Snap { origin, h, dims, values, band, measured, note, cut, open };
        crate::cache::computed();
        crate::cache::store("grids", key, &snap);
        Snap::keep(key, &snap);
        Ok(finish(self, snap))
    }

    /// The context a load-time service finds what it finds of `f` over (a
    /// box, a reach), so that holds wherever `f` is asked: the load's, with
    /// t over the motion of the file whose form it is (not std's, inside a
    /// call) and 0, where load reads it, if `f` reads t, and each ranged
    /// param it reads over its band and its value. What it finds is a
    /// number of the model, so the load holds over these alone: a run
    /// beyond them is refused (`load_t`, `covered`); a param with no range
    /// is pinned (`bake`). Err: `f` reads t and no motion is stated before.
    fn load_scope(&mut self, e: &Sexp, f: NodeId, what: &str) -> R<Ctx> {
        let mut ctx = self.load_ctx_for(&[f]);
        if self.g.value_depends_on(f, |v| v == Var::T) {
            let file = self.running.last().map_or(self.file.clone(), |&i| self.files[i as usize].name.clone());
            let Some((a, b)) = self.span_of(&file) else {
                return Err(self.err_hint(
                    e,
                    format!("{what}: found at load, it holds over the file's motion, and none is stated before it"),
                    Some("state (t-span A B) before it, or bound the shape before it moves (one that does not read t)".into()),
                ));
            };
            let span = (a.min(0.0), b.max(0.0));
            ctx.t_span = Some(span);
            self.load_t.push((span, self.form_at.clone()));
        }
        let ks = self.params_read(&[f]);
        let ranged: Vec<(usize, (f64, f64))> = ks.iter().filter_map(|&k| Some((k, self.load_band(k, &ctx)?))).collect();
        self.pin(&ks.iter().copied().filter(|k| !ranged.iter().any(|r| r.0 == *k)).collect::<Vec<_>>());
        for &(k, band) in &ranged {
            self.covered.push((k, band, self.form_at.clone()));
            // Point evaluation inside the span, as a proof keeps it.
            ctx.params[k] = ink_core::eval::span_mid(band);
        }
        self.baked.borrow_mut().extend(ranged.iter().map(|r| r.0));
        ctx.p_spans = ink_core::eval::ParamSpans::new(ranged.iter().map(|&(k, b)| (k as u32, b)));
        Ok(ctx)
    }

    /// The values of param k in `ctx` a run may use, where it has a band:
    /// its value and its tolerance's band (which may exclude it), or its
    /// (ranged lo hi).
    fn load_band(&self, k: usize, ctx: &Ctx) -> Option<(f64, f64)> {
        let v = ctx.params[k];
        let (lo, hi) = self.ranged.iter().find(|r| r.0 == k).map(|r| r.1).or_else(|| self.params[k].band(v))?;
        Some((lo.min(v), hi.max(v)))
    }

    /// `eval` forms that call kernel services at load time.
    pub(super) fn eval_services(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        match name {
            "grid-bake" | "import-stl" | "redistance" => return self.sampled_grid(name, e, args, env),
            "bounds" => {
                let usage = "usage: (bounds field [default]) gives (vec (vec x0 y0 z0) (vec x1 y1 z1)); default (evaluated only then) where field is empty";
                let (f, default) = match args {
                    [f] => (f, None),
                    [f, d] => (f, Some(d)),
                    _ => return Err(self.err(e, usage)),
                };
                let fe = f;
                let f = self.field(f, env)?;
                // Its numbers would hide a solved value; of sampled data they
                // are measured (measured_box).
                let sampled = self.inexact(&[f], &[]);
                if !sampled.solved.is_empty() {
                    self.exact_at_load(fe, f)?;
                }
                let ctx = self.load_scope(fe, f, "bounds of a field that reads t")?;
                let (clipped, reach) = self.clip_bounded(e, "bounds", f, &ctx)?;
                let tape = std::sync::Arc::new(Tape::new(&self.g, clipped));
                // Every answer, decided here once: what was found of sampled
                // data is measured, so a box of it is (measured_box), and so
                // is finding it empty, which a default (a value of any kind)
                // cannot say: refused.
                let found = crate::query::bounds(&tape, &ctx, reach.lo, reach.hi);
                let v3 = |d: &mut Doc, p: V3| Value::list(p.iter().map(|c| Value::Field(d.g.constant(*c), None)).collect());
                return match (found, sampled.samples.is_empty(), default) {
                    (Some((lo, hi, _)), true, _) => Ok(Value::list(vec![v3(self, lo), v3(self, hi)])),
                    (Some((lo, hi, _)), false, _) => Ok(self.measured_box(e, fe, lo, hi, &sampled.samples)),
                    (None, true, Some(d)) => self.eval(d, env),
                    (None, true, None) => Err(self.err(e, "bounds of an empty shape")),
                    (None, false, _) => Err(self.err_hint(
                        e,
                        format!(
                            "bounds found no interior in sampled data ({}): whether it is empty is measured, so no default stands in for it",
                            crate::inexact::sampled_from(&sampled.samples)
                        ),
                        Some("bound the shape it was sampled from, or claim it empty: (assert none (empty f))".into()),
                    )),
                };
            }
            "profile-reach" => {
                let usage = "usage: (profile-reach profile): how far a 2D profile (a field in x, y) reaches from its origin, proven; \
                             one that reads z is taken over every z in [-1, 1] (sweep-varying's sections)";
                let [f] = args else { return Err(self.err_hint(e, "bad profile-reach", Some(usage.into()))) };
                let fe = f;
                let f = self.field(f, env)?;
                self.exact_at_load(fe, f)?;
                let ctx = self.load_scope(fe, f, "profile-reach of a profile that reads t")?;
                let reads_z = self.g.value_depends_on(f, |v| matches!(v, Var::Z));
                // Bounded in z: the profile over z in [-1, 1].
                let (z, one) = (self.g.var(Var::Z), self.g.constant(1.0));
                let az = self.g.un(UnOp::Abs, z);
                let slab = self.g.bin(BinOp::Sub, az, one);
                let f = self.g.bin(BinOp::Max, f, slab);
                let (clipped, reach) = self.clip_bounded(e, "profile-reach", f, &ctx)?;
                let tape = std::sync::Arc::new(Tape::new(&self.g, clipped));
                let Some((mut lo, mut hi, _)) = crate::query::bounds(&tape, &ctx, reach.lo, reach.hi) else {
                    return Err(self.err_hint(e, "profile-reach of an empty profile", Some(usage.into())));
                };
                // A profile that does not read z is the same at every z: one
                // is enough.
                if !reads_z {
                    (lo[2], hi[2]) = (0.0, 0.0);
                }
                // To a thousandth of its size: far below what a reach is for
                // (where a section's farthest point is, to bound a sweep).
                let size = (hi[0] - lo[0]).max(hi[1] - lo[1]);
                let r = crate::query::radial_extent(&tape, &ctx, lo, hi, reads_z, 1e-3 * size, 1 << 16).map_or(0.0, |e| e.1);
                return Ok(Value::Field(self.g.constant(r), None));
            }
            _ => {}
        }
        unreachable!()
    }
}
