//! Reaction providers and access to whole reactions and support shares.

use crate::model::Num;
use crate::providers::{FIXED_SHORTFALL, Shared, unseparated};
use crate::supports::Support;
use crate::{Ctx, Model, Moment, Serve, V3, command, model, split};
use ink_core::graph::{NodeId, Var};
use ink_lang::{Doc, Form, ParamProvider, Reader, Refused, Steps, Xf};
use std::sync::Arc;

/// A frame compiled, with the params it reads (of `n`).
pub(super) fn compiled(r: &mut Reader, x: &Xf, n: usize) -> (Arc<[Num; 12]>, Vec<usize>) {
    let gr = r.graph();
    let entries: [NodeId; 12] = std::array::from_fn(|k| if k < 9 { x.r[k / 3][k % 3] } else { x.t[k - 9] });
    let reads: Vec<usize> = (0..n).filter(|&k| entries.iter().any(|&e| gr.depends_on(e, |v| v == Var::Param(k as u32)))).collect();
    (Arc::new(entries.map(|e| Num::of(gr, e))), reads)
}

/// The form keeping one reaction's provider: its index among the shared
/// model's reactions, the frame its claim names (compiled, with the params
/// it reads), and the six params it sets; or, split among supports, six for
/// each, and its path's joint's drive.
#[derive(Clone)]
pub(super) struct Reacting {
    pub(super) shared: Arc<Shared>,
    pub(super) index: usize,
    /// Its place among the reactions claims ask for: its steps' key.
    pub(super) claim: usize,
    pub(super) frame: Option<(Arc<[Num; 12]>, Vec<usize>)>,
    pub(super) owns: Vec<usize>,
    pub(super) supports: Vec<Support>,
    pub(super) drive: Option<split::Drive>,
    /// Said by `(supports ...)`: no expression reads its shares, and a located shaft's thrust, which is in
    /// none of them, is its bearings' to carry (split::split_located), not a force nothing carries.
    pub(super) named: bool,
}

impl Form for Reacting {
    fn name(&self) -> &str {
        "reaction"
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for Reacting {
    fn name(&self) -> &str {
        "the reaction on what a path carries along the motion (inverse dynamics of its bodies; (reaction PATH))"
    }

    fn owns(&self) -> &[usize] {
        &self.owns
    }

    fn uses_t(&self) -> bool {
        true
    }

    fn on_demand(&self) -> bool {
        true
    }

    fn warm_starts(&self) -> bool {
        false
    }

    /// Its real force and moment at t, in its frame; none where it is not
    /// this motion's (Model::reaction_refusal). What a stepped claim reads
    /// as its steps read it, round, comes with its band (Reacting::band).
    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, Refused> {
        let model = &self.shared.model()?.0;
        // What a simulation's springs move it carries: not the rigid
        // model's reaction.
        if let Some(sim) = self.shared.springs() {
            let moved = model.carried_by(self.index, &sim?.model().sprung_bodies());
            if !moved.is_empty() {
                return Err(Refused::from(format!(
                    "it carries what the file's simulation moves by its springs ({}): its reaction is not the rigid model's along the written motion, \
                     and is not found",
                    moved.join(", ")
                )));
            }
        }
        let m = model.at(t);
        if let Some(why) = model.reaction_refused(self.index, &m) {
            return Err(why);
        }
        for (k, v) in self.owns.iter().zip(self.values(&m, params, t)?.0) {
            params[*k] = v;
        }
        let unsettled: Vec<&str> = model.unsettled.iter().chain(&model.unsettled_carried).map(String::as_str).collect();
        Ok((!unsettled.is_empty()).then(|| format!("masses it rests on are not measured to 2%: {}", unsettled.join("; "))))
    }

    fn reads(&self) -> Option<Vec<usize>> {
        let (model, _) = self.shared.model().ok()?;
        let mut out: Vec<usize> = model.reads.iter().chain(&model.reads_carried).filter_map(|n| self.shared.names.iter().position(|m| m == n)).collect();
        out.extend(self.frame.iter().chain(self.supports.iter().flat_map(|s| std::iter::once(&s.load).chain(&s.given))).flat_map(|f| f.1.iter().copied()));
        Some(out)
    }

    /// Its own steps (Model::steps_for): no other claim's follow them.
    fn steps(&self, params: &[f64], span: (f64, f64)) -> Option<Result<Steps, Refused>> {
        Some(self.shared.model().map(|(m, _)| self.served(params, |s| m.steps_for(span, s))))
    }

    /// Where its steps read groups it carries at rest, round: their bound
    /// on its force (each component, in any frame) and on its moment about
    /// its frame's origin (the moment about each group's axis point and
    /// that force's arm from there), and its value so read, the band's
    /// centre.
    fn band(&self, params: &[f64], t: f64) -> Option<(Vec<f64>, Vec<f64>)> {
        let (model, _) = self.shared.model().ok()?;
        let (m, shakes) = self.read(model, params, t);
        if shakes.is_empty() {
            return None;
        }
        let w = &m.reactions[self.index];
        // A split's map is from the moment about its path's frame's origin.
        let o: V3 = match &self.frame {
            Some((frame, _)) => std::array::from_fn(|k| frame[9 + k].at(&Ctx::new(t, params.to_vec()))),
            None => std::array::from_fn(|k| w.frame[9 + k]),
        };
        let force: f64 = shakes.iter().map(|s| s.0).sum();
        let moment: f64 = shakes.iter().map(|&(f, mo, at)| mo + f * ink_core::v3::norm(std::array::from_fn(|k| at[k] - o[k]))).sum();
        let (centre, map) = self.values(&m, params, t).ok()?;
        let widths = match map {
            None => [force; 3].into_iter().chain([moment; 3]).collect(),
            // Each share is linear in the reaction: its width is its row's
            // reach over the reaction's band.
            Some(map) => map.iter().flatten().map(|row| (0..6).map(|j| row[j].abs() * if j < 3 { force } else { moment }).sum()).collect(),
        };
        Some((widths, centre))
    }

    /// Where a fixed-effort drive's required effort differs from its rating
    /// within the rounding of its load (its effort claim unknown there), the
    /// motion is that drive's to that rounding: its reaction is known to
    /// the shortfall's share of its own rounding either way (Model::
    /// shortfall), each component's, through its split each share's.
    fn error(&self, params: &[f64], t: f64) -> Option<Vec<f64>> {
        let (model, _) = self.shared.model().ok()?;
        let m = model.at(t);
        let share = model.shortfall(self.index, &m).max(0.0);
        let w = &m.reactions[self.index];
        let (nf, nm) = w.rounding();
        let frame: [f64; 12] = match &self.frame {
            Some((frame, _)) => frame.each_ref().map(|n| n.at(&Ctx::new(t, params.to_vec()))),
            None => w.frame,
        };
        let (force, moment) = (share * nf, share * (nm + nf * ink_core::v3::norm([frame[9], frame[10], frame[11]])));
        let rounding = |j: usize| if j < 3 { force } else { moment };
        // And what its bodies' measured masses, centres and inertias can move each component by (Moment::mass_reach).
        let (_, map) = self.values(&m, params, t).ok()?;
        let widths: Vec<f64> = match map {
            None => m.mass_reach(self.index, &frame).iter().enumerate().map(|(j, x)| rounding(j) + x).collect(),
            Some(map) => {
                let rows: Vec<[f64; 6]> = map.iter().flatten().copied().collect();
                let mass = m.mass_reach_mapped(self.index, &rows);
                rows.iter().zip(mass).map(|(row, x)| (0..6).map(|j| row[j].abs() * rounding(j)).sum::<f64>() + x).collect()
            }
        };
        widths.iter().any(|x| *x != 0.0).then_some(widths)
    }

    /// The bodies moving with its path that it does not carry, once the
    /// model is built (Reaction::moves_with).
    fn note(&self) -> Option<String> {
        // Carriages of one record sharing as a rigid table say so.
        let parts: Vec<(split::Carries, &str)> = self.supports.iter().map(|s| (s.carries, s.rating.as_str())).collect();
        let guides: Vec<&str> = self.supports.iter().filter(|s| s.carries == split::Carries::Guide).map(|s| s.path.as_str()).collect();
        let table = split::table(&parts).then(|| format!("the shares of {} {}", guides.join(", "), split::TABLE));
        let Some(Ok((model, _))) = self.shared.model.get() else { return table };
        let r = &model.reactions[self.index];
        let with = (!r.moves_with.is_empty()).then(|| {
            format!(
                "the reaction of {} sums, unseparated, what these bodies exert on what it carries: they move with it but are not placed on its frames \
                 (moves_with): {}; it is not the joint's load alone",
                r.path,
                r.moves_with.join(", ")
            )
        });
        // A carriage's moments rest on where its slide's drive pulls, which
        // no part named says, nor a coupling pulling along its own line.
        let line = (self.drive.is_some_and(|d| !d.turn)
            && self.supports.iter().any(|s| s.carries == split::Carries::Guide)
            && !self.supports.iter().any(|s| s.carries == split::Carries::Axial)
            && !model.coupled(self.index))
        .then(|| {
            "the slide's drive is taken along its axis through its joint frame's origin: a belt or screw pulling h from there moves h times its force \
                 of the carriages' pitch and yaw moments onto them (name the screw's nut on its line to have it counted); a carriage's moments are about \
                 its load frame, on its rail under its centre, not its ball tracks' height"
                .to_string()
        });
        let pulled = unseparated(model).filter(|_| model.coupled(self.index));
        let fixed = model.fixed_reaches(self.index).then(|| FIXED_SHORTFALL.to_string());
        let loose = model.bears_unpaired(self.index).then(|| UNPAIRED.to_string());
        let notes: Vec<String> =
            with.into_iter().chain(line).chain(pulled).chain(table).chain(fixed).chain(loose).chain(model.spin_centre_note(self.index)).chain(model.moved_by_t_note()).collect();
        (!notes.is_empty()).then(|| notes.join("; "))
    }
}

impl Reacting {
    /// Its moment as its claims' steps over the file's span read it, the
    /// groups they read at rest read round, with each one's bound
    /// (Model::shaken); outside the span, every body as measured.
    fn read(&self, model: &Model, params: &[f64], t: f64) -> (Arc<Moment>, Vec<(f64, f64, V3)>) {
        match model.span {
            Some(s) if s.0 <= t && t <= s.1 => self.served(params, |serve| model.shaken(s, Some(serve), self.index, t)),
            // At one moment its groups are read as they are; what is not known of them is still not.
            _ => (model.at(t), model.unknown_at(self.index, t)),
        }
    }

    /// `f` given what its steps serve: its reaction, itself, and the frame
    /// its components are followed in at each moment (read at `params`):
    /// its claim's, or its first support's.
    fn served<T>(&self, params: &[f64], f: impl FnOnce(&Serve) -> T) -> T {
        let followed = self.frame.as_ref().or(self.supports.first().map(|s| s.given.as_ref().unwrap_or(&s.load)));
        let frame = |t: f64| -> [f64; 12] {
            let ctx = Ctx::new(t, params.to_vec());
            followed.map_or([0.0; 12], |(fr, _)| fr.each_ref().map(|n| n.at(&ctx)))
        };
        let serve = Serve { reaction: self.index, claim: self.claim, frame: followed.map(|_| &frame as &(dyn Fn(f64) -> [f64; 12] + Sync)) };
        f(&serve)
    }

    /// What it sets at moment m: its force and moment (Reacting::wrench), or
    /// each support's share in its frame, with the map from the reaction to
    /// them (split::split); Err where the split is refused.
    fn values(&self, m: &Moment, params: &[f64], t: f64) -> Result<(Vec<f64>, Option<Vec<[[f64; 6]; 6]>>), String> {
        if self.supports.is_empty() {
            let (f, mo) = self.wrench(m, params, t);
            return Ok((f.into_iter().chain(mo).collect(), None));
        }
        let model = &self.shared.model()?.0;
        let ctx = Ctx::new(t, params.to_vec());
        let at = |f: &(Arc<[Num; 12]>, Vec<usize>)| f.0.each_ref().map(|n| n.at(&ctx));
        let parts: Vec<split::Part> = self
            .supports
            .iter()
            .map(|s| {
                let load = at(&s.load);
                split::Part {
                    path: &s.path,
                    carries: s.carries,
                    pulls: s.pulls,
                    load,
                    given: s.given.as_ref().map_or(load, at),
                    coupled: model.coupled_on(&s.path, m),
                    rating: &s.rating,
                }
            })
            .collect();
        let (shares, map, _) = split::split_located(&m.reactions[self.index], &parts, self.drive, self.named).map_err(|why| format!("its split among its supports: {why}"))?;
        Ok((shares.into_iter().flatten().collect(), Some(map)))
    }

    /// Its force and moment at moment m, in its claim's frame (the frame
    /// read at t and `params`), else its path's.
    fn wrench(&self, m: &Moment, params: &[f64], t: f64) -> (V3, V3) {
        let w = &m.reactions[self.index];
        match &self.frame {
            Some((frame, _)) => {
                let ctx = Ctx::new(t, params.to_vec());
                w.in_frame(&frame.each_ref().map(|n| n.at(&ctx)))
            }
            None => w.local(),
        }
    }
}

/// A reaction the design divides among supports, `(reaction PATH (by S ...))`, for a caller that reads what
/// each support carries off a run of its own (ink-effort's facts): the path, and the supports' paths in order.
pub struct Divided {
    pub path: String,
    pub supports: Vec<String>,
    /// Of each support, whether it is a radial part that carries thrust too (its record's `thrust`,
    /// double-direction: a deep-groove ball bearing): one of those that may locate the shaft.
    pub locates: Vec<bool>,
    /// Of each support, whether it is a guide's carriage (moment ratings): its share is a force across its
    /// load frame's x, its travel, and a moment about each of its axes.
    pub guides: Vec<bool>,
    by: Reacting,
}

/// Every reaction the design divides among supports, once each (a path divided twice among the same parts is one).
pub fn divided(doc: &Doc) -> Vec<Divided> {
    let mut out: Vec<Divided> = Vec::new();
    for by in doc.forms::<Reacting>().filter(|r| !r.supports.is_empty()) {
        let (path, supports) = (by.shared.reactions[by.index].clone(), by.supports.iter().map(|s| s.path.clone()).collect::<Vec<_>>());
        if !out.iter().any(|d| d.path == path && d.supports == supports) {
            let locates = by.supports.iter().map(|s| s.carries == split::Carries::Radial && s.pulls).collect();
            let guides = by.supports.iter().map(|s| s.carries == split::Carries::Guide).collect();
            out.push(Divided { path, supports, locates, guides, by: by.clone() });
        }
    }
    out
}

impl Divided {
    /// What each support carries at moment t of `run`, as it is there (no group read round), `(fx fy fz mx my
    /// mz)` in its load frame about its origin, N and N mm, the design's params at `params`: the claim's own
    /// split (Reacting::values), each share in the frame its rating acts in. None where the run does not
    /// carry its path's reaction; Err where the reaction or the split is refused.
    ///
    /// With them, the thrust its radial parts that locate the shaft carry among them, where they are what
    /// carries it (split::split_located): it is in none of their shares, since statics does not say which.
    pub fn shares(&self, run: &command::Read, t: f64, params: &[f64]) -> Option<Result<(Vec<[f64; 6]>, Option<f64>), String>> {
        let index = run.model.reactions.iter().position(|r| r.path == self.path)?;
        let m = run.model.at(t);
        if let Some(why) = run.model.reaction_refused_lossless(index, &m) {
            return Some(Err(why.why));
        }
        let ctx = Ctx::new(t, params.to_vec());
        let parts: Vec<split::Part> = self
            .by
            .supports
            .iter()
            .map(|s| {
                let load = s.load.0.each_ref().map(|n| n.at(&ctx));
                split::Part { path: &s.path, carries: s.carries, pulls: s.pulls, load, given: load, coupled: run.model.coupled_on(&s.path, &m), rating: &s.rating }
            })
            .collect();
        Some(
            split::split_located(&m.reactions[index], &parts, self.by.drive, true)
                .map(|(shares, _, thrust)| (shares, thrust))
                .map_err(|why| format!("its split among its supports: {why}")),
        )
    }

    /// What each change its reaction's masses may make (Moment::mass_w: each body's measured moment
    /// anywhere in its proven bounds) moves each support's share and the locating thrust by, at moment t:
    /// the split of the reaction so changed less the split as read, the split being linear in it. Empty
    /// where `shares` gives none.
    pub fn share_moves(&self, run: &command::Read, t: f64, params: &[f64]) -> Vec<(Vec<[f64; 6]>, f64)> {
        let Some(index) = run.model.reactions.iter().position(|r| r.path == self.path) else { return Vec::new() };
        let m = run.model.at(t);
        let ctx = Ctx::new(t, params.to_vec());
        let parts: Vec<split::Part> = self
            .by
            .supports
            .iter()
            .map(|s| {
                let load = s.load.0.each_ref().map(|n| n.at(&ctx));
                split::Part { path: &s.path, carries: s.carries, pulls: s.pulls, load, given: load, coupled: run.model.coupled_on(&s.path, &m), rating: &s.rating }
            })
            .collect();
        let w = &m.reactions[index];
        let Ok((base, _, thrust)) = split::split_located(w, &parts, self.by.drive, true) else { return Vec::new() };
        m.mass_w[index]
            .iter()
            .filter_map(|(df, dm)| {
                let moved = model::Wrench::of(std::array::from_fn(|k| w.force[k] + df[k]), std::array::from_fn(|k| w.moment[k] + dm[k]), w.frame);
                let (shares, _, t2) = split::split_located(&moved, &parts, self.by.drive, true).ok()?;
                let d = shares.iter().zip(&base).map(|(a, b)| std::array::from_fn(|i| a[i] - b[i])).collect();
                Some((d, t2.unwrap_or(0.0) - thrust.unwrap_or(0.0)))
            })
            .collect()
    }
}

/// What each change its masses may make (Moment::mass_w) moves joint `path`'s reaction by at moment t of
/// `run`, force and moment in `local` (as `reaction_in`), or in its path's own frame with none.
pub fn reaction_moves(run: &command::Read, path: &str, t: f64, local: Option<&[f64; 12]>) -> Vec<([f64; 3], [f64; 3])> {
    let Some(index) = run.model.reactions.iter().position(|r| r.path == path) else { return Vec::new() };
    let m = run.model.at(t);
    let frame = m.reactions[index].frame;
    m.mass_w[index]
        .iter()
        .map(|(df, dm)| {
            let w = model::Wrench::of(*df, *dm, frame);
            match local {
                Some(l) => w.in_local(l),
                None => w.local(),
            }
        })
        .collect()
}

/// How far joint `path`'s reaction at moment t of `run` may be from the one read: over a span the solver
/// reads bodies that spin as if round, at rest, and says the most their departure from round can move the
/// reaction's force (N) and its moment (N mm) at any angle of their spin (`ink effort`'s `bound`). (0, 0)
/// with no span, where nothing is read so.
pub fn reaction_bound(run: &command::Read, path: &str, t: f64) -> (f64, f64) {
    let Some(r) = run.model.reactions.iter().position(|r| r.path == path) else { return (0.0, 0.0) };
    // (With no span, what is not known of the groups it carries, which is so at every moment.)
    let (w, shakes) = match run.span {
        Some(span) => run.model.shaken(span, None, r, t),
        None => (run.model.at(t), run.model.unknown_at(r, t)),
    };
    let w = &w.reactions[r];
    let o = [w.frame[9], w.frame[10], w.frame[11]];
    let far = |p: [f64; 3]| (0..3).map(|i| (p[i] - o[i]).powi(2)).sum::<f64>().sqrt();
    (shakes.iter().map(|s| s.0).sum(), shakes.iter().map(|&(f, mo, p)| mo + f * far(p)).sum())
}

impl Divided {
    /// The books of this division at moment t of `run`: what its supports were given, summed, against what
    /// the joint needs held. Each share turned from its support's load frame into the world; the joint's
    /// reaction as the model has it. (the force given, the force needed, the moment given about the world's
    /// origin, the moment needed, the thrust the locating supports carry among them, which is in no share).
    /// A check of the split's own arithmetic, by adding up what it handed out. None or Err as `shares`.
    ///
    /// With them, the joint's own axis in the world and whether it turns: what the joint's drive carries along
    /// it (a motor's torque about a turn's axis, a slide's force along it) is no support's.
    #[allow(clippy::type_complexity)]
    pub fn books(&self, run: &command::Read, t: f64, params: &[f64]) -> Option<Result<([[f64; 3]; 4], Option<(bool, [f64; 3])>), String>> {
        let index = run.model.reactions.iter().position(|r| r.path == self.path)?;
        let (shares, _) = match self.shares(run, t, params)? {
            Ok(found) => found,
            Err(why) => return Some(Err(why)),
        };
        let m = run.model.at(t);
        let ctx = Ctx::new(t, params.to_vec());
        let (mut force, mut moment) = ([0.0; 3], [0.0; 3]);
        for (s, share) in self.by.supports.iter().zip(&shares) {
            let f = s.load.0.each_ref().map(|n| n.at(&ctx));
            let turn = |v: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| f[3 * i] * v[0] + f[3 * i + 1] * v[1] + f[3 * i + 2] * v[2]) };
            let (fw, mw, o) = (turn([share[0], share[1], share[2]]), turn([share[3], share[4], share[5]]), [f[9], f[10], f[11]]);
            let lever = [o[1] * fw[2] - o[2] * fw[1], o[2] * fw[0] - o[0] * fw[2], o[0] * fw[1] - o[1] * fw[0]];
            for i in 0..3 {
                force[i] += fw[i];
                moment[i] += mw[i] + lever[i];
            }
        }
        let needed = &m.reactions[index];
        let axis = self.by.drive.map(|d| (d.turn, [needed.frame[d.axis], needed.frame[3 + d.axis], needed.frame[6 + d.axis]]));
        Some(Ok(([force, needed.force, moment, needed.moment], axis)))
    }

    /// What support `k` of this division was given at moment t, `(fx fy fz mx my mz)` in its load frame (`shares`).
    pub fn given(&self, run: &command::Read, t: f64, params: &[f64], k: usize) -> Option<[f64; 6]> {
        self.shares(run, t, params)?.ok().map(|(shares, _)| shares[k])
    }
}

/// What the rest of the model exerts on what joint `path` carries at moment t of `run`, as it is there: its
/// force and moment in the joint's frame, N and N mm (`(reaction PATH)`). None where the run carries no
/// reaction of that path; Err where it is refused.
pub fn reaction(run: &command::Read, path: &str, t: f64) -> Option<Result<([f64; 3], [f64; 3]), String>> {
    let index = run.model.reactions.iter().position(|r| r.path == path)?;
    let m = run.model.at(t);
    Some(match run.model.reaction_refused_lossless(index, &m) {
        Some(why) => Err(why.why),
        None => Ok(m.reactions[index].local()),
    })
}

/// A placement's 12 numbers (R row major, then T) at `ctx`: its params' values and its moment.
pub fn placement_at(doc: &Doc, x: &ink_lang::Xf, ctx: &Ctx) -> [f64; 12] {
    std::array::from_fn(|k| Num::of(&doc.g, if k < 9 { x.r[k / 3][k % 3] } else { x.t[k - 9] }).at(ctx))
}

/// `reaction`, about a frame placed in the joint's frame (`local`, R row major then T) and in its axes: a
/// frame a part on the joint declares (`(reaction PATH (frame-of PATH NAME))`).
pub fn reaction_in(run: &command::Read, path: &str, t: f64, local: &[f64; 12]) -> Option<Result<([f64; 3], [f64; 3]), String>> {
    let index = run.model.reactions.iter().position(|r| r.path == path)?;
    let m = run.model.at(t);
    Some(match run.model.reaction_refused_lossless(index, &m) {
        Some(why) => Err(why.why),
        None => Ok(m.reactions[index].in_local(local)),
    })
}

/// What a reaction carrying a routed belt's part says of its friction.
const UNPAIRED: &str = "a routed belt's loss is friction at the parts power enters it by (DYNAMICS.md, Couplings), in the reactions carrying them; \
     what reacts it through the belt's other parts' supports, shares the belt's bending sets, is not found and not in them";
