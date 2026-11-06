//! Interference: every pair of a target's instances whose boxes meet,
//! claimed disjoint outside their declared contact zones (contacts, and
//! relations with a zone), to `ACCURACY`, each answered as a `disjoint`
//! claim is. Still pairs are proven in parallel, each question once
//! (pairs placed alike, in a frame of their own: pairs::framed), and
//! answered by the claims' own step (claims::answer); pairs that move over
//! the motion are checked as `disjoint` claims are (claims::ask), their
//! searches side by side.
//! Results are typed; clients name and present them (`ink interfere`).

use crate::InkError;
use crate::claims::{self, Check, ClaimOptions, Status};
use crate::motion::{self, Asked, Found, Search};
use crate::pairs::{self, Framed};
use crate::view::ModelView;
use crate::{Claim, Inexact, Item, PlacedRelation};
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::tape::Tape;
use ink_kernel::query::{self, Verdict, max_extent};
use ink_kernel::ranged::{self, Proof, ProofOptions};
use std::collections::HashMap;
use std::sync::Arc;
use std::time::{Duration, Instant};

/// The accuracy every pair is asked to, mm, wherever the model lies:
/// twice the claims' default, for its cost (`WHY`). An overlap shallower
/// than parts are made to changes no decision. A stopgap: no part states
/// how exact it is, to take the tighter of a pair's.
pub const ACCURACY: f64 = 0.1;

/// What interfere's accuracy lets pass, as its answer says. Why 0.1, twice
/// a claim's 0.05: a contact costs as its area over the accuracy squared,
/// so a pair at 0.05 costs about 4 times as much (a 193-pair assembly: 19
/// CPU-s at 0.1, 80 at 0.05). The 0.1 is a depth inside both parts, half
/// an overlap's thickness.
pub const WHY: &str = "an overlap up to 0.2 mm thick passes (a shaft 0.4 mm over its bore's diameter). Check a fit with a claim: (disjoint shaft bore) with its (accuracy x)";

/// The pairs to check, found (`plan`): what the answer says besides its
/// pairs, and each pair's question, built.
pub struct Plan {
    /// Items resting on a solve with no value here, with its error: no
    /// pair with one is checked (as assert's claims: ModelView::refusal).
    pub unsolved_items: Vec<(String, InkError)>,
    /// Items not proven bounded: checked against every other item.
    pub unbounded_items: Vec<String>,
    /// Pairs apart by their boxes, proven (to the accuracy).
    pub apart: usize,
    /// Pairs apart by their boxes, these and those of which one reads
    /// inexact inputs (an estimate).
    pub boxed_apart: usize,
    /// The contacts and relations under the target, which pairs name by
    /// index (Pair::zones, Pair::relations).
    pub relations: Vec<PlacedRelation>,
    /// What is asked at one moment though its question reads t
    /// (claims::reads_t: its items' fields and zones), no span covering
    /// it: each pair so asked, and each such instance in none of them
    /// (its pairs apart by their boxes there), with that moment.
    pub one_moment: Vec<(String, f64)>,
    items: Vec<Item>,
    candidates: Vec<(usize, usize)>,
    regions: Vec<(V3, V3)>,
    facts: Vec<Facts>,
    zones: Vec<Option<u32>>,
    unsolved: Vec<Option<String>>,
    riding: Vec<bool>,
    questions: Vec<Question>,
    asks: Vec<(usize, usize, Option<Framed>)>,
}

/// A pair's answer: its items, its contact zones and relations (indices
/// into Plan::relations), the inexact inputs it reads, whether its
/// intersection is proven bounded, and what checking it found.
pub struct Pair {
    pub items: [String; 2],
    /// The relations whose zones its proof leaves out.
    pub zones: Vec<usize>,
    /// The relations between its items (a kind says what relate says).
    pub relations: Vec<usize>,
    /// The inexact inputs it reads (sampled fields, solved values): then
    /// it is an estimate.
    pub inexact: Inexact,
    /// Its intersection is proven bounded. Else it was examined in the
    /// window only: apart there is not apart (`status`).
    pub whole: bool,
    pub answer: Answer,
    /// Wall-clock milliseconds it took: its own proof's (a still pair), or
    /// its own wall time, built, searched and finished (a moving one).
    pub ms: Option<f64>,
    /// Answered from the cache: a still pair whose question another pair
    /// asked first (Some(true); the first, None); a moving pair none of
    /// whose kernel answers was computed.
    pub cached: Option<bool>,
    /// The one moment it was asked at though its question reads t
    /// (Plan::one_moment).
    pub one_moment: Option<f64>,
}

/// What checking a pair found.
pub enum Answer {
    /// A still pair's claim, its proof proven beside the others' (as a
    /// `disjoint` claim answers it: claims::answer), with its overlap's
    /// volume where it overlaps.
    Still { check: Check, volume: Option<Volume> },
    /// A moving pair's claim, checked over the motion.
    Moving { check: Check },
    /// Not checked: the call's time limit had passed.
    NotBegun,
    /// Not checked: a contact zone of its rests on a solve with no value
    /// here, as the error says.
    Unsolved(String),
}

/// An overlap's volume, estimated; measured in a pair's frame, scaled by
/// its |det R|.
pub struct Volume {
    pub estimate: f64,
}

impl Pair {
    /// Whether it was found overlapping.
    pub fn overlaps(&self) -> bool {
        self.status() == Status::Fails
    }

    /// Its status: its claim's.
    pub fn status(&self) -> Status {
        match &self.answer {
            Answer::Still { check, .. } | Answer::Moving { check, .. } => check.status,
            Answer::NotBegun | Answer::Unsolved(_) => Status::Unknown,
        }
    }
}

/// What a pair is besides its answer (Pair's fields).
struct Facts {
    items: [String; 2],
    zones: Vec<usize>,
    relations: Vec<usize>,
    inexact: Inexact,
    whole: bool,
    one_moment: Option<f64>,
    /// A still pair's field as placed, with its zones (Plan::world_question).
    world: Option<u32>,
}

impl Facts {
    fn pair(&self, answer: Answer, ms: Option<f64>, cached: Option<bool>) -> Pair {
        let one_moment = self.one_moment;
        Pair {
            items: self.items.clone(),
            zones: self.zones.clone(),
            relations: self.relations.clone(),
            inexact: self.inexact.clone(),
            whole: self.whole,
            answer,
            ms,
            cached,
            one_moment,
        }
    }
}

/// A still pair's question (`prove_pairs`): its field clipped to its
/// region, and whether it is asked in a shared frame (pairs::framed).
struct Question {
    node: u32,
    region: (V3, V3),
    framed: bool,
}

/// f (an intersection, or b of a pair a, b) where `zone` is not: what must
/// be empty, or what a must not meet, for a pair apart outside its zone.
fn outside(v: &mut ModelView, f: u32, zone: Option<u32>) -> u32 {
    match zone {
        None => f,
        Some(z) => {
            let out = v.doc.g.neg(z);
            v.doc.g.max(f, out)
        }
    }
}

/// What interfere checks without a target: root, else every body, each
/// unplaced at its own origin, so refused once an assembly places a body.
pub fn default_items(v: &mut ModelView) -> Result<Vec<Item>, InkError> {
    let d = &v.doc;
    if d.body("root").is_none() && d.assembly("root").is_none() && d.assemblies.iter().any(|a| a.instances.iter().any(|i| d.body(&i.child).is_some())) {
        return Err(InkError::msg("the file has no root and its bodies sit unplaced: name an assembly, --target NAME, or write (assembly root ...) placing the parts"));
    }
    Ok(v.doc.default_items()?)
}

/// The pairs of `all` (a target's items, `target` what it places, where
/// it names a body or assembly: Doc::target) to check, at `res`, and each
/// still pair's question.
pub fn plan(v: &mut ModelView, all: Vec<Item>, target: Option<&crate::Target>, res: Option<f64>) -> Result<Plan, InkError> {
    if all.len() < 2 {
        return Err(InkError::msg("interfere needs a target with at least two instances"));
    }
    // Every question over its ranged params' whole range, as the claims'
    // entry points: a pointwise view (a library caller's) would box each
    // item at one value and call pairs apart that meet at another.
    v.scope_question()?;
    // As assert's claims (ModelView::refusal): an item resting on a solve
    // with no value here is left out, with the solve's error, and the rest
    // are checked. Over a motion, an item that moves with t is solved at
    // each step. An item's field is placed only where read: over a motion,
    // and for the pairs whose boxes meet asked in the world.
    let stepped = |v: &ModelView, f: u32| v.motion.is_some_and(|span| claims::moves(v, &[f], span));
    let (mut items, mut unsolved_items, mut item_moves) = (Vec::new(), Vec::new(), Vec::new());
    for it in all {
        let moves = v.motion.is_some() && {
            let f = v.doc.placed_field(&it);
            stepped(v, f)
        };
        let reads = v.doc.placed_reads(std::slice::from_ref(&it));
        match v.refusal(&reads).filter(|_| !moves) {
            Some(e) => unsolved_items.push((it.label.clone(), e)),
            None => {
                items.push(it);
                item_moves.push(moves);
            }
        }
    }
    // Contacts and relations: a relation's zone is a contact's. Those
    // enclosing assemblies declare between its parts are the target's too.
    let relations = target.map_or_else(Vec::new, |t| v.doc.relations_of(t));
    // A relation's end covers an item if it names it or one of its parents.
    let covers = |end: &str, label: &str| label == end || label.starts_with(&format!("{end}/"));
    let relating = |p: &str, q: &str| -> Vec<usize> {
        (0..relations.len())
            .filter(|&k| {
                let (a, b) = (&relations[k].a, &relations[k].b);
                (covers(a, p) && covers(b, q)) || (covers(a, q) && covers(b, p))
            })
            .collect()
    };
    let zone = |k: usize| relations[k].zone;
    // Items that read inexact inputs (Doc::inexact: sampled fields, solved
    // values): every pair with one is measured, and so is its box. What
    // reads one of a union's parts reads one, so each field is asked once.
    let certified = v.certified();
    let mut reads_inexact: HashMap<u32, bool> = HashMap::new();
    let measured: Vec<bool> = items
        .iter()
        .map(|it| v.doc.placed_reads(std::slice::from_ref(it)).into_iter().any(|f| *reads_inexact.entry(f).or_insert_with(|| v.doc.inexact(&[f], &certified).any())))
        .collect();
    // Over a motion, boxes holding each item at every moment of it.
    let (mut bounds, mut whole, mut over) = (Vec::with_capacity(items.len()), Vec::with_capacity(items.len()), Vec::with_capacity(items.len()));
    for (b, w, o) in motion::item_boxes(v, &items) {
        bounds.push(b);
        whole.push(w);
        over.push(o);
    }
    // Outside ink's range, no pair is answered, by boxes or otherwise: its
    // parts' numbers and the target's contact zones', over the span asked.
    let mut fields = v.doc.placed_reads(&items);
    fields.extend(relations.iter().filter_map(|r| r.zone));
    let scope = crate::validity::scope(v);
    crate::validity::refuse(&v.doc, &scope, &fields, bounds.iter().flatten().map(|b| (b.0, b.1)), res.unwrap_or(ACCURACY).min(ACCURACY))?;
    let exact: Vec<bool> = measured.iter().map(|m| !m).collect();
    let (candidates, boxed_apart, among) = pairs::candidate_pairs(&bounds, &whole, &over, &[&exact]);
    drop(over);
    let apart = among[0];
    // The reach of each pair neither of whose items is bounded (pair_region),
    // found for all of them at once.
    let open: Vec<u32> = candidates
        .iter()
        .filter(|&&(i, j)| !whole[i] && !whole[j])
        .map(|&(i, j)| {
            let (a, b) = (v.doc.placed_field(&items[i]), v.doc.placed_field(&items[j]));
            v.doc.g.max(a, b)
        })
        .collect();
    v.reach_all(&open);
    let mut regions = Vec::with_capacity(candidates.len());
    let mut pair_whole = Vec::with_capacity(candidates.len());
    for &c in &candidates {
        let (r, w) = pair_region(v, &items, c, &bounds, &whole);
        regions.push(r);
        pair_whole.push(w);
    }
    // Narrow phase: each pair's field as a graph node now (cheap); it is
    // compiled inside the parallel proof and focused on the pair's region at
    // once, so a pair against a huge union (a screen of 10 000 holes) holds
    // only the parts near it, and only one full tape per thread exists.
    // Overlap inside the pair's contact zones is intended: it is left out.
    let pair_relations: Vec<Vec<usize>> = candidates.iter().map(|&(i, j)| relating(&items[i].label, &items[j].label)).collect();
    let pair_zones: Vec<Vec<usize>> = pair_relations.iter().map(|rs| rs.iter().copied().filter(|&k| zone(k).is_some()).collect()).collect();
    // A pair's inputs are what its items' fields read and its zones: the
    // refusal and measured rules apply to them all (a zone may read t or a
    // solve).
    let mut pair_inputs: Vec<Vec<u32>> = Vec::with_capacity(candidates.len());
    let mut zones: Vec<Option<u32>> = Vec::with_capacity(candidates.len());
    // A pair whose zone rests on a solve with no value here: unknown, with
    // the solve's error.
    let mut unsolved: Vec<Option<String>> = vec![None; candidates.len()];
    // Pairs that move over the motion: an item or a zone of theirs does.
    // Of those, the pairs whose items ride one rigid group (the same moving
    // instance, nothing below it moving: motion::rigid_group) and whose
    // zones do not move keep their placement relative to each other: a still
    // question, asked once in that instance's frame.
    let groups: Vec<Option<usize>> = match v.motion {
        Some(span) => items.iter().map(|it| motion::rigid_group(v, it, span)).collect(),
        None => vec![Some(0); items.len()],
    };
    let rides = |i: usize, j: usize| match (groups[i], groups[j]) {
        (Some(m), Some(n)) => m > 0 && m == n && items[i].links[m - 1].path == items[j].links[m - 1].path,
        _ => false,
    };
    let mut moving: Vec<bool> = Vec::with_capacity(candidates.len());
    let mut riding: Vec<bool> = Vec::with_capacity(candidates.len());
    for (k, (&(i, j), zs)) in candidates.iter().zip(&pair_zones).enumerate() {
        let mut inputs = v.doc.placed_reads(&items[i..=i]);
        inputs.extend(v.doc.placed_reads(&items[j..=j]));
        moving.push(item_moves[i] || item_moves[j]);
        riding.push(moving[k] && rides(i, j));
        let z = zs.iter().filter_map(|&k| zone(k)).reduce(|a, b| v.doc.g.min(a, b));
        if let Some(z) = z {
            let zone_moves = stepped(v, z);
            moving[k] |= zone_moves;
            riding[k] &= !zone_moves;
            if let Some(e) = v.refusal(&[z]).filter(|_| !zone_moves) {
                let (a, b) = (&relations[zs[0]].a, &relations[zs[0]].b);
                unsolved[k] = Some(format!("contact zone of '{a}' and '{b}': {}", e.message));
            }
            inputs.push(z);
        }
        zones.push(z);
        pair_inputs.push(inputs);
    }
    // Each pair's own inexact inputs: what its fields read is what each
    // reads, so a pair of items and a zone reading none reads none.
    // Asked at one moment though the question reads t: no span covers
    // it and the run asked for none. Looked for only where the model
    // reads t at all; a pair reads it where one of its parts or its zones
    // does.
    let at = v.unspanned_moment().filter(|_| claims::model_reads_t(v));
    let pair_one_moment: Vec<Option<f64>> = pair_inputs.iter().map(|inputs| at.and_then(|_| claims::one_moment_of(v, inputs))).collect();
    let mut one_moment: Vec<(String, f64)> = Vec::new();
    for (k, &(i, j)) in candidates.iter().enumerate() {
        if let Some(t) = pair_one_moment[k] {
            one_moment.push((format!("pair {} & {}", items[i].label, items[j].label), t));
        }
    }
    if let Some(t) = at {
        for (i, it) in items.iter().enumerate() {
            let named = candidates.iter().enumerate().any(|(k, &(a, b))| (a == i || b == i) && pair_one_moment[k].is_some());
            let reads = v.doc.placed_reads(std::slice::from_ref(it));
            if !named && claims::one_moment_of(v, &reads).is_some() {
                one_moment.push((format!("instance {}", it.label), t));
            }
        }
    }
    let zone_inexact: Vec<bool> = zones.iter().map(|z| z.is_some_and(|z| *reads_inexact.entry(z).or_insert_with(|| v.doc.inexact(&[z], &certified).any()))).collect();
    let facts: Vec<Facts> = (0..candidates.len())
        .map(|k| {
            let (i, j) = candidates[k];
            let inexact = match measured[i] || measured[j] || zone_inexact[k] {
                true => v.doc.inexact(&pair_inputs[k], &certified),
                false => Inexact::default(),
            };
            Facts {
                items: [items[i].label.clone(), items[j].label.clone()],
                zones: pair_zones[k].clone(),
                relations: pair_relations[k].clone(),
                inexact,
                whole: pair_whole[k],
                one_moment: pair_one_moment[k],
                world: None,
            }
        })
        .collect();
    // Pairs that move over the motion are checked as a claim over it is
    // (motion: t ranged), their searches side by side; the rest, and
    // those riding one group, are proven at once, in parallel.
    let still: Vec<usize> = (0..candidates.len()).filter(|&k| (!moving[k] || riding[k]) && unsolved[k].is_none()).collect();
    let mut plan = Plan {
        unsolved_items,
        unbounded_items: items.iter().zip(&whole).filter(|(_, w)| !**w).map(|(it, _)| it.label.clone()).collect(),
        apart,
        boxed_apart,
        relations,
        one_moment,
        items,
        candidates,
        regions,
        facts,
        zones,
        unsolved,
        riding,
        questions: Vec::new(),
        asks: Vec::with_capacity(still.len()),
    };
    // Each still pair's question: asked in a frame of its own (pairs::framed:
    // one of its parts', else the deepest instance holding it), where its
    // fields are bodies placed rigidly and its zones are in its ends'
    // coordinates, so pairs placed alike (every copy of a sub-assembly, a
    // screw on each of its nuts, on a carriage or not) are one question; else
    // where it lies, its fields placed.
    let mut asked: HashMap<u32, usize> = HashMap::new();
    for &k in &still {
        let world = plan.world_question(v, k);
        plan.facts[k].world = Some(world.node);
        let (i, j) = plan.candidates[k];
        let parts: Option<Vec<(u32, &str)>> =
            plan.facts[k].zones.iter().map(|&z| plan.relations[z].local.filter(|l| l.b.is_none()).map(|l| (l.a, plan.relations[z].a.as_str()))).collect();
        let frame = parts.filter(|_| plan.facts[k].whole && !plan.facts[k].inexact.any()).and_then(|p| pairs::framed(v, &plan.items[i], &plan.items[j], &p, ACCURACY));
        let Some(fr) = frame else {
            // Riding a group, it is asked in its frame or over t: the
            // world's question holds at one moment.
            if plan.riding[k] {
                continue;
            }
            plan.questions.push(world);
            plan.asks.push((k, plan.questions.len() - 1, None));
            continue;
        };
        let (fa, fb) = (fr.items[0].0, fr.items[1].0);
        let f = v.doc.g.max(fa, fb);
        let f = outside(v, f, fr.zone);
        let node = v.doc.clip(f, fr.region.0, fr.region.1);
        // Riding a moving group, it is proven once only where its whole
        // question in the group's frame reads no t: a zone still in the
        // world moves in that frame (it counter-moves), and its parts may
        // too. Else over t, as one with no frame.
        if plan.riding[k] && claims::reads_t(v, &[node]) {
            continue;
        }
        let q = *asked.entry(node).or_insert_with(|| {
            plan.questions.push(Question { node, region: fr.region, framed: true });
            plan.questions.len() - 1
        });
        plan.asks.push((k, q, Some(fr)));
    }
    Ok(plan)
}

/// Where a pair of items can meet, and whether that holds all of their
/// intersection: inside the intersection of their boxes (each field is
/// positive outside its own proven bounds), so proofs start there, clipped
/// to it. For a pair of items neither of which is bounded, their
/// intersection's own reach: it may be bounded (or empty) all the same; if
/// not, the window, and only results that need all of it say so.
fn pair_region(v: &mut ModelView, items: &[Item], (i, j): (usize, usize), bounds: &[Option<(V3, V3, f64)>], whole: &[bool]) -> ((V3, V3), bool) {
    let (mut lo, mut hi) = ([f64::NEG_INFINITY; 3], [f64::INFINITY; 3]);
    for k in [i, j].into_iter().filter(|&k| whole[k]) {
        let (a, b, _) = bounds[k].unwrap();
        lo = std::array::from_fn(|d| lo[d].max(a[d]));
        hi = std::array::from_fn(|d| hi[d].min(b[d]));
    }
    let mut w = whole[i] || whole[j];
    if !w {
        let (a, b) = (v.doc.placed_field(&items[i]), v.doc.placed_field(&items[j]));
        let f = v.doc.g.max(a, b);
        let r = v.reach(f);
        (lo, hi, w) = (r.lo, r.hi, r.whole);
    }
    let pad = 1e-6 * max_extent(lo, hi).max(1e-9);
    ((lo.map(|v| v - pad), hi.map(|v| v + pad)), w)
}

/// A pair's row of a caller's answer, reported with it as it is decided
/// (`report`).
pub type Row<'a> = dyn Fn(&Pair) -> Option<serde_json::Value> + Sync + 'a;

/// Pair p, decided, reported (claims::progress) with its row, so a call
/// stopped at its time limit gives the rows of the pairs done.
fn report(p: Pair, row: &Row) -> Pair {
    let mut done = serde_json::json!({"done": p.items, "status": p.status().as_str(), "ms": p.ms});
    if let Some(r) = row(&p) {
        done["row"] = r;
    }
    claims::progress(done);
    p
}

/// What a moving pair's search sends back: the pair, what it found (or its
/// panic), its wall time, and whether the cache did not answer it.
type Sent = (usize, std::thread::Result<Found>, Duration, bool);

/// Pair k's search run on a runner (ink_kernel::runners), beside the view
/// and the other pairs', what it found sent back, a panic too: the view,
/// which keeps a sender itself, would otherwise wait for it forever.
fn beside(k: usize, search: impl FnOnce() -> Found + Send + 'static, tx: &std::sync::mpsc::Sender<Sent>) {
    let tx = tx.clone();
    ink_kernel::runners::spawn(move || {
        let start = Instant::now();
        let found = std::panic::catch_unwind(std::panic::AssertUnwindSafe(search));
        _ = tx.send((k, found, start.elapsed(), !ink_kernel::cache::here().1));
    });
}

/// What `ink interfere` asks at one moment though it reads t, in a view
/// scoped as interfere scopes it: its default target's plan and each
/// assembly's no other places, as `plan` finds them (Plan::one_moment:
/// the pairs it would check, and its moving instances apart by their
/// boxes). A target interfere refuses asks nothing. Named, for `ink check`.
pub fn at_one_moment(v: &mut ModelView) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    if !claims::model_reads_t(v) || v.motion_over(v.doc.spans_hull()).is_err() || v.unspanned_moment().is_none() {
        return out;
    }
    // Placed, or a copy of it built with arguments is.
    let placed: std::collections::HashSet<String> =
        v.doc.assemblies.iter().flat_map(|a| a.instances.iter()).flat_map(|i| [Some(i.child.clone()), v.doc.assembly(&i.child).and_then(|c| c.of.clone())]).flatten().collect();
    let tops: Vec<String> = v.doc.assemblies.iter().filter(|a| a.of.is_none() && !placed.contains(&a.name) && a.name != "root").map(|a| a.name.clone()).collect();
    for name in std::iter::once(None).chain(tops.into_iter().map(Some)) {
        let items = match &name {
            None => default_items(v),
            Some(n) => v.doc.items(n).map_err(InkError::from),
        };
        let target = v.doc.target(name.as_deref().unwrap_or("root"));
        let Ok(items) = items else { continue };
        // Planned only where its parts' fields or its zones read t: else
        // no pair of it, nor part, is asked at one moment.
        let mut reads = v.doc.placed_reads(&items);
        reads.extend(target.as_ref().map_or_else(Vec::new, |t| v.doc.relations_of(t)).iter().filter_map(|r| r.zone));
        if claims::one_moment_of(v, &reads).is_none() {
            continue;
        }
        let Ok(p) = plan(v, items, target.as_ref(), None) else { continue };
        for (what, _) in p.one_moment {
            if !out.contains(&what) {
                out.push(what);
            }
        }
    }
    out
}

impl Plan {
    /// Pair k's question where it lies, its fields placed.
    fn world_question(&self, v: &mut ModelView, k: usize) -> Question {
        let (i, j) = self.candidates[k];
        let (fa, fb) = (v.doc.placed_field(&self.items[i]), v.doc.placed_field(&self.items[j]));
        let f = v.doc.g.max(fa, fb);
        let f = outside(v, f, self.zones[k]);
        let (lo, hi) = self.regions[k];
        Question { node: v.doc.clip(f, lo, hi), region: self.regions[k], framed: false }
    }

    /// Every pair checked, in order, at cell size `res` (None: automatic). Each is reported as it is decided
    /// (claims::progress), with its row of the caller's answer, `row`, so
    /// a call stopped at its time limit gives the rows of the pairs done.
    pub fn run(self, v: &mut ModelView, res: Option<f64>, row: &Row) -> Result<Vec<Pair>, InkError> {
        claims::progress(serde_json::json!({"total": self.candidates.len(), "of": "pairs", "kept": true}));
        let n = self.candidates.len();
        let mut rows: Vec<Option<Pair>> = (0..n).map(|_| None).collect();
        let mut again = Vec::new();
        let opts = ClaimOptions { res, accuracy: Some(ACCURACY), ..Default::default() };
        for (k, pair) in prove_pairs(v, &self.questions, &self.asks, &self.facts, opts, row) {
            match pair {
                Some(p) => rows[k] = Some(p),
                None => again.push(k),
            }
        }
        // What its frame left to the world (a witness not shown inside
        // there): where it lies; riding a group, over t below.
        again.retain(|&k| !self.riding[k]);
        again.sort_unstable();
        let questions: Vec<Question> = again.iter().map(|&k| self.world_question(v, k)).collect();
        let asks: Vec<(usize, usize, Option<Framed>)> = again.iter().enumerate().map(|(q, &k)| (k, q, None)).collect();
        for (k, pair) in prove_pairs(v, &questions, &asks, &self.facts, opts, row) {
            rows[k] = pair;
        }
        // Moving pairs, in order: each one's question is built here, in the
        // graph, one after another; its search (motion::Search: tapes
        // alone) runs beside the others' while the next is built, and what it
        // found is finished here as it comes back, which may ask another. The
        // searches share the runners, their pieces in turn (ink_kernel::runners). A
        // search whose question another in flight asks (its key in the answer
        // cache: every copy of a sub-assembly) waits for that one, then is
        // answered from the cache. At most twice the pool's threads are in
        // flight or waiting, each holding its tapes. A pair's ms is its own
        // wall time (built, searched sharing the threads, finished); cached,
        // that none of it was computed.
        // Moving pairs run beside each other (ClaimOptions::beside).
        let opts = ClaimOptions { beside: true, ..opts };
        struct Waiting {
            then: Box<dyn FnOnce(&mut ModelView, Found) -> Result<Asked<Answer>, InkError>>,
            key: Option<u128>,
            spent: Duration,
            computed: bool,
        }
        let (tx, rx) = std::sync::mpsc::channel::<Sent>();
        let mut waiting: HashMap<usize, Waiting> = HashMap::new();
        // The keys of the searches in flight, each with the pairs whose search
        // asks the same, held back until it returns.
        let mut flying: HashMap<u128, Vec<(usize, Search)>> = HashMap::new();
        let cap = 2 * rayon::current_num_threads();
        let moving: Vec<usize> = (0..n).filter(|&k| self.unsolved[k].is_none() && rows[k].is_none()).collect();
        let mut queue = moving.into_iter();
        // Searches that may replace a pair's answer (motion::Upgrade: t
        // ranged after its sampled answer), begun once every pair has an
        // answer, as runners free, until the call's time limit: what each
        // stands in for (its pair, check, joint's word, time).
        let mut deferred: std::collections::VecDeque<(usize, motion::Upgrade, Check, Duration, bool)> = Default::default();
        let misses = || ink_kernel::cache::here().0;
        loop {
            let next = if waiting.len() < cap { queue.next() } else { None };
            let (k, asked, spent, computed) = match next {
                // Past the call's time limit the rest are not begun, as
                // assert's claims (a pair's setup alone, its shapes and joint,
                // is milliseconds: tens of thousands of them would run ten
                // minutes on).
                Some(k) if ink_kernel::stop::reached() => {
                    rows[k] = Some(report(self.facts[k].pair(Answer::NotBegun, None, None), row));
                    continue;
                }
                Some(k) => {
                    let (start, before) = (Instant::now(), misses());
                    claims::progress(serde_json::json!({"running": self.facts[k].items}));
                    let (i, j) = self.candidates[k];
                    let (a, b) = (v.doc.placed_field(&self.items[i]), v.doc.placed_field(&self.items[j]));
                    let asked = over_motion(v, (a, b), self.zones[k], opts)?;
                    (k, asked, start.elapsed(), misses() > before)
                }
                // Every pair answered: an upgrade begun as a runner frees.
                None if waiting.len() < cap && !deferred.is_empty() && !ink_kernel::stop::reached() => {
                    let (k, up, check, spent, computed) = deferred.pop_front().expect("an upgrade");
                    let start = Instant::now();
                    match (up.start)(v)? {
                        motion::Later::Now(None) => continue,
                        motion::Later::Now(Some(check)) => {
                            let ms = (spent + start.elapsed()).as_secs_f64() * 1e3;
                            rows[k] = Some(report(self.facts[k].pair(Answer::Moving { check }, Some(ms), Some(false)), row));
                            continue;
                        }
                        // Its search waits as a first answer's does: on one
                        // in flight that asks the same (every copy of a
                        // sub-assembly), then from the cache.
                        motion::Later::Wait(search, then) => {
                            let then: Box<dyn FnOnce(&mut ModelView, Found) -> Result<Asked<Answer>, InkError>> =
                                Box::new(move |v, found| Ok(Asked::Done(Answer::Moving { check: then(v, found)?.unwrap_or(check) })));
                            (k, Asked::Wait(search, then), spent + start.elapsed(), computed)
                        }
                    }
                }
                None if waiting.is_empty() => break,
                None => {
                    let (k, found, took, computed) = rx.recv().expect("a search sends what it found");
                    // A search that panicked panics here, as it would searched in turn.
                    let found = found.unwrap_or_else(|p| std::panic::resume_unwind(p));
                    let w = waiting.remove(&k).expect("a pair waits on its search");
                    for (h, search) in w.key.and_then(|key| flying.remove(&key)).unwrap_or_default() {
                        beside(h, move || search.run(), &tx);
                    }
                    let (start, before) = (Instant::now(), misses());
                    let asked = (w.then)(v, found)?;
                    (k, asked, w.spent + took + start.elapsed(), w.computed || computed || misses() > before)
                }
            };
            match asked {
                Asked::Done(answer) => {
                    for up in std::mem::take(&mut v.upgrades) {
                        if let Answer::Moving { check } = &answer {
                            deferred.push_back((k, up, check.clone(), spent, computed));
                        }
                    }
                    rows[k] = Some(report(self.facts[k].pair(answer, Some(spent.as_secs_f64() * 1e3), Some(!computed)), row));
                }
                Asked::Wait(search, then) => {
                    let key = search.key();
                    waiting.insert(k, Waiting { then, key, spent, computed });
                    match key.map(|key| flying.entry(key)) {
                        Some(std::collections::hash_map::Entry::Occupied(mut held)) => held.get_mut().push((k, search)),
                        held => {
                            if let Some(std::collections::hash_map::Entry::Vacant(e)) = held {
                                e.insert(Vec::new());
                            }
                            beside(k, move || search.run(), &tx);
                        }
                    }
                }
            }
        }
        Ok(rows
            .into_iter()
            .enumerate()
            .map(|(k, p)| match &self.unsolved[k] {
                Some(why) => report(self.facts[k].pair(Answer::Unsolved(why.clone()), None, None), row),
                None => p.expect("every pair is checked or not begun"),
            })
            .collect())
    }
}

/// The proofs of still pairs (`asks`: each pair's index into `facts`; its
/// question's into `questions`; its frame), each question proven once, in
/// parallel, and each of its pairs answered as a `disjoint` claim to the
/// accuracy is (claims::answer) and reported as decided (`row`), with the
/// overlap's volume; None where the pair is to be asked in the world. A
/// framed question is asked to the accuracy less the frame's share, at the
/// cells the accuracy reaches, and moved to each pair's world (`framed`).
/// The first pair asking a question gives its proof's time; the others are
/// `cached`.
fn prove_pairs(v: &ModelView, questions: &[Question], asks: &[(usize, usize, Option<Framed>)], facts: &[Facts], opts: ClaimOptions, row: &Row) -> Vec<(usize, Option<Pair>)> {
    use rayon::prelude::*;
    let ctx = v.ctx.clone();
    let (res, tol) = (opts.res, opts.accuracy.unwrap_or(ACCURACY));
    let world = v.proof_options(res, tol);
    let shared = ProofOptions { tol: Framed::asked(tol), ..world };
    let mut of: Vec<Vec<usize>> = vec![Vec::new(); questions.len()];
    for (n, &(_, q, _)) in asks.iter().enumerate() {
        of[q].push(n);
    }
    // Each question's ranged params (a witness names these). A disjoint
    // claim's accuracy is a length.
    let reads: Vec<Vec<(u32, (f64, f64))>> = questions.iter().map(|q| v.ranged_reads(&[q.node])).collect();
    let quantity = Some(crate::units::Quantity::Length);
    let g = &v.doc.g;
    // Each pair's held unsaid where it reads ranged params or t as placed.
    let unsaid: Vec<bool> = facts.iter().map(|f| f.world.is_none_or(|w| claims::held_unsaid(v, w))).collect();
    let why: Vec<Option<String>> =
        asks.iter().map(|(_, q, fr)| crate::validity::outside_region(&v.doc, &ctx, questions[*q].region, fr.as_ref().map(|fr| &fr.at), tol, res)).collect();
    questions
        .par_iter()
        .zip(&of)
        .zip(&reads)
        .flat_map_iter(|((q, of), read)| {
            let (lo, hi) = q.region;
            let asked = Instant::now();
            // The proof focuses it on the region (when not cached); the
            // volume's descent prunes it cell by cell.
            let tape = Arc::new(Tape::new(g, q.node));
            let outside = |k: usize, why: &String| {
                let unknown = Check { accuracy: None, ..claims::checked(claims::Status::Unknown, claims::Evidence::Outside(why.clone())) };
                Answer::Still { check: claims::answer(unknown, tol, facts[k].inexact.clone(), quantity, read.clone()), volume: None }
            };
            if let Some(w) = of.iter().map(|&a| why[a].as_ref()).collect::<Option<Vec<_>>>() {
                return of
                    .iter()
                    .zip(w)
                    .map(|(&a, w)| {
                        let k = asks[a].0;
                        (k, Some(report(facts[k].pair(outside(k, w), None, None), row)))
                    })
                    .collect::<Vec<_>>();
            }
            let p = ranged::prove(&tape, &ctx, lo, hi, if q.framed { &shared } else { &world });
            let proof_ms = asked.elapsed();
            // Overlap volume for colliding pairs, measured once a pair gives
            // it (a frame's witness not shown in the world gives none: the
            // world measures its own).
            let measured = std::cell::OnceCell::new();
            let volume = || *measured.get_or_init(|| p.disproven().then(|| overlap_volume(&tape, &ctx, q.region)).flatten());
            // How it held: unsaid where the pair as placed reads ranged
            // params or t (`unsaid`); outright only where a frame reads back
            // exactly; its near point found in the question (once), moved to
            // each pair's world and confirmed on the pair's world field.
            let candidates = std::cell::OnceCell::new();
            let check = |c: Proof, k: usize, fr: Option<&Framed>| {
                let held = (!unsaid[k]).then(|| match fr {
                    Some(fr) if p.exact && fr.error > 0.0 => claims::Held::ToAccuracy { near: None },
                    _ => {
                        let points: Vec<V3> = match p.exact {
                            true => Vec::new(),
                            false => candidates.get_or_init(|| claims::Held::candidates(&p, &tape, &ctx)).iter().map(|&c| fr.map_or(c, |fr| fr.world(c).0)).collect(),
                        };
                        let world = match (fr, facts[k].world) {
                            (Some(_), Some(w)) if !points.is_empty() => Arc::new(Tape::new(g, w)),
                            _ => tape.clone(),
                        };
                        claims::Held::confirmed(p.exact, points, &world, &ctx)
                    }
                });
                Check { held, ..claims::answer(claims::emptied(c, facts[k].whole), tol, facts[k].inexact.clone(), quantity, read.clone()) }
            };
            let mut first = true;
            of.iter()
                .map(|&a| {
                    let (k, _, frame) = &asks[a];
                    let start = Instant::now();
                    let answer = match frame {
                        _ if let Some(w) = &why[a] => Some(outside(*k, w)),
                        None => Some(Answer::Still { check: check(p.clone(), *k, None), volume: volume().map(|estimate| Volume { estimate }) }),
                        Some(fr) => framed(&p, fr, &tape, &ctx, tol).map(|proof| {
                            // Its volume is the frame's, measured.
                            let volume = volume().map(|x| Volume { estimate: x * fr.det.mid() });
                            Answer::Still { check: check(proof, *k, Some(fr)), volume }
                        }),
                    };
                    let pair = answer.map(|answer| {
                        let ms = if first { proof_ms } else { start.elapsed() };
                        let cached = (!std::mem::take(&mut first)).then_some(true);
                        report(facts[*k].pair(answer, Some(ms.as_secs_f64() * 1e3), cached), row)
                    });
                    (*k, pair)
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// A question's proof in a pair's frame `fr` (pairs::framed), moved to the
/// world: holds where it holds; fails where the cube that puts the cube of
/// half-side the accuracy about the witness's world point inside the
/// world's pair (`Framed::cube`) is shown inside the pair in the frame, at
/// the witness's params; unknown with its undecided cell's box moved.
/// None: a witness not so shown, which the world is asked.
fn framed(p: &Proof, fr: &Framed, tape: &Tape, ctx: &Ctx, tol: f64) -> Option<Proof> {
    use ink_core::interval::I;
    let verdict = match &p.verdict {
        Verdict::Proven => Verdict::Proven,
        Verdict::Disproven { witness, t, value, params } => {
            let (w, off) = fr.world(*witness);
            let h = fr.cube(tol, off);
            let mut at = ctx.clone();
            for &(k, x) in params.iter().filter(|(k, _)| ctx.p_spans.get(*k).is_some() && !ctx.p_spans.is_derived(*k)) {
                at = at.with_param(k, (x, x));
            }
            let cube = std::array::from_fn(|d| I::new((witness[d] - h).next_down(), (witness[d] + h).next_up()));
            if !(ink_core::eval::eval_interval(tape, &at, cube, None).hi < 0.0) {
                return None;
            }
            Verdict::Disproven { witness: w, t: *t, value: *value, params: params.clone() }
        }
        Verdict::Unknown { cells, volume, ended, cell, .. } => {
            let cell = cell.map(|[lo, hi]| {
                let (a, b) = fr.at.aabb(lo, hi);
                [a, b]
            });
            // Its attained accuracy is the frame's, which the frame's
            // rounding moves: none is given.
            Verdict::Unknown { cells: *cells, volume: volume * fr.det.mid(), ended: *ended, cell, attained: None }
        }
    };
    Some(Proof { verdict, res: p.res, tol, exact: p.exact, loose: Vec::new() })
}

/// A moving pair checked over the motion as a claim (claims::ask), a
/// apart from b outside its zones.
fn over_motion(v: &mut ModelView, (a, b): (u32, u32), zone: Option<u32>, opts: ClaimOptions) -> Result<Asked<Answer>, InkError> {
    let b = outside(v, b, zone);
    claims::ask(v, Claim::Disjoint(a, b), None, opts)?.then(v, move |_, check| Ok(Asked::Done(Answer::Moving { check })))
}

/// The volume of a pair's intersection (its tape) in its region,
/// estimated as `ink mass` measures one (query::measured). None: it has no
/// box.
fn overlap_volume(tape: &Arc<Tape>, ctx: &Ctx, (lo, hi): (V3, V3)) -> Option<f64> {
    let (blo, bhi, _) = query::bounds_quick(tape, ctx, lo, hi)?;
    Some(query::measured(tape, &[None], false, ctx, blo, bhi, None).sums[0])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A moving pair's search that panics sends its panic back, for the
    /// view to raise: the view keeps a sender itself, so it would wait for
    /// the search forever (interfere hanging past its time limit).
    #[test]
    fn a_search_that_panics_is_sent_back() {
        let (tx, rx) = std::sync::mpsc::channel::<Sent>();
        beside(7, || panic!("a search's bug"), &tx);
        let (k, found, ..) = rx.recv_timeout(Duration::from_secs(60)).expect("sent back");
        assert!(k == 7 && found.is_err());
    }

    /// A library caller's view is pointwise: plan scopes it over the
    /// ranged params itself, so two cubes apart at u = 0 that meet at
    /// u = 10 are not called apart by their boxes at one value.
    #[test]
    fn a_pointwise_view_is_planned_over_its_ranges() {
        static CORE: crate::Registry = crate::Registry::new(&[]);
        let src = "(import \"std\")\n(def u (ranged 0 10))\n(body cube (shape (box 2 2 2)))\n\
                   (assembly root (instance a cube) (instance b cube (at (move (- 10 u) 0 0))))";
        let doc = crate::Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
        let (ctx, _) = doc.ctx_at(&[], 0.0).unwrap();
        let mut v = ModelView::new(doc, ctx);
        let (items, target) = (v.doc.items("root").unwrap(), v.doc.target("root"));
        let plan = plan(&mut v, items, target.as_ref(), None).unwrap();
        assert_eq!((plan.apart, plan.candidates.len()), (0, 1));
        let pairs = plan.run(&mut v, None, &|_| None).unwrap();
        assert!(pairs[0].overlaps());
    }
}
