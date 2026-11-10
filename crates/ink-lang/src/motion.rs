//! Claims over a motion (`(t-span A B)`, `--t-span`). Time is a range of
//! the question, as a tolerance band is: a claim that reads t is the same
//! claim with t a ranged param over the span (`ranged`), decided by the
//! one search over every moment at once, so what holds or fails is
//! guaranteed as it is over any band. What that search leaves undecided
//! (within the call's time limit), or cannot ask (a solve on t not
//! followed, a package's metric, an integral), is measured at evenly
//! spaced moments (`measured`): holding there is an estimate, failing at
//! one is a witness.

use crate::InkError;
use crate::claims::{self, Check, ClaimOptions, Evidence, Status, checked};
use crate::view::ModelView;
use crate::{Claim, Cmp, Xf};
use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, Node, UnOp, Var};
use ink_core::interval;
use ink_core::tape::Tape;
pub use ink_kernel::motion::{MAX_STEPS, Motion};
use ink_kernel::rigid::Isometry;

/// How a claim over a motion was checked.
#[derive(Clone, Debug)]
pub struct Stepped {
    pub motion: Motion,
    /// Its answer is a proof (t ranged over the span, or asked once);
    /// else it is measured at the steps' middles.
    pub bound: bool,
    /// t was a ranged param of the one search (`ranged`).
    pub ranged: bool,
    /// The pair rides one moving instance: asked once (`riding`).
    pub riding: bool,
    /// The uniform step a measured claim was checked at (its configured
    /// size, however many steps were checked before a disproof).
    pub step: Option<f64>,
    /// Measured at steps a provider of what it reads chose, and how
    /// (ParamProvider::steps), not at uniform ones.
    pub basis: Option<String>,
    /// Measured, why t ranged did not decide it (`ranged`), and where it
    /// searched and was left undecided, its own unknown answer (its hint,
    /// undecided cell, the params it ranged).
    pub left: Option<String>,
    pub left_answer: Option<Box<Check>>,
}

impl Stepped {
    fn of(motion: Motion, bound: bool, step: Option<f64>) -> Stepped {
        Stepped { motion, bound, ranged: false, riding: false, step, basis: None, left: None, left_answer: None }
    }
}

/// A kernel search a check waits on, built: it reads its tapes and its
/// context alone, no view, so the searches of many checks run side by
/// side (`Asked`): an emptiness with t ranged (motion::ranged), as the
/// tape, its context, the region and the proof's options.
pub struct Search(std::sync::Arc<Tape>, Ctx, (V3, V3), ink_kernel::ranged::ProofOptions);

/// What may replace an answer given already (a moving pair's sampled one:
/// `over_span`), begun once every answer is given, with the time left.
pub struct Upgrade {
    pub start: Box<dyn FnOnce(&mut ModelView) -> Result<Later, InkError>>,
}

/// What an Upgrade's start gives: the answer replacing the one given (None:
/// it stands), or the search that decides it and what it then gives.
pub enum Later {
    Now(Option<Check>),
    Wait(Search, Box<dyn FnOnce(&mut ModelView, Found) -> Result<Option<Check>, InkError>>),
}

/// What a search found.
pub type Found = ink_kernel::ranged::Proof;

impl Search {
    /// Its question's key in the answer cache (None: not cached): the one
    /// ranged::prove stores its proof under, so a search another in flight
    /// asks (every copy of a sub-assembly) waits for that one's proof.
    pub fn key(&self) -> Option<u128> {
        let Search(tape, ctx, (lo, hi), opts) = self;
        ink_kernel::ranged::prove_key(tape, None, ctx, *lo, *hi, opts)
    }

    pub fn run(&self) -> Found {
        let Search(tape, ctx, (lo, hi), opts) = self;
        ink_kernel::ranged::prove(tape, ctx, *lo, *hi, opts)
    }
}

/// What a check that may wait on a search gives: its answer, or the
/// search and what is done with what it finds (with the view again: a
/// clash's moment, an exact question, another search). Its question is
/// built in the view (which adds to the graph), one check after another;
/// its search needs no view, so a caller with many (ink interfere's moving
/// pairs) runs their searches side by side while it builds the next.
pub enum Asked<T> {
    Done(T),
    Wait(Search, Box<dyn FnOnce(&mut ModelView, Found) -> Result<Asked<T>, InkError>>),
}

impl<T: 'static> Asked<T> {
    /// Its answer, its searches run here, one after another.
    pub fn wait(mut self, v: &mut ModelView) -> Result<T, InkError> {
        loop {
            match self {
                Asked::Done(x) => return Ok(x),
                Asked::Wait(s, then) => self = then(v, s.run())?,
            }
        }
    }

    /// `f` of its answer, once it has one.
    pub fn then<U: 'static>(self, v: &mut ModelView, f: impl FnOnce(&mut ModelView, T) -> Result<Asked<U>, InkError> + 'static) -> Result<Asked<U>, InkError> {
        match self {
            Asked::Done(x) => f(v, x),
            Asked::Wait(s, then) => Ok(Asked::Wait(s, Box::new(move |v, found| then(v, found)?.then(v, f)))),
        }
    }
}

/// An item's box over the view's motion and whether it is proven bounded
/// (as ModelView::item_bounds and item_whole at a moment): a box holding
/// the item at every moment, for an interfere's broad phase, and where it
/// moves rigidly, its boxes over pieces of the span (`swept`: in order of
/// t, covering it), so that pairs whose boxes meet at no moment are apart
/// (pairs::candidate_pairs).
pub type ItemBox = (Option<(V3, V3, f64)>, bool, Vec<((f64, f64), (V3, V3))>);

/// Pieces of the span an item that moves rigidly is boxed over (`swept`).
const SWEPT_PIECES: usize = 64;

/// Each of `items`' ItemBox. A placed body proven bounded in its frame
/// (over the span, where its own field reads t) has its boxes over
/// SWEPT_PIECES uniform pieces of the span (`swept`), and their hull; any
/// other that moves has its proven bounds with t ranged over the whole span
/// (the one search halves t where it blurs them). One that reads a sampled
/// or solved value has no box: a candidate against every other item.
pub fn item_boxes(v: &mut ModelView, items: &[crate::Item]) -> Vec<ItemBox> {
    items
        .iter()
        .map(|it| {
            let moving = v.motion.map(|span| (span, v.doc.placed_field(it))).filter(|&(span, f)| claims::moves(v, &[f], span));
            let Some((span, f)) = moving else { return (v.item_bounds(it), v.item_whole(it), Vec::new()) };
            if v.doc.inexact(&[f], &[]).any() {
                return (None, false, Vec::new());
            }
            let over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..v.ctx.clone() };
            // Its own box: at rest, or over the span where its own field
            // reads t (a belt's teeth running round its loop).
            let local = |v: &mut ModelView, s: u32| v.at_point(over.clone(), |v| v.whole(s).then(|| v.bounds(s)).flatten());
            // Moved forward, R q + T, only by an isometry: the field reads
            // R^T (p - T), whose solid R^-T moves, not R.
            if let Some(p) = v.doc.placed.get(&f).map(|p| (p.shape, p.at))
                && let Some(iso) = Isometry::over(&v.doc.g, &p.1.r, &p.1.t, &over)
                && let Some((lo, hi, _)) = local(v, p.0)
                && let Some(pieces) = swept(v, &iso, (lo, hi), span, true)
            {
                let (lo, hi) = pieces.iter().map(|p| p.1).reduce(ink_core::v3::box_union).expect("pieces");
                return (Some((lo, hi, 0.0)), true, pieces);
            }
            let (b, w) = v.at_point(over, |v| if v.whole(f) { (v.bounds(f), true) } else { (None, false) });
            (b, w, Vec::new())
        })
        .collect()
}

/// Boxes holding a body (its box (lo, hi) in its own frame) placed by
/// `at` (an isometry over the span: so its solid is where R p + T puts it)
/// over SWEPT_PIECES uniform pieces of `span`, each holding it at every
/// moment of its piece: R p + T with R's and T's entries bounded over the
/// piece. None: an entry not bounded.
fn swept(v: &ModelView, at: &Isometry, (lo, hi): (V3, V3), span: (f64, f64), refine: bool) -> Option<Vec<((f64, f64), (V3, V3))>> {
    use interval::I;
    let entry = |n: u32| Tape::new(&v.doc.g, n);
    let r: Vec<Tape> = at.r().iter().flatten().map(|&n| entry(n)).collect();
    let t: Vec<Tape> = at.t().iter().map(|&n| entry(n)).collect();
    let p: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
    let o = [I::point(0.0); 3];
    let size = |b: &(V3, V3)| (0..3).map(|k| b.1[k] - b.0[k]).fold(0.0, f64::max);
    let own = size(&(lo, hi));
    let boxed = |s: (f64, f64)| -> Option<(V3, V3)> {
        let c = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..v.ctx.clone() };
        let mut b = ([0.0; 3], [0.0; 3]);
        for j in 0..3 {
            let mut x = ink_core::eval::eval_interval(&t[j], &c, o, None);
            for k in 0..3 {
                x = interval::bin(BinOp::Add, x, interval::bin(BinOp::Mul, ink_core::eval::eval_interval(&r[3 * j + k], &c, o, None), p[k]).0).0;
            }
            if !(x.lo.is_finite() && x.hi.is_finite()) {
                return None;
            }
            (b.0[j], b.1[j]) = (x.lo, x.hi);
        }
        Some(b)
    };
    // SWEPT_PIECES uniform pieces, each (`refine`) halved while its box
    // reaches past the body's own size by more than BLUR_SHARE of it, down
    // to a SWEPT_FINEST-th of the span.
    let n = SWEPT_PIECES;
    let finest = (span.1 - span.0) / SWEPT_FINEST as f64;
    let mut todo: Vec<(f64, f64)> = (0..n)
        .rev()
        .map(|i| (span.0 + (span.1 - span.0) * i as f64 / n as f64, if i + 1 == n { span.1 } else { span.0 + (span.1 - span.0) * (i + 1) as f64 / n as f64 }))
        .collect();
    let mut out = Vec::new();
    while let Some(s) = todo.pop() {
        let b = boxed(s)?;
        let m = ink_core::eval::span_mid(s);
        if refine && size(&b) > own * (1.0 + BLUR_SHARE) && s.1 - s.0 > finest && s.0 < m && m < s.1 {
            todo.extend([(m, s.1), (s.0, m)]);
            continue;
        }
        out.push((s, b));
    }
    Some(out)
}

/// The share of a body's size its box over a piece of t may reach past
/// it before the piece is halved (`swept`).
const BLUR_SHARE: f64 = 1.0 / 8.0;

/// The finest piece of the span `swept` halves to: this share of it.
const SWEPT_FINEST: usize = 1 << 10;

/// Whether two items' boxes over pieces of a motion (`swept`, in order of
/// t, covering the span; empty: its box holds it throughout) meet at no
/// moment: each pair of pieces that share a moment, in one pass.
pub fn apart_over(a: &[((f64, f64), (V3, V3))], b: &[((f64, f64), (V3, V3))]) -> bool {
    let meet = |x: &(V3, V3), y: &(V3, V3)| (0..3).all(|k| x.0[k] <= y.1[k] && y.0[k] <= x.1[k]);
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let ((sa, ba), (sb, bb)) = (&a[i], &b[j]);
        if sb.0 <= sa.1 && sa.0 <= sb.1 && meet(ba, bb) {
            return false;
        }
        // The piece ending first is done with (closed: at a tie, the other
        // is still met against the next, which shares its end).
        if sa.1 <= sb.1 { i += 1 } else { j += 1 }
    }
    true
}

/// The claim at every moment of `span` (`bound`: its bound, as a field, if
/// it reads params): asked once where the pair rides one moving instance
/// (`riding`), else with t ranged (`ranged`), else measured.
/// The view is left at its moment.
pub(crate) fn over_span(v: &mut ModelView, claim: Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Asked<Check>, InkError> {
    let t0 = v.ctx.t;
    let bound = bound.filter(|&b| claims::reads_t(v, &[b]));
    Asked::Done(()).then(v, move |v, ()| {
        // Asked once only where the claim was judged over every moment
        // (claims::judged_scope): one reading a solve only measured is
        // judged moment by moment, measured.
        let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
        if bound.is_none() && claims::judged_scope(v, &roots, Some(span)).2 && riding(v, &claim, span) {
            // Asked once, at the view's moment: their relative placement is
            // the same at every moment.
            let motion = v.motion.take();
            let c = claims::decide(v, claim, opts);
            v.motion = motion;
            let m = Motion::none(if c.as_ref().is_ok_and(|c| c.status == Status::Unknown) { ink_kernel::ranged::Status::Unknown } else { ink_kernel::ranged::Status::Proven }, span);
            let c = c?;
            return Ok(Asked::Done(Check { note: c.note.or(Some(RIDING)), motion: Some(Stepped { riding: true, ..Stepped::of(m, true, None) }), ..c }));
        }
        // Where a drive it reads may jump, a claim is unknown there; one
        // reading a solve on t is measured, the jump said beside its
        // estimate (the solve's root is followed only where its equations
        // are shown continuous).
        let solved = v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t);
        let jump = jumps(v, &claim, bound, span);
        if let (Some(c), false) = (&jump, solved) {
            return Ok(Asked::Done(c.clone()));
        }
        // The sampled answer first: cheap, always there, an estimate (a
        // failure at a moment is a witness, guaranteed: nothing to prove,
        // but where it reads a solve on t, found there to its solver's
        // accuracy: a witness with t ranged is confirmed over its root's
        // enclosure).
        let sampled = measured(v, claim.clone(), bound, span, opts);
        v.set_moment(t0)?;
        if sampled.status == Status::Fails && sampled.accuracy.is_some() && !solved {
            return Ok(Asked::Done(sampled));
        }
        let keep = |sampled: Check, why: String, answer: Option<Check>| Check { motion: sampled.motion.map(|s| Stepped { left: Some(why), left_answer: answer.map(Box::new), ..s }), ..sampled };
        let jumped = |why: String| match &jump {
            Some(Check { evidence: Evidence::NotCovered(j), .. }) => format!("{why}, and {j}, so between the moments measured the speed may pass its limits"),
            _ => why,
        };
        if jump.is_some() {
            return Ok(Asked::Done(keep(sampled, jumped("t ranged is not tried".into()), None)));
        }
        // Beside other searches (interfere's moving pairs): the sampled
        // answer is the pair's, and t ranged is tried once every pair has
        // one, with the time left (`Upgrade`), replacing it where it decides.
        if opts.beside && v.doc.claim_roots(&claim).iter().any(|&f| own_field_moves(v, f)) {
            return Ok(Asked::Done(keep(sampled, "a shape's own field reads t (a belt's teeth running round, a chain's links). Beside other pairs, finding its region over the span is a search of its own, so t ranged is not tried".into(), None)));
        }
        if opts.beside {
            let asked = claim.clone();
            v.upgrades.push(Upgrade {
                start: Box::new(move |v| {
                    let r = ranged(v, &asked, bound, span, opts);
                    v.set_moment(t0)?;
                    Ok(match r? {
                        Ranged::Decided(c) => Later::Now(Some(claims::finish_later(v, &asked, bound, c, opts))),
                        Ranged::Left(..) => Later::Now(None),
                        Ranged::Wait(search, whole, read, world, moments, charged) => {
                            let asked_in = (search.0.clone(), search.1.clone());
                            let then = move |v: &mut ModelView, p: Found| {
                                // Over a drive's value only a hold is an
                                // answer: else the sampled one stands.
                                if !moments {
                                    return Ok((p.status() == ink_kernel::ranged::Status::Proven).then(|| {
                                        let held = held_beside(v, &p, &asked_in, None).map(|h| match h.charged(charged) {
                                            claims::Held::Outright => claims::Held::Outright,
                                            claims::Held::ToAccuracy { .. } => claims::Held::ToAccuracy { near: None },
                                            h @ claims::Held::SetAside { .. } => h,
                                        });
                                        let c = Check { ranged: read, held, note: Some(DRIVE), ..claims::in_window(Status::Holds, whole, Evidence::Empty(p)) };
                                        match decided(c, span, false) {
                                            Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                                            _ => None,
                                        }
                                    }).flatten());
                                }
                                let held = held_beside(v, &p, &asked_in, world.as_ref()).map(|h| h.charged(charged));
                                // A witness in the frame shown none in the
                                // world: the sampled answer stands.
                                let e = match &world {
                                    Some(w) if p.status() != ink_kernel::ranged::Status::Unknown => match w.confirmed(v, Evidence::Empty(p.clone())) {
                                        Some(e) => e,
                                        None => return Ok(None),
                                    },
                                    _ => Evidence::Empty(p.clone()),
                                };
                                let c = Check { ranged: read, held, ..claims::in_window(Status::of(p.status()), whole, e) };
                                Ok(match decided(c, span, false) {
                                    Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                                    _ => None,
                                })
                            };
                            Later::Wait(search, Box::new(then))
                        }
                    })
                }),
            });
            return Ok(Asked::Done(keep(sampled, "t ranged is tried once every pair has an answer, with the time left: it did not decide it".into(), None)));
        }
        // Then t ranged, with the time left, to prove it or find a witness;
        // undecided, the sampled answer stands, saying why.
        let r = ranged(v, &claim, bound, span, opts);
        v.set_moment(t0)?;
        Ok(Asked::Done(match r? {
            Ranged::Decided(c) => c,
            Ranged::Left(why, answer) => keep(sampled, why, answer),
            Ranged::Wait(..) => unreachable!("only beside other searches is t ranged searched apart"),
        }))
    })
}

const RIDING: &str =
    "both shapes ride one moving instance, with nothing below it moving. They keep the same relative placement at every moment, so the pair was asked once, at the view's moment";

/// Whether a claim about two shapes (disjoint, a clearance, touches,
/// subset, overlaps: each the same under a rigid motion of both) is about
/// two items of the root that ride one rigid group over `span`
/// (`rigid_group`: the same moving instance, nothing below it moving, and
/// every placement down to it proven a rotation).
fn riding(v: &mut ModelView, claim: &Claim, span: (f64, f64)) -> bool {
    let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else { return false };
    let Ok(items) = v.doc.items("root") else { return false };
    let mut placed = |want: u32| items.iter().position(|i| v.doc.placed_field(i) == want);
    let (Some(ia), Some(ib)) = (placed(a), placed(b)) else { return false };
    let (ia, ib) = (&items[ia], &items[ib]);
    match (rigid_group(v, ia, span), rigid_group(v, ib, span)) {
        (Some(m), Some(n)) => m > 0 && m == n && ia.links[m - 1].path == ib.links[m - 1].path,
        _ => false,
    }
}

/// A speed reads its drive's derivative along t (Graph::derivative's
/// mark), which the formula takes where the drive jumps too: a jump is
/// unbounded speed, so where a drive may jump (`may_jump`) a claim that
/// could be proven over t is unknown, with the moment (its sampled values
/// say nothing there); one only ever measured says so beside its estimate.
fn jumps(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64)) -> Option<Check> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    let marked: Vec<u32> = v.doc.g.reachable(&roots).into_iter().filter_map(|n| if let Node::Deriv(_, d, 3) = v.doc.g.nodes[n as usize] { Some(d) } else { None }).collect();
    let s = marked.into_iter().find_map(|d| may_jump(&Tape::new(&v.doc.g, d), &v.ctx, span))?;
    let why = format!("it reads the speed of a drive whose rate has no bound for t in [{}, {}]. It may jump there, and a jump is unbounded speed", s.0, s.1);
    let m = Motion { undecided: vec![s], t: Some(ink_core::eval::span_mid(s)), ..Motion::none(ink_kernel::ranged::Status::Unknown, span) };
    Some(Check { motion: Some(Stepped { ranged: true, ..Stepped::of(m, true, None) }), ..checked(Status::Unknown, Evidence::NotCovered(why)) })
}

/// Whether a shape's own field (not its placement) reads t: a placed
/// body's, or any part of a union or intersection of them.
fn own_field_moves(v: &ModelView, f: u32) -> bool {
    match v.doc.placement_of(f) {
        Some(p) => claims::reads_t(v, &[p.shape]),
        None => match v.doc.g.nodes[f as usize] {
            Node::Bin(BinOp::Min | BinOp::Max, a, b) => own_field_moves(v, a) || own_field_moves(v, b),
            // Outside a contact zone (interfere's max(b, -zone)).
            Node::Tag(_, a) | Node::Un(UnOp::Neg, a) => own_field_moves(v, a),
            _ => claims::reads_t(v, &[f]),
        },
    }
}

/// What `ranged` settled: an answer over every moment, or why the sampled
/// one stands.
enum Ranged {
    Decided(Check),
    /// Why, and the search's own unknown answer where it searched.
    Left(String, Option<Check>),
    /// A disjoint pair's emptiness to search beside others, whether its
    /// shapes are proven bounded, the ranged params it reads, the frame it
    /// is searched in (`Frame`), whether t's param is a moment (not a
    /// drive's angle: `Drive`, where only a hold is an answer), and what
    /// the shapes asked may be off the world's (Held::charged).
    Wait(Search, bool, Vec<(u32, (f64, f64))>, Option<Frame>, bool, f64),
}

/// The answer with t ranged as param k: holding or failing, the claim's
/// over every moment (its witness's moment, param k's value); undecided,
/// the sampled answer stands (`own`: the time share stopped it).
fn decided(c: Check, span: (f64, f64), own: bool) -> Ranged {
    match c.status {
        Status::Unknown => Ranged::Left(
            match own {
                true => "with t ranged over the span the search was not decided within half the call's time left".into(),
                false => "with t ranged over the span the search left it undecided".into(),
            },
            Some(c),
        ),
        status => {
            let k = c.ranged.iter().map(|r| r.0).max();
            let t = k.and_then(|k| witness_t(&c.evidence, k));
            let m = Motion { t, ..Motion::none(if status == Status::Holds { ink_kernel::ranged::Status::Proven } else { ink_kernel::ranged::Status::Disproven }, span) };
            Ranged::Decided(Check { motion: Some(Stepped { ranged: true, ..Stepped::of(m, true, None) }), ..c })
        }
    }
}

/// The claim with t a ranged param over `span` (the view's param "t",
/// made once): every node of it, and its bound, with t replaced by that
/// param (Graph::subst), decided as any claim over ranged params is
/// (claims::decide), the param's value a witness's moment. Only t's span is
/// added to the view's scope: what the view knows of fields that do not read it
/// holds. What a solve on t sets is derived from that param
/// (solved_over_t). Not ranged: a claim of a kind no band ranges (a
/// package's metric, an integral, a thickness, cavities), one reading a
/// solve on t not followed over the span, one whose bound moves.
/// Asked alone (assert), it stops at half the call's time left (stop::within);
/// beside other searches (interfere), at the call's limit.
fn ranged(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Ranged, InkError> {
    // Through its drive where it has one (`Drive`); not decided there
    // but by a hold, asked over t. A witness in the frame that is none in
    // the world: asked there.
    let (mut framed, mut through) = (true, true);
    loop {
        match ranged_in(v, claim, bound, span, opts, framed, through)? {
            Some(r) => return Ok(r),
            None if through => through = false,
            None if framed => framed = false,
            None => unreachable!("a world search over t needs no confirming"),
        }
    }
}

/// `ranged`, in the first shape's frame where `framed` and it has one
/// (`moving_frame`), through its drive where `through` and it has one
/// (`drive`). None: decided there by a witness that is not shown one in
/// the world (`Frame::confirmed`), or through the drive other than by a
/// hold.
#[allow(clippy::too_many_arguments)]
fn ranged_in(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions, framed: bool, through: bool) -> Result<Option<Ranged>, InkError> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    // The whole scope first (`scope_over`), kept until the answer is given
    // in the world: every rewrite of the claim (its frame, simplified over
    // the scope) and every witness confirmed reads a solved input as
    // ranged, never as the number it is at the view's moment.
    let (scope, k) = match scope_over(v, &roots, span) {
        Ok(s) => s,
        Err(why) => return Ok(Some(Ranged::Left(why, None))),
    };
    let (ctx, motion) = (std::mem::replace(&mut v.ctx, scope), v.motion.take());
    let r = scoped(v, claim, bound, span, opts, (framed, through), k, &roots);
    (v.ctx, v.motion) = (ctx, motion);
    r
}

/// A claim's scope over a motion (`roots` its fields): the view's context
/// with t's param k (made once: ModelView::t_param) ranged over `span`, and
/// what each solve on t they read sets derived from it (solved_over_t:
/// enclosed over each piece of the span by the steps it is followed by).
/// It is what the claim is validated in (claims::ask_in), searched in and
/// confirmed in, in that order. Err: why a solve on t it reads is not
/// derived (it is then measured, moment by moment).
pub(crate) fn scope_over(v: &mut ModelView, roots: &[u32], span: (f64, f64)) -> Result<(Ctx, u32), String> {
    let k = v.t_param();
    let solved = solved_over_t(v, roots, k, span)?;
    let mut c = v.ctx.clone();
    c.params[k as usize] = ink_core::eval::span_mid(span);
    c.p_spans = c.p_spans.with_band(k, span);
    for (d, bands) in solved {
        c.p_spans = c.p_spans.derived(d, &bands);
    }
    c.derive();
    Ok((c, k))
}

/// `ranged_in` in its scope (the view's context, t's param k ranged over
/// `span` and what solves on t set derived from it): in its frame, then
/// through its drive, then searched, then confirmed in the world.
#[allow(clippy::too_many_arguments)]
fn scoped(
    v: &mut ModelView,
    claim: &Claim,
    bound: Option<u32>,
    span: (f64, f64),
    opts: ClaimOptions,
    (framed, through): (bool, bool),
    k: u32,
    roots: &[u32],
) -> Result<Option<Ranged>, InkError> {
    if let Some(why) = sweeps_far(v, roots, span) {
        return Ok(Some(Ranged::Left(why, None)));
    }
    let t = v.doc.g.var(Var::Param(k));
    // Searched where one of the two shapes rests (`frame`): the same claim
    // at every moment, its points moved rigidly.
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let frame = if framed { moving_frame(v, claim, span, Some(tol)) } else { None };
    let charge = frame.as_ref().map_or(0.0, |f| f.2);
    let tn = v.doc.g.var(Var::T);
    // Through the drive all its reads of t pass through, where one is
    // shown periodic and the claim in its frame still moves: t's param is
    // then the drive's value.
    let drive = match through {
        true => {
            let moves = roots.iter().any(|&n| {
                let n = frame.as_ref().map_or(n, |m| v.doc.g.replace(n, &m.0));
                v.doc.g.depends_on(n, |x| x == Var::T)
            });
            moves.then(|| drive(v, roots, span, k)).flatten()
        }
        false => None,
    };
    if through && drive.is_none() {
        return Ok(None);
    }
    let by: std::collections::HashMap<u32, u32> = [(drive.as_ref().map_or(tn, |d| d.node), t)].into();
    // Its witnesses are moments where t is the drive and nothing is
    // charged (`Drive::slack`); else only a hold is an answer.
    let moments = drive.as_ref().is_none_or(|d| d.node == tn && d.slack == 0.0);
    let world = frame.as_ref().filter(|_| moments).map(|f| Frame::new(v, claim, f.1, k, span, opts.accuracy, charge));
    let mut on = |n: u32| {
        let n = match &frame {
            Some(m) => v.doc.g.replace(n, &m.0),
            None => n,
        };
        v.doc.g.replace(n, &by)
    };
    let Some(asked) = with_nodes(claim, &mut on) else {
        return Ok(Some(Ranged::Left("a claim of this kind is not ranged".into(), None)));
    };
    let bound = bound.map(&mut on);
    let asked = match (asked, bound) {
        // A bound that moves: the value checked is the expression minus
        // the bound, against 0.
        (Claim::Value(f, cmp, _), Some(b)) => Claim::Value(v.doc.g.sub(f, b), cmp, 0.0),
        (_, Some(_)) => return Ok(Some(Ranged::Left("its bound moves (reads t), and is not ranged".into(), None))),
        (c, None) => c,
    };
    let read = v.doc.claim_roots(&asked);
    if drive.is_some() && v.doc.g.reachable(&read).contains(&tn) {
        return Ok(None);
    }
    // What the shapes asked may be off the world's, charged to the claim:
    // the drive's slack (Drive::slack) and the pivots' drift
    // (`moving_frame`). An emptiness decided to its accuracy less it, a
    // clearance over its distance plus it, an overlap's witness deeper by
    // it at every moment (not through a drive: its witness is no moment).
    let slack = drive.as_ref().map_or(0.0, |d| d.slack);
    let charged = interval::add_up(charge, slack);
    let (asked, opts) = match asked {
        _ if charged == 0.0 => (asked, opts),
        Claim::Empty(..) | Claim::Disjoint(..) | Claim::Subset(..) => {
            let less = -interval::add_up(-tol, charged);
            if !(less > 0.0) {
                return Ok(None);
            }
            (asked, ClaimOptions { accuracy: Some(less), ..opts })
        }
        Claim::Clearance(a, b, gap) => (Claim::Clearance(a, b, interval::add_up(gap, charged)), opts),
        Claim::Overlaps(..) if slack == 0.0 => (asked, ClaimOptions { accuracy: Some(interval::add_up(tol, charged)), ..opts }),
        _ => return Ok(None),
    };
    match &drive {
        // One period of the drive's values stands for every moment.
        Some(d) => {
            v.ctx.p_spans.narrow(k, d.window);
            v.ctx.params[k as usize] = ink_core::eval::span_mid(d.window);
        }
        None => {
            // Over one period of a motion longer than one
            // (ink_core::period): the rest repeats it.
            let tapes: Vec<Tape> = read.iter().map(|&n| Tape::new(&v.doc.g, n)).collect();
            if let Some(one) = ink_core::period::quotient(&tapes.iter().collect::<Vec<_>>(), &v.ctx) {
                v.ctx = one;
            }
        }
    }
    // t's param is reported over the span (through the drive its window
    // is the drive's values, not moments).
    let at_span = |mut r: Vec<(u32, (f64, f64))>| {
        r.iter_mut().filter(|x| x.0 == k && !moments).for_each(|x| x.1 = span);
        r
    };
    match (opts.beside, &asked) {
        (true, Claim::Disjoint(a, b)) => {
            let tol = claims::claim_accuracy(v, &asked, opts.accuracy);
            let f = v.doc.g.max(*a, *b);
            let tape = v.clipped(f);
            let crate::Reach { lo, hi, whole } = v.reach(f);
            Ok(Some(Ranged::Wait(Search(tape, v.ctx.clone(), (lo, hi), v.proof_options(opts.res, tol)), whole, at_span(v.ranged_reads(&read)), world, moments, charged)))
        }
        _ => {
            let (c, own) = ink_kernel::stop::within(if opts.beside { 1.0 } else { 0.5 }, || claims::decide(v, asked, opts));
            // t's param among them where the world reads it, though the
            // frame cancels it (a witness, a near point, is at a moment).
            let mut ranged = at_span(v.ranged_reads(&read));
            if (world.is_some() || !moments) && !ranged.iter().any(|r| r.0 == k) {
                ranged.push((k, span));
            }
            let c = c?;
            // An overlap fails only proven apart outright, which shapes off
            // the world's do not show of it: asked there. The answer is to
            // the claim's own accuracy, and holds outright only of shapes
            // asked as they are (Held::charged).
            if charged > 0.0 && matches!(claim, Claim::Overlaps(..)) && c.status == Status::Fails {
                return Ok(None);
            }
            let c = Check { accuracy: c.accuracy.map(|_| tol), held: c.held.map(|h| h.charged(charged)), ..c };
            // Through the drive, a hold is the answer (its near point is
            // at a drive's value, not a moment: none given); else asked
            // over t.
            let c = match (moments, c.status) {
                (true, _) => c,
                (false, Status::Holds) => Check {
                    held: c.held.map(|h| match h {
                        claims::Held::Outright => claims::Held::Outright,
                        claims::Held::ToAccuracy { .. } => claims::Held::ToAccuracy { near: None },
                        h @ claims::Held::SetAside { .. } => h,
                    }),
                    note: c.note.or(Some(DRIVE)),
                    ..c
                },
                (false, Status::Fails) => return Ok(None),
                (false, _) => c,
            };
            // Its witnesses in the world, each shown one there (undecided,
            // there are none).
            let c = match &world {
                Some(w) if c.status != Status::Unknown => match w.confirmed(v, c.evidence.clone()) {
                    Some(evidence) => Check { evidence, held: c.held.clone().map(|h| w.held(v, h)), ..c },
                    None => return Ok(None),
                },
                _ => c,
            };
            Ok(Some(decided(Check { ranged, ..c }, span, own)))
        }
    }
}

/// What each solve on t that `roots` read sets, derived from t ranged as
/// param k over `span` (follow::OverT), every one followed over all of it.
/// Err: why not: one not followed or bounded there (a dead point, more
/// equations than unknowns), resting on a value measured along t (a
/// simulation: Follow::certifies), or on a ranged param (it is followed at
/// one value of each).
fn solved_over_t(v: &ModelView, roots: &[u32], k: u32, span: (f64, f64)) -> Result<Vec<(std::sync::Arc<dyn ink_core::eval::Derive>, Vec<(u32, (f64, f64))>)>, String> {
    let on_t: Vec<usize> = v.doc.solves_read(roots).into_iter().filter(|&i| v.doc.solves[i].on_t).collect();
    if on_t.is_empty() {
        return Ok(Vec::new());
    }
    let follow = v.doc.follow();
    let mut out = Vec::new();
    for i in on_t {
        if let Some(&r) = follow.rests_on(i).iter().find(|&&r| v.ctx.p_spans.get(r as u32).is_some()) {
            return Err(format!("it reads {}, which rests on ranged '{}' and is followed along t at one value of it", follow.followed(i).2, v.doc.params[r].name));
        }
        if !follow.certifies(i) {
            return Err(format!("it reads {}, which rests on a value measured along t (a simulation's), not bounded", follow.followed(i).2));
        }
        let d = crate::follow::OverT::new(follow.clone(), i, k, &v.ctx.params, span)?;
        let bands = d.bands();
        out.push((std::sync::Arc::new(d) as std::sync::Arc<dyn ink_core::eval::Derive>, bands));
    }
    Ok(out)
}

/// How a search beside others (`Search`'s tape and context: `asked`)
/// held, as `shallow` says it: its near point found in the scope, moved
/// to the world where it was asked in a frame. None: not proven.
fn held_beside(v: &ModelView, p: &Found, (tape, ctx): &(std::sync::Arc<Tape>, Ctx), world: Option<&Frame>) -> Option<claims::Held> {
    if p.status() != ink_kernel::ranged::Status::Proven {
        return None;
    }
    let points = if p.exact { Vec::new() } else { claims::Held::scoped(p, tape, ctx) };
    let h = claims::Held::confirmed_at(p.exact, points, tape, ctx);
    Some(match world {
        Some(w) => w.held(v, h),
        None => h,
    })
}

const DRIVE: &str = "searched over one period of the drive every read of t passes through (each shape turned with it repeating by its own turn symmetry): that period's angles stand for every moment";

/// A claim over a motion read through its drive: `node`, a value of t
/// (and params) that every read of t the claim makes passes through (a
/// motor's angle spinning up as t^2), `window`, one period of the claim
/// in it from the least value it takes over the span, and `slack`, how far
/// (a length) a shape may sit from where the window places it.
struct Drive {
    node: u32,
    window: (f64, f64),
    slack: f64,
    /// How far the drive travels over the span (its tests weigh the
    /// window by it).
    #[cfg_attr(not(test), expect(dead_code))]
    travel: f64,
}

/// The claim's drive (`Drive`): the first of the values `tries` proposes
/// that `through` shows periodic. None: none shown, or the call's time
/// limit passed. It costs a reading of the claim per try.
fn drive(v: &mut ModelView, roots: &[u32], span: (f64, f64), k: u32) -> Option<Drive> {
    use ink_core::symmetry::{Order, turn_symmetry};
    let tries = tries(v, roots, k)?;
    // The shapes placed moving with t whose own fields do not, with their
    // turn symmetries: read once for every try.
    let g = &v.doc.g;
    let reached: std::collections::HashSet<u32> = g.reachable(roots).into_iter().collect();
    let moving = |n: u32| g.depends_on(n, |x| x == Var::T);
    let placed: Vec<(u32, crate::Placed)> = v
        .doc
        .placements()
        .into_iter()
        .filter(|(w, p)| reached.contains(w) && p.at.r.iter().flatten().chain(&p.at.t).any(|&n| moving(n)) && !g.depends_on(p.shape, |x| x == Var::T))
        .collect();
    let mut shapes = Vec::new();
    for (w, p) in placed {
        if ink_kernel::stop::reached() {
            return None;
        }
        let sym = turn_symmetry(&v.doc.g, p.shape);
        if sym.order != Order::Of(1) {
            shapes.push((w, p, sym));
        }
    }
    for d in tries {
        if ink_kernel::stop::reached() {
            return None;
        }
        if let Some(x) = through(v, roots, &shapes, d, span, k) {
            return (!ink_kernel::stop::reached()).then_some(x);
        }
    }
    None
}

/// The values `drive` reads the claim through, in turn, each once: t (a
/// witness's value is a moment), then of the values every read of t passes
/// through, those the walk down from the claim proposes, deepest first. One
/// per run of the walk (each but the nearest begins at a sin, cos or tan of
/// the value below it), found in a pass over the claim's readers of t.
/// None: the claim reads no drive.
fn tries(v: &ModelView, roots: &[u32], k: u32) -> Option<Vec<u32>> {
    let g = &v.doc.g;
    let tn = g.nodes.iter().position(|n| matches!(n, Node::Var(Var::T)))? as u32;
    let reach = g.reachable(roots);
    // What a solve on t sets is derived from t's param in the scope
    // (scope_over), not read through t's node: a claim reading one is not
    // a function of a drive found among t's readers, and narrowing t's
    // param to a drive's values would derive it from values that are no
    // moments. It is ranged over the span as it is (ink_core::period's
    // quotient leaves t whole where such a param rests on it).
    if reach.iter().any(|&n| matches!(g.nodes[n as usize], Node::Var(Var::Param(j)) if v.ctx.p_spans.is_derived(j) && v.ctx.p_spans.sources(j).contains(&k))) {
        return None;
    }
    // In node order (each after what it reads): what reads t, and what
    // reads a point. ts: the readers of t, in order; at: each one's place.
    let (mut on_t, mut on_p) = (std::collections::HashSet::new(), std::collections::HashSet::new());
    let mut ts: Vec<u32> = Vec::new();
    for &n in &reach {
        let (c, m) = g.nodes[n as usize].children();
        let node = g.nodes[n as usize];
        if n == tn || c[..m].iter().any(|a| on_t.contains(a)) {
            on_t.insert(n);
            ts.push(n);
        }
        if matches!(node, Node::Var(Var::X | Var::Y | Var::Z)) || c[..m].iter().any(|a| on_p.contains(a)) {
            on_p.insert(n);
        }
    }
    let at: std::collections::HashMap<u32, usize> = ts.iter().enumerate().map(|(i, &n)| (n, i)).collect();
    let first_root = roots.iter().filter_map(|r| at.get(r)).min().copied()?;
    // Every reader of t lies on a path from t to a root, so a reader d is
    // on every one exactly where no root comes before it and no edge
    // between readers passes over it (from before it to after it): one
    // sweep of the edges' spans.
    let mut over = vec![0i64; ts.len() + 1];
    for (i, &n) in ts.iter().enumerate() {
        let (c, m) = g.nodes[n as usize].children();
        for a in &c[..m] {
            if let Some(&j) = at.get(a)
                && j + 1 < i
            {
                over[j + 1] += 1;
                over[i] -= 1;
            }
        }
    }
    let mut covered = 0;
    let mut chain: Vec<u32> = Vec::new();
    for (i, &n) in ts.iter().enumerate().take(first_root + 1) {
        covered += over[i];
        // A value of no point (not a field): past one, every later reads it.
        if on_p.contains(&n) {
            break;
        }
        if covered == 0 {
            chain.push(n);
        }
    }
    // The walk proposes, `through` decides. Down the chain from its
    // nearest value, each stretch between two adjacent values is read once
    // (every reader in it, in node order) as a function of the lower
    // (`Read`, the period module's rules, its coefficients not computed):
    // the upper linear in it, repeating in it (sin and cos of it, however
    // their branches rejoin), or neither, where the walk stops (below, the
    // period module reads a use of the drive it cannot). Of a run of linear
    // stretches only the top is proposed: the claim is read through each
    // value of it alike, its coefficients scaled, with no less shown.
    #[derive(Clone, Copy, PartialEq)]
    enum Read {
        Lin,
        Kept,
        Other,
    }
    let read_of = |n: u32, read: &std::collections::HashMap<u32, Read>| {
        let (c, m) = g.nodes[n as usize].children();
        let arg = |a: u32| if on_t.contains(&a) { read.get(&a).copied().unwrap_or(Read::Other) } else { Read::Kept };
        let free = |a: u32| !on_t.contains(&a);
        let args: Vec<Read> = c[..m].iter().map(|&a| arg(a)).collect();
        if args.contains(&Read::Other) {
            return Read::Other;
        }
        if !args.contains(&Read::Lin) {
            return Read::Kept;
        }
        match g.nodes[n as usize] {
            Node::Un(UnOp::Neg | UnOp::Rad | UnOp::Deg, _) | Node::Tag(..) | Node::Real(_) => Read::Lin,
            Node::Un(UnOp::Sin | UnOp::Cos | UnOp::Tan, _) => Read::Kept,
            Node::Bin(BinOp::Add | BinOp::Sub, a, b) if (arg(a) == Read::Lin || free(a)) && (arg(b) == Read::Lin || free(b)) => Read::Lin,
            Node::Bin(BinOp::Mul, a, b) if free(a) || free(b) => Read::Lin,
            Node::Bin(BinOp::Div, _, b) if free(b) => Read::Lin,
            _ => Read::Other,
        }
    };
    let last = *chain.last()?;
    let mut read: std::collections::HashMap<u32, Read> = std::collections::HashMap::new();
    // keeps[j]: the stretch from chain[j] to chain[j + 1], read in chain[j].
    let mut keeps: Vec<Read> = Vec::new();
    let mut next = 1;
    for &n in &ts[at[&chain[0]]..=at[&last]] {
        if n == chain[0] {
            read.insert(n, Read::Lin);
            continue;
        }
        let r = read_of(n, &read);
        if chain.get(next) == Some(&n) {
            keeps.push(r);
            read.insert(n, Read::Lin);
            next += 1;
        } else {
            read.insert(n, r);
        }
    }
    // The tops of the runs the walk passes, nearest the claim first.
    let mut proposed = vec![last];
    for j in (0..keeps.len()).rev() {
        match keeps[j] {
            Read::Other => break,
            Read::Kept => proposed.push(chain[j]),
            Read::Lin => {}
        }
    }
    let mut tries = vec![tn];
    for d in proposed.into_iter().rev() {
        if !tries.contains(&d) {
            tries.push(d);
        }
    }
    Some(tries)
}

/// The claim read through d (`Drive`), with t's param k standing for d:
/// each placed shape (`shapes`: moving, n-fold by symmetry::turn_symmetry)
/// that turns with d about its own z (rigid::own_turn) stands for a value
/// repeating whenever d turns it by 2 pi / n (period::period_kept), so the
/// claim's period is the shapes' and the rest's together. A shape's turn
/// symmetry holds within its slack (radians of its fold), which moves its
/// points by at most that times twice their distance from its axis: the
/// drive's slack is their sum. None: no period shorter than d's travel.
fn through(v: &mut ModelView, roots: &[u32], shapes: &[(u32, crate::Placed, ink_core::symmetry::Symmetry)], d: u32, span: (f64, f64), k: u32) -> Option<Drive> {
    use ink_core::period::{Shown, form, period_kept};
    use ink_core::symmetry::Order;
    let tn = v.doc.g.var(Var::T);
    let a = v.doc.g.var(Var::Param(k));
    let by: std::collections::HashMap<u32, u32> = [(d, a)].into();
    // d's values over the span.
    let over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..v.ctx.clone() };
    let range = ink_core::eval::eval_interval_scope(&Tape::new(&v.doc.g, d), &over, [interval::I::point(0.0); 3], None);
    if !(range.lo.is_finite() && range.hi.is_finite() && range.lo < range.hi) {
        return None;
    }
    let mut ctx = v.ctx.clone();
    if ctx.params.len() <= k as usize {
        ctx.params.resize(k as usize + 1, f64::NAN);
    }
    ctx.params[k as usize] = ink_core::eval::span_mid((range.lo, range.hi));
    let ctx = ctx.with_band(k, (range.lo, range.hi));
    let mut kept_ctx = ctx.clone();
    let (mut map, mut kept, mut also, mut slack) = (by.clone(), Vec::new(), Vec::new(), 0.0);
    let is_a = |x: Var| x == Var::Param(k);
    for &(w, p, sym) in shapes {
        if ink_kernel::stop::reached() {
            return None;
        }
        let r = p.at.r.map(|row| row.map(|n| v.doc.g.replace(n, &by)));
        let g = &v.doc.g;
        let Some(wound) = ink_kernel::rigid::own_turn(g, &r, &ctx, |n| g.depends_on(n, is_a)) else { continue };
        // Its turn's rate: sum m_i c_i, each angle c_i a plus a constant.
        let rate = wound.iter().try_fold(None, |sum: Option<ink_core::period::Coef>, &(n, m)| match form(&Tape::new(g, n), Var::Param(k), &ctx)? {
            None => Some(sum),
            Some(c) => {
                let c = c.times(m as i128)?;
                match sum {
                    None => Some(Some(c)),
                    Some(s) => s.plus(c),
                }
            }
        });
        let Some(rate) = rate else { continue };
        let coef = match (sym.order, rate) {
            (Order::Of(n), Some(c)) => match c.times(n as i128) {
                Some(c) => Some(c),
                None => continue,
            },
            _ => None,
        };
        let moved = match sym.slack > 0.0 {
            false => 0.0,
            true => {
                let Some((lo, hi, _)) = v.whole(p.shape).then(|| v.bounds(p.shape)).flatten() else { continue };
                let r = [lo[0].abs().max(hi[0].abs()), lo[1].abs().max(hi[1].abs())];
                interval::mul_up(interval::mul_up(2.0, sym.slack), r[0].hypot(r[1]).next_up())
            }
        };
        let kd = kept_ctx.params.len().max(v.doc.params.len()) as u32;
        kept_ctx.params.resize(kd as usize + 1, 0.0);
        kept_ctx = kept_ctx.with_band(kd, (0.0, 0.0));
        let stand = v.doc.g.var(Var::Param(kd));
        map.insert(w, stand);
        kept.push((Var::Param(kd), coef));
        let t = p.at.t.map(|n| v.doc.g.replace(n, &by));
        also.extend(ink_kernel::rigid::simplified(&mut v.doc.g, &t, &ctx).unwrap_or(t.to_vec()));
        slack = interval::add_up(slack, moved);
    }
    let tops: Vec<u32> = roots.iter().map(|&n| v.doc.g.replace(n, &map)).chain(also).collect();
    if v.doc.g.reachable(&tops).contains(&tn) {
        return None;
    }
    let tapes: Vec<Tape> = tops.iter().map(|&n| Tape::new(&v.doc.g, n)).collect();
    let Shown::Base(b) = period_kept(&tapes.iter().collect::<Vec<_>>(), Var::Param(k), &kept_ctx, &kept)? else { return None };
    let p = b.shift().hi;
    (p > 0.0 && p.is_finite() && range.hi - range.lo > p).then(|| Drive { node: d, window: (range.lo, (range.lo + p).next_up()), slack, travel: range.hi - range.lo })
}

/// How far, in its own sizes, a moving body may travel over the span
/// before t ranged is not tried (`sweeps_far`): the search halves a ranged
/// param to 2^-PARAM_SPLITS of its band, and a body that travels more than
/// that many sizes still sweeps over more than its own size in its finest
/// piece of t, so near another shape no cell of it settles (a chain's pin
/// along its loop, a belt's tooth; a carriage 9 mm wide sliding 580 mm is
/// tried, and proven).
const FAR: f64 = (1u32 << ink_kernel::scope::PARAM_SPLITS) as f64;

/// Why t ranged is not tried: a placed body the claim reads that travels
/// over the span (its boxes over SWEPT_PIECES pieces: `swept`) more than
/// FAR times its own size. Cheap: interval bounds of its placement alone.
fn sweeps_far(v: &mut ModelView, roots: &[u32], span: (f64, f64)) -> Option<String> {
    let placed: Vec<(u32, Xf)> = v.doc.g.reachable(roots).into_iter().filter_map(|n| v.doc.placed.get(&n).map(|p| (p.shape, p.at))).collect();
    let size = |(lo, hi): (V3, V3)| (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
    for (shape, at) in placed {
        if !claims::reads_t(v, &at.r.iter().flatten().chain(&at.t).copied().collect::<Vec<u32>>()) {
            continue;
        }
        let local = match claims::reads_t(v, &[shape]) {
            false => v.whole(shape).then(|| v.bounds(shape)).flatten(),
            true => {
                let over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..v.ctx.clone() };
                v.at_point(over, |v| v.whole(shape).then(|| v.bounds(shape)).flatten())
            }
        };
        let Some((lo, hi, _)) = local else { continue };
        // A cost, not a proof: a placement no isometry is not measured.
        let Some(iso) = Isometry::over(&v.doc.g, &at.r, &at.t, &Ctx { t_span: Some(span), ..v.ctx.clone() }) else { continue };
        let Some(pieces) = swept(v, &iso, (lo, hi), span, false) else { continue };
        let hull = pieces.iter().map(|p| p.1).reduce(ink_core::v3::box_union).expect("pieces");
        let (own, travel) = (size((lo, hi)), size(hull));
        if own > 0.0 && travel > FAR * own {
            return Some(format!(
                "a body it reads travels {travel:.1} over the span, over {FAR} times its size ({own:.2}). Even the finest pieces of t would sweep it by more than its size, so t ranged was not tried"
            ));
        }
    }
    None
}

/// The frame a shape of a pair claim rests in: the placement of a shape
/// placed in `a` (an instance's, a construction's in it, a zone's: the
/// first in node order), else in `b`, that moves over `span`. A claim about two shapes is the same claim of them
/// both moved rigidly at each moment (disjoint, a clearance, a subset,
/// overlaps, touches: distances and depths are kept), so it is asked of
/// them in that frame: `a` rests there and does not move with t, `b` moves
/// only as the two move apart (a pin turning in its hole moves nothing),
/// and the search covers `a`'s own box, not the space the pair sweeps.
/// In world bounds over a piece of the motion each shape blurs by how far
/// its own points move, and their max does not see that they move
/// together; and the space the pair sweeps is searched cell by cell, the
/// contact passing through each (benchmark 32's pin-b-seated: 4.2 million
/// boxes and left to samples; in the crank's frame, proven in seconds).
///
/// Each shape placed (Doc::placements: a field f at R_x^T (p - t_x)) is f
/// placed relative to the frame, at R_a^T R_x, R_a^T (t_x - t_a), the
/// frame's own f itself; anything else read at x, y, z is read at the
/// frame's point R_a p + t_a. Exact in real numbers where R_a is a
/// rotation (its inverse is R_a^T), so a placement is taken only where it
/// is one over the span (rigid::turn_bounds: det exactly 1, r^T r exactly
/// I). A shape turning about z relative to the frame whose move is not one
/// point is placed about its pivot (`pivot`), shown to drift less than the
/// claim's accuracy `tol`, where the claim can be charged its drift
/// (`charge`). Gives the node map, the frame's point (x, y, z of the world
/// as fields of the frame's p and t: `to_world`) and that charge (0: the
/// frame is exact). None: no shape placed in `a` or `b` moves (a shape
/// placed at a moving frame of another, a belt's tooth at its shift, is no
/// instance's: it is searched where it is).
fn moving_frame(v: &mut ModelView, claim: &Claim, span: (f64, f64), tol: Option<f64>) -> Option<(std::collections::HashMap<u32, u32>, [u32; 3], f64)> {
    let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else { return None };
    let over = Ctx { t_span: Some(span), ..v.ctx.clone() };
    // The placements of the shapes placed in `a` (Doc::placed), in node
    // order, then in `b` (a still shape against a moving one): the first
    // that moves and is an isometry over the scope (rigid::Isometry: a
    // rotation, its translation finite). Either is exact: every shape is
    // placed relative to it.
    let placements = v.doc.placements();
    let moves = |xf: &Xf| xf.r.iter().flatten().chain(&xf.t).any(|&n| v.doc.g.depends_on(n, |x| matches!(x, Var::T)));
    let first = |s: u32| {
        let within: std::collections::HashSet<u32> = v.doc.g.reachable(&[s]).into_iter().collect();
        placements.iter().filter(|(w, _)| within.contains(w)).map(|(_, p)| p.at).filter(moves).find_map(|xf| Isometry::over(&v.doc.g, &xf.r, &xf.t, &over))
    };
    let mut frames: Vec<Xf> = [a, b].into_iter().filter_map(first).map(|iso| Xf { r: *iso.r(), t: *iso.t() }).collect();
    frames.dedup();
    if frames.is_empty() {
        return None;
    }
    let read: std::collections::HashSet<u32> = v.doc.g.reachable(&[a, b]).into_iter().collect();
    // In node order, so the nodes made are the same each run.
    let mut placed: Vec<((u32, Xf), u32)> = placements.iter().filter(|(w, _)| read.contains(w)).map(|(w, p)| ((p.shape, p.at), *w)).collect();
    // Each placed relative to the frame, its entries finite over the scope:
    // one infinite at a moment places nothing there, which no relative
    // placement need say alike.
    if !placed.iter().all(|((_, xf), _)| ink_kernel::rigid::finite_over(&v.doc.g, xf.r.iter().flatten().chain(&xf.t), &over)) {
        return None;
    }
    // The world's frame too, where a shape rests in it: a shaft turning
    // about its own axis in a still block's bore moves nothing of the
    // shaft there, where in the shaft's frame the block's screws sweep
    // round it. Asked in the first frame where nothing of the claim moves
    // (reads t), else in the first.
    let g = &mut v.doc.g;
    let (z, o) = (g.constant(0.0), g.constant(1.0));
    let identity = Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: [z, z, z] };
    // A shape no instance places (a still part's own field) is placed
    // where it is, so that in another's frame it is one placed shape too:
    // a still bush round about a turning shaft is unchanged by its turn.
    for s in [a, b] {
        let within: std::collections::HashSet<u32> = g.reachable(&[s]).into_iter().collect();
        if !placed.iter().any(|(_, w)| within.contains(w)) && g.depends_on(s, |x| matches!(x, Var::X | Var::Y | Var::Z)) {
            placed.push(((s, identity), s));
        }
    }
    frames.push(identity);
    // Nothing of the claim rests in any frame where a shape's own field
    // reads t (a belt's teeth), nor in one where a shape turns about no
    // axis of its own fixed on the frame's (`Axes`): those frames, but the
    // first, are not built.
    let own = placed.iter().any(|((f, _), _)| v.doc.g.depends_on(*f, |x| matches!(x, Var::T)));
    let mut first = None;
    for (i, fa) in frames.into_iter().enumerate() {
        let rels = relative(v, fa, &placed, &over);
        let may = !own && rels.iter().all(|(_, _, rel)| rel.as_ref().is_none_or(|rel| !moving(&v.doc.g, rel) || Axes::of(&v.doc.g, rel, &over).is_some()));
        if i > 0 && !may {
            continue;
        }
        let m = in_frame(v, claim, fa, rels, span, tol, &over);
        let g = &mut v.doc.g;
        let still = may
            && [a, b].into_iter().all(|s| {
                let s = g.replace(s, &m.0);
                !g.depends_on(s, |x| matches!(x, Var::T))
            });
        if still {
            return Some(m);
        }
        first.get_or_insert(m);
    }
    first
}

/// Each placed shape's placement relative to the frame `fa`: R_a^T R_x,
/// R_a^T (t_x - t_a), as (node placing it, shape, placement); the frame's
/// own shape where it rests (None).
fn relative(v: &mut ModelView, fa: Xf, placed: &[((u32, Xf), u32)], over: &Ctx) -> Vec<(u32, u32, Option<Xf>)> {
    let over = over.clone();
    let mut rels: Vec<(u32, u32, Option<Xf>)> = Vec::with_capacity(placed.len());
    for &((f, xf), w) in placed {
        if xf == fa {
            rels.push((w, f, None));
            continue;
        }
        let g = &mut v.doc.g;
        let dot = |g: &mut ink_core::graph::Graph, u: [u32; 3], w: [u32; 3]| {
            let p: Vec<u32> = (0..3).map(|k| g.mul(u[k], w[k])).collect();
            let s = g.add(p[0], p[1]);
            g.add(s, p[2])
        };
        let col = |m: &[[u32; 3]; 3], j: usize| [m[0][j], m[1][j], m[2][j]];
        let d: [u32; 3] = std::array::from_fn(|k| g.sub(xf.t[k], fa.t[k]));
        let r: [[u32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| dot(g, col(&fa.r, i), col(&xf.r, j))));
        let t: [u32; 3] = std::array::from_fn(|i| dot(g, col(&fa.r, i), d));
        // What the frame's turn cancels, cancelled (R_a^T R_a = I).
        let entries: Vec<u32> = r.iter().flatten().copied().chain(t).collect();
        let rel = match ink_kernel::rigid::simplified(g, &entries, &over) {
            Some(e) => Xf { r: std::array::from_fn(|i| std::array::from_fn(|j| e[3 * i + j])), t: [e[9], e[10], e[11]] },
            None => Xf { r, t },
        };
        rels.push((w, f, Some(rel)));
    }
    rels
}

/// Whether a placement reads t.
fn moving(g: &ink_core::graph::Graph, xf: &Xf) -> bool {
    xf.r.iter().flatten().chain(&xf.t).any(|&n| g.depends_on(n, |x| matches!(x, Var::T)))
}

/// `moving_frame` in frame `fa` (an isometry over the scope `over`), of
/// the shapes placed in the claim (`placed`: each shape and its placement,
/// with the node placing it).
fn in_frame(
    v: &mut ModelView,
    claim: &Claim,
    fa: Xf,
    rels: Vec<(u32, u32, Option<Xf>)>,
    span: (f64, f64),
    tol: Option<f64>,
    over: &Ctx,
) -> (std::collections::HashMap<u32, u32>, [u32; 3], f64) {
    let over = over.clone();
    let moves = |v: &ModelView, xf: &Xf| moving(&v.doc.g, xf);
    // A proper turn about the shape's own z (an isometry, r_zz exactly 1):
    // its parts every such turn leaves unchanged (a pin, a bore about the
    // axis) placed by the move alone, equal in real numbers
    // (`split_about`). Where its move is not one point over the scope,
    // the turn about a pivot it holds still to within a certified drift
    // (`pivot`), where the claim can be charged it (`charge`).
    // A turn about the shape's own x or y fixed on one of the frame's axes
    // (a shaft along y turned by a rot-y) is the same turn about z seen
    // through relabelled axes (`Axes`): each shape so turning is asked as
    // its own (`seen`: the shape and its placement relabelled), and what is
    // made of it relabelled back.
    let about_z = |v: &ModelView, rel: &Xf| Isometry::over(&v.doc.g, &rel.r, &rel.t, &over).filter(|_| ink_kernel::rigid::constant_over(&v.doc.g, rel.r[2][2], &over) == Some(1.0));
    let mut seen: std::collections::HashMap<u32, (Axes, u32, Xf, Xf)> = std::collections::HashMap::new();
    for (w, f, rel) in &rels {
        let Some(rel) = rel.filter(|rel| moves(v, rel)) else { continue };
        let Some(axes) = Axes::of(&v.doc.g, &rel, &over) else { continue };
        let (f, rel) = axes.seen(&mut v.doc.g, *f, &rel);
        if let Some(turn) = about_z(v, &rel) {
            seen.insert(*w, (axes, f, rel, Xf { r: *turn.r(), t: *turn.t() }));
        }
    }
    let mut pivots: std::collections::HashMap<u32, Pivot> = std::collections::HashMap::new();
    for (w, (_, f, rel, _)) in &seen {
        if !rel.t.iter().all(|&n| ink_kernel::rigid::constant_over(&v.doc.g, n, &over).is_some())
            && let Some(p) = tol.and_then(|tol| pivot(v, *f, rel, span, tol))
        {
            pivots.insert(*w, p);
        }
    }
    let charge = match pivots.is_empty() {
        true => 0.0,
        false => {
            let moved: std::collections::HashMap<u32, f64> = pivots.iter().map(|(w, p)| (*w, p.delta)).collect();
            match self::charge(&v.doc.g, claim, &moved) {
                Some(c) => c,
                None => {
                    pivots.clear();
                    0.0
                }
            }
        }
    };
    let mut map = std::collections::HashMap::new();
    for (w, f, rel) in rels {
        let to = match rel {
            None => f,
            Some(rel) => match (seen.get(&w).cloned(), pivots.get(&w)) {
                (Some((axes, f, rel, _)), Some(p)) => {
                    let g = &mut v.doc.g;
                    // The turn about the pivot: p -> R (q - l) + c.
                    let l = [g.constant(p.l[0]), g.constant(p.l[1])];
                    let t: [u32; 3] = std::array::from_fn(|k| {
                        let (a, b) = (g.mul(rel.r[k][0], l[0]), g.mul(rel.r[k][1], l[1]));
                        let rl = g.add(a, b);
                        let c = g.constant(p.c[k]);
                        g.sub(c, rl)
                    });
                    let turn = Xf { r: rel.r, t };
                    // c - l as a real number (Graph::real), not folded to
                    // its float.
                    let moved: [u32; 3] = std::array::from_fn(|k| {
                        let c = g.constant(p.c[k]);
                        if k < 2 {
                            let c = g.real(c);
                            g.sub(c, l[k])
                        } else {
                            c
                        }
                    });
                    let placed = split_about(v, f, &turn, p.l, moved).unwrap_or_else(|| v.doc.warped(f, &turn));
                    axes.back(&mut v.doc.g, placed)
                }
                (Some((axes, seen_f, _, turn)), None) => match split_about(v, seen_f, &turn, [0.0; 2], turn.t) {
                    Some(placed) => axes.back(&mut v.doc.g, placed),
                    None => v.doc.warped(f, &rel),
                },
                (None, _) => v.doc.warped(f, &rel),
            },
        };
        map.insert(w, to);
    }
    // Anything else in the world at the frame's point: R_a p + t_a.
    let g = &mut v.doc.g;
    let p = [Var::X, Var::Y, Var::Z].map(|x| g.var(x));
    let world: [u32; 3] = std::array::from_fn(|i| {
        let terms: Vec<u32> = (0..3).map(|k| g.mul(fa.r[i][k], p[k])).collect();
        let s = g.add(terms[0], terms[1]);
        let s = g.add(s, terms[2]);
        g.add(s, fa.t[i])
    });
    for i in 0..3 {
        map.entry(p[i]).or_insert(world[i]);
    }
    (map, world, charge)
}

/// A placement relative to a frame (R, T) whose turn keeps an axis of the
/// shape's own frame, k, on an axis of the frame's, j (R e_k = s e_j, s
/// +-1, exactly over the scope), seen as a turn about z: the shape's axes
/// relabelled by P (a cycle, P e_z = e_k) and the frame's by M (a signed
/// cycle, M e_z = s e_j, det 1). The shape f placed, f(R^T (p - T)), is f'
/// = f o P placed by (M^T R P, M^T T), read at M^T p: a turn about z,
/// exactly (its r_zz is s s), and relabelling is a substitution of x, y
/// and z by signed x, y and z, so every value is the same in real numbers.
#[derive(Clone, Copy)]
struct Axes {
    own: [[f64; 3]; 3],
    frame: [[f64; 3]; 3],
}

impl Axes {
    /// The relabelling of `rel`, its own z first: None where no column of
    /// its turn is a signed axis of the frame's over the scope.
    fn of(g: &ink_core::graph::Graph, rel: &Xf, over: &Ctx) -> Option<Axes> {
        // The cycle taking z to axis k: (P q)_i = q_((i + 2 - k) mod 3).
        let cycle = |k: usize| -> [[f64; 3]; 3] { std::array::from_fn(|i| std::array::from_fn(|j| if j == (i + 2 - k) % 3 { 1.0 } else { 0.0 })) };
        [2, 0, 1].into_iter().find_map(|k| {
            let col: Vec<f64> = (0..3).map(|i| ink_kernel::rigid::constant_over(g, rel.r[i][k], over)).collect::<Option<_>>()?;
            let j = (0..3).find(|&j| col[j].abs() == 1.0 && (0..3).all(|i| i == j || col[i] == 0.0))?;
            let s = col[j];
            // M = P_j diag(1, s, s): M e_z = s e_j, det s^2.
            let frame = cycle(j).map(|row| std::array::from_fn(|c| row[c] * if c == 0 { 1.0 } else { s }));
            Some(Axes { own: cycle(k), frame })
        })
    }

    fn identity(&self) -> bool {
        let i: [[f64; 3]; 3] = std::array::from_fn(|a| std::array::from_fn(|b| if a == b { 1.0 } else { 0.0 }));
        self.own == i && self.frame == i
    }

    /// The shape and its placement seen with z the axis: f o P and
    /// (M^T R P, M^T T).
    fn seen(&self, g: &mut ink_core::graph::Graph, f: u32, rel: &Xf) -> (u32, Xf) {
        if self.identity() {
            return (f, *rel);
        }
        let (p, m) = (self.own, self.frame);
        // The one term of a signed permutation's product: sign times node.
        let term = |g: &mut ink_core::graph::Graph, terms: Vec<(f64, u32)>| {
            let (s, n) = terms.into_iter().find(|(s, _)| *s != 0.0).expect("a signed permutation has one term a row");
            if s > 0.0 { n } else { g.neg(n) }
        };
        // (M^T R P)_ab = sum_i,c M_ia R_ic P_cb.
        let r = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let terms = (0..3).flat_map(|i| (0..3).map(move |c| (i, c))).map(|(i, c)| (m[i][a] * p[c][b], rel.r[i][c])).collect();
                term(g, terms)
            })
        });
        let t = std::array::from_fn(|a| term(g, (0..3).map(|i| (m[i][a], rel.t[i])).collect()));
        (relabel(g, f, &p), Xf { r, t })
    }

    /// A field made at the relabelled point read at the frame's: q = M^T p.
    fn back(&self, g: &mut ink_core::graph::Graph, f: u32) -> u32 {
        if self.identity() {
            return f;
        }
        let mt = std::array::from_fn(|a| std::array::from_fn(|b| self.frame[b][a]));
        relabel(g, f, &mt)
    }
}

/// Field `f` read at m q (m a signed permutation): x, y and z replaced by
/// signed x, y and z.
fn relabel(g: &mut ink_core::graph::Graph, f: u32, m: &[[f64; 3]; 3]) -> u32 {
    let xyz = [Var::X, Var::Y, Var::Z];
    let map: Vec<(Var, u32)> = (0..3)
        .map(|i| {
            let (j, s) = (0..3).find_map(|j| (m[i][j] != 0.0).then_some((j, m[i][j]))).expect("a signed permutation");
            let q = g.var(xyz[j]);
            (xyz[i], if s > 0.0 { q } else { g.neg(q) })
        })
        .collect();
    g.subst(f, &map)
}

/// Field `f` with its fixed params read as the floats they are, as
/// `carried` reads centres.
fn fixed_params(v: &mut ModelView, f: u32) -> u32 {
    let ks: Vec<u32> = v
        .doc
        .g
        .reachable(&[f])
        .into_iter()
        .filter_map(|n| match v.doc.g.nodes[n as usize] {
            Node::Var(Var::Param(k)) if v.ctx.p_spans.get(k).is_none() => Some(k),
            _ => None,
        })
        .collect();
    let fixed: Vec<(Var, u32)> = ks.into_iter().map(|k| (Var::Param(k), v.doc.g.constant(v.ctx.params[k as usize]))).collect();
    if fixed.is_empty() { f } else { v.doc.g.subst(f, &fixed) }
}

/// Field `f` placed by `turn` (an isometry, a proper turn about the
/// vertical axis through `centre` of its own frame, then a move), its parts
/// every such turn leaves unchanged (symmetry::carried: a pin, a bore about
/// the axis) placed by the move alone, by `moved`: for such a part F,
/// F(R^T (p - T)) = F(p - (T + R l - l)) exactly, l the centre, and
/// `moved` is T + R l - l (T itself about the origin). A placement is a
/// substitution of x, y and z, which passes through min, max and negation,
/// so each part may be placed on its own. Its fixed params read as the
/// floats they are (`fixed_params`). None: no part is unchanged.
fn split_about(v: &mut ModelView, f: u32, turn: &Xf, centre: [f64; 2], moved: [u32; 3]) -> Option<u32> {
    let f = fixed_params(v, f);
    let invariant = ink_core::symmetry::carried(&v.doc.g, f, centre, None).invariant;
    if invariant.is_empty() {
        return None;
    }
    let (z, o) = (v.doc.g.constant(0.0), v.doc.g.constant(1.0));
    let moved = Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: moved };
    let (tm, mm) = (v.doc.inverse_map(turn), v.doc.inverse_map(&moved));
    fn go(g: &mut ink_core::graph::Graph, n: u32, inv: &std::collections::HashSet<u32>, tm: &[(Var, u32)], mm: &[(Var, u32)]) -> u32 {
        if inv.contains(&n) {
            return g.subst(n, mm);
        }
        match g.nodes[n as usize] {
            Node::Bin(op @ (BinOp::Min | BinOp::Max), a, b) => {
                let (a, b) = (go(g, a, inv, tm, mm), go(g, b, inv, tm, mm));
                g.bin(op, a, b)
            }
            Node::Un(UnOp::Neg, a) => {
                let a = go(g, a, inv, tm, mm);
                g.neg(a)
            }
            // The same value without its name or promise.
            Node::Tag(_, a) | Node::Lipschitz(a, ..) => go(g, a, inv, tm, mm),
            _ => g.subst(n, tm),
        }
    }
    Some(go(&mut v.doc.g, f, &invariant, &tm, &mm))
}

/// A shape's placement relative to a frame, a turn about z, read as a
/// turn about a pivot: the point `l` of the shape's own frame (z 0) that
/// the motion holds within `delta` of the frame's point `c` at every
/// moment and param value of the scope. The shape so placed, p -> R (q -
/// l) + c, is the shape as placed moved by R l + T - c, one translation no
/// longer than `delta` (`charge`).
struct Pivot {
    l: [f64; 2],
    c: [f64; 3],
    delta: f64,
}

/// The pivot of `rel` (an isometry turning about z, its move not one
/// point over `span`), shown to drift less than `tol`: a linkage's
/// coupler seen from its rocker turns about the pin they share, though
/// the coupler is placed through its own angle and the pin through the
/// rocker's, so they meet only by the loop's closure, which no bound of
/// the two placements over a piece of the motion sees. The point is
/// fitted at the moments an estimate reads (claims::SAMPLES: the one
/// point two of their turns carry alike), taken at the centres the
/// shape's parts are written about (symmetry::centres, nearest on each
/// axis), so its parts round it are unchanged by the turn
/// (`split_about`); `c` is its image at the span's middle. Its drift
/// R l + T - c is bounded over pieces of the span, each t's band in the
/// view's scope (mean value bounds in it and the ranged params, what
/// solves on t set derived from it over the piece), the worst piece
/// halved until MAX_STEPS pieces: the largest bound, the sum of its axes'
/// magnitudes, rounded up. None: no turn over the span, no point fits,
/// or its drift is not shown under `tol`.
fn pivot(v: &mut ModelView, f: u32, rel: &Xf, span: (f64, f64), tol: f64) -> Option<Pivot> {
    // Read with t as its param k, each piece of the span k's band in the
    // view's scope, so what solves on t set derived from it follows the
    // piece (Ctx::with_band).
    let k = v.t_param();
    let tk = v.doc.g.var(Var::Param(k));
    let rel = Xf { r: rel.r.map(|row| row.map(|n| v.doc.g.subst(n, &[(Var::T, tk)]))), t: rel.t.map(|n| v.doc.g.subst(n, &[(Var::T, tk)])) };
    let scope = v.ctx.clone();
    let over = |s: (f64, f64)| Ctx { t: ink_core::eval::span_mid(s), t_span: None, ..scope.with_band(k, s) };
    let tape = |v: &ModelView, n: u32| Tape::new(&v.doc.g, n);
    let (r, t) = (rel.r.map(|row| row.map(|n| tape(v, n))), rel.t.map(|n| tape(v, n)));
    let at = |s: f64| {
        let c = over((s, s));
        let e = |x: &Tape| ink_core::eval::eval_interval_scope(x, &c, [interval::I::point(0.0); 3], None).mid();
        (r.each_ref().map(|row| row.each_ref().map(e)), t.each_ref().map(e))
    };
    let n = claims::SAMPLES - 1;
    let moments: Vec<_> = (0..=n).map(|i| at(span.0 + (span.1 - span.0) * i as f64 / n as f64)).collect();
    // The two moments turned furthest apart: (R_0 - R_1) l = T_1 - T_0.
    let (r0, t0) = moments[0];
    let apart = |m: &([[f64; 3]; 3], [f64; 3])| (m.0[0][0] - r0[0][0]).hypot(m.0[1][0] - r0[1][0]);
    let (r1, t1) = *moments.iter().max_by(|x, y| apart(x).total_cmp(&apart(y)))?;
    let (a, b, c, d) = (r0[0][0] - r1[0][0], r0[0][1] - r1[0][1], r0[1][0] - r1[1][0], r0[1][1] - r1[1][1]);
    let det = a * d - b * c;
    if !(det.abs() > 0.0 && det.is_finite()) {
        return None;
    }
    let (u, w) = (t1[0] - t0[0], t1[1] - t0[1]);
    let fit = [(d * u - b * w) / det, (a * w - c * u) / det];
    let f = fixed_params(v, f);
    let centres = ink_core::symmetry::centres(&v.doc.g, f);
    let l: [f64; 2] = std::array::from_fn(|i| *centres[i].iter().min_by(|x, y| (*x - fit[i]).abs().total_cmp(&(*y - fit[i]).abs())).expect("0 is a centre"));
    let image = |m: &([[f64; 3]; 3], [f64; 3])| -> [f64; 3] { std::array::from_fn(|k| m.0[k][0] * l[0] + m.0[k][1] * l[1] + m.1[k]) };
    let c = image(&at(ink_core::eval::span_mid(span)));
    // Not near one point at the moments: no pivot, and no bound sought.
    if moments.iter().any(|m| (0..3).map(|k| (image(m)[k] - c[k]).abs()).sum::<f64>() >= tol) {
        return None;
    }
    let g = &mut v.doc.g;
    let drift: Vec<Tape> = (0..3)
        .map(|k| {
            let (lx, ly, ck) = (g.constant(l[0]), g.constant(l[1]), g.constant(c[k]));
            let (x, y) = (g.mul(rel.r[k][0], lx), g.mul(rel.r[k][1], ly));
            let s = g.add(x, y);
            let s = g.add(s, rel.t[k]);
            let e = g.sub(s, ck);
            Tape::new(g, e)
        })
        .collect();
    let bound = |s: (f64, f64)| -> f64 {
        let at = over(s);
        let b = drift.iter().fold(0.0, |b, e| interval::add_up(b, interval::mag(ink_core::eval::eval_interval_scope(e, &at, [interval::I::point(0.0); 3], None))));
        if b.is_nan() { f64::INFINITY } else { b }
    };
    // Bounds are at least 0: their bits order as they do.
    let mut pieces = std::collections::BinaryHeap::from([(bound(span).to_bits(), span.0.to_bits(), span.1.to_bits())]);
    while pieces.len() < MAX_STEPS {
        let &(b, s0, s1) = pieces.peek().expect("a piece");
        let (s0, s1) = (f64::from_bits(s0), f64::from_bits(s1));
        let m = ink_core::eval::span_mid((s0, s1));
        if f64::from_bits(b) == 0.0 || !(s0 < m && m < s1) {
            break;
        }
        pieces.pop();
        for h in [(s0, m), (m, s1)] {
            pieces.push((bound(h).to_bits(), h.0.to_bits(), h.1.to_bits()));
        }
    }
    let delta = f64::from_bits(pieces.peek().expect("a piece").0);
    (delta < tol).then_some(Pivot { l, c, delta })
}

/// What a pair claim asked of shapes moved (each node of `moved` its
/// world shape moved by one translation no longer than its bound, the
/// rest exact) is charged to say the world's (`ranged_in` asks it so):
/// None where the moved shapes do not keep it within that. A clearance's
/// shapes each moved whole (an operand that is one moved shape or reads
/// none): the gap is off by at most the sum. Disjoint, subset and
/// overlaps ask whether a region of their shapes' points (a ∩ b; a less
/// b) holds a cube of half-side the accuracy: where it is an intersection
/// of parts each moved whole (`literal`), a cube of half-side s in the
/// world's region holds one of s - d in the region asked, d the largest
/// move, and a cube asked of s + d one of s in the world's at each moment.
/// Not where a moved part is one of a union: its members moving apart can
/// open a gap that closes deep inside the cube.
fn charge(g: &ink_core::graph::Graph, claim: &Claim, moved: &std::collections::HashMap<u32, f64>) -> Option<f64> {
    let reads = |n: u32| g.reachable(&[n]).into_iter().any(|m| moved.contains_key(&m));
    // A name or a promise wraps the same shape.
    let whole = |mut n: u32| loop {
        match (moved.get(&n), g.nodes[n as usize]) {
            (Some(&d), _) => return Some(d),
            (None, Node::Tag(_, a) | Node::Lipschitz(a, ..)) => n = a,
            (None, _) => return (!reads(n)).then_some(0.0),
        }
    };
    // The largest move of the parts `n` (a region where `inside`, else its
    // complement) is the intersection of, each moved whole.
    fn literal(g: &ink_core::graph::Graph, n: u32, inside: bool, whole: &dyn Fn(u32) -> Option<f64>) -> Option<f64> {
        if let Some(d) = whole(n) {
            return Some(d);
        }
        match g.nodes[n as usize] {
            Node::Bin(BinOp::Max, a, b) if inside => Some(literal(g, a, inside, whole)?.max(literal(g, b, inside, whole)?)),
            Node::Bin(BinOp::Min, a, b) if !inside => Some(literal(g, a, inside, whole)?.max(literal(g, b, inside, whole)?)),
            Node::Un(UnOp::Neg, a) => literal(g, a, !inside, whole),
            Node::Tag(_, a) | Node::Lipschitz(a, ..) => literal(g, a, inside, whole),
            _ => None,
        }
    }
    match *claim {
        Claim::Clearance(a, b, _) => Some(interval::add_up(whole(a)?, whole(b)?)),
        Claim::Disjoint(a, b) | Claim::Overlaps(a, b) => Some(literal(g, a, true, &whole)?.max(literal(g, b, true, &whole)?)),
        Claim::Subset(a, b) => Some(literal(g, a, true, &whole)?.max(literal(g, b, false, &whole)?)),
        _ => None,
    }
}

/// A pair claim searched in a moving frame (`moving_frame`), for its
/// answer in the world: the frame's point as fields of x, y, z and t's
/// param k (`world`), and the claim's shapes in the world (t read as k).
/// A witness found in the frame is a real point of the world claim's
/// violation (the frame is exact in real numbers), but the point reported
/// is its image in floats, which may be none: a frame 1e16 from the origin
/// rounds a witness off both shapes. So every witness reported is shown
/// one again, of the world's own fields at its moment, by bounds rounded
/// outward (`confirmed`); one that is not leaves the claim to the world's
/// search.
#[derive(Clone, Debug)]
pub(crate) struct Frame {
    world: [u32; 3],
    k: u32,
    span: (f64, f64),
    /// The world's fields: what an emptiness proof's witness must be inside
    /// (max(a, b); a subset's max(a, -b), whose witness is inside a and not
    /// inside b: `subset`), and the two shapes.
    empty: u32,
    subset: bool,
    overlap: u32,
    a: u32,
    b: u32,
    /// The distance a pair of points must be within: a clearance's, or a
    /// contact's (the claim's accuracy).
    d: f64,
    /// The scope the frame was searched in (t's param over the span, what
    /// solves on t set derived from it): a witness's point of it is its
    /// world's (`ctx`).
    scope: Ctx,
    /// How far the frame's shapes may be off the world's (`charge`).
    charge: f64,
}

/// Halvings of a witness's cube `Frame::inside` may make on each axis.
const CONFIRM_SPLITS: u32 = 3;

impl Frame {
    fn new(v: &mut ModelView, claim: &Claim, world: [u32; 3], k: u32, span: (f64, f64), accuracy: Option<f64>, charge: f64) -> Frame {
        let t = v.doc.g.var(Var::Param(k));
        let world = world.map(|n| v.doc.g.subst(n, &[(Var::T, t)]));
        let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else {
            unreachable!("a pair claim has a frame")
        };
        let (a, b) = (v.doc.g.subst(a, &[(Var::T, t)]), v.doc.g.subst(b, &[(Var::T, t)]));
        let overlap = v.doc.g.max(a, b);
        let empty = match claim {
            Claim::Subset(..) => claims::subset_field(&mut v.doc.g, a, b),
            _ => overlap,
        };
        let d = match *claim {
            Claim::Clearance(_, _, d) => d,
            _ => claims::claim_accuracy(v, claim, accuracy),
        };
        let scope = if v.ctx.p_spans.get(k).is_some() { v.ctx.clone() } else { v.ctx.with_band(k, span) };
        Frame { world, k, span, empty, subset: matches!(claim, Claim::Subset(..)), overlap, a, b, d, scope, charge }
    }

    /// The witness's point of the question's scope (`scope`): the view's
    /// ranged params and t's, each at its witness value (`params`) as a
    /// span of that one value, the params derived from them following
    /// (Ctx::point): what a solve on t sets, at the witness's moment. A ranged param keeps its
    /// meaning there, its expressions' values in real numbers bounded
    /// outward: pinned as a float instead, t at 0.5 reads 1e16 + 0.001 t
    /// as its rounded float, and a pair 0.9 off both shapes shows inside.
    fn ctx(&self, params: &[(u32, f64)]) -> Ctx {
        let mut c = self.scope.clone();
        for &(k, x) in params {
            c.params[k as usize] = x;
        }
        c.point()
    }

    /// A witness's params with t's: where the frame took t out of the
    /// question, the witness holds at every moment of the span in the
    /// frame, and is given at its start.
    fn with_t(&self, params: &[(u32, f64)]) -> Vec<(u32, f64)> {
        let mut ps = params.to_vec();
        if !ps.iter().any(|p| p.0 == self.k) {
            ps.push((self.k, self.span.0));
        }
        ps
    }

    /// p of the frame moved into the world at ctx's moment (floats).
    fn point(&self, v: &ModelView, c: &Ctx, p: V3) -> V3 {
        self.world.map(|n| ink_core::eval::eval_point(&Tape::new(&v.doc.g, n), c, p))
    }

    /// Whether field f is below 0 all over box b at ctx's moment, by bounds
    /// rounded outward (halved up to CONFIRM_SPLITS times where undecided).
    fn inside(tape: &Tape, c: &Ctx, b: [interval::I; 3], splits: u32) -> bool {
        let r = ink_core::eval::eval_interval(tape, c, b, None);
        if r.all_inside() {
            return true;
        }
        if splits == 0 || r.none_inside() || b.iter().all(|x| x.lo == x.hi) {
            return false;
        }
        (0..8).all(|m: usize| {
            let h: [interval::I; 3] = std::array::from_fn(|j| {
                let x = b[j];
                let mid = x.mid();
                if m >> j & 1 == 0 { interval::I::new(x.lo, mid) } else { interval::I::new(mid, x.hi) }
            });
            Frame::inside(tape, c, h, splits - 1)
        })
    }

    /// An emptiness proof's witness (a point every point within `tol` of
    /// which is inside the solid: its cube) in the world, shown so of
    /// field f there; None where it is not. Any other outcome as it is (a
    /// proof holds in real numbers; an undecided box is not one of the
    /// world's: dropped).
    fn proof(&self, v: &ModelView, f: u32, p: ink_kernel::ranged::Proof) -> Option<ink_kernel::ranged::Proof> {
        use ink_kernel::{query::Verdict, ranged::Proof};
        let verdict = match p.verdict {
            Verdict::Disproven { witness, params, .. } => {
                let params = self.with_t(&params);
                let c = self.ctx(&params);
                let w = self.point(v, &c, witness);
                // Asked of shapes off the world's (`charge`), a witness is
                // shown one to the claim's own accuracy.
                let tol = if self.charge > 0.0 { p.tol.max(self.d) } else { p.tol };
                let cube = w.map(|x| if tol > 0.0 { interval::I::new((x - tol).next_down(), (x + tol).next_up()) } else { interval::I::point(x) });
                let tape = Tape::new(&v.doc.g, f);
                let shown = match self.subset && f == self.empty {
                    true => {
                        Frame::inside(&Tape::new(&v.doc.g, self.a), &c, cube, CONFIRM_SPLITS)
                            && ink_core::eval::eval_interval(&Tape::new(&v.doc.g, self.b), &c, cube, None).none_inside()
                    }
                    false => Frame::inside(&tape, &c, cube, CONFIRM_SPLITS),
                };
                if !shown {
                    return None;
                }
                let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                Verdict::Disproven { witness: w, t, value: ink_core::eval::eval_point(&tape, &c, w), params }
            }
            // The accuracy it holds to, off by the charge.
            Verdict::Unknown { cells, volume, ended, attained, hunt, .. } => {
                Verdict::Unknown { cells, volume, ended, cell: None, attained: attained.map(|a| if self.charge > 0.0 { interval::add_up(a, self.charge) } else { a }), hunt }
            }
            o => o,
        };
        Some(Proof { verdict, ..p })
    }

    /// How the claim held, its near point (found in the frame) moved into
    /// the world at its moment and confirmed there on the world's field,
    /// as a witness is (`proof`): outright in the frame, it holds outright
    /// in the world where the frame is exact in real numbers (a charged
    /// one's is not: Held::charged first); a near point the
    /// world does not show is not measured.
    fn held(&self, v: &ModelView, h: claims::Held) -> claims::Held {
        let claims::Held::ToAccuracy { near: Some((_, Some((p, at)))) } = h else { return h };
        let params = self.with_t(&at.params);
        let c = self.ctx(&params);
        let w = self.point(v, &c, p);
        let at = ink_kernel::scope::At { t: params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1), params };
        match claims::Held::confirmed_at(false, [(w, at)], &Tape::new(&v.doc.g, self.empty), &c) {
            claims::Held::ToAccuracy { near: Some((_, None)) } => claims::Held::ToAccuracy { near: None },
            h => h,
        }
    }

    /// The evidence with its witnesses in the world, each shown one of the
    /// world's fields at its moment: an emptiness's (or an overlap's) point
    /// inside to its tolerance; a clearance's (or a contact's) pair, each
    /// point in its shape and nearer each other than its distance. None
    /// where one is not.
    fn confirmed(&self, v: &ModelView, e: Evidence) -> Option<Evidence> {
        use ink_kernel::query::Verdict;
        Some(match e {
            Evidence::Empty(p) => Evidence::Empty(self.proof(v, self.empty, p)?),
            Evidence::Overlaps(p, x) => Evidence::Overlaps(self.proof(v, self.overlap, p)?, x),
            Evidence::ClearanceSearch { verdict, pair, res } => {
                let (verdict, pair) = match verdict {
                    Verdict::Disproven { params, .. } => {
                        let params = self.with_t(&params);
                        let c = self.ctx(&params);
                        let (p, q) = pair?;
                        // One in each shape (the claim's way round).
                        let (p, q) = (self.point(v, &c, p), self.point(v, &c, q));
                        let at = |f: u32, x: V3| Frame::inside(&Tape::new(&v.doc.g, f), &c, x.map(interval::I::point), 0);
                        let each = at(self.b, p) && at(self.a, q) || at(self.a, p) && at(self.b, q);
                        if !(each && ink_kernel::query::distance_up(p, q) < self.d) {
                            return None;
                        }
                        let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                        let value = ink_kernel::query::distance_up(p, q);
                        (Verdict::Disproven { witness: p, t, value, params }, Some((p, q)))
                    }
                    // The clearance it holds to, less the charge.
                    Verdict::Unknown { cells, volume, ended, attained, hunt, .. } => {
                        let attained = attained.map(|a| if self.charge > 0.0 { -interval::add_up(-a, self.charge) } else { a });
                        (Verdict::Unknown { cells, volume, ended, cell: None, attained, hunt }, None)
                    }
                    verdict => (verdict, pair),
                };
                Evidence::ClearanceSearch { verdict, pair, res }
            }
            Evidence::Touches { overlap, gap } => Evidence::Touches { overlap: Box::new(self.confirmed(v, *overlap)?), gap: Box::new(self.confirmed(v, *gap)?) },
            e => e,
        })
    }
}

/// Where drive `tape` (read at the origin, every ranged param at every
/// value) may jump in `span`: the first piece, halved down to a
/// MAX_STEPS-th of the span where its rate along t has no bound, whose rate
/// still has none (a select or floor switching in it, whatever its
/// branches' values there). None: bounded over pieces covering the span, so
/// it is continuous.
fn may_jump(tape: &Tape, ctx: &Ctx, span: (f64, f64)) -> Option<(f64, f64)> {
    let least = (span.1 - span.0) / MAX_STEPS as f64;
    let mut todo = vec![span];
    while let Some(s) = todo.pop() {
        if ink_kernel::motion::rate(tape, ctx, [interval::I::point(0.0); 3], s).is_some() {
            continue;
        }
        let m = ink_core::eval::span_mid(s);
        if s.1 - s.0 <= least || !(s.0 < m && m < s.1) {
            return Some(s);
        }
        todo.extend([(m, s.1), (s.0, m)]);
    }
    None
}

/// A witness's value of param k (t ranged: its moment), where it has one.
fn witness_t(e: &Evidence, k: u32) -> Option<f64> {
    use ink_kernel::query::Verdict;
    let of = |ps: &[(u32, f64)]| ps.iter().find(|p| p.0 == k).map(|p| p.1);
    match e {
        Evidence::Empty(p) | Evidence::Overlaps(p, _) => match &p.verdict {
            Verdict::Disproven { params, .. } => of(params),
            _ => None,
        },
        Evidence::ValueSpan(v) => v.witness.as_ref().and_then(|(at, _)| of(&at.params)),
        Evidence::ClearanceSearch { verdict: ink_kernel::query::Verdict::Disproven { params, .. }, .. } => of(params),
        Evidence::Touches { overlap, gap } => witness_t(overlap, k).or_else(|| witness_t(gap, k)),
        _ => None,
    }
}

/// The claim with each node it reads mapped by `f`, for the kinds a band
/// ranges (None: a package's metric, an integral or a mean: measured).
fn with_nodes(c: &Claim, f: &mut impl FnMut(u32) -> u32) -> Option<Claim> {
    Some(match *c {
        Claim::Empty(a) => Claim::Empty(f(a)),
        Claim::Range(a) => Claim::Range(f(a)),
        Claim::Disjoint(a, b) => Claim::Disjoint(f(a), f(b)),
        Claim::Subset(a, b) => Claim::Subset(f(a), f(b)),
        Claim::Overlaps(a, b) => Claim::Overlaps(f(a), f(b)),
        Claim::Touches(a, b) => Claim::Touches(f(a), f(b)),
        Claim::Point { f: a, at, inside } => Claim::Point { f: f(a), at: at.map(&mut *f), inside },
        Claim::Clearance(a, b, d) => Claim::Clearance(f(a), f(b), d),
        Claim::Value(a, cmp, x) => Claim::Value(f(a), cmp, x),
        Claim::Metric { .. } | Claim::Integral { .. } | Claim::Mean { .. } => return None,
    })
}

/// The claim checked at claims::SAMPLES moments of `span`, evenly spaced,
/// its ends included (where a motion starts and stops: a part landing on
/// its stop at the last moment), in order, with its bound evaluated
/// there: an estimate (a failure at a moment is still a witness: the
/// first stops it). A moment what it rests on has no value at (a solve
/// that does not solve there, a simulation stopped at a dead point, a
/// rate not defined where the motion starts) is undecided, with why.
/// Where a provider of what it reads chooses the steps
/// (ParamProvider::steps), at those.
fn measured(v: &mut ModelView, claim: Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Check {
    use ink_kernel::ranged::Status as S;
    let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
    let chosen: Vec<Result<crate::Steps, String>> =
        v.doc.solves_read(&roots).into_iter().filter_map(|i| v.doc.ext[v.doc.solves[i].ext].form.provider(&v.doc)?.steps(&v.ctx.params, span)).collect();
    match chosen.into_iter().collect::<Result<Vec<_>, _>>().map(merged).transpose() {
        Some(Ok(steps)) => return measured_at(v, &claim, bound, &roots, steps, span, opts),
        Some(Err(why)) => {
            let m = Motion { undecided: vec![span], ..Motion::none(S::Unknown, span) };
            let c = finish(Some(checked(Status::Unknown, Evidence::Unsolved(InkError::msg(why)))), m, false, None);
            return Check { motion: c.motion.map(|s| Stepped { basis: Some("not stepped: the steps rest on what has no value (unavailable)".into()), ..s }), ..c };
        }
        None => {}
    }
    let moments = claims::samples(span);
    let mut m = Motion::none(S::Proven, span);
    let mut deciding: Option<Check> = None;
    for &t in &moments {
        if ink_kernel::stop::reached() {
            past_limit(&mut m, (t, t));
            continue;
        }
        let refused = match v.set_moment(t) {
            Err(e) => Some(e),
            Ok(()) => v.refusal(&roots),
        };
        let solved = refused.is_none();
        m.steps += 1;
        let c = match refused {
            Some(e) => checked(Status::Unknown, Evidence::Unsolved(e)),
            None => match validate_at(v, &claim, bound, opts) {
                Some(o) => o,
                None => claims::decide(v, v.doc.claim_at(&claim, bound, &v.ctx), opts).unwrap_or_else(|e| checked(Status::Unknown, Evidence::Unsolved(e))),
            },
        };
        match c.status {
            Status::Fails => {
                (m.status, m.t, deciding) = (S::Disproven, Some(t), Some(certified(v, &claim, bound, c, opts)));
                break;
            }
            Status::Holds => {
                deciding.get_or_insert(c);
            }
            Status::Unknown => {
                m.undecided.push((t, t));
                m.status = S::Unknown;
                m.stopped += usize::from(solved && ink_kernel::stop::reached());
                deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(c));
            }
        }
    }
    let step = (moments.len() > 1).then(|| (span.1 - span.0) / (moments.len() - 1) as f64);
    finish(deciding, m, false, step)
}

/// `claim` (with its `bound`, a field) judged at the view's moment
/// (claims::validate) where it, its bound included, reads a solve on t,
/// before anything is put in its place (the bound's number there:
/// Doc::claim_at): over the motion it was judged before it was asked
/// (claims::ask_in) with what such a solve sets ranged where it could be;
/// where not (one only measured), it is judged here, at each moment
/// measured, as it is asked.
fn validate_at(v: &mut ModelView, claim: &Claim, bound: Option<u32>, opts: ClaimOptions) -> Option<Check> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    if !v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t) {
        return None;
    }
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let at = Ctx { t_span: None, ..v.ctx.clone() };
    claims::validate(v, claim, &roots, tol, opts.res, at, None)
}

/// A moment (or step) of a measured claim not checked: the call's time
/// limit had passed (ink_kernel::stop). A moment checked and left unknown
/// once it has is the limit's too (its searches stopped); one a solve
/// refuses is not.
fn past_limit(m: &mut Motion, s: (f64, f64)) {
    m.undecided.push(s);
    m.status = ink_kernel::ranged::Status::Unknown;
    m.stopped += 1;
}

/// The steps several providers chose over one span, each a claim reads
/// (an effort and a reaction, two reactions): every end any chose, a step
/// unchecked or loose where it lies in one some provider left so, their
/// bases each once. None: none chose.
fn merged(all: Vec<crate::Steps>) -> Option<crate::Steps> {
    let mut it = all.into_iter();
    let first = it.next()?;
    let rest: Vec<crate::Steps> = it.collect();
    if rest.is_empty() {
        return Some(first);
    }
    let every: Vec<&crate::Steps> = std::iter::once(&first).chain(&rest).collect();
    let mut ends: Vec<f64> = every.iter().flat_map(|s| s.ends.iter().copied()).collect();
    ends.sort_by(f64::total_cmp);
    ends.dedup();
    let within = |s: (f64, f64), of: &[(f64, f64)]| of.iter().any(|o| o.0 <= s.0 && s.1 <= o.1);
    let steps: Vec<(f64, f64)> = ends.windows(2).map(|w| (w[0], w[1])).collect();
    let unchecked = steps.iter().copied().filter(|&s| every.iter().any(|x| within(s, &x.unchecked))).collect();
    let loose = steps.iter().copied().filter(|&s| every.iter().any(|x| within(s, &x.loose))).collect();
    let mut bases: Vec<String> = Vec::new();
    for b in every.iter().map(|s| s.basis.clone()) {
        if !bases.contains(&b) {
            bases.push(b);
        }
    }
    Some(crate::Steps { ends, unchecked, loose, basis: bases.join("; and at the steps of what else it reads: ") })
}

/// The claim measured at the span's ends and the middles of the steps a
/// provider chose (`steps`), in order: the first failing is the witness;
/// a moment it rests on a solve with no value at, undecided, with the
/// solve's error; the steps the provider could not choose as it says,
/// unchecked or loose. Holding, its evidence is the moment where it came
/// nearest to failing (a value claim's value of least margin: its greatest
/// under an upper bound, its least over a lower).
fn measured_at(v: &mut ModelView, claim: &Claim, bound: Option<u32>, roots: &[u32], steps: crate::Steps, span: (f64, f64), opts: ClaimOptions) -> Check {
    let mut m = Motion { unchecked: steps.unchecked.clone(), loose: steps.loose.clone(), ..Motion::none(ink_kernel::ranged::Status::Proven, span) };
    let skipped = |s: (f64, f64)| steps.unchecked.contains(&s) || steps.loose.contains(&s);
    let (mut deciding, mut least): (Option<Check>, f64) = (None, f64::INFINITY);
    let middles = steps.ends.windows(2).map(|w| (w[0], w[1])).filter(|&s| !skipped(s)).map(|s| (ink_core::eval::span_mid(s), s));
    let moments: Vec<(f64, (f64, f64))> = std::iter::once((span.0, (span.0, span.0))).chain(middles).chain([(span.1, (span.1, span.1))]).collect();
    for (t, s) in moments {
        if ink_kernel::stop::reached() {
            past_limit(&mut m, s);
            continue;
        }
        m.steps += 1;
        let refused = match v.set_moment(t) {
            Err(e) => Some(e),
            Ok(()) => v.refusal(roots),
        };
        if let Some(e) = refused {
            m.undecided.push(s);
            deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(checked(Status::Unknown, Evidence::Unsolved(e))));
            continue;
        }
        if let Some(o) = validate_at(v, claim, bound, opts) {
            m.undecided.push(s);
            deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(o));
            continue;
        }
        let c = v.doc.claim_at(claim, bound, &v.ctx);
        let decided = match banded(v, &c, roots, t, opts) {
            Some(c) => Ok(c),
            None => claims::decide(v, c, opts),
        };
        let Ok(c) = decided else {
            m.undecided.push(s);
            continue;
        };
        match c.status {
            Status::Fails => {
                (m.status, m.t) = (ink_kernel::ranged::Status::Disproven, Some(t));
                let c = certified(v, claim, bound, c, opts);
                return Check { motion: Some(Stepped { basis: Some(steps.basis.clone()), ..Stepped::of(m, false, None) }), ..c };
            }
            Status::Unknown => {
                m.undecided.push(s);
                m.stopped += usize::from(ink_kernel::stop::reached());
                deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(c));
            }
            Status::Holds => {
                let x = match (&c.evidence, claim) {
                    (Evidence::Value(x), Claim::Value(_, Cmp::Lt | Cmp::Le, bound)) => bound - x,
                    (Evidence::Value(x), Claim::Value(_, _, bound)) => x - bound,
                    (Evidence::Value(x), _) => *x,
                    _ => f64::INFINITY,
                };
                if deciding.as_ref().is_none_or(|d| d.status == Status::Holds && x < least) {
                    (deciding, least) = (Some(c), x);
                }
            }
        }
    }
    if !(m.unchecked.is_empty() && m.undecided.is_empty() && m.loose.is_empty()) {
        m.status = ink_kernel::ranged::Status::Unknown;
    }
    let c = finish(deciding, m, false, None);
    Check { motion: c.motion.map(|s| Stepped { basis: Some(steps.basis), ..s }), ..c }
}

/// A failure found at a sampled moment (the view at it), a guarantee only
/// where the same claim, built from its own nodes with t the param `ranged`
/// asks over the span, registered at the moment's one value (every other
/// ranged input as it is), fails too: then with that answer's witness. A
/// sample reads t as a number, and a check's fast path may read what its
/// nodes give there in floats ((1e16 + t) - 1e16 is t, 0 in floats; a
/// probe moved by it frozen where it rounds): asked so, its points and
/// values are the original expressions' in real numbers, bounded outward
/// (Frame::ctx's range before its point). Otherwise an estimate. A claim
/// whose nodes read no t froze nothing: its sample is the question.
fn certified(v: &mut ModelView, claim: &Claim, bound: Option<u32>, c: Check, opts: ClaimOptions) -> Check {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    if c.accuracy.is_none() || !roots.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)) {
        return c;
    }
    let estimate = |c: Check| Check { accuracy: None, ..c };
    let k = v.t_param();
    let tk = v.doc.g.var(Var::Param(k));
    let g = &mut v.doc.g;
    let mut on = |n: u32| g.subst(n, &[(Var::T, tk)]);
    let Some(asked) = with_nodes(claim, &mut on) else { return estimate(c) };
    let asked = match (asked, bound.map(&mut on)) {
        (Claim::Value(f, cmp, _), Some(b)) => Claim::Value(v.doc.g.sub(f, b), cmp, 0.0),
        (_, Some(_)) => return estimate(c),
        (a, None) => a,
    };
    let (ctx, motion) = (v.ctx.clone(), v.motion.take());
    let t = v.ctx.t;
    v.ctx.params[k as usize] = t;
    v.ctx.p_spans = v.ctx.p_spans.clone().with_band(k, (t, t));
    let r = claims::decide(v, asked, opts);
    (v.ctx, v.motion) = (ctx, motion);
    match r {
        Ok(r) if r.status == Status::Fails && r.accuracy.is_some() => Check { evidence: r.evidence, ..c },
        _ => estimate(c),
    }
}

/// A claim at moment t reading params whose provider reads them over the
/// step holding t within a band (ParamProvider::band): a round reading as
/// the band's centre, the real values what it set. A value claim holds
/// where it holds for every value in the band (interval bounds; its value
/// then the band's worst), fails where it fails at the real values, and is
/// undecided between, saying so; another claim fails where it fails at the
/// real values, else is undecided. A band of no finite width decides
/// nothing by holding. None: no band.
fn banded(v: &mut ModelView, claim: &Claim, roots: &[u32], t: f64, opts: ClaimOptions) -> Option<Check> {
    let mut band: Vec<(usize, f64, f64)> = Vec::new();
    for i in v.doc.solves_read(roots) {
        let s = &v.doc.solves[i];
        let Some(p) = v.doc.ext[s.ext].form.provider(&v.doc) else { continue };
        if let Some((w, centre)) = p.band(&v.ctx.params, t) {
            band.extend(s.owns.iter().zip(w.into_iter().zip(centre)).map(|(&k, (w, c))| (k, w, c)));
        }
    }
    let read = v.doc.params_read(roots);
    band.retain(|b| read.contains(&b.0) && b.1 != 0.0);
    if band.is_empty() {
        return None;
    }
    let tol = claims::claim_accuracy(v, claim, opts.accuracy);
    let x = match claim {
        Claim::Value(f, cmp, bound) => {
            let comparison = ink_kernel::ranged::Comparison { cmp: *cmp, bound: *bound, tol };
            let tape = Tape::new(&v.doc.g, *f);
            if band.iter().all(|b| b.1.is_finite()) {
                let spans: Vec<(u32, (f64, f64))> = band.iter().map(|&(k, w, c)| (k as u32, (c - w, c + w))).collect();
                let mut centre = v.ctx.clone();
                for &(k, _, c) in &band {
                    centre.params[k] = c;
                }
                let over = Ctx { p_spans: v.ctx.p_spans.clone().with_derived(&spans), ..centre };
                let iv = ink_core::eval::eval_interval(&tape, &over, [interval::I::point(0.0); 3], None);
                if comparison.holds_on(iv) {
                    let worst = if matches!(cmp, crate::Cmp::Lt | crate::Cmp::Le) { iv.hi } else { iv.lo };
                    return Some(Check { accuracy: None, ..checked(Status::Holds, Evidence::Value(worst)) });
                }
            }
            let x = ink_core::eval::eval_point(&tape, &v.ctx, [0.0; 3]);
            if comparison.fails_on(interval::I::point(x)) || x.is_nan() {
                return Some(Check { accuracy: None, ..checked(Status::Fails, Evidence::Value(x)) });
            }
            Some(x)
        }
        _ => {
            let at = v.doc.claim_at(claim, None, &v.ctx);
            match claims::decide(v, at, opts) {
                Ok(c) if c.status == Status::Fails => return Some(c),
                _ => None,
            }
        }
    };
    let widths: Vec<String> = band.iter().map(|&(k, w, _)| format!("{} within {w:.6}", v.doc.params[k].name)).collect();
    let why = format!(
        "what it reads is set at the moment's real value, and over its step read with parts that spin at rest, round, the real value within a band of \
         that ({}, in their units: motion.basis): neither holding over that band nor failing at this moment's real value{}",
        widths.join(", "),
        x.map_or(String::new(), |x| format!(" ({x:.6})"))
    );
    Some(checked(Status::Unknown, Evidence::Unsolved(InkError::msg(why))))
}

/// The stepped answer: its status from the steps, with the evidence of the
/// step that decided it (a disproof, the first undecided, else the
/// tightest proven). Measured between steps (not `bound`), it is an
/// estimate, except a failure: the failing step's own answer, a witness.
fn finish(deciding: Option<Check>, m: Motion, bound: bool, step: Option<f64>) -> Check {
    let status = Status::of(m.status);
    let c = deciding.unwrap_or_else(|| checked(status, Evidence::Stepped));
    let accuracy = c.accuracy.filter(|_| bound || status == Status::Fails);
    Check { status, accuracy, in_window: c.in_window.filter(|_| status == Status::Unknown), motion: Some(Stepped::of(m, bound, step)), ..c }
}

/// The rigid group item `it` rides over the motion `span`: how many of its
/// links, from the top, reach the last whose placement moves (0: none
/// does, the world's). Items of one group (the same instance at that link)
/// keep their placements relative to each other at every moment: nothing
/// below it moves, and every placement down to it is proven an isometry
/// over the span (rigid::Isometry: a rotation exactly, its translation
/// finite), so what they share moves them both rigidly. Sharing a moving link is not
/// enough: a link whose matrix scales with t (a turn about an axis given
/// too small to normalize exactly) stretches the distance between its
/// parts as it turns, and one placement above it that is no rotation
/// stretches the turning distances by direction, and one whose
/// translation is infinite at a moment (a move by ln((t - 0.5)^2)) places
/// nothing there. None: its own field moves (it deforms), or a placement
/// down to its group is not so proven.
pub fn rigid_group(v: &ModelView, it: &crate::Item, span: (f64, f64)) -> Option<usize> {
    if claims::moves(v, &[it.shape], span) {
        return None;
    }
    let moves = |l: &crate::Link| claims::moves(v, &l.at.r.iter().flatten().chain(&l.at.t).copied().collect::<Vec<u32>>(), span);
    let k = it.links.iter().rposition(moves).map_or(0, |i| i + 1);
    let over = Ctx { t_span: Some(span), ..v.ctx.clone() };
    it.links[..k].iter().all(|l| Isometry::over(&v.doc.g, &l.at.r, &l.at.t, &over).is_some()).then_some(k)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A view of `src` over the span, and its def `name`'s field.
    fn view(src: &str, span: (f64, f64)) -> ModelView {
        static CORE: crate::Registry = crate::Registry::new(&[]);
        let doc = crate::Doc::load_str_in(&CORE, src, "t", std::path::Path::new("."), &[]).unwrap();
        let (ctx, _) = doc.ctx_at(&[], span.0).unwrap();
        ModelView::new(doc, ctx)
    }

    fn field(v: &ModelView, name: &str) -> u32 {
        let crate::Value::Field(f, _) = v.doc.def(name).unwrap().clone() else { panic!("{name} is no field") };
        f
    }

    /// The window `drive` gives the claim that the def `sa` is disjoint of
    /// `sb`, in turns of the drive.
    fn tooth(src: &str) -> Option<f64> {
        let mut v = view(src, (0.0, 1.0));
        let (a, b) = (field(&v, "sa"), field(&v, "sb"));
        let k = v.t_param();
        let d = drive(&mut v, &[a, b], (0.0, 1.0), k)?;
        // The drive turns 10 times over the span.
        Some(10.0 * (d.window.1 - d.window.0) / d.travel)
    }

    /// A claim reads a spin-up (t^2) only through the turns it drives: a
    /// gear and a post repeat by one of the gear's teeth, a 15:45 pair by
    /// one of the pinion's (the wheel turns a whole tooth with it), and an
    /// arm fixed to the gear only by a whole turn. A gear tilted, turning
    /// about its own z, repeats by a tooth; turned about its x, by a turn.
    /// Waves of one drive summed repeat every turn of it; waves with no
    /// common period leave the turn they drive repeating.
    #[test]
    fn a_drive_repeats_by_the_teeth_it_turns() {
        let head = "(import \"std\")\n(t-span 0 1)\n(def drive (* 3600 (square t)))\n\
                    (fn gear (n) (extrude-z (involute-gear-2d n 1 20 0.05) 0 5))\n\
                    (body g12 (shape (gear 12)))\n(body g15 (shape (gear 15)))\n(body g45 (shape (gear 45)))\n(body post (shape (box 2 2 5)))\n";
        let near = |got: Option<f64>, want: f64| matches!(got, Some(w) if (w - want).abs() < 1e-9 * want);
        let post =
            format!("{head}(assembly root (instance g g12 (at) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))");
        assert!(near(tooth(&post), 1.0 / 12.0), "{:?}", tooth(&post));
        let pair = format!(
            "{head}(assembly root (instance p g15 (at) (turn z drive)) (instance w g45 (at (move (mesh-distance 15 45 1) 0 0)) (turn z (mesh-turn 15 45 0 drive))))\n\
             (def sa (shape-of root/p))\n(def sb (shape-of root/w))"
        );
        assert!(near(tooth(&pair), 1.0 / 15.0), "{:?}", tooth(&pair));
        let arm = format!(
            "{head}(assembly geared (instance g g12) (instance a post (at (move 0 8 2.5))))\n\
             (assembly root (instance g geared (at) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n\
             (def sa (shape-of root/g))\n(def sb (shape-of root/p))"
        );
        assert!(near(tooth(&arm), 1.0), "{:?}", tooth(&arm));
        let tilted = format!(
            "{head}(assembly root (instance g g12 (at (rot-x 10)) (turn z drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))"
        );
        assert!(near(tooth(&tilted), 1.0 / 12.0), "{:?}", tooth(&tilted));
        let off =
            format!("{head}(assembly root (instance g g12 (at) (turn x drive)) (instance p post (at (move 9 0 2.5))))\n(def sa (shape-of root/g))\n(def sb (shape-of root/p))");
        assert!(near(tooth(&off), 1.0), "{:?}", tooth(&off));
        // Two waves of one spin-up, summed: their branches rejoin above the
        // drive, which still repeats every turn.
        let mut v = view("(t-span 0 1)\n(def a (* 3600000000 (square t)))\n(def f (+ (sin (deg->rad a)) (cos (deg->rad a))))", (0.0, 1.0));
        let (f, k) = (field(&v, "f"), v.t_param());
        let d = drive(&mut v, &[f], (0.0, 1.0), k).expect("a drive");
        assert!((1e7 * (d.window.1 - d.window.0) / d.travel - 1.0).abs() < 1e-9, "{:?}", d.window);
        // Waves of it with no common period (sin of a in radians, cos of it
        // in degrees) drive an outer turn b: tried deepest first, the
        // waves' drive is shown no period and b's turn is.
        let mut v = view(
            "(t-span 0 1)\n(def a (* 3600000000 (square t)))\n(def b (* 1000000000 (+ (sin a) (cos (deg->rad a)))))\n\
             (def f (sin (deg->rad b)))",
            (0.0, 1.0),
        );
        let (f, k) = (field(&v, "f"), v.t_param());
        let d = drive(&mut v, &[f], (0.0, 1.0), k).expect("a drive");
        assert!(d.window.1 - d.window.0 < 1e-6 * d.travel, "{:?}", d.window);
    }

    /// Finding a drive costs a few readings of the claim, however many
    /// values of t it reads in a row: a value 4,096 additions from t that
    /// repeats in nothing is read through t and the top of the run (the
    /// claim), not again through each addition.
    #[test]
    fn a_drive_is_found_in_a_few_readings() {
        let n = 4096;
        let src: String = std::iter::once("(t-span 0 1)\n(def a0 t)\n".to_string()).chain((1..=n).map(|i| format!("(def a{i} (+ a{} 1))\n", i - 1))).collect();
        let mut v = view(&src, (0.0, 1.0));
        let (f, k) = (field(&v, &format!("a{n}")), v.t_param());
        let tries = tries(&v, &[f], k).expect("t is read");
        assert!(tries.len() <= 2, "{} readings", tries.len());
        assert!(drive(&mut v, &[f], (0.0, 1.0), k).is_none());
    }

    /// A frame's witness and pair, as the search in it gives them (frame
    /// coordinates; t's param absent where the frame cancels t), checked
    /// in the world by `Frame::confirmed`.
    fn witnesses(v: &mut ModelView, claim: Claim, span: (f64, f64), y: V3) -> (Option<Evidence>, Option<Evidence>) {
        use ink_kernel::query::Verdict;
        use ink_kernel::ranged::Proof;
        let (_, world, _) = moving_frame(v, &claim, span, None).expect("a moving frame");
        let k = v.t_param();
        let fr = Frame::new(v, &claim, world, k, span, None, 0.0);
        let proof = Proof { verdict: Verdict::Disproven { witness: y, t: 0.0, value: -0.1, params: vec![] }, res: 0.05, tol: 0.05, exact: false, loose: Vec::new() };
        let pair = Evidence::ClearanceSearch { verdict: Verdict::Disproven { witness: y, t: 0.0, value: 0.0, params: vec![] }, pair: Some((y, y)), res: 0.05 };
        (fr.confirmed(v, Evidence::Empty(proof)), fr.confirmed(v, pair))
    }

    /// Two balls 0.1 apart along x carried 1e16 out: in the first's frame
    /// their overlap's point (1.05, 0, 0) is a witness, but moved into the
    /// world in floats it rounds to x = 1e16, outside both: neither it nor a
    /// pair there is reported. Turned together by two instances instead, the frame
    /// cancels t: the witness is given at the span's start, in the world
    /// there (R(90 deg) (10.25, 0, 0)), and the pair too.
    #[test]
    fn a_frame_witness_is_reported_only_where_the_world_shows_it() {
        let far = "(import \"std\")\n(t-span 0 1)\n(body a (shape (move (sphere 0.25) 1 0 0)))\n(body b (shape (move (sphere 0.25) 1.1 0 0)))\n\
                   (assembly root (instance a a (at (move (+ 1e16 (* 0.001 t)) 0 0))) (instance b b (at (move (+ 1e16 (* 0.001 t)) 0 0))))\n\
                   (def sa (shape-of root/a))\n(def sb (shape-of root/b))";
        let mut v = view(far, (0.0, 1.0));
        let (a, b) = (field(&v, "sa"), field(&v, "sb"));
        let (e, p) = witnesses(&mut v, Claim::Disjoint(a, b), (0.0, 1.0), [1.05, 0.0, 0.0]);
        assert!(e.is_none() && p.is_none(), "{e:?} {p:?}");
        let (_, p) = witnesses(&mut v, Claim::Clearance(a, b, 0.5), (0.0, 1.0), [1.05, 0.0, 0.0]);
        assert!(p.is_none(), "{p:?}");

        let turned = "(import \"std\")\n(t-span 0.25 1)\n(body a (shape (move (sphere 1) 10 0 0)))\n(body b (shape (move (sphere 1) 10.5 0 0)))\n\
                      (assembly root (instance a a (at (rot-z (* 360 t)))) (instance b b (at (rot-z (* 360 t)))))\n\
                      (def sa (shape-of root/a))\n(def sb (shape-of root/b))";
        let mut v = view(turned, (0.25, 1.0));
        let (a, b) = (field(&v, "sa"), field(&v, "sb"));
        let k = v.t_param();
        let (e, p) = witnesses(&mut v, Claim::Disjoint(a, b), (0.25, 1.0), [10.25, 0.0, 0.0]);
        let Some(Evidence::Empty(ink_kernel::ranged::Proof { verdict: ink_kernel::query::Verdict::Disproven { witness, t, params, .. }, .. })) = e else { panic!("{e:?}") };
        assert!(t == 0.25 && params.contains(&(k, 0.25)) && witness[0].abs() < 1e-9 && (witness[1] - 10.25).abs() < 1e-9, "{witness:?} at {t}, {params:?}");
        let Some(Evidence::ClearanceSearch { pair: Some((p, q)), verdict: ink_kernel::query::Verdict::Disproven { t, .. }, .. }) = p else { panic!("{p:?}") };
        assert!(t == 0.25 && p == q && (p[1] - 10.25).abs() < 1e-9, "{p:?} {q:?} at {t}");
    }

    /// A sampled failure stays a guarantee only where the claim asked of its
    /// own nodes, t a param of the moment's one value, fails too (`certified`):
    /// one a pointwise check found in floats, ((1e16 + t) - 1e16 > 0.001)
    /// failing at its value 0, is an estimate (the real value is t, 0.0156
    /// here); (t > 0.02) failing at t = 0.0156 is certified.
    #[test]
    fn a_sampled_failure_is_a_guarantee_only_once_its_own_claim_fails_at_the_moment() {
        let mut v = view("(def far (- (+ 1e16 t) 1e16))\n(def near t)", (0.01, 0.02));
        v.set_moment(0.015625).unwrap();
        let sampled = |x: f64| Check { accuracy: Some(1e-6), ..checked(Status::Fails, Evidence::Value(x)) };
        let opts = ClaimOptions::default();
        let (far, near) = (field(&v, "far"), field(&v, "near"));
        let c = certified(&mut v, &Claim::Value(far, Cmp::Gt, 0.001), None, sampled(0.0), opts);
        assert!(c.status == Status::Fails && c.accuracy.is_none(), "{c:?}");
        let c = certified(&mut v, &Claim::Value(near, Cmp::Gt, 0.02), None, sampled(0.015625), opts);
        assert!(c.accuracy == Some(1e-6), "{c:?}");
    }

    /// tests/motion/loop-closed-pin.ink's coupler, seen from its rocker,
    /// turns about the pin they share, though the eye is placed through the
    /// coupler's angle and the pin through the rocker's: the eye is placed
    /// about that pivot, its drift bounded far under the accuracy and
    /// charged. Seen from a coupler whose eye is 0.3 off the pin, the pin
    /// turns about the point the rocker shares with that coupler: asked
    /// there, where nothing of the claim moves, its drift charged. A union
    /// of the eye and a part at rest is not one shape moved, so it is not
    /// charged one move and is placed as it is.
    #[test]
    fn a_loop_closed_joint_is_searched_about_its_pivot() {
        let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/motion/loop-closed-pin.ink")).unwrap();
        let src = format!(
            "{src}\n(def pin (max (shape-of root/rocker rocker-pin) (slab-z 11 15)))\n(def in-eye (shape-of root/coupler eye))\n(def in-far (shape-of root/long far-eye))\n\
             (def coupler-w (shape-of root/coupler))\n(def rocker-w (shape-of root/rocker))\n(def with-wall (min (shape-of root/coupler eye) (move (box 1 1 1) 100 0 0)))"
        );
        let mut v = view(&src, (0.0, 1.0));
        let [pin, in_eye, in_far, coupler, rocker, with_wall] = ["pin", "in-eye", "in-far", "coupler-w", "rocker-w", "with-wall"].map(|n| field(&v, n));
        let charge = |v: &mut ModelView, c: Claim| moving_frame(v, &c, (0.0, 1.0), Some(0.05)).expect("a moving frame").2;
        let seated = charge(&mut v, Claim::Subset(pin, in_eye));
        assert!(seated > 0.0 && seated < 1e-3, "{seated}");
        let clear = charge(&mut v, Claim::Clearance(coupler, rocker, 0.05));
        assert!(clear > 0.0 && clear < 2e-3, "{clear}");
        let far = charge(&mut v, Claim::Subset(pin, in_far));
        assert!(far > 0.0 && far < 1e-3, "{far}");
        assert_eq!(charge(&mut v, Claim::Disjoint(pin, with_wall)), 0.0);
        assert_eq!(charge(&mut v, Claim::Subset(pin, in_eye)), seated);
    }

    /// A hold of shapes asked off the world's (a pivot's drift, a drive's
    /// slack, charged to the accuracy) is to the accuracy, never outright:
    /// the world's may overlap by up to the charge. Two cylinders 0.01 into
    /// each other at t = 0, one turning about the other's axis on a centre
    /// that drifts 0.04 (its pivot), as asserted and as interfere's pair;
    /// and a 15:45 gear pair apart, through a drive whose wheel turns at a
    /// float rate (its slack).
    #[test]
    fn a_charged_hold_is_never_outright() {
        let src = "(import \"std\")\n(t-span 0 1)\n(body a (shape (cylinder-z 1 2)))\n(body b (shape (cylinder-z 1 2)))\n\
                   (def shift (- 2.01 (* 0.02 (cos (deg->rad (* 360 t))))))\n\
                   (assembly root (instance a a (at (rot-z (* 360 t)))) (instance b b (at (rot-z (* 720 t)) (move shift 0 0) (rot-z (* 360 t)))))\n\
                   (assert apart (disjoint (shape-of root/a) (shape-of root/b)))";
        let gears = "(import \"std\")\n(t-span 0 1)\n(def drive (* 3600 (square t)))\n(fn gear (n) (extrude-z (involute-gear-2d n 1 20 0.05) 0 5))\n\
                     (body g15 (shape (gear 15)))\n(body g45 (shape (gear 45)))\n\
                     (assembly root (instance p g15 (at) (turn z drive)) (instance w g45 (at (move (+ 2 (mesh-distance 15 45 1)) 0 0)) (turn z (mesh-turn 15 45 0 drive))))\n\
                     (assert apart (disjoint (shape-of root/p) (shape-of root/w)))";
        let held = |src: &str| {
            let mut v = view(src, (0.0, 1.0));
            v.set_motion(v.doc.t_span).unwrap();
            let r = claims::run_asserts(&mut v, None, false, ClaimOptions::default()).unwrap();
            r.asserts[0].check.clone()
        };
        for src in [src, gears] {
            let c = held(src);
            assert!(c.holds_guaranteed() && matches!(c.held, Some(claims::Held::ToAccuracy { .. })), "{c:?}");
        }
        let mut v = view(src, (0.0, 1.0));
        v.set_motion(v.doc.t_span).unwrap();
        let (items, target) = (v.doc.items("root").unwrap(), v.doc.target("root"));
        let plan = crate::interfere::plan(&mut v, items, target.as_ref(), None, None).unwrap();
        let pairs = plan.run(&mut v, None, &|_| None).unwrap();
        let crate::interfere::Answer::Moving { check } = &pairs[0].answer else { panic!("a moving pair") };
        assert!(check.holds_guaranteed() && matches!(check.held, Some(claims::Held::ToAccuracy { .. })), "{check:?}");
    }
}
