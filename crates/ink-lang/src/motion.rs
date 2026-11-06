//! Claims over a motion (`(t-span A B)`, `--t-span`). Time is a range of
//! the question, as a tolerance band is: a claim that reads t is the same
//! claim with t a ranged param over the span (`ranged`), decided by the
//! one search over every moment at once, so what holds or fails is
//! guaranteed as it is over any band. What that search leaves undecided
//! (within the call's time limit), or cannot ask (what reads a solve on t,
//! a package's metric, an integral), is measured at the middles of
//! uniform steps (`measured`): holding there is an estimate, failing at
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

/// Uniform steps a measured claim is checked at.
pub const MEASURED_STEPS: usize = 32;

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
        ink_kernel::ranged::prove_key(tape, ctx, *lo, *hi, opts)
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
        if bound.is_none() && riding(v, &claim, span) {
            // Asked once, at the view's moment: their relative placement is
            // the same at every moment.
            let motion = v.motion.take();
            let c = claims::decide(v, claim, opts);
            v.motion = motion;
            let m = Motion::none(if c.as_ref().is_ok_and(|c| c.status == Status::Unknown) { ink_kernel::ranged::Status::Unknown } else { ink_kernel::ranged::Status::Proven }, span);
            let c = c?;
            return Ok(Asked::Done(Check { note: c.note.or(Some(RIDING)), motion: Some(Stepped { riding: true, ..Stepped::of(m, true, None) }), ..c }));
        }
        // A claim reading a solve on t (a simulated motion's effort) is
        // measured only, never proven over t: a drive that may jump is said
        // beside its estimate. Any other is unknown there.
        let roots: Vec<u32> = v.doc.claim_roots(&claim).into_iter().chain(bound).collect();
        let measured_only = v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t);
        let jump = jumps(v, &claim, bound, span);
        if let (Some(c), false) = (&jump, measured_only) {
            return Ok(Asked::Done(c.clone()));
        }
        // The sampled answer first: cheap, always there, an estimate (a
        // failure at a moment is a witness, guaranteed: nothing to prove).
        let sampled = measured(v, claim.clone(), bound, span, opts);
        v.set_moment(t0)?;
        if sampled.status == Status::Fails && sampled.accuracy.is_some() {
            return Ok(Asked::Done(sampled));
        }
        let keep = |sampled: Check, why: String, answer: Option<Check>| Check { motion: sampled.motion.map(|s| Stepped { left: Some(why), left_answer: answer.map(Box::new), ..s }), ..sampled };
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
                        Ranged::Wait(search, whole, read, world) => Later::Wait(
                            search,
                            Box::new(move |v, found| {
                                let p = found;
                                // A witness in the frame shown none in the
                                // world: the sampled answer stands.
                                let e = match &world {
                                    Some(w) if p.status() != ink_kernel::ranged::Status::Unknown => match w.confirmed(v, Evidence::Empty(p.clone())) {
                                        Some(e) => e,
                                        None => return Ok(None),
                                    },
                                    _ => Evidence::Empty(p.clone()),
                                };
                                let c = Check { ranged: read, ..claims::in_window(Status::of(p.status()), whole, e) };
                                Ok(match decided(c, span, false) {
                                    Ranged::Decided(c) => Some(claims::finish_later(v, &asked, bound, c, opts)),
                                    _ => None,
                                })
                            }),
                        ),
                    })
                }),
            });
            return Ok(Asked::Done(keep(sampled, "t ranged is tried once every pair has an answer, with the time left: it did not decide it".into(), None)));
        }
        // Then t ranged, with the time left, to prove it or find a witness;
        // undecided, the sampled answer stands, saying why.
        let r = ranged(v, &claim, bound, span, opts);
        v.set_moment(t0)?;
        let jumped = |why: String| match &jump {
            Some(Check { evidence: Evidence::NotCovered(j), .. }) => format!("{why}, and {j}, so between the moments measured the speed may pass its limits"),
            _ => why,
        };
        Ok(Asked::Done(match r? {
            Ranged::Decided(c) => c,
            Ranged::Left(why, answer) => keep(sampled, jumped(why), answer),
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
    /// shapes are proven bounded, the ranged params it reads, and the
    /// frame it is searched in (`Frame`).
    Wait(Search, bool, Vec<(u32, (f64, f64))>, Option<Frame>),
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
/// holds. Not ranged: a claim of a kind no band ranges (a package's metric, an
/// integral, a thickness, cavities), one reading a solve on t (its root is
/// found at each moment, not over a range of them), one whose bound moves.
/// Asked alone (assert), it stops at half the call's time left (stop::within);
/// beside other searches (interfere), at the call's limit.
fn ranged(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions) -> Result<Ranged, InkError> {
    match ranged_in(v, claim, bound, span, opts, true)? {
        Some(r) => Ok(r),
        // Its witness in the frame is none in the world: asked there.
        None => Ok(ranged_in(v, claim, bound, span, opts, false)?.expect("a world search's answer needs no confirming")),
    }
}

/// `ranged`, in the first shape's frame where `framed` and it has one
/// (`moving_frame`). None: decided there by a witness that is not shown
/// one in the world (`Frame::confirmed`).
fn ranged_in(v: &mut ModelView, claim: &Claim, bound: Option<u32>, span: (f64, f64), opts: ClaimOptions, framed: bool) -> Result<Option<Ranged>, InkError> {
    let roots: Vec<u32> = v.doc.claim_roots(claim).into_iter().chain(bound).collect();
    if v.doc.solves_read(&roots).iter().any(|&i| v.doc.solves[i].on_t) {
        return Ok(Some(Ranged::Left("it reads a solve on t, found at each moment, not over a range of them".into(), None)));
    }
    if let Some(why) = sweeps_far(v, &roots, span) {
        return Ok(Some(Ranged::Left(why, None)));
    }
    let k = v.t_param();
    let t = v.doc.g.var(Var::Param(k));
    // Searched where one of the two shapes rests (`frame`): the same claim
    // at every moment, its points moved rigidly.
    let frame = if framed { moving_frame(v, claim, span) } else { None };
    let world = frame.as_ref().map(|f| Frame::new(v, claim, f.1, k, span, opts.accuracy));
    let mut on = |n: u32| {
        let n = match &frame {
            Some(m) => v.doc.g.replace(n, &m.0),
            None => n,
        };
        v.doc.g.subst(n, &[(Var::T, t)])
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
    let (ctx, motion) = (v.ctx.clone(), v.motion.take());
    v.ctx.params[k as usize] = ink_core::eval::span_mid(span);
    v.ctx.p_spans = v.ctx.p_spans.clone().with_band(k, span);
    let read = v.doc.claim_roots(&asked);
    let r = match (opts.beside, &asked) {
        (true, Claim::Disjoint(a, b)) => {
            let tol = claims::claim_accuracy(v, &asked, opts.accuracy);
            let f = v.doc.g.max(*a, *b);
            let tape = v.clipped(f);
            let crate::Reach { lo, hi, whole } = v.reach(f);
            Ok(Some(Ranged::Wait(Search(tape, v.ctx.clone(), (lo, hi), v.proof_options(opts.res, tol)), whole, v.ranged_reads(&read), world)))
        }
        _ => {
            let (c, own) = ink_kernel::stop::within(if opts.beside { 1.0 } else { 0.5 }, || claims::decide(v, asked, opts));
            let ranged = v.ranged_reads(&read);
            (v.ctx, v.motion) = (ctx, motion);
            let c = c?;
            // Its witnesses in the world, each shown one there (undecided,
            // there are none).
            let c = match &world {
                Some(w) if c.status != Status::Unknown => match w.confirmed(v, c.evidence.clone()) {
                    Some(evidence) => Check { evidence, ..c },
                    None => return Ok(None),
                },
                _ => c,
            };
            return Ok(Some(decided(Check { ranged, ..c }, span, own)));
        }
    };
    (v.ctx, v.motion) = (ctx, motion);
    r
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
/// I). Gives the node map and the frame's point (x, y, z of the world as
/// fields of the frame's p and t: `to_world`). None: no shape placed in `a`
/// or `b` moves (a shape placed at a moving frame of another, a belt's
/// tooth at its shift, is no instance's: it is searched where it is).
fn moving_frame(v: &mut ModelView, claim: &Claim, span: (f64, f64)) -> Option<(std::collections::HashMap<u32, u32>, [u32; 3])> {
    let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else { return None };
    let over = Ctx { t_span: Some(span), ..v.ctx.clone() };
    // The placements of the shapes placed in `a` (Doc::placed), in node
    // order, then in `b` (a still shape against a moving one): the first
    // that moves and is an isometry over the scope (rigid::Isometry: a
    // rotation, its translation finite). Either is exact: every shape is
    // placed relative to it.
    let placements = v.doc.placements();
    let moves = |xf: &Xf| xf.r.iter().flatten().chain(&xf.t).any(|&n| v.doc.g.depends_on(n, |x| matches!(x, Var::T)));
    let iso = [a, b].into_iter().find_map(|s| {
        let within: std::collections::HashSet<u32> = v.doc.g.reachable(&[s]).into_iter().collect();
        placements.iter().filter(|(w, _)| within.contains(w)).map(|(_, p)| p.at).filter(moves).find_map(|xf| Isometry::over(&v.doc.g, &xf.r, &xf.t, &over))
    })?;
    let fa = Xf { r: *iso.r(), t: *iso.t() };
    let read: std::collections::HashSet<u32> = v.doc.g.reachable(&[a, b]).into_iter().collect();
    // In node order, so the nodes made are the same each run.
    let placed: Vec<((u32, Xf), u32)> = placements.iter().filter(|(w, _)| read.contains(w)).map(|(w, p)| ((p.shape, p.at), *w)).collect();
    // Each placed relative to the frame, its entries finite over the scope:
    // one infinite at a moment places nothing there, which no relative
    // placement need say alike.
    if !placed.iter().all(|((_, xf), _)| ink_kernel::rigid::finite_over(&v.doc.g, xf.r.iter().flatten().chain(&xf.t), &over)) {
        return None;
    }
    let mut map = std::collections::HashMap::new();
    // Each placed shape at its placement relative to the frame's: R_a^T
    // R_x, R_a^T (t_x - t_a); the frame's own shape where it rests.
    for ((f, xf), w) in placed {
        let to = if xf == fa {
            f
        } else {
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
            // A proper turn about the shape's own z (an isometry, r_zz
            // exactly 1): its parts every such turn leaves unchanged (a pin,
            // a bore about the axis) placed by the move alone, equal in real
            // numbers (`split_about`).
            let about_z = Isometry::over(&v.doc.g, &rel.r, &rel.t, &over).filter(|_| ink_kernel::rigid::constant_over(&v.doc.g, rel.r[2][2], &over) == Some(1.0));
            match about_z {
                Some(turn) => split_about(v, f, &turn).unwrap_or_else(|| v.doc.warped(f, &rel)),
                None => v.doc.warped(f, &rel),
            }
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
    Some((map, world))
}

/// Field `f` placed by `turn` (an isometry, a proper turn about its own z
/// axis through its origin, then a move), its parts every such turn leaves unchanged
/// (symmetry::carried: a pin, a bore about the axis) placed by the move
/// alone: for such a part F, F(R^T (p - T)) = F(p - T) exactly. A placement
/// is a substitution of x, y and z, which passes through min, max and
/// negation, so each part may be placed on its own. Its fixed params read
/// as the floats they are, as `carried` reads centres. None: no part is
/// unchanged.
fn split_about(v: &mut ModelView, f: u32, turn: &Isometry) -> Option<u32> {
    let turn = &Xf { r: *turn.r(), t: *turn.t() };
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
    let f = if fixed.is_empty() { f } else { v.doc.g.subst(f, &fixed) };
    let invariant = ink_core::symmetry::carried(&v.doc.g, f, [0.0, 0.0], None).invariant;
    if invariant.is_empty() {
        return None;
    }
    let (z, o) = (v.doc.g.constant(0.0), v.doc.g.constant(1.0));
    let moved = Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: turn.t };
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
    /// (max(a, b); a subset's max(a, -b)), and the two shapes.
    empty: u32,
    overlap: u32,
    a: u32,
    b: u32,
    /// The distance a pair of points must be within: a clearance's, or a
    /// contact's (the claim's accuracy).
    d: f64,
}

/// Halvings of a witness's cube `Frame::inside` may make on each axis.
const CONFIRM_SPLITS: u32 = 3;

impl Frame {
    fn new(v: &mut ModelView, claim: &Claim, world: [u32; 3], k: u32, span: (f64, f64), accuracy: Option<f64>) -> Frame {
        let t = v.doc.g.var(Var::Param(k));
        let world = world.map(|n| v.doc.g.subst(n, &[(Var::T, t)]));
        let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else {
            unreachable!("a pair claim has a frame")
        };
        let (a, b) = (v.doc.g.subst(a, &[(Var::T, t)]), v.doc.g.subst(b, &[(Var::T, t)]));
        let overlap = v.doc.g.max(a, b);
        let empty = match claim {
            Claim::Subset(..) => {
                let nb = v.doc.g.neg(b);
                v.doc.g.max(a, nb)
            }
            _ => overlap,
        };
        let d = match *claim {
            Claim::Clearance(_, _, d) => d,
            _ => claims::claim_accuracy(v, claim, accuracy),
        };
        Frame { world, k, span, empty, overlap, a, b, d }
    }

    /// The witness's point of the question's scope: the view's ranged
    /// params and t's (over the span, as the search asked it), each at its
    /// witness value (`params`) as a span of that one value, the params
    /// derived from them following (Ctx::point). A ranged param keeps its
    /// meaning there, its expressions' values in real numbers bounded
    /// outward: pinned as a float instead, t at 0.5 reads 1e16 + 0.001 t
    /// as its rounded float, and a pair 0.9 off both shapes shows inside.
    fn ctx(&self, v: &ModelView, params: &[(u32, f64)]) -> Ctx {
        let mut c = v.ctx.with_band(self.k, self.span);
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
    /// rounded outward (halved up to CONFIRM_SPLITS times where undecided);
    /// `closed`: at most 0 (a pair's points may lie on their surfaces).
    fn inside(tape: &Tape, c: &Ctx, b: [interval::I; 3], closed: bool, splits: u32) -> bool {
        let r = ink_core::eval::eval_interval(tape, c, b, None);
        if r.hi < 0.0 || closed && r.hi <= 0.0 {
            return true;
        }
        if splits == 0 || r.lo >= 0.0 || b.iter().all(|x| x.lo == x.hi) {
            return false;
        }
        (0..8).all(|m: usize| {
            let h: [interval::I; 3] = std::array::from_fn(|j| {
                let x = b[j];
                let mid = x.mid();
                if m >> j & 1 == 0 { interval::I::new(x.lo, mid) } else { interval::I::new(mid, x.hi) }
            });
            Frame::inside(tape, c, h, closed, splits - 1)
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
                let c = self.ctx(v, &params);
                let w = self.point(v, &c, witness);
                let tol = p.tol;
                let cube = w.map(|x| if tol > 0.0 { interval::I::new((x - tol).next_down(), (x + tol).next_up()) } else { interval::I::point(x) });
                let tape = Tape::new(&v.doc.g, f);
                if !Frame::inside(&tape, &c, cube, false, CONFIRM_SPLITS) {
                    return None;
                }
                let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                Verdict::Disproven { witness: w, t, value: ink_core::eval::eval_point(&tape, &c, w), params }
            }
            Verdict::Unknown { cells, volume, ended, attained, .. } => Verdict::Unknown { cells, volume, ended, cell: None, attained },
            o => o,
        };
        Some(Proof { verdict, ..p })
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
                        let c = self.ctx(v, &params);
                        let (p, q) = pair?;
                        // One in each shape (the claim's way round).
                        let (p, q) = (self.point(v, &c, p), self.point(v, &c, q));
                        let at = |f: u32, x: V3| Frame::inside(&Tape::new(&v.doc.g, f), &c, x.map(interval::I::point), true, 0);
                        let each = at(self.b, p) && at(self.a, q) || at(self.a, p) && at(self.b, q);
                        if !(each && ink_kernel::query::distance_up(p, q) < self.d) {
                            return None;
                        }
                        let t = params.iter().find(|x| x.0 == self.k).map_or(self.span.0, |x| x.1);
                        let value = ink_kernel::query::distance_up(p, q);
                        (Verdict::Disproven { witness: p, t, value, params }, Some((p, q)))
                    }
                    Verdict::Unknown { cells, volume, ended, attained, .. } => (Verdict::Unknown { cells, volume, ended, cell: None, attained }, None),
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

/// The claim checked at the middles of MEASURED_STEPS uniform steps of
/// `span` and at its two ends, with its bound evaluated there: an estimate
/// (a failure at a moment is still a witness). Where a provider of what it reads chooses the steps
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
    let n = MEASURED_STEPS;
    let width = (span.1 - span.0) / n as f64;
    let mut m = Motion::none(S::Proven, span);
    let mut deciding: Option<Check> = None;
    for i in 0..n {
        let s = (span.0 + width * i as f64, if i + 1 == n { span.1 } else { span.0 + width * (i + 1) as f64 });
        let t = ink_core::eval::span_mid(s);
        if ink_kernel::stop::reached() {
            past_limit(&mut m, s);
            continue;
        }
        // A moment what it rests on has no value at (a solve that does not
        // solve there, a simulation stopped at a dead point) is undecided,
        // with why.
        let refused = match v.set_moment(t) {
            Err(e) => Some(e),
            Ok(()) => v.refusal(&roots),
        };
        let c = match refused {
            Some(e) => Some(checked(Status::Unknown, Evidence::Unsolved(e))),
            None => claims::decide(v, v.doc.claim_at(&claim, bound, &v.ctx), opts).ok(),
        };
        m.steps += c.is_some() as usize;
        match c {
            Some(c) if c.status == Status::Fails => {
                (m.status, m.t, deciding) = (S::Disproven, Some(t), Some(certified(v, &claim, bound, c, opts)));
                break;
            }
            Some(c) if c.status == Status::Holds => {
                deciding.get_or_insert(c);
            }
            c => {
                m.undecided.push(s);
                m.status = S::Unknown;
                m.stopped += usize::from(ink_kernel::stop::reached());
                if let Some(c) = c.filter(|_| deciding.as_ref().is_none_or(|d| d.status != Status::Unknown)) {
                    deciding = Some(c);
                }
            }
        }
    }
    // The span's ends too: where a motion starts and stops (a part landing
    // on its stop at the last moment, which no step's middle reaches), and
    // a rate not defined at one (a dead point the motion starts or ends at)
    // leaves the claim unknown there; a failure is a witness.
    if m.status != S::Disproven {
        let ends = if span.0 < span.1 { vec![span.0, span.1] } else { vec![span.0] };
        for t in ends {
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
                None => match claims::decide(v, v.doc.claim_at(&claim, bound, &v.ctx), opts) {
                    Ok(c) => c,
                    Err(e) => checked(Status::Unknown, Evidence::Unsolved(e)),
                },
            };
            match c.status {
                Status::Fails => {
                    (m.status, m.t, deciding) = (S::Disproven, Some(t), Some(certified(v, &claim, bound, c, opts)));
                    break;
                }
                Status::Unknown => {
                    m.undecided.push((t, t));
                    m.status = S::Unknown;
                    m.stopped += usize::from(solved && ink_kernel::stop::reached());
                    deciding = deciding.filter(|d| d.status == Status::Unknown).or(Some(c));
                }
                Status::Holds => {}
            }
        }
    }
    finish(deciding, m, false, Some(width))
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
        let crate::Value::Field(f) = v.doc.def(name).unwrap().clone() else { panic!("{name} is no field") };
        f
    }

    /// A frame's witness and pair, as the search in it gives them (frame
    /// coordinates; t's param absent where the frame cancels t), checked
    /// in the world by `Frame::confirmed`.
    fn witnesses(v: &mut ModelView, claim: Claim, span: (f64, f64), y: V3) -> (Option<Evidence>, Option<Evidence>) {
        use ink_kernel::query::Verdict;
        use ink_kernel::ranged::Proof;
        let (_, world) = moving_frame(v, &claim, span).expect("a moving frame");
        let k = v.t_param();
        let fr = Frame::new(v, &claim, world, k, span, None);
        let proof = Proof { verdict: Verdict::Disproven { witness: y, t: 0.0, value: -0.1, params: vec![] }, res: 0.05, tol: 0.05, exact: false };
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
}
