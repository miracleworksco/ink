//! Solves on t followed along t (docs/IR.md, Solving constraints). A
//! solve whose params are one root among others of as many equations
//! (ParamProvider::follows) sets, at moment t, the root reached from its
//! root at load time (t = 0) by moving along t, not whichever one a
//! numerical method finds from a guess there. Each step of t has a box (its
//! tube) holding one root of the equations at every moment of the step,
//! and holding the root followed to the step's start (Roots::enclose on
//! its equations: Krawczyk's test with t ranged over the step), the
//! equations' rate along t bounded over it (Roots::continuous: they do not
//! jump in it), so inside the step the root is that one, moved along t; it
//! cannot pass to another. At a moment of the step it is the provider's solve from the
//! root followed, where that lands in the box, else the box narrowed
//! (Roots::inside). One with more equations than unknowns has no root
//! certified to exist (roots.rs), so no tube: it is not followed, but
//! solved at each moment from the params' values there, as one not on t
//! is; and so is every solve that reads what it sets, through the solves
//! on t between, none of them bounded over a step.
//!
//! A step is twice the last one, or as long where that one was halved (the
//! first 1/16 s; as long as certifying allows, however slowly the
//! equations change), and is halved where it is not certified.
//! Where one under 2^-40 s is not, two roots meet or the root followed ends
//! (a dead point), or the equations switch there and are not shown
//! continuous (Roots::enclose and Roots::continuous say which): which root follows
//! past there is not known, so every moment past it is refused, saying
//! why. So, in this call alone, is every moment past where following had
//! got when the call's time limit passed (ink_kernel::stop). A solve with
//! no regular root at t = 0 has none to follow: every other moment is
//! refused. One reading what an earlier solve on t sets is
//! certified over that one's enclosures over the step: a followed one's
//! narrowed from its steps, one with one value at each moment (a solve's
//! rates) bounded over the step as a root is, its equations linear in it;
//! and along their path, moving at their rates (roots::Moving).
//!
//! A solve whose equations read t only through one expression of t (a
//! clock: t itself, or a crank's angle 10 t + t^2 / 1000) and repeat in it
//! (ink_core::period), and the solves on t they read repeat in it too, is
//! followed one period, where its root is shown to be the one at 0 moved
//! by exact amounts that leave the equations unchanged (Follow::fold): a
//! moment past it or before 0, the clock not jumping between, is a moment
//! of that period where the clock is where it is, moved by whole periods
//! (Follow::folded); where the period's enclosure leaves the count of
//! periods to it uncertain, any of the first period's moments, certified
//! at the moment itself.
//!
//! The steps depend only on where following starts, never on the moments
//! asked for, and are kept (Paths) by the solve and the values of what its
//! root rests on that does not vary along t: a moment's root is the same
//! whichever moments were asked for first, in a document's views, its
//! packages' runs and its imports' alike. They are kept on disk too
//! (ink_kernel::cache), CHUNK at a time, keyed by all they rest on
//! (Follow::key): a later command reads the steps an earlier one
//! certified, each the step it would certify, its box included.

use crate::eval::{Ctx, ParamSpans, span_mid};
use crate::ext::ParamProvider;
use crate::roots::{Moving, Roots};
use crate::tape::Tape;
use ink_core::eval;
use ink_core::graph::{BinOp, Fingerprint, Var};
use ink_core::interval::{self, I, add_o, mul_o, sub_o};
use ink_core::period::{Clock, Coef, Shown, clocks, common_base, keeps, period, period_in};
use ink_kernel::cache;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};

/// The first step tried from t = 0, and the shortest: one not certified
/// there is a dead point (s).
const FIRST: f64 = 1.0 / 16.0;
const SHORTEST: f64 = 1.0 / (1u64 << 40) as f64;

/// Steps kept on disk together: a way's steps n CHUNK to (n + 1) CHUNK
/// are its chunk n, stored once full.
pub const CHUNK: usize = 64;

/// The kind of entry the steps are kept under (ink_kernel::cache).
const KIND: &str = "follow";

thread_local! {
    /// Steps this thread certified (Follow::step), not read from disk.
    pub static CERTIFIED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// The steps each solve was followed by, by its index in Doc::ext and the
/// bits of what its root rests on that does not vary along t: the params it
/// and the solves on t under it read that no solve on t sets, and those
/// solves' params at load time (where each starts from at t = 0).
pub type Paths = Arc<Mutex<HashMap<(usize, Vec<u64>), Arc<Mutex<Path>>>>>;

/// A document's solves in their order (Doc::solves), compiled, each on t
/// followed along t from load time.
pub struct Follow {
    solves: Vec<Solve>,
    /// The params' values at load time, before the solves on t are solved
    /// at t = 0 from them.
    load: Vec<f64>,
    paths: Paths,
}

struct Solve {
    ext: usize,
    p: Arc<dyn ParamProvider>,
    /// Its equations (Form::equations), whose roots in its params it sets.
    equations: Option<Vec<Tape>>,
    /// What it is, for messages ("solve 'loop'").
    label: String,
    /// It is on t (ext::Solve::on_t), and followed: on t, a root among
    /// others, and not measured.
    on_t: bool,
    followed: bool,
    /// Its root is certified nowhere, only measured at each moment: it has
    /// more equations than unknowns, or reads a solve on t that is so.
    measured: bool,
    /// The solves on t it reads from (indices in `solves`).
    upstream: Vec<usize>,
    /// The params its path rests on that no solve on t sets.
    fixed: Vec<usize>,
    /// The params of the solves on t its path rests on, itself included.
    chain: Vec<usize>,
}

/// One solve's path from t = 0: its root there, and the steps away from it
/// each way.
pub struct Path {
    /// The root at t = 0 and its enclosure, once found; Err: why it has
    /// none to follow.
    anchor: Option<Result<(Vec<f64>, Vec<I>), String>>,
    /// Forward (t > 0) and backward (t < 0).
    ways: [Way; 2],
    /// Its period, once looked for, and once it is passed, whether it holds
    /// (Follow::fold).
    period: Option<Option<Period>>,
    fold: Option<Option<Fold>>,
}

/// A solve's period (Follow::period): its equations read t through `clock`
/// and are unchanged by the clock moved by p (exactly, and enclosed), each
/// solve on t they read (`upstream`: its index, and how many of its own
/// periods p is) moved by its own periods. The clock is g0 at t = 0, and
/// has moved by more than p from there by t = `reach`, up (`way` 1) or down
/// (-1).
#[derive(Clone)]
struct Period {
    clock: Clock,
    base: Coef,
    p: Coef,
    pe: I,
    upstream: Vec<(usize, i128)>,
    g0: I,
    reach: f64,
    way: f64,
}

/// A period that holds (Follow::fold): each unknown's move over it,
/// exactly (None: 0), and enclosed.
#[derive(Clone)]
struct Fold {
    exact: Vec<Option<Coef>>,
    shift: Vec<I>,
}

struct Way {
    /// Its key on disk (from Follow::key's), None: kept in memory only.
    disk: Option<u128>,
    /// The unknowns of each step.
    n: usize,
    /// The steps found, in order away from 0, each starting where the last
    /// ends, in chunks of CHUNK but the last: each chunk's far end, and its
    /// steps once read (a full one is on disk).
    chunks: Vec<(f64, Option<Vec<Step>>)>,
    /// Why it is not followed past the last step's end, once found.
    end: Option<String>,
}

/// A span of t folded into a path's first period (Follow::folded): pieces
/// of the period from 0 to p (enclosed), each with the whole numbers m of
/// periods its moments may be moved back by (a range, where the period's
/// enclosure leaves the count uncertain). The root followed to a moment of
/// the span is the one followed to a moment of a piece, its unknowns moved
/// by m `shift`s, for some m of that piece's.
struct Folded {
    pieces: Vec<(I, I)>,
    shift: Vec<I>,
}

impl Folded {
    /// Each unknown's move by a piece's m shifts, enclosed.
    fn moved(&self, m: I) -> Vec<I> {
        self.shift.iter().map(|d| mul_o(m, *d)).collect()
    }
}

/// A step of t, from `a` (nearer 0) to `b`: the box holding one root at
/// every moment of it, the root at each end (with its enclosure at b,
/// where the next step starts), and whether a longer one was tried first
/// (then the next is as long, else twice as long).
#[derive(Clone)]
struct Step {
    a: f64,
    b: f64,
    tube: Vec<I>,
    ra: Vec<f64>,
    rb: Vec<f64>,
    eb: Vec<I>,
    halved: bool,
}

impl Follow {
    /// The solves in their order: the form each is (index in Doc::ext), its
    /// provider, its equations compiled (Form::equations), the params it
    /// reads (Solve::reads), whether it is on t (Solve::on_t), and what it
    /// is; `load` the params at load time; `paths` the steps found so far.
    pub fn new(solves: Vec<(usize, Arc<dyn ParamProvider>, Option<Vec<Tape>>, Vec<usize>, bool, String)>, load: Vec<f64>, paths: Paths) -> Follow {
        let on_t: Vec<bool> = solves.iter().map(|s| s.4).collect();
        let owner = |k: usize| solves.iter().position(|s| s.1.owns().contains(&k));
        let mut out: Vec<Solve> = Vec::with_capacity(solves.len());
        for (i, (ext, p, equations, reads, t, label)) in solves.iter().enumerate() {
            let mut upstream: Vec<usize> = reads.iter().filter_map(|&k| owner(k)).filter(|&j| j < i && on_t[j]).collect();
            upstream.sort_unstable();
            upstream.dedup();
            let mut fixed: Vec<usize> = reads.iter().copied().filter(|&k| owner(k).is_none_or(|j| !on_t[j])).collect();
            let mut chain: Vec<usize> = p.owns().to_vec();
            for &j in &upstream {
                fixed.extend(&out[j].fixed);
                chain.extend(&out[j].chain);
            }
            fixed.sort_unstable();
            fixed.dedup();
            chain.sort_unstable();
            chain.dedup();
            let measured = equations.as_ref().is_some_and(|e| e.len() != p.owns().len()) || upstream.iter().any(|&j| out[j].measured);
            let followed = *t && p.follows() && equations.is_some() && !measured;
            out.push(Solve { ext: *ext, p: p.clone(), equations: equations.clone(), label: label.clone(), on_t: *t, followed, measured, upstream, fixed, chain });
        }
        Follow { solves: out, load, paths }
    }

    /// Whether solve `i` is followed along t (Follow::provide), the params
    /// it sets, and what it is ("solve 'loop'").
    pub fn followed(&self, i: usize) -> (bool, &[usize], &str) {
        (self.solves[i].followed, self.solves[i].p.owns(), &self.solves[i].label)
    }

    /// Set solve `i`'s params at moment t by its provider alone, from their
    /// values in `params` (its root near them: a simulation's state from
    /// the last one's), not followed along t.
    pub fn solve_near(&self, i: usize, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        self.solves[i].p.provide(params, t)
    }

    /// Set solve `i`'s params at moment t: on t, the root followed there
    /// (at t = 0, its provider's from their values at load time, as at
    /// load); else as its provider sets them from params' values. What
    /// falls short of the accuracy its provider states, why; Err: why it
    /// sets nothing.
    pub fn provide(&self, i: usize, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        let s = &self.solves[i];
        if s.on_t && !t.is_finite() {
            return Err(format!("{} reads t, and t = {t} is not a moment", s.label));
        }
        if s.followed && t == 0.0 {
            for &k in s.p.owns() {
                params[k] = self.load[k];
            }
        }
        if !s.followed || t == 0.0 {
            return s.p.provide(params, t);
        }
        let not = |m: String| format!("{} is not followed along t from t = 0 to {t}: {m}", s.label);
        if let Some(f) = self.folded(i, params, (t, t)) {
            return self.provide_folded(i, params, t, &f).map_err(not);
        }
        // Its provider's solve from the root followed, between the step's
        // ends: in the step's box, it is the one root there.
        let step = self.steps(i, params, (t, t)).map_err(not)?.swap_remove(0);
        let w = (t - step.a) / (step.b - step.a);
        for ((&k, a), b) in s.p.owns().iter().zip(&step.ra).zip(&step.rb) {
            params[k] = a + w * (b - a);
        }
        let solved = s.p.provide(params, t);
        if solved.is_ok() && s.p.owns().iter().zip(&step.tube).all(|(&k, x)| x.lo <= params[k] && params[k] <= x.hi) {
            return solved;
        }
        let (r, _) = self.over(i, params, (t, t)).map_err(not)?;
        for (&k, v) in s.p.owns().iter().zip(r) {
            params[k] = v;
        }
        Ok(None)
    }

    /// Solve `i`'s path for what its root rests on in `params`, with the
    /// steps kept on disk for it.
    fn path(&self, i: usize, params: &[f64]) -> Arc<Mutex<Path>> {
        let s = &self.solves[i];
        let bits: Vec<u64> = s.fixed.iter().map(|&k| params[k]).chain(s.chain.iter().map(|&k| self.load[k])).map(f64::to_bits).collect();
        if let Some(p) = self.paths.lock().unwrap().get(&(s.ext, bits.clone())) {
            return p.clone();
        }
        let key = self.key(i, params);
        let ways = [false, true].map(|back| Way::kept(key.map(|k| Fingerprint::default().key(k).word(back as u64).value()), s.p.owns().len()));
        self.paths.lock().unwrap().entry((s.ext, bits)).or_insert_with(|| Arc::new(Mutex::new(Path { anchor: None, ways, period: None, fold: None }))).clone()
    }

    /// Solve `i`'s key on disk, for what its root rests on in `params`:
    /// whether it is followed; its equations as the kernel sees them, with
    /// the values of the params they read (ink_kernel::cache::question),
    /// those of the solves on t at load time (where following starts: its
    /// load-time root); which of those params are its unknowns and which
    /// each solve on t it reads sets, with that one's key. None: it or one
    /// it reads has no equations.
    fn key(&self, i: usize, params: &[f64]) -> Option<u128> {
        let s = &self.solves[i];
        let tapes: Vec<&Tape> = s.equations.as_ref()?.iter().collect();
        let mut at = params.to_vec();
        s.chain.iter().for_each(|&k| at[k] = self.load[k]);
        let mut fp = Fingerprint::default();
        fp.word(s.followed as u64);
        let names = cache::question(&mut fp, &tapes, &Ctx::new(0.0, at))?;
        let place = |fp: &mut Fingerprint, j: usize| self.solves[j].p.owns().iter().for_each(|&k| _ = fp.word(names.param(k as u32) as u64));
        place(&mut fp, i);
        for &j in &s.upstream {
            fp.key(self.key(j, params)?);
            place(&mut fp, j);
        }
        Some(cache::key(fp.value()))
    }

    /// Solve `i`'s roots followed over `span` (a moment, or a span on one
    /// side of 0): its root at the span's middle, and an enclosure of
    /// every one over the span (each step's tube narrowed to it by K).
    fn over(&self, i: usize, params: &[f64], span: (f64, f64)) -> Result<(Vec<f64>, Vec<I>), String> {
        if span == (0.0, 0.0) {
            return self.anchor(i, params, &mut self.path(i, params).lock().unwrap());
        }
        if let Some(f) = self.folded(i, params, span) {
            // Each piece's roots moved by its whole periods; the root at the
            // middle, as one of the first piece's.
            let mut x: Option<Vec<I>> = None;
            let mut r = Vec::new();
            for &(near, m) in &f.pieces {
                let (at, e) = self.over(i, params, (near.lo, near.hi))?;
                let e: Vec<I> = e.iter().zip(f.moved(m)).map(|(e, d)| add_o(*e, d)).collect();
                if r.is_empty() {
                    r = at.iter().zip(f.moved(I::point(m.lo))).map(|(v, d)| v + d.mid()).collect();
                }
                x = Some(x.map_or(e.clone(), |x| hull(&x, &e)));
            }
            return Ok((r, x.expect("a folded span has pieces")));
        }
        let steps = self.steps(i, params, span)?;
        let mid = span_mid(span);
        let (mut root, mut all): (Option<Vec<f64>>, Option<Vec<I>>) = (None, None);
        for s in steps {
            let sub = (span.0.max(s.a.min(s.b)), span.1.min(s.a.max(s.b)));
            let w = (span_mid(sub) - s.a) / (s.b - s.a);
            let (mut ctx, _) = self.context(i, params, sub, false)?;
            for ((&k, a), b) in self.solves[i].p.owns().iter().zip(&s.ra).zip(&s.rb) {
                ctx.params[k] = a + w * (b - a);
            }
            let (r, x) = self.solves[i].within(&ctx, &s.tube)?;
            if root.is_none() && sub.0 <= mid && mid <= sub.1 {
                root = Some(r);
            }
            all = Some(all.map_or(x.clone(), |a| hull(&a, &x)));
        }
        Ok((root.expect("a step holds the middle"), all.expect("steps hold the span")))
    }

    /// The steps of solve `i`'s path holding `span` (not at 0), in order
    /// away from 0, following it there first.
    fn steps(&self, i: usize, params: &[f64], span: (f64, f64)) -> Result<Vec<Step>, String> {
        let path = self.path(i, params);
        let mut path = path.lock().unwrap();
        let far = if span.1 > 0.0 { span.1 } else { span.0 };
        let sign = if far < 0.0 { -1.0 } else { 1.0 };
        let (near, to) = if far < 0.0 { (-span.1, -span.0) } else { span };
        loop {
            self.reach(i, params, &mut path, far)?;
            // A chunk gone from disk is followed again.
            if let Some(steps) = path.ways[(sign < 0.0) as usize].holding(sign, near, to) {
                return Ok(steps);
            }
        }
    }

    /// Solve `i`'s root at t = 0 (from its params at load time, as at load)
    /// and its enclosure: where following starts.
    fn anchor(&self, i: usize, params: &[f64], path: &mut Path) -> Result<(Vec<f64>, Vec<I>), String> {
        if let Some(a) = &path.anchor {
            return a.clone();
        }
        let s = &self.solves[i];
        let found = (|| {
            let (mut ctx, _) = self.context(i, params, (0.0, 0.0), false)?;
            for &k in s.p.owns() {
                ctx.params[k] = self.load[k];
            }
            s.p.provide(&mut ctx.params, 0.0).map_err(|e| format!("it has no root at t = 0 to start from ({e}). Make t = 0 a moment where it solves"))?;
            s.tube(&ctx, None, &[]).map_err(|e| format!("its root at t = 0 is not regular ({e}). Start the motion where it is"))
        })();
        path.anchor = Some(found.clone());
        found
    }

    /// Solve `i`'s period, looked for once a path: its equations read t
    /// through one expression of t alone (period::Clock: t itself, or a
    /// crank's angle 10 t + t^2 / 1000) and repeat in it, as ink_core::period
    /// shows, with the unknowns held; or the clock of the solves on t they
    /// read, each of which repeats in it (Follow::repeats). The least whole
    /// multiple of their periods; the clock at t = 0, and a moment by
    /// which it has moved more than that from there, and which way.
    /// None: none shown.
    fn period(&self, i: usize, params: &[f64], path: &mut Path) -> Option<Period> {
        if let Some(p) = &path.period {
            return p.clone();
        }
        let s = &self.solves[i];
        let found = (|| {
            let equations: Vec<&Tape> = s.equations.as_ref()?.iter().collect();
            let ctx = self.at_rest(i, params);
            let mut ups: Vec<(usize, Period)> = Vec::new();
            for &j in &s.upstream {
                ups.push((j, self.repeats(j, params)?.0));
            }
            let (clock, own) = match ups.first() {
                Some((_, u)) if ups.iter().all(|(_, v)| v.clock.same(&u.clock)) => {
                    let c = u.clock.in_tapes(&equations, &ctx)?;
                    let own = period_in(&equations, &c, &ctx)?;
                    (c, own)
                }
                Some(_) => return None,
                None => clocks(&equations, &ctx).into_iter().find_map(|c| match period_in(&equations, &c, &ctx)? {
                    Shown::Base(b) => Some((c, Shown::Base(b))),
                    Shown::Every => None,
                })?,
            };
            let mut bases: Vec<Coef> = ups.iter().map(|(_, u)| u.base).collect();
            if let Shown::Base(b) = own {
                bases.push(b);
            }
            let base = common_base(&bases)?;
            let p = base.period()?;
            let pe = p.value();
            let upstream = ups.iter().map(|(j, u)| Some((*j, p.whole_of(u.p)?))).collect::<Option<Vec<_>>>()?;
            // The first moment, doubling, then halving the gap, by which the
            // clock has moved past p.
            let g0 = self.clock(&clock, params, (0.0, 0.0));
            let past = |t: f64| {
                let g = self.clock(&clock, params, (t, t));
                if g.lo > interval::add_up(g0.hi, pe.hi) {
                    Some(1.0)
                } else if -g.hi > interval::add_up(-g0.lo, pe.hi) {
                    Some(-1.0)
                } else {
                    None
                }
            };
            let mut reach = FIRST;
            while past(reach).is_none() {
                reach *= 2.0;
                if !(pe.lo > 0.0 && pe.hi.is_finite() && reach < f64::MAX / 4.0) {
                    return None;
                }
            }
            let way = past(reach)?;
            let mut below = 0.0;
            for _ in 0..64 {
                let m = 0.5 * (below + reach);
                match past(m) == Some(way) {
                    true => reach = m,
                    false => below = m,
                }
            }
            // Its value over it does not jump.
            self.continuous_clock(&clock, params, (0.0, reach)).then_some(Period { clock, base, p, pe, upstream, g0, reach, way })
        })();
        // Not kept where a solve it reads stopped at the time limit: a
        // later call may find it.
        if found.is_some() || !ink_kernel::stop::reached() {
            path.period = Some(found.clone());
        }
        found
    }

    /// Solve `j`'s period and fold, where it is followed and both hold.
    fn repeats(&self, j: usize, params: &[f64]) -> Option<(Period, Fold)> {
        if !self.solves[j].followed {
            return None;
        }
        let path = self.path(j, params);
        let mut path = path.lock().unwrap();
        let per = self.period(j, params, &mut path)?;
        let fold = self.fold(j, params, &mut path, &per)?;
        Some((per, fold))
    }

    /// A clock's value over `span` of t, enclosed.
    fn clock(&self, c: &Clock, params: &[f64], span: (f64, f64)) -> I {
        let mut ctx = Ctx::new(span_mid(span), params.to_vec());
        ctx.t_span = (span.0 < span.1).then_some(span);
        eval::eval_interval(&c.cone, &ctx, [I::point(0.0); 3], None)
    }

    /// Whether a clock's rate along t over `span` is bounded: it does not
    /// jump there.
    fn continuous_clock(&self, c: &Clock, params: &[f64], span: (f64, f64)) -> bool {
        let mut ctx = Ctx::new(span_mid(span), params.to_vec());
        ctx.t_span = (span.0 < span.1).then_some(span);
        let (v, d) = eval::eval_interval_rates(&c.cone, &ctx, [I::point(0.0); 3], &eval::Inputs { t: true, params: Vec::new() });
        [v.lo, v.hi].iter().chain(d.iter().flat_map(|d| [&d.lo, &d.hi])).all(|x| x.is_finite())
    }

    /// The moments from 0 to `reach` at which a clock may be in `target`,
    /// enclosed (t: those of the target; else halving the span, keeping the
    /// halves it may be in there). None: none.
    fn preimage(&self, c: &Clock, params: &[f64], target: I, reach: f64) -> Option<I> {
        if c.is_t() {
            let at = interval::meet(target, I::new(0.0, reach));
            return (at.lo <= at.hi).then_some(at);
        }
        let meets = |span: (f64, f64)| {
            let g = self.clock(c, params, span);
            g.lo <= target.hi && target.lo <= g.hi
        };
        let mut front = vec![(0.0, reach)];
        for _ in 0..64 {
            let mut next = Vec::new();
            for &(a, b) in &front {
                let m = 0.5 * (a + b);
                if !(a < m && m < b) {
                    next.push((a, b));
                    continue;
                }
                next.extend([(a, m), (m, b)].into_iter().filter(|&h| meets(h)));
            }
            if next.is_empty() {
                return None;
            }
            if next.len() > 64 || next == front {
                break;
            }
            front = next;
        }
        Some(front.iter().fold(I::new(front[0].0, front[0].1), |h, &(a, b)| interval::hull(h, I::new(a, b))))
    }

    /// The context solve `i`'s equations are evaluated in over a step, for
    /// ink_core::period: its params' values, its unknowns and the params of
    /// the solves on t it reads ranged (as Krawczyk's test ranges them:
    /// Roots::enclose), so none is taken as a fixed float.
    fn at_rest(&self, i: usize, params: &[f64]) -> Ctx {
        let s = &self.solves[i];
        let mut ctx = Ctx::new(0.0, params.to_vec());
        let ranged = s.p.owns().iter().chain(s.upstream.iter().flat_map(|&j| self.solves[j].p.owns()));
        ctx.p_spans = ParamSpans::new(ranged.map(|&k| (k as u32, (params[k], params[k]))));
        ctx
    }

    /// Whether solve `i`'s period holds, looked for once a path: the path is
    /// followed past where the clock has moved by p, and there each unknown
    /// has moved by an exact amount d (a whole number of its own period, of
    /// p, or of what a solve it reads moves by: w = 2 th moves by twice th's
    /// turn; 0), its move rounded to one, such that the equations are
    /// unchanged by the clock moved by p, each solve's they read by its
    /// whole periods' moves, and the unknowns by d, at once
    /// (period::keeps). Where the root at 0 so moved (its enclosure, d
    /// enclosed) lies in the tube of every step that may hold that
    /// moment, which holds one root there, it is the root followed there;
    /// then the root followed to where the clock is g + p is the one where it
    /// is g moved, for every g, both ways: the moved path is a path of roots
    /// through it, and a regular root is followed one way. The unknowns'
    /// shifts; None: it does not hold.
    fn fold(&self, i: usize, params: &[f64], path: &mut Path, per: &Period) -> Option<Fold> {
        if let Some(f) = &path.fold {
            return f.clone();
        }
        let s = &self.solves[i];
        // Not kept where following stopped at the time limit or a chunk
        // is gone from disk (Err): a later call may find it.
        let found = (|| -> Result<Option<Fold>, ()> {
            let Some(equations) = s.equations.as_ref() else { return Ok(None) };
            let equations: Vec<&Tape> = equations.iter().collect();
            let Ok((r0, e0)) = self.anchor(i, params, path) else { return Ok(None) };
            if self.reach(i, params, path, per.reach).is_err() {
                return if ink_kernel::stop::reached() { Err(()) } else { Ok(None) };
            }
            let target = add_o(per.g0, mul_o(I::point(per.way), per.pe));
            let Some(near) = self.preimage(&per.clock, params, target, per.reach) else { return Ok(None) };
            let steps = path.ways[0].holding(1.0, near.lo, near.hi).ok_or(())?;
            let Some(at) = steps.iter().find(|st| st.a <= near.mid() && near.mid() <= st.b) else { return Ok(None) };
            let w = (near.mid() - at.a) / (at.b - at.a);
            let ctx = self.at_rest(i, params);
            let Some(by) = per.p.times(per.way as i128) else { return Ok(None) };
            // What the solves on t read move by over p.
            let mut ups: Vec<(u32, Coef)> = Vec::new();
            for &(j, n) in &per.upstream {
                let Some((_, f)) = self.repeats(j, params) else { return Ok(None) };
                for (&k, d) in self.solves[j].p.owns().iter().zip(&f.exact) {
                    if let Some(d) = d {
                        let Some(d) = d.times(n) else { return Ok(None) };
                        ups.push((k as u32, d));
                    }
                }
            }
            // Each unknown's candidate moves: 0, and its move rounded to a
            // whole number of each unit.
            let mut choices: Vec<Vec<Option<Coef>>> = Vec::new();
            for (j, &k) in s.p.owns().iter().enumerate() {
                let moved = at.ra[j] + w * (at.rb[j] - at.ra[j]) - r0[j];
                let own = match period(&equations, Var::Param(k as u32), &ctx) {
                    Some(Shown::Base(c)) => c.period(),
                    _ => None,
                };
                let mut c: Vec<Option<Coef>> = vec![None];
                for u in own.into_iter().chain([per.p]).chain(ups.iter().map(|u| u.1)) {
                    let n = (moved / u.value().mid()).round();
                    if n != 0.0
                        && n.abs() < 1e15
                        && let Some(d) = u.times(n as i128)
                        && !c.contains(&Some(d))
                    {
                        c.push(Some(d));
                    }
                }
                choices.push(c);
            }
            let lands = |d: &[I], tube: &[I]| e0.iter().zip(d).zip(tube).all(|((e, d), x)| x.lo <= add_o(*e, *d).lo && add_o(*e, *d).hi <= x.hi);
            let mut pick = vec![0usize; choices.len()];
            for _ in 0..256 {
                let exact: Vec<Option<Coef>> = pick.iter().zip(&choices).map(|(&n, c)| c[n]).collect();
                let shift: Vec<I> = exact.iter().map(|d| d.map_or(I::point(0.0), |d| d.value())).collect();
                if steps.iter().all(|st| lands(&shift, &st.tube)) {
                    let mut moves = ups.clone();
                    moves.extend(s.p.owns().iter().zip(&exact).filter_map(|(&k, d)| Some((k as u32, (*d)?))));
                    if keeps(&equations, Some((&per.clock, by)), &moves, &ctx) {
                        return Ok(Some(Fold { exact, shift }));
                    }
                }
                // The next choice of each unknown's, in turn.
                let Some(j) = (0..pick.len()).find(|&j| pick[j] + 1 < choices[j].len()) else { return Ok(None) };
                pick[j] += 1;
                pick[..j].iter_mut().for_each(|n| *n = 0);
            }
            Ok(None)
        })();
        let found = found.ok()?;
        path.fold = Some(found.clone());
        found
    }

    /// `span` (not 0) folded into solve `i`'s first period, where it has
    /// one that holds and the span lies farther than it from 0, either
    /// way, its clock not jumping between: the moments from 0 to the
    /// period's reach where the clock is where it is over the span, moved
    /// by whole periods. None: not folded.
    fn folded(&self, i: usize, params: &[f64], span: (f64, f64)) -> Option<Folded> {
        if !self.solves[i].followed {
            return None;
        }
        let path = self.path(i, params);
        let mut path = path.lock().unwrap();
        let per = self.period(i, params, &mut path)?;
        // Folding is a shortcut: it follows the path to the period's
        // reach, so only a span farther from 0 than that, which following
        // straight there would pass it for, is folded (a moment near 0
        // never pays for a long period, t = -0.1 against a period of 1e17
        // s).
        if span.0.abs().min(span.1.abs()) <= per.reach || (span.0 < 0.0) != (span.1 < 0.0) {
            return None;
        }
        let fold = self.fold(i, params, &mut path, &per)?;
        if !self.continuous_clock(&per.clock, params, (span.0.min(0.0), span.1.max(0.0))) {
            return None;
        }
        // The whole numbers m with the clock moved back by m periods in the
        // first period, for some moment of the span and p in its enclosure.
        let g = self.clock(&per.clock, params, span);
        let along = |x: I| mul_o(I::point(per.way), sub_o(x, per.g0));
        let q = interval::bin(BinOp::Div, along(g), per.pe).0;
        let (m0, m1) = (q.lo.floor(), q.hi.floor());
        if !(m0.is_finite() && m1.is_finite()) {
            return None;
        }
        let first = interval::hull(per.g0, add_o(per.g0, mul_o(I::point(per.way), per.pe)));
        let back = |m: I| {
            let at = interval::meet(sub_o(g, mul_o(I::point(per.way), mul_o(m, per.pe))), first);
            (at.lo <= at.hi).then(|| self.preimage(&per.clock, params, at, per.reach)).flatten()
        };
        // Two counts, a piece each; more, and one between them covers the
        // whole period: one piece, every count.
        let pieces: Vec<(Option<I>, I)> = match m1 - m0 {
            0.0 => vec![(back(I::point(m0)), I::point(m0))],
            1.0 => vec![(back(I::point(m0)), I::point(m0)), (back(I::point(m1)), I::point(m1))],
            _ => vec![(Some(I::new(0.0, per.reach)), I::new(m0, m1))],
        };
        let pieces: Vec<(I, I)> = pieces.into_iter().filter_map(|(near, m)| Some((near?, m))).collect();
        (!pieces.is_empty()).then_some(Folded { pieces, shift: fold.shift })
    }

    /// Solve `i`'s params at moment t, folded (Follow::folded): its
    /// provider's solve from a piece's root moved, where it lands in every
    /// tube that may hold the moment t is folded to, moved by its piece's
    /// shifts (rounded in: whichever holds it holds one root there); else
    /// the one root at t Krawczyk's test shows in a box holding them all,
    /// moved (rounded out). A count of periods too uncertain to tell which
    /// turn an unknown is on leaves that box more than one root: refused.
    fn provide_folded(&self, i: usize, params: &mut [f64], t: f64, f: &Folded) -> Result<Option<String>, String> {
        let s = &self.solves[i];
        let mut held: Option<Vec<I>> = None;
        let mut tubes: Vec<(Vec<I>, Vec<I>)> = Vec::new();
        for &(near, m) in &f.pieces {
            let moved = f.moved(m);
            let steps = self.steps(i, params, (near.lo, near.hi))?;
            if tubes.is_empty() {
                let at = &steps[0];
                let w = ((near.mid() - at.a) / (at.b - at.a)).clamp(0.0, 1.0);
                for (j, &k) in s.p.owns().iter().enumerate() {
                    params[k] = at.ra[j] + w * (at.rb[j] - at.ra[j]) + moved[j].mid();
                }
            }
            for st in steps {
                let e: Vec<I> = st.tube.iter().zip(&moved).map(|(x, d)| add_o(*x, *d)).collect();
                held = Some(held.map_or(e.clone(), |h| hull(&h, &e)));
                tubes.push((st.tube, moved.clone()));
            }
        }
        let solved = s.p.provide(params, t);
        let down = |a: f64, b: f64| -interval::add_up(-a, -b);
        let within =
            |tube: &[I], moved: &[I], p: &[f64]| s.p.owns().iter().zip(tube).zip(moved).all(|((&k, x), d)| interval::add_up(x.lo, d.hi) <= p[k] && p[k] <= down(x.hi, d.lo));
        if solved.is_ok() && tubes.iter().all(|(tube, moved)| within(tube, moved, params)) {
            return solved;
        }
        let held = held.expect("a folded moment has steps");
        let (mut ctx, _) = self.context(i, params, (t, t), false)?;
        for (&k, x) in s.p.owns().iter().zip(&held) {
            ctx.params[k] = x.mid();
        }
        let counts = f.pieces.iter().fold(I::point(f.pieces[0].1.lo), |h, (_, m)| interval::hull(h, *m));
        let (r, _) = Roots::enclose(s.p.owns(), s.equations.clone().expect("a followed solve has equations"), &ctx, Some(&held), &[]).map_err(|e| match counts.lo < counts.hi {
            true => format!("its path repeats, but t is {} to {} periods on, and over those the root is not shown to be one ({e})", counts.lo, counts.hi),
            false => e,
        })?;
        for (&k, v) in s.p.owns().iter().zip(r) {
            params[k] = v;
        }
        Ok(None)
    }

    /// Follow solve `i`'s path to moment t (not 0). Err: why it is not
    /// followed there.
    fn reach(&self, i: usize, params: &[f64], path: &mut Path, t: f64) -> Result<(), String> {
        let (r0, e0) = self.anchor(i, params, path)?;
        let w = &mut path.ways[(t < 0.0) as usize];
        let sign = if t < 0.0 { -1.0 } else { 1.0 };
        let mut filled = false;
        let reached = loop {
            if w.chunks.last().is_some_and(|c| sign * t <= sign * c.0) {
                break Ok(());
            }
            if let Some(end) = &w.end {
                break Err(end.clone());
            }
            if ink_kernel::stop::reached() {
                let at = w.chunks.last().map_or(0.0, |c| c.0);
                break Err(format!("following stopped at t = {at}, where the call's time limit passed"));
            }
            let (a, ra, ea, mut h) = match w.last() {
                Some(s) => (s.b, s.rb.clone(), s.eb.clone(), ((if s.halved { 1.0 } else { 2.0 }) * (s.b - s.a).abs())),
                None => (0.0, r0.clone(), e0.clone(), FIRST),
            };
            let mut halved = false;
            loop {
                match self.step(i, params, a, a + sign * h, &ra, &ea) {
                    Ok(s) => {
                        filled |= w.push(Step { halved, ..s });
                        break;
                    }
                    Err(why) if h <= SHORTEST => {
                        w.end = Some(format!("past t = {a}, {why}. Which root it moves to after that is not known"));
                        break;
                    }
                    Err(_) => (h, halved) = (h / 2.0, true),
                }
            }
        };
        if filled {
            w.keep();
        }
        reached
    }

    /// Solve `i` over the step from a to b, from its root ra at a
    /// (enclosed by ea). Err: why it is not certified.
    fn step(&self, i: usize, params: &[f64], a: f64, b: f64, ra: &[f64], ea: &[I]) -> Result<Step, String> {
        let s = &self.solves[i];
        let (mut ctx, moving) = self.context(i, params, (a.min(b), a.max(b)), true)?;
        for (&k, &v) in s.p.owns().iter().zip(ra) {
            ctx.params[k] = v;
        }
        let (_, tube) = s.tube(&ctx, Some(ea), &moving)?;
        // One root in the box at every moment is the root followed, moved
        // along t, only where the equations do not jump.
        let equations = s.equations.as_deref().expect("a followed solve has equations");
        if !Roots::continuous(s.p.owns(), equations, &ctx, &tube, &moving) {
            return Err("its equations switch there (a select, floor or atan2 reading t), and their rates along t over the step are not bounded, \
                        so they are not shown continuous (a jump, or a seam t moves even where its sides agree)"
                .into());
        }
        let (mut at, _) = self.context(i, params, (b, b), false)?;
        for (&k, &v) in s.p.owns().iter().zip(ra) {
            at.params[k] = v;
        }
        let (rb, eb) = s.within(&at, &tube)?;
        CERTIFIED.with(|n| n.set(n.get() + 1));
        Ok(Step { a, b, tube, ra: ra.to_vec(), rb, eb, halved: false })
    }

    /// What solve `j`, on t, sets over `span`: its value at the span's
    /// middle and its enclosure. Followed, narrowed from its steps
    /// (Follow::over); else (a solve's rates: its equations linear in them,
    /// so one solution where their coefficients are regular) Krawczyk's
    /// test as for a root, from its params' values at load time. Where
    /// `moving` and over a span of t, also its params moving along it
    /// (roots::Moving): each over its enclosure at rate ENTIRE where its
    /// rate is not bounded, so a rate that reads it is not bounded either
    /// (Roots::rate) and a path bound over it falls back to the plain one.
    /// Err: why not, as what reads it says.
    fn sets(&self, j: usize, params: &[f64], span: (f64, f64), moving: bool) -> Result<(Vec<f64>, Vec<I>, Vec<Moving>), String> {
        let s = &self.solves[j];
        let moving = moving && span.0 < span.1;
        let unbounded = |e: &[I]| -> Vec<Moving> { s.p.owns().iter().zip(e).map(|(&k, &x)| (k as u32, x, I::ENTIRE)).collect() };
        let (r, e, mut ctx, along) = match s.followed {
            true => {
                let (r, e) = self.over(j, params, span).map_err(|m| format!("it reads what {} sets, which is not followed there: {m}", s.label))?;
                if !moving {
                    return Ok((r, e, Vec::new()));
                }
                match self.context(j, params, span, true) {
                    Ok((ctx, along)) => (r, e, ctx, along),
                    Err(_) => {
                        let m = unbounded(&e);
                        return Ok((r, e, m));
                    }
                }
            }
            false => {
                let not = |m: String| format!("it reads what {} sets, which is not bounded there: {m}", s.label);
                // A provider without equations that bounds what it sets.
                if s.equations.is_none()
                    && let Some(over) = s.p.over_span(params, span)
                {
                    let (r, e, rates) = over.map_err(not)?;
                    let m = if moving { s.p.owns().iter().zip(&e).zip(rates).map(|((&k, &x), d)| (k as u32, x, d)).collect() } else { Vec::new() };
                    return Ok((r, e, m));
                }
                let (mut ctx, along) = self.context(j, params, span, moving).map_err(not)?;
                for &k in s.p.owns() {
                    ctx.params[k] = self.load[k];
                }
                let (r, e) = s.tube(&ctx, None, &along).map_err(not)?;
                if !moving {
                    return Ok((r, e, Vec::new()));
                }
                (r, e, ctx, along)
            }
        };
        let Ok((_, at, _)) = self.sets(j, params, (span_mid(span), span_mid(span)), false) else {
            let m = unbounded(&e);
            return Ok((r, e, m));
        };
        let spans: Vec<(u32, (f64, f64))> = s.p.owns().iter().zip(&e).map(|(&k, x)| (k as u32, (x.lo, x.hi))).collect();
        ctx.p_spans = ctx.p_spans.with_derived(&spans);
        s.p.owns().iter().zip(&r).for_each(|(&k, &v)| ctx.params[k] = v);
        let rate = s.equations.as_deref().and_then(|eq| Roots::rate(s.p.owns(), eq, &ctx, &along));
        let moved = match rate {
            Some(d) => s.p.owns().iter().zip(at).zip(d).map(|((&k, x), d)| (k as u32, x, d)).collect(),
            None => unbounded(&at),
        };
        Ok((r, e, moved))
    }

    /// What solve `i`, on t, sets over `span` of t, either side of 0
    /// (Follow::sets, the two sides' hull where it holds 0): its value at
    /// the span's middle and an enclosure over the span. Err: why not.
    pub fn over_t(&self, i: usize, params: &[f64], span: (f64, f64)) -> Result<(Vec<f64>, Vec<I>), String> {
        if !(span.0 < 0.0 && 0.0 < span.1) {
            return self.sets(i, params, span, false).map(|(r, e, _)| (r, e));
        }
        let (_, a, _) = self.sets(i, params, (span.0, 0.0), false)?;
        let (_, b, _) = self.sets(i, params, (0.0, span.1), false)?;
        let (r, _, _) = self.sets(i, params, (span_mid(span), span_mid(span)), false)?;
        Ok((r, hull(&a, &b)))
    }

    /// What solve `i`, on t, sets over `span` of t, as it moves along it:
    /// for each param, an enclosure at the span's middle moment and bounds
    /// on its rate along t over the span (Follow::sets; either side of 0,
    /// the two sides' rates' hull). Err: why not.
    pub fn moving_over_t(&self, i: usize, params: &[f64], span: (f64, f64)) -> Result<Vec<(I, I)>, String> {
        let m = span_mid(span);
        let (_, at, _) = self.sets(i, params, (m, m), false)?;
        let sides: Vec<(f64, f64)> = if span.0 < 0.0 && 0.0 < span.1 { vec![(span.0, 0.0), (0.0, span.1)] } else { vec![span] };
        let mut rates: Option<Vec<I>> = None;
        for side in sides {
            let (_, _, moving) = self.sets(i, params, side, true)?;
            let r: Vec<I> = moving.iter().map(|m| m.2).collect();
            rates = Some(rates.map_or(r.clone(), |h| hull(&h, &r)));
        }
        let rates = rates.expect("a side");
        if rates.len() != at.len() {
            return Err("its rates along t are not bounded".into());
        }
        Ok(at.into_iter().zip(rates).collect())
    }

    /// The params solve `i`'s path rests on that no solve on t sets.
    pub fn rests_on(&self, i: usize) -> &[usize] {
        &self.solves[i].fixed
    }

    /// Whether what solve `i` sets over a span of t (Follow::over_t) is
    /// certified: it and every solve on t under it have equations, so each
    /// enclosure is a root's or a rate's. A provider's own bounds (a
    /// simulation's: its interpolant's) are measured.
    pub fn certifies(&self, i: usize) -> bool {
        let s = &self.solves[i];
        s.equations.is_some() && s.upstream.iter().all(|&j| self.certifies(j))
    }

    /// A context for solve `i` over `span`: t over it (at it, a moment),
    /// and the solves on t `i` reads followed over it, each at its root at
    /// the span's middle and ranged over its enclosure (Follow::sets); with
    /// `moving`, their params moving along t over it.
    fn context(&self, i: usize, params: &[f64], span: (f64, f64), moving: bool) -> Result<(Ctx, Vec<Moving>), String> {
        let mut ctx = Ctx::new(span_mid(span), params.to_vec());
        ctx.t_span = (span.0 < span.1).then_some(span);
        let (mut spans, mut along) = (Vec::new(), Vec::new());
        for &j in &self.solves[i].upstream {
            let (r, e, m) = self.sets(j, params, span, moving)?;
            for ((&k, v), x) in self.solves[j].p.owns().iter().zip(r).zip(e) {
                ctx.params[k] = v;
                spans.push((k as u32, (x.lo, x.hi)));
            }
            along.extend(m);
        }
        ctx.p_spans = ParamSpans::new(spans);
        Ok((ctx, along))
    }
}

impl Solve {
    /// Its root about ctx's values and a box holding one root of its
    /// equations at every moment and value of ctx's spans, and `from`
    /// (Roots::enclose), reading params `moving` along t.
    fn tube(&self, ctx: &Ctx, from: Option<&[I]>, moving: &[Moving]) -> Result<(Vec<f64>, Vec<I>), String> {
        let equations = self.equations.clone().ok_or("it is measured at each moment, with no bound over a span of t")?;
        Roots::enclose(self.p.owns(), equations, ctx, from, moving)
    }

    /// The one root in `tube` (from `tube` over spans holding ctx's) at
    /// ctx: its value and enclosure (Roots::inside).
    fn within(&self, ctx: &Ctx, tube: &[I]) -> Result<(Vec<f64>, Vec<I>), String> {
        Roots::inside(self.p.owns(), self.equations.clone().expect("a followed solve has equations"), ctx, tube)
    }
}

impl Way {
    /// The way keyed `disk` with n unknowns, holding the chunks kept on
    /// disk: their far ends, their steps read when needed.
    fn kept(disk: Option<u128>, n: usize) -> Way {
        let ends: Vec<u64> = disk.and_then(|d| cache::load(KIND, d)).unwrap_or_default();
        Way { disk, n, chunks: ends.into_iter().map(|b| (f64::from_bits(b), None)).collect(), end: None }
    }

    /// Chunk c's key on disk.
    fn chunk_key(disk: u128, c: usize) -> u128 {
        Fingerprint::default().key(disk).word(c as u64).value()
    }

    /// Chunk c's steps, read from disk where needed; None where it is not
    /// there (evicted): then it and the chunks after it are dropped, to be
    /// followed again.
    fn read(&mut self, c: usize) -> Option<&[Step]> {
        if self.chunks[c].1.is_none() {
            let n = self.n;
            let words = self.disk.and_then(|d| cache::load::<Vec<u64>>(KIND, Way::chunk_key(d, c)));
            match words.filter(|w| w.len() == CHUNK * Step::words(n)) {
                Some(w) => self.chunks[c].1 = Some(w.chunks(Step::words(n)).map(|s| Step::of(s, n)).collect()),
                None => {
                    self.chunks.truncate(c);
                    return None;
                }
            }
        }
        self.chunks[c].1.as_deref()
    }

    /// The last step, read from disk where needed.
    fn last(&mut self) -> Option<&Step> {
        let mut c = self.chunks.len();
        while c > 0 && self.read(c - 1).is_none() {
            c -= 1;
        }
        self.chunks.last()?.1.as_ref()?.last()
    }

    /// Add `s` after the last step (read: Way::last); whether it filled a
    /// chunk, which is then kept on disk.
    fn push(&mut self, s: Step) -> bool {
        let b = s.b;
        match self.chunks.last_mut() {
            Some((end, Some(steps))) if steps.len() < CHUNK => {
                *end = b;
                steps.push(s);
            }
            _ => self.chunks.push((b, Some(vec![s]))),
        }
        let c = self.chunks.len() - 1;
        let steps = self.chunks[c].1.as_ref().expect("just added");
        if steps.len() < CHUNK {
            return false;
        }
        if let Some(d) = self.disk {
            let mut words = Vec::with_capacity(CHUNK * Step::words(self.n));
            steps.iter().for_each(|s| s.bits(&mut words));
            cache::store(KIND, Way::chunk_key(d, c), &words);
        }
        true
    }

    /// Keep on disk the far ends of the full chunks, which are there.
    fn keep(&self) {
        let full = self.chunks.iter().take_while(|c| c.1.as_ref().is_none_or(|s| s.len() == CHUNK));
        if let Some(d) = self.disk {
            cache::store(KIND, d, &full.map(|c| c.0.to_bits()).collect::<Vec<u64>>());
        }
    }

    /// Its steps holding the moments from `near` to `far` away from 0
    /// (`sign` the way's): the first ending at or past `near`, then those
    /// starting before `far`. None where a chunk is gone from disk.
    fn holding(&mut self, sign: f64, near: f64, far: f64) -> Option<Vec<Step>> {
        let first = self.chunks.partition_point(|c| sign * c.0 < near);
        let mut out = Vec::new();
        for c in first..self.chunks.len() {
            if c > 0 && sign * self.chunks[c - 1].0 > far {
                break;
            }
            let steps = self.read(c)?;
            let i = steps.partition_point(|s| sign * s.b < near);
            out.extend(steps[i..].iter().take_while(|s| sign * s.a <= far).cloned());
        }
        Some(out)
    }
}

impl Step {
    /// The words a step of n unknowns takes on disk.
    const fn words(n: usize) -> usize {
        3 + 6 * n
    }

    /// Its words on disk: its numbers' bits, exactly.
    fn bits(&self, out: &mut Vec<u64>) {
        out.extend([self.a.to_bits(), self.b.to_bits(), self.halved as u64]);
        out.extend(self.tube.iter().chain(&self.eb).flat_map(|x| [x.lo.to_bits(), x.hi.to_bits()]));
        out.extend(self.ra.iter().chain(&self.rb).map(|v| v.to_bits()));
    }

    /// The step of n unknowns in `w` (Step::bits).
    fn of(w: &[u64], n: usize) -> Step {
        let f = |i: usize| f64::from_bits(w[i]);
        let boxed = |at: usize| (0..n).map(|j| I { lo: f(at + 2 * j), hi: f(at + 2 * j + 1) }).collect();
        let point = |at: usize| (0..n).map(|j| f(at + j)).collect();
        Step { a: f(0), b: f(1), halved: w[2] != 0, tube: boxed(3), eb: boxed(3 + 2 * n), ra: point(3 + 4 * n), rb: point(3 + 5 * n) }
    }
}

/// The smallest box holding a and b.
fn hull(a: &[I], b: &[I]) -> Vec<I> {
    a.iter().zip(b).map(|(a, b)| I { lo: a.lo.min(b.lo), hi: a.hi.max(b.hi) }).collect()
}

/// What a solve on t sets, derived over t ranged as param `k` (a claim
/// over a motion: motion::ranged): over a span of k, its params' enclosure
/// over it from the steps followed (Follow::over_t), and at k's value its
/// root there. Made over the whole span `whole` (Err: why not), so every
/// piece of it is followed already; one whose enclosure cannot be narrowed
/// to a piece keeps the whole span's.
pub struct OverT {
    follow: Arc<Follow>,
    i: usize,
    k: [u32; 1],
    owns: Vec<u32>,
    params: Vec<f64>,
    whole: (f64, f64),
    all: Vec<I>,
    cache: Mutex<HashMap<(u64, u64), (Vec<f64>, Vec<I>)>>,
    /// Rate bounds found (`rate_bounds`), by span.
    moving: Mutex<HashMap<(u64, u64), Option<Vec<(I, Vec<I>)>>>>,
}

impl OverT {
    pub fn new(follow: Arc<Follow>, i: usize, k: u32, params: &[f64], whole: (f64, f64)) -> Result<OverT, String> {
        let (_, all) = follow.over_t(i, params, whole)?;
        let owns = follow.solves[i].p.owns().iter().map(|&u| u as u32).collect();
        Ok(OverT { follow, i, k: [k], owns, params: params.to_vec(), whole, all, cache: Mutex::new(HashMap::new()), moving: Mutex::new(HashMap::new()) })
    }

    /// Its params' enclosures over the whole span.
    pub fn bands(&self) -> Vec<(u32, (f64, f64))> {
        self.owns.iter().zip(&self.all).map(|(&u, x)| (u, (x.lo, x.hi))).collect()
    }

    /// Its value at the middle of `span` and enclosure over it (inside the
    /// whole span).
    fn at(&self, span: (f64, f64)) -> (Vec<f64>, Vec<I>) {
        let key = (span.0.to_bits(), span.1.to_bits());
        if let Some(hit) = self.cache.lock().unwrap().get(&key) {
            return hit.clone();
        }
        let inside = self.whole.0 <= span.0 && span.1 <= self.whole.1;
        let r = match inside.then(|| self.follow.over_t(self.i, &self.params, span)) {
            Some(Ok(r)) => r,
            Some(Err(_)) if span.0 < span.1 => (self.at((span_mid(span), span_mid(span))).0, self.all.clone()),
            Some(Err(_)) => (self.all.iter().map(|x| span_mid((x.lo, x.hi))).collect(), self.all.clone()),
            None => (vec![f64::NAN; self.owns.len()], vec![I::ENTIRE; self.owns.len()]),
        };
        let mut cache = self.cache.lock().unwrap();
        if cache.len() > 1 << 16 {
            cache.clear();
        }
        cache.insert(key, r.clone());
        r
    }
}

impl ink_core::eval::Derive for OverT {
    fn owns(&self) -> &[u32] {
        &self.owns
    }

    /// Its solve's key on disk (Follow::key: its equations and what its
    /// root rests on, which the steps it reads are kept by), t's param, and
    /// the span it is made over with the enclosure there. Its own caches
    /// hold what it derived, by the span: they change no answer.
    fn fingerprint(&self, fp: &mut Fingerprint, names: &ink_core::tape::Names) -> bool {
        // Every field, named: one added must be hashed too.
        let OverT { follow, i, k, owns, params, whole, all, cache: _, moving: _ } = self;
        let Some(key) = follow.key(*i, params) else { return false };
        fp.key(key).word(names.param(k[0]) as u64).word(owns.len() as u64);
        owns.iter().for_each(|&u| _ = fp.word(names.param(u) as u64));
        fp.f64(whole.0).f64(whole.1).word(all.len() as u64);
        all.iter().for_each(|x| _ = fp.f64(x.lo).f64(x.hi));
        true
    }

    fn reads(&self) -> &[u32] {
        &self.k
    }

    fn derive(&self, ctx: &mut Ctx) {
        let v = ctx.params[self.k[0] as usize];
        let span = ctx.p_spans.get(self.k[0]).unwrap_or((v, v));
        let (mut r, e) = self.at(span);
        if v != span_mid(span) {
            r = self.at((v, v)).0;
        }
        for ((&u, x), c) in self.owns.iter().zip(e).zip(r) {
            ctx.p_spans.set_span(u, (x.lo, x.hi));
            ctx.params[u as usize] = c.clamp(x.lo, x.hi);
        }
    }

    /// Its params read through t's param: their rates along t are bounded
    /// over a span (Follow::moving_over_t).
    fn chains(&self) -> bool {
        true
    }

    /// Over param k's span in ctx (inside the whole span): each param's
    /// enclosure at its middle moment and its rate along t over it, kept
    /// by the span.
    fn rate_bounds(&self, ctx: &Ctx) -> Option<Vec<(I, Vec<I>)>> {
        let span = ctx.p_spans.get(self.k[0])?;
        if !(self.whole.0 <= span.0 && span.1 <= self.whole.1 && span.0 < span.1) {
            return None;
        }
        let key = (span.0.to_bits(), span.1.to_bits());
        if let Some(hit) = self.moving.lock().unwrap().get(&key) {
            return hit.clone();
        }
        let r = self.follow.moving_over_t(self.i, &self.params, span).ok().map(|m| m.into_iter().map(|(x, d)| (x, vec![d])).collect());
        let mut cache = self.moving.lock().unwrap();
        if cache.len() > 1 << 16 {
            cache.clear();
        }
        cache.insert(key, r.clone());
        r
    }

    /// Along t, by the roots a little either side (floats: for witness
    /// searches only).
    fn rates(&self, ctx: &Ctx) -> Vec<Vec<f64>> {
        let v = ctx.params[self.k[0] as usize];
        let h = 1e-6 * (self.whole.1 - self.whole.0).max(1e-9);
        let (a, b) = ((v - h).max(self.whole.0), (v + h).min(self.whole.1));
        if b <= a {
            return vec![vec![0.0]; self.owns.len()];
        }
        let (ra, rb) = (self.at((a, a)).0, self.at((b, b)).0);
        ra.iter().zip(&rb).map(|(x, y)| vec![(y - x) / (b - a)]).collect()
    }
}
