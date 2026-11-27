use super::*;

pub(super) const DRIVE: &str = "searched over one period of the drive every read of t passes through (each shape turned with it repeating by its own turn symmetry): that period's angles stand for every moment";

/// A claim over a motion read through its drive: `node`, a value of t
/// (and params) that every read of t the claim makes passes through (a
/// motor's angle spinning up as t^2), `window`, one period of the claim
/// in it from the least value it takes over the span, and `slack`, how far
/// (a length) a shape may sit from where the window places it.
pub(super) struct Drive {
    pub(super) node: u32,
    pub(super) window: (f64, f64),
    pub(super) slack: f64,
    /// How far the drive travels over the span (its tests weigh the
    /// window by it).
    #[cfg_attr(not(test), expect(dead_code))]
    pub(super) travel: f64,
}

/// The claim's drive (`Drive`): the first of the values `tries` proposes
/// that `through` shows periodic. None: none shown, or the call's time
/// limit passed. It costs a reading of the claim per try.
pub(super) fn drive(v: &mut ModelView, roots: &[u32], span: (f64, f64), k: u32) -> Option<Drive> {
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
    // Every shape placed, moving or still: what a frame carried round with them turns (`carried_round`).
    let g = &v.doc.g;
    let every: Vec<(u32, crate::Placed)> = v.doc.placements().into_iter().filter(|(w, p)| reached.contains(w) && !g.depends_on(p.shape, |x| x == Var::T)).collect();
    let every: Vec<(u32, crate::Placed, ink_core::symmetry::Symmetry)> = every.into_iter().map(|(w, p)| (w, p, turn_symmetry(&v.doc.g, p.shape))).collect();
    for d in tries {
        if ink_kernel::stop::reached() {
            return None;
        }
        let found = match through(v, roots, &shapes, d, span, k) {
            Some(x) => Some(x),
            None => carried_round(v, roots, &every, d, span, k),
        };
        if let Some(x) = found {
            return (!ink_kernel::stop::reached()).then_some(x);
        }
    }
    None
}

/// `through`, where the shapes do not each turn in place but are carried round one axis together
/// (planets on their carrier about its sun and ring): the claim is about placed shapes alone, so it is
/// the same claim with all of them moved rigidly, and in the frame that goes round with them each only
/// turns about its own axis, less the frame's turn. Shown exactly: every turn is about the world's z
/// (its angle read off its first column, `rigid::carried`), each shape's place less the first's is 0
/// in x and y or a constant length carried round z by one angle th, the same rate for all (read the
/// same way), and no height moves. Each shape then stands for a value repeating whenever its turn less
/// th passes 2 pi / n of its own symmetry, and the claim's period in d is theirs together. None: not
/// so shown, or the claim reads a point of space besides the shapes placed.
fn carried_round(v: &mut ModelView, roots: &[u32], placed: &[(u32, crate::Placed, ink_core::symmetry::Symmetry)], d: u32, span: (f64, f64), k: u32) -> Option<Drive> {
    use ink_core::period::{Coef, Shown, form, period_kept};
    use ink_core::symmetry::Order;
    if placed.len() < 2 {
        return None;
    }
    let tn = v.doc.g.var(Var::T);
    let a = v.doc.g.var(Var::Param(k));
    let by: std::collections::HashMap<u32, u32> = [(d, a)].into();
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
    let is_a = |x: Var| x == Var::Param(k);
    // The rate of an angle read as windings: sum m_i c_i (None: it does not turn with d).
    let rate = |g: &ink_core::graph::Graph, wound: &[(u32, i32)]| -> Option<Option<Coef>> {
        wound.iter().try_fold(None, |sum: Option<Coef>, &(n, m)| match form(&Tape::new(g, n), Var::Param(k), &ctx)? {
            None => Some(sum),
            Some(c) => {
                let c = c.times(m as i128)?;
                match sum {
                    None => Some(Some(c)),
                    Some(s) => s.plus(c),
                }
            }
        })
    };
    // Each shape's turn about the world's z, and its place, read through d.
    let mut each: Vec<(u32, ink_core::symmetry::Symmetry, Option<Coef>, [u32; 3], u32)> = Vec::new();
    for &(w, p, sym) in placed {
        if ink_kernel::stop::reached() {
            return None;
        }
        let r = p.at.r.map(|row| row.map(|n| v.doc.g.replace(n, &by)));
        let t = p.at.t.map(|n| v.doc.g.replace(n, &by));
        let g = &v.doc.g;
        let is = |n: u32, x: f64| {
            let i = ink_core::eval::eval_interval(&Tape::new(g, n), &ctx, [interval::I::point(0.0); 3], None);
            i.lo == x && i.hi == x
        };
        if !(is(r[0][2], 0.0) && is(r[1][2], 0.0) && is(r[2][0], 0.0) && is(r[2][1], 0.0) && (is(r[2][2], 1.0) || is(r[2][2], -1.0))) {
            return None;
        }
        let moving = |n: u32| g.depends_on(n, is_a);
        // A turn about z, as own_turn reads one (facing up or down the axis: turned over, its heading
        // still goes round z, and less the frame's it is its own turn either way); its heading's
        // windings off its first column.
        ink_kernel::rigid::own_turn(g, &r, &ctx, moving)?;
        let turned = rate(g, &ink_kernel::rigid::carried(g, r[0][0], r[1][0], &ctx, moving)?)?;
        each.push((w, sym, turned, t, p.shape));
    }
    // The angle the others' places go round the first's by: one rate for all that are off its axis.
    let first = each[0].3;
    let mut round: Option<Coef> = None;
    for e in &each[1..] {
        let off: Vec<u32> = (0..3).map(|j| v.doc.g.sub(e.3[j], first[j])).collect();
        let off = ink_kernel::rigid::simplified(&mut v.doc.g, &off, &ctx).unwrap_or(off);
        let g = &v.doc.g;
        let moving = |n: u32| g.depends_on(n, is_a);
        if moving(off[2]) {
            return None;
        }
        if !moving(off[0]) && !moving(off[1]) {
            // Still beside the first: carried round with the rest only on its axis.
            let zero = |n: u32| {
                let i = ink_core::eval::eval_interval(&Tape::new(g, n), &ctx, [interval::I::point(0.0); 3], None);
                i.lo == 0.0 && i.hi == 0.0
            };
            if !(zero(off[0]) && zero(off[1])) {
                return None;
            }
            continue;
        }
        let th = rate(g, &ink_kernel::rigid::carried(g, off[0], off[1], &ctx, moving)?)??;
        if round.is_some_and(|r| r != th) {
            return None;
        }
        round = Some(th);
    }
    let round = round?;
    let (mut map, mut kept, mut slack) = (by.clone(), Vec::new(), 0.0);
    let mut kept_ctx = ctx.clone();
    for (w, sym, turned, _, shape) in each {
        // Its turn in the frame carried round: its own less the frame's.
        let less = round.times(-1)?;
        let own = match turned {
            None => Some(less),
            Some(c) => c.plus(less)?,
        };
        let coef = match (sym.order, own) {
            (Order::Of(n), Some(c)) => Some(c.times(n as i128)?),
            // It does not turn there: one value. A turn its symmetry keeps whole: every shift keeps it.
            (Order::Of(_), None) | (Order::Every, _) => None,
        };
        if sym.slack > 0.0 {
            let (lo, hi, _) = v.whole(shape).then(|| v.bounds(shape)).flatten()?;
            let r = [lo[0].abs().max(hi[0].abs()), lo[1].abs().max(hi[1].abs())];
            slack = interval::add_up(slack, interval::mul_up(interval::mul_up(2.0, sym.slack), r[0].hypot(r[1]).next_up()));
        }
        let kd = kept_ctx.params.len().max(v.doc.params.len()) as u32;
        kept_ctx.params.resize(kd as usize + 1, 0.0);
        kept_ctx = kept_ctx.with_band(kd, (0.0, 0.0));
        map.insert(w, v.doc.g.var(Var::Param(kd)));
        kept.push((Var::Param(kd), coef));
    }
    let tops: Vec<u32> = roots.iter().map(|&n| v.doc.g.replace(n, &map)).collect();
    // The claim is of the shapes placed alone: nothing else of space, nothing else of t.
    if v.doc.g.reachable(&tops).iter().any(|&n| n == tn || matches!(v.doc.g.nodes[n as usize], Node::Var(Var::X | Var::Y | Var::Z))) {
        return None;
    }
    let tapes: Vec<Tape> = tops.iter().map(|&n| Tape::new(&v.doc.g, n)).collect();
    let Shown::Base(b) = period_kept(&tapes.iter().collect::<Vec<_>>(), Var::Param(k), &kept_ctx, &kept)? else { return None };
    Some(Drive { node: d, window: b.window((range.lo, range.hi))?, slack, travel: range.hi - range.lo })
}

/// The values `drive` reads the claim through, in turn, each once: t (a
/// witness's value is a moment), then of the values every read of t passes
/// through, those the walk down from the claim proposes, deepest first. One
/// per run of the walk (each but the nearest begins at a sin, cos or tan of
/// the value below it), found in a pass over the claim's readers of t.
/// None: the claim reads no drive.
pub(super) fn tries(v: &ModelView, roots: &[u32], k: u32) -> Option<Vec<u32>> {
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
    let mut rates: Vec<(f64, f64)> = Vec::new();
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
        // How fast it turns with d (radians a unit, at most), and how far its points lie from its axis.
        if let Some(c) = rate {
            let reach = v
                .whole(p.shape)
                .then(|| v.bounds(p.shape))
                .flatten()
                .map_or(f64::INFINITY, |(lo, hi, _)| lo[0].abs().max(hi[0].abs()).hypot(lo[1].abs().max(hi[1].abs())).next_up());
            let c = c.value();
            rates.push((c.lo.abs().max(c.hi.abs()), reach));
        }
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
    if let Some(window) = b.window((range.lo, range.hi)) {
        return Some(Drive { node: d, window, slack, travel: range.hi - range.lo });
    }
    // The shapes repeat by their turns and the rest by its own repeats (a belt's teeth along its path),
    // the two periods written apart: a pitch radius is a rounded number, and no whole number of the
    // one is exactly a whole number of the other, though they are meant to be.
    let none: Vec<(Var, Option<ink_core::period::Coef>)> = kept.iter().map(|k| (k.0, None)).collect();
    let (Shown::Base(turns), Shown::Base(rest)) =
        (period_kept(&[], Var::Param(k), &kept_ctx, &kept)?, period_kept(&tapes.iter().collect::<Vec<_>>(), Var::Param(k), &kept_ctx, &none)?)
    else {
        return None;
    };
    nearly(turns.shift(), rest.shift(), (range.lo, range.hi)).map(|(window, off)| {
        // Each shape may sit off where the window places it by its turn over `off` of d.
        let moved = rates.iter().fold(0.0, |sum, &(rate, reach)| interval::add_up(sum, interval::mul_up(interval::mul_up(2.0, interval::mul_up(off, rate)), reach)));
        Drive { node: d, window, slack: interval::add_up(slack, moved), travel: range.hi - range.lo }
    })
}

/// Where the shapes of a claim repeat every `turns` of its drive and the rest every `rest`, each an
/// enclosure, and i of the one are j of the other but for a rounding (to a billionth; i and j to
/// NEARLY): the window that stands for every value of the drive in `range`, and how far (of the drive)
/// what is read there may be off. A value d is k whole i `turns` past one in the window's first i
/// turns, where the shapes are as at d and the rest as it is k (i turns - j rests) further on: under
/// `off`, the most that comes to over the range, so the window is that much longer at each end and
/// the shapes are placed within their turn over `off`. None: no such i and j, or the window is no
/// shorter than the range.
fn nearly(turns: interval::I, rest: interval::I, range: (f64, f64)) -> Option<((f64, f64), f64)> {
    const NEARLY: u32 = 64;
    if !(turns.lo > 0.0 && rest.lo > 0.0 && turns.hi.is_finite() && rest.hi.is_finite()) {
        return None;
    }
    let (i, j, apart) = (1..=NEARLY)
        .flat_map(|i| (1..=NEARLY).map(move |j| (i, j)))
        .map(|(i, j)| {
            let (a, b) = (interval::mul_o(interval::I::point(i as f64), turns), interval::mul_o(interval::I::point(j as f64), rest));
            (i, j, (a.hi - b.lo).abs().max((b.hi - a.lo).abs()).next_up())
        })
        .find(|&(i, _, apart)| apart <= 1e-9 * i as f64 * turns.lo)?;
    let _ = j;
    let period = interval::mul_up(i as f64, turns.hi);
    let whole = ((range.1 - range.0) / (i as f64 * turns.lo)).ceil() + 1.0;
    let off = interval::mul_up(whole, apart);
    let end = interval::add_up(interval::add_up(range.0, period), interval::mul_up(2.0, off));
    (off.is_finite() && end < range.1).then_some(((range.0, end), off))
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
pub(super) fn sweeps_far(v: &mut ModelView, roots: &[u32], span: (f64, f64)) -> Option<String> {
    let (travel, own) = sweeps(v, roots, span, FAR)?;
    Some(format!(
        "a body it reads travels {travel:.1} over the span, over {FAR} times its size ({own:.2}). Even the finest pieces of t would sweep it by more than its size, so t ranged was not tried"
    ))
}

/// A placed body `roots` read that travels over `span` more than `far` times its own size: how far,
/// and its size. Cheap: interval bounds of its placement alone.
pub(super) fn sweeps(v: &mut ModelView, roots: &[u32], span: (f64, f64), far: f64) -> Option<(f64, f64)> {
    // Each body's own placement, the outermost: one placed inside another (a member of a union asked
    // alone, placed as its body is) travels with it, and is no body of its own size.
    let within: Vec<u32> = v.doc.g.reachable(roots).into_iter().filter(|n| v.doc.placed.contains_key(n)).collect();
    let placed: Vec<(u32, Xf)> =
        within.iter().filter(|&&n| !within.iter().any(|&o| o != n && v.doc.g.reachable(&[o]).contains(&n))).filter_map(|n| v.doc.placed.get(n).map(|p| (p.shape, p.at))).collect();
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
        let Some(pieces) = swept(&v.doc.g, &v.ctx, &iso, (lo, hi), span, false) else { continue };
        let hull = pieces.iter().map(|p| p.1).reduce(ink_core::v3::box_union).expect("pieces");
        let (own, travel) = (size((lo, hi)), size(hull));
        if own > 0.0 && travel > far * own {
            return Some((travel, own));
        }
    }
    None
}
