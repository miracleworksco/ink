use super::*;

/// A number as a claim's words say it: to 4 significant figures (a
/// displacement of 1.055e-5 mm is not 0), "unbounded" where it is not finite.
pub fn sig(x: f64) -> String {
    if !x.is_finite() {
        return "unbounded".into();
    }
    if x == 0.0 {
        return "0".into();
    }
    let digits = 3 - x.abs().log10().floor() as i32;
    let r = if digits > 0 { format!("{:.*}", digits.min(17) as usize, x).parse::<f64>().unwrap_or(x) } else { (x / 10f64.powi(-digits)).round() * 10f64.powi(-digits) };
    if r.abs() < 1e-3 || r.abs() >= 1e9 { format!("{r:e}") } else { format!("{r}") }
}

/// The grid error of value `x` read on cells of h, from the same read on
/// cells of 2h (`x2`) and 4h (`x4`, if it has one): Roache's grid
/// convergence index. Where the three converge monotonically, at the
/// observed order p = log2((x4 - x2) / (x2 - x)), 1.25 |x - x2| / (2^p - 1);
/// with two grids (a read the call's time limit stopped) 3 |x - x2|, the
/// index of an assumed first order; three that oscillate (differences of
/// opposite signs) 3 times the larger difference, but none where the newer
/// is more than 3 times the older: the error stated is three times the
/// larger difference, and a swing growing by more than three would leave
/// it on the very next grid. Three
/// whose differences grow the same way by more than `noise`, what the solves
/// themselves resolve, give none, infinite: a value growing as the cells shrink may go on growing (100 N on a 0.1 mm
/// patch read 1.022e-4, 3.262e-4 and 7.736e-4 mm on cells of 4, 2 and 1 mm:
/// 3 |x - x2| held `< 0.0022 mm`, and cells of 0.125 mm read 6.336e-3). A
/// read's error or a band is no such noise: it widens the answer after the
/// grid error (`spread`), and says nothing of whether three grids approach
/// anything. Never under |x - x2| (1.25 |x - x2| where they converge: Roache's factor), the optimizer's rule (ink-design
/// `finer`). With how it was found.
pub fn grid_error(x: f64, x2: f64, x4: Option<f64>, noise: f64) -> (f64, &'static str) {
    // Never under what the solves themselves resolve: grids alike to four figures know the value no closer
    // (a conductor's layer read 49.8472 +- 0.0003 where it is 49.8475).
    let (e, how) = grid_error_read(x, x2, x4, noise);
    (if noise.is_finite() { e.max(noise) } else { e }, how)
}

fn grid_error_read(x: f64, x2: f64, x4: Option<f64>, noise: f64) -> (f64, &'static str) {
    // An infinite value (a buckling factor where nothing is compressed: no
    // buckling mode) is no number to difference: grids alike have none,
    // and one infinite where another is not leaves any value possible.
    if !(x.is_finite() && x2.is_finite() && x4.is_none_or(f64::is_finite)) {
        let alike = x == x2 && x4.is_none_or(|x4| x4 == x);
        return match (alike, x4) {
            (true, Some(_)) => (0.0, "three grids alike"),
            (true, None) => (0.0, "two grids alike"),
            (false, _) => (f64::INFINITY, "a grid's value infinite where another's is not"),
        };
    }
    let d1 = x - x2;
    match x4.map(|x4| x2 - x4) {
        Some(d2) if d1.abs() <= noise && d2.abs() <= noise => (3.0 * d1.abs().max(d2.abs()), "three grids alike within what the solves resolve"),
        // A difference no larger than what the solves resolve is no evidence of a trend, and no ratio is
        // taken against it: the grids are judged as two, by the one real difference (a motion's 0.06278 and
        // 0.06278 on cells of 4 and 2 mm, then 0.06285 on 1).
        Some(d2) if d1.abs() <= noise || d2.abs() <= noise => (3.0 * d1.abs().max(d2.abs()), "three grids, one difference within what the solves resolve"),
        Some(d2) if d1 * d2 > 0.0 && d1.abs() < d2.abs() => {
            let p = (d2 / d1).log2();
            // Roache's grid convergence index, 1.25 |x - x2| / (2^p - 1), never under 1.25 |x - x2|: the
            // standard's safety factor on the last change too, which an order over 1.17 left out (a bar's
            // hold read 0.1860121 +- 1.11e-5 where it is 0.186: 1.21e-5 off).
            (1.25 * d1.abs() * (1.0f64 / (2f64.powf(p) - 1.0)).max(1.0), "three grids converging")
        }
        Some(d2) if d1 == 0.0 && d2 == 0.0 => (0.0, "three grids alike"),
        // Differences of opposite signs: no order can be read from them, and the larger stands for the
        // error, Roache's factor for an order not observed (a notch's root read 2.774, 2.799 and 1.046 MPa:
        // 3 |x - x2|, 0.075, held a claim its 4.5 broke). Swings that grow are not told from the cells'
        // placement jitter with three readings and no measure of that jitter, but the error stated is three
        // times the larger difference: a swing growing by more than three per grid would leave it on the very
        // next grid, so none is stated (100, 99, 140: swings of 1 then 41). Two coarse grids that agree and a
        // finest that leaves them is this too, whatever the sign of their near-nothing difference: as 100,
        // 100.0001, 103 moves apart, so 100, 99.9999, 103 swings wider (a sliver's 0.1486, 0.1486, 0.1479).
        Some(d2) if d1 * d2 < 0.0 && d1.abs() > 3.0 * d2.abs() => (f64::INFINITY, "three grids swinging wider"),
        Some(d2) if d1 * d2 < 0.0 => (3.0 * d1.abs().max(d2.abs()), "three grids oscillating"),
        // Differences that grow the same way, past what the solves resolve: no order can be read from them,
        // and no error (100 N on a 0.1 mm patch read 1.022e-4, 3.262e-4 and 7.736e-4 mm on cells of 4, 2 and
        // 1 mm: 3 |x - x2| held `< 0.0022 mm`, and cells of 0.125 mm read 6.336e-3).
        Some(_) => (f64::INFINITY, "three grids moving apart"),
        None => (3.0 * d1.abs(), "two grids"),
    }
}

/// `x` widened by a grid error `e` on each side, `(x - e, x + e)`: where
/// both are infinite (grids infinite and not, `grid_error`) any value is
/// possible, not NaN.
pub fn widened(x: f64, e: f64) -> (f64, f64) {
    let (lo, hi) = (x - e, x + e);
    (if lo.is_nan() { f64::NEG_INFINITY } else { lo }, if hi.is_nan() { f64::INFINITY } else { hi })
}

/// `x` widened by a grid error `e`, and by what its read may miss where
/// that is measured (`ReadError`, not proven: a proven one decides alone).
pub(super) fn spread(x: f64, e: f64, read: Option<ReadError>) -> (f64, f64) {
    let (lo, hi) = widened(x, e);
    match read.filter(|r| !r.proven && x.is_finite()) {
        Some(r) => {
            let (a, b) = r.span(x);
            (lo - (x - a), hi + (b - x))
        }
        None => (lo, hi),
    }
}

/// The tolerance a ladder's grids are solved to, relative: values that
/// differ by no more are alike (`grid_error`), as are those within what the
/// solve itself does not resolve, where its package states it
/// (`Measured::resolves`: a study's stiffness kept in void; `noise`).
pub(super) const LADDER_TOL: f64 = 1e-8;

/// The grid error of the value on the finest of `grids` (two or three,
/// finest first, each twice as coarse as the last: `grid_error`), the one
/// way a claim's grids are compared. Where a grid splits cells finer than
/// its own (a study's where a feature is thinner than two of them), its finest cells are
/// what its answer rests on there (`Measured::finest`), and grids whose
/// finest cells are not each at least twice the last's give none: grids of
/// 1, 0.5 and 0.25 mm all split to 0.25 at a ligament read it alike, and that
/// agreement is no convergence (such a ligament held 6.1e-5 mm, grid error
/// 1.3e-7, against 1.1e-4 on 0.05 mm cells). Three grids converging faster
/// than their elements can (an observed order above `order`, the metric's:
/// Package::metric_order) are not
/// yet where their difference measures the error: a coarsest grid all but
/// blind to a ligament reads it rigid, and the finer two's difference, a
/// quarter of the error left, held a claim 23% short. Nor does their finest
/// pair's agreement measure it: an order above what the elements can give
/// is luck, and a bar held at x < 1 read 0.18616 +- 0.00046 (3 |x - x2|) on
/// 0.3 mm cells where 0.075 read 0.18904. Their error is three times the
/// whole movement they show, finest against coarsest, 3 |x - x4|. Where the
/// coarsest could not see the part, that is wide, and a grid finer shows
/// whether the answer was luck (a flexure's travel 0.4722, 0.4726, then
/// 0.4736 on cells of 0.125 mm: moving apart). Reading that finer grid
/// whenever three converge so (tried, 2026-11-24) answered no more: the
/// flexure's eight claims all moved apart on it, its first such solve took
/// 167 s and the process 20 GB, and a pressure vessel's would have needed
/// 31 GB; only where the wide band leaves a claim undecided would it be
/// worth it, and then with the memory of the solves held before it counted.
/// Each grid must have a value. Differences within what the solves resolve
/// (`noise`) are no evidence of a trend.
pub fn grids_error(grids: &[&Measured], order: f64) -> Result<(f64, &'static str), String> {
    let finest: Vec<f64> = grids.iter().map(|g| g.finest.or(g.res).expect("a grid's measurement")).collect();
    if !finest.windows(2).all(|w| w[1] >= 2.0 * w[0] * (1.0 - 1e-9)) {
        return Err(format!(
            "its grids read the part on cells of {} mm at their finest (split where a feature is thinner than two of their cells), not each at least twice the last: \
             they do not show how the answer moves with the cells. A grid stops splitting at cells half the feature's thickness, so the feature is at least {} mm \
             thick, and at cells an eighth of that none of the three grids splits: refine (res r) to {} mm or finer (halving once may only repeat this), or thicken the feature",
            finest.iter().map(|v| format!("{v}")).collect::<Vec<_>>().join(", "),
            2.0 * finest.iter().copied().fold(f64::INFINITY, f64::min),
            finest.iter().copied().fold(f64::INFINITY, f64::min) / 4.0
        ));
    }
    let x: Vec<f64> = grids.iter().map(|g| g.value.clone().expect("a grid's value")).collect();
    // Differences within what the solves resolve are no evidence of an order, faster or slower: they are
    // asked of the noise first (`grid_error`). A reading 3.6e-15 then 1.6e-14 apart under a floor of 3.2e-4
    // was "converging faster than its elements can".
    let floor = noise(grids);
    if let [x, x2, x4] = x[..]
        && (x - x2) * (x2 - x4) > 0.0
        && (x - x2).abs() > floor
        && (x2 - x4).abs() > floor
        && ((x2 - x4) / (x - x2)).log2() > order
    {
        return Ok((3.0 * (x - x4).abs(), "three grids converging faster than their elements can"));
    }

    let (e, how) = grid_error(x[0], x[1], x.get(2).copied(), floor);
    // Close but not in order: the refusal guards a value that grows without bound as the cells shrink (100 N on a
    // 0.1 mm patch read 1.022e-4, 3.262e-4, 7.736e-4, then 6.336e-3 on cells of 0.125 mm). Where no grid's solve
    // has found anything that could (`Measured::unbounded`), three times the larger change covers it: on every
    // recorded ladder with a truth, 2026-11-25 (a heat hold's 129.156, 129.203, 129.273 for 129.28; columns'
    // buckling factors on cells off their faces).
    if let [x, x2, x4] = x[..]
        && matches!(how, "three grids moving apart" | "three grids swinging wider")
        && grids.iter().all(|g| g.unbounded.is_none())
        && !shrinks_with_the_cells(grids)
    {
        let larger = (x - x2).abs().max((x2 - x4).abs());
        return Ok(((3.0 * larger).max(floor), "three grids close, not in order"));
    }
    Ok((e, how))
}

/// Whether the mode a value is of shrinks with the cells (`Measured::extent`):
/// on the finest grid it reaches under three quarters of what it reaches on
/// the next. A mode spread over the part keeps its reach (a ratio near 1);
/// one held to a few cells halves it with them.
fn shrinks_with_the_cells(grids: &[&Measured]) -> bool {
    matches!(grids, [a, b, ..] if a.extent.zip(b.extent).is_some_and(|(a, b)| a < 0.75 * b))
}

/// What grids' values may differ by and be alike: the ladder's solve
/// tolerance, of the largest value (an iterative solve's error is a share of
/// the whole solution, offset and all: a block carried 1000 mm read its
/// deflection 1e-6 apart from the same block held still), or what any grid's
/// solve does not resolve (`Measured::resolves`), whichever is more.
fn noise(grids: &[&Measured]) -> f64 {
    let tol = LADDER_TOL * grids.iter().filter_map(|g| g.value.as_ref().ok()).fold(0.0_f64, |m, v| m.max(v.abs()));
    grids.iter().fold(tol, |m, g| m.max(g.resolves))
}

/// The grids a value on cells of h is checked on: three, each twice as
/// coarse as the last, whose coarsest is the coarsest of 4h, 2h and h that
/// reads the part (solves it, and sees what it holds). A part thinner than
/// the coarser cells hold (a printed wall, a sheet) is so checked on cells
/// finer than its own, and answered on the finest. A claim, holding or
/// failing, is checked so (`metric_claim`), and an optimize's
/// search measures each design so, holding or not, so that the design it
/// returns keeps the margin its claim needs.
pub struct Ladder {
    /// The finest grid's measurement: the answer.
    pub finest: Measured,
    /// The two coarser grids' measurements, finer first.
    pub coarser: [Measured; 2],
    /// The value on each grid, finest first, and its cells.
    pub values: [f64; 3],
    pub cells: [f64; 3],
    /// Each grid's cells in words: "2 mm", or "2 mm split to 1 mm" where
    /// it split them (`Measured::finest`).
    pub named: [String; 3],
    /// The finest value's grid error, and how it was found (`grids_error`).
    error: (f64, &'static str),
}

impl Ladder {
    /// Of three grids, finest first, each twice as coarse as the last, if
    /// they give a grid error (`grids_error`).
    fn of(grids: [&Measured; 3], order: f64) -> Result<Ladder, String> {
        let error = grids_error(&grids, order)?;
        let at = |i: usize| -> f64 { grids[i].res.expect("a grid's measurement") };
        let values = grids.map(|g| g.value.clone().expect("a grid's value"));
        let named = grids.map(|g| match (g.res, g.finest) {
            (Some(h), Some(f)) if f < h => format!("{} mm split to {} mm", sig(h), sig(f)),
            (h, _) => format!("{} mm", sig(h.expect("a grid's measurement"))),
        });
        Ok(Ladder { finest: grids[0].clone(), coarser: [grids[1].clone(), grids[2].clone()], values, cells: [at(0), at(1), at(2)], named, error })
    }

    /// The finest value's grid error, and how it was found.
    pub fn error(&self) -> (f64, &'static str) {
        self.error
    }

    /// Why a grid's value may grow without bound (`Measured::unbounded`): the first grid, finest first, to say;
    /// else a mode that shrinks with the cells.
    pub fn unbounded(&self) -> Option<String> {
        let grids = [&self.finest, &self.coarser[0], &self.coarser[1]];
        grids.iter().find_map(|g| g.unbounded.clone()).or_else(|| {
            shrinks_with_the_cells(&grids).then(|| {
                format!(
                    "its mode shrinks with the cells, reaching {} mm on cells of {} and {} mm on cells of {}: it lives in a few of them, where nothing bounds it",
                    sig(self.finest.extent.unwrap_or(0.0)),
                    self.named[0],
                    sig(self.coarser[0].extent.unwrap_or(0.0)),
                    self.named[1]
                )
            })
        })
    }

    /// The grid error of what the finest grid leaves out under `key`
    /// (`Measured::excluded`), as the value's is found: from its value
    /// counted on each grid (a grid that leaves nothing out under it counts
    /// its own value); 0 where it is unbounded. With the side of the value
    /// counting it falls on (+1 above, -1 below, 0 on none): the first grid
    /// whose count differs from its value says.
    pub fn excluded_error(&self, key: &str) -> (f64, f64) {
        let grids = [&self.finest, &self.coarser[0], &self.coarser[1]];
        let counted: Vec<f64> = grids.iter().zip(self.values).map(|(m, x)| m.excluded.iter().find(|e| e.key == key).map_or(x, |e| e.value)).collect();
        let side = counted.iter().zip(self.values).map(|(c, x)| c - x).find(|d| *d != 0.0).map_or(0.0, f64::signum);
        // What each grid's solve does not resolve, as a share of its value, of what it counts.
        let share = grids.iter().zip(self.values).fold(LADDER_TOL, |m, (g, x)| if x != 0.0 { m.max(g.resolves / x.abs()) } else { m });
        let noise = share * counted.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let e = if counted.iter().all(|v| v.is_finite()) { grid_error(counted[0], counted[1], Some(counted[2]), noise).0 } else { 0.0 };
        (e, side)
    }
}

/// The ladder of the measurement `base` (`Ladder`), or the cells it was
/// asked on and why no ladder reads the part. Each other grid splits
/// cells where the claim's own would, as deep (`RunOpts::split_as`), so
/// that its cells differ from those everywhere as its size does, and the
/// grids can be compared (`grids_error`).
pub fn ladder(v: &mut ModelView, m: &FormMetric, accuracy: Option<f64>, base: &Measured) -> Result<Ladder, (f64, String)> {
    let h = base.res.expect("a grid's measurement");
    let order = v.metric_order(m);
    // Each grid splits cells where h's grid would, as deep (`RunOpts::split_as`):
    // its cells there differ from h's as its own size does.
    // A grid that cannot see the part says whether a finer one would
    // (`Measured::unseen_finer`): one blind to a peak's own curve buys no
    // finer grids, and leaves the claim unknown.
    let mut read = |k: f64| -> Result<Measured, (String, bool)> {
        let r = v.measure(m, RunOpts { res: Some(k * h), tol: LADDER_TOL, accuracy, split_as: Some(h), both_ways: false }).map_err(|e| (e.message, true))?;
        match (&r.unseen, &r.value) {
            (Some(why), _) => Err((format!("cells of {} mm cannot see the part ({why})", k * h), r.unseen_finer)),
            (None, Err(why)) => Err((format!("cells of {} mm read nothing ({why})", k * h), true)),
            (None, Ok(_)) => Ok(r),
        }
    };
    let of = |grids: [&Measured; 3]| Ladder::of(grids, order).map_err(|why| (h, why));
    let up = read(2.0).and_then(|x2| read(4.0).map(|x4| (x2, x4)));
    let up = match up {
        Ok((x2, x4)) => return of([base, &x2, &x4]),
        Err((why, false)) => return Err((h, why)),
        Err((why, true)) => why,
    };
    // 2h reads the part (4h does not): h/2, h, 2h; else h/4, h/2, h.
    let fail = |(why, _): (String, bool)| (h, format!("{up}; {why}"));
    match read(2.0) {
        Ok(x2) => {
            let f = read(0.5).map_err(fail)?;
            of([&f, base, &x2])
        }
        Err(_) => {
            let f2 = read(0.5).map_err(fail)?;
            let f4 = read(0.25).map_err(fail)?;
            of([&f4, &f2, base])
        }
    }
}
