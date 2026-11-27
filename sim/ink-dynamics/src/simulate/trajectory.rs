use super::{DEVIATION, EQUATIONS, INTEGRATOR, PER_PERIOD, RESOLVE, SPRINGS_METHOD, TOLERANCE, UNRESOLVED};

/// A simulated motion: its knots' moments, and at each its components'
/// positions, speeds and accelerations (the simulated coordinates in their
/// units, then the deviations of the released joints carrying inertia,
/// radians a turn, mm a slide) and each spring line's tension (N); where it
/// stopped short of the span's end, when and why.
pub struct Trajectory {
    pub t: Vec<f64>,
    pub x: Vec<Vec<[f64; 3]>>,
    /// Position increment of each accepted step, before adding it to the
    /// starting position. Taking the difference of rounded positions can
    /// erase the displacement that fixes the higher derivatives.
    pub(super) displacement: Vec<Vec<f64>>,
    /// Each step, whether it was taken across a switch in what drives the
    /// motion (at the finest): an acceleration jumps inside it, from its
    /// start's to its end's, and is read as those (`at`), not off the
    /// quintic, which holds the step's mean speed by bulging past both (a
    /// force pulse's 1.0e6 read 1.113e6 there).
    pub(super) across: Vec<bool>,
    pub tensions: Vec<Vec<f64>>,
    /// At each knot, the largest deviation so far at the parts the springs
    /// move (mm: Simulation::scales), and each stepped component's radius
    /// of gyration (the deviations' tolerance and band read them).
    pub scale: Vec<f64>,
    pub gyration: Vec<f64>,
    /// Whether the springs are stepped (a joint released).
    pub springs: bool,
    pub stopped: Option<(f64, String)>,
    pub rejected: usize,
    pub evaluations: usize,
}

/// Interior sign changes of a polynomial on [0, 1]. Its derivative's
/// sign changes partition it into monotone pieces. A zero at one of those
/// cuts is harmless to include even when the polynomial only touches zero.
fn crossings(c: &[f64]) -> Vec<f64> {
    let Some(degree) = c.iter().rposition(|&v| v != 0.0) else { return Vec::new() };
    if degree == 0 {
        return Vec::new();
    }
    if degree == 1 {
        let root = -c[0] / c[1];
        return (0.0 < root && root < 1.0).then_some(root).into_iter().collect();
    }
    let derivative: Vec<_> = (1..=degree).map(|k| k as f64 * c[k]).collect();
    let cuts = crossings(&derivative);
    let value = |s: f64| c.iter().rev().fold(0.0, |y, a| y * s + a);
    let mut roots = Vec::new();
    let mut lo = 0.0;
    for hi in cuts.into_iter().chain([1.0]) {
        let (a, b) = (value(lo), value(hi));
        if lo > 0.0 && a == 0.0 {
            roots.push(lo);
        }
        if (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) {
            roots.push(crate::steps::turn_back(lo, hi, a < 0.0, value));
        }
        lo = hi;
    }
    roots
}

impl Trajectory {
    /// Piece j's quintic Hermite for component i, through its ends'
    /// positions, speeds and accelerations: its coefficients in s, the
    /// share of the piece (0 to 1), and the piece's length.
    fn piece(&self, j: usize, i: usize) -> ([f64; 6], f64) {
        const H: [[f64; 6]; 6] = [
            [1.0, 0.0, 0.0, -10.0, 15.0, -6.0],
            [0.0, 1.0, 0.0, -6.0, 8.0, -3.0],
            [0.0, 0.0, 0.5, -1.5, 1.5, -0.5],
            [0.0, 0.0, 0.0, 0.5, -1.0, 0.5],
            [0.0, 0.0, 0.0, -4.0, 7.0, -3.0],
            [0.0, 0.0, 0.0, 10.0, -15.0, 6.0],
        ];
        let h = self.t[j + 1] - self.t[j];
        let ([p0, v0, a0], [_, v1, a1]) = (self.x[j][i], self.x[j + 1][i]);
        let c = [0.0, h * v0, h * h * a0, h * h * a1, h * v1, self.displacement[j][i]];
        let mut m = std::array::from_fn(|k| (0..6).map(|n| c[n] * H[n][k]).sum());
        m[0] = p0;
        (m, h)
    }

    /// The coefficients in s of a piece's derivative of order d along t.
    fn derived(m: &[f64; 6], h: f64, d: usize) -> Vec<f64> {
        (d..6).map(|k| ((k + 1 - d)..=k).map(|f| f as f64).product::<f64>() * m[k] / h.powi(d as i32)).collect()
    }

    fn no_piece(&self) -> String {
        self.stopped.as_ref().map_or_else(|| "the simulation has no step from which to read jerk or snap".into(), |(_, why)| format!("the simulation has no first step: {why}"))
    }

    /// Component i's position, speed, acceleration, jerk and snap at t:
    /// between its knots the quintic Hermite through theirs (at a knot, the
    /// piece it starts). Err: t outside what was simulated.
    pub fn at(&self, i: usize, t: f64) -> Result<[f64; 5], String> {
        let (a, b) = (self.t[0], self.t[self.t.len() - 1]);
        if !(t >= a && t <= b) {
            return Err(match &self.stopped {
                Some((at, why)) if t > b => format!("the simulation stopped at t = {at}: {why}"),
                _ => format!("the simulated motion is over the file's span from {a} to {b}: t = {t} is outside it"),
            });
        }
        if self.t.len() == 1 {
            return Err(self.no_piece());
        }
        let j = self.t.partition_point(|&x| x <= t).clamp(1, self.t.len() - 1) - 1;
        let (m, h) = self.piece(j, i);
        let s = (t - self.t[j]) / h;
        let horner = |c: &[f64]| c.iter().rev().fold(0.0, |acc, x| acc * s + x);
        let mut out: [f64; 5] = std::array::from_fn(|d| horner(&Self::derived(&m, h, d)));
        if self.across[j] {
            out[2] = self.x[if s < 0.5 { j } else { j + 1 }][i][2];
        }
        Ok(out)
    }

    /// Position extrema used by the spring and ringing reports.
    pub fn moments(&self, i: usize, span: (f64, f64)) -> Vec<f64> {
        self.extremes(i, span, 1)
    }

    /// Position (rate_order = 1) or speed (2) extrema of the quintic
    /// interpolant. Search each monotone piece of the rate polynomial:
    /// it can change sign twice between two like-signed knots.
    pub(crate) fn extremes(&self, i: usize, span: (f64, f64), rate_order: usize) -> Vec<f64> {
        let mut times = Vec::new();
        for (j, w) in self.t.windows(2).enumerate().filter(|(_, w)| w[0] < span.1 && w[1] > span.0) {
            let (m, h) = self.piece(j, i);
            let rate = Self::derived(&m, h, rate_order);
            if rate.iter().all(|&a| a == 0.0) {
                continue;
            }
            if self.x[j][i][rate_order] == 0.0 {
                times.push(w[0]);
            }
            times.extend(crossings(&rate).into_iter().map(|s| w[0] + h * s));
            if self.x[j + 1][i][rate_order] == 0.0 {
                times.push(w[1]);
            }
        }
        times.retain(|&t| span.0 < t && t < span.1);
        times.dedup();
        times
    }

    /// A scalar derivative limit needs the minimum and maximum of the
    /// interpolant on the requested span. Read its stationary points and
    /// both ends of each piece: jerk and snap can jump at a knot. Retain
    /// just the winning times so geometric claims do not visit every knot.
    pub(crate) fn derivative_peaks(&self, i: usize, span: (f64, f64)) -> Vec<f64> {
        let mut peaks = [[(f64::INFINITY, span.0), (f64::NEG_INFINITY, span.0)]; 3];
        for (j, w) in self.t.windows(2).enumerate().filter(|(_, w)| w[0] <= span.1 && w[1] >= span.0) {
            let (m, h) = self.piece(j, i);
            let lo = w[0].max(span.0);
            let hi = (if j + 2 == self.t.len() { w[1] } else { w[1].next_down() }).min(span.1);
            if lo > hi {
                continue;
            }
            for (d, [min, max]) in (2..=4).zip(&mut peaks) {
                let coeff = Self::derived(&m, h, d);
                // (Across a switch the acceleration is its ends', as `at`
                // reads it.)
                let jumps = d == 2 && self.across[j];
                let roots = if jumps { Vec::new() } else { crossings(&Self::derived(&m, h, d + 1)) };
                for t in roots.into_iter().map(|r| w[0] + h * r).chain([lo, hi]).filter(|&t| lo <= t && t <= hi) {
                    let u = (t - w[0]) / h;
                    let x = if jumps { self.x[if u < 0.5 { j } else { j + 1 }][i][2] } else { coeff.iter().rev().fold(0.0, |y, a| y * u + a) };
                    if x < min.0 {
                        *min = (x, t);
                    }
                    if x > max.0 {
                        *max = (x, t);
                    }
                }
            }
        }
        peaks.into_iter().flatten().map(|(_, t)| t).collect()
    }

    /// What answers reading it say of it. It is read before a claim steps
    /// the simulation, so it says how the motion is found, not how far it
    /// has got (a claim reading past where it stops is refused, saying why).
    pub fn note(&self) -> String {
        if self.springs {
            format!(
                "the motion of each (simulate ...) coordinate and each spring's joint is an estimate, {EQUATIONS}; {SPRINGS_METHOD}; its relative and \
                 absolute tolerance {TOLERANCE:e} a step (a deviation's absolute {RESOLVE:e} of the largest deviation the springs have made so far \
                 at the parts they move, at least {DEVIATION:e} mm; a speed's by what it moves over the step), read between its steps by the \
                 quintic Hermite through each end's position, speed and acceleration; a following error is known either way to within its band, \
                 {UNRESOLVED} of a deviation's tolerance and {PER_PERIOD} more each period of the springs' fastest mode stepped (Radau IIA damps an \
                 undamped ringing near the tolerance by about 0.2 of them a period), a claim within it of its bound undecided; the band is \
                 measured against closed forms over up to 1000 periods (a 335 Hz ringing over 3 s, a 50 Hz one over 15 s), no bound: a run of \
                 more may pass it; {MASSES}"
            )
        } else {
            format!(
                "the motion of each (simulate ...) coordinate is an estimate, {EQUATIONS}; {INTEGRATOR}, its relative and absolute tolerance {TOLERANCE:e} a step (a \
                 speed's its position's times its own frequency, or, once stiff, by what it moves over the step), read between its steps by the \
                 quintic Hermite through each end's position, speed and acceleration; {MASSES}"
            )
        }
    }
}

/// What a simulated motion reads of the bodies' masses: what they are measured to moves it, by no bound found.
const MASSES: &str = "the bodies' masses, centres and inertias read as measured, exactly: how far what their measures are known to \
    moves the motion is not bounded, so a claim near its limit may be overturned by it";
