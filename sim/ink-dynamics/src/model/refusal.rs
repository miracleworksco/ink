use super::drives::UNCERTAIN_DRIVES;
use super::{Model, Moment, NO_MASS, Num, ROUNDING};
use ink_core::eval::Ctx;
use ink_lang::Refused;

impl Model {
    /// The moving bodies with no mass (their labels) that actuator a's effort would rest on: those moving
    /// along a coordinate coupled with one it rests on. Empty: its effort is as if none were left out.
    pub fn blocked(&self, a: usize) -> Vec<&str> {
        self.massless_in(&self.actuators[a].coupled)
    }

    /// The same of reaction r: those moving along a coordinate coupled with one moving what it carries (a
    /// coupling's pull on it rests on those coordinates' efforts).
    pub fn blocks_reaction(&self, r: usize) -> Vec<&str> {
        self.massless_in(&self.reactions[r].coords)
    }

    /// The same of coupling c's tensions (they rest on the efforts along the coordinates it relates), and of
    /// the effort along coordinate k.
    pub fn blocks_coupling(&self, c: usize) -> Vec<&str> {
        // The coordinates its ends move along: their own placements', and the joints' they ride.
        let ends = self.couplings[c].ends.iter();
        let coords: Vec<usize> = ends.flat_map(|e| e.dq.iter().map(|(k, _)| *k).chain(e.input.iter().flat_map(|i| i.along.iter().map(|(k, _)| *k)))).collect();
        self.massless_in(&coords)
    }

    pub fn blocks_coord(&self, k: usize) -> Vec<&str> {
        self.massless_in(&[k])
    }

    fn massless_in(&self, coords: &[usize]) -> Vec<&str> {
        let sets: Vec<usize> = coords.iter().map(|&k| self.sets[k]).collect();
        self.massless.iter().filter(|(_, _, along)| along.iter().any(|k| sets.contains(&self.sets[*k]))).map(|(label, ..)| label.as_str()).collect()
    }

    /// A passive joint may be fixed in place. A prescribed motion or a
    /// coordinate held by a declared drive must balance even at rest.
    pub(super) fn needs_drive(&self, k: usize) -> bool {
        !self.simulated.iter().any(|s| s.0 == k) && (self.coords[k].moves || self.actuators.iter().any(|a| !matches!(a.along[k], Num::C(0.0))))
    }

    /// Why actuator a's effort is not this motion's at any moment, if it
    /// leaves out deforming bodies moving with its coordinates (Model::
    /// refusal, its part that reads no moment).
    pub fn leaves_out(&self, a: usize) -> Option<String> {
        self.left_out(a).map(String::from)
    }

    /// `leaves_out`, with its kind: outside what the model covers.
    pub fn left_out(&self, a: usize) -> Option<Refused> {
        let out = &self.actuators[a].leaves_out;
        (!out.is_empty()).then(|| Refused::outside(format!("its effort leaves out bodies whose shape moves with what moves it, whose inertia is not counted: {}", out.join(", "))))
    }

    /// Why actuator a's effort at moment m is not this motion's, if it is
    /// not: it leaves out deforming bodies moving with its coordinates, or
    /// a coordinate it rests on needs effort no actuator supplies (the
    /// motion as written cannot happen).
    pub fn refusal(&self, a: usize, m: &Moment) -> Option<String> {
        self.refusal_words(a, m, false)
    }

    /// `refusal`, but that a coupling passing power with no efficiency given is taken as lossless (`unknown_losses`
    /// names them): for one who says that it is assumed.
    pub fn refusal_lossless(&self, a: usize, m: &Moment) -> Option<String> {
        self.refusal_words(a, m, true)
    }

    /// Moving bodies with no mass that actuator a's effort rests on (`blocked`),
    /// refused as something the design lacks (`Open::Lacks`).
    pub fn missing_mass(&self, a: usize) -> Option<Refused> {
        let without = self.blocked(a);
        (!without.is_empty()).then(|| Refused::lacks(format!("these moving bodies have no mass: {}: {NO_MASS}", self.massless_said(&without))))
    }

    /// `refusal`, with its kind where it has one (`missing_mass`, `left_out`).
    pub fn refused(&self, a: usize, m: &Moment) -> Option<Refused> {
        if let Some(why) = &m.unsolved {
            return Some(Refused::from(why.clone()));
        }
        if let Some(why) = self.missing_mass(a) {
            return Some(why);
        }
        self.refusal_words(a, m, false).map(|why| match self.left_out(a) {
            Some(left) if left.why == why => left,
            _ => Refused::from(why),
        })
    }

    fn refusal_words(&self, a: usize, m: &Moment, lossless: bool) -> Option<String> {
        if let Some(why) = &m.unsolved {
            return Some(why.clone());
        }
        let act = &self.actuators[a];
        // Built without the moving bodies that have no mass: an effort one of them would change is not found.
        if let Some(why) = self.missing_mass(a) {
            return Some(why.why);
        }
        if m.pm[a].iter().any(|p| !p.is_finite()) {
            return Some(UNCERTAIN_DRIVES.into());
        }
        if let Some(why) = self.coupling_refusal(&act.coupled, m, lossless) {
            return Some(why);
        }
        // Its couplings' loss rests on tensions that leave the unseparated
        // bodies out: within that of a limit, undecided.
        let (b, x) = (m.loss_bound[a], m.margin(a));
        if b > 0.0 && -b <= x && x <= b {
            return Some(format!(
                "its couplings' loss rests on tensions that leave out what {} exert, moving with its coordinates through nothing a joint or \
                 coupling declares: it may be off by up to {b:.6}, and its margin is {x:.6}: undecided (place them on the joint they ride, or \
                 relate them)",
                self.unseparated.join(", ")
            ));
        }
        if let Some(why) = self.leaves_out(a) {
            return Some(why);
        }
        // Another actuator's limits crossing leave this effort as it is, but where the two drive one coordinate
        // between them (their shares rest on each other's limits).
        let own = |b: usize| (0..act.along.len()).filter(|&k| !matches!(self.actuators[b].along[k], Num::C(x) if x == 0.0)).collect::<Vec<usize>>();
        if let Some(why) = (0..self.actuators.len()).filter(|&b| b == a || own(b).iter().any(|k| own(a).contains(k))).find_map(|b| self.crossed_at(b, m)) {
            return Some(why);
        }
        // coupled is the transitive closure of all coordinates joined by
        // drives and bodies, including held axes outside this motor's row.
        // Read over a span, an unpowered effort is within bound_u of the
        // one read, at any angle of the spins read at rest: past that it
        // is one, under it it may be.
        act.coupled.iter().find_map(|&k| {
            let (r, b) = (m.unpowered[k].unwrap_or(0.0), m.bound_u[k]);
            if m.unpowered[k].is_some() && (!r.is_finite() || r.abs() > b) {
                return Some(self.unsupplied(m, k, r));
            }
            (self.needs_drive(k) && b > m.noise_f[k]).then(|| {
                format!(
                    "the motion as written may need up to {:.6} N mm per unit of {} that no actuator supplies, at some angle of the bodies read with their spin at rest (motion.basis): undecided",
                    r.abs() + b,
                    self.coords[k].name
                )
            })
        })
    }

    /// Why the couplings passing power along `coords` at moment m leave
    /// what moves them unknown, if they do: one whose efficiency is not
    /// given, or a lossy one whose tensions are not found.
    fn coupling_refusal(&self, coords: &[usize], m: &Moment, lossless: bool) -> Option<String> {
        let ctx = Ctx::new(m.t, self.params.clone());
        for (ci, c) in self.couplings.iter().enumerate() {
            let moves = c.moves(coords);
            let lossy = c.efficiency.as_ref().is_none_or(|(e, b)| e.at(&ctx) != 1.0 || b.at(&ctx) != 1.0);
            if moves && m.lossless.contains(&ci) {
                if lossless {
                    continue;
                }
                if c.efficiency.is_some() {
                    return Some(format!("{}, so its loss is unknown", self.loss_on_a_carrier(ci)));
                }
                return Some(format!(
                    "{} passes power along its coordinates and gives no efficiency, so its loss is unknown (belt-drive, chain-drive, belt-span, \
                     the gear meshes and lead-screw take a published one last)",
                    c.name
                ));
            }
            if moves
                && lossy
                && let Some(why) = &m.coupling
            {
                return Some(format!("{} passes power along its coordinates, and its loss rests on its tensions, which are not found: {why}", c.name));
            }
        }
        None
    }

    /// Why a simulation cannot go on from state m (simulate.rs), if it
    /// cannot: what inverse dynamics refuses of the efforts along the
    /// simulated coordinates. A coupling passing power along them whose
    /// loss is unknown, or rests on tensions not found or that leave out
    /// what unseparated bodies exert; a body moving with them whose inertia
    /// is left out (its shape moves).
    pub(crate) fn forward_refusal(&self, m: &Moment) -> Option<String> {
        let sims: Vec<usize> = self.simulated.iter().map(|s| s.0).collect();
        if let Some(why) = self.coupling_refusal(&sims, m, false) {
            return Some(why);
        }
        let ctx = Ctx::new(m.t, self.params.clone());
        if !self.unseparated.is_empty() && self.couplings.iter().any(|c| c.moves(&sims) && c.efficiency.as_ref().is_some_and(|(e, b)| e.at(&ctx) != 1.0 || b.at(&ctx) != 1.0)) {
            return Some(format!(
                "a coupling's loss along the simulated coordinates rests on tensions that leave out what {} exert (moving with the joints whose balances \
                 give them, carried by no actuated or balancing joint): its forces are not found",
                self.unseparated.join(", ")
            ));
        }
        let left: Vec<&str> = self.leaves.iter().filter(|(_, reads)| reads.iter().any(|k| sims.contains(k))).map(|(l, _)| l.as_str()).collect();
        (!left.is_empty()).then(|| format!("bodies whose shape moves move with the simulated coordinates, their inertia not counted: {}", left.join(", ")))
    }

    /// Why reaction r at moment m is not this motion's, if it is not: a
    /// body it carries has no mass, or its shape moves (its inertia left
    /// out), or a coordinate moving what it carries needs effort no
    /// actuator supplies.
    pub fn reaction_refusal(&self, r: usize, m: &Moment) -> Option<String> {
        self.reaction_refused(r, m).map(String::from)
    }

    /// `reaction_refusal`, with its kind where it has one: a body with no
    /// mass or no bound is the design's to give; a body left out, or loads
    /// the model does not separate, are outside what it covers.
    pub fn reaction_refused(&self, r: usize, m: &Moment) -> Option<Refused> {
        self.reaction_refused_as(r, m, false)
    }

    /// `reaction_refused`, but that a coupling passing power with no
    /// efficiency given is taken as lossless (`unknown_losses` names them):
    /// for one who says that it is assumed (a requirement of ratings).
    pub fn reaction_refused_lossless(&self, r: usize, m: &Moment) -> Option<Refused> {
        self.reaction_refused_as(r, m, true)
    }

    fn reaction_refused_as(&self, r: usize, m: &Moment, lossless: bool) -> Option<Refused> {
        // A moment not solved (a simulation stopped short of it) has no
        // reaction to read: its numbers are NaN, which a claim would read.
        if let Some(why) = &m.unsolved {
            return Some(Refused::from(why.clone()));
        }
        if let Some(why) = self.static_refusal(r) {
            return Some(why);
        }
        let (re, why) = (&self.reactions[r], self.reaction_words(r, m, lossless)?);
        // The sum of what bodies moving with it exert is not its joint's load: the model does not cover theirs.
        Some(if re.moved_by_t.is_empty() && !re.moves_with.is_empty() { Refused::outside(why) } else { Refused::from(why) })
    }

    fn reaction_words(&self, r: usize, m: &Moment, lossless: bool) -> Option<String> {
        // One summing what bodies moving with it exert on what it carries,
        // where no joint the model declares meets them (a linkage's pin), is
        // not its joint's load, and the model does not determine theirs.
        let re = &self.reactions[r];
        // One carrying a body t moves through no coordinate: nothing the
        // model declares holds that body, so neither is what holds it.
        if !re.moved_by_t.is_empty() {
            return Some(format!(
                "it carries bodies placed by a formula of t, through no coordinate or joint, which nothing the model declares moves or holds \
                 (driven_by_nothing): {}. What holds them, and so the load where {} meets the rest, is not determined: place them through a \
                 coordinate (an argument of their assembly), or on a joint that moves them",
                re.moved_by_t.join(", "),
                re.path
            ));
        }
        if !re.moves_with.is_empty() {
            return Some(format!(
                "it is not the load where {} meets the rest: {} moves with it and is not on it, and nothing says what joins them, so the push between them \
                 is not known (moves_with). Put what it moves on its frame, (instance NAME CHILD (frame-of ITS-PART FRAME)): a drum on a motor's shaft is \
                 (frame-of motor shaft), not a joint of its own. Or relate the two by a coupling, one line: (relate A B KIND), KIND std's (cord R) for \
                 what hangs from a drum, (belt-span ...), (lead-screw LEAD), a gear mesh",
                re.path,
                re.moves_with.join(", ")
            ));
        }
        // One a coupling reaches rests on the balances giving its tensions,
        // those of its couplings' trains.
        let trains: Vec<usize> = (0..self.couplings.len()).filter(|&c| self.reaches(c, r)).map(|c| self.trains[c]).collect();
        if !trains.is_empty()
            && let Some(why) = self.tension_refusal_in(m, &trains)
        {
            return Some(format!("its couplings' tensions are not found: {why}"));
        }
        // A coupling of those trains passing power with no efficiency given:
        // its loss's friction, in the reactions and the tensions, is unknown.
        if !lossless && let Some(why) = self.loss_unknown(m, &trains) {
            return Some(why);
        }
        self.unpowered(&self.reactions[r].coords, m)
    }

    /// What a claim on an effort or a reaction says of bodies t moves through
    /// no coordinate, where there are any: no actuator is charged for them,
    /// and where they touch what it reads, their push is not in it.
    /// What a reaction's claim says of the spinning bodies it carries whose centre is not exact.
    pub fn spin_centre_note(&self, r: usize) -> Option<String> {
        let each: Vec<String> =
            self.spin_centres.iter().filter(|m| m.within > 0.0 && m.of.contains(&r)).map(|m| format!("{} (within {} mm)", m.body, ink_lang::claims::sig(m.within))).collect();
        (!each.is_empty()).then(|| {
            format!(
                "spinning bodies whose centre is known only as far as its measure goes: {}. What so far off it would add at their spin is in this \
                 reaction's band. To narrow it: slow the spin, or centre the body on its spin (a shape centred on its axis, or repeating about \
                 it, has its centre there exactly)",
                each.join(", ")
            )
        })
    }

    pub fn moved_by_t_note(&self) -> Option<String> {
        (!self.held_by_nothing.is_empty()).then(|| {
            format!(
                "bodies placed by a formula of t, through no coordinate or joint, move by nothing the model declares (driven_by_nothing): {}. No \
                 actuator is charged for their motion, and where they touch what this claim reads (a link of a linkage placed by t), their push \
                 is not in it: place them through a coordinate to have them counted",
                self.held_by_nothing.join(", ")
            )
        })
    }

    /// Why the couplings' tensions at moment m are not this motion's, if
    /// they are not: they are not found, or a balance they rest on is not
    /// (a body with no mass, a coordinate moving it that no actuator
    /// supplies).
    pub fn tension_refusal(&self, m: &Moment, c: usize) -> Option<String> {
        self.tension_refusal_in(m, &[self.trains[c]])
    }

    /// `tension_refusal`, and that a coupling of its train passing power
    /// gives no efficiency (its loss, in the tensions of the stages driving
    /// it, unknown): a tension's claim.
    pub fn tension_refusal_strict(&self, m: &Moment, c: usize) -> Option<String> {
        self.tension_refusal(m, c).or_else(|| self.loss_unknown(m, &[self.trains[c]]))
    }

    /// Why what rests on the couplings of `trains` at moment m rests on a
    /// loss not given: one of them passes power and gives no efficiency.
    fn loss_unknown(&self, m: &Moment, trains: &[usize]) -> Option<String> {
        let c = m.lossless.iter().find(|&&c| trains.contains(&self.trains[c]))?;
        if self.couplings[*c].efficiency.is_some() {
            return Some(format!("{}, so its loss, whose friction the reactions and tensions of its train carry, is unknown", self.loss_on_a_carrier(*c)));
        }
        Some(format!(
            "{} passes power here and gives no efficiency, so its loss, whose friction the reactions and tensions of its train carry, is unknown \
             (belt-drive, chain-drive, belt-span, the gear meshes and lead-screw take a published one last)",
            self.couplings[*c].name
        ))
    }

    /// `tension_refusal` of the couplings of `trains`: the balances giving
    /// their tensions, not every coupling's (machines sharing no coupling
    /// stand apart).
    fn tension_refusal_in(&self, m: &Moment, trains: &[usize]) -> Option<String> {
        if let Some(why) = m.unsolved.as_ref().or(m.coupling.as_ref()) {
            return Some(why.clone());
        }
        let mut without: Vec<&str> = Vec::new();
        for b in self.balances.iter().filter(|b| (0..self.couplings.len()).any(|c| self.couplings[c].side[b.reaction] != 0 && trains.contains(&self.trains[c]))) {
            for label in &self.reactions[b.reaction].without_mass {
                if !without.contains(&label.as_str()) {
                    without.push(label.as_str());
                }
            }
        }
        if !without.is_empty() {
            return Some(format!("bodies carried by its balances have no mass: {} ({NO_MASS})", without.join(", ")));
        }
        let mut blocked: Vec<&str> = Vec::new();
        for c in (0..self.couplings.len()).filter(|&c| trains.contains(&self.trains[c])) {
            for label in self.blocks_coupling(c) {
                if !blocked.contains(&label) {
                    blocked.push(label);
                }
            }
        }
        if !blocked.is_empty() {
            return Some(format!("these moving bodies have no mass: {}: {NO_MASS}", self.massless_said(&blocked)));
        }
        let coords: Vec<usize> = self
            .balances
            .iter()
            .filter(|b| (0..self.couplings.len()).any(|c| self.couplings[c].side[b.reaction] != 0 && trains.contains(&self.trains[c])))
            .flat_map(|b| self.reactions[b.reaction].coords.iter().copied())
            .collect();
        self.unpowered(&coords, m)
    }

    /// Where the motion reaction r rests on is a fixed-effort drive's whose
    /// required effort differs from its rating within the rounding of its
    /// load (its effort claim unknown there: Moment::fixed_margin_unresolved),
    /// the shortfall's share of that rounding, summed over such drives: the
    /// real motion may differ from the one read by that share of what the
    /// drive's terms round to, so the reaction by that share of what its own
    /// do. 0 where there is none.
    pub fn shortfall(&self, r: usize, m: &Moment) -> f64 {
        let trains: Vec<usize> = (0..self.couplings.len()).filter(|&c| self.reaches(c, r)).map(|c| self.trains[c]).collect();
        self.shortfall_in(m, &self.reactions[r].coords, &trains)
    }

    /// `shortfall` of the balances of `trains` and coordinates `coords`.
    fn shortfall_in(&self, m: &Moment, coords: &[usize], trains: &[usize]) -> f64 {
        let mut coords = coords.to_vec();
        coords.extend(
            self.balances
                .iter()
                .filter(|b| (0..self.couplings.len()).any(|c| self.couplings[c].side[b.reaction] != 0 && trains.contains(&self.trains[c])))
                .flat_map(|b| self.reactions[b.reaction].coords.iter().copied()),
        );
        (0..self.actuators.len()).filter(|&a| m.fixed_margin_unresolved(a) && self.actuators[a].coupled.iter().any(|k| coords.contains(k))).map(|a| -m.margin(a) / m.noise[a]).sum()
    }

    /// Coupling k's tensions' band at moment m: where a fixed drive's
    /// shortfall reaches them (`shortfall`), its share of their rounding,
    /// the balances giving them; and what the bodies' measured masses can
    /// move each by. None where neither does.
    pub fn tension_shortfall(&self, m: &Moment, k: usize) -> Option<Vec<f64>> {
        let share = self.shortfall_in(m, &[], &[self.trains[k]]).max(0.0);
        let c = &self.couplings[k];
        let lines = c.first..c.first + c.lines.len();
        // And what its bodies' measured masses can move each line by (Moment::mass_t).
        let mass = |l: usize| m.mass_t.get(l).copied().unwrap_or(0.0);
        let band: Vec<f64> = lines.clone().map(|l| share * (ROUNDING * m.tensions[l].abs() + m.tension_noise[l]) + mass(l)).collect();
        band.iter().any(|x| *x != 0.0).then_some(band)
    }

    /// Whether a fixed-effort drive (its limits one number) moves what
    /// reaction r rests on: its read may carry a shortfall's band.
    pub fn fixed_reaches(&self, r: usize) -> bool {
        let coords = &self.reactions[r].coords;
        self.actuators.iter().any(|a| matches!((&a.lo, &a.hi), (Num::C(lo), Num::C(hi)) if lo == hi) && a.coupled.iter().any(|k| coords.contains(k)))
            || (0..self.couplings.len()).any(|c| self.reaches(c, r)) && self.actuators.iter().any(|a| matches!((&a.lo, &a.hi), (Num::C(lo), Num::C(hi)) if lo == hi))
    }

    /// The first actuator whose effort rests on one of `coords` and whose
    /// limits cross at moment m, saying so: no effort lies within them (a
    /// drive past its published curve's last speed), so
    /// the motion as written is past what it can do, not unpowered.
    fn crossed(&self, coords: &[usize], m: &Moment) -> Option<String> {
        let b = (0..self.actuators.len()).find(|&b| m.lo[b] > m.hi[b] && self.actuators[b].coupled.iter().any(|k| coords.contains(k)))?;
        self.crossed_at(b, m)
    }

    /// Actuator b's limits crossing at moment m, said, if they do.
    fn crossed_at(&self, b: usize, m: &Moment) -> Option<String> {
        if !(m.lo[b] > m.hi[b]) {
            return None;
        }
        let a = &self.actuators[b];
        let speed = |rate: f64| if a.turn { format!("{:.1} rpm ({rate:.6} deg/s)", rate.abs() / 6.0) } else { format!("{rate:.6} mm/s") };
        let ends = match self.crossing(b, m) {
            Some(r) => format!("ends at {}", speed(r)),
            None => "ends below that speed".into(),
        };
        Some(format!(
            "the joint of {} moves at {} here, where its effort limits (at {}) cross, no effort lying within them: its published \
             torque-speed curve {ends} (a drive past its curve's last speed, tools/drive.ink), so the motion as written is faster than its \
             data goes",
            a.path,
            speed(m.rate[b]),
            a.source
        ))
    }

    /// The rate of actuator b at which its limits begin to cross, the
    /// mechanism at moment m's position and acceleration and every speed
    /// scaled alike (by bisection, in its model of states: `twin`). None
    /// where there is none to read it in, or its limits cross at rest.
    fn crossing(&self, b: usize, m: &Moment) -> Option<f64> {
        let twin = match &self.states {
            Some(_) => self,
            None => {
                let (reload, target, built) = self.twin.as_ref()?;
                ink_kernel::once::get_or_compute(built, || crate::build_states(reload, target, self.gravity).ok()).as_ref()?
            }
        };
        let tb = twin.actuators.iter().position(|x| x.path == self.actuators[b].path)?;
        let names = |m: &Model| m.coords.iter().map(|c| c.name.clone()).collect::<Vec<_>>();
        if names(twin) != names(self) {
            return None;
        }
        let crosses = |f: f64| -> Option<bool> {
            let state: Vec<[f64; 3]> = m.q.iter().map(|c| [c[0], f * c[1], c[2]]).collect();
            let (limits, _) = twin.limits_at(m.t, &state).ok()?;
            Some(limits[tb].0 > limits[tb].1)
        };
        if crosses(0.0)? || !crosses(1.0)? {
            return None;
        }
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if crosses(mid)? { hi = mid } else { lo = mid }
        }
        Some(lo * m.rate[b])
    }

    /// The first of `coords` where no actuator supplies the effort the
    /// motion needs along it, saying so.
    fn unpowered(&self, coords: &[usize], m: &Moment) -> Option<String> {
        if let Some(why) = self.crossed(coords, m) {
            return Some(why);
        }
        // The inverse solve finds what a fixed drive would have to
        // supply, even when impossible. A tension or reaction premised
        // on that motion cannot use the hypothetical effort as supplied;
        // nor where the two differ within the rounding of the drive's load,
        // where whether it supplies the motion is undecided (its effort
        // claim is unknown there).
        // Where they differ within the rounding of the drive's load, the
        // motion is the drive's to that rounding: what rests on it is known
        // to its share of its own (Model::shortfall).
        if let Some((a, act)) = self.actuators.iter().enumerate().find(|(a, act)| {
            m.lo[*a] == m.hi[*a] && m.effort[*a].is_finite() && m.margin(*a) < 0.0 && !m.fixed_margin_unresolved(*a) && act.coupled.iter().any(|k| coords.contains(k))
        }) {
            return Some(format!("the required effort of {} is {}, different from its fixed rating of {}: the motion's balance is undetermined", act.path, m.effort[a], m.lo[a]));
        }
        coords.iter().find_map(|&k| m.unpowered[k].map(|u| self.unsupplied(m, k, u)))
    }

    /// A non-finite balance has no numerical residual to report.
    pub fn unsupplied(&self, m: &Moment, k: usize, residual: f64) -> String {
        if residual.is_finite() {
            format!("the prescribed pose or motion needs {residual:.6} N mm per unit of {} that no actuator supplies", self.coords[k].name)
        } else if m.pm.iter().any(|row| !row[k].is_finite()) {
            UNCERTAIN_DRIVES.into()
        } else {
            format!("the force balance along {} is not finite", self.coords[k].name)
        }
    }

    /// Bodies with no mass, by their labels, as they are said: each definition once, with where it is placed.
    pub fn massless_said(&self, labels: &[&str]) -> String {
        let mut by: Vec<(&str, Vec<&str>)> = Vec::new();
        for (label, def, _) in self.massless.iter().filter(|(label, ..)| labels.contains(&label.as_str())) {
            match by.iter_mut().find(|(d, _)| *d == def.as_str()) {
                Some((_, at)) => at.push(label),
                None => by.push((def, vec![label])),
            }
        }
        by.iter().map(|(def, at)| format!("{def} ({})", at.join(", "))).collect::<Vec<_>>().join("; ")
    }

    /// Why reaction r is found at no moment, if it is not: a body it
    /// carries has no mass, or is not proven bounded, or its shape moves.
    pub(super) fn static_refusal(&self, r: usize) -> Option<Refused> {
        let re = &self.reactions[r];
        if !re.without_mass.is_empty() {
            return Some(Refused::lacks(format!("bodies it carries have no mass: {} ({NO_MASS})", re.without_mass.join(", "))));
        }
        let blocked = self.blocks_reaction(r);
        if !blocked.is_empty() {
            return Some(Refused::lacks(format!("these moving bodies have no mass: {}: {NO_MASS}", self.massless_said(&blocked))));
        }
        if !re.unbounded.is_empty() {
            return Some(Refused::lacks(format!(
                "bodies it carries are not proven bounded, so their mass would be the window's: {} (bound each, e.g. intersect it with a box, or place it outside what the path carries)",
                re.unbounded.join(", ")
            )));
        }
        if !re.leaves_out.is_empty() {
            return Some(Refused::outside(format!("it leaves out bodies it carries whose shape moves, whose inertia is not counted: {}", re.leaves_out.join(", "))));
        }
        None
    }
}
