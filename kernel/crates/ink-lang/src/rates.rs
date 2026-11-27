//! Derivatives through what a form provides. Graph::diff holds every
//! param (a partial); a param some form sets (a simulated motion's, a
//! measurement's) has no rate here: no form gives equations for one. A
//! derivative taken of a param no form sets yet remembers that it held it
//! (`held`), and a form read later that makes it vary is refused
//! (`varies_held`): the derivative would have been wrong.

use crate::Doc;
use crate::graph::{BinOp, NodeId, Var};

impl Doc {
    /// The derivative of `root` along `var`: its partial (every param held,
    /// Graph::derivative) plus, for each solved param it reads that varies
    /// along `var`, its partial along that param times the param's rate
    /// (`rate`): the chain rule through solves.
    pub fn derivative(&mut self, root: NodeId, var: Var) -> NodeId {
        if var == Var::T {
            return self.speed(root, &[root]);
        }
        let d = self.g.derivative(root, var);
        let d = if matches!(var, Var::X | Var::Y | Var::Z) { d } else { self.through_solves(root, var, d) };
        crate::estimate::keep_conditions(self, root, d)
    }

    /// A turn's speed, d deg / dt (`derivative` along t), marked as the
    /// speed of the turn's sine and cosine (Graph::speed_read): an angle
    /// wrapped at a whole turn moves the turn by none, so its speed is the
    /// formula's, and only a jump the turn sees is unbounded speed.
    pub fn turn_rate(&mut self, deg: NodeId) -> NodeId {
        let read = crate::turn_sin_cos(&mut self.g, deg);
        self.speed(deg, &read)
    }

    /// Mark the total rate, before a zero partial can erase the only
    /// evidence that the drive switches on a moving solved coordinate.
    fn speed(&mut self, drive: NodeId, read: &[NodeId]) -> NodeId {
        let partial = self.g.diff(drive, Var::T);
        let rate = self.through_solves(drive, Var::T, partial);
        let rate = if self.g.depends_on(drive, |v| v == Var::T) || self.solves_read(&[drive]).iter().any(|&i| self.solves[i].on_t) { self.g.mark_speed(rate, read) } else { rate };
        crate::estimate::keep_conditions(self, drive, rate)
    }

    /// `d`, root's partial along `var`, every param held: each param root
    /// reads is noted as held there (`rate`).
    fn through_solves(&mut self, root: NodeId, var: Var, d: NodeId) -> NodeId {
        for k in self.params_read(&[root]) {
            self.rate(k, var);
        }
        d
    }

    /// The joints' speed-limit claims made with their rates through every
    /// solve their drives read (`rated`): what a value no solve sets as
    /// they were read, a solve read after them sets.
    pub(super) fn rate_joints(&mut self) {
        for (i, drive, lo, hi, turn) in std::mem::take(&mut self.rated) {
            let rate = if turn { self.turn_rate(drive) } else { self.derivative(drive, Var::T) };
            let (a, b) = (self.g.bin_real(BinOp::Sub, rate, lo), self.g.bin_real(BinOp::Sub, hi, rate));
            let claim = self.g.min(a, b);
            let x = &mut self.asserts[i];
            x.claim = crate::Claim::Value(claim, crate::Cmp::Ge, 0.0);
            x.terms[0].1 = rate;
        }
    }

    /// Param k's rate along `var`: none, since no form gives equations for
    /// what it sets. Where no form sets k (yet), the derivative holds it,
    /// and that is kept with where (`held`): a form read later that makes
    /// it vary along `var` is refused (`keep_form`).
    pub fn rate(&mut self, k: usize, var: Var) -> Option<NodeId> {
        if var != Var::Param(k as u32) && !self.solves.iter().any(|s| s.owns.contains(&k)) {
            let site = self.form_at.clone();
            self.held.entry((k, var)).or_insert(site);
        }
        None
    }

    /// Refuse what the solves read so far make vary that a derivative
    /// held fixed (`held`): a held param a solve now sets, varying along
    /// the derivative's input through everything that solve rests on, so a
    /// solve that leaves it fixed passes the hold on to what it reads (u =
    /// 3 m, then m = 100 t). Before the solve just read none varied, so
    /// one that does now is that solve's doing. Err: which, and where.
    pub(super) fn varies_held(&self) -> Result<(), String> {
        let mut held: Vec<(&(usize, Var), &String)> = self.held.iter().collect();
        held.sort_by(|a, b| (a.0.0, a.1).cmp(&(b.0.0, b.1)));
        for (&(k, var), site) in held {
            let Some(i) = self.solves.iter().position(|s| s.owns.contains(&k)) else { continue };
            let varies = match var {
                Var::T => self.solves[i].on_t,
                Var::Param(p) => self.solves_under(&[k]).iter().any(|&j| self.solves[j].reads.contains(&(p as usize))),
                _ => false,
            };
            if varies {
                let along = match var {
                    Var::Param(p) => format!("'{}'", self.params[p as usize].name),
                    _ => "t".into(),
                };
                return Err(format!("it makes '{}' vary along {along}, which the derivative at {site} held fixed (no solve set it there)", self.params[k].name));
            }
        }
        Ok(())
    }
}
