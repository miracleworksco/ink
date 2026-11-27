//! Placements as transforms (Xf): composing them, rotations, `at` and
//! axes, directions that must have a length (Condition), and values moved
//! into the world.

use super::*;

impl Doc {
    pub fn identity(&mut self) -> Xf {
        let (o, z) = (self.g.constant(1.0), self.g.constant(0.0));
        Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: [z, z, z] }
    }

    /// `a` then `b`: p -> Rb (Ra p + ta) + tb.
    pub fn compose(&mut self, a: &Xf, b: &Xf) -> Xf {
        let identity = self.identity();
        if *a == identity {
            return *b;
        }
        if *b == identity {
            return *a;
        }
        let mut r = a.r;
        let mut t = a.t;
        for i in 0..3 {
            for j in 0..3 {
                let mut acc = self.g.constant(0.0);
                for k in 0..3 {
                    let p = self.g.frame_bin(BinOp::Mul, b.r[i][k], a.r[k][j]);
                    acc = self.g.frame_bin(BinOp::Add, acc, p);
                }
                r[i][j] = acc;
            }
            let mut acc = b.t[i];
            for k in 0..3 {
                let p = self.g.frame_bin(BinOp::Mul, b.r[i][k], a.t[k]);
                acc = self.g.frame_bin(BinOp::Add, acc, p);
            }
            t[i] = acc;
        }
        Xf { r, t }
    }

    /// The unit vector along `v`: v / |v|. A direction must have a length,
    /// or the division by zero makes it 0 (a rotation about it a uniform
    /// scale). Where `v` reads no x, y, z (a placement's axis, a frame's
    /// direction), a zero length is refused: a constant one here, at
    /// `at`; otherwise its positive length is a condition (`conditions`),
    /// checked once the document is read at the params' values, at each
    /// moment a view is made at, and proven positive over every scope it
    /// reads (Doc::check_conditions). A `v` over space is a field: its
    /// zero is an undefined case, reported at samples like any division.
    pub(super) fn unit(&mut self, v: &[NodeId], at: &Sexp, what: &str, hint: &str, constructed: bool) -> R<Vec<NodeId>> {
        let g = &mut self.g;
        // For an axis-aligned frame the norm is the absolute component,
        // even when its square would underflow. An authored unit follows
        // the usual evaluation rules for fixed scalar expressions.
        let un = if constructed { Graph::un_real } else { Graph::un };
        let bin = if constructed { Graph::bin_real } else { Graph::bin };
        let mut nonzero = v.iter().copied().filter(|&a| g.constant_value(a) != Some(0.0));
        let single = nonzero.next().filter(|_| constructed && nonzero.next().is_none());
        let n2 = if let Some(a) = single {
            let length = un(g, UnOp::Abs, a);
            // Zero components do not contribute to the length, but its
            // use still requires the drives read to compute all of them.
            v.iter().fold(length, |n, &input| g.with_reads_of(n, input))
        } else {
            let sq: Vec<_> = v.iter().map(|&a| un(g, UnOp::Square, a)).collect();
            sq[1..].iter().fold(sq[0], |acc, &s| bin(g, BinOp::Add, acc, s))
        };
        if !self.g.depends_on(n2, |x| matches!(x, Var::X | Var::Y | Var::Z)) {
            match self.g.as_const(n2) {
                Some(c) if c > 0.0 && c.is_finite() => {}
                Some(c) => {
                    let why = if c > 0.0 { "is too long to square in doubles" } else { "has zero length" };
                    return Err(self.err_hint(at, format!("{what} {why}: it has no direction"), Some(hint.into())));
                }
                None if self.conditions.iter().any(|c| c.len2 == n2) => {}
                None => {
                    let c = Condition { len2: n2, what: what.to_string(), site: self.form_at.clone(), at: self.err_hint(at, "", Some(hint.into())) };
                    self.conditions.push(c);
                }
            }
        }
        let g = &mut self.g;
        let n = if single.is_some() { n2 } else { un(g, UnOp::Sqrt, n2) };
        Ok(v.iter().map(|&a| bin(g, BinOp::Div, a, n)).collect())
    }

    /// Whether the directions (Doc::unit) have a length at ctx: at its
    /// moment and params' values, and, where they read them, for every
    /// value of its ranged params (ranged::value_span, as a value claim)
    /// and every moment of `motion` (value_span with t ranged over it). Err: the first that does not, naming it and where it is
    /// zero (a witness), or that its length could not be proven nonzero
    /// over the scope. `at_load`: only those not depending on t (load-time
    /// values are at t = 0; a moment is checked when a view is made at it).
    /// Solves that failed at ctx (`failed`) are skipped: what rests on them
    /// is refused already (ModelView::refusal), and over a motion so is
    /// what rests on a solve depending on t that ctx does not derive (solved
    /// again at each step); one it derives from t's param over the span
    /// (motion::scope_over) is checked over it. `of`: only the directions
    /// these fields read that rest on a solve on t (a claim's, in its
    /// scope over a motion: claims::ask_in).
    pub fn check_conditions(&self, ctx: &Ctx, motion: Option<(f64, f64)>, failed: &[usize], at_load: bool, of: Option<&[u32]>) -> R<()> {
        use ink_kernel::ranged::{Cmp, Comparison, ProofOptions, Status, value_span};
        if of.is_some() || (!self.claim_materials && motion.is_none() && ctx.p_spans.is_empty()) {
            self.check_magnitudes(ctx, motion, failed, at_load, of)?;
        }
        if self.conditions.is_empty() {
            return Ok(());
        }
        let read: Option<HashSet<u32>> = of.map(|r| self.g.reachable(r).into_iter().collect());
        for c in &self.conditions {
            let solves = self.solves_read(&[c.len2]);
            let on_t: Vec<usize> = solves.iter().copied().filter(|&i| self.solves[i].on_t).collect();
            let solved_on_t = !on_t.is_empty();
            let derived = solved_on_t && on_t.iter().all(|&i| self.solves[i].owns.iter().all(|&k| ctx.p_spans.is_derived(k as u32)));
            if read.as_ref().is_some_and(|r| !solved_on_t || !r.contains(&c.len2)) {
                continue;
            }
            let reads_t = solved_on_t || self.g.depends_on(c.len2, |v| matches!(v, Var::T));
            if solves.iter().any(|i| failed.contains(i)) || (at_load && reads_t) || (motion.is_some() && solved_on_t && !derived) {
                continue;
            }
            let ranged = self.g.depends_on(c.len2, |v| matches!(v, Var::Param(k) if ctx.p_spans.get(k).is_some()));
            let tape = Tape::new(&self.g, c.len2);
            let err = |message: String| InkError { message, ..c.at.clone() };
            let named = |t: f64, values: &[(usize, f64)]| -> String {
                let mut at: Vec<String> = values.iter().map(|&(k, v)| format!("{} = {v}", self.params[k].name)).collect();
                // t's own param (a motion's scope) names the moment.
                if reads_t && !values.iter().any(|&(k, _)| self.params[k].name == "t") {
                    at.insert(0, format!("t = {t}"));
                }
                if at.is_empty() { String::new() } else { format!(" at {}", at.join(", ")) }
            };
            let over = motion.filter(|_| reads_t);
            if let Some(span) = over {
                // t over the span, an axis of the one search as a ranged
                // param is (with the params ranged, where they are).
                let opts = ProofOptions { res: None, finest: 0.0, depth: ink_kernel::scope::PARAM_SPLITS, tol: 0.0, boxes: ink_kernel::query::BOXES };
                let over = Ctx { t: ink_core::eval::span_mid(span), t_span: Some(span), ..ctx.clone() };
                let v = value_span(&tape, &over, [0.0; 3], Comparison::exact(Cmp::Gt, 0.0), &opts);
                // A length bounded above too: over the span at once, else over it in pieces (a turn of a
                // whole revolution bounds no better than its circle's box at once, and may overflow).
                let bounded = || {
                    v.range.1.is_finite()
                        || (0..256).all(|i| {
                            let piece = (span.0 + (span.1 - span.0) * i as f64 / 256.0, span.0 + (span.1 - span.0) * (i + 1) as f64 / 256.0);
                            let at = Ctx { t: ink_core::eval::span_mid(piece), t_span: Some(piece), ..over.clone() };
                            eval::eval_interval_scope(&tape, &at, [ink_core::interval::I::point(0.0); 3], None).hi.is_finite()
                        })
                };
                match (v.status, v.witness) {
                    (Status::Proven, _) if bounded() => {}
                    (Status::Disproven, Some((at, len2))) => {
                        let values: Vec<(usize, f64)> = at.params.iter().map(|&(k, x)| (k as usize, x)).collect();
                        return Err(err(format!(
                            "{} at {} has zero length{} (squared length {len2}), within t in [{}, {}]: it has no direction there",
                            c.what,
                            c.site,
                            named(at.t, &values),
                            span.0,
                            span.1
                        )));
                    }
                    _ => {
                        return Err(err(format!(
                            "{} at {} is not proven to have a length at every moment of t in [{}, {}]: a direction must have one at every moment",
                            c.what, c.site, span.0, span.1
                        )));
                    }
                }
            } else if ranged {
                let opts = ProofOptions { res: None, finest: 0.0, depth: ink_kernel::scope::PARAM_SPLITS, tol: 0.0, boxes: ink_kernel::query::BOXES };
                let v = value_span(&tape, ctx, [0.0; 3], Comparison::exact(Cmp::Gt, 0.0), &opts);
                let scope: Vec<String> = ctx
                    .p_spans
                    .iter()
                    .filter(|s| self.params_read(&[c.len2]).contains(&(s.k as usize)))
                    .map(|s| format!("{} in [{}, {}]", self.params[s.k as usize].name, s.span.0, s.span.1))
                    .collect();
                let scope = scope.join(", ");
                match v.status {
                    Status::Proven if v.range.1.is_finite() => {}
                    Status::Disproven => {
                        let (at, len2) = v.witness.expect("a disproof has a witness");
                        let values: Vec<(usize, f64)> = at.params.iter().map(|&(k, x)| (k as usize, x)).collect();
                        return Err(err(format!(
                            "{} at {} has zero length{} (squared length {len2}), within {scope}: it has no direction there",
                            c.what,
                            c.site,
                            named(at.t, &values)
                        )));
                    }
                    // Where the bounds reach 0 but no value is found to
                    // (a zero within rounding), not proven, and refused.
                    _ => {
                        return Err(err(format!(
                            "{} at {} is not proven to have a length over {scope} (its positive-length check is bounded within [{}, {}] there): a direction must have one at every value",
                            c.what, c.site, v.range.0, v.range.1
                        )));
                    }
                }
            } else {
                let len2 = ink_core::eval::eval_interval(&tape, ctx, [ink_core::interval::I::point(0.0); 3], None);
                if !(len2.lo > 0.0 && len2.hi.is_finite()) {
                    let values: Vec<(usize, f64)> = self.params_read(&[c.len2]).into_iter().map(|k| (k, ctx.params[k])).collect();
                    let why = if len2.hi == 0.0 {
                        "has zero length"
                    } else if !len2.hi.is_finite() {
                        "is too long to square in doubles"
                    } else {
                        "is not proven to have a nonzero length"
                    };
                    return Err(err(format!("{} at {} {why}{}: it has no direction", c.what, c.site, named(ctx.t, &values))));
                }
            }
        }
        Ok(())
    }

    pub(super) fn rotation(&mut self, axis: [NodeId; 3], deg: NodeId, at: &Sexp) -> R<Xf> {
        let hint = "(rot ax ay az deg) turns about the axis (ax, ay, az), which needs a nonzero length. If it reads a param or t, keep it nonzero over their ranges";
        let u = self.unit(&axis, at, &format!("the axis of {at}"), hint, true)?;
        let g = &mut self.g;
        let [s, c] = turn_sin_cos(g, deg);
        let one = g.constant(1.0);
        let omc = g.frame_bin(BinOp::Sub, one, c);
        // Rodrigues: R = c I + s [u]x + (1 - c) u u^T. A diagonal entry is
        // written u_i^2 + c (1 - u_i^2), c read once: along a constant axis
        // it is exactly 1, where (1 - c) + c would read the angle (and t,
        // over a motion) in an entry that never changes.
        let mut r = [[0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let uu = g.frame_bin(BinOp::Mul, u[i], u[j]);
                let mut e;
                if i == j {
                    let rest = g.frame_bin(BinOp::Sub, one, uu);
                    let cr = g.frame_bin(BinOp::Mul, c, rest);
                    e = g.frame_bin(BinOp::Add, uu, cr);
                } else {
                    e = g.frame_bin(BinOp::Mul, omc, uu);
                    // [u]x[i][j] = -eps_ijk u_k
                    let kk = 3 - i - j;
                    let sign = if (i + 1) % 3 == j { -1.0 } else { 1.0 };
                    let sgn = g.constant(sign);
                    let su = g.frame_bin(BinOp::Mul, s, u[kk]);
                    let term = g.frame_bin(BinOp::Mul, sgn, su);
                    e = g.frame_bin(BinOp::Add, e, term);
                }
                r[i][j] = e;
            }
        }
        let z = g.constant(0.0);
        Ok(Xf { r, t: [z, z, z] })
    }

    /// `(at STEP ...)`: steps apply in order; a step is (move dx dy dz),
    /// (rot-x deg), (rot-y deg), (rot-z deg), (rot ax ay az deg) or any
    /// frame. A placement reads no x, y or z (frame_value).
    pub(super) fn eval_at(&mut self, steps: &[Sexp], env: &mut Env) -> R<Xf> {
        let usage = "usage: (at step ...) with steps (move dx dy dz) (rot-x deg) (rot-y deg) (rot-z deg) (rot ax ay az deg) or a frame, applied in order";
        let mut xf = self.identity();
        for step in steps {
            let head = step.list().and_then(|l| l.first()).and_then(|h| h.sym()).unwrap_or("");
            let s = if matches!(head, "move" | "rot-x" | "rot-y" | "rot-z" | "rot") {
                let sp = step.list().unwrap();
                let ps = match head {
                    "move" => "dx dy dz",
                    "rot" => "ax ay az deg",
                    _ => "deg",
                };
                if sp.len() - 1 != ps.split(' ').count() {
                    return Err(self.err_hint(step, format!("the step '{head}' takes {} number(s) ({ps}), got {}", ps.split(' ').count(), sp.len() - 1), Some(usage.into())));
                }
                // A move's numbers are lengths, a turn's an angle (bare: mm, degrees): a time or a force is neither.
                let mut nums = Vec::new();
                for (k, a) in sp[1..].iter().enumerate() {
                    let mut v = self.eval(a, env)?;
                    match head {
                        "move" => self.fits(&mut v, a, Quantity::Length, "move")?,
                        "rot" if k < 3 => {}
                        _ => self.fits(&mut v, a, Quantity::Angle, head)?,
                    }
                    nums.push(self.as_field(v, a)?);
                }
                let (z, o) = (self.g.constant(0.0), self.g.constant(1.0));
                match (head, nums.as_slice()) {
                    ("move", [dx, dy, dz]) => {
                        let mut m = self.identity();
                        m.t = [*dx, *dy, *dz];
                        m
                    }
                    ("rot-x", [d]) => self.rotation([o, z, z], *d, step)?,
                    ("rot-y", [d]) => self.rotation([z, o, z], *d, step)?,
                    ("rot-z", [d]) => self.rotation([z, z, o], *d, step)?,
                    ("rot", [ax, ay, az, d]) => self.rotation([*ax, *ay, *az], *d, step)?,
                    _ => unreachable!("counted above"),
                }
            } else {
                match self.eval(step, env).map_err(|err| self.misplaced(err, step, usage))? {
                    Value::Frame(f) => f,
                    other => return Err(self.err_hint(step, format!("bad placement step: a {}, not a step or a frame", other.kind()), Some(usage.into()))),
                }
            };
            xf = self.compose(&xf, &s);
        }
        Ok(xf)
    }

    /// `(axes (origin ox oy oz) (x xx xy xz) (z zx zy zz))`: a frame as ink
    /// prints one, its clauses in any order. z keeps its direction (a
    /// face's normal); x is made square to it, its part along z dropped
    /// (refused where nothing is left: x along z); y is z cross x. A
    /// placement reads no x, y or z (frame_value).
    pub(super) fn eval_axes(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Xf> {
        let usage = "usage: (axes (origin ox oy oz) (x xx xy xz) (z zx zy zz)): the frame at origin with those directions, as ink check prints one \
                     (each clause 3 numbers or a vec of 3); z keeps its direction, x is made square to it, y is z cross x";
        let mut got: [Option<[NodeId; 3]>; 3] = [None; 3];
        for a in args {
            let (k, nums) = match a.list() {
                Some([k, nums @ ..]) => (k.sym().unwrap_or(""), nums),
                _ => ("", &[][..]),
            };
            let Some(slot) = ["origin", "x", "z"].iter().position(|&n| n == k) else {
                let why = if k == "y" { "(axes ...) takes no y: it is z cross x".to_string() } else { format!("bad clause {a} of (axes ...)") };
                return Err(self.err_hint(a, why, Some(usage.into())));
            };
            let v = match nums {
                [v] => self.fields(v, env)?,
                _ => nums.iter().map(|n| self.field(n, env)).collect::<R<Vec<_>>>()?,
            };
            if got[slot].is_some() || v.len() != 3 {
                let why = if v.len() != 3 { format!("({k} ...) takes 3 numbers or a vec of 3") } else { format!("({k} ...) is given twice") };
                return Err(self.err_hint(a, why, Some(usage.into())));
            }
            got[slot] = Some([v[0], v[1], v[2]]);
        }
        let [Some(o), Some(x), Some(z)] = got else { return Err(self.err_hint(e, "(axes ...) needs origin, x and z", Some(usage.into()))) };
        // z cross x is the normal to their plane. Form it before
        // normalizing z: a rounded projection along z can swamp a thin
        // perpendicular component. Enclose each product and difference.
        let cross = |g: &mut Graph, a: &[NodeId], b: &[NodeId]| -> [NodeId; 3] {
            std::array::from_fn(|i| {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                let (p, q) = (g.bin_real(BinOp::Mul, a[j], b[k]), g.bin_real(BinOp::Mul, a[k], b[j]));
                g.bin_real(BinOp::Sub, p, q)
            })
        };
        // Canonical-axis inputs can be scaled exactly before the product.
        // This also keeps the product of two tiny axial directions nonzero.
        let axial = |g: &mut Graph, v: [NodeId; 3]| {
            let mut components = v.iter().filter(|&&a| g.as_const(a) != Some(0.0));
            match components.next().and_then(|&a| g.as_const(a)) {
                Some(c) if c.is_finite() && components.next().is_none() => v.map(|a| if g.as_const(a) == Some(0.0) { a } else { g.constant(c.signum()) }),
                _ => v,
            }
        };
        let x = axial(&mut self.g, x);
        let z = axial(&mut self.g, z);
        let normal = cross(&mut self.g, &z, &x);
        let z =
            self.unit(&z, e, &format!("the z of {e}"), "z is the frame's normal, which needs a nonzero length. If it reads a param or t, keep it nonzero over their ranges", true)?;
        let y =
            self.unit(&normal, e, &format!("the x of {e}, square to its z,"), "x must not be parallel to z: x is made square to z, and of an x along z nothing is left", true)?;
        let x = cross(&mut self.g, &y, &z);
        Ok(Xf { r: std::array::from_fn(|i| [x[i], y[i], z[i]]), t: o })
    }

    /// Express a local-space value in world space: substitute p_local = R^T (p - t).
    pub fn to_world(&mut self, v: &Value, xf: &Xf) -> Value {
        let map = self.inverse_map(xf);
        let w = self.placed_value(v, xf, &map);
        if let (Value::Field(local, _), Value::Field(world, _)) = (v, &w) {
            self.record_placed(*world, *local, xf);
        }
        w
    }

    /// Record that `world` is `local` placed at `at` (`placed`). Equal
    /// fields are one node, so an entry is one true decomposition of the
    /// field, not the item's: a reader takes the shape and the placement
    /// from the same entry. One whose placement changes nothing (world is
    /// local) says nothing and would hide one that does.
    pub(super) fn record_placed(&mut self, world: NodeId, local: NodeId, at: &Xf) {
        if world != local {
            self.placed.insert(world, Placed { shape: local, at: *at });
        }
    }

    /// `v` put at `xf`: a field warped by its inverse (`map`), a point
    /// moved by it (R p + T), a vec each of its items.
    fn placed_value(&mut self, v: &Value, xf: &Xf, map: &[(Var, NodeId)]) -> Value {
        match v {
            Value::Point(p) => Value::Point(std::array::from_fn(|i| {
                let mut acc = xf.t[i];
                for j in 0..3 {
                    let m = self.g.frame_bin(BinOp::Mul, xf.r[i][j], p[j]);
                    acc = self.g.frame_bin(BinOp::Add, acc, m);
                }
                acc
            })),
            Value::Vec(items) => Value::list(items.iter().map(|i| self.placed_value(i, xf, map)).collect()),
            other => self.warp_value(other, map),
        }
    }

    pub(crate) fn inverse_map(&mut self, xf: &Xf) -> Vec<(Var, NodeId)> {
        let p = [self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z)];
        let d: Vec<NodeId> = (0..3).map(|j| self.g.sub(p[j], xf.t[j])).collect();
        let mut out = Vec::new();
        for (i, var) in [Var::X, Var::Y, Var::Z].into_iter().enumerate() {
            let mut acc = self.g.constant(0.0);
            for j in 0..3 {
                let m = self.g.mul(xf.r[j][i], d[j]);
                acc = self.g.add(acc, m);
            }
            out.push((var, acc));
        }
        out
    }

    /// The subexpressions of `f` named `name` (by a def or a tag), among
    /// the pieces its value is built of (Graph::pieces).
    pub(super) fn named_in(&mut self, f: NodeId, name: &str) -> Vec<NodeId> {
        self.g
            .pieces(&[f])
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Tag(t, a) if self.g.tag_names[t as usize] == name => Some(self.g.with_reads_of(a, f)),
                _ => None,
            })
            .collect()
    }

    /// Every construction name inside `f`.
    pub(super) fn named_in_all(&self, f: NodeId) -> Vec<&str> {
        let mut names: Vec<&str> = self
            .g
            .pieces(&[f])
            .into_iter()
            .filter_map(|n| match self.g.nodes[n as usize] {
                crate::graph::Node::Tag(t, _) => Some(self.g.tag_names[t as usize].as_str()),
                _ => None,
            })
            .collect();
        names.sort_unstable();
        names.dedup();
        names
    }

    pub(super) fn warp_value(&mut self, v: &Value, map: &[(Var, NodeId)]) -> Value {
        match v {
            Value::Field(f, said) => Value::Field(self.g.subst(*f, map), said.clone()),
            Value::Vec(items) => Value::list(items.iter().map(|i| self.warp_value(i, map)).collect()),
            other => other.clone(),
        }
    }
}
