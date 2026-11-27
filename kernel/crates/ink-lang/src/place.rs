//! Where instances are: an instance's placement and joint, frames, and the
//! relations and clauses of an assembly placed in the world.

use super::*;

impl Doc {
    /// An instance's placement from its clauses after the child: WHERE,
    /// then optionally CFRAME (a frame the child declares, met face to face
    /// with WHERE), MOTION and limits. Steps in order: without CFRAME
    /// (MOTION, WHERE); with it (C^-1, FLIP, MOTION, WHERE), FLIP the exact
    /// half-turn about x. Limits add the claim `limits:ASSEMBLY/NAME`.
    pub(super) fn instance_placement(&mut self, c: &Sexp, iname: &str, child: &str, rest: &[Sexp], facts: &mut Vec<Stated>) -> R<(Xf, Option<Joint>, Vec<String>, Option<String>)> {
        let Some((w, after)) = rest.split_first() else { return Ok((self.identity(), None, Vec::new(), None)) };
        // `(fact ...)` is the language's clause unless the file defines a `fact` of its own (a function giving a frame).
        let fact_clause = |d: &Doc, x: &Sexp| x.head_args().is_some_and(|(head, _)| head == "fact") && d.scoped("fact").is_none();
        // An instance that states a fact before any WHERE is placed as one with no clauses is.
        let (anchor, on, rest) = if fact_clause(self, w) {
            (self.identity(), Vec::new(), rest)
        } else {
            self.building.as_mut().unwrap().read.clear();
            let anchor = self.frame_value(w, &mut self.arg_env(), "WHERE")?;
            (anchor, std::mem::take(&mut self.building.as_mut().unwrap().read), after)
        };
        // A `for` among the clauses repeats them, its list read at load time
        // with the assembly's arguments: a clause one placement has and
        // another lacks (an actuator's limits only at a supply its record
        // publishes them for).
        let rest = &self.expand_for(rest)?;
        let rest = &self.optional_clauses(rest, &mut self.arg_env(), |p| p.joint_clauses(), INSTANCE, Some(iname))?;
        let (mut cframe, mut motion) = (None, None);
        // (limits LO HI) and (speed-limits LO HI): bounds on the drive and
        // on its rate, each with the accuracy it states, (accuracy x).
        const BOUNDS: [&str; 2] = ["limits", "speed-limits"];
        let mut bounds: [Option<(NodeId, NodeId, String, String, Option<&Sexp>)>; 2] = [None, None];
        let mut written: [Option<[(Value, &Sexp); 2]>; 2] = [None, None];
        // Packages' joint clauses, read once the motion is known.
        let mut owned: Vec<(&'static dyn ext::Package, &'static str, &Sexp)> = Vec::new();
        for x in rest {
            let head = x.list().and_then(|l| l.first()).and_then(|h| h.sym());
            match (x.sym(), head, x.list().unwrap_or(&[])) {
                (Some(n), ..) if cframe.is_none() && motion.is_none() && bounds.iter().all(Option::is_none) && owned.is_empty() => cframe = Some((n.to_string(), x)),
                (None, Some("fact"), _) if fact_clause(self, x) => facts.push(self.stated(x, child)?),
                (None, Some("turn" | "slide"), [m, a, v]) => {
                    if motion.is_some() {
                        return Err(self.err(x, "one motion per instance; put fixed offsets in WHERE, e.g. (at (move 0 0 -8) (frame-of m shaft))"));
                    }
                    let axis = match a.sym() {
                        Some("x") => 0,
                        Some("y") => 1,
                        Some("z") => 2,
                        _ => return Err(self.err(a, format!("the axis of {} must be x, y or z: an axis of the frame WHERE places the instance at", m.sym().unwrap()))),
                    };
                    let mut value = self.eval(v, &mut self.arg_env())?;
                    self.drives(&mut value, v, head == Some("turn"))?;
                    let drive = self.as_field(value, v)?;
                    if self.g.depends_on(drive, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                        return Err(self.err(v, "placements cannot depend on x, y or z"));
                    }
                    motion = Some((head == Some("turn"), axis, drive, v.to_string()));
                }
                (None, Some(k @ ("limits" | "speed-limits")), [_, lo, hi, tol @ ..]) if tol.len() <= 1 => {
                    let slot = BOUNDS.iter().position(|b| *b == k).unwrap();
                    if bounds[slot].is_some() {
                        return Err(self.err(x, format!("one ({k} LO HI) per instance")));
                    }
                    let (vl, vh) = (self.eval(lo, &mut self.arg_env())?, self.eval(hi, &mut self.arg_env())?);
                    let (l, h) = (self.as_field(vl.clone(), lo)?, self.as_field(vh.clone(), hi)?);
                    let fixed = |v: Var| matches!(v, Var::X | Var::Y | Var::Z | Var::T);
                    if self.g.depends_on(l, fixed) || self.g.depends_on(h, fixed) {
                        return Err(self.err(x, format!("{k} are numbers or params, not fields of x y z or t")));
                    }
                    if let (Some(a), Some(b)) = (self.g.as_const(l), self.g.as_const(h))
                        && a > b
                    {
                        return Err(self.err(x, format!("{k} {a} > {b}: write ({k} LO HI) with LO <= HI")));
                    }
                    bounds[slot] = Some((l, h, lo.to_string(), hi.to_string(), tol.first()));
                    written[slot] = Some([(vl, lo), (vh, hi)]);
                }
                (None, Some(h), _) if let Some((package, h)) = self.registry.joint_clause(h) => {
                    if owned.iter().any(|o| o.1 == h) {
                        return Err(self.err(x, format!("one ({h} ...) per instance")));
                    }
                    owned.push((package, h, x));
                }
                (None, Some(h), _) if !["turn", "slide", "limits", "speed-limits"].contains(&h) => {
                    let msg = format!(
                        "bad instance clause: '{h}' is not the kernel's, nor from a package this ink loads. Check its name; a package's clause needs an ink that loads that package"
                    );
                    return Err(self.err_hint(x, msg, Some(INSTANCE.into())));
                }
                _ => return Err(self.err_hint(x, "bad instance clause", Some(INSTANCE.into()))),
            }
        }
        if motion.is_none() {
            if let Some(k) = bounds.iter().position(Option::is_some) {
                return Err(self.err(c, format!("({} LO HI) bounds a motion: add (turn AXIS deg) or (slide AXIS mm)", BOUNDS[k])));
            }
            if let Some((_, h, x)) = owned.first() {
                return Err(self.err(x, format!("({h} ...) is a joint's: add (turn AXIS deg) or (slide AXIS mm)")));
            }
        }
        // A slide's limits are lengths, its speed limits speeds; a turn's
        // limits are angles (a bare number: degrees), its speed limits turning rates.
        // Their accuracies are in their units.
        let (mut tols, mut quantities) = ([None; 2], [None; 2]);
        if let Some((turn, ..)) = motion {
            for (k, b) in written.iter_mut().enumerate() {
                let q = match (k, turn) {
                    (0, true) => Some(Quantity::Angle),
                    (0, false) => Some(Quantity::Length),
                    (_, true) => Some(Quantity::AngularSpeed),
                    (_, false) => Some(Quantity::Speed),
                };
                quantities[k] = q;
                for (end, (v, e)) in b.iter_mut().flatten().enumerate() {
                    if let Some(q) = q {
                        self.fits(v, e, q, BOUNDS[k])?;
                        let read = self.as_field(v.clone(), e)?;
                        let bound = bounds[k].as_mut().expect("written with its bounds");
                        *(if end == 0 { &mut bound.0 } else { &mut bound.1 }) = read;
                    }
                }
                if let Some(t) = bounds[k].as_ref().and_then(|b| b.4) {
                    tols[k] = Some(self.accuracy_clause(t, q, &mut self.arg_env())?);
                }
            }
        }
        let mut at = self.identity();
        if let Some((n, x)) = &cframe {
            let frames = self.frames_of(child);
            let Some((_, f)) = frames.iter().find(|(m, _)| m == n) else { return Err(self.err(x, Self::no_frame(child, n, &frames))) };
            let f = *f;
            let inv = self.inverse(&f);
            let m = self.g.constant(-1.0);
            let mut flip = self.identity();
            (flip.r[1][1], flip.r[2][2]) = (m, m);
            at = self.compose(&inv, &flip);
        }
        let below = at;
        let joint = match motion {
            Some((turn, axis, drive, text)) => {
                let z = self.g.constant(0.0);
                let m = if turn {
                    let mut u = [z; 3];
                    u[axis] = self.g.constant(1.0);
                    self.rotation(u, drive, c)?
                } else {
                    let mut m = self.identity();
                    m.t[axis] = drive;
                    m
                };
                at = self.compose(&at, &m);
                // Each bound is a value claim, LO <= v <= HI, on the drive
                // or on its rate, stated where the motion is set: where
                // the arguments were passed, else this clause.
                let asm = self.building.as_ref().unwrap().asm.name.clone();
                let (file, line, scope) = self.building.as_ref().unwrap().passed.clone().unwrap_or((self.file.clone(), c.line, self.scope));
                // The rate as the drive reads so far: taken again once the
                // file is read (rate_joints).
                let rate = if turn {
                    let read = turn_sin_cos(&mut self.g, drive);
                    self.g.speed_read(drive, &read)
                } else {
                    self.g.derivative(drive, Var::T)
                };
                let what = [("limits", drive, text.clone(), "drive"), ("speed-limits", rate, format!("(diff {text} t)"), "rate")];
                for ((((lo, hi, lt, ht, _), tol), quantity), (kind, v, vt, term)) in
                    bounds.iter().zip(tols).zip(quantities).zip(what).filter_map(|(((b, t), q), w)| Some((((b.as_ref()?, t), q), w)))
                {
                    if kind == "speed-limits" {
                        self.rated.push((self.asserts.len(), drive, *lo, *hi, turn));
                    }
                    let (a, b) = (self.g.bin_real(BinOp::Sub, v, *lo), self.g.bin_real(BinOp::Sub, *hi, v));
                    let claim = Claim::Value(self.g.min(a, b), Cmp::Ge, 0.0);
                    let text = format!("{lt} <= {vt} <= {ht}, as value = min({vt} - {lt}, {ht} - {vt}) >= 0");
                    let terms = vec![(term.into(), v), ("lo".into(), *lo), ("hi".into(), *hi)];
                    self.asserts.push(AssertDef {
                        name: format!("{kind}:{asm}/{iname}"),
                        claim,
                        text,
                        bound: None,
                        res: None,
                        accuracy: tol,
                        quantity,
                        file: file.clone(),
                        scope,
                        line,
                        terms,
                        fixed: Vec::new(),
                    });
                }
                let mut clauses = Vec::new();
                for (package, head, x) in owned {
                    let mut r = Reader::in_joint(self, self.arg_env(), turn);
                    let kept = package.clause(&mut r, x, head, &x.list().unwrap()[1..])?;
                    let reads = r.into_reads();
                    clauses.push(OwnedClause { package: package.name(), head, kept, reads, file: file.clone(), scope, line });
                }
                let [limits, speed_limits] = bounds.map(|b| b.map(|l| (l.0, l.1)));
                Some(Joint { turn, axis, anchor, below, drive, text, limits, speed_limits, clauses })
            }
            None => None,
        };
        let met = cframe.map(|(n, _)| n);
        Ok((self.compose(&at, &anchor), joint, on, met))
    }

    /// `(frame NAME F)` in a body or an assembly: F in its coordinates.
    pub(super) fn frame_clause(&mut self, c: &Sexp, frames: &[(String, Xf)]) -> R<(String, Xf)> {
        let [_, n, f] = c.list().unwrap_or(&[]) else { return Err(self.err(c, "usage: (frame name F), F a frame such as (at (move 0 0 5))")) };
        let n = self.name_of(n)?;
        if frames.iter().any(|(m, _)| *m == n) {
            return Err(self.err(c, format!("duplicate frame '{n}'")));
        }
        Ok((n, self.frame_value(f, &mut self.arg_env(), "F")?))
    }

    /// `e` as a placement (`what` names the slot in the error): a frame
    /// whose numbers read params, t and local names, not x, y or z, however
    /// it was made. Every placement is read here: an instance's WHERE, a
    /// frame clause, `place`.
    pub(super) fn frame_value(&mut self, e: &Sexp, env: &mut Env, what: &str) -> R<Xf> {
        let how = "a placement is a frame: (at STEP ...) with steps (move dx dy dz) (rot-x deg) (rot-y deg) (rot-z deg) (rot ax ay az deg), or (frame-of inst name)";
        match self.eval(e, env).map_err(|err| self.misplaced(err, e, how))? {
            Value::Frame(f) => {
                let numbers: Vec<NodeId> = f.r.iter().flatten().chain(&f.t).copied().collect();
                if self.g.reachable(&numbers).iter().any(|&n| matches!(self.g.nodes[n as usize], graph::Node::Var(Var::X | Var::Y | Var::Z))) {
                    return Err(self.err(e, "placements cannot depend on x, y or z"));
                }
                Ok(f)
            }
            other => Err(self.err_hint(e, format!("{what} must be a frame, not a {}", other.kind()), Some("e.g. (at (move 0 0 5)) or (frame-of inst name)".into()))),
        }
    }

    /// `err` from evaluating `e`, written where a placement or a step of
    /// one goes: a call of std's (a shape's move or turn) or of no function
    /// that failed as written (`(move 0 0 5)`, `(rotate-z 90)`, `(rotz
    /// 90)`), with no hint of its own, says how one is written (`how`). A
    /// call of the model's own function keeps its error as it is.
    pub(super) fn misplaced(&self, mut err: InkError, e: &Sexp, how: &str) -> InkError {
        let Some(head) = e.list().and_then(|l| l.first()) else { return err };
        let call = head.sym().is_some_and(|h| match self.function(h) {
            Some(f) => &*f.file == "std",
            None => self.def(h).is_none() && !SPECIAL.contains(&h) && builtin(h).is_none() && !COMPARISONS.contains(&h),
        });
        if call && err.hint.is_none() && err.file == *self.file && [e, head].iter().any(|s| (s.line, s.col) == (err.line, err.col)) {
            err.hint = Some(how.into());
        }
        err
    }

    /// `(frame-of PATH [NAME])`: the frame NAME declared by what PATH ends
    /// at, placed; NAME omitted, the placement itself. In an assembly's
    /// clauses a path starting at one of its instances is in its
    /// coordinates; otherwise PATH is a body (its coordinates) or a path
    /// from a top-level assembly (world coordinates).
    pub(super) fn frame_of(&mut self, path: &str, name: Option<&str>) -> Result<Xf, String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or_default();
        let mut xf = self.identity();
        let bound = self.files[self.scope.file as usize].objects.get(first);
        let mut child = bound.cloned().unwrap_or_else(|| first.to_owned());
        let mut on = false;
        let mut visible = bound.is_some();
        // Old files could read the caller's instance here. Both readings
        // are meaningful, and choosing the file's frame can move a part.
        if bound.is_some()
            && let Some(b) = &self.building
            && self.local_assembly.as_deref() != Some(b.asm.name.as_str())
            && (b.asm.instances.iter().any(|i| i.name == first) || b.later.iter().any(|n| n == first))
        {
            return Err(format!(
                "ambiguous frame: '{first}' names both '{child}' in this file and an instance of '{}'. Pass the frame as an argument: read (frame-of {path}{}) in the assembly for the instance, or in a top-level def for the file's frame, or rename instance '{first}'",
                b.asm.name,
                name.map_or(String::new(), |n| format!(" {n}"))
            ));
        }
        if let Some(owner) = self.local_assembly.as_deref() {
            let building = self.building.as_ref().filter(|b| b.asm.name == owner);
            let asm = building.map(|b| &b.asm).or_else(|| self.assembly(owner)).ok_or_else(|| format!("the local frames of '{owner}' are still being built"))?;
            if let Some(i) = asm.instances.iter().find(|i| i.name == first) {
                if self.own_frame {
                    return Err(format!(
                        "a relation end's FRAME is in its part's own coordinates, but (frame-of {path} ...) is in assembly '{}''s. Name a frame the part declares, (PART NAME), or write one, (PART (at ...))",
                        asm.name
                    ));
                }
                (xf, child) = (i.at, i.child.clone());
                on = building.is_some();
                visible = true;
            } else if building.is_some_and(|b| b.later.iter().any(|l| l == first)) {
                return Err(format!("'{first}' is placed later in '{owner}': move its clause above this one"));
            }
        }
        if !visible || (self.body(&child).is_none() && self.assembly(&child).is_none()) {
            let mut names: Vec<&str> = self.bodies.iter().map(|b| b.name.as_str()).chain(self.assemblies.iter().map(|a| a.name.as_str())).collect();
            if let Some(b) = &self.building {
                names.extend(b.asm.instances.iter().map(|i| i.name.as_str()));
            }
            let hint = closest(first, names.into_iter()).unwrap_or_default();
            return Err(format!("no instance, body or assembly named '{first}'. {hint}").trim().to_string());
        }
        for seg in segs {
            let Some(sub) = self.assembly(&child) else { return Err(format!("'{child}' is a body; '{path}' goes through it")) };
            let Some(i) = sub.instances.iter().find(|i| i.name == seg) else {
                let hint = choices(seg, "Its instances", sub.instances.iter().map(|i| i.name.as_str()));
                return Err(format!("'{child}' has no instance '{seg}'. {hint}"));
            };
            let (at, c) = (i.at, i.child.clone());
            xf = self.compose(&at, &xf);
            child = c;
        }
        // What the instance being built is placed on: the instance under
        // PATH the frame is of (an assembly's frame of one of its own
        // instances: a servo's output frame is its output's).
        let under: Vec<String> = name.and_then(|n| self.assembly(&child)?.frame_on.iter().find(|(f, ..)| f == n).map(|(_, on, _)| on.clone())).unwrap_or_default();
        if on {
            let read = &mut self.building.as_mut().unwrap().read;
            match under.is_empty() {
                true => read.push(path.to_string()),
                false => read.extend(under.iter().map(|u| format!("{path}/{u}"))),
            }
        }
        let Some(name) = name else { return Ok(xf) };
        let frames = self.frames_of(&child);
        match frames.iter().find(|(n, _)| n == name) {
            Some((_, f)) => Ok(self.compose(f, &xf)),
            None => Err(Self::no_frame(&child, name, &frames)),
        }
    }

    /// The frames a body or an assembly declares.
    pub(super) fn frames_of(&self, child: &str) -> Vec<(String, Xf)> {
        self.body(child).map(|b| b.frames.clone()).or(self.assembly(child).map(|a| a.frames.clone())).unwrap_or_default()
    }

    pub(super) fn no_frame(child: &str, name: &str, frames: &[(String, Xf)]) -> String {
        if frames.is_empty() {
            return format!("'{child}' has no frame '{name}' (it declares none). Add (frame {name} F) to it, or place it by WHERE alone");
        }
        let names: Vec<&str> = frames.iter().map(|(n, _)| n.as_str()).collect();
        format!("'{child}' has no frame '{name}'. Its frames: {}", names.join(", "))
    }

    /// The inverse of a rigid placement: (R^T, -R^T t).
    pub fn inverse(&mut self, a: &Xf) -> Xf {
        let mut r = a.r;
        let mut t = a.t;
        for i in 0..3 {
            for j in 0..3 {
                r[i][j] = a.r[j][i];
            }
        }
        for i in 0..3 {
            let mut acc = self.g.constant(0.0);
            for k in 0..3 {
                let p = self.g.frame_bin(BinOp::Mul, r[i][k], a.t[k]);
                acc = self.g.frame_bin(BinOp::Add, acc, p);
            }
            t[i] = self.g.un(UnOp::Neg, acc);
        }
        Xf { r, t }
    }

    /// Contact ends are instance names, or paths into sub-assemblies
    /// (`inst/sub-inst/...`), relative to the assembly being defined.
    pub(super) fn check_contact_path(&self, asm: &Assembly, path: &str) -> Result<(), String> {
        let mut segs = path.split('/');
        let first = segs.next().unwrap_or("");
        let Some(inst) = asm.instances.iter().find(|i| i.name == first) else {
            let hint = closest(first, asm.instances.iter().map(|i| i.name.as_str())).unwrap_or_default();
            return Err(format!("no instance '{first}' (read so far). {hint}").trim().to_string());
        };
        let mut child = inst.child.clone();
        for seg in segs {
            let sub = self.assembly(&child).ok_or_else(|| format!("'{child}' is a body; '{path}' goes too deep"))?;
            let i = sub
                .instances
                .iter()
                .find(|i| i.name == seg)
                .ok_or_else(|| format!("'{child}' has no instance '{seg}'. {}", choices(seg, "Its instances", sub.instances.iter().map(|i| i.name.as_str()))))?;
            child = i.child.clone();
        }
        Ok(())
    }

    /// Every relation of what `target` places (Doc::target): those
    /// declared in it and under it, and those an enclosing assembly
    /// declares between parts of it. Its ends are full paths (matching
    /// `Item::label`) and its zone in world coordinates: A's part placed
    /// where A is, B's where B is.
    pub fn relations_under(&mut self, target: &str) -> Vec<PlacedRelation> {
        match self.target(target) {
            Some(t) => self.relations_of(&t),
            None => Vec::new(),
        }
    }

    /// The relations of a resolved target (Doc::relations_under).
    pub fn relations_of(&mut self, t: &Target) -> Vec<PlacedRelation> {
        let mut out = Vec::new();
        // An end is the target's where it is at or under it, or above it
        // (an enclosing instance holds all of it).
        let within = |end: &str| end == t.path || end.starts_with(&format!("{}/", t.path)) || t.path.starts_with(&format!("{end}/"));
        for (path, asm, xf) in &t.above {
            let mut declared = Vec::new();
            self.relations_here(asm, path, xf, &mut declared);
            out.extend(declared.into_iter().filter(|r| within(&r.a) && within(&r.b)));
        }
        self.collect_relations(&t.def, &t.path, &t.at, &mut out);
        out
    }
    /// Every clause `package` owns written in an assembly placed under
    /// `target` (a top-level assembly, or a path into one), each with the
    /// path of its assembly's instance and where that is placed in the world.
    pub fn clauses_under(&mut self, target: &str, package: &str) -> Vec<PlacedClause> {
        let mut out = Vec::new();
        let Some((xf, asm_name)) = self.placement(target) else { return out };
        self.collect_clauses(&asm_name, target, &xf, package, &mut out);
        out
    }

    fn collect_clauses(&mut self, asm: &str, prefix: &str, xf: &Xf, package: &str, out: &mut Vec<PlacedClause>) {
        let Some(a) = self.assembly(asm) else { return };
        out.extend(a.clauses.iter().filter(|c| c.package == package).map(|c| PlacedClause { path: prefix.to_string(), at: *xf, clause: c.clone() }));
        let insts: Vec<(String, String, Xf)> = a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at)).collect();
        for (name, child, at) in insts {
            let xf = self.compose(&at, xf);
            self.collect_clauses(&child, &format!("{prefix}/{name}"), &xf, package, out);
        }
    }

    /// Where an instance path (`asm/inst/...`) is placed in its top
    /// assembly's coordinates, its instances composed from the top (as
    /// `items` places what is under it), and what it places. A top-level
    /// name: the identity and itself. None: no such path.
    pub fn placement(&mut self, path: &str) -> Option<(Xf, String)> {
        self.target(path).map(|t| (t.at, t.def))
    }

    /// The body or assembly `path` names (a definition name, or an
    /// instance path from a top-level assembly), where it is placed, the
    /// file that defines it and the assemblies enclosing it. None: no
    /// such body, assembly or path.
    pub fn target(&mut self, path: &str) -> Option<Target> {
        let mut segs = path.split('/');
        let mut def = segs.next().unwrap_or_default().to_string();
        let (mut at, mut above, mut here) = (self.identity(), Vec::new(), def.clone());
        for seg in segs {
            let (i, next) = self.assembly(&def).and_then(|a| a.instances.iter().find(|i| i.name == seg)).map(|i| (i.at, i.child.clone()))?;
            above.push((here.clone(), def, at));
            at = self.compose(&i, &at);
            here = format!("{here}/{seg}");
            def = next;
        }
        let file = self.assembly(&def).map(|a| a.file.clone()).or_else(|| self.body(&def).map(|b| b.file.clone()))?;
        Some(Target { path: path.to_string(), def, file, at, above })
    }

    fn collect_relations(&mut self, asm: &str, prefix: &str, xf: &Xf, out: &mut Vec<PlacedRelation>) {
        self.relations_here(asm, prefix, xf, out);
        let Some(a) = self.assembly(asm) else { return };
        let insts: Vec<(String, String, Xf)> = a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at)).collect();
        for (name, child, at) in insts {
            let xf = self.compose(&at, xf);
            self.collect_relations(&child, &format!("{prefix}/{name}"), &xf, out);
        }
    }

    /// The relations assembly `asm` itself declares, placed at `xf` under
    /// path `prefix`.
    fn relations_here(&mut self, asm: &str, prefix: &str, xf: &Xf, out: &mut Vec<PlacedRelation>) {
        let Some(a) = self.assembly(asm) else { return };
        let relations: Vec<(String, String, Option<Zone>, Option<Rc<Kind>>)> = a.relations.iter().map(|r| (r.a.clone(), r.b.clone(), r.zone, r.kind.clone())).collect();
        for (pa, pb, zone, kind) in relations {
            // A field in an end's coordinates, placed: its path's
            // instances, in order.
            let placed = |doc: &mut Self, path: &str, f: NodeId| {
                let (mut at, mut child) = (*xf, asm.to_string());
                for seg in path.split('/') {
                    let i = doc.assembly(&child).and_then(|a| a.instances.iter().find(|i| i.name == seg)).map(|i| (i.at, i.child.clone())).unwrap();
                    at = doc.compose(&i.0, &at);
                    child = i.1;
                }
                let Value::Field(z, _) = doc.to_world(&Value::Field(f, None), &at) else { unreachable!() };
                z
            };
            let placed_zone = zone.map(|z| {
                let w = placed(self, &pa, z.a);
                match z.b {
                    Some(b) => {
                        let wb = placed(self, &pb, b);
                        self.g.max(w, wb)
                    }
                    None => w,
                }
            });
            out.push(PlacedRelation { a: format!("{prefix}/{pa}"), b: format!("{prefix}/{pb}"), zone: placed_zone, local: zone, kind, at: *xf });
        }
    }
}
