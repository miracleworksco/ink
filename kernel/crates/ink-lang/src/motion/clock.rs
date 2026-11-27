use super::*;

pub(super) fn required_drives(doc: &crate::Doc, roots: Vec<(u32, [interval::I; 3])>) -> Vec<(DriveRead, [interval::I; 3], bool)> {
    let mut todo: Vec<_> = roots.into_iter().map(|(f, b)| (f, b, false)).collect();
    let mut providers = std::collections::HashSet::new();
    let mut seen = std::collections::HashSet::new();
    let mut drives = Vec::new();
    while let Some((f, b, integrated)) = todo.pop() {
        for &read in doc.g.motion_reads(f) {
            match read {
                MotionRead::Param(k) => {
                    if let Some(s) = doc.solves.iter().find(|s| s.owns.contains(&(k as usize)))
                        && providers.insert((s.ext, integrated))
                    {
                        let form = &doc.ext[s.ext];
                        todo.extend(form.form.integration_inputs().iter().map(|&f| (f, [interval::I::ENTIRE; 3], true)));
                        todo.extend(doc.form_motion_inputs(form).into_iter().map(|f| (f, [interval::I::ENTIRE; 3], integrated)));
                    }
                }
                MotionRead::Drive(d) if seen.insert((d, b.map(|i| (i.lo.to_bits(), i.hi.to_bits())), integrated)) => drives.push((d, b, integrated)),
                MotionRead::Drive(_) => {}
            }
        }
    }
    drives
}

/// Load-time consumers that replace a field with numbers or samples must
/// retain its domain, or decline to materialize it.
pub(crate) fn needs_speed(doc: &crate::Doc, f: u32) -> bool {
    doc.g.has_speed_reads() && !doc.g.motion_reads(f).is_empty() && !required_drives(doc, vec![(f, [interval::I::ENTIRE; 3])]).is_empty()
}

/// The derivative at an endpoint needs its two-sided neighbourhood. In
/// particular, a clock fixed exactly at a switch cannot see just one side.
fn neighbourhood((a, b): (f64, f64)) -> (f64, f64) {
    (a.next_down(), b.next_up())
}

pub(super) fn drive_jump(v: &mut ModelView, read: DriveRead, b: [interval::I; 3], span: (f64, f64)) -> Option<((f64, f64), Evidence, bool)> {
    let (drive, at) = (read.field, read.at);
    let check = LocalDrive::new(&v.doc, read);
    let moving: Vec<_> = v.doc.solves_read(&[drive]).into_iter().filter(|&i| v.doc.solves[i].on_t).collect();
    let roots: Vec<_> = std::iter::once(drive).chain(at.into_iter().flatten()).collect();
    // A refused input is unavailable, not a jump: say so as its refusal.
    if let Err(e) = v.want(&roots) {
        return Some((span, Evidence::Unsolved(e), false));
    }
    if let Some(e) = v.refusal(&roots) {
        return Some((span, Evidence::Unsolved(e), false));
    }
    // At its own clock, a simulated coordinate has the position and
    // speed shared by consecutive integration steps. Composing it with a
    // globally Lipschitz unary operation (e.g. a turn's sine) preserves
    // that property. Commands are checked separately in required_drives.
    if !moving.is_empty() && at.is_none() && !read.near_params && !read.near_space.contains(&true) && v.ctx.p_spans.iter().next().is_none() {
        // The speed of an integrated coordinate as a claim writes it,
        // (diff q t), is a marked node, not the form's own: match it
        // exactly, and nothing else that carries a speed mark.
        let coordinates: Vec<(u32, usize)> = moving.iter().filter_map(|&i| v.doc.ext[v.doc.solves[i].ext].form.integrated_coordinate().map(|q| (q, i))).collect();
        let speeds: Vec<(u32, u32, usize)> = coordinates.into_iter().map(|(q, i)| (v.doc.derivative(q, Var::T), q, i)).collect();
        let mut f = drive;
        loop {
            if let Some(r) = moving.iter().find_map(|&i| v.doc.ext[v.doc.solves[i].ext].form.integrated_on(f, span)) {
                return r.err().map(|why| (span, Evidence::NotCovered(why), false));
            }
            if let Some(&(_, q, i)) = speeds.iter().find(|s| s.0 == f)
                && let Some(r) = v.doc.ext[v.doc.solves[i].ext].form.integrated_on(q, span)
            {
                return r.err().map(|why| (span, Evidence::NotCovered(why), false));
            }
            match v.doc.g.nodes[f as usize] {
                Node::Tag(_, a) | Node::Real(a) | Node::Un(UnOp::Neg | UnOp::Abs | UnOp::Rad | UnOp::Deg | UnOp::Sin | UnOp::Cos, a) => f = a,
                _ => break,
            }
        }
    }
    let solved = v.doc.solves_read(&roots).into_iter().any(|i| v.doc.solves[i].on_t);
    if read.near_params && solved {
        return Some((span, Evidence::NotCovered("a mixed derivative needs the solve's path over a parameter neighbourhood as well as time".into()), false));
    }
    let scope = if solved {
        match scope_over(v, &roots, neighbourhood(span)) {
            Ok(scope) => Some(scope),
            Err(why) => return Some((span, Evidence::NotCovered(format!("the drive's rate cannot be bounded: {why}")), false)),
        }
    } else {
        None
    };
    let s = may_jump(span, |piece| {
        let s = neighbourhood(piece);
        let mut ctx = scope.as_ref().map_or_else(|| v.ctx.clone(), |(c, k)| c.with_param(*k, s));
        ctx.t = ink_core::eval::span_mid(s);
        ctx.t_span = Some(s);
        // The local clock and coordinates may have been warped. Provider
        // values and rates are always enclosed on the model clock.
        ctx.p_spans = ink_core::eval::ParamSpans::new(ctx.p_spans.iter().map(|p| (p.k, p.span)));
        check.bounded(ctx, b)
    })?;
    Some((
        s,
        Evidence::NotCovered(format!("it reads the speed of a drive whose rate has no bound for t in [{}, {}]. It may jump there, and a jump is unbounded speed", s.0, s.1)),
        true,
    ))
}

/// A derivative before any warp, and the map at which it is evaluated.
struct LocalDrive {
    tape: Tape,
    coords: Option<[Tape; 4]>,
    near_space: [bool; 3],
    near_params: Vec<usize>,
}

impl LocalDrive {
    fn new(doc: &crate::Doc, read: DriveRead) -> Self {
        Self {
            tape: Tape::new(&doc.g, read.field),
            coords: read.at.map(|at| at.map(|n| Tape::new(&doc.g, n))),
            near_space: read.near_space,
            near_params: if read.near_params { doc.params_read(&[read.field]) } else { Vec::new() },
        }
    }

    fn bounded(&self, mut ctx: Ctx, mut b: [interval::I; 3]) -> bool {
        if let Some(coords) = &self.coords {
            let p = coords.each_ref().map(|t| ink_core::eval::eval_interval_scope(t, &ctx, b, None));
            b = [p[0], p[1], p[2]];
            let t = neighbourhood((p[3].lo, p[3].hi));
            if !t.0.is_finite() || !t.1.is_finite() {
                return false;
            }
            ctx.t = ink_core::eval::span_mid(t);
            ctx.t_span = Some(t);
        }
        for (interval, near) in b.iter_mut().zip(self.near_space) {
            if near {
                *interval = interval::I::new(interval.lo.next_down(), interval.hi.next_up());
            }
        }
        for &k in &self.near_params {
            let point = ctx.params[k];
            let band = ctx.p_spans.get(k as u32).unwrap_or((point, point));
            ctx.p_spans = ctx.p_spans.with_band(k as u32, neighbourhood(band));
        }
        ink_kernel::motion::rate(&self.tape, &ctx, b, ctx.t_span.expect("time scope")).is_some()
    }
}

/// A scalar load on a package's own clock. Its model parameters are held
/// while the package advances local time; model-clock providers are
/// checked separately by `jumps`. Do this before using a cached solution.
pub fn check_local_clock(doc: &crate::Doc, ctx: &Ctx, f: u32, span: (f64, f64)) -> Result<(), InkError> {
    for &read in doc.g.motion_reads(f) {
        if let MotionRead::Drive(read) = read {
            // Local time advances with the supplied parameter values held.
            // A derivative through a model-time provider also requires its
            // path, which this local interval does not enclose.
            let inputs: Vec<_> = std::iter::once(read.field).chain(read.at.into_iter().flatten()).collect();
            if doc.solves_read(&inputs).iter().any(|&i| doc.solves[i].on_t) {
                return Err(InkError::msg("the load's speed reads a moving model parameter on its own clock; its rate is not bounded over both clocks"));
            }
            let drive = LocalDrive::new(doc, read);
            if let Some((a, b)) = may_jump(span, |piece| {
                let s = neighbourhood(piece);
                let ctx = Ctx { t: ink_core::eval::span_mid(s), t_span: Some(s), ..ctx.clone() };
                drive.bounded(ctx, [interval::I::point(0.0); 3])
            }) {
                return Err(InkError::msg(format!("the load reads a speed whose drive may jump on its own clock in [{a}, {b}]")));
            }
        }
    }
    Ok(())
}

/// Whether a shape's own field (not its placement) reads t: a placed
/// body's, or any part of a union or intersection of them.
pub(super) fn own_field_moves(v: &ModelView, f: u32) -> bool {
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
