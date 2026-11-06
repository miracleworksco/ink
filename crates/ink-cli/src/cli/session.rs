//! The loaded model a command works on (ink::view::ModelView), and what
//! commands need of it beyond the library: targets, values, colors.

use crate::*;
use ink::lang::command::{claim_scope, nominal, param_args};

/// The model a command works on (`ink::view`), loaded from the command's
/// file.
pub type Session = ModelView;

/// Commands that prove claims over ranges of params (tolerances, and each
/// (ranged lo hi), which no --param sets there: session); the others show,
/// measure or solve at one value of each param.
fn proves_over_ranges(cmd: &str) -> bool {
    matches!(cmd, "assert" | "interfere")
}

pub(crate) fn session(reg: &'static Registry, a: &Args) -> Result<Session, InkError> {
    let file = a.file.clone().ok_or_else(|| InkError::msg(format!("'{}' needs a file: ink {} <file.ink>", a.cmd, a.cmd)))?;
    let params = param_args(a)?;
    let nominal = nominal(a)?;
    if !proves_over_ranges(&a.cmd) {
        if let Some((k, _)) = params.ranges.first() {
            return Err(InkError::msg(format!("--param {k}=lo:hi: '{}' works at one value of each param. Ranges are for assert and interfere", a.cmd)));
        }
        if nominal {
            return Err(InkError::msg(format!("--nominal: '{}' already uses the params' values. It is for assert and interfere", a.cmd)));
        }
    }
    let mut s = ModelView::load_in(reg, Path::new(&file), &params.values, a.moment()?)?;
    if proves_over_ranges(&a.cmd)
        && let Some((k, (lo, hi))) = params.values.iter().find_map(|(k, _)| s.doc.ranged_named(k).map(|r| (k, r.1)))
    {
        return Err(InkError {
            hint: Some(format!("look at one value with ink eval, ink check (its conditions), ink render --param {k}=V")),
            ..InkError::msg(format!("--param {k}: '{k}' is (ranged {lo} {hi}), a range of the question. '{}' answers each claim over all of it", a.cmd))
        });
    }
    if proves_over_ranges(&a.cmd) {
        let scope = claim_scope(&s, a)?;
        s.set_claim_scope(&scope)?;
    }
    // assert and interfere refuse what rests on a solve that does not
    // solve at this moment, one by one (ModelView::refusal), and check
    // warns of it; the other commands need every solve the shapes rest on
    // (Doc::needed: rates no shape reads are read by claims and eval's
    // derivatives alone, each refusing or giving null where one is not
    // defined). A package's command asks itself (command::solved).
    if !proves_over_ranges(&a.cmd) && a.cmd != "check" && reg.command(&a.cmd).is_none() {
        ink::lang::command::solved(&s)?;
    }
    Ok(s)
}

/// The scope a claim command covered, for its report: each ranged param's
/// span, or that declared tolerances were dropped (--nominal).
pub(crate) fn scope_fields(out: &mut J, s: &Session) {
    let spans = |derived: bool| -> serde_json::Map<String, J> {
        s.ctx
            .p_spans
            .iter()
            .filter(|p| p.derived == derived && !ink::lang::unnamed_range(&s.doc.params[p.k as usize].name))
            .map(|p| (s.doc.params[p.k as usize].name.clone(), bounds(p.span.0, p.span.1)))
            .collect()
    };
    let (spans, solved) = (spans(false), spans(true));
    if !solved.is_empty() {
        out["solved_over_ranges"] = J::Object(solved);
    }
    let dropped = s.doc.params.iter().enumerate().any(|(k, p)| p.tol.is_some() && s.ctx.p_spans.get(k as u32).is_none());
    if !spans.is_empty() || dropped {
        out["param_spans"] = J::Object(spans);
        if dropped {
            out["tolerances_dropped"] = json!(true);
        }
    }
}

/// What commands need of a model beyond the library's view: targets named
/// on the command line, and presentation.
pub trait CliView {
    fn items(&mut self, a: &Args) -> Result<Vec<Item>, InkError>;
    /// The file as named on the command line.
    fn file(&self) -> String;
    fn stem(&self) -> String;
    fn value_at(&self, v: &Value, p: V3) -> J;
    fn describe(&self, v: &Value) -> J;
    fn frame_json(&self, xf: &Xf) -> J;
    fn kind_json(&self, k: &ink::lang::Kind) -> J;
    fn colors(&self, item: &Item, verts: &[V3], k: usize) -> Vec<[f32; 3]>;
}

/// A published mass's density not measured here (Doc::unmeasured), by
/// the params it is derived as: it has no value to show.
fn unmeasured(doc: &ink::lang::Doc, f: ink::graph::NodeId) -> Option<J> {
    let ks = doc.unmeasured(&[f]);
    (!ks.is_empty()).then(|| json!({"derived": ks.iter().map(|&k| &doc.params[k].name).collect::<Vec<_>>()}))
}

impl CliView for Session {
    fn items(&mut self, a: &Args) -> Result<Vec<Item>, InkError> {
        ink::lang::command::items(self, a)
    }

    fn file(&self) -> String {
        self.path.as_ref().map_or(String::new(), |p| p.to_string_lossy().into_owned())
    }

    fn stem(&self) -> String {
        Path::new(&self.file()).file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or("model".into())
    }

    fn value_at(&self, v: &Value, p: V3) -> J {
        match v {
            Value::Field(f) => match self.doc.g.as_const(*f) {
                Some(c) => n(c),
                None if let Some(j) = unmeasured(&self.doc, *f) => j,
                None => n(ink::eval::eval_point(&Tape::new(&self.doc.g, *f), &self.ctx, p)),
            },
            Value::Str(s) => json!(s),
            Value::Vec(items) => J::Array(items.iter().map(|i| self.value_at(i, p)).collect()),
            Value::Data(d) => d.value().clone(),
            Value::Fn(_) => json!("function"),
            Value::Frame(xf) => self.frame_json(xf),
            Value::Point(at) => json!({"point": at.map(|e| n(self.doc.scalar(e, &self.ctx)))}),
            Value::Relation(r) => json!({"relation": r.kind}),
        }
    }

    /// A frame at the current params and t: its origin and unit axes.
    fn frame_json(&self, xf: &Xf) -> J {
        let v = |e: ink::graph::NodeId| n(self.doc.scalar(e, &self.ctx));
        let axis = |j: usize| json!([v(xf.r[0][j]), v(xf.r[1][j]), v(xf.r[2][j])]);
        json!({"origin": xf.t.map(v), "x": axis(0), "y": axis(1), "z": axis(2)})
    }

    fn describe(&self, v: &Value) -> J {
        match v {
            Value::Field(f) => match self.doc.g.as_const(*f) {
                Some(c) => n(c),
                None if let Some(j) = unmeasured(&self.doc, *f) => j,
                // One resting on a solve that fails here has none.
                None if self.refusal(&[*f]).is_some() => J::Null,
                // Expressions of params (and t) have a current value.
                None if !self.doc.g.depends_on(*f, |x| matches!(x, ink::graph::Var::X | ink::graph::Var::Y | ink::graph::Var::Z)) => n(self.doc.scalar(*f, &self.ctx)),
                None => json!({"field": {"nodes": self.doc.g.reachable(&[*f]).len()}}),
            },
            Value::Str(s) => json!(s),
            Value::Vec(items) => J::Array(items.iter().map(|i| self.describe(i)).collect()),
            Value::Data(d) => json!({"data": d.value()}),
            Value::Fn(_) => json!("function"),
            Value::Frame(xf) => json!({"frame": self.frame_json(xf)}),
            Value::Point(at) => json!({"point": at.map(|e| n(self.doc.scalar(e, &self.ctx)))}),
            Value::Relation(r) => json!({"relation": r.kind}),
        }
    }

    /// A relation's kind (`relate`): its name, its data at the current
    /// params and t, what its zone leaves unchecked and its claims.
    fn kind_json(&self, k: &ink::lang::Kind) -> J {
        let data: serde_json::Map<String, J> = k.data.iter().map(|(n, v)| (n.clone(), self.describe(v))).collect();
        let mut j = json!({"kind": k.name, "data": data, "claims": k.claims});
        if let Some(what) = &k.skipped {
            j["skipped"] = json!(what);
        }
        // The heads of its packages' clauses, with their packages.
        if !k.clauses.is_empty() {
            j["clauses"] = k.clauses.iter().map(|c| json!({"head": c.head, "package": c.package})).collect();
        }
        j
    }

    /// Per-vertex colors from the body's look (ink::lang::look), the
    /// palette's k-th where it gives none.
    fn colors(&self, item: &Item, verts: &[V3], k: usize) -> Vec<[f32; 3]> {
        use ink::lang::look::{Look, PALETTE};
        let rgb = Look::of(&item.attrs).color.map_or([None; 3], |c| c.map(Some));
        paint(self, rgb, PALETTE[k % PALETTE.len()], verts)
    }
}

/// Each paint channel's field at each vertex, clamped to 0..1, `default`'s
/// where a channel has none.
pub fn paint<const N: usize>(s: &Session, fields: [Option<ink::graph::NodeId>; N], default: [f32; N], verts: &[V3]) -> Vec<[f32; N]> {
    let mut out = vec![default; verts.len()];
    for (j, f) in fields.into_iter().enumerate() {
        let Some(f) = f else { continue };
        let set = |o: &mut [f32; N], v: f64| o[j] = v.clamp(0.0, 1.0) as f32;
        match s.doc.g.as_const(f) {
            Some(c) => out.iter_mut().for_each(|o| set(o, c)),
            None => out.iter_mut().zip(eval_points_par(&Tape::new(&s.doc.g, f), &s.ctx, verts)).for_each(|(o, v)| set(o, v)),
        }
    }
    out
}

pub fn union_bounds(bs: &[(V3, V3, f64)]) -> Option<(V3, V3)> {
    bs.iter().map(|&(lo, hi, _)| (lo, hi)).reduce(ink::v3::box_union)
}
