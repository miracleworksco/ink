//! Results as fields: `(study-field form result)` reads any package form's
//! result (a study's stress, say) as a sampled field, through the
//! `Form` trait; and the helper packages build such results with.

use super::*;

impl Doc {
    /// `(study-field form result [(mode k)])`.
    pub(super) fn eval_result(&mut self, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        let usage = "usage: (study-field form result [(mode k)]): a package form's result as a field (a study's stress, a flow's speed), \
                     the results its package names; (mode k): of its mode k, where the form gives modes";
        let (st, q, opts) = match args {
            [st, q, opts @ ..] if opts.len() <= 1 => (st, q, opts),
            _ => return Err(self.err_hint(e, "bad study-field", Some(usage.into()))),
        };
        let (st, q) = (self.name_of(st)?, self.name_of(q)?);
        // (mode k): a whole number from 1, known at load time.
        let mode = match opts.first() {
            None => None,
            Some(o) => match o.list() {
                Some([h, k]) if h.sym() == Some("mode") => {
                    let k = self.const_num(k, env)?;
                    if !(k >= 1.0 && k == k.floor() && k < 1e9) {
                        return Err(self.err_hint(o, format!("(mode {k}): a mode is a whole number from 1"), Some(usage.into())));
                    }
                    Some(k as usize)
                }
                _ => return Err(self.err_hint(o, "bad option", Some(usage.into()))),
            },
        };
        // A solved result, sampled: what reads it is measured.
        let st = self.resolve(&st);
        let Some(x) = self.ext.iter().find(|x| x.name == st && !x.form.results().is_empty()) else {
            return Err(self.err(e, format!("no form with results named '{st}'. Define the study or package form before using it")));
        };
        let (head, form, inputs) = (x.head, x.form.clone(), self.form_inputs(x));
        let ctx = self.load_ctx_for(&inputs);
        let results: Vec<&str> = form.results().iter().map(|r| r.name).collect();
        if !results.contains(&q.as_str()) {
            let hint = format!("{head} '{st}' gives {}", results.join(", "));
            return Err(self.err_hint(e, format!("'{q}' is not available for {head} '{st}'"), Some(hint)));
        }
        // Solved at the params' values: every param it read is now
        // one number in the model.
        self.bake(&inputs);
        let out = form.run(self, &ctx, &RunOpts::default())?;
        if let Some(why) = out.failure() {
            return Err(self.err(e, format!("{head} '{st}' failed, so it has no '{q}': {why}")));
        }
        // Its mode k's outcome, and what the mode is.
        let (out, of_mode) = match mode {
            None => (out, String::new()),
            Some(k) => {
                let (m, what) = out.mode(k).map_err(|why| self.err(&opts[0], format!("{head} '{st}' has no mode {k}: {why}")))?;
                (m, format!(", {what}"))
            }
        };
        let data = out.field(&q).ok_or_else(|| self.err(e, format!("{head} '{st}' gave no '{q}'")))?;
        // A quantity holds its edge value past its grid, a shape's
        // field grows by the distance (ResultSpec::shape); its unit,
        // and its least and most value where it has one.
        let specs = form.results();
        let spec = specs.iter().find(|r| r.name == q).expect("listed above");
        let (unit, data) = (spec.unit, if spec.shape { data } else { data.holding() });
        // A vector's component (displacement-x): the vector, and the
        // component's axis.
        let component = specs.iter().find_map(|r| Some((format!("{head} {st}: {}{of_mode}", r.name), r.vector?.iter().position(|c| *c == q)?)));
        let span = out.samples(&q).and_then(|s| {
            let v = s.values.iter().filter(|v| !v.is_nan()).map(|&v| v as f64);
            let (lo, hi) = v.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)));
            (lo <= hi).then_some((lo, hi))
        });
        // A sample of the package's solution, which is itself numerical.
        let source = crate::grid::Source {
            what: format!("{head} {st}: {q}{of_mode}"),
            at: format!("{}:{}", self.file, e.line),
            kind: "result",
            res: data.h(),
            band: None,
            measured: None,
            note: Some("a package's solved result, sampled: numerical, with no error bound".into()),
            unit,
            span,
            component,
            undecided: out.undecided(&q),
        };
        let at = [self.g.var(Var::X), self.g.var(Var::Y), self.g.var(Var::Z)];
        Ok(Value::Field(self.g.grid(std::sync::Arc::new(data.with_source(source)), at), None))
    }
}

/// A result sampled on a grid, `NaN` where it has no value (outside the
/// material), as a field's data: values are carried a few cells past where
/// they end (averaging known neighbours), so interpolation near surfaces
/// keeps its level, and `fill` beyond.
pub fn sampled_result(mut vals: Vec<f32>, dims: [usize; 3], origin: V3, h: f64, fill: f32) -> crate::grid::GridData {
    let at = |i: usize, j: usize, k: usize| crate::grid::index([i, j, k], dims);
    for _ in 0..3 {
        let prev = vals.clone();
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    if !prev[at(i, j, k)].is_nan() {
                        continue;
                    }
                    let mut acc = (0.0f32, 0);
                    let mut add = |ii: isize, jj: isize, kk: isize| {
                        if ii >= 0 && jj >= 0 && kk >= 0 && (ii as usize) < dims[0] && (jj as usize) < dims[1] && (kk as usize) < dims[2] {
                            let v = prev[at(ii as usize, jj as usize, kk as usize)];
                            if !v.is_nan() {
                                acc = (acc.0 + v, acc.1 + 1);
                            }
                        }
                    };
                    let (ii, jj, kk) = (i as isize, j as isize, k as isize);
                    for (a, b, c) in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)] {
                        add(ii + a, jj + b, kk + c);
                    }
                    if acc.1 > 0 {
                        vals[at(i, j, k)] = acc.0 / acc.1 as f32;
                    }
                }
            }
        }
    }
    for v in vals.iter_mut() {
        if v.is_nan() {
            *v = fill;
        }
    }
    crate::grid::GridData::new(origin, h, dims, vals)
}
