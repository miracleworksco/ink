//! Load-time programming: functions, lists, `for`, conditionals, JSON,
//! derivatives.

use super::*;

/// `e` with the symbol `var` replaced by `val`, and `$var` inside other
/// symbols by `text`.
fn substitute(e: &Sexp, var: &str, val: &Sexp, text: &str) -> Sexp {
    let pat = format!("${var}");
    match &e.k {
        SK::Sym(s) if s == var => Sexp { k: val.k.clone(), line: e.line, col: e.col },
        SK::Sym(s) if s.contains(&pat) => Sexp { k: SK::Sym(s.replace(&pat, text)), line: e.line, col: e.col },
        SK::Str(s) if s.contains(&pat) => Sexp { k: SK::Str(s.replace(&pat, text)), line: e.line, col: e.col },
        SK::List(items) => Sexp { k: SK::List(items.iter().map(|x| substitute(x, var, val, text)).collect()), line: e.line, col: e.col },
        _ => e.clone(),
    }
}

/// The parts of `(for var list form ...)`, if `it` is one.
pub(super) fn for_parts(it: &Sexp) -> Option<(&Sexp, &Sexp, &[Sexp])> {
    match it.list() {
        Some([head, var, list, body @ ..]) if head.sym() == Some("for") => Some((var, list, body)),
        _ => None,
    }
}

/// The forms a `for` makes for one element: `body` with `var` replaced by
/// `val`, and `$var` inside names by `text`.
pub(super) fn for_made(body: &[Sexp], var: &str, val: &Sexp, text: &str) -> Vec<Sexp> {
    body.iter().map(|b| substitute(b, var, val, text)).collect()
}

impl Doc {
    /// Expand `(for var list form ...)` in a sequence of clauses: the list
    /// is evaluated at load time, and for each element the forms are
    /// repeated with `var` replaced by the element and `$var` inside
    /// symbols replaced by its text (for names: `bolt-$i`). Works inside
    /// study, assembly, body and optimize clauses (any package form's);
    /// nests (a file's
    /// top-level forms are expanded as they are reached: Doc::run_forms).
    pub(super) fn expand_for(&mut self, items: &[Sexp]) -> R<Vec<Sexp>> {
        self.expand_for_in(items, &mut self.arg_env())
    }

    /// `expand_for` with the lists read in `env`: inside an expression (an
    /// `except`'s paths), where a function's arguments and `let` names are.
    pub(super) fn expand_for_in(&mut self, items: &[Sexp], env: &mut Env) -> R<Vec<Sexp>> {
        let mut out = Vec::with_capacity(items.len());
        for it in items {
            let Some((var, list, body)) = for_parts(it) else {
                out.push(it.clone());
                continue;
            };
            let (var, values) = self.for_values_in(var, list, env)?;
            for (val, text) in values {
                let made = for_made(body, &var, &val, &text);
                out.extend(self.expand_for_in(&made, env)?);
            }
        }
        Ok(out)
    }

    /// A `for`'s variable, and its list's elements as source text and
    /// display.
    pub(super) fn for_values(&mut self, var: &Sexp, list: &Sexp) -> R<(String, Vec<(Sexp, String)>)> {
        self.for_values_in(var, list, &mut self.arg_env())
    }

    fn for_values_in(&mut self, var: &Sexp, list: &Sexp, env: &mut Env) -> R<(String, Vec<(Sexp, String)>)> {
        let usage = "usage: (for var list form ...): repeats the forms for each element, with var and $var in names replaced";
        let var = var.sym().ok_or_else(|| self.err_hint(var, "for needs a variable name", Some(usage.into())))?.to_string();
        let values = match self.eval(list, env)? {
            Value::Vec(v) => v,
            v => vec![v],
        };
        let literals = values.iter().map(|v| self.literal(v, list)).collect::<R<Vec<_>>>()?;
        Ok((var, literals))
    }

    /// A load-time value as source text (for `for`), and its display.
    fn literal(&self, v: &Value, at: &Sexp) -> R<(Sexp, String)> {
        let mk = |k: SK| Sexp { k, line: at.line, col: at.col };
        Ok(match v {
            Value::Field(f, said) => {
                let x = match self.g.as_const(*f) {
                    Some(x) => x,
                    None => self.load_number(at, *f)?.ok_or_else(|| self.err(at, "for values must be known at load time"))?,
                };
                let text = if x.fract() == 0.0 && x.abs() < 1e15 { format!("{}", x as i64) } else { format!("{x}") };
                // The number keeps the unit its value was said in.
                (mk(SK::Num(x, said.as_ref().map(|s| crate::sexp::Written::From(Box::new((**s).clone()))))), text)
            }
            Value::Str(t) => (mk(SK::Str(t.clone())), t.clone()),
            Value::Vec(items) => {
                let mut parts = vec![mk(SK::Sym("vec".into()))];
                let mut texts = Vec::new();
                for i in items {
                    let (s, t) = self.literal(i, at)?;
                    parts.push(s);
                    texts.push(t);
                }
                (mk(SK::List(parts)), texts.join("-"))
            }
            _ => return Err(self.err(at, "for iterates over numbers, strings and vectors")),
        })
    }

    /// Balanced-tree reduction (shallow graphs; interval pruning works per branch).
    fn reduce_values(&mut self, f: &Callable, items: &[Value], at: &Sexp) -> R<Value> {
        if items.len() == 1 {
            return Ok(items[0].clone());
        }
        // Fast path: an associative kernel op over plain fields.
        if let Callable::Builtin(name) = f
            && let Some(Builtin::Variadic(op)) = builtin(name)
        {
            let ids: Option<Vec<NodeId>> = items.iter().map(|v| if let Value::Field(n, _) = v { Some(*n) } else { None }).collect();
            if let Some(ids) = ids {
                let said = self.sum_said(name, items, at)?;
                return Ok(Value::Field(self.reduce(op, &ids), said));
            }
        }
        let m = items.len() / 2;
        let a = self.reduce_values(f, &items[..m], at)?;
        let b = self.reduce_values(f, &items[m..], at)?;
        self.call(f, vec![a, b], at)
    }

    /// Load-time programming forms of `eval`.
    pub(super) fn eval_author(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        match name {
            "lambda" => {
                let [ps, body] = args else { return Err(self.err(e, "usage: (lambda (arg ...) body)")) };
                let (params, rest) = self.arg_names(ps, "lambda")?;
                return Ok(Value::Fn(Callable::Closure(Rc::new(Closure { params, rest, body: body.clone(), env: env.clone(), file: self.file.clone(), scope: self.scope }))));
            }
            "map" => {
                let usage = "usage: (map f list ...) applies f element-wise (lists of equal length)";
                let Some((f, lists)) = args.split_first() else { return Err(self.err(e, usage)) };
                let f = self.callable(f, env)?;
                let lists = lists.iter().map(|l| self.list(l, env)).collect::<R<Vec<_>>>()?;
                let Some(n) = lists.first().map(|l| l.len()) else { return Err(self.err(e, usage)) };
                if lists.iter().any(|l| l.len() != n) {
                    let lens: Vec<usize> = lists.iter().map(|l| l.len()).collect();
                    return Err(self.err(e, format!("map over lists of different lengths {lens:?}")));
                }
                let mut out = Vec::with_capacity(n);
                for i in 0..n {
                    out.push(self.call(&f, lists.iter().map(|l| l[i].clone()).collect(), e)?);
                }
                return Ok(Value::Vec(out));
            }
            "reduce" => {
                let [f, l] = args else { return Err(self.err(e, "usage: (reduce f list); f must be associative (reduced as a balanced tree)")) };
                let f = self.callable(f, env)?;
                let items = self.list(l, env)?;
                if items.is_empty() {
                    return Err(self.err(e, "reduce of an empty list"));
                }
                return self.reduce_values(&f, &items, e);
            }
            "fold" => {
                let [f, init, l] = args else { return Err(self.err(e, "usage: (fold f init list) left fold: (f (f init a) b) ...")) };
                let f = self.callable(f, env)?;
                let mut acc = self.eval(init, env)?;
                for item in self.list(l, env)? {
                    acc = self.call(&f, vec![acc, item], e)?;
                }
                return Ok(acc);
            }
            "flatten" => {
                let [l] = args else { return Err(self.err(e, "usage: (flatten list-of-lists)")) };
                let mut out = Vec::new();
                for item in self.list(l, env)? {
                    match item {
                        Value::Vec(v) => out.extend(v),
                        other => out.push(other),
                    }
                }
                return Ok(Value::Vec(out));
            }
            "apply" => {
                let [f, l] = args else { return Err(self.err(e, "usage: (apply f list)")) };
                let f = self.callable(f, env)?;
                let items = self.list(l, env)?;
                return self.call(&f, items, e);
            }
            "range" => {
                let nums = args.iter().map(|a| self.const_num(a, env)).collect::<R<Vec<_>>>()?;
                if nums.iter().any(|v| !v.is_finite()) {
                    return Err(self.err(e, format!("{name} needs finite arguments")));
                }
                let usage = "usage: (range n) (range a b) (range a b step), constant arguments";
                // Counts are whole numbers (none is rounded), and are
                // checked against the limit before anything is allocated.
                if let Some(&v) = match (name, nums.as_slice()) {
                    ("range", [_] | [_, _]) => nums.iter().find(|v| v.fract() != 0.0),
                    _ => None,
                } {
                    return Err(self.err_hint(e, format!("{name} counts in whole numbers, got {v}"), Some(usage.into())));
                }
                let (start, step, count) = match (name, nums.as_slice()) {
                    ("range", [n]) => (0.0, 1.0, n.max(0.0)),
                    ("range", [a, b]) => (*a, 1.0, (b - a).max(0.0)),
                    ("range", [a, b, st]) if *st != 0.0 => (*a, *st, ((b - a) / st).ceil().max(0.0)),
                    _ => return Err(self.err(e, usage)),
                };
                if !count.is_finite() || !step.is_finite() {
                    return Err(self.err(e, format!("{name} needs a finite count and spacing")));
                }
                if count > 1e6 {
                    return Err(self.err(e, format!("{name} would produce {} elements (limit 1e6)", if count < 1e12 { count.to_string() } else { format!("{count:e}") })));
                }
                let vals: Vec<f64> = (0..count as usize).map(|i| start + step * i as f64).collect();
                return Ok(Value::Vec(vals.into_iter().map(|v| Value::Field(self.g.constant(v), None)).collect()));
            }
            "len" => {
                let [l] = args else { return Err(self.err(e, "usage: (len list)")) };
                let n = self.list(l, env)?.len();
                return Ok(Value::Field(self.g.constant(n as f64), None));
            }
            "concat" => {
                // Lists end to end, or strings (names for (tag ...)).
                let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
                if !vals.is_empty() && vals.iter().all(|v| matches!(v, Value::Str(_))) {
                    return Ok(Value::Str(vals.iter().map(|v| if let Value::Str(s) = v { s.as_str() } else { "" }).collect()));
                }
                let mut out = Vec::new();
                for (v, a) in vals.into_iter().zip(args) {
                    match v {
                        Value::Vec(v) => out.extend(v),
                        other => return Err(self.err(a, format!("concat takes lists, or strings: got a {}", other.kind()))),
                    }
                }
                return Ok(Value::Vec(out));
            }
            "char-codes" => {
                // A string's characters as numbers (Unicode code points), so
                // load-time code can read text: index tables, decode data.
                let [a] = args else { return Err(self.err(e, "usage: (char-codes \"text\")")) };
                let Value::Str(s) = self.eval(a, env)? else { return Err(self.err(a, "char-codes takes a string")) };
                return Ok(Value::Vec(s.chars().map(|c| Value::Field(self.g.constant(c as u32 as f64), None)).collect()));
            }
            "slice" => {
                let usage = "usage: (slice list start [end]); negative indices count from the end";
                let (l, rest) = args.split_first().ok_or_else(|| self.err(e, usage))?;
                let items = self.list(l, env)?;
                let n = items.len() as i64;
                let idx = |d: &mut Doc, a: &Sexp, env: &mut Env| -> R<usize> {
                    let v = d.whole(a, env, "a slice index")? as i64;
                    Ok((if v < 0 { n + v } else { v }).clamp(0, n) as usize)
                };
                let (a, b) = match rest {
                    [a] => (idx(self, a, env)?, n as usize),
                    [a, b] => (idx(self, a, env)?, idx(self, b, env)?),
                    _ => return Err(self.err(e, usage)),
                };
                return Ok(Value::Vec(if a < b { items[a..b].to_vec() } else { Vec::new() }));
            }
            "if" => {
                let [c, a, b] = args else { return Err(self.err(e, "usage: (if constant-condition then else)")) };
                let cf = self.field(c, env)?;
                let Some(cv) = self.load_number(c, cf)? else {
                    return Err(self.err_hint(c, "if needs a condition known at load time", Some("for a condition that varies over space, use (select c a b)".into())));
                };
                return self.eval(if cv != 0.0 { a } else { b }, env);
            }
            "read-json" => {
                let [p] = args else { return Err(self.err(e, "usage: (read-json \"file.json\")")) };
                let SK::Str(rel) = &p.k else { return Err(self.err(p, "read-json path must be a string")) };
                let path = self.dir().join(rel);
                let key = path.canonicalize().unwrap_or(path.clone());
                let data = match self.json.get(&key).cloned() {
                    Some(d) => {
                        self.read_again(&path);
                        d
                    }
                    None => {
                        let text = std::fs::read_to_string(&path)
                            .inspect_err(|_| self.files_read.push((path.clone(), None)))
                            .map_err(|er| self.err(p, format!("cannot read {}: {er}", shown(&path))))?;
                        self.read(&shown(&path), Some(&path), text.as_bytes());
                        let json: serde_json::Value = serde_json::from_str(&text).map_err(|er| self.err(p, format!("invalid JSON in {}: {er}", shown(&path))))?;
                        let d = crate::Data::new(json);
                        self.json.insert(key, d.clone());
                        d
                    }
                };
                return self.json_value(&data).map_err(|m| self.err(p, m));
            }
            "get" => {
                let (d, k, default) = match args {
                    [d, k] => (d, k, None),
                    [d, k, x] => (d, k, Some(x)),
                    _ => return Err(self.err(e, "usage: (get data \"key\" [default])")),
                };
                let obj = match self.eval(d, env)? {
                    Value::Data(obj) => obj,
                    Value::Frame(xf) => {
                        let Value::Str(key) = self.eval(k, env)? else { return Err(self.err(k, "get key must be a string")) };
                        let v = match key.as_str() {
                            "origin" => xf.t,
                            "x" | "y" | "z" => {
                                let j = (key.as_bytes()[0] - b'x') as usize;
                                [xf.r[0][j], xf.r[1][j], xf.r[2][j]]
                            }
                            _ => return Err(self.err_hint(k, format!("a frame has no '{key}'"), Some("its keys: \"origin\", \"x\", \"y\", \"z\"".into()))),
                        };
                        return Ok(Value::Vec(v.map(|f| Value::Field(f, None)).to_vec()));
                    }
                    _ => return Err(self.err(d, "get expects data (a JSON object) or a frame")),
                };
                let Value::Str(key) = self.eval(k, env)? else { return Err(self.err(k, "get key must be a string")) };
                let Some(v) = obj.get(&key) else {
                    if let Some(x) = default {
                        return self.eval(x, env);
                    }
                    let keys: Vec<&str> = obj.value().as_object().map(|o| o.keys().map(String::as_str).collect()).unwrap_or_default();
                    return Err(self.err_hint(k, format!("no key '{key}'"), Some(format!("its keys: {}", keys.join(", ")))));
                };
                return self.json_value(&v).map_err(|m| self.err(e, m));
            }
            "diff" => {
                let usage = "usage: (diff f var): the derivative of a field with respect to x, y, z, t or a param";
                let [f, v] = args else { return Err(self.err_hint(e, "bad diff", Some(usage.into()))) };
                let f = self.field(f, env)?;
                let name = self.name_of(v)?;
                let var = match name.as_str() {
                    "x" => Var::X,
                    "y" => Var::Y,
                    "z" => Var::Z,
                    "t" => Var::T,
                    // A param by name, or a value that is a param itself (a (ranged lo hi)).
                    p => match self.param_index(p).map(|k| k as u32).or_else(|| match self.eval(v, env) {
                        Ok(Value::Field(g, _)) => match self.g.nodes[g as usize] {
                            crate::graph::Node::Var(Var::Param(k)) => Some(k),
                            _ => None,
                        },
                        _ => None,
                    }) {
                        Some(k) => Var::Param(k),
                        None => return Err(self.err_hint(v, format!("'{p}' is not x, y, z, t or a param"), Some(usage.into()))),
                    },
                };
                return Ok(Value::Field(self.derivative(f, var), None));
            }
            _ => {}
        }
        unreachable!()
    }
}

impl Doc {
    /// A whole number known at load time (`what` it is, for the error):
    /// an index or a count is never rounded.
    pub(super) fn whole(&mut self, e: &Sexp, env: &mut Env, what: &str) -> R<f64> {
        let v = self.const_num(e, env)?;
        if v.fract() != 0.0 {
            return Err(self.err(e, format!("{what} must be a whole number, got {v}")));
        }
        Ok(v)
    }
}
