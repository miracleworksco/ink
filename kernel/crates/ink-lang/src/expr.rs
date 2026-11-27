//! Expressions: evaluating S-expressions into values (fields, points,
//! lists, callables), names, calls and the built-in forms.

use super::*;

impl Doc {
    /// A context at moment t with every solve solved there (`ctx_at`; Err:
    /// the first that fails, or a direction with no length there:
    /// Doc::check_conditions).
    pub fn ctx(&self, overrides: &[(String, Given)], t: f64) -> Result<Ctx, String> {
        let (ctx, failed) = self.ctx_at(overrides, t)?;
        match failed.first() {
            Some((_, e)) => Err(e.message.clone()),
            None => {
                self.check_conditions(&ctx, None, &[], false, None).map_err(|e| e.message)?;
                Ok(ctx)
            }
        }
    }

    /// Evaluate a space-independent node to a number.
    pub fn scalar(&self, n: NodeId, ctx: &Ctx) -> f64 {
        eval_point(&Tape::new(&self.g, n), ctx, [0.0; 3])
    }

    /// The fields combined by `op` (min or max), as a balanced tree. Of
    /// none: the union is the empty shape (as `integrands` reads one), the
    /// intersection all of space.
    pub fn reduce(&mut self, op: BinOp, ids: &[NodeId]) -> NodeId {
        if ids.is_empty() {
            assert!(matches!(op, BinOp::Min | BinOp::Max), "{op:?} of no fields");
            return self.g.constant(if op == BinOp::Min { 1.0 } else { -1.0 });
        }
        if ids.len() == 1 {
            return ids[0];
        }
        let m = ids.len() / 2;
        let a = self.reduce(op, &ids[..m]);
        let b = self.reduce(op, &ids[m..]);
        self.g.bin(op, a, b)
    }

    /// A fn's or lambda's argument names. The last may be followed by
    /// `...`: it then takes the remaining arguments (none or more), as a vec.
    pub(super) fn arg_names(&self, ps: &Sexp, what: &str) -> R<(Vec<String>, bool)> {
        let pl = ps.list().ok_or_else(|| self.err(ps, format!("{what} arguments must be a list of names")))?;
        let rest = pl.last().and_then(|p| p.sym()) == Some("...");
        let names = if rest { &pl[..pl.len() - 1] } else { pl };
        if rest && names.is_empty() {
            return Err(self.err(ps, "'...' follows the argument that takes the rest: (a rest ...)"));
        }
        let mut params: Vec<String> = Vec::new();
        for p in names {
            let pn = self.name_of(p)?;
            if pn == "..." {
                return Err(self.err(p, "'...' goes only at the end, after the argument that takes the rest: (a rest ...)"));
            }
            if params.contains(&pn) {
                return Err(self.err(p, format!("duplicate argument '{pn}'")));
            }
            params.push(pn);
        }
        Ok((params, rest))
    }

    /// A name: a bare word. Words quoted or in parentheses say to write
    /// them bare (quoted words joined by '-').
    pub(super) fn name_of(&self, e: &Sexp) -> R<String> {
        e.sym().map(str::to_string).ok_or_else(|| {
            let word = match &e.k {
                SK::Str(s) => parse(s).ok().and_then(|f| f.iter().map(|w| w.sym()).collect::<Option<Vec<&str>>>().filter(|w| !w.is_empty()).map(|w| w.join("-"))),
                SK::List(l) => l.first().and_then(|w| w.sym()).filter(|_| l.len() == 1).map(str::to_string),
                _ => None,
            };
            self.err_hint(e, format!("expected a name, got {e}"), word.map(|w| format!("a name is a bare word: {w}")))
        })
    }

    pub(super) fn const_v3(&mut self, e: &Sexp) -> R<V3> {
        let mut out = [0.0; 3];
        for (i, (f, at)) in self.point(e)?.into_iter().enumerate() {
            out[i] = self.g.as_const(f).ok_or_else(|| self.err(&at, "expected a constant number"))?;
        }
        Ok(out)
    }

    /// A point's coordinates, each with what it was written in: (x y z), or
    /// an expression whose value is a vec of three numbers ((vec x y z),
    /// (nth box 0), a def holding one, a call of any function). A list is a
    /// call where its head's value is a function (or it is a form): decided
    /// by that value, not by how the head is written.
    fn point(&mut self, e: &Sexp) -> R<Vec<(NodeId, Sexp)>> {
        let whole = match &e.k {
            SK::Sym(_) => true,
            SK::List(l) => match l.first() {
                Some(h) if h.sym().is_some_and(|h| SPECIAL.contains(&h)) => true,
                Some(h) => matches!(self.eval(h, &mut Vec::new()), Ok(Value::Fn(_))),
                None => false,
            },
            _ => false,
        };
        if whole {
            return match self.eval(e, &mut Vec::new())? {
                Value::Vec(v) if v.len() == 3 => v.iter().map(|c| Ok((self.as_field(c.clone(), e)?, e.clone()))).collect(),
                Value::Vec(v) => Err(self.err(e, format!("expected a point, (x y z) or a vec of three numbers, got a vec of {}", v.len()))),
                other => Err(self.err(e, format!("expected a point, (x y z) or a vec of three numbers, got a {}", other.kind()))),
            };
        }
        let parts = e.list().filter(|p| p.len() == 3).ok_or_else(|| self.err_hint(e, "expected a point (x y z)", Some("or a vec of three numbers, (vec x y z)".into())))?;
        parts.iter().map(|p| Ok((self.field(p, &mut Vec::new())?, p.clone()))).collect()
    }

    /// A point (x y z) of numbers: constants and params, kept as fields
    /// (not evaluated, so its params stay symbolic).
    /// A claim's operand: a shape (Ok), or a point (Err: `probe`).
    pub(super) fn shape_or_point(&mut self, e: &Sexp, env: &mut Env) -> R<Result<NodeId, [NodeId; 3]>> {
        match self.eval(e, env)? {
            Value::Field(f, _) => Ok(Ok(f)),
            Value::Point(p) => Ok(Err(p)),
            other => Err(self.err(e, format!("expected a shape (a field) or a point (probe), got a {}: {e}", other.kind()))),
        }
    }

    pub(super) fn corner(&mut self, e: &Sexp) -> R<[NodeId; 3]> {
        let mut out = [0; 3];
        for (i, (f, at)) in self.point(e)?.into_iter().enumerate() {
            if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T)) {
                return Err(self.err(&at, "a corner is a number (constants and params), not a field of x y z or t"));
            }
            out[i] = f;
        }
        Ok(out)
    }

    pub(super) fn field(&mut self, e: &Sexp, env: &mut Env) -> R<NodeId> {
        let v = self.eval(e, env)?;
        self.as_field(v, e)
    }

    /// Value `v` of `e` as a field.
    pub(super) fn as_field(&self, v: Value, e: &Sexp) -> R<NodeId> {
        match v {
            Value::Field(f, _) => Ok(f),
            other => Err(self.err(e, format!("expected a field (number expression), got a {}: {e}", other.kind()))),
        }
    }

    fn lookup(&mut self, s: &str, at: &Sexp, env: &Env) -> R<Value> {
        if let Some((_, v)) = env.iter().rev().find(|(k, _)| k == s) {
            return Ok(v.clone());
        }
        // Lexical: what is visible where it is read (Doc::scoped).
        self.shadows(s, at)?;
        if let Some(e) = self.later(s, at) {
            return Err(e);
        }
        if let Some(i) = self.param_index(s) {
            return Ok(Value::Field(self.g.var(Var::Param(i as u32)), self.param_said.get(&i).cloned()));
        }
        if let Some(v) = self.def(s) {
            return Ok(v.clone());
        }
        if s == "t" {
            // The moment asked, a time in seconds: a motion's formula carries it through as any quantity.
            let said = units::Said::worked(units::unit("s").expect("the second is a unit"), None, None, &self.file, at.line);
            return Ok(Value::Field(self.g.var(Var::T), Some(Rc::new(said))));
        }
        let v = match s {
            "x" => self.g.var(Var::X),
            "y" => self.g.var(Var::Y),
            "z" => self.g.var(Var::Z),
            "pi" => {
                // pi lies between its float and the next.
                let band = ink_core::interval::I::new(std::f64::consts::PI.next_down(), std::f64::consts::PI.next_up());
                let pi = self.g.constant(std::f64::consts::PI);
                self.g.real_constant(pi, ink_core::exact::Exact { ratio: None, band })
            }
            "inf" => self.g.constant(f64::INFINITY),
            _ => {
                if let Some(f) = self.function(s) {
                    return Ok(Value::Fn(Callable::User(s.to_string(), f.clone())));
                }
                if let Some(hint) = self.hidden_hint(s) {
                    return Err(self.err_hint(at, format!("unknown name '{s}'"), Some(hint)));
                }
                let o = self.resolve(s).map_err(|m| self.err(at, m))?;
                if let Some(kind) = self.body(&o).map(|_| "body").or(self.assembly(&o).map(|_| "assembly")) {
                    return Err(self.err_hint(
                        at,
                        format!("'{s}' is {} {kind}, not a value", if kind == "body" { "a" } else { "an" }),
                        Some(format!("its field is (shape-of {s})")),
                    ));
                }
                if builtin(s).is_some() || COMPARISONS.contains(&s) || FORM_FUNCTIONS.contains(&s) {
                    return Ok(Value::Fn(Callable::Builtin(s.to_string())));
                }
                if SPECIAL.contains(&s) {
                    return Err(self.err_hint(
                        at,
                        format!("'{s}' is a form, not a function: it cannot be passed as a value"),
                        Some(format!("wrap it in a function of its arguments, e.g. (lambda (a b) ({s} a b))")),
                    ));
                }
                let mut names = self.value_names();
                names.extend(env.iter().map(|(k, _)| k.clone()));
                names.extend(CONSTANTS.iter().map(|c| c.to_string()));
                // Declared as something else (a body or an assembly further on,
                // an instance) or an import's, it is not to be declared again.
                let instance = self.building.as_ref().is_some_and(|b| b.asm.instances.iter().any(|i| i.name == s) || b.later.iter().any(|l| l == s));
                let undeclared = !s.contains('.') && !instance && !self.files[self.scope.file as usize].declared.contains(s);
                // A CAD word ((reduce union shapes)) names a function as one
                // called would.
                let word = CAD_WORDS.iter().any(|w| w.0 == s).then(|| self.unknown_function_hint(s)).flatten();
                let hint = word.or_else(|| closest(s, names.iter().map(|n| n.as_str()))).or_else(|| undeclared.then(|| format!("declare it: (param {s} VALUE) or (def {s} EXPR)")));
                return Err(self.err_hint(at, format!("unknown name '{s}'"), hint));
            }
        };
        Ok(Value::Field(v, None))
    }

    pub(super) fn callable(&mut self, e: &Sexp, env: &mut Env) -> R<Callable> {
        match self.eval(e, env)? {
            Value::Fn(f) => Ok(f),
            other => Err(self.err(e, format!("expected a function, got a {}", other.kind()))),
        }
    }

    pub(super) fn list(&mut self, e: &Sexp, env: &mut Env) -> R<Vec<Value>> {
        match self.eval(e, env)? {
            Value::Vec(v) => Ok(Rc::unwrap_or_clone(v)),
            other => Err(self.err(e, format!("expected a list (vec), got a {}", other.kind()))),
        }
    }

    pub(super) fn const_num(&mut self, e: &Sexp, env: &mut Env) -> R<f64> {
        let f = self.field(e, env)?;
        self.const_of(e, f)
    }

    /// Field `f` of `e` as a number known at load time.
    pub(super) fn const_of(&mut self, e: &Sexp, f: NodeId) -> R<f64> {
        self.load_number(e, f)?.ok_or_else(|| self.err_hint(e, "expected a number known at load time", Some("constants and params work here; x, y, z and t do not".into())))
    }

    /// The value of data `d`: numbers as constants, strings, arrays as
    /// vecs, an object as data (the place, not a copy).
    pub(super) fn json_value(&mut self, d: &Data) -> Result<Value, String> {
        use serde_json::Value as J;
        let at = |v: &J, place: &dyn Fn() -> Data, me: &mut Self| -> Result<Value, String> {
            Ok(match v {
                J::Number(n) => {
                    // A record's number is written in digits, as a literal is: in a value claim its real
                    // value is theirs, so a record's 3.2 and a written 3.2 are one number, not a double
                    // and a decimal 1e-16 apart.
                    let plain = me.g.constant(n.as_f64().ok_or("number out of range")?);
                    let f = match ink_core::exact::Ratio::of_decimal(&n.to_string()) {
                        Some(r) => me.g.real_constant(plain, ink_core::exact::Exact::of(r)),
                        None => plain,
                    };
                    Value::Field(f, None)
                }
                J::Bool(b) => Value::Field(me.g.constant(if *b { 1.0 } else { 0.0 }), None),
                J::String(s) => Value::Str(s.clone()),
                J::Array(_) | J::Object(_) => return me.json_value(&place()),
                J::Null => return Err("null is not a value".into()),
            })
        };
        match d.value() {
            J::Object(_) => Ok(Value::Data(d.clone())),
            J::Array(a) => Ok(Value::list(a.iter().enumerate().map(|(i, v)| at(v, &|| d.nth(i), self)).collect::<Result<_, _>>()?)),
            v => at(v, &|| d.clone(), self),
        }
    }

    /// Whether a promise made here (its k and field, `roots`) reads an input
    /// of the call it is made in that varies over space. A promise in a
    /// function is about its body as a function of x, y, z with its inputs
    /// fixed (a sphere's radius): true of every fixed radius, it says
    /// nothing of a radius that reads x, whose slope adds to the field's
    /// (the derivative of r(x) - |p| has r'(x) in it). Inputs are the
    /// fields among what the call received (arguments and captured values,
    /// vecs included; t and params may vary, a promise holds for each of
    /// their values). Functions received are code, not inputs: a promise
    /// over what one returns is its author's.
    fn reads_moving_input(&self, roots: &[NodeId], env: &Env) -> bool {
        let n = self.inputs.last().map_or(0, |&n| n.min(env.len()));
        let mut inputs = Vec::new();
        env[..n].iter().for_each(|(_, v)| v.fields_into(&mut inputs, false));
        inputs.retain(|&i| self.g.as_const(i).is_none());
        if inputs.is_empty() {
            return false;
        }
        let read = self.g.reachable(roots);
        inputs.iter().any(|&i| read.binary_search(&i).is_ok() && self.g.depends_on(i, |v| matches!(v, Var::X | Var::Y | Var::Z)))
    }

    /// Call a function value with evaluated arguments.
    pub(super) fn call(&mut self, f: &Callable, args: Vec<Value>, at: &Sexp) -> R<Value> {
        let (label, params, rest, body, file, mut env) = match f {
            Callable::Builtin(name) if FORM_FUNCTIONS.contains(&name.as_str()) => return self.call_form(name, args, at),
            Callable::Builtin(name) => return self.apply_builtin(name, args, at),
            Callable::User(name, d) => (name.clone(), d.params.clone(), d.rest, &d.body, (d.file.clone(), d.scope), Vec::new()),
            Callable::Closure(c) => ("lambda".to_string(), c.params.clone(), c.rest, &c.body, (c.file.clone(), c.scope), c.env.clone()),
        };
        let mut args = args;
        if rest && args.len() + 1 >= params.len() {
            let more = args.split_off(params.len() - 1);
            args.push(Value::list(more));
        }
        if args.len() != params.len() {
            let (n, ps) = if rest { (format!("{} or more", params.len() - 1), format!("{} ...", params.join(" "))) } else { (params.len().to_string(), params.join(" ")) };
            return Err(self.err(at, format!("'{label}' takes {n} argument(s) ({ps}), got {}", args.len())));
        }
        env.extend(params.into_iter().zip(args));
        if self.depth >= MAX_CALL_DEPTH {
            return Err(self.err(at, format!("call depth exceeded {MAX_CALL_DEPTH} in '{label}' (recursion is not allowed)")));
        }
        self.depth += 1;
        self.inputs.push(env.len());
        self.calls.push((label, self.file.clone(), at.line));
        // The body reads the names of where it was written (Scope).
        let saved = (std::mem::replace(&mut self.file, file.0), std::mem::replace(&mut self.scope, file.1));
        let local = std::mem::replace(
            &mut self.local_assembly,
            match f {
                Callable::Closure(c) => c.assembly.clone(),
                _ => None,
            },
        );
        let r = self.eval(body, &mut env);
        self.local_assembly = local;
        (self.file, self.scope) = saved;
        let (label, ..) = self.calls.pop().expect("pushed above");
        self.inputs.pop();
        self.depth -= 1;
        // An error inside the body: record this call, so it can be traced
        // back to the author's line.
        r.map_err(|mut e| {
            e.called_from.push((label, self.file.to_string(), at.line, at.col));
            e
        })
    }

    /// Call a form that is a function (FORM_FUNCTIONS) with evaluated
    /// arguments: the form, written at the call site, reading each argument
    /// by a name bound to its value, as a call of `(lambda (a ...) (form a
    /// ...))` would.
    fn call_form(&mut self, name: &str, args: Vec<Value>, at: &Sexp) -> R<Value> {
        let sx = |k| Sexp { k, line: at.line, col: at.col };
        let mut items = vec![sx(SK::BuiltIn(name.to_string()))];
        let mut env: Env = Vec::with_capacity(args.len());
        for (i, v) in args.into_iter().enumerate() {
            let a = format!("arg-{}", i + 1);
            items.push(sx(SK::Sym(a.clone())));
            env.push((a, v));
        }
        self.inputs.push(env.len());
        let r = self.eval(&sx(SK::List(items.into())), &mut env);
        self.inputs.pop();
        r
    }

    /// Apply a kernel op or comparison. Vec arguments broadcast element-wise
    /// (scalars repeat), so vector maths reads naturally: (- p a), (* v 2).
    pub(super) fn apply_builtin(&mut self, name: &str, args: Vec<Value>, at: &Sexp) -> R<Value> {
        if let Some(n) = args.iter().find_map(|v| if let Value::Vec(v) = v { Some(v.len()) } else { None }) {
            // Every vec the same length (an empty one too), and as many
            // arguments as the op takes.
            let takes = match builtin(name) {
                Some(Builtin::Un(_)) => Some(1),
                Some(Builtin::Bin(_)) => Some(2),
                Some(Builtin::Select) => Some(3),
                Some(Builtin::Minus) if !(1..=2).contains(&args.len()) => Some(2),
                _ if COMPARISONS.contains(&name) => Some(2),
                _ => None,
            };
            if let Some(k) = takes.filter(|&k| k != args.len()) {
                return Err(self.err(at, format!("'{name}' takes {k} argument(s), got {}", args.len())));
            }
            if let Some(m) = args.iter().find_map(|v| if let Value::Vec(v) = v { Some(v.len()).filter(|&m| m != n) } else { None }) {
                return Err(self.err(at, format!("'{name}' on vecs of lengths {n} and {m}")));
            }
            let mut out = Vec::with_capacity(n);
            for i in 0..n {
                let row = args
                    .iter()
                    .map(|v| match v {
                        Value::Vec(v) if v.len() == n => Ok(v[i].clone()),
                        Value::Vec(v) => Err(self.err(at, format!("'{name}' on vecs of lengths {n} and {}", v.len()))),
                        other => Ok(other.clone()),
                    })
                    .collect::<R<Vec<_>>>()?;
                out.push(self.apply_builtin(name, row, at)?);
            }
            return Ok(Value::list(out));
        }
        if let ("==", [Value::Str(a), Value::Str(b)]) = (name, args.as_slice()) {
            return Ok(Value::Field(self.g.constant(if a == b { 1.0 } else { 0.0 }), None));
        }
        // A choice of frames is one of their numbers: a part that changes
        // hands (picked up at a moment) is placed by (select c F G).
        if let (Some(Builtin::Select), [Value::Field(c, _), Value::Frame(a), Value::Frame(b)]) = (builtin(name), args.as_slice()) {
            let mut pick = |p: NodeId, q: NodeId| self.g.select(*c, p, q);
            let r = std::array::from_fn(|i| std::array::from_fn(|k| pick(a.r[i][k], b.r[i][k])));
            return Ok(Value::Frame(Xf { r, t: std::array::from_fn(|i| pick(a.t[i], b.t[i])) }));
        }
        let mut fs = Vec::with_capacity(args.len());
        for v in &args {
            match v {
                Value::Field(f, _) => fs.push(*f),
                other => {
                    let frames = if name == "select" { ", or a number and two frames" } else { "" };
                    let got = match other {
                        Value::Str(s) => format!("the string \"{s}\""),
                        o => format!("a {}", o.kind()),
                    };
                    return Err(self.err(at, format!("'{name}' expects numbers or vecs{frames}, got {got}")));
                }
            }
        }
        let arity = |n: usize, d: &Doc| {
            if fs.len() != n { Err(d.err(at, format!("'{name}' takes {n} argument(s), got {}", fs.len()))) } else { Ok(()) }
        };
        if COMPARISONS.contains(&name) {
            arity(2, self)?;
            let (one, zero) = (self.g.constant(1.0), self.g.constant(0.0));
            let (a, b) = (fs[0], fs[1]);
            let (ab, ba) = (self.g.sub(a, b), self.g.sub(b, a));
            let id = match name {
                "<" => self.g.select(ab, one, zero),
                ">=" => self.g.select(ab, zero, one),
                ">" => self.g.select(ba, one, zero),
                "<=" => self.g.select(ba, zero, one),
                _ => {
                    let ne = self.g.select(ba, zero, one);
                    self.g.select(ab, zero, ne)
                }
            };
            return Ok(Value::Field(id, None));
        }
        let id = match builtin(name).expect("builtin") {
            Builtin::Un(op) => {
                arity(1, self)?;
                // An angle said with its unit is in degrees. A sine, cosine or tangent takes radians, so is given it in
                // them; `rad->deg` of it is its number of degrees, which it already is.
                let angle = matches!(&args[0], Value::Field(_, Some(s)) if s.dims().is_some_and(|d| d.0 == [0, 0, 0, 0, 1]));
                match op {
                    UnOp::Deg if angle => fs[0],
                    UnOp::Sin | UnOp::Cos | UnOp::Tan if angle => {
                        let radians = self.radians(fs[0]);
                        let radians = radians.unwrap_or_else(|| self.g.un(UnOp::Rad, fs[0]));
                        self.g.un(op, radians)
                    }
                    _ => self.g.un(op, fs[0]),
                }
            }
            Builtin::Bin(op) => {
                arity(2, self)?;
                self.g.bin(op, fs[0], fs[1])
            }
            Builtin::Variadic(op) => {
                if fs.is_empty() {
                    return Err(self.err(at, format!("'{name}' needs at least one argument")));
                }
                self.reduce(op, &fs)
            }
            Builtin::Minus => match fs.as_slice() {
                [] => return Err(self.err(at, "'-' needs at least one argument")),
                [a] => self.g.neg(*a),
                // (- a b c ...) = a - b - c - ...
                [a, rest @ ..] => rest.iter().fold(*a, |acc, b| self.g.sub(acc, *b)),
            },
            Builtin::Select => {
                arity(3, self)?;
                self.g.select(fs[0], fs[1], fs[2])
            }
        };
        let (said, times) = self.said_of(name, &args, at)?;
        // The number as worked out, in the unit of the quantity it is.
        let id = self.scaled(id, times);
        Ok(Value::Field(id, said))
    }

    /// An angle in degrees `f`, in radians, where a constant of it can take the degree: a rate times t, (* 180deg/s
    /// t), over a number or times one, ((* 3 (/ (* 360deg/s t) period))). The sine is then of the number a bare
    /// (* pi t) gave, to the last bit, however the rest is grouped. None: no such constant.
    fn radians(&mut self, f: NodeId) -> Option<NodeId> {
        use crate::graph::Node;
        match self.g.nodes[f as usize] {
            Node::Bin(BinOp::Mul, a, b) => {
                let (k, y) = match (self.g.as_const(a), self.g.as_const(b)) {
                    (Some(_), None) => (a, b),
                    (None, Some(_)) => (b, a),
                    _ => return None,
                };
                Some(match self.radians(y) {
                    Some(inner) => self.g.mul(k, inner),
                    None => {
                        // k degrees in radians: pi's interval, not its float, carried (a value claim reads it).
                        let kv = self.g.as_const(k).expect("a constant");
                        let pi = ink_core::interval::I::new(std::f64::consts::PI.next_down(), std::f64::consts::PI.next_up());
                        let kr = self.g.real_value(k).map_or(ink_core::interval::I::point(kv), |e| e.band);
                        let band = ink_core::interval::bin(BinOp::Mul, kr, ink_core::interval::bin(BinOp::Div, pi, ink_core::interval::I::point(180.0)).0).0;
                        let turned = self.g.constant(kv * (std::f64::consts::PI / 180.0));
                        let turned = self.g.real_constant(turned, ink_core::exact::Exact { ratio: None, band });
                        self.g.mul(turned, y)
                    }
                })
            }
            Node::Bin(BinOp::Div, a, b) => {
                let over = self.radians(a)?;
                Some(self.g.bin(BinOp::Div, over, b))
            }
            _ => None,
        }
    }

    /// `id` times the scale a product's units come to (millimetres a metre, newtons a gram millimetre a second
    /// squared). A power of ten under one divides by its whole number: 0.001 is not a thousandth to the last
    /// bit, and a motion's joins are judged in real numbers, where a profile said in mm/s2 must still meet.
    pub(crate) fn scaled(&mut self, id: NodeId, times: f64) -> NodeId {
        if (times - 1.0).abs() <= 1e-12 {
            return id;
        }
        let over = (1.0 / times).round();
        let whole = over >= 10.0 && (over * times - 1.0).abs() < 1e-12 && (over.log10() - over.log10().round()).abs() < 1e-9;
        if whole {
            let by = self.g.constant(over);
            self.g.bin(BinOp::Div, id, by)
        } else {
            let by = self.g.constant(times);
            self.g.mul(id, by)
        }
    }

    pub(super) fn eval(&mut self, e: &Sexp, env: &mut Env) -> R<Value> {
        let items = match &e.k {
            SK::Num(v, unit, exact) => {
                if self.scope != Scope::STD && self.wild_written.is_none() && crate::graph::wild(*v) {
                    self.wild_written = Some((*v, format!("{}:{}", self.file, e.line)));
                }
                // What it is in real numbers, as written (sexp: a literal's); its float everywhere but a value
                // claim on constants. A number the program made is its float.
                let plain = if self.scope == Scope::STD { self.g.own(*v) } else { self.g.constant(*v) };
                let f = match exact {
                    Some(real) => self.g.real_constant(plain, *real),
                    None => plain,
                };
                let said = unit.as_ref().map(|w| match w {
                    sexp::Written::Here(unit, written) => {
                        Rc::new(units::Said { explicit_when: None, unit, written: *written, file: std::sync::Arc::from(&*self.file), line: e.line, powers: None })
                    }
                    sexp::Written::From(s) => Rc::new((**s).clone()),
                });
                return Ok(Value::Field(f, said));
            }
            SK::Str(s) => return Ok(Value::Str(s.clone())),
            SK::Sym(s) => return self.lookup(s, e, env),
            SK::For(_) => return Err(self.err(e, "for is a form, not an expression")),
            SK::BuiltIn(s) => return Err(self.err(e, format!("the built-in {s} is written only at the head of a form"))),
            SK::List(items) => items,
        };
        let Some(head) = items.first() else { return Err(self.err(e, "empty list () is not an expression")) };
        let args = &items[1..];
        let (name, shadowed) = match &head.k {
            SK::BuiltIn(n) => (n.as_str(), false),
            SK::Sym(n) => (n.as_str(), self.shadows(n, head)?),
            // ((lambda (x) ...) 3), or any expression that yields a function.
            _ => {
                let Value::Fn(f) = self.eval(head, env)? else {
                    return Err(self.err(head, format!("expected a function at the head of {e}")));
                };
                let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
                return self.call(&f, vals, e);
            }
        };
        match name {
            // The file's own definition, called below.
            _ if shadowed => {}
            "let" => {
                let [binds, body] = args else { return Err(self.err(e, "usage: (let ((name expr) ...) body)")) };
                let bl = binds.list().ok_or_else(|| self.err(binds, "let bindings must be a list of (name expr)"))?;
                let mark = env.len();
                for b in bl {
                    let [n, v] = b.list().unwrap_or(&[]) else { return Err(self.err(b, "let binding must be (name expr)")) };
                    let n = self.name_of(n)?;
                    if is_reserved(&n) {
                        return Err(self.err_hint(b, format!("'{n}' is a built-in name"), Some(reserved_hint())));
                    }
                    let v = self.eval(v, env)?;
                    env.push((n, v));
                }
                let r = self.eval(body, env);
                env.truncate(mark);
                return r;
            }
            "vec" => {
                let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
                return Ok(Value::list(vals));
            }
            "unit" => {
                let [v] = args else { return Err(self.err(e, "usage: (unit vec): the unit vector along a vec of fields")) };
                let Value::Vec(items) = self.eval(v, env)? else { return Err(self.err(v, "unit expects a vec")) };
                let fs = items.iter().map(|i| if let Value::Field(f, _) = i { Some(*f) } else { None }).collect::<Option<Vec<_>>>();
                let fs = fs.filter(|f| !f.is_empty()).ok_or_else(|| self.err(v, "unit expects a vec of fields (numbers)"))?;
                let hint = "(unit v) is v / |v|, which needs v nonzero (a path's tangent is zero between two equal points, and a sweep's first normal where `up` lies along its first segment); where it reads a param or t, keep it nonzero over their ranges";
                let u = self.unit(&fs, e, &format!("the vector of {e}"), hint, false)?;
                return Ok(Value::list(u.into_iter().map(|f| Value::Field(f, None)).collect()));
            }
            "lambda" | "map" | "reduce" | "fold" | "flatten" | "apply" | "range" | "len" | "concat" | "char-codes" | "slice" | "if" | "read-json" | "get" | "diff" => {
                return self.eval_author(name, e, args, env);
            }
            "nth" => {
                let [v, i] = args else { return Err(self.err(e, "usage: (nth vec index)")) };
                let Value::Vec(items) = self.eval(v, env)? else { return Err(self.err(v, "nth expects a vec")) };
                let i = self.whole(i, env, "an nth index")?;
                // Negative indices count from the end: -1 is the last element.
                let k = if i < 0.0 { items.len() as f64 + i } else { i };
                return (k >= 0.0)
                    .then(|| items.get(k as usize).cloned())
                    .flatten()
                    .ok_or_else(|| self.err(e, format!("index {i} out of range for vec of length {}", items.len())));
            }
            "data" => {
                let v = match args {
                    [single] => single.to_json(),
                    _ => serde_json::Value::Array(args.iter().map(|a| a.to_json()).collect()),
                };
                return Ok(Value::Data(Data::new(v)));
            }
            "warp" => {
                if args.len() != 4 && args.len() != 5 {
                    return Err(self.err(e, "usage: (warp expr fx fy fz [ft])"));
                }
                let v = self.eval(&args[0], env)?;
                let mut map = Vec::new();
                for (var, a) in [Var::X, Var::Y, Var::Z, Var::T].into_iter().zip(&args[1..]) {
                    map.push((var, self.field(a, env)?));
                }
                // A warp says where a shape's points are read from, not
                // where a point goes: that is its inverse, which rounding
                // can put anywhere (a matrix singular in real numbers whose
                // rounded determinant is not 0).
                if has_point(&v) {
                    return Err(self.err_hint(
                        e,
                        "a warp moves a shape, not a point (probe): it says where the shape is read from, and a point would need its inverse",
                        Some("put a point at a frame, which moves it forward: (place (probe x y z) (at (move dx dy dz) (rot-z deg)))".into()),
                    ));
                }
                if args.len() == 5 && !self.ext.is_empty() {
                    estimate::retime(self, e, &v, &mut map)?;
                }
                let warped = self.warp_value(&v, &map);
                // A field read elsewhere is the same values: what checks the original (a density's physical
                // range) checks its copy, though the copy is new arithmetic.
                if let (Value::Field(f, _), Value::Field(w, _)) = (&v, &warped)
                    && f != w
                {
                    self.warped.push((*f, *w));
                    if let [(_, fx), (_, fy), (_, fz)] = map[..] {
                        self.warps.entry(*w).or_insert((*f, [fx, fy, fz]));
                    }
                }
                return Ok(warped);
            }
            "probe" => {
                let [px, py, pz] = args else {
                    return Err(self.err_hint(e, "wrong arguments", Some("usage: (probe x y z): a point, for claims that it is in a shape or not".into())));
                };
                let mut at = [0; 3];
                for (k, a) in [px, py, pz].into_iter().enumerate() {
                    let f = self.field(a, env)?;
                    if self.g.depends_on(f, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                        return Err(self.err(a, "a point's coordinates are numbers (constants and params), not fields of x y z"));
                    }
                    at[k] = f;
                }
                return Ok(Value::Point(at));
            }
            "tag" => {
                let usage = "usage: (tag name field [opaque]) names a construction for provenance (faces); the name is a symbol, a string, \
                             or an expression giving a string; opaque: points under it report this name only, not the names inside";
                let (n, f, opaque) = match args {
                    [n, f] => (n, f, false),
                    [n, f, o] if o.sym() == Some("opaque") => (n, f, true),
                    _ => return Err(self.err_hint(e, "bad tag", Some(usage.into()))),
                };
                let name = match &n.k {
                    SK::Sym(s) | SK::Str(s) => s.clone(),
                    SK::List(_) => match self.eval(n, env)? {
                        Value::Str(s) => s,
                        other => return Err(self.err_hint(n, format!("a tag name must be a string, not a {}", other.kind()), Some(usage.into()))),
                    },
                    _ => return Err(self.err_hint(n, "tag name must be a name or string", Some(usage.into()))),
                };
                // A name on a value: the value, and its unit, as they were.
                let v = self.eval(f, env)?;
                let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
                let f = self.as_field(v, f)?;
                return Ok(Value::Field(if opaque { self.g.tag_opaque(&name, f) } else { self.g.tag(&name, f) }, said));
            }
            "real" => {
                let [a] = args else { return Err(self.err(e, "usage: (real value) uses the value in real arithmetic")) };
                let v = self.eval(a, env)?;
                let Value::Field(n, unit) = v else { return Err(self.err(a, "real expects a scalar field")) };
                return Ok(Value::Field(self.g.real(n), unit));
            }
            "lipschitz" => {
                let [k, f] = args else { return Err(self.err(e, "usage: (lipschitz k field) promises |grad field| <= k in x y z (exact distance: 1); proofs trust it")) };
                // A promise about a value: the value, and its unit, as they were.
                let v = self.eval(f, env)?;
                let said = if let Value::Field(_, s) = &v { s.clone() } else { None };
                let (k, f) = (self.field(k, env)?, self.as_field(v, f)?);
                if self.reads_moving_input(&[k, f], env) {
                    return Ok(Value::Field(f, said));
                }
                return Ok(Value::Field(self.g.lipschitz(f, k), said));
            }
            "corner-reach" => {
                let (a, b, r, guard) = match args {
                    [a, b, r] => (self.field(a, env)?, self.field(b, env)?, self.field(r, env)?, self.g.constant(-1.0)),
                    [a, b, r, g] => (self.field(a, env)?, self.field(b, env)?, self.field(r, env)?, self.field(g, env)?),
                    _ => return Err(self.err(e, "usage: (corner-reach a b r [guard])")),
                };
                return Ok(Value::Field(self.g.corner_reach(a, b, r, guard), None));
            }
            "attr-of" => {
                let usage = "usage: (attr-of PATH NAME): the attribute NAME of the body placed at PATH (or the body PATH), as the body writes it: a record, string or number";
                let [p, k] = args else { return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))) };
                let k = self.name_of(k)?;
                let (path, items) = self.scoped_items(&self.name_of(p)?).map_err(|m| self.err(p, m))?;
                let [it] = items.as_slice() else {
                    return Err(self.err_hint(p, format!("'{path}' places {} bodies; (attr-of PATH NAME) reads one", items.len()), Some(usage.into())));
                };
                let body = it.instance.as_ref().map_or(path.clone(), |(b, _)| b.clone());
                // Placed as the part is, as (attr NAME) reads it: a field attribute where the part lies.
                let attrs = self.placed_attrs(it);
                return match attrs.iter().find(|(n, _)| *n == k) {
                    Some((_, v)) => Ok(v.clone()),
                    None => {
                        let names: Vec<&str> = attrs.iter().map(|(n, _)| n.as_str()).collect();
                        Err(self.err_hint(p, format!("body '{body}' (at '{path}') has no attribute '{k}'"), Some(format!("its attributes: {}", names.join(", ")))))
                    }
                };
            }
            "attr" => {
                let usage = "(attr NAME) is an attribute of the body being integrated (its own, else its material's), placed with it: it is read in the fields of an \
                             integral or mean claim over a body, an assembly or an instance path, e.g. (mean root z (weight (attr density)) < 80)";
                let [k] = args else { return Err(self.err_hint(e, "usage: (attr NAME)", Some(usage.into()))) };
                let k = self.name_of(k)?;
                let Some((label, attrs)) = self.integrand_of.clone() else { return Err(self.err_hint(e, "(attr NAME) outside an integrand", Some(usage.into()))) };
                let attrs = match attrs {
                    Ok(attrs) => attrs,
                    Err(why) => return Err(self.err_hint(e, format!("the target {label} {why}"), Some(usage.into()))),
                };
                return match materials::required(&mut self.g, &label, &attrs, &k) {
                    Ok(f) => {
                        self.want(&[f]);
                        Ok(Value::Field(f, None))
                    }
                    Err(m) => Err(self.err_hint(e, m.message, m.hint)),
                };
            }
            "grid-bake" | "import-stl" | "redistance" | "bounds" | "profile-reach" => return self.eval_services(name, e, args, env),
            "print-overhang-mask" if self.scope == Scope::STD => {
                let [f, bed, angle, h, s] = args else { return Err(self.err(e, "expected field, bed, angle, height and ranged depth")) };
                let (f, bed, angle, h, s) = (self.field(f, env)?, self.field(bed, env)?, self.field(angle, env)?, self.field(h, env)?, self.field(s, env)?);
                let graph::Node::Var(Var::Param(k)) = self.g.nodes[s as usize] else { return Err(self.err(e, "expected a ranged section position")) };
                self.support_depths.insert(k);
                let unsupported = self.overhang_support(f, bed, angle, h);
                let g = &mut self.g;
                let z = g.var(Var::Z);
                let above = g.sub(bed, z);
                let material = g.max(f, above);
                let possible = g.max(material, unsupported);
                let depth = g.bin_real(BinOp::Mul, h, s);
                let slope = bridge::print_slope(g, angle);
                let width = g.bin_real(BinOp::Mul, depth, slope);
                let section = support::section(g, f, width, depth, &self.params[k as usize].name);
                let z_at = g.sub(z, depth);
                let above = g.sub(bed, z_at);
                let section = g.max(section, above);
                let certain = g.max(f, section);
                let mask = g.max(possible, certain);
                self.print_bounds.insert(mask, support::Overhang { f, bed, h, slope, k, certain, possible });
                return Ok(Value::Field(mask, None));
            }
            "material" | "by-body" => return self.eval_materials(name, e, args, env),
            "check" => return self.eval_check(e, args, env),
            "offset-shape" => return self.eval_offset(e, args, env),
            // Only std can certify that its two fields enclose one region.
            "bridge-region" if self.scope == Scope::STD => {
                let [certain, possible] = args else { return Err(self.err(e, "usage: (bridge-region CERTAIN POSSIBLE)")) };
                let certain = self.field(certain, env)?;
                let possible = self.field(possible, env)?;
                return Ok(Value::Region { certain, possible });
            }
            "bridge-support" if self.scope == Scope::STD => {
                let [f, bed, angle, height, depth, limit, step] = args else { return Err(self.err(e, "usage: (bridge-support FIELD BED ANGLE HEIGHT DEPTH LIMIT STEP)")) };
                let (f, bed) = (self.field(f, env)?, self.field(bed, env)?);
                let (angle, height, depth) = (self.field(angle, env)?, self.field(height, env)?, self.field(depth, env)?);
                let (limit, step) = (self.field(limit, env)?, self.field(step, env)?);
                let flags = self.bridge_support(f, bed, angle, height, depth, limit, step);
                return Ok(Value::list(flags.map(|f| Value::Field(f, None)).into()));
            }
            "ranged" => {
                let usage = "usage: (ranged lo hi [at]): a new variable, which every claim reading it ranges over [lo, hi] (a range of the question: \
                             how far a tool is withdrawn, how far a view looks), and no run pins (--nominal, --param); at (lo if not given) where nothing ranges (eval, mesh, view)";
                if !(2..=3).contains(&args.len()) {
                    return Err(self.err_hint(e, "bad ranged", Some(usage.into())));
                }
                let ends = args.iter().map(|x| self.eval(x, env)).collect::<R<Vec<_>>>()?;
                let said = self.one_said("ranged", &ends, e)?;
                let end = |i: usize| -> R<f64> {
                    let f = self.as_field(ends[i].clone(), &args[i])?;
                    self.g.as_const(f).filter(|v| v.is_finite()).ok_or_else(|| self.err_hint(&args[i], "a range's ends are constant finite numbers", Some(usage.into())))
                };
                let (lo, hi) = (end(0)?, end(1)?);
                let at = if args.len() == 3 { end(2)? } else { lo };
                if !(lo <= at && at <= hi) {
                    return Err(self.err_hint(e, format!("(ranged {lo} {hi} {at}) needs lo <= at <= hi"), Some(usage.into())));
                }
                // Named by its form, "t-slot-reaches 0", the first number free:
                // no declared name has a space (a def that is one renames it).
                let base = &self.form_name;
                let name = (0..).map(|i| format!("{base} {i}")).find(|n| !self.ranged.iter().any(|r| self.params[r.0].name == *n)).unwrap();
                // A value given it by that name (--param "b 0=2mm") says
                // its unit as one given a declared param does.
                let said = self.param_given(&name, at, said, e)?;
                self.add_param(Param { name, default: at, tol: None });
                self.ranged.push((self.params.len() - 1, (lo, hi)));
                if let Some(s) = said.clone() {
                    self.param_said.insert(self.params.len() - 1, s);
                }
                return Ok(Value::Field(self.g.var(Var::Param(self.params.len() as u32 - 1)), said));
            }
            "estimate-when" => return estimate::read(self, env, e, args),
            "error" => {
                let usage = "usage: (error \"text\" [VALUE ...]): a load error saying text and the values after it (strings as written, numbers, records), where it is evaluated (a library refusing an input no part has)";
                let [m, rest @ ..] = args else { return Err(self.err_hint(e, "bad error", Some(usage.into()))) };
                let Value::Str(mut m) = self.eval(m, env)? else { return Err(self.err_hint(e, "an error's text is a string", Some(usage.into()))) };
                for x in rest {
                    let v = self.eval(x, env)?;
                    m += &self.value_text(&v);
                }
                return Err(self.err(e, m));
            }
            "met-by" => {
                let usage = "usage: (met-by PATH): the name of the frame the instance PATH was placed by (its CFRAME, met face to face with WHERE)";
                let [p] = args else { return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))) };
                let path = self.path_of(p)?;
                // Resolved as attr-of resolves a path, with its errors.
                self.items(&path).map_err(|m| self.err(p, m))?;
                let (parent, name) = path.rsplit_once('/').ok_or_else(|| self.err_hint(p, format!("'{path}' is no instance path"), Some(usage.into())))?;
                let parent = self.placement(parent).map(|(_, child)| child);
                let met = parent.as_deref().and_then(|a| self.assembly(a)).and_then(|a| a.instances.iter().find(|i| i.name == name)).map(|i| i.met.clone());
                return match met {
                    Some(Some(n)) => Ok(Value::Str(n)),
                    Some(None) => Err(self.err_hint(p, format!("'{path}' is placed by no frame of its own (no CFRAME)"), Some(usage.into()))),
                    None => Err(self.err_hint(p, format!("'{path}' is no instance path"), Some(usage.into()))),
                };
            }
            "relation" => return self.eval_relation(e, args, env),
            "polygon" => return self.eval_polygon(e, args, env),
            "font" | "glyph-2d" | "glyph-metrics" | "glyph-run" | "glyph-reach" => return self.eval_font(name, e, args, env),
            "revolve-z" => return self.eval_revolve(e, args, env),
            "study-field" => return self.eval_result(e, args, env),
            "at" => return Ok(Value::Frame(self.eval_at(args, env)?)),
            "sin-cos-deg" => {
                let [arg] = args else { return Err(self.err(e, "sin-cos-deg takes one angle in degrees")) };
                let deg = self.field(arg, env)?;
                return Ok(Value::list(turn_sin_cos(&mut self.g, deg).map(|n| Value::Field(n, None)).to_vec()));
            }
            "axes" => return Ok(Value::Frame(self.eval_axes(e, args, env)?)),
            "frame-of" => {
                let usage = "usage: (frame-of PATH [NAME]): the frame NAME declared by the body or assembly PATH ends at, placed; NAME omitted, PATH's placement";
                let (p, n) = match args {
                    [p] => (p, None),
                    [p, n] => (p, Some(self.name_of(n)?)),
                    _ => return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))),
                };
                let path = self.name_of(p)?;
                return self.frame_of(&path, n.as_deref()).map(Value::Frame).map_err(|m| self.err(e, m));
            }
            "place" => {
                let [v, f] = args else { return Err(self.err_hint(e, "wrong arguments", Some("usage: (place SHAPE F): SHAPE, written in F's coordinates, put at F".into()))) };
                let v = self.eval(v, env)?;
                let f = self.frame_value(f, env, "place's second argument")?;
                return Ok(self.to_world(&v, &f));
            }
            "shape-of" => {
                let usage = "usage: (shape-of name-or-path [construction | (except PATH ...)]): an instance's world field, that of a named construction inside it, or its parts but those under each PATH; a (for var list PATH ...) among them repeats its paths";
                // (except PATH ...): the parts under the path but those at or under each PATH,
                // a `for` among them repeating its paths for a list known at load time.
                let except = match args {
                    [_, n] => n.list().filter(|l| l.first().and_then(|w| w.sym()) == Some("except")).map(|l| l[1..].to_vec()),
                    _ => None,
                };
                let except = match except {
                    Some(left) => Some(self.expand_for_in(&left, env)?),
                    None => None,
                };
                let (p, part) = match args {
                    [p] => (p, None),
                    [p, _] if except.is_some() => (p, None),
                    [p, n] => (p, Some(self.name_of(n)?)),
                    _ => return Err(self.err_hint(e, "wrong arguments", Some(usage.into()))),
                };
                let key = (self.scope, self.name_of(p)?);
                if part.is_none()
                    && except.is_none()
                    && let Some(&f) = self.shape_of.get(&key)
                {
                    return Ok(Value::Field(f, None));
                }
                let (path, mut items) = self.scoped_items(&key.1).map_err(|m| self.err(p, m))?;
                if let Some(left) = &except {
                    if left.is_empty() {
                        return Err(self.err_hint(e, "(except) names no path", Some(usage.into())));
                    }
                    // At or under `out`: the label itself, or one below it past a '/'.
                    let under = |label: &str, out: &str| label == out || label.strip_prefix(out).is_some_and(|r| r.starts_with('/'));
                    let mut outs: Vec<String> = Vec::new();
                    for x in left {
                        let out = self.path_of(x)?;
                        if let Some(o) = outs.iter().find(|o| under(&out, o) || under(o, &out)) {
                            let what = if *o == out {
                                "is named twice".to_string()
                            } else if under(&out, o) {
                                format!("is under '{o}', already left out")
                            } else {
                                format!("holds '{o}', named too")
                            };
                            return Err(self.err_hint(x, format!("(except ...): '{out}' {what}"), Some("name each part, or the assembly holding it, once".into())));
                        }
                        if !items.iter().any(|it| under(&it.label, &out)) {
                            let hint = format!("a path from '{path}' down to a part it places, as {}", items.first().map_or("root/part", |it| it.label.as_str()));
                            return Err(self.err_hint(x, format!("'{out}' names no part of '{path}'"), Some(hint)));
                        }
                        outs.push(out);
                    }
                    items.retain(|it| !outs.iter().any(|o| under(&it.label, o)));
                    if items.is_empty() {
                        return Err(self.err_hint(
                            e,
                            format!("(except ...) leaves no part of '{path}'"),
                            Some("a claim against nothing holds whatever it asks: leave out less".into()),
                        ));
                    }
                }
                let mut shapes = Vec::new();
                for it in &items {
                    let Some(part) = &part else {
                        shapes.push(self.placed_field(it));
                        continue;
                    };
                    let Some((body, xf)) = it.instance.clone() else { return Err(self.err(p, format!("'{path}' is not a body or an instance"))) };
                    let local = self.body(&body).unwrap().shape;
                    let found = self.construction_names(part).iter().map(|n| self.named_in(local, n)).find(|f| !f.is_empty()).unwrap_or_default();
                    if found.is_empty() {
                        let names = self.named_in_all(local);
                        return Err(self.err_hint(p, format!("body '{body}' has no piece built as (tag {part} ...)"), Some(format!("its pieces: {}", names.join(", ")))));
                    }
                    // Every occurrence (a patterned def appears once per copy), placed like the instance,
                    // under the body's own tag: a piece of a part is of that part (Doc::parts_of).
                    let f = self.reduce(BinOp::Min, &found);
                    let f = match self.g.nodes[local as usize] {
                        crate::graph::Node::Tag(t, _) => {
                            let name = self.g.tag_names[t as usize].clone();
                            self.g.tag(&name, f)
                        }
                        _ => f,
                    };
                    let map = self.inverse_map(&xf);
                    let world = self.g.subst(f, &map);
                    self.record_placed(world, f, &xf);
                    self.placed_labels.entry(world).or_insert_with(|| it.label.clone());
                    shapes.push(world);
                }
                // Of no items (an empty assembly), the empty shape.
                let f = self.reduce(BinOp::Min, &shapes);
                if part.is_none() && except.is_none() {
                    self.shape_of.insert(key, f);
                }
                return Ok(Value::Field(f, None));
            }
            _ => {}
        }
        // A local name bound to a function value (e.g. an argument of a lambda).
        if let Some((_, Value::Fn(f))) = env.iter().rev().find(|(k, _)| k == name) {
            let f = f.clone();
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.call(&f, vals, e);
        }
        if !shadowed && (builtin(name).is_some() || COMPARISONS.contains(&name)) {
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            return self.apply_builtin(name, vals, e);
        }
        // A package's value (`(reaction PATH)`), read in this scope.
        if let Some((p, h)) = self.registry.value(name) {
            let mut r = Reader::in_env(self, env.clone());
            return p.value(&mut r, e, h, args);
        }
        // Only here is this head a user call. Resolve it in the file that
        // wrote it: syntax readers above (and placements, claims, packages)
        // must see their own words, even if a module defines one as a fn.
        if let Some(mut err) = self.later(name, head) {
            // (def t-span (t-span 0 2)): the def's own name, read as a call.
            err.hint = self.top_form_hint(e, "an expression").or(err.hint);
            return Err(err);
        }
        // A def whose value is a function is called as a local name holding one is.
        let f = self.function(name).map(|f| Callable::User(name.to_string(), f.clone())).or_else(|| match self.def(name) {
            Some(Value::Fn(f)) => Some(f.clone()),
            _ => None,
        });
        if let Some(f) = f {
            let vals = args.iter().map(|a| self.eval(a, env)).collect::<R<Vec<_>>>()?;
            // Called with the wrong number of arguments: the built-in may
            // have been meant.
            return self.call(&f, vals, e).map_err(|mut err| {
                if shadowed && err.called_from.is_empty() && (err.line, err.col) == (e.line, e.col) {
                    err.hint = Some(self.shadow_hint(name));
                }
                err
            });
        }
        let shadow = shadowed.then(|| self.shadow_hint(name));
        if let Some(v) = self.def(name) {
            return Err(self.err_hint(head, format!("'{name}' is a def holding a {}, not a function", v.kind()), shadow));
        }
        if shadowed && self.param_index(name).is_some() {
            return Err(self.err_hint(head, format!("'{name}' is a param, not a function"), shadow));
        }
        let hint = if let Some(h) = self.hidden_hint(name) {
            Some(h)
        } else if self.registry.top_forms().split(", ").any(|f| f == name) {
            // A top-level form's head inside an expression.
            let here = match name {
                "fn" => "; a function value here is (lambda (args) body)",
                "def" | "param" => "; a local name here is (let ((name value)) body)",
                _ => "",
            };
            Some(format!("'{name}' is a top-level form, written only at a file's top level (`ink spec --section \"Top-level forms\"`){here}"))
        } else if let Some(axis) = name.strip_prefix("rot-").filter(|a| ["x", "y", "z"].contains(a)) {
            Some(format!("'{name}' is a placement step, written inside (at ...); a shape turns by std's (rotate-{axis} f deg)"))
        } else if !self.std_imported && in_std(name) {
            Some(format!("'{name}' is in the standard library: add (import \"std\")"))
        } else {
            self.unknown_function_hint(name)
        };
        Err(self.err_hint(head, format!("unknown function '{name}'"), hint))
    }

    /// `unknown_function_hint` over the functions in scope.
    fn unknown_function_hint(&self, name: &str) -> Option<String> {
        let mut names: Vec<String> = self.fns.keys().cloned().collect();
        names.extend(SPECIAL.iter().map(|s| s.to_string()));
        names.extend(BUILTINS.iter().map(|(n, _)| n.to_string()));
        names.sort();
        unknown_function_hint(name, &names, self.std_imported)
    }
}
