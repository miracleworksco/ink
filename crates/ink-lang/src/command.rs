//! Commands: what a package adds to a front end (`ink`'s command line and
//! its MCP tool), as its forms are what it adds to a document. A package
//! lists its commands (`Package::commands`); a front end finds them in the
//! registry it loads documents with (`Registry::command`), checks their
//! options, loads the command's file into a view and calls `run`, which
//! answers in JSON. Nothing in a front end is any one package's: a third
//! party's package plugs in as the standard distribution's do.
//!
//! What every command shares lives here: its arguments (`Args`), numbers
//! as answers give them (`n`, `v3`), in the language's units, records
//! shaped by `--keys` (`Records`), the inexact inputs an answer reads
//! (`inexact_reads`), and the call's time limit (`limit`).

use crate::InkError;
use crate::units::Quantity;
use crate::view::ModelView;
use crate::{Form, Inexact, Item};
use ink_core::V3;
use ink_core::graph::NodeId;
use serde_json::{Value as J, json};
use std::collections::HashMap;

pub use crate::records::{Records, to_compact, to_pretty};

/// A command a package adds: `ink NAME FILE [--option ...]`.
pub struct Command {
    pub name: &'static str,
    /// The options it reads besides those every command on a file reads
    /// (`--param`, `--t`, `--nominal`, `--compact`: they write nothing),
    /// each with what it writes: any other is refused before it runs.
    pub options: &'static [Opt],
    /// What `ink help` says of it: a line of what it does, and lines of its
    /// options (each indented 16 spaces).
    pub help: &'static str,
    /// What the MCP tool's description says of it, on one line: what it
    /// answers, and an example of its options as JSON (`{"name": "duct"}`).
    pub tool: &'static str,
    /// What it runs on, and how: its answer.
    pub run: Run,
}

/// An option a command reads, and the file it has the command write
/// besides its answer (which decides the MCP tool a call is served by).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Opt {
    pub name: &'static str,
    pub writes: Writes,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Writes {
    /// None.
    No,
    /// The file it names, when given.
    File,
    /// The file it names, else one of the command's own: the command
    /// writes a file either way.
    FileOrOwn,
}

impl Opt {
    /// An option that writes nothing.
    pub const fn reads(name: &'static str) -> Opt {
        Opt { name, writes: Writes::No }
    }
    /// One naming a file it writes, when given.
    pub const fn writes(name: &'static str) -> Opt {
        Opt { name, writes: Writes::File }
    }
    /// One naming the file it writes, else the command writes its own.
    pub const fn writes_or_own(name: &'static str) -> Opt {
        Opt { name, writes: Writes::FileOrOwn }
    }
}

/// What a command runs on.
pub enum Run {
    /// The model its file loads, at the params and moment the arguments
    /// give.
    Model(fn(&mut ModelView, &Args) -> Result<J, InkError>),
    /// No file: its words (`Args::file`, then `Args::extra`) and options.
    /// It loads what it reads with `registry` (the distribution's), and
    /// stops once `stop` says so (the call's time limit; an MCP call
    /// cancelled). `ink mcp` answers it in the server's own process where
    /// `here` says so (what it keeps between calls is kept: a catalog
    /// read), else in a process of its own, as every other command.
    Words { run: fn(&'static crate::Registry, &Args, &(dyn Fn() -> bool + Sync)) -> Result<J, InkError>, here: fn(&Args) -> bool },
}

/// A command's arguments: `ink CMD FILE --key value --flag ...`.
pub struct Args {
    pub cmd: String,
    pub file: Option<String>,
    pub opts: HashMap<String, Vec<String>>,
    /// Positional arguments after the file (none is taken).
    pub extra: Vec<String>,
}

impl Args {
    pub fn parse(v: &[String]) -> Args {
        if v.len() == 1 && matches!(v[0].as_str(), "--help" | "-h") {
            return Args { cmd: "help".into(), file: None, opts: HashMap::new(), extra: Vec::new() };
        }
        let mut pos = Vec::new();
        let mut opts: HashMap<String, Vec<String>> = HashMap::new();
        let mut i = 0;
        while i < v.len() {
            // After `--`, every word is positional.
            if v[i] == "--" {
                pos.extend(v[i + 1..].iter().cloned());
                break;
            }
            if let Some(k) = v[i].strip_prefix("--") {
                let (k, val) = match k.split_once('=') {
                    Some((k, val)) => (k.to_string(), val.to_string()),
                    None if i + 1 < v.len() && !v[i + 1].starts_with("--") => {
                        i += 1;
                        (k.to_string(), v[i].clone())
                    }
                    None => (k.to_string(), "true".to_string()),
                };
                opts.entry(k).or_default().push(val);
            } else {
                pos.push(v[i].clone());
            }
            i += 1;
        }
        let mut pos = pos.into_iter();
        Args { cmd: pos.next().unwrap_or_else(|| "help".into()), file: pos.next(), opts, extra: pos.collect() }
    }

    pub fn get(&self, k: &str) -> Option<&str> {
        self.opts.get(k).and_then(|v| v.last()).map(|s| s.as_str())
    }

    /// A flag: whether `--k` was given; it takes no value.
    pub fn flag(&self, k: &str) -> Result<bool, InkError> {
        match self.get(k) {
            None => Ok(false),
            Some("true") => Ok(true),
            Some(v) => Err(InkError::msg(format!("--{k} takes no value, got '{v}'"))),
        }
    }

    /// A number option. Those of a quantity (`--res`, `--t`) may say their
    /// unit (`--res 0.01in`), read in the language's.
    pub fn f64(&self, k: &str) -> Result<Option<f64>, InkError> {
        let q = match k {
            "res" | "max-res" | "accuracy" | "at" => Some(Quantity::Length),
            "t" => Some(Quantity::Time),
            _ => None,
        };
        self.get(k)
            .map(|s| match (s.parse::<f64>(), q) {
                (Ok(v), _) => Ok(v),
                (Err(_), Some(q)) => crate::sexp::number_of(s, q, &format!("--{k}")).map_err(InkError::msg),
                (Err(_), None) => Err(InkError::msg(format!("--{k} expects a number, got '{s}'"))),
            })
            .transpose()
    }

    /// The moment (`--t`, 0 by default): a finite number of seconds.
    pub fn moment(&self) -> Result<f64, InkError> {
        match self.f64("t")? {
            Some(t) if !t.is_finite() => Err(InkError::msg(format!("--t expects a moment, a finite number of seconds, got '{}'", self.get("t").unwrap_or_default()))),
            t => Ok(t.unwrap_or(0.0)),
        }
    }

    /// A cell size (`--res`, `--max-res`): a number above 0 and finite.
    /// No command takes cells of no size, and refining toward one never ends.
    pub fn cell(&self, k: &str) -> Result<Option<f64>, InkError> {
        match self.f64(k)? {
            Some(h) if !(h > 0.0 && h.is_finite()) => Err(InkError::msg(format!("--{k} expects a cell size above 0, got '{}'", self.get(k).unwrap_or_default()))),
            h => Ok(h),
        }
    }
}

/// The call's time limit, Unix seconds, if it has one: a front end that
/// limits a call (`ink mcp`) sets INK_DEADLINE in its process. Searches
/// stop before it (ink_kernel::stop); a command that runs long otherwise
/// stops in time to answer whole by it.
pub fn limit() -> Option<f64> {
    std::env::var("INK_DEADLINE").ok().and_then(|v| v.parse::<f64>().ok())
}

/// Err: the first solve the model's shapes rest on (Doc::needed) that does
/// not solve at this moment. A command that needs every value its shapes
/// read refuses the model then; one that solves again at each moment it
/// reads (effort), or refuses claims one by one (assert), does not ask.
pub fn solved(s: &ModelView) -> Result<(), InkError> {
    match s.unsolved.iter().find(|(i, _)| s.doc.needed(*i)) {
        Some((_, e)) => Err(e.clone()),
        None => Ok(()),
    }
}

/// The items `--target` names, by default the root, else every body.
pub fn items(s: &mut ModelView, a: &Args) -> Result<Vec<Item>, InkError> {
    Ok(match a.get("target") {
        Some(t) => s.doc.items(t)?,
        None => s.doc.default_items()?,
    })
}

/// The public name of the form of type `T` that option `--flag` names, or
/// of the document's only one (Doc::picked): the one name that finds it
/// again (an imported one's `ns.name`, never its own short name, which a
/// local form may share). An error says which option picks one.
pub fn form_named<T: Form>(s: &ModelView, kind: &str, a: &Args, flag: &str) -> Result<String, InkError> {
    s.doc.picked::<T>(kind, a.get(flag)).map(|e| e.name.clone()).map_err(|e| match s.doc.forms::<T>().next() {
        Some(_) if e.hint.is_none() => InkError { hint: Some(format!("--{flag} NAME picks one")), ..e },
        _ => e,
    })
}

/// A param's name as a file names it: plain words (an ASCII letter, then
/// letters, digits, '-' and '_') joined by '.', a namespace's, the grammar
/// the reader holds every param's name to (Doc::check_param_name). So a
/// name printed (an optimization's `apply`) sets that param and no other
/// through `--param`, with no quoting. A name a file cannot write (a
/// claim's range, "travel 0", which `--param "travel 0=v"` looks at) is
/// not one.
pub fn param_name(name: &str) -> bool {
    name.split('.').all(param_word)
}

/// One word of a param's name (`param_name`).
pub fn param_word(w: &str) -> bool {
    let mut c = w.chars();
    c.next().is_some_and(|c| c.is_ascii_alphabetic()) && c.all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
}

/// --param values: name=v sets a value (a declared tolerance band follows
/// it); name=lo:hi a range (the value at its middle), which assert and
/// interfere prove over; name=v:v pins a value (no band). A text param
/// takes the text as written (name=Grace); written as a JSON string
/// (name="Smith, Ada"), it may hold commas, colons and escapes.
pub struct ParamArgs {
    /// Each value, with the unit it was given in (d=0.25in), if any: a
    /// load reads the param as said in it (crate::Given), and the text it
    /// was written as, which a text param takes.
    pub values: Vec<(String, crate::Given)>,
    pub ranges: Vec<(String, (f64, f64))>,
    pub pinned: Vec<String>,
}

/// The `name=value` specs of `--param` options: each option split at
/// commas outside a quoted text ("...", escapes as JSON's).
fn param_specs(raw: &str) -> Vec<String> {
    let (mut out, mut cur, mut quoted, mut escaped) = (Vec::new(), String::new(), false, false);
    for c in raw.chars() {
        match c {
            _ if escaped => escaped = false,
            '\\' if quoted => escaped = true,
            '"' => quoted = !quoted,
            ',' if !quoted => {
                out.push(std::mem::take(&mut cur));
                continue;
            }
            _ => {}
        }
        cur.push(c);
    }
    out.push(cur);
    out
}

/// Why `--param k=v`, text no number param takes, is no number or range:
/// what a number param given it says.
pub fn not_a_number(k: &str, v: &str) -> String {
    let num = |v: &str| crate::sexp::number(v.trim()).map(|_| ());
    let e = match v.split_once(':') {
        Some((lo, hi)) => num(lo).and_then(|_| num(hi)).err(),
        None => num(v).err(),
    };
    format!("--param {k}: {}", e.unwrap_or_else(|| format!("'{v}' is no number")))
}

pub fn param_args(a: &Args) -> Result<ParamArgs, InkError> {
    let mut p = ParamArgs { values: Vec::new(), ranges: Vec::new(), pinned: Vec::new() };
    // A value may say its unit (d=0.25in), read in the language's.
    let num = |k: &str, v: &str| crate::sexp::number(v.trim()).map_err(|e| InkError::msg(format!("--param {k}: {e}")));
    for spec in a.opts.get("param").into_iter().flatten().flat_map(|s| param_specs(s)) {
        let (k, v) = spec.split_once('=').ok_or_else(|| InkError::msg(format!("--param expects name=value or name=lo:hi, got '{spec}'")))?;
        let k = k.trim().to_string();
        // One meaning per name: given twice, which one a command used would
        // depend on the order options are applied in.
        if p.values.iter().any(|(n, _)| *n == k) {
            return Err(InkError::msg(format!("--param {k} is given twice: give each param one value (name=v), range (name=lo:hi) or pin (name=v:v)")));
        }
        // Quoted: text, whole.
        if v.trim_start().starts_with('"') {
            let t: String =
                serde_json::from_str(v.trim()).map_err(|e| InkError::msg(format!("--param {k}: a quoted text is a JSON string, \"...\" with escapes as JSON's ({e}), got {v}")))?;
            p.values.push((k, crate::Given::text(&t)));
            continue;
        }
        let text = Some(std::sync::Arc::<str>::from(v));
        let range = v.split_once(':').map(|(lo, hi)| (num(&k, lo), num(&k, hi)));
        match range {
            None | Some((Err(_), _) | (_, Err(_))) => match num(&k, v) {
                Ok((value, unit)) => p.values.push((k.clone(), crate::Given { value, unit, text, case: None })),
                // No number: text, refused by a number param with this
                // error (Doc's param form).
                Err(_) => p.values.push((k.clone(), crate::Given { value: f64::NAN, unit: None, text, case: None })),
            },
            Some((Ok((lo, ul)), Ok((hi, uh)))) => {
                // A range's ends are one quantity, said as a sum's terms are
                // (units::kept).
                let unit = match (ul, uh) {
                    (Some(a), Some(b)) => Some(crate::units::kept(a, b).map_err(|m| InkError::msg(format!("--param {k}: a range's ends {m}")))?),
                    (a, b) => a.or(b),
                };
                if !(lo <= hi && lo.is_finite() && hi.is_finite()) {
                    return Err(InkError::msg(format!("--param {k}={lo}:{hi}: a range needs finite lo <= hi")));
                }
                p.values.push((k.clone(), crate::Given { value: ink_core::eval::span_mid((lo, hi)), unit, text: None, case: None }));
                if lo == hi {
                    p.pinned.push(k);
                } else {
                    p.ranges.push((k, (lo, hi)));
                }
            }
        }
    }
    Ok(p)
}

/// --nominal: claims at the params' values, their declared tolerances
/// dropped.
pub fn nominal(a: &Args) -> Result<bool, InkError> {
    match a.get("nominal") {
        None => Ok(false),
        Some("true") => Ok(true),
        Some(v) => Err(InkError::msg(format!("--nominal takes no value, got '{v}'"))),
    }
}

/// The scope assert and interfere answer claims over (ClaimScope): the motion,
/// --t-span, else unless a moment is given each claim's file's (t-span), the
/// run file's first; and the
/// params' spans: declared tolerances (around the values, after
/// --param name=v), unless --nominal; a pinned param has none; an explicit
/// range replaces a declared one. A (ranged lo hi) is always ranged
/// (ModelView::scope_spans). The one setup of a claim's question from a
/// run's options: `ink assert`'s, and an optimization's of the design it
/// returns (ink-design verify), on the view of that design.
pub fn claim_scope(s: &ModelView, a: &Args) -> Result<crate::view::ClaimScope, InkError> {
    let (params, nominal) = (param_args(a)?, nominal(a)?);
    let given = t_span(a)?;
    let by_file = given.is_none() && a.get("t").is_none();
    let motion = given.or(s.doc.t_span.filter(|_| by_file));
    let mut spans: Vec<(usize, (f64, f64))> = if nominal { Vec::new() } else { s.declared_spans()? };
    let index = |k: &str| s.doc.param_index(k).ok_or_else(|| InkError::msg(format!("--param: unknown param '{k}'")));
    for k in &params.pinned {
        let i = index(k)?;
        spans.retain(|x| x.0 != i);
    }
    for (k, r) in &params.ranges {
        let i = index(k)?;
        spans.retain(|x| x.0 != i);
        spans.push((i, *r));
    }
    Ok(crate::view::ClaimScope { motion, by_file, spans })
}

/// --t-span a:b (claims that read t are checked at every moment in it).
pub fn t_span(a: &Args) -> Result<Option<(f64, f64)>, InkError> {
    let Some(raw) = a.get("t-span") else { return Ok(None) };
    let (x, y) = raw.split_once(':').ok_or_else(|| InkError::msg("--t-span expects a:b"))?;
    let time = |v: &str| crate::sexp::number_of(v.trim(), Quantity::Time, "--t-span").map_err(|e| InkError::msg(format!("--t-span a:b needs times: {e}")));
    let (x, y) = (time(x)?, time(y)?);
    if !(x < y && x.is_finite() && y.is_finite()) {
        return Err(InkError::msg("--t-span needs finite a < b"));
    }
    Ok(Some((x, y)))
}

/// The JSON value of `--key`: given inline, or, since large probe sets
/// exceed the argument limit, `@file` or `-` (stdin). None: not given.
pub fn json_arg(a: &Args, key: &str) -> Result<Option<J>, InkError> {
    let Some(raw) = a.get(key) else { return Ok(None) };
    let text = match raw {
        "-" => std::io::read_to_string(std::io::stdin()).map_err(|e| InkError::msg(format!("--{key} -: {e}")))?,
        _ => match raw.strip_prefix('@') {
            Some(path) => std::fs::read_to_string(path).map_err(|e| InkError::msg(format!("--{key} @{path}: {e}")))?,
            None => raw.to_string(),
        },
    };
    serde_json::from_str(&text).map(Some).map_err(|e| InkError::msg(format!("--{key} is not valid JSON: {e}")))
}

/// [x, y, z] as a point.
pub fn json_v3(p: &J) -> Option<V3> {
    let a = p.as_array().filter(|a| a.len() == 3)?;
    Some([a[0].as_f64()?, a[1].as_f64()?, a[2].as_f64()?])
}

/// `--points`: one point [x, y, z] or a list of them.
pub fn parse_points(a: &Args) -> Result<Vec<V3>, InkError> {
    let v = json_arg(a, "points")?.ok_or_else(|| InkError::msg("eval needs --points '[[x,y,z],...]'"))?;
    if let Some(p) = json_v3(&v) {
        return Ok(vec![p]);
    }
    v.as_array().and_then(|arr| arr.iter().map(json_v3).collect::<Option<Vec<_>>>()).ok_or_else(|| InkError::msg("--points must be [x,y,z] or [[x,y,z],...]"))
}

/// A number in an answer, written exactly: rounded, a bound (an
/// enclosure, a box, a band, an accuracy) could move past what it bounds.
/// Only `measured` rounds.
pub fn n(v: f64) -> J {
    if v.is_finite() { json!(v) } else { J::Null }
}

pub fn v3(p: V3) -> J {
    json!([n(p[0]), n(p[1]), n(p[2])])
}

/// An interval, `n` at both ends.
pub fn bounds(lo: f64, hi: f64) -> J {
    json!([n(lo), n(hi)])
}

/// A measured value (a sample, an estimate), rounded to 10 significant
/// digits for reading. Never a bound or an accuracy: those are `n`.
pub fn measured(v: f64) -> J {
    if v.is_finite() { json!(format!("{v:.9e}").parse::<f64>().unwrap()) } else { J::Null }
}

/// `more`'s keys added to `out` (both objects).
pub fn extend(out: &mut J, more: J) {
    if let (Some(o), J::Object(m)) = (out.as_object_mut(), more) {
        o.extend(m);
    }
}

/// Say in `out` which inexact inputs (Doc::inexact: an imported mesh, a
/// baked field, a solved value not certified) the fields `roots` an answer
/// read are: the answer is measured.
pub fn inexact_reads(s: &ModelView, out: &mut J, roots: &[NodeId]) {
    let x = s.doc.inexact(roots, &s.certified());
    if x.any() && out.is_object() {
        let names: Vec<String> = s.doc.params.iter().map(|p| p.name.clone()).collect();
        inexact_fields(out, &x, &names);
        // Any proof in the output holds for the inputs as evaluated (the
        // samples, the values the solve found), not for their sources: the
        // answer is measured.
        out["measured"] = json!(true);
    }
    outside_band_reads(s, out, roots);
}

/// Say in `out` which toleranced params the fields `roots` read at their
/// declared value where their band excludes it (Param::outside_band): the
/// answer measured a size no real part has (`outside_band`: each one's
/// value and band). One given by --param is a size asked for.
fn outside_band_reads(s: &ModelView, out: &mut J, roots: &[NodeId]) {
    let mut ks = s.doc.params_read(roots);
    ks.sort_unstable();
    ks.dedup();
    let o: serde_json::Map<String, J> = ks
        .into_iter()
        .filter_map(|k| {
            let p = &s.doc.params[k];
            let at_declared = s.ctx.params.get(k) == Some(&p.default) && !s.overrides.iter().any(|(name, _)| *name == p.name);
            p.outside_band().filter(|_| at_declared).and(p.band(p.default)).map(|(a, b)| (p.name.clone(), json!({"value": n(p.default), "band": bounds(a, b)})))
        })
        .collect();
    if !o.is_empty() && out.is_object() {
        out["outside_band"] = J::Object(o);
    }
}

/// The inexact inputs an answer reads (Doc::inexact), which make it an
/// estimate: the sampled fields, and the solved values not certified as
/// roots (it holds for the values the solve found), with why any of them
/// fall short of their accuracy (`unsettled`).
pub fn inexact_fields(r: &mut J, x: &Inexact, names: &[String]) {
    if !x.samples.is_empty() {
        r["samples"] = sources_json(&x.samples);
    }
    if !x.solved.is_empty() {
        r["solved_inputs"] = json!(x.solved.iter().map(|&k| &names[k as usize]).collect::<Vec<_>>());
    }
    if !x.unsettled.is_empty() {
        r["unsettled"] = J::Object(x.unsettled.iter().map(|(k, why)| (names[*k as usize].clone(), json!(why))).collect());
    }
    if !x.notes.is_empty() {
        r["notes"] = json!(x.notes);
    }
}

/// Sampled fields an answer reads: what each was sampled from, where,
/// at what spacing, and how far it may be from it (a sound bound, and the
/// deviation measured at probe points).
pub fn sources_json(sources: &[std::sync::Arc<ink_core::grid::Source>]) -> J {
    J::Array(
        sources
            .iter()
            .map(|s| {
                let mut o = json!({"source": s.what, "kind": s.kind, "at": s.at, "res": n(s.res),
                    "error_bound": s.band.map_or(J::Null, |b| n(b.err))});
                // The largest |samples - source| found at probe points: an
                // estimate, not a bound.
                if let Some((m, k)) = s.measured {
                    o["measured_error"] = json!({"max": n(m), "points": k});
                }
                if let Some(note) = &s.note {
                    o["note"] = json!(note);
                }
                o
            })
            .collect(),
    )
}
