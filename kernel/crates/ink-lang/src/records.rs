//! Per-sample output of batch commands (eval's points, ray's rays, a solve's
//! probes, a slice's loops), a table: the columns named once and
//! a row per sample; a mass's targets, records. `--keys` picks what each
//! record holds. A record's bulk (a slice's points, a mass's inertia) is
//! built only when `--keys` names it (`names`).

use crate::InkError;
use crate::command::Args;
use serde_json::{Value as J, json};

pub struct Records {
    /// Paths into a record: keys joined by dots (a key may hold dots
    /// itself, as a param's does). None: everything.
    keys: Option<Vec<String>>,
    /// Keys a command asked about (`wants`) and may have left out.
    asked: std::cell::RefCell<Vec<String>>,
}

impl Records {
    pub fn new(a: &Args) -> Records {
        let keys = a.opts.get("keys").map(|v| v.iter().flat_map(|s| s.split(',')).map(|k| k.trim().to_string()).filter(|k| !k.is_empty()).collect());
        Records { keys, asked: Default::default() }
    }

    /// Whether the output holds anything of the top-level key `key` (a
    /// command can skip computing what is not wanted).
    pub fn wants(&self, key: &str) -> bool {
        self.ask(key);
        self.keys.as_ref().is_none_or(|ks| ks.iter().any(|k| k == key || k.strip_prefix(key).is_some_and(|r| r.starts_with('.'))))
    }

    /// Whether `--keys` names `path` or a path inside it (not one
    /// around it): a part of a record left out unless asked for.
    pub fn names(&self, path: &str) -> bool {
        self.ask(path);
        self.keys.as_ref().is_some_and(|ks| ks.iter().any(|k| k == path || k.strip_prefix(path).is_some_and(|r| r.starts_with('.'))))
    }

    /// Whether `--keys` names `path` itself, whole.
    pub fn names_exactly(&self, path: &str) -> bool {
        self.keys.as_ref().is_some_and(|ks| ks.iter().any(|k| k == path))
    }

    /// `key` is one a record may hold: a typo's hint lists it.
    fn ask(&self, key: &str) {
        if !self.asked.borrow().iter().any(|k| k == key) {
            self.asked.borrow_mut().push(key.to_string());
        }
    }

    /// The records as a table: {columns, rows}. Keys in `keep` say how far
    /// a record can be trusted, or why it is empty (an undefined case the
    /// value passed through, a probe outside the grid): they stay whatever
    /// is picked, listed apart, as [row, value] pairs for the rows that
    /// have them.
    pub fn shape(&self, records: Vec<J>, keep: &[&str]) -> Result<J, InkError> {
        let paths = self.paths(&records, keep)?;
        // Columns: the leaves of what is picked, objects opened (a null
        // where others hold an object stands for its leaves, all null).
        let flat: Vec<Vec<(Vec<String>, J)>> = records
            .iter()
            .map(|r| {
                let mut cells = Vec::new();
                for p in &paths {
                    if let Some(x) = pick(r, p) {
                        leaves(&x, &mut Vec::new(), &mut cells);
                    }
                }
                cells
            })
            .collect();
        let mut cols = union(flat.iter().flatten().map(|(c, _)| c));
        for p in &paths {
            if records.iter().all(|r| pick(r, p).is_none()) {
                cols.push(vec![p.to_string()]);
            }
        }
        let opened = |c: &Vec<String>| cols.iter().any(|d| d.len() > c.len() && d.starts_with(c));
        let cols: Vec<Vec<String>> = cols.iter().filter(|c| !opened(c)).cloned().collect();
        let rows: Vec<J> = flat.iter().map(|cells| J::Array(cols.iter().map(|c| cells.iter().find(|(d, _)| d == c).map_or(J::Null, |(_, v)| v.clone())).collect())).collect();
        let mut out = json!({"columns": cols.iter().map(|c| c.join(".")).collect::<Vec<_>>(), "rows": rows});
        for &k in keep {
            let at: Vec<J> = records.iter().enumerate().filter_map(|(i, r)| Some(json!([i, r.get(k)?]))).collect();
            if !at.is_empty() {
                out[k] = J::Array(at);
            }
        }
        Ok(out)
    }

    /// The records, each holding its keys in `ids` (what it is of), what
    /// --keys picks (all of it, without) and its keys in `keep` (how far it
    /// can be trusted, or why it lacks what is picked): a mass's targets.
    pub fn picked(&self, records: Vec<J>, ids: &[&str], keep: &[&str]) -> Result<Vec<J>, InkError> {
        let paths: Vec<String> = self.paths(&records, keep)?.into_iter().filter(|p| !ids.contains(&p.as_str())).collect();
        Ok(records
            .iter()
            .map(|r| {
                let mut o = J::Object(ids.iter().filter_map(|&k| Some((k.to_string(), r.get(k)?.clone()))).collect());
                match &self.keys {
                    None => merge(&mut o, r.clone()),
                    Some(_) => paths.iter().filter_map(|p| pick(r, p)).for_each(|x| merge(&mut o, x)),
                }
                for &k in keep {
                    if let Some(v) = r.get(k) {
                        o[k] = v.clone();
                    }
                }
                o
            })
            .collect())
    }

    /// What is picked of each record, `keep` aside. A key no record has is
    /// a typo, unless the records say why they lack it (a probe outside the
    /// grid).
    fn paths(&self, records: &[J], keep: &[&str]) -> Result<Vec<String>, InkError> {
        let keys = || records.iter().filter_map(|r| r.as_object()).flat_map(|m| m.keys());
        let paths: Vec<String> = self.keys.clone().unwrap_or_else(|| union(keys())).into_iter().filter(|p| !keep.contains(&p.as_str())).collect();
        let whole: Vec<&J> = records.iter().filter(|r| keep.iter().all(|&k| r.get(k).is_none())).collect();
        for p in &paths {
            if !whole.is_empty() && whole.iter().all(|r| pick(r, p).is_none()) {
                let seen = union(keys().chain(self.asked.borrow().iter()));
                return Err(InkError {
                    hint: Some(format!("keys here: {}. A dot steps into one (KEY.SUBKEY)", seen.join(", "))),
                    ..InkError::msg(format!("--keys: no record has '{p}'"))
                });
            }
        }
        Ok(paths)
    }
}

/// The distinct items, in order of first appearance.
fn union<'a, T: PartialEq + Clone + 'a>(items: impl Iterator<Item = &'a T>) -> Vec<T> {
    let mut out: Vec<T> = Vec::new();
    for x in items {
        if !out.contains(x) {
            out.push(x.clone());
        }
    }
    out
}

/// The part of `v` at `path`, nested as in `v`; passing through a null
/// gives the null. None: `path` names nothing in `v`.
fn pick(v: &J, path: &str) -> Option<J> {
    let J::Object(m) = v else { return v.is_null().then_some(J::Null) };
    if let Some(x) = m.get(path) {
        return Some(json!({path: x}));
    }
    path.match_indices('.').find_map(|(i, _)| {
        let k = &path[..i];
        Some(json!({k: pick(m.get(k)?, &path[i + 1..])?}))
    })
}

fn merge(into: &mut J, from: J) {
    match (into, from) {
        (J::Object(a), J::Object(b)) => {
            for (k, v) in b {
                match a.get_mut(&k) {
                    Some(x) => merge(x, v),
                    None => {
                        a.insert(k, v);
                    }
                }
            }
        }
        (a, b) => *a = b,
    }
}

/// The leaves of `v` (arrays and scalars; objects are opened), by path.
fn leaves(v: &J, at: &mut Vec<String>, out: &mut Vec<(Vec<String>, J)>) {
    match v {
        J::Object(m) => {
            for (k, x) in m {
                at.push(k.clone());
                leaves(x, at, out);
                at.pop();
            }
        }
        _ => out.push((at.clone(), v.clone())),
    }
}

/// `--compact` printing: no whitespace, except a line per table row.
pub fn to_compact(v: &J) -> String {
    fn write(v: &J, s: &mut String) {
        match v {
            J::Object(m) => {
                s.push('{');
                for (i, (k, x)) in m.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    s.push_str(&J::String(k.clone()).to_string());
                    s.push(':');
                    match (k.as_str(), x) {
                        ("rows", J::Array(rows)) => {
                            s.push('[');
                            for (j, r) in rows.iter().enumerate() {
                                s.push_str(if j > 0 { ",\n" } else { "\n" });
                                s.push_str(&r.to_string());
                            }
                            s.push_str("\n]");
                        }
                        _ => write(x, s),
                    }
                }
                s.push('}');
            }
            J::Array(a) => {
                s.push('[');
                for (i, x) in a.iter().enumerate() {
                    if i > 0 {
                        s.push(',');
                    }
                    write(x, s);
                }
                s.push(']');
            }
            _ => s.push_str(&v.to_string()),
        }
    }
    let mut s = String::new();
    write(v, &mut s);
    s
}

/// Default printing: indented, but a value that fits on its line (a
/// vector, a matrix, a small object) printed on it, so a frame is a line
/// and not twelve.
pub fn to_pretty(v: &J) -> String {
    const WIDTH: usize = 100;
    fn flat(v: &J, s: &mut String) {
        match v {
            J::Object(m) => {
                s.push('{');
                for (i, (k, x)) in m.iter().enumerate() {
                    s.push_str(if i > 0 { ", " } else { "" });
                    s.push_str(&J::String(k.clone()).to_string());
                    s.push_str(": ");
                    flat(x, s);
                }
                s.push('}');
            }
            J::Array(a) => {
                s.push('[');
                for (i, x) in a.iter().enumerate() {
                    s.push_str(if i > 0 { ", " } else { "" });
                    flat(x, s);
                }
                s.push(']');
            }
            _ => s.push_str(&v.to_string()),
        }
    }
    // `v` starting at column `col`, its lines indented by `indent`.
    fn write(v: &J, indent: usize, col: usize, s: &mut String) {
        let mut one = String::new();
        flat(v, &mut one);
        let pad = |s: &mut String, n: usize| s.extend(std::iter::repeat_n(' ', n));
        match v {
            J::Object(m) if !m.is_empty() && col + one.len() > WIDTH => {
                s.push_str("{\n");
                for (i, (k, x)) in m.iter().enumerate() {
                    pad(s, indent + 2);
                    let key = J::String(k.clone()).to_string();
                    s.push_str(&key);
                    s.push_str(": ");
                    write(x, indent + 2, indent + 4 + key.len(), s);
                    s.push_str(if i + 1 < m.len() { ",\n" } else { "\n" });
                }
                pad(s, indent);
                s.push('}');
            }
            J::Array(a) if !a.is_empty() && col + one.len() > WIDTH => {
                s.push_str("[\n");
                for (i, x) in a.iter().enumerate() {
                    pad(s, indent + 2);
                    write(x, indent + 2, indent + 2, s);
                    s.push_str(if i + 1 < a.len() { ",\n" } else { "\n" });
                }
                pad(s, indent);
                s.push(']');
            }
            _ => s.push_str(&one),
        }
    }
    let mut s = String::new();
    write(v, 0, 0, &mut s);
    s
}
