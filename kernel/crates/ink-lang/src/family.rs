//! Families: one file that draws many parts, each a member with a name of
//! its own (docs/IR.md, Top-level forms: Families). A directory's
//! `families.json` holds the table of each family file beside it, by the
//! file's name less `.ink`: each member's name and record. `DIR/NAME.ink`,
//! where NAME is a member there, is that member: one lookup, in that one
//! file, which is also every tool's list of the directory's members.

use serde_json::{Map, Value as J};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// The file in a directory holding its families' tables.
pub const INDEX: &str = "families.json";

/// A table's members in order: each its name and its record (its groups'
/// records, then its own, merged: `merge`), and where each name is.
#[derive(Debug, Default)]
pub struct Members {
    pub list: Vec<(String, J)>,
    at: HashMap<String, usize>,
}

impl Members {
    /// The record of member `name`.
    pub fn get(&self, name: &str) -> Option<&J> {
        self.at.get(name).map(|&i| &self.list[i].1)
    }

    /// The first member: a family file loaded as itself is it.
    pub fn first(&self) -> Option<&(String, J)> {
        self.list.first()
    }
}

/// A directory's families (`INDEX`): each family file's name less `.ink`
/// and its members, in the file's order, and where each member is.
#[derive(Debug, Default)]
pub struct Families {
    pub list: Vec<(String, Members)>,
    at: HashMap<String, usize>,
}

impl Families {
    /// The members of family file `stem`.ink.
    pub fn family(&self, stem: &str) -> Option<&Members> {
        self.list.iter().find(|(f, _)| f == stem).map(|(_, m)| m)
    }

    /// Member `name`: its family's name and its record.
    pub fn member(&self, name: &str) -> Option<(&str, &J)> {
        let (f, m) = &self.list[*self.at.get(name)?];
        Some((f.as_str(), m.get(name)?))
    }

    /// Every member's name, family by family.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.list.iter().flat_map(|(_, m)| m.list.iter().map(|(k, _)| k.as_str()))
    }
}

/// What a member's or a family's name is: a file's name, less `.ink`, in
/// its directory (a tool reading families checks it alike).
const NAMES: &str = "not empty, no '/' or '\\', not hidden";
fn is_name(name: &str) -> bool {
    !name.is_empty() && !name.contains(['/', '\\']) && !name.starts_with('.')
}

/// A table's members (`Members`). A table is a group: an object of
/// `record` (what each of its members' records has), `members` (each
/// name's own record) and `groups` (groups within it, each its record over
/// this one's), all optional. A member's record is its groups' records,
/// outermost first, then its own, merged (`merge`). Its members come in
/// order: a group's own, then its groups'.
pub fn members(table: &J) -> Result<Members, String> {
    fn group(g: &J, at: &str, outer: &J, out: &mut Members) -> Result<(), String> {
        let J::Object(m) = g else { return Err(format!("{at}: a group is an object of \"record\", \"members\" and \"groups\"")) };
        if let Some(k) = m.keys().find(|k| !matches!(k.as_str(), "record" | "members" | "groups")) {
            return Err(format!("{at}: a group has \"record\", \"members\" and \"groups\", not \"{k}\""));
        }
        let record = match m.get("record") {
            None => outer.clone(),
            Some(r @ J::Object(_)) => merge(outer, r),
            Some(_) => return Err(format!("{at}/record: a record is an object")),
        };
        match m.get("members") {
            None => {}
            Some(J::Object(ms)) => {
                for (name, own) in ms {
                    if !matches!(own, J::Object(_)) {
                        return Err(format!("{at}/members/{name}: a member's record is an object"));
                    }
                    if !is_name(name) {
                        return Err(format!("{at}/members/{name}: a member's name is a file's name less .ink: {NAMES}"));
                    }
                    if out.at.insert(name.clone(), out.list.len()).is_some() {
                        return Err(format!("{at}/members/{name}: the table has member {name} twice"));
                    }
                    out.list.push((name.clone(), merge(&record, own)));
                }
            }
            Some(_) => return Err(format!("{at}/members: an object of each member's name and record")),
        }
        match m.get("groups") {
            None => Ok(()),
            Some(J::Array(gs)) => gs.iter().enumerate().try_for_each(|(i, g)| group(g, &format!("{at}/groups/{i}"), &record, out)),
            Some(_) => Err(format!("{at}/groups: a list of groups")),
        }
    }
    let mut out = Members::default();
    group(table, "#", &J::Object(Map::new()), &mut out)?;
    Ok(out)
}

/// A directory's families from its `INDEX`'s JSON: an object of each
/// family file's name less `.ink` and its table (`members`). A member's
/// name is once in the directory, and no family's.
pub fn families(index: &J) -> Result<Families, String> {
    let J::Object(fs) = index else { return Err("an object of each family file's name (less .ink) and its table".into()) };
    let mut out = Families::default();
    for (stem, table) in fs {
        if !is_name(stem) {
            return Err(format!("{stem}: a family's name is its file's name less .ink: {NAMES}"));
        }
        let m = members(table).map_err(|e| format!("{stem}{}", e.trim_start_matches('#')))?;
        for (name, _) in &m.list {
            if let Some(&i) = out.at.get(name) {
                return Err(format!("member {name} is in both {} and {stem}", out.list[i].0));
            }
            if fs.contains_key(name) {
                return Err(format!("{stem}: member {name} has a family's name"));
            }
            out.at.insert(name.clone(), out.list.len());
        }
        out.list.push((stem.clone(), m));
    }
    Ok(out)
}

/// `inner` over `outer`: two objects key by key (`outer`'s keys in its
/// order, merged where `inner` has them too, then `inner`'s others in its
/// order); else `inner`.
pub fn merge(outer: &J, inner: &J) -> J {
    match (outer, inner) {
        (J::Object(o), J::Object(i)) => {
            let mut m: Map<String, J> = o.iter().map(|(k, v)| (k.clone(), i.get(k).map_or_else(|| v.clone(), |x| merge(v, x)))).collect();
            m.extend(i.iter().filter(|(k, _)| !o.contains_key(*k)).map(|(k, v)| (k.clone(), v.clone())));
            J::Object(m)
        }
        _ => inner.clone(),
    }
}

/// The families of directory `dir` (its `INDEX`) and the index's content
/// key (`Doc::content_key`); None where it has none. Read once while the
/// file's size and time are unchanged (the process's: a catalog's members
/// are asked for many times, and an index is megabytes).
pub fn in_dir(dir: &Path) -> Result<Option<(Arc<Families>, u128)>, String> {
    type Kept = HashMap<PathBuf, ((u64, std::time::SystemTime), u128, Arc<Families>)>;
    static KEPT: Mutex<Option<Kept>> = Mutex::new(None);
    let path = dir.join(INDEX);
    let meta = match std::fs::metadata(&path) {
        Ok(m) if m.is_file() => m,
        Ok(_) => return Ok(None),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(format!("cannot read {}: {e}", crate::shown(&path))),
    };
    let stamp = (meta.len(), meta.modified().unwrap_or(std::time::UNIX_EPOCH));
    let canon = path.canonicalize().unwrap_or_else(|_| path.clone());
    if let Some((s, k, f)) = KEPT.lock().unwrap_or_else(|e| e.into_inner()).get_or_insert_with(HashMap::new).get(&canon)
        && *s == stamp
    {
        return Ok(Some((f.clone(), *k)));
    }
    let bytes = std::fs::read(&path).map_err(|e| format!("cannot read {}: {e}", crate::shown(&path)))?;
    let key = crate::Doc::content_key(&bytes);
    let j: J = serde_json::from_slice(&bytes).map_err(|e| format!("invalid JSON in {}: {e}", crate::shown(&path)))?;
    let f = Arc::new(families(&j).map_err(|e| format!("{}: {e}", crate::shown(&path)))?);
    KEPT.lock().unwrap_or_else(|e| e.into_inner()).get_or_insert_with(HashMap::new).insert(canon, (stamp, key, f.clone()));
    Ok(Some((f, key)))
}

/// A member, as `member_of` finds it.
#[derive(Debug, Clone)]
pub struct Member {
    /// Its family's file, beside it.
    pub family: PathBuf,
    /// Its name: its path's file name less `.ink`.
    pub name: String,
    /// The index that names it, and the index's content key.
    pub index: PathBuf,
    pub content: u128,
    pub record: J,
}

/// The member `DIR/NAME.ink` names: NAME, where `DIR/families.json` lists
/// it. None for a path no index names (a file, or nothing). An error where
/// the index cannot be read, or where a file has a member's name too
/// (`NAME.ink` or its record's `NAME.json`): a name is one part.
pub fn member_of(path: &Path) -> Result<Option<Member>, String> {
    if path.extension().is_none_or(|x| x != "ink") {
        return Ok(None);
    }
    let Some(name) = path.file_stem().and_then(|s| s.to_str()) else { return Ok(None) };
    let dir = path.parent().unwrap_or(Path::new(""));
    let Some((fams, content)) = in_dir(dir)? else { return Ok(None) };
    let Some((family, record)) = fams.member(name) else { return Ok(None) };
    let index = dir.join(INDEX);
    if let Some(f) = [path.to_path_buf(), path.with_extension("json")].into_iter().find(|f| f.exists()) {
        return Err(format!(
            "{} is a file, and {} names member {name} of {family}.ink too: remove the file or the member (one name, one part)",
            crate::shown(&f),
            crate::shown(&index)
        ));
    }
    Ok(Some(Member { family: dir.join(format!("{family}.ink")), name: name.to_string(), index, content, record: record.clone() }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    /// family_cases.json, which a tool reading families tests itself by too, so
    /// the two read a table alike: each index's members in order, each its
    /// family and its record (its keys in order too), or an error.
    #[test]
    fn a_member_is_its_groups_records_then_its_own() {
        let cases: J = serde_json::from_str(include_str!("family_cases.json")).unwrap();
        for case in cases.as_array().unwrap() {
            let why = &case["why"];
            match families(&case["index"]) {
                Err(e) if case["error"] == J::Bool(true) => assert!(!e.is_empty(), "{why}"),
                Err(e) => panic!("{why}: {e}"),
                Ok(_) if case["error"] == J::Bool(true) => panic!("{why}: read, not refused"),
                Ok(fs) => {
                    let got: Vec<J> = fs.list.iter().flat_map(|(f, m)| m.list.iter().map(move |(n, r)| json!([n, f, r]))).collect();
                    let want = case["members"].as_array().unwrap();
                    assert_eq!(serde_json::to_string(&got).unwrap(), serde_json::to_string(want).unwrap(), "{why}");
                    for w in want {
                        assert_eq!(fs.member(w[0].as_str().unwrap()).map(|m| m.0), w[1].as_str(), "{why}: {}", w[0]);
                    }
                }
            }
        }
        assert!(members(&json!({"members": {"A": {}}, "groups": [{"members": {"A": {}}}]})).unwrap_err().contains("twice"));
        assert!(families(&json!({"F": {"members": {"X": {}}}, "G": {"members": {"X": {}}}})).unwrap_err().contains("in both F and G"));
        assert!(families(&json!({"../f/f": {}})).unwrap_err().contains("a family's name"));
    }
}
