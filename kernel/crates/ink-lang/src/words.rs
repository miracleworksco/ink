//! The words the language keeps (its built-in forms, operators and
//! constants), and the hints an error gives about a word: a typo's nearest
//! name, a word from another language or CAD tool.

use super::*;

#[derive(Clone, Copy)]
pub(super) enum Builtin {
    Un(UnOp),
    Bin(BinOp),
    Variadic(BinOp),
    Minus,
    Select,
}

/// The arithmetic forms, by name.
pub(super) const BUILTINS: &[(&str, Builtin)] = {
    use Builtin::*;
    &[
        ("+", Variadic(BinOp::Add)),
        ("*", Variadic(BinOp::Mul)),
        ("min", Variadic(BinOp::Min)),
        ("max", Variadic(BinOp::Max)),
        ("emax", Variadic(BinOp::Emax)),
        ("-", Minus),
        ("/", Bin(BinOp::Div)),
        ("atan2", Bin(BinOp::Atan2)),
        ("pow", Bin(BinOp::Pow)),
        ("neg", Un(UnOp::Neg)),
        ("abs", Un(UnOp::Abs)),
        ("sqrt", Un(UnOp::Sqrt)),
        ("square", Un(UnOp::Square)),
        ("exp", Un(UnOp::Exp)),
        ("ln", Un(UnOp::Ln)),
        ("sin", Un(UnOp::Sin)),
        ("cos", Un(UnOp::Cos)),
        ("tan", Un(UnOp::Tan)),
        ("floor", Un(UnOp::Floor)),
        ("deg->rad", Un(UnOp::Rad)),
        ("rad->deg", Un(UnOp::Deg)),
        ("select", Select),
    ]
};

pub(super) fn builtin(name: &str) -> Option<Builtin> {
    BUILTINS.iter().find(|(n, _)| *n == name).map(|&(_, b)| b)
}

pub(super) const SPECIAL: &[&str] = &[
    "let",
    "vec",
    "unit",
    "sin-cos-deg",
    "nth",
    "data",
    "warp",
    "estimate-when",
    "shape-of",
    "lambda",
    "map",
    "reduce",
    "apply",
    "range",
    "len",
    "concat",
    "char-codes",
    "slice",
    "if",
    "read-json",
    "get",
    "at",
    "axes",
    "frame-of",
    "place",
    "fold",
    "flatten",
    "bounds",
    "tag",
    "real",
    "lipschitz",
    "corner-reach",
    "grid-bake",
    "import-stl",
    "redistance",
    "study-field",
    "diff",
    "material",
    "by-body",
    "attr",
    "attr-of",
    "met-by",
    "check",
    "offset-shape",
    "error",
    "profile-reach",
    "polygon",
    "revolve-z",
    "probe",
    "relation",
    "ranged",
    "font",
    "glyph-2d",
    "glyph-metrics",
    "glyph-run",
    "glyph-reach",
];
/// The forms that are functions: each argument an expression, evaluated
/// once (no names, clauses or literal text), so they are values like a
/// user `fn`: `(map polygon outlines)`. The other forms are syntax.
pub(super) const FORM_FUNCTIONS: &[&str] = &[
    "vec",
    "unit",
    "nth",
    "map",
    "reduce",
    "fold",
    "flatten",
    "apply",
    "range",
    "len",
    "concat",
    "char-codes",
    "slice",
    "get",
    "warp",
    "lipschitz",
    "corner-reach",
    "offset-shape",
    "bounds",
    "profile-reach",
    "polygon",
    "revolve-z",
    "probe",
    "ranged",
    "glyph-2d",
    "glyph-metrics",
    "glyph-run",
    "glyph-reach",
];
/// Builtins beyond the kernel ops: comparisons yielding 1 or 0 (as fields).
pub(super) const COMPARISONS: &[&str] = &["<", "<=", ">", ">=", "=="];
pub(super) const CONSTANTS: &[&str] = &["x", "y", "z", "t", "pi", "inf"];

/// Whether a value is, or holds, a point.
pub(super) fn has_point(v: &Value) -> bool {
    match v {
        Value::Point(..) => true,
        Value::Vec(items) => items.iter().any(has_point),
        _ => false,
    }
}

pub(super) fn is_reserved(name: &str) -> bool {
    builtin(name).is_some() || SPECIAL.contains(&name) || CONSTANTS.contains(&name) || COMPARISONS.contains(&name)
}

/// The forms built in that a word names (the operators' symbols and the
/// constants aside), sorted: the kernel's operations and the language's
/// own forms, which no std definition takes (a file's own may: Doc::shadows).
pub fn builtin_words() -> Vec<&'static str> {
    let mut words: Vec<&str> = BUILTINS.iter().map(|(n, _)| *n).chain(SPECIAL.iter().copied()).filter(|n| n.starts_with(char::is_alphabetic)).collect();
    words.sort();
    words.dedup();
    words
}

/// An error's hint naming what is_reserved keeps, words first (the
/// operators' symbols aside): every name refused is in it.
pub(super) fn reserved_hint() -> String {
    format!("the built-in names: {}; {}", CONSTANTS.join(" "), builtin_words().join(" "))
}

/// `(fn NAME ARGS ...)` for a function written as a def: `(def (f x)
/// body)`, `(def f (x) body)` or `(def f (fn (x) body))`.
pub(super) fn fn_written_as_def(args: &[Sexp]) -> Option<String> {
    let names = |e: &Sexp| e.list().filter(|l| l.iter().all(|a| a.sym().is_some())).map(|l| l.iter().map(|a| a.to_string()).collect::<Vec<_>>().join(" "));
    let (name, ps) = match args {
        [head, _] if head.list().is_some_and(|l| l.len() > 1) => {
            let all = names(head)?;
            let (n, ps) = all.split_once(' ')?;
            (n.to_string(), ps.to_string())
        }
        [n, ps, _] => (n.sym()?.to_string(), names(ps)?),
        [n, e] => match e.list()? {
            [h, ps, _] if h.sym() == Some("fn") => (n.sym()?.to_string(), names(ps)?),
            _ => return None,
        },
        _ => return None,
    };
    Some(format!("(fn {name} ({ps}) body)"))
}

/// Typos between two names: letters inserted, deleted, changed, or two
/// swapped ('stduy' is one from 'study').
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let (mut before, mut prev): (Vec<usize>, Vec<usize>) = (Vec::new(), (0..=b.len()).collect());
    for i in 1..=a.len() {
        let mut cur = vec![i; b.len() + 1];
        for j in 1..=b.len() {
            let c = if a[i - 1] == b[j - 1] { 0 } else { 1 };
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + c);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                cur[j] = cur[j].min(before[j - 2] + 1);
            }
        }
        before = std::mem::replace(&mut prev, cur);
    }
    prev[b.len()]
}

/// std's forms, parsed once.
pub(super) fn std_forms() -> &'static std::sync::Arc<[Sexp]> {
    static STD: std::sync::OnceLock<std::sync::Arc<[Sexp]>> = std::sync::OnceLock::new();
    STD.get_or_init(|| parse(STD_SRC).expect("std parses").into())
}

/// Whether std defines a function `name`.
pub(super) fn in_std(name: &str) -> bool {
    static FNS: std::sync::OnceLock<std::collections::HashSet<String>> = std::sync::OnceLock::new();
    FNS.get_or_init(|| {
        std_forms()
            .iter()
            .filter_map(|f| match f.list() {
                Some([h, n, ..]) if h.sym() == Some("fn") => n.sym().map(String::from),
                _ => None,
            })
            .collect()
    })
    .contains(name)
}

/// The candidate nearest `name` by edit distance, as a hint, if a typo
/// could make one of the other: at most a third of the name's letters
/// differ (one, in a short name: 'hole' is not 'fold').
pub fn closest<'a>(name: &str, candidates: impl Iterator<Item = &'a str>) -> Option<String> {
    near(name, candidates).first().map(|c| format!("did you mean '{c}'?"))
}

/// What to pick instead of `name`, not one of `names` (`what`, e.g.
/// "its instances"): the nearest by `closest`, else all of them.
pub fn choices<'a>(name: &str, what: &str, names: impl Iterator<Item = &'a str> + Clone) -> String {
    let all = names.clone().collect::<Vec<_>>().join(", ");
    closest(name, names).unwrap_or_else(|| format!("{what}: {}", if all.is_empty() { "none" } else { &all }))
}

/// The candidates a typo could make of `name` (closest), nearest first:
/// each keeps a letter of it in place ('w' is no typo of 't').
pub(super) fn near<'a>(name: &str, candidates: impl Iterator<Item = &'a str>) -> Vec<&'a str> {
    let n = name.chars().count();
    let most = (n / 3).max(1);
    // An edit changes the length by one at most: a candidate longer or
    // shorter by more than `most` is further (not measured).
    let mut c: Vec<(usize, &str)> = candidates
        .map(|c| (c, c.chars().count()))
        .filter(|(_, k)| k.abs_diff(n) <= most)
        .filter_map(|(c, k)| Some((edit_distance(name, c), c)).filter(|&(d, _)| d <= most && d < n.max(k)))
        .collect();
    c.sort();
    c.into_iter().map(|(_, c)| c).collect()
}

/// An error's hint where what it is at is, or directly holds, a word
/// another language opens a comment with (`//`, `/*`): read as a name, the
/// comment is what is refused. Such a word is a name elsewhere (data's
/// paths and globs), so only an error says it.
pub(super) fn comment_hint(at: &Sexp) -> Option<String> {
    let opens = |e: &Sexp| e.sym().is_some_and(|w| w.starts_with("//") || w.starts_with("/*"));
    (opens(at) || at.list().is_some_and(|l| l.iter().any(opens))).then(|| "a comment runs from ';' to the line's end".to_string())
}

/// Words agents bring from other CAD tools and languages, with what ink
/// calls the same thing, and how it is written where that is not a call
/// by another name.
pub(super) const CAD_WORDS: &[(&str, &[&str], &str)] = &[
    ("list", &["vec"], "a list is (vec a b ...)"),
    (
        "for",
        &["map", "reduce"],
        "in an expression, map over a list: (reduce min (map (lambda (i) (move s (* 10 i) 0 0)) (range n))); (for ...) repeats forms, at top level and among clauses",
    ),
    ("hole", &["subtract", "cylinder-z"], "a hole is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("bore", &["subtract", "cylinder-z"], "a bore is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("drill", &["subtract", "cylinder-z"], "a drilled hole is the part minus a cylinder: (subtract part (move (cylinder-z r h) x y z))"),
    ("cut", &["subtract"], ""),
    ("difference", &["subtract"], ""),
    ("minus", &["subtract"], ""),
    ("fillet", &["fillet-union", "corner-fillet", "smooth-union", "rounded-box"], ""),
    ("round", &["floor", "fillet-union", "corner-fillet", "smooth-union", "rounded-box"], "a number rounds as (floor (+ v 0.5)); an edge, with a fillet"),
    ("chamfer", &["chamfer-union", "chamfer-intersect"], ""),
    ("bevel", &["chamfer-union", "chamfer-intersect"], ""),
    ("extrude", &["extrude-z", "sweep", "loft-z"], ""),
    ("prism", &["extrude-z"], ""),
    ("pad", &["extrude-z"], ""),
    ("revolve", &["revolve-z"], ""),
    ("lathe", &["revolve-z"], ""),
    ("translate", &["move"], ""),
    ("rotate", &["rotate-x", "rotate-y", "rotate-z"], ""),
    ("mirror", &["mirror-x", "mirror-y", "mirror-z"], ""),
    ("cube", &["box"], ""),
    ("cuboid", &["box"], ""),
    ("block", &["box", "rounded-box"], ""),
    ("ball", &["sphere"], ""),
    ("union", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("join", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("combine", &["min"], "a union is the min of its shapes: (min a b ...)"),
    ("intersect", &["max"], "an intersection is the max of its shapes: (max a b ...)"),
    ("intersection", &["max"], "an intersection is the max of its shapes: (max a b ...)"),
    ("or", &["max"], "a condition is 1 or 0: (max a b) is a or b, (min a b) a and b"),
    ("and", &["min"], "a condition is 1 or 0: (min a b) is a and b, (max a b) a or b"),
    ("array", &["pattern-linear", "pattern-polar", "repeat-linear", "repeat-x"], ""),
    ("hollow", &["shell"], ""),
    ("thicken", &["shell", "offset-shape"], ""),
    ("offset", &["offset-shape"], "a shape grown by d is (offset-shape f d), shrunk for d < 0"),
    ("plane", &["half-space"], ""),
    ("ring", &["ring-2d", "torus-z"], ""),
    ("annulus", &["ring-2d"], ""),
    ("hex", &["hexagon-2d"], ""),
    ("hexagon", &["hexagon-2d"], ""),
    ("cone", &["cone-z"], ""),
    ("frustum", &["cone-z"], ""),
    ("taper", &["cone-z"], ""),
    ("countersink", &["cone-z"], ""),
    ("rectangle", &["rect-2d"], ""),
    ("screw", &["thread-iso"], ""),
    ("spring", &["helix-points", "sweep"], ""),
];

/// A hint for a call of an unknown function `name`: what std calls the
/// thing a CAD word names (CAD_WORDS: the whole name or a part of it
/// between dashes), the names one part of which a typo could make of it
/// ('cylnder': 'cylinder-z'), and the names a typo could make of it; `names` the
/// functions in scope, std's included where it is imported.
pub(super) fn unknown_function_hint(name: &str, names: &[String], std_imported: bool) -> Option<String> {
    let known = |c: &str| names.iter().any(|n| n == c) || in_std(c);
    let mut pick: Vec<String> = Vec::new();
    let mut add = |c: &str| {
        if !pick.iter().any(|p| p == c) {
            pick.push(c.to_string());
        }
    };
    let words: Vec<&(&str, &[&str], &str)> = CAD_WORDS.iter().filter(|(w, ..)| name == *w || name.split('-').any(|s| s == *w || s == format!("{w}s"))).collect();
    for (_, std_names, _) in &words {
        std_names.iter().filter(|c| known(c)).for_each(|c| add(c));
    }
    if name.len() >= 3 {
        names.iter().filter(|c| c.contains('-') && !near(name, c.split('-')).is_empty()).for_each(|c| add(c));
    }
    near(name, names.iter().map(String::as_str)).into_iter().take(3).for_each(&mut add);
    if pick.is_empty() {
        return None;
    }
    pick.truncate(6);
    let quoted: Vec<String> = pick.iter().map(|c| format!("'{c}'")).collect();
    let mut hint = match quoted.as_slice() {
        [one] => format!("did you mean {one}?"),
        all => format!("did you mean one of {}?", all.join(", ")),
    };
    if let Some((.., how)) = words.iter().find(|(.., how)| !how.is_empty()) {
        hint += &format!(" ({how})");
    }
    if !std_imported && pick.iter().any(|c| in_std(c) && !names.contains(c)) {
        hint += " (std's: add (import \"std\"))";
    }
    Some(hint + " `ink std --section NAME` shows one.")
}
