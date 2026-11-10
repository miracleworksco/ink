//! `ink`: command-line interface. Every command on a file prints JSON. A library, so
//! a distribution with more packages builds its own `ink` on it
//! ([`main`] with its registry): this crate's binary is the kernel's.

use ink::eval::{eval_dual, eval_points_par};
use ink::lang::command::{Args, Command, Opt, Records, Run, bounds, extend, inexact_fields, json_arg, json_v3, measured, n, parse_points, sources_json, to_compact, to_pretty, v3};
use ink::lang::{Attrs, Item, Reach, Registry, Value, Xf};
use ink::mesh::{self, Mesh};
use ink::query::{self, max_extent};
use ink::rigid::Rigid;
use ink::tape::Tape;
use ink::view::ModelView;
use ink::{InkError, V3, export, render};
use serde_json::{Value as J, json};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

mod cli;
mod mcp;
use cli::*;

/// Meshing for seeing and its GLB, for a package that writes a GLB of its
/// own on them (a platform's, for its viewer): the plan of what to mesh
/// (`plan_meshes`), which it may change before the meshes are made
/// (`make_meshes`), and `ink mesh`'s answer.
pub mod meshes {
    pub use crate::cli::queries::{
        Budget, MeshHow, Meshed, Plan, Quality, REPORTED, SEEN, drawn, glb_attrs, looks, make_meshes, mesh_answer, motion, paint_tapes, plan_meshes, reported, span_times, tape_of,
    };
    pub use crate::cli::session::{CliView, Session, union_bounds};
}

/// Each core command: its name, what it answers (`ink help` and the MCP
/// tool's description both print it), the options it reads besides
/// those every command on a file reads (MODEL_OPTIONS, PROCESS_OPTIONS),
/// and the file it writes (so the MCP tool it is served under).
/// `ink help` and the tool's description are printed from this table and
/// OPTION_HELP, so they cannot disagree.
pub const COMMANDS: &[(&str, &str, &[Opt])] = &[
    (
        "check",
        "Load the model: its own params, defs, fns, bodies, assemblies with frames, claims and warnings, and each import in brief.",
        &[Opt::reads("load"), Opt::reads("imports"), Opt::reads("name")],
    ),
    (
        "eval",
        "The field's value, inside or not, gradient, and rates of change with t and each param, at points. Attributes where keys names them (attrs.spec).",
        &[Opt::reads("target"), Opt::reads("points"), Opt::reads("keys")],
    ),
    ("bounds", "Each target's proven bounding box.", &[Opt::reads("target")]),
    (
        "mass",
        "Mass, volume and centre of mass, estimated to 2%. Density is the target's, its material's or its published mass's; with none, mass and total mass are null. Inertia and principal axes where keys names them.",
        &[Opt::reads("target"), Opt::reads("density"), Opt::reads("res"), Opt::reads("keys")],
    ),
    (
        "slice",
        "Cross-section loops in a plane: kind, area, perimeter, bounds. Their points where keys names them.",
        &[Opt::reads("target"), Opt::reads("axis"), Opt::reads("at"), Opt::reads("res"), Opt::reads("accuracy"), Opt::reads("keys"), Opt::writes("out")],
    ),
    (
        "mesh",
        "Write STL or 3MF (to print) or GLB (to view). Reports triangles, cells where a feature under res went missing, and (STL, 3MF) sampled deviation from the surface.",
        &[Opt::reads("target"), Opt::writes_or_own("out"), Opt::reads("res"), Opt::reads("max-res"), Opt::reads("accuracy")],
    ),
    ("components", "Count cavities by proof, as (cavities f < n) does, and list a mesh's shells.", &[Opt::reads("target"), Opt::reads("res")]),
    (
        "interfere",
        "Overlaps between an assembly's instances outside contact zones, at rest or over the motion.",
        &[Opt::reads("target"), Opt::reads("res"), Opt::reads("accuracy"), Opt::reads("t-span")],
    ),
    (
        "assert",
        "Check the file's claims and joint limits: each holds, fails (with a witness point where one shows it) or is unknown (with an undecided_cell, and a hint where it can). accuracy says how far off each answer may be; null means unknown or estimated. Exits 3 if any fails or is unknown.",
        &[Opt::reads("name"), Opt::reads("own"), Opt::reads("res"), Opt::reads("t-span"), Opt::reads("keys")],
    ),
    ("ray", "Where rays first hit each target: point, normal, name.", &[Opt::reads("target"), Opt::reads("rays"), Opt::reads("from"), Opt::reads("dir"), Opt::reads("keys")]),
    ("distance", "Minimum distance between two targets: a proven lower bound and a witness pair of points.", &[Opt::reads("a"), Opt::reads("b"), Opt::reads("res")]),
    (
        "render",
        "Labelled PNG from iso, top, front and right, or one view, with the res each target was meshed at.",
        &[Opt::reads("target"), Opt::writes_or_own("out"), Opt::reads("view"), Opt::reads("cut"), Opt::reads("size"), Opt::reads("res")],
    ),
    ("spec", "The language reference: its headings, a section, or all of it.", &[Opt::reads("contents"), Opt::reads("section"), Opt::reads("all")]),
    ("std", "The standard library: sections and signatures, a section or named definitions, or all of it.", &[Opt::reads("contents"), Opt::reads("section"), Opt::reads("all")]),
    ("mcp", "Serve these commands over MCP (stdio): claude mcp add ink -- ink mcp", &[]),
    ("help", "This text.", &[]),
];

/// Each option: its name, an example of its value as the MCP tool takes it
/// (JSON), and what it does.
pub const OPTION_HELP: &[(&str, &str, &str)] = &[
    (
        "param",
        r#"["w=40", "h=5"]"#,
        "Set a param (repeatable). assert and interfere also take a range, w=38:42. A text param takes text, name=Grace; one with a comma is a JSON string, name=\"Lovelace, Ada\".",
    ),
    ("t", "1.5", "The moment, s (default 0)."),
    ("nominal", "true", "assert, interfere: ignore declared tolerances."),
    ("compact", "true", "No whitespace (the MCP tool's default)."),
    ("target", r#""root/arm""#, "A body, assembly, def or instance path. Default root, else every body (render: else the first assembly no other places)."),
    ("res", "0.5", "Cell size, mm. Commands that sample use it as given (default from the target's size); assert, interfere and distance refine undecided cells down to it."),
    ("keys", r#""value,inside""#, "Keep only these keys of each record (a dot steps into an object: nearest.point), including ones left out by default."),
    ("load", "true", "Only what loading decides: skips body bounds, surface sampling, published-mass checks and solves over tolerance bands."),
    ("imports", "true", "The imports in full."),
    ("points", "[[0, 0, 0]]", "Points [[x,y,z],...], or @file.json, or - (stdin)."),
    ("density", r#""1.24g/cm3""#, "A density (a number, g/mm3 unless it gives a unit) or the attribute naming one."),
    ("axis", r#""z""#, "The plane's normal: x, y or z."),
    ("at", "5", "The plane's position along axis, mm (default the middle)."),
    (
        "accuracy",
        "0.02",
        "How far simplifying may move the outline or mesh from the samples, mm (default: slice res / 2, STL and 3MF res / 30, GLB res / 4; 0 keeps it uniform). interfere: how deep an overlap may pass unseen, mm (default: the model's accuracy, as a claim's).",
    ),
    ("out", r#""part.stl""#, "The file written: mesh .stl/.3mf/.glb, slice .dxf (cut at the model's accuracy), render .png."),
    ("max-res", "1", "mesh: caps the cell size ink picks for a target (not a given --res)."),
    ("t-span", r#""0:2""#, "A:B: check claims that read t at every moment in [A, B] (default: the file's t-span)."),
    ("name", r#""fit""#, "assert: only this claim. check: only this def's value, whole (to see a value, def it)."),
    ("own", "true", "Only this file's own claims, not its imports'."),
    ("rays", "[[[0, 0, 50], [0, 0, -1]]]", "Rays [[[ox,oy,oz],[dx,dy,dz]],...], or @file.json, or - (stdin)."),
    ("from", "[0, 0, 50]", "One ray's origin, with dir."),
    ("dir", "[0, 0, -1]", "One ray's direction, with from."),
    ("a", r#""root/arm""#, "The first target."),
    ("b", r#""root/base""#, "The second target."),
    ("view", r#""back""#, "One view: iso top front right back left bottom, or [x,y,z] to look from."),
    ("cut", r#""z=25""#, "AXIS[=AT]: cut away the half toward the camera."),
    ("size", "1024", "The image's side, px, 64 to 4096."),
    ("contents", "true", "The headings (default)."),
    ("section", r#""Motion""#, "A section by title or path (spec Tour: a short tour), a form by name (polygon), else those containing these words."),
    ("all", "true", "All of it."),
];

/// Examples of an option for one command, where OPTION_HELP's is not one
/// it takes.
const EXAMPLES: &[(&str, &str, &str)] = &[
    ("interfere", "target", r#""root""#),
    ("check", "name", r#""reach""#),
    ("assert", "keys", r#""status,accuracy""#),
    ("ray", "keys", r#""nearest.point""#),
    ("mass", "keys", r#""mass,inertia""#),
    ("slice", "keys", r#""kind,points""#),
    ("slice", "out", r#""part.dxf""#),
    ("render", "out", r#""part.png""#),
    ("std", "section", r#""belts and chains""#),
];

/// An option's example for a command, and its help (OPTION_HELP).
pub(crate) fn option(cmd: &str, k: &str) -> (&'static str, &'static str) {
    let (ex, help) = OPTION_HELP.iter().find(|o| o.0 == k).map_or(("", ""), |o| (o.1, o.2));
    (EXAMPLES.iter().find(|e| (e.0, e.1) == (cmd, k)).map_or(ex, |e| e.2), help)
}

/// `ink help COMMAND`: what it does, and what each option it reads does
/// (those every command on a file reads too).
fn usage_of(cmd: &str) -> Result<String, InkError> {
    if let Some(c) = registry().commands().find(|c| c.name == cmd) {
        return Ok(format!("ink {cmd}: {}\n", c.help));
    }
    let Some((_, help, own)) = COMMANDS.iter().find(|c| c.0 == cmd && !c.1.is_empty()) else {
        let names = COMMANDS.iter().filter(|c| !c.1.is_empty()).map(|c| c.0).chain(registry().commands().map(|c| c.name));
        return Err(InkError { hint: ink::lang::closest(cmd, names), ..InkError::msg(format!("no command '{cmd}'")) });
    };
    let no_file = NO_FILE.contains(&cmd);
    let mut s = format!("usage: ink {cmd}{} [options]\n{help}\n", if no_file { "" } else { " <file.ink>" });
    let shared = MODEL_OPTIONS.iter().filter(|_| !no_file).chain(PROCESS_OPTIONS.iter().filter(|_| !TEXT.contains(&cmd)));
    let opts: Vec<&str> = own.iter().map(|o| o.name).chain(shared.copied()).collect();
    if !opts.is_empty() {
        s += "\noptions:\n";
    }
    for k in opts {
        s += &format!("  --{k:<10} {}\n", option(cmd, k).1);
    }
    Ok(s)
}

/// `ink help`: the commands (the core's, then the packages'), each with its
/// options, then what each option does.
fn usage() -> String {
    let mut s = String::from(
        "ink - an implicit geometry kernel for AI agents\n\nusage: ink <command> <file.ink> [options]; ink help <command>: one command and its options.\nAnswers are JSON; help, spec and std print text.\n\ncommands:\n",
    );
    let line = |name: &str, help: &str, opts: &[Opt]| {
        let opts: String = opts.iter().map(|o| format!(" --{}", o.name)).collect();
        format!("  {name:<11} {help}\n{}", if opts.is_empty() { String::new() } else { format!("              {}\n", opts.trim_start()) })
    };
    for (name, help, opts) in COMMANDS.iter().filter(|c| !c.1.is_empty()) {
        s += &line(name, help, opts);
    }
    s += "\npackage commands:\n";
    for c in registry().commands() {
        s += &format!("  {:<11} {}\n", c.name, c.help);
    }
    s += "\noptions (every command on a file takes --param and --t, every JSON answer --compact; others only where listed):\n";
    for (k, _, help) in OPTION_HELP {
        s += &format!("  --{k:<10} {help}\n");
    }
    s += "\nUnits in options and answers: mm, g, N, MPa, N mm, degC, s, deg (ink spec --section Units).\nAn option may give its own: --res 0.01in, --param wall=0.25in.\n\nNew to ink: read docs/GUIDE.md.\n";
    s
}

/// The core's commands that read no file, and of them those that print
/// text (mcp serves the others) rather than a JSON answer.
pub const NO_FILE: &[&str] = &["help", "std", "spec", "mcp"];
const TEXT: &[&str] = &["help", "std", "spec", "mcp"];
/// The options every command answering in JSON reads (the answer's layout:
/// `publish`), and those every command on a file reads besides (the
/// model's loading).
const PROCESS_OPTIONS: &[&str] = &["compact"];
pub const MODEL_OPTIONS: &[&str] = &["param", "t", "nominal"];
/// The options read as many times as given (each --param sets one, each
/// --keys adds keys); any other reads one value, so giving it twice is
/// refused rather than answered for the last alone.
const REPEATED: &[&str] = &["param", "keys"];

/// A package's command that reads words, not a file (Run::Words).
pub fn words(reg: &Registry, cmd: &str) -> Option<&'static ink::lang::command::Command> {
    reg.command(cmd).filter(|c| matches!(c.run, Run::Words { .. }))
}

/// Whether `cmd` reads no file: a core command's (NO_FILE) or a package's
/// that reads words.
pub fn no_file(reg: &Registry, cmd: &str) -> bool {
    NO_FILE.contains(&cmd) || words(reg, cmd).is_some()
}

/// The core's commands and those of the registry's packages, each with the
/// options it reads (besides those every command on a file reads).
pub fn commands(reg: &Registry) -> Vec<(&'static str, &'static [Opt])> {
    COMMANDS.iter().map(|c| (c.0, c.2)).chain(reg.commands().map(|c| (c.name, c.options))).collect()
}

/// Err: an option the command does not read (COMMANDS) or given more
/// times than it reads (REPEATED), or an argument after the file.
fn check(reg: &Registry, a: &Args) -> Result<(), InkError> {
    let all = commands(reg);
    // An unknown command, before its file is read (a package's this build
    // lacks).
    let Some((_, own)) = all.iter().find(|(c, _)| *c == a.cmd) else {
        if matches!(a.cmd.as_str(), "--help" | "-h") {
            return Ok(());
        }
        let names: Vec<&str> = all.iter().map(|(c, _)| *c).collect();
        let hint = ink::lang::closest(&a.cmd, names.iter().copied()).unwrap_or_else(|| format!("commands: {}", names.join(" ")));
        return Err(InkError { hint: Some(hint), ..InkError::msg(format!("unknown command '{}'", a.cmd)) });
    };
    let no_file = no_file(reg, &a.cmd);
    // A package's command may read words instead of a file.
    let x = match no_file {
        _ if words(reg, &a.cmd).is_some() => None,
        // ink help COMMAND: that command's help.
        true if a.cmd == "help" => a.extra.first(),
        true => a.file.as_ref(),
        false => a.extra.first(),
    };
    if let Some(x) = x {
        let hint = match a.cmd.as_str() {
            "spec" => "ink spec describes the language and takes no file; use ink spec (its headings) or ink spec --section NAME; use ink check FILE to inspect a model".into(),
            _ if no_file => format!("usage: ink {} [options]; options start with --", a.cmd),
            _ => format!("usage: ink {} <file.ink> [options]; options start with --", a.cmd),
        };
        return Err(InkError {
            hint: Some(hint),
            ..InkError::msg(format!("unexpected argument '{x}': '{}' takes {}", a.cmd, if no_file { "no positional arguments" } else { "one file" }))
        });
    }
    let shared = PROCESS_OPTIONS.iter().filter(|_| !TEXT.contains(&a.cmd.as_str())).chain(MODEL_OPTIONS.iter().filter(|_| !no_file));
    let known: Vec<String> = shared.copied().chain(own.iter().map(|o| o.name)).map(|k| format!("--{k}")).collect();
    let mut keys: Vec<&String> = a.opts.keys().collect();
    keys.sort();
    for k in keys {
        if !known.iter().any(|o| o[2..] == **k) {
            let elsewhere: Vec<&str> = all.iter().filter(|(_, o)| o.iter().any(|o| o.name == k)).map(|(c, _)| *c).collect();
            let hint = match elsewhere.is_empty() {
                false => format!("--{k} is an option of {}, not of '{}'", elsewhere.join(", "), a.cmd),
                true if known.is_empty() => format!("'{}' takes no options", a.cmd),
                true => ink::lang::closest(&format!("--{k}"), known.iter().map(String::as_str)).unwrap_or_else(|| format!("'{}' takes {}", a.cmd, known.join(" "))),
            };
            return Err(InkError { hint: Some(hint), ..InkError::msg(format!("unknown option --{k} for '{}'", a.cmd)) });
        }
        let given = &a.opts[k];
        if given.len() > 1 && !REPEATED.contains(&k.as_str()) {
            let values: Vec<String> = given.iter().map(|v| format!("'{v}'")).collect();
            return Err(InkError {
                hint: Some(format!("give it once: a call per value (only --{} repeat)", REPEATED.join(" and --"))),
                ..InkError::msg(format!("--{k} takes one value, given {}: {}", given.len(), values.join(", ")))
            });
        }
    }
    Ok(())
}

/// The registry the command line runs with and its language reference
/// (`ink spec`): the kernel's (`ink::registry()`, docs/IR.md), or a
/// distribution's, given to [`main`].
static DISTRIBUTION: std::sync::OnceLock<(&'static Registry, &'static str)> = std::sync::OnceLock::new();

pub fn registry() -> &'static Registry {
    DISTRIBUTION.get().map_or_else(ink::registry, |d| d.0)
}

pub fn reference() -> &'static str {
    DISTRIBUTION.get().map_or(ink::SPEC, |d| d.1)
}

/// The `ink` command line, with the packages `reg` holds and `spec` its
/// language reference: the kernel's binary passes `ink::registry()` and
/// `ink::SPEC`, a distribution its own.
pub fn main(reg: &'static Registry, spec: &'static str) {
    DISTRIBUTION.set((reg, spec)).ok();
    let args: Vec<String> = std::env::args().skip(1).collect();
    // `ink mcp` runs on this thread: its calls compute in processes of their
    // own, but a package's command it answers itself (Run::Words), on a
    // thread of its own (mcp::call), and it
    // blocks signals before any other thread starts (mcp::serve).
    if args.first().is_some_and(|c| c == "mcp") {
        std::process::exit(run(args));
    }
    let worker = std::thread::Builder::new().stack_size(256 << 20).spawn(move || run(args)).expect("spawn");
    std::process::exit(worker.join().unwrap_or(2));
}

pub fn run(raw: Vec<String>) -> i32 {
    let t0 = Instant::now();
    let a = Args::parse(&raw);
    if let Err(e) = check(registry(), &a) {
        return fail(e, &a);
    }
    match a.cmd.as_str() {
        "help" | "--help" | "-h" => {
            return match a.file.as_deref().map(usage_of) {
                None => {
                    print!("{}", usage());
                    0
                }
                Some(Ok(text)) => {
                    print!("{text}");
                    0
                }
                Some(Err(e)) => fail(e, &a),
            };
        }
        "std" => {
            return match std_lib(&a) {
                Ok(text) => {
                    print!("{text}");
                    0
                }
                Err(e) => fail(e, &a),
            };
        }
        "spec" => {
            return match spec(&a) {
                Ok(text) => {
                    print!("{text}");
                    0
                }
                Err(e) => fail(e, &a),
            };
        }
        "mcp" => return mcp::serve(),
        _ => {}
    }
    // The call's time limit (INK_DEADLINE, Unix seconds: `ink mcp` sets
    // it): searches stop in time to answer whole, with what they decided
    // (ink::stop), as a package's long solve may. Read here, in the call's own process,
    // not in `answer`: ink mcp answers some commands in its own (`printed`).
    if let Some(at) = ink::lang::command::limit() {
        ink::stop::at(at - ANSWER_S);
    }
    let (text, code) = answer(&a, t0, &|| false);
    print!("{text}");
    code
}

/// The exit status of an answer that succeeded: [`NOT_ALL_HOLD`] where
/// `ink assert` answered with a claim that does not hold (fails or
/// unknown), else 0. An error exits 1.
fn status(a: &Args, v: &J) -> i32 {
    if a.cmd == "assert" && v["all_hold"] == json!(false) { NOT_ALL_HOLD } else { 0 }
}

/// `ink assert`'s exit status when its answer has a claim that does not
/// hold: not 1 (an error, no answer), so a script tells them apart.
pub const NOT_ALL_HOLD: i32 = 3;

/// What a command answering in JSON prints (`run`), and its exit status
/// (`status`; 1 an error); a command reading words stopped once `stop`
/// says so.
fn answer(a: &Args, t0: Instant, stop: &(dyn Fn() -> bool + Sync)) -> (String, i32) {
    let answer = match (a.cmd.as_str(), words(registry(), &a.cmd)) {
        (_, Some(Command { run: Run::Words { run, .. }, .. })) => run(registry(), a, stop),
        _ => dispatch(registry(), a),
    };
    match answer {
        Ok(v) => {
            let code = status(a, &v);
            (format!("{}\n", publish(v, a, t0.elapsed().as_secs_f64() * 1e3)), code)
        }
        Err(e) => (failed(e, a), 1),
    }
}

/// Whether ink mcp answers `ink ARGV` in its own process (`printed`): a
/// package's command reading words, where it says so for these arguments
/// (Run::Words `here`).
pub fn answered_here(argv: &[String]) -> bool {
    let a = Args::parse(argv);
    matches!(words(registry(), &a.cmd), Some(Command { run: Run::Words { here, .. }, .. }) if here(&a))
}

/// What `ink ARGV` (`answered_here`) prints, and whether it succeeded,
/// answered in this process; it stops once `stop` says so.
pub fn printed(argv: &[String], stop: &(dyn Fn() -> bool + Sync)) -> (String, bool) {
    let t0 = Instant::now();
    let a = Args::parse(argv);
    match check(registry(), &a) {
        Ok(()) => {
            let (text, code) = answer(&a, t0, stop);
            (text, code != 1)
        }
        Err(e) => (failed(e, &a), false),
    }
}

/// Seconds before the call's time limit at which searches stop: what
/// follows (the cells being searched finishing, the claims after it not
/// begun, the answer written) took 0.7 to 0.95 s on a wrist's 4-param
/// clearance (6 threads, Ryzen AI MAX+ 395).
const ANSWER_S: f64 = 3.0;

/// An answer as printed: `ok` first, then its fields, marked `stopped`
/// where the call's time limit stopped its searches, with its time; compact
/// or not.
fn publish(mut v: J, a: &Args, ms: f64) -> String {
    if let J::Object(m) = &mut v {
        let mut out = serde_json::Map::new();
        out.insert("ok".into(), json!(true));
        out.append(m);
        // Stopped by the call's time limit (ink::stop): what it had then, as
        // ink mcp marks a call it stops.
        if ink::stop::reached() {
            out.insert("stopped".into(), json!(true));
        }
        out.insert("ms".into(), n(ms));
        v = J::Object(out);
    }
    match a.get("compact") {
        Some(_) => to_compact(&v),
        None => to_pretty(&v),
    }
}

fn fail(e: InkError, a: &Args) -> i32 {
    print!("{}", failed(e, a));
    1
}

/// An error as printed.
fn failed(e: InkError, a: &Args) -> String {
    let mut v = json!({"ok": false, "error": e.to_json()});
    // Stopped by the call's time limit (a load it cut short), as `publish`
    // marks an answer.
    if ink::stop::reached() {
        v["stopped"] = json!(true);
    }
    match a.get("compact") {
        Some(_) => format!("{}\n", to_compact(&v)),
        None => format!("{}\n", to_pretty(&v)),
    }
}

/// A command on a file, its model loaded with the packages of `reg`.
fn dispatch(reg: &'static Registry, a: &Args) -> Result<J, InkError> {
    let mut s = session(reg, a)?;
    // A package's command is its own (ink::lang::command).
    if let Some(Command { run: Run::Model(run), .. }) = reg.command(&a.cmd) {
        return run(&mut s, a);
    }
    // A published mass's density is measured only where a command reads
    // it: the attributes it weighs or evaluates (attrs_read).
    // eval measures one only where --keys names it: its default record
    // shows an unmeasured density by name.
    if matches!(a.cmd.as_str(), "eval" | "mass") {
        let r = Records::new(a);
        let named = |attrs: &Attrs| -> Vec<u32> {
            match a.cmd.as_str() {
                "eval" => attrs_read(a, attrs)
                    .into_iter()
                    .filter(|_| r.names("attrs"))
                    .filter(|&f| attrs.iter().any(|(k, v)| matches!(v, Value::Field(g, _) if *g == f) && (r.names(&format!("attrs.{k}")) || r.names_exactly("attrs"))))
                    .collect(),
                _ => attrs_read(a, attrs),
            }
        };
        let mut fields: Vec<u32> = Vec::new();
        for it in s.items(a)? {
            fields.extend(named(&it.attrs).into_iter().map(|f| s.doc.place(f, &it)));
        }
        s.want(&fields)?;
    }
    let mut out = run_command(&mut s, a)?;
    // Queries reading inexact inputs (Doc::inexact: an imported mesh, a
    // baked field, a solved value not certified) are measured: say which.
    let names: Vec<String> = match a.cmd.as_str() {
        // As Session::items: a whole root, else the bodies (one each).
        "eval" | "bounds" | "slice" | "mesh" | "components" | "render" | "ray" | "mass" => {
            vec![a.get("target").unwrap_or("root").to_string()]
        }
        "distance" => ["a", "b"].iter().filter_map(|k| a.get(k)).map(str::to_string).collect(),
        _ => vec![],
    };
    let targets: Vec<Result<Vec<Item>, String>> = match a.cmd.as_str() {
        "distance" => names.iter().map(|t| s.doc.items(t)).collect(),
        _ if names.is_empty() => vec![],
        _ => vec![s.items(a).map_err(|e| e.message)],
    };
    let items: Vec<Item> = targets.into_iter().flatten().flatten().collect();
    // The fields the command reads: the items' where they are placed
    // (Doc::placed_reads), the attributes it evaluates, integrates or
    // weighs (placed: Doc::place). (A published density shown by name,
    // unmeasured, is not read.)
    let mut roots = s.doc.placed_reads(&items);
    for it in &items {
        roots.extend(attrs_read(a, &it.attrs).into_iter().map(|f| s.doc.place(f, it)));
    }
    let mut seen = std::collections::HashSet::new();
    roots.retain(|&f| seen.insert(f) && s.doc.unmeasured(&[f]).is_empty());
    ink::lang::command::inexact_reads(&s, &mut out, &roots);
    Ok(out)
}

/// The fields of an item's attributes the command reads: eval's (when
/// its records hold them), mass's density (unless a number is given).
fn attrs_read(a: &Args, attrs: &Attrs) -> Vec<u32> {
    let names: Vec<&str> = match a.cmd.as_str() {
        "eval" if Records::new(a).wants("attrs") => attrs.iter().map(|(k, _)| k.as_str()).collect(),
        "mass" => match a.get("density") {
            Some(v) if ink::lang::materials::density_option(v).is_some() => vec![],
            v => vec![v.unwrap_or("density")],
        },
        _ => vec![],
    };
    fn fields(v: &Value, out: &mut Vec<u32>) {
        match v {
            Value::Field(f, _) => out.push(*f),
            Value::Vec(items) => items.iter().for_each(|x| fields(x, out)),
            _ => {}
        }
    }
    let mut out = Vec::new();
    attrs.iter().filter(|(k, _)| names.contains(&k.as_str())).for_each(|(_, v)| fields(v, &mut out));
    out
}

fn run_command(s: &mut Session, a: &Args) -> Result<J, InkError> {
    match a.cmd.as_str() {
        "check" => cmd_check(s, a),
        "eval" => cmd_eval(s, a),
        "bounds" => cmd_bounds(s, a),
        "slice" => cmd_slice(s, a),
        "mesh" => cmd_mesh(s, a),
        "components" => cmd_components(s, a),
        "interfere" => cmd_interfere(s, a),
        "assert" => cmd_assert(s, a),
        "render" => cmd_render(s, a),
        "ray" => cmd_ray(s, a),
        "mass" => cmd_mass(s, a),
        "distance" => cmd_distance(s, a),
        other => Err(InkError {
            hint: Some(format!("commands: {}", commands(s.doc.registry).iter().map(|(c, _)| *c).collect::<Vec<_>>().join(" "))),
            ..InkError::msg(format!("unknown command '{other}'"))
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ink::lang::{Form, Package, Reader};
    use ink::sexp::Sexp;

    /// A third party's package: one command, `echo`, and no forms.
    struct Echo;

    fn echo(s: &mut ModelView, a: &Args) -> Result<J, InkError> {
        Ok(json!({"word": a.get("word"), "has_b": s.doc.body("b").is_some(), "t": s.ctx.t}))
    }

    static ECHO: [Command; 1] = [Command { name: "echo", options: &[Opt::reads("word")], help: "say --word", tool: r#"says the word. {"word": "hi"}"#, run: Run::Model(echo) }];

    impl Package for Echo {
        fn name(&self) -> &'static str {
            "echo"
        }
        fn forms(&self) -> &'static [&'static str] {
            &[]
        }
        fn read(&self, _: &mut Reader, _: &Sexp, head: &str, _: &[Sexp]) -> Result<Arc<dyn Form>, InkError> {
            Err(InkError::msg(format!("echo reads no form '{head}'")))
        }
        fn commands(&self) -> &'static [Command] {
            &ECHO
        }
        fn check(&self, _: &mut ModelView) -> serde_json::Map<String, J> {
            serde_json::Map::from_iter([("echo_checked".to_string(), json!(true))])
        }
    }

    static REGISTRY: Registry = Registry::new(&[&Echo]);

    /// A package's command reaches the command line through the registry
    /// its model is loaded with, as the standard packages' do: its options
    /// are checked as its own, it runs on the loaded model at the moment
    /// asked, and a registry without the package has no such command.
    #[test]
    fn a_package_command_is_dispatched_through_the_registry() {
        let dir = std::env::temp_dir().join(format!("ink-echo-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("m.ink");
        std::fs::write(&file, "(import \"std\")\n(body b (shape (sphere 1)))\n").unwrap();
        let argv = |more: &[&str]| Args::parse(&["echo", &file.to_string_lossy()].iter().chain(more).map(|s| s.to_string()).collect::<Vec<_>>());
        let a = argv(&["--word", "hi", "--t", "0.5"]);
        check(&REGISTRY, &a).unwrap();
        assert_eq!(dispatch(&REGISTRY, &a).unwrap(), json!({"word": "hi", "has_b": true, "t": 0.5}));
        assert!(commands(&REGISTRY).contains(&("echo", &[Opt::reads("word")][..])));
        let typo = check(&REGISTRY, &argv(&["--wrod", "hi"])).unwrap_err();
        assert_eq!(typo.message, "unknown option --wrod for 'echo'");
        // What the package adds to `ink check`'s answer (Package::check).
        let checked = dispatch(&REGISTRY, &Args::parse(&["check".to_string(), file.to_string_lossy().to_string()])).unwrap();
        assert_eq!(checked["echo_checked"], json!(true), "{checked}");
        let absent = dispatch(registry(), &a).unwrap_err();
        assert_eq!(absent.message, "unknown command 'echo'");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
