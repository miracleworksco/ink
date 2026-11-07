//! Refusals of bad arguments, and the plain queries: eval, bounds, mass's volume,
//! slice (its DXF too) and components.

use super::*;
use serde_json::json;
use std::collections::HashMap;

#[test]
fn an_error_in_a_model_says_its_line() {
    let s = Scratch::new("arity");
    let bad = s.model("bad.ink", "(import \"std\")\n(body b (shape (sphere 1 2)))\n");
    let d = ink_err(["check", &bad]);
    assert!(d["error"]["at"]["line"] == 2, "{}", pretty(&d));
}

#[test]
fn a_cell_size_or_moment_that_is_not_a_finite_positive_number_is_refused() {
    let sphere = sphere();
    for v in ["0", "-1", "inf", "nan"] {
        let e = &ink_err(["assert", &sphere, "--res", v])["error"];
        assert!(e["message"].as_str().unwrap().contains(&format!("--res expects a cell size above 0, got '{v}'")), "--res {v}: {}", pretty(e));
    }
    // A moment is a finite number of seconds: --t nan is refused, not a
    // solve on t giving its load-time values, as if solved.
    for v in ["nan", "inf"] {
        let e = &ink_err(["check", &sphere, "--t", v])["error"];
        assert!(e["message"].as_str().unwrap().contains(&format!("--t expects a moment, a finite number of seconds, got '{v}'")), "--t {v}: {}", pretty(e));
    }
}

/// `ink ARGS` is refused, its message containing `says` and its hint `hint`.
fn refused(args: &[&str], says: &str, hint: Option<&str>) {
    let e = &ink_err(args)["error"];
    assert!(
        e["message"].as_str().unwrap().contains(says) && hint.is_none_or(|h| e["hint"].as_str().unwrap_or("").contains(h)),
        "ink {}: wanted {says:?} (hint {hint:?}):\n{}",
        args.join(" "),
        pretty(e)
    );
}

#[test]
fn an_option_a_command_does_not_read_or_a_stray_argument_is_refused_with_a_hint() {
    // A mistyped --t-span would leave a claim proven at t = 0 only.
    let s = Scratch::new("options");
    let moving = s.model("moving-value.ink", "(assert late (value t < 0.5))\n");
    let sphere = sphere();
    let d = ink(["assert", &moving, "--t-span", "0:1"]);
    assert_eq!(verdict(&d["asserts"][0]), "fails", "--t-span is read: {}", pretty(&d));
    let jpg = s.path("x.jpg");
    let cases: &[(&[&str], &str, Option<&str>)] = &[
        (&["assert", &moving, "--tspan", "0:1"], "unknown option --tspan for 'assert'", Some("did you mean '--t-span'?")),
        (&["assert", &moving, "--params", "r=2"], "unknown option --params for 'assert'", Some("did you mean '--param'?")),
        (&["std", "extra"], "unexpected argument 'extra'", None),
        (&["spec", "model.ink"], "unexpected argument 'model.ink'", Some("ink spec --section NAME; use ink check FILE")),
        (&["help", "--typo"], "unknown option --typo for 'help'", None),
        (&["mesh", &sphere, "--format", "stl"], "unknown option --format for 'mesh'", Some("'mesh' takes --compact --param --t --nominal --target --out --res")),
        (&["spec", "--section", "nosuch"], "no section of the reference is called 'nosuch'", Some("ink spec --contents")),
        (&["check", &moving, "--res", "1"], "unknown option --res for 'check'", Some("--res is an option of")),
        (&["assert", &moving, &sphere], &format!("unexpected argument '{sphere}'"), Some("usage: ink assert")),
        (&["render", &sphere, "--out", &jpg], "render writes a PNG", None),
        (&["mesh", &sphere, "--accuracy", "nan"], "--accuracy must be 0 (uniform) or positive", None),
    ];
    for (args, says, hint) in cases {
        refused(args, says, *hint);
    }
}

#[test]
fn eval_reads_points_from_a_file_beyond_what_argv_could_hold() {
    let s = Scratch::new("eval-file");
    let pts: Vec<Value> = (0..50000).map(|i| json!([i as f64 * 0.01, 0, 0])).collect();
    let path = s.model("points.json", &Value::Array(pts).to_string());
    let d = ink(["eval", &sphere(), "--points", &format!("@{path}")]);
    let n = d["targets"][0]["results"]["rows"].as_array().map_or(0, Vec::len);
    assert_eq!(n, 50000, "results for 50000 points");
}

#[test]
fn eval_keys_pick_nested_fields_and_compact_makes_a_table() {
    // Batch output for agents: a table, columns once and a row per point;
    // --keys picks keys (a dot steps into an object); --compact, a row a
    // line. The values are the sphere's field and its derivative in r.
    let sphere = sphere();
    let base = ["eval", &sphere, "--points", "[[3,4,0],[0,0,0]]", "--param", "r=10", "--keys"];
    let want = json!({"columns": ["value", "d_dparams.r"], "rows": [[-5.0, -1.0], [-10.0, -1.0]]});
    let out = Ink::new([&base[..], &["value,d_dparams", "--compact"][..]].concat()).text();
    let t: Value = serde_json::from_str(&out).unwrap();
    let t = &t["targets"][0]["results"];
    assert!(*t == want && out.matches('\n').count() == 4, "--keys picks keys, nested; --compact: a row per line: {out}");
    let d = ink([&base[..], &["value,d_dparams.r"][..]].concat());
    assert!(d["targets"][0]["results"] == want, "the same table without --compact: {}", pretty(&d["targets"][0]["results"]));
    let e = &ink_err(["eval", &sphere, "--points", "[[0,0,0]]", "--keys", "value,d_dparams.q"])["error"];
    assert!(e["message"].as_str().unwrap().contains("d_dparams.q") && e["hint"].as_str().unwrap_or("").contains("attrs"), "an unknown key lists the keys: {}", pretty(e));
}

#[test]
fn bounds_of_a_target_enclose_its_shape() {
    let d = ink(["bounds", &shelf(), "--target", "bracket"]);
    let b = &d["targets"][0];
    assert!(num(&b["min"][0]) <= -17.0 && num(&b["max"][2]) >= 45.0 && num(&b["max"][2]) < 46.0, "{}", pretty(b));
}

#[test]
fn a_curved_bodys_volume_estimate_is_close() {
    // `ink mass`'s volume: an estimate.
    let d = ink(["mass", &sphere(), "--res", "0.2"]);
    let v = &d["targets"][0]["volume"];
    let exact = 4.0 / 3.0 * std::f64::consts::PI * 20f64.powi(3);
    assert!((num(&v["estimate"]) - exact).abs() / exact < 1e-3, "{}", pretty(v));
}

#[test]
fn an_empty_instance_measured_once_gives_no_spread() {
    // A spread is shown only by a halving: `--res` measures once.
    let dir = Scratch::new("empty-instance");
    let m = dir.model("empty.ink", "(import \"std\")\n(body gone (shape (max (sphere 1) (move (sphere 1) 5 0 0))))\n(assembly root (instance g gone))\n");
    let d = ink(["mass", &m, "--res", ".1"]);
    assert!(d["targets"].as_array().unwrap().iter().all(|t| t.get("spread").is_none()), "{}", pretty(&d));
}

#[test]
fn a_mass_key_no_target_has_is_refused() {
    // Every record has its target, so a typo kept only that, silently.
    let sphere = sphere();
    refused(&["mass", &sphere, "--res", "2", "--keys", "mass,absent"], "--keys: no record has 'absent'", Some("center_of_mass"));
    let d = ink(["mass", &sphere, "--res", "2", "--keys", "volume"]);
    let t = &d["targets"][0];
    assert!(t["target"].is_string() && t["volume"].is_object() && t.get("mass").is_none(), "{}", pretty(&d));
}

const ARM: [&str; 8] = ["slice", "", "--target", "bracket", "--axis", "z", "--at", "42"];

fn arm(extra: &[&str]) -> Value {
    let shelf = shelf();
    let mut args = ARM.to_vec();
    args[1] = &shelf;
    args.extend_from_slice(extra);
    ink(args)
}

#[test]
fn a_slice_gives_its_loops_areas_and_bounds_and_points_only_when_asked() {
    let d = arm(&[]);
    let want = 34.0 * 40.0 - 2.0 * std::f64::consts::PI * 2.75f64.powi(2);
    assert!((num(&d["area"]) - want).abs() < 2.0, "arm slice area {} vs {want}", d["area"]);
    let loops = records(&d["loops"]);
    let kinds: Vec<&str> = loops.iter().map(|l| l["kind"].as_str().unwrap()).collect();
    assert!(kinds.iter().filter(|k| **k == "outer").count() == 1 && kinds.iter().filter(|k| **k == "hole").count() == 2 && kinds.len() == 3, "outer and 2 holes: {kinds:?}");
    // A loop's points only when --keys names them (70 loops of a keyboard
    // plate are 22 KB).
    assert!(loops.iter().all(|l| l.get("points").is_none() && l["min"].as_array().map(Vec::len) == Some(2)), "{}", pretty(&d["loops"]));
    let t = arm(&["--keys", "kind,points", "--compact"]);
    let t = &t["loops"];
    assert!(
        t["columns"] == json!(["kind", "points"])
            && t["rows"].as_array().unwrap().iter().all(|r| r[1].as_array().is_some_and(|p| p.len() >= 4 && p[0].as_array().map(Vec::len) == Some(2))),
        "{}",
        pretty(t)
    );
}

#[test]
fn a_dxf_slice_is_one_closed_polyline_per_loop_at_the_models_accuracy() {
    // --out FILE.dxf: the outline to cut, at the model's accuracy.
    let s = Scratch::new("dxf");
    let dxf = s.path("bracket.dxf");
    let d = arm(&["--out", &dxf]);
    assert!(d["res"] == 0.05 && d["accuracy"] == 0.025, "{}", pretty(&d));
    let text = std::fs::read_to_string(&dxf).unwrap();
    let lines: Vec<&str> = text.split('\n').collect();
    let pairs: Vec<(&str, &str)> = lines.as_chunks::<2>().0.iter().map(|&[a, b]| (a, b)).collect();
    let (mut polys, mut verts) = (vec![], vec![]);
    for (k, &p) in pairs.iter().enumerate() {
        if p == ("0", "POLYLINE") {
            polys.push(pairs[k + 1..k + 6].iter().copied().collect::<std::collections::HashMap<_, _>>());
        }
        if p == ("0", "VERTEX") {
            verts.push((pairs[k + 2].1.trim().parse::<f64>().unwrap(), pairs[k + 3].1.trim().parse::<f64>().unwrap()));
        }
    }
    let mm = pairs.windows(2).any(|w| w[0].1 == "$INSUNITS" && w[1].0 == "70");
    assert!(
        polys.len() == 3 && polys.iter().all(|p| p.get("8") == Some(&"bracket") && p.get("70") == Some(&"1")) && mm && pairs.last() == Some(&("0", "EOF")),
        "polylines {polys:?}, $INSUNITS {mm}, last {:?}",
        pairs.last()
    );
    let loops = records(&d["loops"]);
    let (b0, b1) = (&loops[0]["min"], &loops[0]["max"]);
    let inside = |x: f64, y: f64| num(&b0[0]) - 1e-9 <= x && x <= num(&b1[0]) + 1e-9 && num(&b0[1]) - 1e-9 <= y && y <= num(&b1[1]) + 1e-9;
    assert!(!verts.is_empty() && verts.iter().all(|&(x, y)| inside(x, y)), "vertices outside the outline's bounds {b0} {b1}");
}

/// A DXF's entities (table entries included), in order: each its type
/// and its group codes, the first of each.
fn dxf_entities(path: &str) -> Vec<(String, HashMap<String, String>)> {
    let text = std::fs::read_to_string(path).unwrap();
    let lines: Vec<&str> = text.split('\n').collect();
    let mut ents: Vec<(String, HashMap<String, String>)> = vec![];
    for &[c, v] in lines.as_chunks::<2>().0 {
        if c == "0" {
            ents.push((v.to_string(), HashMap::new()));
        } else if let Some(e) = ents.last_mut() {
            e.1.entry(c.to_string()).or_insert_with(|| v.to_string());
        }
    }
    ents
}

/// A model of one 40 x 20 x 2 box (x -20..20, y -10..10) marked by MARK.
fn marked_box(s: &Scratch, name: &str, mark: &str) -> String {
    s.model(name, &format!("(import \"std\")\n(body flat (shape (box 40 20 2)) {mark})\n"))
}

#[test]
fn a_mark_sharing_a_cut_layer_is_refused() {
    // A cut layer and a mark's are told apart by name, as DXF reads them
    // (any case, its forbidden characters as _): sharing one would cut a
    // mark or dash an outline.
    let s = Scratch::new("dxf-clash");
    let out = s.path("clash.dxf");
    for mark in ["(mark flat (vec x -1))", "(mark FLAT (vec x -1))"] {
        let m = marked_box(&s, "clash.ink", mark);
        refused(&["slice", &m, "--target", "flat", "--out", &out], "would share DXF layer FLAT", Some("rename the mark"));
    }
}

#[test]
fn a_mark_ends_at_its_parts_bounds() {
    // Its region is everywhere: the line x = 0 runs across the part, y -10..10.
    let s = Scratch::new("dxf-mark-bounds");
    let m = marked_box(&s, "scribe.ink", "(mark scribe (vec x -1))");
    let d = ink(["slice", &m, "--target", "flat", "--keys", "kind,points"]);
    let marks: Vec<Value> = records(&d["loops"]).into_iter().filter(|l| l["kind"] == "mark").collect();
    let ys: Vec<f64> = marks.iter().flat_map(|l| l["points"].as_array().unwrap().iter().map(|p| num(&p[1])).collect::<Vec<_>>()).collect();
    let (lo, hi) = ys.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &y| (a.min(y), b.max(y)));
    assert!(marks.len() == 1 && (lo + 10.0).abs() < 0.1 && (hi - 10.0).abs() < 0.1, "the mark from y {lo} to {hi}: {}", pretty(&d));
}

#[test]
fn a_curved_marks_dashes_run_on_across_its_vertices() {
    // A polyline's dash pattern restarts at each vertex unless its flags
    // carry 128 (continuous linetype generation): short chords draw solid.
    let s = Scratch::new("dxf-mark-circle");
    let m = s.model("circle.ink", "(import \"std\")\n(body flat (shape (box 60 60 2)) (mark circle (vec (- (len2 x y) 20) -1)))\n");
    let out = s.path("circle.dxf");
    ink(["slice", &m, "--target", "flat", "--out", &out]);
    let ents = dxf_entities(&out);
    let flags: Vec<(&str, &str)> = ents.iter().filter(|(t, _)| t == "POLYLINE").map(|(_, e)| (e["8"].as_str(), e["70"].as_str())).collect();
    assert_eq!(flags, [("flat", "1"), ("circle", "129")], "the outline closed; the mark closed, its linetype continuous");
}

/// x = C lies on a column of the slicing grid's corners for a 40 x 20 x 2
/// box sliced at z = 0, res 0.125: fields 0 there are 0 at the corners.
const GRID_X: &str = "-0.019531249999989342";

#[test]
fn a_field_touching_zero_on_grid_corners_is_reported_and_one_crossing_there_draws_one() {
    // -(x - C)^2 is 0 on the corners at x = C and negative beside them: it
    // touches 0 there, never crossing. Nothing shows its 0 has no width,
    // so the corners read as air: lines drawn there are reported. x - C
    // crosses there: one line through the corners, nothing unresolved.
    let s = Scratch::new("dxf-mark-corner");
    for (g, lines) in [(format!("(neg (square (- x {GRID_X})))"), 0), (format!("(- x {GRID_X})"), 1)] {
        let m = marked_box(&s, "corner.ink", &format!("(mark m (vec {g} -1))"));
        let out = s.path("corner.dxf");
        let d = ink(["slice", &m, "--target", "flat", "--at", "0", "--res", "0.125", "--out", &out]);
        let on: Vec<(String, String)> = dxf_entities(&out).into_iter().filter(|(t, _)| t == "POLYLINE" || t == "LINE").map(|(t, e)| (t, e["8"].clone())).collect();
        let marks = on.iter().filter(|(_, l)| l == "m").count();
        let reported = !d["unresolved"].is_null();
        assert!(on.iter().all(|(t, l)| l == "m" && t == "LINE" || l == "flat" && t == "POLYLINE"), "{g}: {on:?}");
        assert!(if lines == 0 { marks == 0 || reported } else { marks == lines && !reported }, "{g}: {on:?} {}", pretty(&d));
    }
}

#[test]
fn an_outline_touching_zero_on_grid_corners_is_one_loop_or_reported() {
    // The plate's field is 0 on the corners at x = C, negative beside
    // them. Nothing shows that 0 has no width: if the outline is cut there
    // into two halves, the slice says so.
    let s = Scratch::new("slice-corner");
    let m = s.model("slit.ink", &format!("(import \"std\")\n(body flat (shape (max (box 40 20 2) (neg (square (- x {GRID_X}))))))\n"));
    let d = ink(["slice", &m, "--target", "flat", "--at", "0", "--res", "0.125"]);
    let kinds: Vec<Value> = records(&d["loops"]).iter().map(|l| l["kind"].clone()).collect();
    assert!(kinds == [json!("outer")] && (num(&d["area"]) - 800.0).abs() < 1.0 || !d["unresolved"].is_null(), "{}", pretty(&d));
}

#[test]
fn a_mark_in_the_slice_plane_draws_nothing_and_leaves_the_cut() {
    // g = z is 0 over the slice at z = 0: no sign changes there, so no
    // curve, as a field touching 0 makes no material; the outline is cut.
    let s = Scratch::new("dxf-mark-plane");
    let m = marked_box(&s, "plane.ink", "(mark level (vec z -1))");
    let out = s.path("plane.dxf");
    let d = ink(["slice", &m, "--target", "flat", "--at", "0", "--out", &out]);
    let ents = dxf_entities(&out);
    let on: Vec<(&str, &str)> = ents.iter().filter(|(t, _)| t == "POLYLINE" || t == "LINE").map(|(t, e)| (t.as_str(), e["8"].as_str())).collect();
    assert!(on == [("POLYLINE", "flat")] && d["loop_count"] == 1, "{on:?}: {}", pretty(&d));
}

#[test]
fn a_sheets_dxf_has_its_bend_lines_dashed_on_a_layer_of_their_own() {
    // A flat 60 wide, its top edge from (0, 30) up to (60, 40), bent along
    // x = 30 and x = 45. Each bend's line is marked from y = 0 to the
    // slanted edge, at y = 30 + x / 6: 35 and 37.5.
    let s = Scratch::new("dxf-bends");
    let m = s.model(
        "u.ink",
        "(import \"std\")\n\
         (def S (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 30))) (vec (move (circle-2d 2) 10 15 0)) 2 0.4\n\
           (vec (sheet-bend (vec 30 50) (vec 30 0) 90 2) (sheet-bend (vec 45 50) (vec 45 0) 90 2))))\n\
         (body flat (shape (sheet-flat S)) (mark bend (sheet-bend-lines S)))\n",
    );
    let dxf = s.path("u.dxf");
    let d = ink(["slice", &m, "--target", "flat", "--out", &dxf]);
    let ents = dxf_entities(&dxf);
    let f = |e: &HashMap<String, String>, c: &str| e[c].trim().parse::<f64>().unwrap();
    let layer = |n: &str| ents.iter().find(|(t, e)| t == "LAYER" && e.get("2").map(String::as_str) == Some(n)).map(|(_, e)| e["6"].as_str());
    let polys: Vec<&str> = ents.iter().filter(|(t, _)| t == "POLYLINE").map(|(_, e)| e["8"].as_str()).collect();
    let mut bends: Vec<[f64; 4]> = ents.iter().filter(|(t, e)| t == "LINE" && e["8"] == "bend").map(|(_, e)| [f(e, "10"), f(e, "20"), f(e, "11"), f(e, "21")]).collect();
    bends.iter_mut().for_each(|b| {
        if b[1] > b[3] {
            *b = [b[2], b[3], b[0], b[1]]
        }
    });
    bends.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let near = |a: f64, b: f64| (a - b).abs() < 0.05;
    let at = |b: &[f64; 4], x, top| near(b[0], x) && near(b[2], x) && near(b[1], 0.0) && near(b[3], top);
    assert!(
        polys == ["flat", "flat"]
            && layer("flat") == Some("CONTINUOUS")
            && layer("bend") == Some("DASHED")
            && bends.len() == 2
            && at(&bends[0], 30.0, 35.0)
            && at(&bends[1], 45.0, 37.5),
        "outline and hole on flat, a dashed line per bend: polylines {polys:?}, bends {bends:?}\n{}",
        pretty(&d)
    );
    let marks: Vec<_> = records(&d["loops"]).into_iter().filter(|l| l["kind"] == "mark").collect();
    assert!(marks.len() == 2 && marks.iter().all(|l| l["mark"] == "bend" && l["closed"] == false), "{}", pretty(&d));
}

#[test]
fn components_count_solids_even_of_a_sheet_thinner_than_its_sizes_res() {
    let d = ink(["components", &shelf(), "--target", "bracket", "--res", "0.3"]);
    assert!(d["targets"][0]["solids"] == 1, "the bracket is one solid: {}", pretty(&d));
    // A sheet thinner than its size's res (1.48): meshed again finer, not 0 solids.
    let s = Scratch::new("components-sheet");
    let sheet = s.model("sheet.ink", "(import \"std\")\n(body sheet (shape (box 235 235 0.5)))\n");
    let d = ink(["components", &sheet]);
    let t = &d["targets"][0];
    assert!(t["solids"] == 1 && num(&t["res"]) < 1.0, "{}", pretty(t));
}

#[test]
fn faces_meeting_flush_are_no_surface_to_a_void_count() {
    // A shape is where its field is negative; its surface is where the sign
    // changes. A pocket cut ending flush with the top face reads 0 across
    // its mouth, air either side: no surface, so it seals no void (it ran
    // the count to its cell budget, unknown). Pieces meeting on a face (a
    // turned hub under a flange with a toothed ring below it) read 0
    // across it inside the metal: air no wider than 0, which no bound
    // shows airless; counted to the accuracy, it is set aside: no void.
    let s = Scratch::new("flush-voids");
    let m = s.model(
        "flush.ink",
        "(import \"std\")\n\
         (body pocket (shape (subtract (move (box 40 40 10) 20 20 5) (move (box 20 20 5) 20 20 7.5))))\n\
         (body turned (shape (min (revolve-z (polygon (vec (vec 0 0) (vec 4.5 0) (vec 4.5 18) (vec 0 18)))) \
           (revolve-z (polygon (vec (vec 4.5 16.5) (vec 7.5 16.5) (vec 7.5 17.5) (vec 4.5 17.5)))) (extrude-z (circle-2d 6) 9.5 16.5))))\n",
    );
    let d = ink(["components", &m, "--target", "pocket"]);
    assert!(d["targets"][0]["cavities"] == 0, "{}", pretty(&d["targets"][0]["proven"]));
    let d = ink(["components", &m, "--target", "turned"]);
    let c = &d["targets"][0]["proven"]["cavities"];
    assert!(d["targets"][0]["cavities"] == 0 && c["set_aside"]["cells"].as_u64().is_some_and(|n| n > 0), "{}", pretty(c));
}

#[test]
fn a_name_whose_only_zero_is_inside_the_metal_is_not_seen() {
    // A rod filling its sleeve's bore exactly, both ends capped: its field
    // and the sleeve's are 0 together across the bore, a zero the field
    // keeps negative around, so no surface. It read as the rod's face.
    let s = Scratch::new("buried-rod");
    let m = s.model(
        "rod.ink",
        "(import \"std\")\n\
         (body part (shape (min (tag sleeve (subtract (cylinder-z 5 20) (cylinder-z 2 30))) (tag rod (cylinder-z 2 20)) \
           (tag cap (move (cylinder-z 5 2) 0 0 10)) (tag cap (move (cylinder-z 5 2) 0 0 -10)))))\n",
    );
    let d = ink(["check", &m]);
    let w = d["warnings"].to_string();
    assert!(w.contains("'rod' was not seen"), "{}", pretty(&d["warnings"]));
}

#[test]
fn ray_gives_the_fields_its_help_lists() {
    // Help said `face`, which ray answered as `name` since the rename.
    let help = Ink::new(["help"]).text();
    let line = help.lines().find(|l| l.trim_start().starts_with("ray ")).expect("ray in help");
    let fields: Vec<&str> = line.rsplit(':').next().unwrap().trim_end().trim_end_matches('.').split(',').map(str::trim).collect();
    let d = ink(["ray", &sphere(), "--from", "[0,0,50]", "--dir", "[0,0,-1]"]);
    assert!(fields.len() == 3 && fields.iter().all(|f| d["nearest"].get(*f).is_some()), "{line}: {}", pretty(&d["nearest"]));
}

#[test]
fn a_ray_from_inside_says_so_at_its_origin() {
    // It answered its origin as if a surface, with a zero normal.
    let s = Scratch::new("ray-inside");
    let m = s.model("box.ink", "(import \"std\")\n(body b (shape (box 4 6 8)))\n");
    let d = ink(["ray", &m, "--from", "[0,0,0]", "--dir", "[1,0,0]"]);
    let h = &d["nearest"];
    assert!(h["starts_inside"] == true && num(&h["t"]) == 0.0 && h["normal"].is_null(), "{}", pretty(&d));
    // From outside: where it enters, with its normal, not starting inside.
    let d = ink(["ray", &m, "--from", "[-10,0,0]", "--dir", "[1,0,0]"]);
    let h = &d["nearest"];
    assert!((num(&h["t"]) - 8.0).abs() < 1e-6 && (num(&h["normal"][0]) + 1.0).abs() < 1e-6 && h.get("starts_inside").is_none(), "{}", pretty(&d));
}

#[test]
fn an_error_follows_compact_as_an_answer_does() {
    // An error is printed on one line under --compact, so over ink mcp
    // (compact by default) too.
    let s = Scratch::new("compact-error");
    let bad = s.model("bad.ink", "(body b (shape (sphre 2)))\n");
    for args in [vec!["check", bad.as_str(), "--compact"], vec!["check", "no-such-file.ink", "--compact"], vec!["check", bad.as_str(), "--bogus", "1", "--compact"]] {
        let text = Ink::new(&args).run().stdout;
        assert!(text.trim_end().lines().count() == 1 && text.starts_with("{\"ok\":false"), "{args:?}: {text}");
    }
}

#[test]
fn a_mark_s_undecided_samples_reach_the_slice_s_report() {
    // The mark's field touches 0 from below on x = y, a line no face the
    // field writes: samples on it are undecided.
    let s = Scratch::new("mark-unresolved");
    let m = s.model("marks.ink", "(import \"std\")\n(body p (shape (box 12 12 12)) (mark line (vec (neg (abs (- x y))) -1)))\n");
    let d = ink(["slice", &m, "--axis", "z", "--at", "0.3", "--res", "1"]);
    assert!(d["unresolved"]["samples"].as_u64().is_some_and(|n| n > 0), "{}", pretty(&d));
}
