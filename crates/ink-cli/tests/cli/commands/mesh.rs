//! `ink mesh`: STL, 3MF and GLB files, read here from their formats, and
//! what the answer says of each mesh (its res, deviation).

use super::*;
use serde_json::json;
use std::collections::HashMap;
use std::f64::consts::PI;

const BALL: f64 = 4.0 / 3.0 * PI * 8000.0;

/// A float field of a JSON object, by key.
fn f(v: &Value, k: &str) -> f64 {
    num(&v[k])
}

#[test]
fn an_stl_is_closed_says_how_it_was_made_and_matches_the_answer() {
    let s = Scratch::new("stl");
    let stl = s.path("sphere.stl");
    let d = ink(["mesh", &sphere(), "--res", "1", "--out", &stl]);
    let o = &d["objects"][0];
    let (n, bad_edges, vol) = stl_edges(&stl);
    assert_eq!(n as u64, o["triangles"].as_u64().unwrap(), "STL triangles vs the answer's");
    assert_eq!(bad_edges, 0, "a closed 2-manifold (edges not shared by two triangles)");
    assert!((vol - BALL).abs() / BALL < 5e-3, "STL volume {vol} vs {BALL}: its coordinates are the model's");
    let head = String::from_utf8_lossy(&std::fs::read(&stl).unwrap()[..80]).into_owned();
    assert!(head.starts_with("ink ball: res 1 dev ") && head.contains("(sampled)"), "header: {head}");
    // Simplified within res / 30 by default: fewer triangles than at res,
    // its deviation sampled and stated, and what simplification moved.
    let dev = &o["deviation"];
    assert!((f(dev, "simplified_within") - 1.0 / 30.0).abs() < 1e-9 && f(o, "triangles") * 2.0 < f(o, "uniform_triangles"), "{}", pretty(o));
    assert!(f(dev, "simplification") <= f(dev, "simplified_within") && 0.0 < f(dev, "max") && f(dev, "max") < 0.1 && dev["beyond"] == 0, "{}", pretty(dev));
    assert!(o["unresolved"]["cells"] == 0, "no features thinner than res: {}", pretty(&o["unresolved"]));
}

#[test]
fn accuracy_0_meshes_uniformly() {
    // (Its deviation's bound is kernel.rs mesh_deviation_bounds_its_samples'.)
    let s = Scratch::new("stl-uniform");
    let stl = s.path("sphere.stl");
    let d = ink(["mesh", &sphere(), "--res", "1", "--accuracy", "0", "--out", &stl]);
    let o = &d["objects"][0];
    assert!(o["triangles"] == o["uniform_triangles"] && f(&o["deviation"], "simplification") == 0.0, "{}", pretty(o));
}

#[test]
fn a_deviation_past_the_accuracy_asked_is_warned_of_and_the_answer_claims_no_accuracy() {
    // A phone's camera plateau ends 0.44 mm short of its rounded side, under
    // a cell at res 0.6: the mesh at res cuts that step's edge, 0.12 mm off
    // the surface. --accuracy bounds only simplifying, so the answer states
    // no accuracy of its own, and warns that the one asked for is missed.
    let s = Scratch::new("phone");
    let phone = s.model(
        "phone.ink",
        "(import \"std\")\n(body phone (shape (min (extrude-z (rounded-rect-2d 71.45 149.61 12) 0 8.25)\n\
  (move (extrude-z (rounded-rect-2d 44.42 45.74 9) -2.05 0.05) 12.71 51.13 0))))\n",
    );
    let d = ink(["mesh", &phone, "--res", "0.6", "--accuracy", "0.02", "--out", &s.path("phone.stl")]);
    let o = &d["objects"][0];
    let w: Vec<&str> = d["warnings"].as_array().map_or(vec![], |w| w.iter().filter_map(Value::as_str).filter(|x| x.contains("more than --accuracy 0.02")).collect());
    assert!(
        o.get("accuracy").is_none() && f(&o["deviation"], "max") > 0.1 && f(&o["deviation"], "simplified_within") == 0.02 && w.len() == 1 && w[0].contains("bounds only"),
        "{}",
        pretty(&d)
    );
}

#[test]
fn a_smooth_union_is_meshed_whole() {
    // Its blend meets min(a, b) continuously where a and b are k apart:
    // no cell there is left undecided.
    let s = Scratch::new("smooth-union");
    let coin = s.model("coin.ink", "(import \"std\")\n(body coin (shape (extrude-z (smooth-union (circle-2d 20) (move (circle-2d 6) 0 22 0) 3) 0 3)))\n");
    let d = ink(["mesh", &coin, "--out", &s.path("coin.stl")]);
    let quiet = d.get("warnings").is_none_or(|w| w.as_array().is_some_and(Vec::is_empty));
    assert!(d["objects"][0]["unresolved"]["cells"] == 0 && quiet, "{}", pretty(&d));
}

#[test]
fn a_smooth_union_is_its_polynomial_with_its_slopes() {
    // As written with a switch at |a - b| = k (h clamped, mixed): its
    // value, gradient and rates over t and a param, where the blend meets
    // min(a, b), where a = b, and at points about.
    let s = Scratch::new("smooth-union-slopes");
    let m = s.model(
        "both.ink",
        "(import \"std\")\n(param k 1)\n(fn switched (a b k)\n  (let ((h (clamp (+ 0.5 (/ (* 0.5 (- b a)) k)) 0 1)) (m (min a b)))\n\
         (select (- (abs (- a b)) k) (min m (max (- m (* 0.25 k)) (- (mix b a h) (* k h (- 1 h))))) m)))\n\
         (def a (- (len3 x y z) 2))\n(def b (* (+ 1 (* 0.1 t)) (- (len3 (- x 2.5) y z) 1)))\n\
         (body new (shape (min (smooth-union x y k) (move (smooth-union a b (* 0.7 k)) 0 0 10))))\n\
         (body old (shape (min (switched x y k) (move (switched a b (* 0.7 k)) 0 0 10))))\n",
    );
    let mut seed = 7u64;
    let mut r = |lo: f64, hi: f64| {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        lo + (hi - lo) * (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut pts = vec![json!([1, 0, 0]), json!([0, 1, 0]), json!([0, 0, 0]), json!([0.5, 0.5, 0]), json!([-1, 0, 0])];
    pts.extend((0..400).map(|_| json!([r(-3.0, 5.0), r(-3.0, 3.0), r(-2.0, 12.0)])));
    let pts = Value::Array(pts).to_string();
    let rows = |t: &str| ink(["eval", &m, "--target", t, "--points", &pts, "--t", "0.5", "--keys", "value,grad,d_dt,d_dparams"])["targets"][0]["results"]["rows"].clone();
    let (new, old) = (rows("new"), rows("old"));
    let leaves = |v: &Value| v.as_array().unwrap().iter().flat_map(|c| c.as_array().cloned().unwrap_or_else(|| vec![c.clone()])).map(|x| num(&x)).collect::<Vec<f64>>();
    for (i, (a, b)) in new.as_array().unwrap().iter().zip(old.as_array().unwrap()).enumerate() {
        let (a, b) = (leaves(a), leaves(b));
        assert!(a.len() == b.len() && a.iter().zip(&b).all(|(x, y)| (x - y).abs() < 1e-12), "point {i}: {a:?} vs {b:?}");
        // No steeper than its steepest input (b's 1.05 at t = 0.5): its
        // gradient is a mix of a's and b's.
        let g = (a[1] * a[1] + a[2] * a[2] + a[3] * a[3]).sqrt();
        assert!(g <= 1.05 + 1e-12, "point {i}: |grad| {g}");
    }
}

#[test]
fn mesh_warns_of_an_empty_target_it_leaves_out() {
    let s = Scratch::new("mixed-empty");
    let mixed =
        s.model("mixed-empty.ink", "(import \"std\") (body ball (shape (sphere 1))) (body absent (shape 1))\n(assembly root (instance ball ball) (instance absent absent))\n");
    let d = ink(["mesh", &mixed, "--res", "0.5", "--out", &s.path("mixed.stl")]);
    let warned = d["warnings"].as_array().is_some_and(|w| w.iter().any(|w| w.as_str().unwrap().contains("root/absent: empty target")));
    assert!(d["objects"].as_array().unwrap().len() == 1 && warned, "{}", pretty(&d));
}

/// Each `<NAME ...>` start tag in `xml`: its attributes.
fn tags<'a>(xml: &'a str, name: &str) -> Vec<HashMap<&'a str, &'a str>> {
    let open = format!("<{name}");
    let mut out = vec![];
    let mut rest = xml;
    while let Some(i) = rest.find(&open) {
        rest = &rest[i + open.len()..];
        if !rest.starts_with([' ', '>', '/']) {
            continue;
        }
        let end = rest.find('>').unwrap();
        let mut attrs = HashMap::new();
        let mut a = &rest[..end];
        while let Some(eq) = a.find("=\"") {
            let key = a[..eq].trim();
            let close = a[eq + 2..].find('"').unwrap();
            attrs.insert(key, &a[eq + 2..eq + 2 + close]);
            a = &a[eq + 3 + close..];
        }
        out.push(attrs);
        rest = &rest[end..];
    }
    out
}

/// The text between `<NAME ...>` and `</NAME>` for each such element.
fn bodies<'a>(xml: &'a str, name: &str) -> Vec<(&'a str, &'a str)> {
    let (open, close) = (format!("<{name} "), format!("</{name}>"));
    let mut out = vec![];
    let mut rest = xml;
    while let Some(i) = rest.find(&open) {
        let tag_end = i + rest[i..].find('>').unwrap();
        let end = tag_end + rest[tag_end..].find(&close).unwrap();
        out.push((&rest[i..=tag_end], &rest[tag_end + 1..end]));
        rest = &rest[end..];
    }
    out
}

#[test]
fn a_3mf_has_a_mesh_per_distinct_body_and_an_item_per_instance() {
    // The shelf has four distinct meshes and eleven instances.
    let s = Scratch::new("3mf");
    let tmf = s.path("shelf.3mf");
    let d = ink(["mesh", &shelf(), "--res", "1", "--accuracy", "0", "--out", &tmf]);
    let files = super::zip::entries(&std::fs::read(&tmf).unwrap());
    let mut names: Vec<&str> = files.iter().map(|(n, _)| n.as_str()).collect();
    names.sort();
    assert_eq!(names, ["3D/3dmodel.model", "[Content_Types].xml", "_rels/.rels"], "3MF parts");
    let model = std::str::from_utf8(&files.iter().find(|(n, _)| n == "3D/3dmodel.model").unwrap().1).unwrap();
    let objects = bodies(model, "object");
    let items = tags(model, "item");
    assert!(
        objects.len() as u64 == d["distinct_meshes"].as_u64().unwrap() && items.len() as u64 == d["instances"].as_u64().unwrap() && items.len() == 11,
        "objects {}, items {}: {}",
        objects.len(),
        items.len(),
        pretty(&d)
    );
    // Each object's indices are in range and it has a colour per vertex.
    let groups: HashMap<&str, usize> = bodies(model, "m:colorgroup").iter().map(|(tag, body)| (tags(tag, "m:colorgroup")[0]["id"], tags(body, "m:color").len())).collect();
    for (tag, body) in &objects {
        let pid = tags(tag, "object")[0]["pid"];
        let nv = tags(body, "vertex").len();
        assert_eq!(groups.get(pid), Some(&nv), "object {tag}: a colour per vertex");
        for t in tags(body, "triangle") {
            for k in ["v1", "v2", "v3", "p1", "p2", "p3"] {
                let i: usize = t[k].parse().unwrap();
                assert!(i < nv, "{tag}: {k} = {i} of {nv} vertices");
            }
        }
    }
    let moved = items.iter().filter(|i| i.contains_key("transform")).count();
    assert!(moved >= 8, "instances carry transforms: {moved}");
    let desc = model.split("<metadata name=\"Description\">").nth(1).map(|r| &r[..r.find('<').unwrap()]);
    assert!(desc.is_some_and(|t| t.contains(" res ") && t.contains(" dev ")), "Description metadata: {desc:?}");
    assert!(d["objects"].as_array().unwrap().iter().all(|o| o["watertight"] == true), "{}", pretty(&d["objects"]));
}

/// A GLB's animation channels: (node, path) -> sampler.
fn tracks(g: &Glb) -> HashMap<(u64, String), Value> {
    let a = &g.json["animations"][0];
    a["channels"]
        .as_array()
        .unwrap()
        .iter()
        .map(|ch| ((ch["target"]["node"].as_u64().unwrap(), ch["target"]["path"].as_str().unwrap().to_string()), a["samplers"][ch["sampler"].as_u64().unwrap() as usize].clone()))
        .collect()
}

fn acc(v: &Value) -> usize {
    v.as_u64().unwrap() as usize
}

#[test]
fn a_glb_moves_each_instance_over_the_span_and_is_made_again_from_the_cache() {
    let s = Scratch::new("glb-motion");
    let glb = s.path("slider-crank.glb");
    let crank = root("benchmarks/10-slider-crank.ink");
    let d = ink(["mesh", &crank, "--res", "1", "--out", &glb]);
    let data = std::fs::read(&glb).unwrap();
    let u = |o: usize| u32::from_le_bytes(data[o..o + 4].try_into().unwrap());
    assert!(&data[..4] == b"glTF" && u(4) == 2 && u(8) as usize == data.len() && &data[16..20] == b"JSON", "GLB header {:?}", &data[..20]);
    let g = Glb::read(&glb);
    let nodes = g.json["nodes"].as_array().unwrap();
    let top = nodes.iter().find(|n| n["name"] == "ink").unwrap();
    let parts: HashMap<&str, u64> = top["children"].as_array().unwrap().iter().map(|k| (nodes[acc(k)]["name"].as_str().unwrap(), k.as_u64().unwrap())).collect();
    let mut names: Vec<_> = parts.keys().copied().collect();
    names.sort();
    let meshes = g.json["meshes"].as_array().unwrap();
    assert!(names == ["root/crank", "root/frame", "root/rod", "root/slider"] && meshes.len() == 4, "a node per instance, a mesh per body: {names:?}, {} meshes", meshes.len());
    for m in meshes {
        for pr in m["primitives"].as_array().unwrap() {
            let nv = g.json["accessors"][acc(&pr["attributes"]["POSITION"])]["count"].as_u64().unwrap() as u32;
            let ix = g.indices(acc(&pr["indices"]));
            assert!(ix.len().is_multiple_of(3) && ix.iter().all(|&i| i < nv), "{}: indices in range", m["name"]);
        }
    }
    // The file's reference values: the slider pin at x = 90 (t = 0) and 50
    // (t = 0.5); the crank turned a quarter at t = 0.25.
    let track = tracks(&g);
    let slider = &track[&(parts["root/slider"], "translation".to_string())];
    let times: Vec<f32> = g.floats(acc(&slider["input"])).iter().map(|t| t[0]).collect();
    let at = |t: f32| (0..times.len()).min_by(|&a, &b| (times[a] - t).abs().total_cmp(&(times[b] - t).abs())).unwrap();
    let slide = g.floats(acc(&slider["output"]));
    assert!((slide[at(0.0)][0] - 90.0).abs() < 1e-3 && (slide[at(0.5)][0] - 50.0).abs() < 1e-3, "slider at 0 and 0.5: {:?} {:?}", slide[at(0.0)], slide[at(0.5)]);
    let q = &g.floats(acc(&track[&(parts["root/crank"], "rotation".to_string())]["output"]))[at(0.25)];
    let h = 0.5f32.sqrt();
    assert!((q[2].abs() - h).abs() < 1e-3 && (q[3].abs() - h).abs() < 1e-3, "crank a quarter turn about z at 0.25: {q:?}");
    assert!(!track.contains_key(&(parts["root/frame"], "translation".to_string())), "the frame does not move");
    if !no_cache() {
        let again = s.path("again.glb");
        let d2 = ink(["mesh", &crank, "--res", "1", "--out", &again]);
        assert!(d2["meshes_made"] == 0 && std::fs::read(&again).unwrap() == data, "every mesh from the cache, the same file: {}", pretty(&d2));
    }
    _ = d;
}

#[test]
fn a_finish_varying_over_a_body_is_a_few_primitives_not_one_per_value() {
    let s = Scratch::new("glb-rough");
    let rough = s.model("rough.ink", "(import \"std\")\n(body ball (shape (sphere 10)) (attr roughness (/ (+ x 10) 20)))\n(assembly root (instance b ball))\n");
    let glb = s.path("rough.glb");
    ink(["mesh", &rough, "--out", &glb]);
    let g = Glb::read(&glb);
    let prims: usize = g.json["meshes"].as_array().unwrap().iter().map(|m| m["primitives"].as_array().unwrap().len()).sum();
    assert!((2..=17).contains(&prims), "roughness 0 to 1 over a ball is at most 17 primitives (to 1/16): {prims}");
}

#[test]
fn a_clear_body_is_a_transmission_material_and_an_opaque_one_has_no_extension() {
    let s = Scratch::new("glb-clear");
    let glb = |name: &str, attr: &str| {
        let m = s.model(&format!("{name}.ink"), &format!("(import \"std\")\n(body box (shape (box 10 10 10)) {attr})\n(assembly root (instance b box))\n"));
        let out = s.path(&format!("{name}.glb"));
        ink(["mesh", &m, "--res", "1", "--out", &out]);
        Glb::read(&out).json
    };
    let tr = |j: &Value| -> Vec<f64> {
        j["materials"].as_array().unwrap().iter().map(|m| m["extensions"]["KHR_materials_transmission"]["transmissionFactor"].as_f64().unwrap_or(0.0)).collect()
    };
    let opaque = glb("opaque", "");
    assert!(opaque.get("extensionsUsed").is_none() && tr(&opaque) == [0.0], "opaque: no extension: {}", opaque["materials"]);
    // 0.8 to the nearest 1/16.
    let clear = glb("clear", "(attr transmission 0.8)");
    assert!(clear["extensionsUsed"] == serde_json::json!(["KHR_materials_transmission"]) && tr(&clear) == [0.8125], "clear: {} {}", clear["extensionsUsed"], clear["materials"]);
    // Clear where x > 0, opaque elsewhere: two primitives, one clear.
    let half = glb("half", "(attr transmission (select x 0 0.8))");
    let mut t = tr(&half);
    t.sort_by(f64::total_cmp);
    assert!(t == [0.0, 0.8125] && half["meshes"][0]["primitives"].as_array().unwrap().len() == 2, "half clear: {}", half["materials"]);
}

#[test]
fn max_res_meshes_each_part_no_coarser_and_a_small_one_at_its_own_finer_res() {
    // The big plate is at --max-res (not halved for its thin fin, as with
    // --res), the small ball at its own finer res. (--max-res 1, not 0.5:
    // the plate's own res is 1.19, and at 0.5 it cost 14 CPU s.)
    let s = Scratch::new("glb-max-res");
    let sizes = s.model("sizes.ink", "(import \"std\")\n(body big (shape (min (box 200 200 20) (move (box 0.3 50 10) 0 0 15))))\n(body small (shape (move (sphere 2) 0 0 30)))\n(assembly root (instance a big) (instance b small))\n");
    let d = ink(["mesh", &sizes, "--out", &s.path("sizes.glb"), "--max-res", "1"]);
    let r: HashMap<&str, f64> = d["objects"].as_array().unwrap().iter().map(|o| (o["target"].as_str().unwrap(), f(o, "res"))).collect();
    assert!(r["root/a"] == 1.0 && r["root/b"] < 0.1, "{r:?}");
}

/// The res ladder: powers of 2^(1/4).
fn on_ladder(r: f64) -> bool {
    let l = 4.0 * r.log2();
    (l - l.round()).abs() < 1e-6
}

/// `ink ARGS` with XDG_CACHE_HOME at `cache`: its answer, which must be ok.
fn ink_in<S: AsRef<std::ffi::OsStr>>(s: &Scratch, args: impl IntoIterator<Item = S>, cache: &str) -> Value {
    s.ink(args).env("XDG_CACHE_HOME", cache).ok()
}

#[test]
fn a_thin_sheet_far_from_a_speck_is_drawn_at_a_res_on_the_ladder() {
    // A 40 x 30 x 0.5 mm sheet with a speck 500 mm off: its own res, by its
    // box, is 2.83, where at 9 of 12 offsets it had no triangles.
    let s = Scratch::new("glb-sheet");
    for dz in ["0.3", "0.6", "1.1"] {
        let sheet = s.model(&format!("sheet{dz}.ink"), &format!("(import \"std\")\n(body s (shape (min (move (box 40 30 0.5) 0 0 {dz}) (move (sphere 1) 300 300 300))))\n"));
        let d = ink(["mesh", &sheet, "--out", &s.path("sheet.glb")]);
        let o = &d["objects"][0];
        assert!(f(o, "mesh_area") > 2400.0 && on_ladder(f(o, "res")), "{dz} mm up: {}", pretty(o));
    }
}

#[test]
fn a_face_of_no_thickness_or_a_gap_that_thins_to_nothing_is_not_drawn_finer() {
    // Beside a cube, the shell of another (a subtract of it from itself: a
    // face of no thickness, which no res draws): not made again for it.
    let s = Scratch::new("glb-shell");
    let shell = s.model("shell.ink", "(import \"std\")\n(def ghost (move (box 10 10 10) 30 0 0))\n(body b (shape (min (box 20 20 20) (subtract ghost ghost))))\n");
    let d = ink_in(&s, ["mesh", &shell, "--out", &s.path("shell.glb")], &s.path("cache"));
    assert!(d["meshes_made"] == 1 && on_ladder(f(&d["objects"][0], "res")), "{}", pretty(&d));
    // A ball resting in a pocket 0.1 mm larger: the gap thins to nothing
    // where they touch (a bearing's), and no res resolves it.
    let cup = s.model("cup.ink", "(import \"std\")\n(body cup (shape (min (subtract (box 10 10 6) (move (sphere 2.1) 0 0 0.1)) (sphere 2))))\n");
    let d = ink(["mesh", &cup, "--out", &s.path("cup.glb")]);
    let o = &d["objects"][0];
    assert!(f(o, "res") > 0.07 && f(&o["unresolved"], "cells") > 0.0, "kept at its own res, its unresolved cells said: {}", pretty(o));
}

#[test]
fn a_toothed_belt_looped_round_a_large_box_is_drawn_whole() {
    // A 0.63 mm back with 0.75 mm teeth, looped 300 mm across: by its box
    // its res is 1.41, where its back falls between samples (it was drawn
    // in 21 pieces, none reported). By its surface's area it is one piece.
    let s = Scratch::new("glb-belt");
    let belt = s.model(
        "belt.ink",
        "(import \"std\")\n(def back (subtract (cylinder-z 150 6) (cylinder-z 149.37 7)))\n\
         (body belt (shape (min back (repeat-polar (move (box 0.75 1.2 6) 148.995 0 0) 471))))\n\
         (body block (shape (box 400 400 20)))\n(assembly root (instance a belt (at (move 0 0 30))) (instance b block))\n",
    );
    let d = ink(["mesh", &belt, "--out", &s.path("belt.glb")]);
    let o = d["objects"].as_array().unwrap().iter().find(|o| o["target"] == "root/a").unwrap();
    assert!(o["shells"] == 1, "{}", pretty(o));
}

#[test]
fn a_film_thinner_than_the_budget_affords_is_reported_and_its_glb_stays_within_the_budget() {
    // A 0.02 mm film over a block, one body: drawing it takes res 0.02,
    // which the budget cannot afford (making it finer and finer reached 8
    // GB in 6 s). The GLB samples within the budget and says where the film
    // may be.
    let s = Scratch::new("glb-film");
    let film = s.model("film.ink", "(import \"std\")\n(body b (shape (min (box 30 30 30) (move (box 40 40 0.02) 0 0 40))))\n");
    let d = ink(["mesh", &film, "--out", &s.path("film.glb")]);
    let o = &d["objects"][0];
    let u = &o["unresolved"];
    let warned = d["warnings"].as_array().unwrap().iter().any(|w| w.as_str().unwrap().contains("could not rule out material"));
    assert!(f(o, "uniform_triangles") < 4.6e6 && f(u, "undecided_at_least") > 0.0 && warned, "{}", pretty(&d));
}

#[test]
fn a_res_finer_than_the_budget_pays_for_is_refused_and_reported() {
    // A ball 0.2 across beside a field whose bounds cannot prove it
    // positive, at --res 0.0001: every brick of its box would be sampled,
    // its search sized for that res (it reached 5.4 GB). Made where its
    // grid fits the budget, its search a quarter as large each time.
    let s = Scratch::new("glb-refused");
    let m = s.model("loose.ink", "(import \"std\")\n(body b (shape (min (sphere 0.1) (max (box 1 1 1) (+ 0.000001 (- (* 2 x) (+ x x)))))))\n");
    let d = ink(["mesh", &m, "--res", "0.0001", "--out", &s.path("loose.glb")]);
    let u = &d["objects"][0]["unresolved"];
    assert!(f(u, "refused") == 0.0001 && f(&d["objects"][0], "res") > 0.001 && !u["regions"].as_array().unwrap().is_empty(), "{}", pretty(&d["objects"][0]));
}

#[test]
fn a_parts_mesh_made_in_one_model_is_found_in_another_and_progress_says_so() {
    let s = Scratch::new("glb-progress");
    s.model("ring.ink", "(import \"std\")\n(def band (subtract (cylinder-z 8 4) (cylinder-z 6 5)))\n(body ring (shape (min band (tag spoke (box 14 1 2)))))\n");
    let one = s.model("one.ink", "(import \"ring.ink\" as r)\n(assembly root (instance a r.ring))\n");
    // Its own body and tag before the import: the ring's names are numbered after them.
    let two = s.model(
        "two.ink",
        "(import \"std\")\n(body own (shape (tag lump (sphere 3))))\n(import \"ring.ink\" as r)\n(assembly root (instance o own (at (move 20 0 0))) (instance a r.ring))\n",
    );
    // A cache of its own: in the shared one, a mesh a failing run stored
    // would be found by every run after it.
    let cache = s.path("cache");
    let mut lines = vec![];
    for m in [&one, &two] {
        let progress = format!("{m}.progress");
        s.ink(["mesh", m, "--out", &m.replace(".ink", ".glb")]).env("INK_PROGRESS", &progress).env("XDG_CACHE_HOME", &cache).ok();
        lines.push(std::fs::read_to_string(&progress).unwrap().lines().map(|l| serde_json::from_str::<Value>(l).unwrap()).collect::<Vec<_>>());
    }
    let done: HashMap<String, String> =
        lines[1].iter().filter_map(|l| Some((l["done"].as_str()?.split(" at res ").next().unwrap().to_string(), l["status"].as_str().unwrap_or("").to_string()))).collect();
    let total = lines[1].iter().any(|l| l["of"] == "meshes" && l["total"] == 2);
    let mut parts: Vec<&str> = done.keys().map(String::as_str).collect();
    parts.sort();
    assert!(total && parts == ["own", "r.ring"], "how many meshes, then each as it is made or found (INK_PROGRESS): {:?}", lines[1]);
    if !no_cache() {
        assert!(done["r.ring"] == "cached", "the ring's view mesh, made in one.ink, is found in two.ink, which numbers its names otherwise: {done:?}");
    }
}

#[test]
fn a_body_built_with_arguments_says_which_it_is_an_instance_of() {
    let s = Scratch::new("glb-rods");
    let rods = s.model("rods.ink", "(import \"std\")\n(body (rod (h 10)) (shape (box 2 2 h)))\n(assembly root (instance a rod) (instance b (rod 20) (at (move 5 0 0))))\n");
    let glb = s.path("rods.glb");
    ink(["mesh", &rods, "--out", &glb]);
    let g = Glb::read(&glb);
    let ex: HashMap<&str, &Value> =
        g.json["nodes"].as_array().unwrap().iter().filter_map(|n| Some((n["name"].as_str().filter(|s| s.starts_with("root/"))?, &n["extras"]))).collect();
    assert!(ex["root/a"].get("of").is_none() && ex["root/b"]["of"] == "rod" && ex["root/b"]["body"] != ex["root/a"]["body"], "{ex:?}");
}

#[test]
fn a_chain_on_a_path_is_a_node_per_link_moving_smoothly_round_its_sprockets() {
    // Thirty-two simple links at path-link-at on a chordal pitch path,
    // running 20 mm/s. Each link a node moving smoothly (no jump where it
    // changes piece or crosses the loop's start); link 0 starts at A's top
    // and after 1 s is 20 mm round A's pitch circle, at 16 x 6.35 / 2 pi mm
    // of chain per radian.
    let s = Scratch::new("glb-chain");
    let glb = s.path("chain.glb");
    ink(["mesh", &fixture("chain.ink"), "--res", "0.5", "--out", &glb]);
    let g = Glb::read(&glb);
    let nodes = g.json["nodes"].as_array().unwrap();
    let links: Vec<u64> = (0..nodes.len() as u64).filter(|&k| nodes[k as usize]["name"].as_str().is_some_and(|n| n.starts_with("root/link-"))).collect();
    let nmesh = g.json["meshes"].as_array().unwrap().len();
    assert!(links.len() == 32 && nmesh == 3, "a node per link: {} links, {nmesh} meshes", links.len());
    let track = tracks(&g);
    let dist = |a: &[f32], b: &[f32]| (0..3).map(|i| ((a[i] - b[i]) as f64).powi(2)).sum::<f64>().sqrt();
    let mut step = 0.0f64;
    for &k in &links {
        let out = g.floats(acc(&track[&(k, "translation".to_string())]["output"]));
        step = out.windows(2).map(|w| dist(&w[0], &w[1])).fold(step, f64::max);
    }
    assert!(step < 0.7, "every link moves at most 20 mm/s / 30 between samples: {step}");
    let l0 = nodes.iter().position(|n| n["name"] == "root/link-0").unwrap() as u64;
    let tr = &track[&(l0, "translation".to_string())];
    let times: Vec<f32> = g.floats(acc(&tr["input"])).iter().map(|t| t[0]).collect();
    let i = (0..times.len()).min_by(|&a, &b| (times[a] - 1.0).abs().total_cmp(&(times[b] - 1.0).abs())).unwrap();
    let p0 = &g.floats(acc(&tr["output"]))[i];
    let (rp, rho) = (6.35 / 2.0 / (PI / 16.0).sin(), 16.0 * 6.35 / (2.0 * PI));
    let a = PI / 2.0 + 20.0 / rho;
    let off = ((p0[0] as f64 - rp * a.cos()).powi(2) + (p0[1] as f64 - rp * a.sin()).powi(2)).sqrt();
    assert!(off < 1e-3, "link 0 20 mm round sprocket A at t 1: {p0:?}, {off} off");
}
