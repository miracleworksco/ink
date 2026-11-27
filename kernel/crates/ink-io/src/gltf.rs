//! Binary glTF (GLB): an assembly as the web draws it. One mesh per
//! distinct body (in its own frame, shared by its instances), a node per
//! instance with its placement, and the motion: each moving node's
//! placement sampled over the span (LINEAR: translation lerped, rotation
//! slerped), and a body whose shape itself changes as a flipbook (a mesh
//! per moment, each shown from its moment to the next by STEP scale keys).
//! Millimetres and Z up in the file's nodes; a root node scales to metres
//! and turns Z up to glTF's Y up, so other tools read it the right size and
//! way up.

use ink_kernel::mesh::Mesh;
use std::fmt::Write as _;

/// A placement: world = r p + t (rows of r).
pub type Place = ([[f64; 3]; 3], [f64; 3]);

pub struct GlbMesh<'a> {
    pub name: &'a str,
    pub mesh: &'a Mesh,
    /// One RGB colour (0..1, sRGB) per vertex.
    pub colors: &'a [[f32; 3]],
    /// The surface's unit normal at each triangle's corners (three per
    /// triangle, in order), as the field gives it; empty: from the faces.
    pub normals: &'a [[f32; 3]],
    /// Metallic, roughness and transmission (0..1) per vertex.
    pub looks: &'a [[f32; 3]],
}

pub enum Shape {
    /// An index into the meshes.
    Mesh(usize),
    /// Meshes shown one after another, at the moments in `GlbScene::flip_times`.
    Flipbook(Vec<usize>),
    Empty,
}

pub struct GlbNode<'a> {
    pub name: &'a str,
    pub shape: Shape,
    /// Fixed placement (None: identity), or one per moment of `GlbScene::times`.
    pub place: Option<Place>,
    pub motion: Option<Vec<Place>>,
    /// JSON object written as the node's `extras` ("" for none).
    pub extras: String,
}

pub struct GlbScene<'a> {
    pub meshes: Vec<GlbMesh<'a>>,
    pub nodes: Vec<GlbNode<'a>>,
    /// Seconds from the span's start of each motion sample and each flipbook frame.
    pub times: Vec<f64>,
    pub flip_times: Vec<f64>,
    /// JSON object written as the scene's `extras` ("" for none).
    pub extras: String,
}

/// Crease angle: faces meeting at more than this keep a sharp edge; given
/// normals (exact, so the same wherever the surface is smooth), corners
/// whose normals differ by more than a degree.
const CREASE_COS: f64 = 0.766; // 40 degrees

/// The glTF extension a material lets light through by.
const TRANSMISSION: &str = "KHR_materials_transmission";
const GIVEN_COS: f64 = 0.99985; // 1 degree

#[derive(Default)]
struct Bin {
    data: Vec<u8>,
    views: Vec<String>,
    accessors: Vec<String>,
}

impl Bin {
    /// A buffer view of `bytes` (target: 34962 vertices, 34963 indices, 0 none).
    fn view(&mut self, bytes: &[u8], target: u32) -> usize {
        while !self.data.len().is_multiple_of(4) {
            self.data.push(0);
        }
        let t = if target > 0 { format!(",\"target\":{target}") } else { String::new() };
        self.views.push(format!("{{\"buffer\":0,\"byteOffset\":{},\"byteLength\":{}{t}}}", self.data.len(), bytes.len()));
        self.data.extend_from_slice(bytes);
        self.views.len() - 1
    }

    fn floats(&mut self, v: &[f32], width: usize, target: u32, bounds: bool) -> usize {
        let kind = ["", "SCALAR", "VEC2", "VEC3", "VEC4"][width];
        let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
        let view = self.view(&bytes, target);
        let mut mm = String::new();
        if bounds {
            let (mut lo, mut hi) = (vec![f32::MAX; width], vec![f32::MIN; width]);
            for c in v.chunks(width) {
                for k in 0..width {
                    lo[k] = lo[k].min(c[k]);
                    hi[k] = hi[k].max(c[k]);
                }
            }
            let list = |x: &[f32]| x.iter().map(|f| num(*f as f64)).collect::<Vec<_>>().join(",");
            mm = format!(",\"min\":[{}],\"max\":[{}]", list(&lo), list(&hi));
        }
        self.accessors.push(format!("{{\"bufferView\":{view},\"componentType\":5126,\"count\":{},\"type\":\"{kind}\"{mm}}}", v.len() / width));
        self.accessors.len() - 1
    }

    fn indices(&mut self, v: &[u32]) -> usize {
        let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
        let view = self.view(&bytes, 34963);
        self.accessors.push(format!("{{\"bufferView\":{view},\"componentType\":5125,\"count\":{},\"type\":\"SCALAR\"}}", v.len()));
        self.accessors.len() - 1
    }
}

/// A number as JSON (finite; shortest round trip).
fn num(x: f64) -> String {
    if x.is_finite() { format!("{x}") } else { "0".into() }
}

pub fn json_str(s: &str) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            c if (c as u32) < 0x20 => {
                let _ = write!(o, "\\u{:04x}", c as u32);
            }
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

fn srgb_to_linear(c: f32) -> f32 {
    if c <= 0.04045 { c / 12.92 } else { ((c + 0.055) / 1.055).powf(2.4) }
}

/// Vertices split where faces meet at a crease, each with its normal (the
/// mean of the corners' `normals` around it that agree, else the faces'
/// area weighted on its side), and the source vertex it came from.
fn shade(m: &Mesh, normals: &[[f32; 3]]) -> (Vec<[f32; 3]>, Vec<[f32; 3]>, Vec<u32>, Vec<u32>) {
    let face_n: Vec<[f64; 3]> = m.tris.iter().map(|t| ink_core::v3::normal(t.map(|i| m.verts[i as usize]))).collect();
    let corner_n = |f: usize, c: usize| match normals.get(3 * f + c) {
        Some(n) => n.map(|x| x as f64),
        None => face_n[f],
    };
    let unit = |v: [f64; 3]| {
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if l > 0.0 { v.map(|c| c / l) } else { v }
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    // Faces around each vertex.
    let mut start = vec![0u32; m.verts.len() + 1];
    for t in &m.tris {
        for &i in t {
            start[i as usize + 1] += 1;
        }
    }
    for i in 0..m.verts.len() {
        start[i + 1] += start[i];
    }
    let mut fill = start.clone();
    let mut around = vec![0u32; 3 * m.tris.len()];
    for (f, t) in m.tris.iter().enumerate() {
        for &i in t {
            around[fill[i as usize] as usize] = f as u32;
            fill[i as usize] += 1;
        }
    }
    // Per vertex, its corners in groups by direction: each joins the first
    // group whose summed normal is within the crease angle.
    let (mut pos, mut nrm, mut src) = (Vec::new(), Vec::new(), Vec::new());
    let mut corner = vec![0u32; 3 * m.tris.len()];
    for v in 0..m.verts.len() {
        let mut groups: Vec<([f64; 3], u32)> = Vec::new();
        for &f in &around[start[v] as usize..start[v + 1] as usize] {
            let t = m.tris[f as usize];
            let c = t.iter().position(|&i| i as usize == v).unwrap();
            let n = corner_n(f as usize, c);
            let nu = unit(n);
            let crease = if normals.is_empty() { CREASE_COS } else { GIVEN_COS };
            let g = match groups.iter().position(|(s, _)| dot(unit(*s), nu) >= crease) {
                Some(g) => g,
                None => {
                    groups.push(([0.0; 3], pos.len() as u32 + groups.len() as u32));
                    groups.len() - 1
                }
            };
            groups[g].0 = std::array::from_fn(|k| groups[g].0[k] + n[k]);
            corner[3 * f as usize + c] = groups[g].1;
        }
        for (s, _) in &groups {
            pos.push(m.verts[v].map(|c| c as f32));
            nrm.push(unit(*s).map(|c| c as f32));
            src.push(v as u32);
        }
    }
    (pos, nrm, corner, src)
}

/// Rotation (unit quaternion x, y, z, w) and scale of a placement's r,
/// a mirror taken out as scale x = -1.
fn rotation(r: &[[f64; 3]; 3]) -> ([f64; 4], f64) {
    let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
    let sx = if det < 0.0 { -1.0 } else { 1.0 };
    let m = [[r[0][0] * sx, r[0][1], r[0][2]], [r[1][0] * sx, r[1][1], r[1][2]], [r[2][0] * sx, r[2][1], r[2][2]]];
    let tr = m[0][0] + m[1][1] + m[2][2];
    let q = if tr > 0.0 {
        let s = (tr + 1.0).sqrt() * 2.0;
        [(m[2][1] - m[1][2]) / s, (m[0][2] - m[2][0]) / s, (m[1][0] - m[0][1]) / s, 0.25 * s]
    } else if m[0][0] > m[1][1] && m[0][0] > m[2][2] {
        let s = (1.0 + m[0][0] - m[1][1] - m[2][2]).sqrt() * 2.0;
        [0.25 * s, (m[0][1] + m[1][0]) / s, (m[0][2] + m[2][0]) / s, (m[2][1] - m[1][2]) / s]
    } else if m[1][1] > m[2][2] {
        let s = (1.0 + m[1][1] - m[0][0] - m[2][2]).sqrt() * 2.0;
        [(m[0][1] + m[1][0]) / s, 0.25 * s, (m[1][2] + m[2][1]) / s, (m[0][2] - m[2][0]) / s]
    } else {
        let s = (1.0 + m[2][2] - m[0][0] - m[1][1]).sqrt() * 2.0;
        [(m[0][2] + m[2][0]) / s, (m[1][2] + m[2][1]) / s, 0.25 * s, (m[1][0] - m[0][1]) / s]
    };
    let l = q.iter().map(|c| c * c).sum::<f64>().sqrt();
    (q.map(|c| c / l), sx)
}

/// A glTF matrix (column major) of a placement.
fn matrix(p: &Place) -> String {
    let (r, t) = p;
    let m = [r[0][0], r[1][0], r[2][0], 0.0, r[0][1], r[1][1], r[2][1], 0.0, r[0][2], r[1][2], r[2][2], 0.0, t[0], t[1], t[2], 1.0];
    m.iter().map(|x| num(*x)).collect::<Vec<_>>().join(",")
}

pub fn glb(scene: &GlbScene) -> Vec<u8> {
    let mut bin = Bin::default();
    // Materials: one per distinct look. A mesh has a primitive per
    // metallic, roughness and transmission its triangles have (a
    // triangle's: that of two of its corners, else its first's), each to
    // the nearest 1/16 (so a finish varying smoothly is a few primitives,
    // not one per value); a primitive of one colour has it as its
    // material's, else its vertices carry colours. Transmission above 0 is
    // KHR_materials_transmission's factor.
    let mut materials: Vec<String> = Vec::new();
    let mut gmeshes: Vec<String> = Vec::new();
    for m in &scene.meshes {
        let (pos, nrm, idx, src) = shade(m.mesh, m.normals);
        if pos.is_empty() {
            gmeshes.push(format!("{{\"name\":{},\"primitives\":[]}}", json_str(m.name)));
            continue;
        }
        let look = |v: u32| m.looks.get(src[v as usize] as usize).copied().unwrap_or([0.0, 0.4, 0.0]).map(|x| (x * 16.0).round() / 16.0);
        let mut groups: Vec<([f32; 3], Vec<u32>)> = Vec::new();
        let mut group_of: std::collections::HashMap<[u32; 3], usize> = std::collections::HashMap::new();
        for t in idx.chunks(3) {
            let [a, b, c] = [t[0], t[1], t[2]].map(look);
            let l = if b == c { b } else { a };
            let g = *group_of.entry(l.map(f32::to_bits)).or_insert_with(|| {
                groups.push((l, Vec::new()));
                groups.len() - 1
            });
            groups[g].1.extend_from_slice(t);
        }
        let pa = bin.floats(&pos.concat(), 3, 34962, true);
        let na = bin.floats(&nrm.concat(), 3, 34962, false);
        let mut colors_at = None;
        let mut prims = Vec::new();
        for (l, ix) in &groups {
            let colour = |v: &u32| m.colors.get(src[*v as usize] as usize).copied().unwrap_or([1.0; 3]);
            let uniform = ix.first().map(colour).filter(|c| ix.iter().all(|v| colour(v) == *c));
            let base = uniform.unwrap_or([1.0; 3]).map(srgb_to_linear);
            let clear = match l[2] > 0.0 {
                true => format!(",\"extensions\":{{\"{TRANSMISSION}\":{{\"transmissionFactor\":{}}}}}", num(l[2] as f64)),
                false => String::new(),
            };
            let mat = format!(
                "{{\"pbrMetallicRoughness\":{{\"baseColorFactor\":[{},{},{},1],\"metallicFactor\":{},\"roughnessFactor\":{}}}{clear}}}",
                num(base[0] as f64),
                num(base[1] as f64),
                num(base[2] as f64),
                num(l[0] as f64),
                num(l[1] as f64)
            );
            let mi = materials.iter().position(|x| *x == mat).unwrap_or_else(|| {
                materials.push(mat);
                materials.len() - 1
            });
            let mut attrs = format!("\"POSITION\":{pa},\"NORMAL\":{na}");
            if uniform.is_none() {
                let ca = *colors_at.get_or_insert_with(|| {
                    let cols: Vec<f32> = src.iter().flat_map(|&v| m.colors.get(v as usize).copied().unwrap_or([1.0; 3]).map(srgb_to_linear)).collect();
                    bin.floats(&cols, 3, 34962, false)
                });
                let _ = write!(attrs, ",\"COLOR_0\":{ca}");
            }
            let ia = bin.indices(ix);
            prims.push(format!("{{\"attributes\":{{{attrs}}},\"indices\":{ia},\"material\":{mi}}}"));
        }
        gmeshes.push(format!("{{\"name\":{},\"primitives\":[{}]}}", json_str(m.name), prims.join(",")));
    }

    // Nodes: 0 is the root (mm, Z up -> m, Y up); then each node, and a
    // flipbook's frames after it as its children.
    let mut nodes: Vec<String> = vec![String::new()];
    let mut roots: Vec<usize> = Vec::new();
    let mut channels: Vec<String> = Vec::new();
    let mut samplers: Vec<String> = Vec::new();
    let mut tracks: std::collections::HashMap<Vec<u8>, usize> = std::collections::HashMap::new();
    let time_in = (!scene.times.is_empty()).then(|| bin.floats(&scene.times.iter().map(|&t| t as f32).collect::<Vec<_>>(), 1, 0, true));
    // A track's values, shared by every node that moves the same way.
    let mut sampler = |bin: &mut Bin, input: usize, v: Vec<f32>, width: usize, interp: &str| -> usize {
        let key: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).chain((input as u32).to_le_bytes()).chain([width as u8]).collect();
        let out = *tracks.entry(key).or_insert_with(|| bin.floats(&v, width, 0, false));
        samplers.push(format!("{{\"input\":{input},\"output\":{out},\"interpolation\":\"{interp}\"}}"));
        samplers.len() - 1
    };
    for n in &scene.nodes {
        let k = nodes.len();
        roots.push(k);
        nodes.push(String::new());
        let mut fields = vec![format!("\"name\":{}", json_str(n.name))];
        if !n.extras.is_empty() {
            fields.push(format!("\"extras\":{}", n.extras));
        }
        match (&n.motion, time_in) {
            (Some(places), Some(input)) => {
                let (mut tr, mut rot) = (Vec::new(), Vec::new());
                let mut prev = [0.0, 0.0, 0.0, 1.0];
                let mut sx = 1.0;
                for (r, t) in places {
                    let (mut q, s) = rotation(r);
                    // The shorter way round from the last sample.
                    if q.iter().zip(&prev).map(|(a, b)| a * b).sum::<f64>() < 0.0 {
                        q = q.map(|c| -c);
                    }
                    prev = q;
                    sx = s;
                    tr.extend(t.map(|c| c as f32));
                    rot.extend(q.map(|c| c as f32));
                }
                let q0 = &rot[0..4];
                fields.push(format!("\"translation\":[{},{},{}]", tr[0], tr[1], tr[2]));
                fields.push(format!("\"rotation\":[{},{},{},{}]", q0[0], q0[1], q0[2], q0[3]));
                if sx < 0.0 {
                    fields.push("\"scale\":[-1,1,1]".into());
                }
                let moves = |v: &[f32], w: usize| v.chunks(w).any(|c| c != &v[..w]);
                if moves(&tr, 3) {
                    let s = sampler(&mut bin, input, tr, 3, "LINEAR");
                    channels.push(format!("{{\"sampler\":{s},\"target\":{{\"node\":{k},\"path\":\"translation\"}}}}"));
                }
                if moves(&rot, 4) {
                    let s = sampler(&mut bin, input, rot, 4, "LINEAR");
                    channels.push(format!("{{\"sampler\":{s},\"target\":{{\"node\":{k},\"path\":\"rotation\"}}}}"));
                }
            }
            _ => {
                if let Some(p) = &n.place {
                    fields.push(format!("\"matrix\":[{}]", matrix(p)));
                }
            }
        }
        match &n.shape {
            Shape::Mesh(m) => fields.push(format!("\"mesh\":{m}")),
            Shape::Flipbook(frames) => {
                // Frame i is scaled 1 from its moment to the next, else 0.
                let ft = &scene.flip_times;
                let mut kids = Vec::new();
                for (i, &m) in frames.iter().enumerate() {
                    let c = nodes.len();
                    kids.push(c.to_string());
                    let shown = if i == 0 { "" } else { ",\"scale\":[0,0,0]" };
                    nodes.push(format!("{{\"name\":{},\"mesh\":{m}{shown}}}", json_str(&format!("{} {i}", n.name))));
                    let mut keys: Vec<(f64, f32)> = Vec::new();
                    if i > 0 {
                        keys.push((ft[0], 0.0));
                    }
                    keys.push((ft[i], 1.0));
                    if i + 1 < frames.len() {
                        keys.push((ft[i + 1], 0.0));
                    }
                    let times: Vec<f32> = keys.iter().map(|k| k.0 as f32).collect();
                    let input = bin.floats(&times, 1, 0, true);
                    let vals: Vec<f32> = keys.iter().flat_map(|k| [k.1; 3]).collect();
                    let s = sampler(&mut bin, input, vals, 3, "STEP");
                    channels.push(format!("{{\"sampler\":{s},\"target\":{{\"node\":{c},\"path\":\"scale\"}}}}"));
                }
                fields.push(format!("\"children\":[{}]", kids.join(",")));
            }
            Shape::Empty => {}
        }
        nodes[k] = format!("{{{}}}", fields.join(","));
    }
    let up = "0.001,0,0,0, 0,0,-0.001,0, 0,0.001,0,0, 0,0,0,1";
    nodes[0] = format!("{{\"name\":\"ink\",\"matrix\":[{}],\"children\":[{}]}}", up.replace(' ', ""), roots.iter().map(|r| r.to_string()).collect::<Vec<_>>().join(","));

    let mut j = String::new();
    let _ = write!(j, "{{\"asset\":{{\"version\":\"2.0\",\"generator\":\"ink\"}},\"scene\":0,\"scenes\":[{{\"nodes\":[0]");
    if !scene.extras.is_empty() {
        let _ = write!(j, ",\"extras\":{}", scene.extras);
    }
    let _ = write!(j, "}}],\"nodes\":[{}]", nodes.join(","));
    if materials.iter().any(|m| m.contains(TRANSMISSION)) {
        let _ = write!(j, ",\"extensionsUsed\":[\"{TRANSMISSION}\"]");
    }
    let _ = write!(j, ",\"meshes\":[{}],\"materials\":[{}]", gmeshes.join(","), materials.join(","));
    if !channels.is_empty() {
        let _ = write!(j, ",\"animations\":[{{\"name\":\"motion\",\"channels\":[{}],\"samplers\":[{}]}}]", channels.join(","), samplers.join(","));
    }
    while bin.data.len() % 4 != 0 {
        bin.data.push(0);
    }
    let _ = write!(j, ",\"accessors\":[{}],\"bufferViews\":[{}]", bin.accessors.join(","), bin.views.join(","));
    let _ = write!(j, ",\"buffers\":[{{\"byteLength\":{}}}]}}", bin.data.len());
    let mut jb = j.into_bytes();
    while !jb.len().is_multiple_of(4) {
        jb.push(b' ');
    }
    let mut out = Vec::with_capacity(28 + jb.len() + bin.data.len());
    out.extend(b"glTF");
    out.extend(2u32.to_le_bytes());
    out.extend(((12 + 8 + jb.len() + 8 + bin.data.len()) as u32).to_le_bytes());
    out.extend((jb.len() as u32).to_le_bytes());
    out.extend(b"JSON");
    out.extend(jb);
    out.extend((bin.data.len() as u32).to_le_bytes());
    out.extend(b"BIN\0");
    out.extend(bin.data);
    out
}
