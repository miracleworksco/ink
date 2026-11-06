//! Orthographic renderer: labelled views of meshes, so a multimodal agent
//! can glance at a result. Surfaces are shaded by their normals, smooth
//! across a smooth surface's facets and turning at a crease (`CREASE`).
//! Outlines mark silhouettes, steps in depth, creases and where one object
//! meets another, so a view straight down an open tray shows its pockets
//! though every face it sees is lit alike.

use ink_core::V3;
use ink_core::v3::{add, cross, dot, norm, normal, scale as times, unit};
use ink_kernel::mesh::Mesh;
use rayon::prelude::*;

/// Where a view looks from.
#[derive(Clone, Debug)]
pub struct Camera {
    pub name: String,
    /// Direction from the model toward the camera.
    pub eye: V3,
    /// The world direction drawn upward (made square to `eye`).
    pub up: V3,
}

/// The named views: iso, and each face's (front looks from -y; bottom
/// keeps x to the right). The first four are the sheet.
pub const NAMED: [(&str, V3, V3); 7] = [
    ("iso", [1.0, -1.0, 0.8], [0.0, 0.0, 1.0]),
    ("top", [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]),
    ("front", [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]),
    ("right", [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    ("back", [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
    ("left", [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    ("bottom", [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]),
];

impl Camera {
    pub fn named(name: &str) -> Option<Camera> {
        NAMED.iter().find(|v| v.0 == name).map(|&(n, eye, up)| Camera { name: n.into(), eye: unit(eye), up })
    }

    /// Looking from direction `from`, z up on screen (y when it looks
    /// along z), named by where it looks from.
    pub fn looking_from(from: V3) -> Option<Camera> {
        let l = norm(from);
        if !(l > 0.0 && l.is_finite()) {
            return None;
        }
        let eye = times(from, 1.0 / l);
        let up = if eye[2].abs() > 0.999 { [0.0, eye[2].signum(), 0.0] } else { [0.0, 0.0, 1.0] };
        Some(Camera { name: format!("from {}", from_text(eye)), eye, up })
    }

    /// iso, top, front and right.
    pub fn sheet() -> Vec<Camera> {
        NAMED[..4].iter().map(|v| Camera::named(v.0).unwrap()).collect()
    }
}

/// The signs of a direction's clear components: "+x -y +z".
fn from_text(eye: V3) -> String {
    let parts: Vec<String> = (0..3).filter(|&k| eye[k].abs() > 0.05).map(|k| format!("{}{}", if eye[k] > 0.0 { '+' } else { '-' }, ['x', 'y', 'z'][k])).collect();
    parts.join(" ")
}

/// A mesh in world coordinates and how to draw it. Where two objects'
/// surfaces are closer than their errors (two bodies sharing a face), their
/// order is unknown, and the object drawn first shows: deciding by rounding
/// would alternate between them pixel by pixel.
#[derive(Clone, Copy)]
pub struct Object<'a> {
    pub mesh: &'a Mesh,
    /// Each vertex's colour.
    pub colors: &'a [[f32; 3]],
    /// The surface's unit normal at each triangle's corners, three a
    /// triangle (the field's: seen::corner_normals); empty: the mesh's own
    /// (`corner_normals`).
    pub normals: &'a [[f32; 3]],
    /// Each vertex's transmission; empty: opaque.
    pub transmission: &'a [f32],
    /// How far each vertex may lie from the surface it stands for; empty:
    /// none.
    pub error: &'a [f32],
}

impl<'a> Object<'a> {
    /// Opaque, exact, shaded by the mesh's own normals.
    pub fn plain(mesh: &'a Mesh, colors: &'a [[f32; 3]]) -> Object<'a> {
        Object { mesh, colors, normals: &[], transmission: &[], error: &[] }
    }
}

/// Triangles whose normals differ by more than this angle (its cosine, 30
/// degrees) meet at a crease: shading turns there, and an outline marks it.
/// A smooth surface's facets differ by far less: a mesh simplified to its
/// tolerance keeps them within 20 degrees (seen::SMOOTH_COS).
const CREASE: f64 = 0.866;

/// The mesh's own normal at each triangle's corners, for a mesh the field
/// does not give normals for (a solution's deformed shape): the
/// area-weighted mean of the normals of the triangles round the corner's
/// vertex within CREASE of the triangle's own. Across a smooth surface's
/// facets shading is then continuous; at a crease each side keeps its own.
pub fn corner_normals(m: &Mesh) -> Vec<[V3; 3]> {
    // Each triangle's normal, as long as twice its area.
    let faces: Vec<V3> = m.tris.iter().map(|t| normal(t.map(|i| m.verts[i as usize]))).collect();
    let mut round = vec![Vec::new(); m.verts.len()];
    for (f, t) in m.tris.iter().enumerate() {
        t.iter().for_each(|&i| round[i as usize].push(f as u32));
    }
    m.tris
        .par_iter()
        .zip(&faces)
        .map(|(t, &face)| {
            let own = unit(face);
            t.map(|i| unit(round[i as usize].iter().map(|&g| faces[g as usize]).filter(|&n| dot(unit(n), own) >= CREASE).fold([0.0; 3], add)))
        })
        .collect()
}

/// A surface seen more nearly edge-on than this (the cosine of its normal
/// with the view) is taken as at this: steps in depth across it, and its
/// error in depth, at most 10 times what they are face-on.
const EDGE_ON: f64 = 0.1;

/// Depths computed from coordinates of magnitude m agree to m times this:
/// far above their rounding, far below a part's size.
const ROUNDING: f64 = 1e-9;

pub struct View {
    pub name: String,
    pub eye: V3,
    pub right: V3,
    pub up: V3,
    pub mm_per_px: f64,
    /// Top-left pixel of this view in the image.
    pub origin: [usize; 2],
}

/// The views in a grid, row by row (one view: the whole image), the image
/// `size` pixels across. Transmission composites surface layers in depth
/// order; it does not model refraction. `section`
/// (axis, at): faces lying on that plane, as a cut leaves them, are hatched.
pub fn render(objects: &[Object], cameras: &[Camera], size: usize, section: Option<(usize, f64)>) -> (Vec<u8>, usize, usize, Vec<View>) {
    let (cols, rows) = grid(cameras.len());
    let panel = size / cols;
    let (w, h) = (panel * cols, panel * rows);
    let drawn: Vec<(Vec<u8>, View)> = cameras.par_iter().enumerate().map(|(i, c)| draw(objects, c, panel, [i % cols * panel, i / cols * panel], section)).collect();
    let mut img = vec![245u8; w * h * 3];
    for (p, v) in &drawn {
        for y in 0..panel {
            let dst = ((v.origin[1] + y) * w + v.origin[0]) * 3;
            img[dst..dst + panel * 3].copy_from_slice(&p[y * panel * 3..(y + 1) * panel * 3]);
        }
    }
    // Panel separators.
    for c in 1..cols {
        for x in [c * panel - 1, c * panel] {
            (0..h).for_each(|y| img[(y * w + x) * 3..][..3].copy_from_slice(&[150; 3]));
        }
    }
    for r in 1..rows {
        for y in [r * panel - 1, r * panel] {
            (0..w).for_each(|x| img[(y * w + x) * 3..][..3].copy_from_slice(&[150; 3]));
        }
    }
    // Each view's name and where it looks from, top left.
    let k = (panel / 256).max(1);
    for (_, v) in &drawn {
        let label = match NAMED.iter().any(|n| n.0 == v.name) {
            true => format!("{}: from {}", v.name, from_text(v.eye)),
            false => v.name.clone(),
        };
        text(&mut img, w, [v.origin[0] + 4 * k, v.origin[1] + 4 * k], v.origin[0] + panel - 2, &label, k, [70; 3]);
    }
    (img, w, h, drawn.into_iter().map(|d| d.1).collect())
}

/// `img` (`w` wide) with a colour scale under it: `label` (the quantity and
/// its unit), then `lo`, the colours from `lo` to `hi` (`color` of 0 to 1),
/// `hi`, and `note`. The image and its height.
pub fn with_scale(mut img: Vec<u8>, w: usize, label: &str, (lo, hi): (f64, f64), color: &dyn Fn(f64) -> [f32; 3], note: &str) -> (Vec<u8>, usize) {
    let k = (w / 512).max(1);
    let (top, strip) = (img.len() / 3 / w, 16 * k);
    img.extend(std::iter::repeat_n(245, w * strip * 3));
    (0..w).for_each(|x| img[(top * w + x) * 3..][..3].copy_from_slice(&[150; 3]));
    let (y, ink) = (top + 5 * k, [70; 3]);
    let width = |s: &str| s.chars().count() * 6 * k;
    let (lo, hi) = (significant(lo), significant(hi));
    let mut x = 4 * k;
    text(&mut img, w, [x, y], w, label, k, ink);
    x += width(label) + 6 * k;
    text(&mut img, w, [x, y], w, &lo, k, ink);
    x += width(&lo) + k;
    let bar = w.saturating_sub(x + width(&hi) + width(note) + 24 * k).max(16 * k);
    for i in 0..bar {
        let c = color(i as f64 / (bar - 1) as f64).map(|v| (v.clamp(0.0, 1.0) * 255.0) as u8);
        for yy in top + 4 * k..top + 12 * k {
            if x + i < w {
                img[(yy * w + x + i) * 3..][..3].copy_from_slice(&c);
            }
        }
    }
    x += bar + 2 * k;
    text(&mut img, w, [x, y], w, &hi, k, ink);
    text(&mut img, w, [x + width(&hi) + 18 * k, y], w, note, k, ink);
    (img, top + strip)
}

/// `v` to 3 significant digits, no exponent: 34.7, 0, 1250, 0.0123.
fn significant(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v}");
    }
    let places = (2 - v.abs().log10().floor() as i32).max(0) as usize;
    format!("{v:.places$}")
}

/// mm per pixel of the finest of `render`'s views of the box lo..hi in an
/// image `size` pixels across (or less: what the box holds may span less).
pub fn mm_per_px(lo: V3, hi: V3, cameras: &[Camera], size: usize) -> f64 {
    let panel = (size / grid(cameras.len()).0) as f64;
    let corners: Vec<V3> = (0..8).map(|k| std::array::from_fn(|a| if k >> a & 1 == 0 { lo[a] } else { hi[a] })).collect();
    let per_px = |c: &Camera| {
        let (right, upv) = screen(c.eye, c.up);
        let (x0, x1, y0, y1) = extent(&corners, right, upv);
        (x1 - x0).max(y1 - y0).max(1e-9) / (FILL * panel)
    };
    cameras.iter().map(per_px).fold(f64::INFINITY, f64::min)
}

/// The columns and rows of the grid of n views.
fn grid(n: usize) -> (usize, usize) {
    let cols = (n as f64).sqrt().ceil() as usize;
    (cols, n.div_ceil(cols))
}

/// The share of a panel the model's extent fills.
const FILL: f64 = 0.9;

/// A view's directions on screen: right, and up (made square to `eye`).
fn screen(eye: V3, up: V3) -> (V3, V3) {
    let f = [-eye[0], -eye[1], -eye[2]];
    let right = unit(cross(f, up));
    (right, cross(right, f))
}

/// Points' extent on screen: (x0, x1, y0, y1).
fn extent<'a>(pts: impl IntoIterator<Item = &'a V3>, right: V3, upv: V3) -> (f64, f64, f64, f64) {
    pts.into_iter().fold((f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY), |(x0, x1, y0, y1), v| {
        let (sx, sy) = (dot(*v, right), dot(*v, upv));
        (x0.min(sx), x1.max(sx), y0.min(sy), y1.max(sy))
    })
}

type Basis = ((f64, f64, f64, f64), f64, V3, V3);

fn draw_basis(objects: &[Object], eye: V3, up: V3, panel: usize) -> Basis {
    let (right, upv) = screen(eye, up);
    let (x0, x1, y0, y1) = extent(objects.iter().flat_map(|o| &o.mesh.verts), right, upv);
    let span = (x1 - x0).max(y1 - y0).max(1e-9);
    let scale = panel as f64 * FILL / span;
    ((x0, x1, y0, y1), scale, right, upv)
}

/// No object at a pixel.
const EMPTY: u32 = u32::MAX;

fn draw(objects: &[Object], cam: &Camera, panel: usize, origin: [usize; 2], section: Option<(usize, f64)>) -> (Vec<u8>, View) {
    let eye = cam.eye;
    let ((x0, x1, y0, y1), scale, right, upv) = draw_basis(objects, eye, cam.up, panel);
    let (cx, cy) = ((x0 + x1) / 2.0, (y0 + y1) / 2.0);
    let half = panel as f64 / 2.0;
    let light = unit([eye[0] + 0.4 * upv[0] - 0.3 * right[0], eye[1] + 0.4 * upv[1] - 0.3 * right[1], eye[2] + 0.4 * upv[2] - 0.3 * right[2]]);
    let mut img = vec![245u8; panel * panel * 3];
    let mut zbuf = vec![f64::NEG_INFINITY; panel * panel];
    // The nearest opaque surface's normal (toward the camera) as shaded,
    // its triangle's, its object, and whether it lies on the section plane.
    let mut nbuf = vec![[0.0; 3]; panel * panel];
    let mut fbuf = vec![[0.0; 3]; panel * panel];
    let mut ibuf = vec![EMPTY; panel * panel];
    let mut cut = vec![false; panel * panel];
    // How far that surface may lie from where it is drawn.
    let mut ebuf = vec![0.0; panel * panel];
    let rounding = ROUNDING * objects.iter().flat_map(|o| &o.mesh.verts).fold(0.0, |m: f64, v| m.max(v[0].abs()).max(v[1].abs()).max(v[2].abs()));
    let on_plane = 0.01 / scale;
    struct Fragment {
        pixel: usize,
        depth: f64,
        color: [f32; 3],
        transmission: f32,
    }
    let mut fragments = Vec::new();
    for (o, ob) in objects.iter().enumerate() {
        let (m, colors, transmission, error) = (ob.mesh, ob.colors, ob.transmission, ob.error);
        let corners: Vec<[V3; 3]> = match ob.normals.len() == 3 * m.tris.len() {
            true => ob.normals.chunks(3).map(|c| std::array::from_fn(|k| c[k].map(|x| x as f64))).collect(),
            false => corner_normals(m),
        };
        for (t, corner) in m.tris.iter().zip(&corners) {
            let v = t.map(|i| m.verts[i as usize]);
            let n = normal(v);
            let nl = norm(n);
            if nl == 0.0 {
                continue;
            }
            let u = times(n, if dot(n, eye) < 0.0 { -1.0 / nl } else { 1.0 / nl });
            // The corners' normals on the triangle's side toward the camera
            // (a corner with none: the triangle's).
            let corner = corner.map(|c| {
                if c == [0.0; 3] {
                    u
                } else if dot(c, u) < 0.0 {
                    times(c, -1.0)
                } else {
                    c
                }
            });
            let sectioned = section.is_some_and(|(a, at)| u[a].abs() > 0.999 && v.iter().all(|p| (p[a] - at).abs() <= on_plane));
            let col: [f32; 3] = std::array::from_fn(|k| (colors[t[0] as usize][k] + colors[t[1] as usize][k] + colors[t[2] as usize][k]) / 3.0);
            let tr = t.map(|i| transmission.get(i as usize).copied().unwrap_or(0.0) as f64);
            let err = t.map(|i| error.get(i as usize).copied().unwrap_or(0.0) as f64);
            let px = v.map(|p| [half + (dot(p, right) - cx) * scale, half - (dot(p, upv) - cy) * scale, dot(p, eye)]);
            let bx0 = px.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min).floor().max(0.0) as usize;
            let bx1 = (px.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max).ceil() as usize).min(panel - 1);
            let by0 = px.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min).floor().max(0.0) as usize;
            let by1 = (px.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max).ceil() as usize).min(panel - 1);
            let area = (px[1][0] - px[0][0]) * (px[2][1] - px[0][1]) - (px[2][0] - px[0][0]) * (px[1][1] - px[0][1]);
            if area.abs() < 1e-12 {
                continue;
            }
            // Own a shared edge on exactly one of its triangles. Blending
            // both would make a transparent face darker along every seam.
            let owns = |mut a: [f64; 3], mut b: [f64; 3]| {
                if area < 0.0 {
                    std::mem::swap(&mut a, &mut b);
                }
                b[1] < a[1] || b[1] == a[1] && b[0] > a[0]
            };
            for yy in by0..=by1 {
                for xx in bx0..=bx1 {
                    let (sx, sy) = (xx as f64 + 0.5, yy as f64 + 0.5);
                    let e = |a: [f64; 3], b: [f64; 3]| (b[0] - a[0]) * (sy - a[1]) - (b[1] - a[1]) * (sx - a[0]);
                    let (w0, w1, w2) = (e(px[1], px[2]) / area, e(px[2], px[0]) / area, e(px[0], px[1]) / area);
                    if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 || w0 == 0.0 && !owns(px[1], px[2]) || w1 == 0.0 && !owns(px[2], px[0]) || w2 == 0.0 && !owns(px[0], px[1]) {
                        continue;
                    }
                    let depth = w0 * px[0][2] + w1 * px[1][2] + w2 * px[2][2];
                    let k = yy * panel + xx;
                    if depth > zbuf[k] {
                        let transmission = (w0 * tr[0] + w1 * tr[1] + w2 * tr[2]).clamp(0.0, 1.0);
                        let shaded = unit(add(add(times(corner[0], w0), times(corner[1], w1)), times(corner[2], w2)));
                        let shaded = if shaded == [0.0; 3] { u } else { shaded };
                        let shade = 0.25 + 0.75 * dot(shaded, light).abs();
                        let color = col.map(|c| (c as f64 * shade).clamp(0.0, 1.0) as f32 * 255.0);
                        let e = w0 * err[0] + w1 * err[1] + w2 * err[2];
                        // Another object's surface within both errors, in
                        // depth along the less face-on of the two: the one
                        // drawn first keeps the pixel (Object).
                        let tie = ibuf[k] != EMPTY && ibuf[k] != o as u32 && {
                            let face_on = dot(u, eye).min(dot(fbuf[k], eye)).max(EDGE_ON);
                            depth <= zbuf[k] + (e + ebuf[k]) / face_on + rounding
                        };
                        if transmission == 0.0 && !tie {
                            zbuf[k] = depth;
                            (nbuf[k], fbuf[k], ibuf[k], cut[k], ebuf[k]) = (shaded, u, o as u32, sectioned, e);
                            for c in 0..3 {
                                img[k * 3 + c] = color[c] as u8;
                            }
                        } else if transmission > 0.0 && transmission < 1.0 {
                            fragments.push(Fragment { pixel: k, depth, color, transmission: transmission as f32 });
                        }
                    }
                }
            }
        }
    }
    // Section faces hatched at 45 degrees.
    let pitch = (panel / 64).max(4);
    for k in 0..panel * panel {
        if cut[k] && (k % panel + k / panel).is_multiple_of(pitch) {
            mark(&mut img[k * 3..k * 3 + 3], 0.55);
        }
    }
    outline(&mut img, &zbuf, &nbuf, &fbuf, &ibuf, panel, eye, [times(right, 1.0 / scale), times(upv, -1.0 / scale)]);
    // Opaque surfaces have established the background. Composite only the
    // layers in front of it, farthest first, regardless of object order.
    fragments.sort_unstable_by(|a, b| a.pixel.cmp(&b.pixel).then_with(|| a.depth.total_cmp(&b.depth)));
    for f in fragments {
        if f.depth > zbuf[f.pixel] {
            for c in 0..3 {
                let behind = img[f.pixel * 3 + c] as f32;
                img[f.pixel * 3 + c] = (f.color[c] * (1.0 - f.transmission) + behind * f.transmission) as u8;
            }
        }
    }
    (img, View { name: cam.name.clone(), eye, right, up: upv, mm_per_px: 1.0 / scale, origin })
}

/// Feature lines from the depth and normal buffers (Saito and Takahashi's
/// G-buffer edges): between two neighbouring pixels, a silhouette, another
/// object, a crease (their shading normals, `nbuf`, apart by more than
/// CREASE), or a step in depth that neither pixel's plane (its triangle's,
/// `fbuf`) accounts for by more than 1.5 pixels; drawn dark on the nearer
/// pixel. `step`: the world move of one pixel right and one down.
#[allow(clippy::too_many_arguments)]
fn outline(img: &mut [u8], zbuf: &[f64], nbuf: &[V3], fbuf: &[V3], ibuf: &[u32], panel: usize, eye: V3, step: [V3; 2]) {
    let gap = 1.5 * norm(step[0]);
    // The depth a pixel's plane has one step `m` away.
    let rise = |k: usize, m: V3| -dot(fbuf[k], m) / dot(fbuf[k], eye).max(EDGE_ON);
    let mut line = vec![false; panel * panel];
    for y in 0..panel {
        for x in 0..panel {
            let k = y * panel + x;
            for (j, m) in [(x + 1 < panel).then(|| (k + 1, step[0])), (y + 1 < panel).then(|| (k + panel, step[1]))].into_iter().flatten() {
                if ibuf[k] == EMPTY && ibuf[j] == EMPTY {
                    continue;
                }
                let edge = ibuf[k] != ibuf[j]
                    || dot(nbuf[k], nbuf[j]) < CREASE
                    || (zbuf[j] - zbuf[k] - rise(k, m)).abs() > gap && (zbuf[k] - zbuf[j] - rise(j, times(m, -1.0))).abs() > gap;
                if edge {
                    line[if zbuf[k] >= zbuf[j] { k } else { j }] = true;
                }
            }
        }
    }
    for (k, _) in line.iter().enumerate().filter(|l| *l.1) {
        mark(&mut img[k * 3..k * 3 + 3], 0.7);
    }
}

/// A line drawn over a pixel: `by` of the way to black, on every surface
/// alike, so a line reads as one line along its length.
fn mark(px: &mut [u8], by: f32) {
    px.iter_mut().for_each(|c| *c = (*c as f32 * (1.0 - by)) as u8);
}

/// `s` in the 5x7 font, `k` pixels a dot, its top left at `at`, cut off at
/// column `end`.
fn text(img: &mut [u8], w: usize, at: [usize; 2], end: usize, s: &str, k: usize, color: [u8; 3]) {
    let h = img.len() / 3 / w;
    for (n, ch) in s.chars().enumerate() {
        for (r, row) in glyph(ch.to_ascii_uppercase()).iter().enumerate() {
            for (c, _) in row.bytes().enumerate().filter(|b| b.1 == b'#') {
                let (x0, y0) = (at[0] + (n * 6 + c) * k, at[1] + r * k);
                for y in y0..(y0 + k).min(h) {
                    for x in x0..(x0 + k).min(end.min(w)) {
                        img[(y * w + x) * 3..][..3].copy_from_slice(&color);
                    }
                }
            }
        }
    }
}

/// A character's 5x7 bitmap (upper case; unknown characters blank).
pub fn glyph(c: char) -> [&'static str; 7] {
    match c {
        'A' => [".###.", "#...#", "#...#", "#####", "#...#", "#...#", "#...#"],
        'B' => ["####.", "#...#", "#...#", "####.", "#...#", "#...#", "####."],
        'C' => [".###.", "#...#", "#....", "#....", "#....", "#...#", ".###."],
        'D' => ["####.", "#...#", "#...#", "#...#", "#...#", "#...#", "####."],
        'E' => ["#####", "#....", "#....", "####.", "#....", "#....", "#####"],
        'F' => ["#####", "#....", "#....", "####.", "#....", "#....", "#...."],
        'G' => [".###.", "#...#", "#....", "#.###", "#...#", "#...#", ".####"],
        'H' => ["#...#", "#...#", "#...#", "#####", "#...#", "#...#", "#...#"],
        'I' => [".###.", "..#..", "..#..", "..#..", "..#..", "..#..", ".###."],
        'J' => ["..###", "...#.", "...#.", "...#.", "...#.", "#..#.", ".##.."],
        'K' => ["#...#", "#..#.", "#.#..", "##...", "#.#..", "#..#.", "#...#"],
        'L' => ["#....", "#....", "#....", "#....", "#....", "#....", "#####"],
        'M' => ["#...#", "##.##", "#.#.#", "#.#.#", "#...#", "#...#", "#...#"],
        'N' => ["#...#", "#...#", "##..#", "#.#.#", "#..##", "#...#", "#...#"],
        'O' => [".###.", "#...#", "#...#", "#...#", "#...#", "#...#", ".###."],
        'P' => ["####.", "#...#", "#...#", "####.", "#....", "#....", "#...."],
        'Q' => [".###.", "#...#", "#...#", "#...#", "#.#.#", "#..#.", ".##.#"],
        'R' => ["####.", "#...#", "#...#", "####.", "#.#..", "#..#.", "#...#"],
        'S' => [".####", "#....", "#....", ".###.", "....#", "....#", "####."],
        'T' => ["#####", "..#..", "..#..", "..#..", "..#..", "..#..", "..#.."],
        'U' => ["#...#", "#...#", "#...#", "#...#", "#...#", "#...#", ".###."],
        'V' => ["#...#", "#...#", "#...#", "#...#", "#...#", ".#.#.", "..#.."],
        'W' => ["#...#", "#...#", "#...#", "#.#.#", "#.#.#", "#.#.#", ".#.#."],
        'X' => ["#...#", "#...#", ".#.#.", "..#..", ".#.#.", "#...#", "#...#"],
        'Y' => ["#...#", "#...#", ".#.#.", "..#..", "..#..", "..#..", "..#.."],
        'Z' => ["#####", "....#", "...#.", "..#..", ".#...", "#....", "#####"],
        '0' => [".###.", "#...#", "#..##", "#.#.#", "##..#", "#...#", ".###."],
        '1' => ["..#..", ".##..", "..#..", "..#..", "..#..", "..#..", ".###."],
        '2' => [".###.", "#...#", "....#", "...#.", "..#..", ".#...", "#####"],
        '3' => ["#####", "...#.", "..#..", "...#.", "....#", "#...#", ".###."],
        '4' => ["...#.", "..##.", ".#.#.", "#..#.", "#####", "...#.", "...#."],
        '5' => ["#####", "#....", "####.", "....#", "....#", "#...#", ".###."],
        '6' => ["..##.", ".#...", "#....", "####.", "#...#", "#...#", ".###."],
        '7' => ["#####", "....#", "...#.", "..#..", ".#...", ".#...", ".#..."],
        '8' => [".###.", "#...#", "#...#", ".###.", "#...#", "#...#", ".###."],
        '9' => [".###.", "#...#", "#...#", ".####", "....#", "...#.", ".##.."],
        '.' => [".....", ".....", ".....", ".....", ".....", ".##..", ".##.."],
        ':' => [".....", ".##..", ".##..", ".....", ".##..", ".##..", "....."],
        '(' => ["...#.", "..#..", ".#...", ".#...", ".#...", "..#..", "...#."],
        ')' => [".#...", "..#..", "...#.", "...#.", "...#.", "..#..", ".#..."],
        '-' => [".....", ".....", ".....", "#####", ".....", ".....", "....."],
        '+' => [".....", "..#..", "..#..", "#####", "..#..", "..#..", "....."],
        '/' => [".....", "....#", "...#.", "..#..", ".#...", "#....", "....."],
        '·' => [".....", ".....", ".....", "..#..", ".....", ".....", "....."],
        '°' => [".##..", "#..#.", ".##..", ".....", ".....", ".....", "....."],
        ',' => [".....", ".....", ".....", ".....", ".##..", "..#..", ".#..."],
        '²' => [".##..", "...#.", "..#..", ".###.", ".....", ".....", "....."],
        _ => [".....", ".....", ".....", ".....", ".....", ".....", "....."],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transmission_composites_in_depth_order_once_per_surface() {
        let square = |z| Mesh { verts: vec![[-1.0, -1.0, z], [1.0, -1.0, z], [1.0, 1.0, z], [-1.0, 1.0, z]], tris: vec![[0, 1, 2], [0, 2, 3]] };
        let (back, front, behind) = (square(0.0), square(1.0), square(-1.0));
        let (red, blue, green) = ([[1.0, 0.0, 0.0]; 4], [[0.0, 0.0, 1.0]; 4], [[0.0, 1.0, 0.0]; 4]);
        let glass = [0.78; 4];
        let top = Camera::named("top").unwrap();
        let draw = |objects: &[Object]| super::draw(objects, &top, 64, [0, 0], None).0;
        let opaque = draw(&[Object::plain(&back, &red)]);
        let shown =
            draw(&[Object::plain(&back, &red), Object { transmission: &glass, ..Object::plain(&front, &blue) }, Object { transmission: &glass, ..Object::plain(&behind, &green) }]);
        let reversed =
            draw(&[Object { transmission: &glass, ..Object::plain(&behind, &green) }, Object { transmission: &glass, ..Object::plain(&front, &blue) }, Object::plain(&back, &red)]);
        assert_eq!(shown, reversed, "submission order must not change transparency or opaque occlusion");
        // The centre is on the square's shared triangle edge: one layer,
        // not two. Also check away from it.
        for (x, y) in [(31, 32), (32, 32), (25, 30)] {
            let k = (y * 64 + x) * 3;
            assert!((shown[k] as f64 / opaque[k] as f64 - 0.78).abs() < 0.01, "{:?}", &shown[k..k + 3]);
            assert_eq!(shown[k + 1], 0, "the green layer behind the opaque one is hidden");
            assert!(shown[k + 2] > 0);
        }
        assert_eq!(draw(&[Object::plain(&back, &red), Object { transmission: &[1.0; 4], ..Object::plain(&front, &blue) }]), opaque);
    }

    /// Squares of half-side `c` at height `z`, one mesh.
    fn squares(sq: &[(f64, f64)]) -> Mesh {
        let mut m = Mesh { verts: vec![], tris: vec![] };
        for &(c, z) in sq {
            let b = m.verts.len() as u32;
            m.verts.extend([[-c, -c, z], [c, -c, z], [c, c, z], [-c, c, z]]);
            m.tris.extend([[b, b + 1, b + 2], [b, b + 2, b + 3]]);
        }
        m
    }

    /// A plate raised over a floor, seen straight down: both face the camera
    /// and are lit alike (an open tray's rim and pockets), so only the step
    /// in depth tells them apart. Row 32 crosses the floor (x 3..17), the
    /// plate (18..45, on z = 0.5) and the floor again.
    #[test]
    fn steps_in_depth_are_outlined_and_section_faces_hatched() {
        let m = squares(&[(1.0, 0.0), (0.5, 0.5)]);
        let grey = [[0.6; 3]; 8];
        let top = Camera::named("top").unwrap();
        let px = |img: &[u8], x: usize| img[(32 * 64 + x) * 3];
        let plain = super::draw(&[Object::plain(&m, &grey)], &top, 64, [0, 0], None).0;
        assert_eq!(px(&plain, 10), px(&plain, 25), "floor and plate are lit alike");
        assert!(px(&plain, 18) < px(&plain, 25) / 2 && px(&plain, 3) < px(&plain, 10) / 2, "the plate's rim and the silhouette are lines");
        assert_eq!(px(&plain, 17), px(&plain, 10), "on the nearer side only");
        let cut = super::draw(&[Object::plain(&m, &grey)], &top, 64, [0, 0], Some((2, 0.5))).0;
        assert!(px(&cut, 24) < px(&plain, 24) && px(&cut, 25) == px(&plain, 25), "the plate on the plane is hatched");
        assert_eq!(px(&cut, 8), px(&plain, 8), "the floor off it is not");
    }

    /// Two bodies sharing a face, one's mesh off it by its own error; then a
    /// face as far proud of another, both meshes exact (a cover on a block):
    /// whatever the scene's size, it shows.
    #[test]
    fn coincident_faces_of_two_objects_show_one() {
        let (a, b) = (squares(&[(1.0, 0.0)]), squares(&[(1.0, 0.001)]));
        let (red, blue) = ([[1.0, 0.0, 0.0]; 4], [[0.0, 0.0, 1.0]; 4]);
        let top = Camera::named("top").unwrap();
        let img = super::draw(&[Object::plain(&a, &red), Object { error: &[0.001; 4], ..Object::plain(&b, &blue) }], &top, 64, [0, 0], None).0;
        assert!(img.chunks(3).all(|p| p[2] <= p[0] || p == [245; 3]), "the blue square never shows through the red one");
        let img = super::draw(&[Object::plain(&a, &red), Object::plain(&b, &blue)], &top, 64, [0, 0], None).0;
        assert!(img.chunks(3).all(|p| p[0] <= p[2]), "the blue square, nearer, hides the red one");
    }

    /// A sphere of 24 facets round, `n` across.
    fn sphere(n: usize) -> Mesh {
        let mut m = Mesh { verts: vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]], tris: vec![] };
        let (rings, at) = (n / 2, |r: usize, s: usize| (2 + (r - 1) * n + s % n) as u32);
        for r in 1..rings {
            let th = std::f64::consts::PI * r as f64 / rings as f64;
            m.verts.extend((0..n).map(|s| {
                let ph = std::f64::consts::TAU * s as f64 / n as f64;
                [th.sin() * ph.cos(), th.sin() * ph.sin(), th.cos()]
            }));
        }
        for s in 0..n {
            m.tris.push([0, at(1, s), at(1, s + 1)]);
            m.tris.push([1, at(rings - 1, s + 1), at(rings - 1, s)]);
            for r in 1..rings - 1 {
                m.tris.extend([[at(r, s), at(r + 1, s), at(r + 1, s + 1)], [at(r, s), at(r + 1, s + 1), at(r, s + 1)]]);
            }
        }
        m
    }

    /// A sphere's facets (15 degrees apart) are shaded as one smooth
    /// surface, by the normals given at its corners or by the mesh's own:
    /// no line inside its outline, and no step in shade between
    /// neighbouring pixels. Shaded facet by facet, neighbours differed by up
    /// to 30 levels.
    #[test]
    fn a_faceted_smooth_surface_is_shaded_smooth() {
        let m = sphere(24);
        let grey = vec![[0.6; 3]; m.verts.len()];
        let exact: Vec<[f32; 3]> = m.tris.iter().flat_map(|t| t.map(|i| m.verts[i as usize].map(|x| x as f32))).collect();
        for normals in [&[][..], &exact[..]] {
            let img = super::draw(&[Object { normals, ..Object::plain(&m, &grey) }], &Camera::named("iso").unwrap(), 128, [0, 0], None).0;
            // The outline's radius is 0.9 * 64 pixels.
            let inside = |x: usize, y: usize| ((x as f64 - 63.5).powi(2) + (y as f64 - 63.5).powi(2)).sqrt() < 0.8 * 57.6;
            let mut step = 0;
            for y in 0..127 {
                for x in (0..127).filter(|&x| inside(x, y) && inside(x + 1, y + 1)) {
                    let k = (y * 128 + x) * 3;
                    step = step.max(img[k].abs_diff(img[k + 3])).max(img[k].abs_diff(img[k + 128 * 3]));
                }
            }
            assert!(step <= 6, "neighbouring pixels differ by {step} levels, normals given: {}", !normals.is_empty());
        }
    }

    /// A ridge (two faces meeting square, seen from above) is a crease:
    /// a line along it, darker than the faces either side, on a dark
    /// surface as on a light one.
    #[test]
    fn a_crease_is_a_dark_line_on_any_surface() {
        let roof = Mesh {
            verts: vec![[-1.0, -1.0, 0.0], [0.0, -1.0, 1.0], [0.0, 1.0, 1.0], [-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0]],
            tris: vec![[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]],
        };
        for shade in [0.9, 0.1] {
            let paint = [[shade; 3]; 6];
            let img = super::draw(&[Object::plain(&roof, &paint)], &Camera::named("top").unwrap(), 64, [0, 0], None).0;
            let px = |x: usize| img[(32 * 64 + x) * 3];
            let ridge = (30..34).map(px).min().unwrap();
            assert!(ridge < px(20) / 2 && ridge < px(44) / 2, "the ridge {ridge}, the faces {} and {}", px(20), px(44));
        }
    }

    #[test]
    fn a_colour_scale_reads_its_range_under_the_image() {
        let (img, h) = with_scale(vec![245; 512 * 100 * 3], 512, "stress (MPa)", (0.0, 34.66), &|t| [t as f32, 0.0, 1.0 - t as f32], "");
        assert_eq!((img.len(), h), (512 * 116 * 3, 116));
        let row: Vec<&[u8]> = img[(108 * 512) * 3..(109 * 512) * 3].chunks(3).collect();
        let (blue, red) = (row.iter().position(|p| p == &[0, 0, 255]).unwrap(), row.iter().rposition(|p| p == &[255, 0, 0]).unwrap());
        assert!(blue < red && red - blue > 200, "the bar runs low to high, left to right");
        assert_eq!([significant(34.66), significant(0.0), significant(1250.4), significant(0.012345)], ["34.7", "0", "1250", "0.0123"]);
    }

    #[test]
    fn each_view_is_labelled() {
        let m = squares(&[(1.0, 0.0)]);
        let (img, w, _, views) = render(&[Object::plain(&m, &[[0.6; 3]; 4])], &Camera::sheet(), 512, None);
        for v in &views {
            let label = (v.origin[1] + 4..v.origin[1] + 18).any(|y| (v.origin[0] + 4..v.origin[0] + 60).any(|x| img[(y * w + x) * 3..][..3] == [70; 3]));
            assert!(label, "{}", v.name);
        }
        assert_eq!((w, views.len()), (512, 4));
        let (_, w, h, one) = render(&[Object::plain(&m, &[[0.6; 3]; 4])], &[Camera::named("back").unwrap()], 300, None);
        assert_eq!((w, h, one[0].right), (300, 300, [-1.0, 0.0, 0.0]));
    }
}
