//! File writers: binary STL, 3MF (multi-object, per-vertex color), PNG.

use ink_core::v3::{norm, normal};
use ink_kernel::mesh::Mesh;
use std::io::Write;

pub struct Object<'a> {
    pub name: &'a str,
    pub mesh: &'a Mesh,
    /// One RGB color (0..1) per vertex.
    pub colors: &'a [[f32; 3]],
}

/// Binary STL; `note` goes in the 80-byte header after "ink", cut to fit
/// (a header must not start with "solid", which marks ASCII STL).
pub fn stl(objects: &[Object], note: &str) -> Vec<u8> {
    let n: usize = objects.iter().map(|o| o.mesh.tris.len()).sum();
    let mut out = Vec::with_capacity(84 + 50 * n);
    let mut header = [b' '; 80];
    let tag = format!("ink {note}");
    let tag = &tag.as_bytes()[..tag.len().min(80)];
    header[..tag.len()].copy_from_slice(tag);
    out.extend_from_slice(&header);
    out.extend_from_slice(&(n as u32).to_le_bytes());
    for o in objects {
        for t in &o.mesh.tris {
            let v = t.map(|i| o.mesh.verts[i as usize]);
            let mut nrm = normal(v);
            let l = norm(nrm);
            if l > 0.0 {
                nrm = nrm.map(|c| c / l);
            }
            for c in nrm.iter().chain(v.iter().flatten()) {
                out.extend_from_slice(&(*c as f32).to_le_bytes());
            }
            out.extend_from_slice(&[0, 0]);
        }
    }
    out
}

fn hex(c: [f32; 3]) -> String {
    let b = c.map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8);
    format!("#{:02X}{:02X}{:02X}", b[0], b[1], b[2])
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// A placed copy of a mesh object: 3MF build items reference one mesh
/// resource with a transform (rotation rows r, translation t), so thousands of
/// instances cost one mesh plus a line each.
pub struct Placed<'a> {
    pub object: usize,
    pub name: &'a str,
    pub transform: Option<([[f64; 3]; 3], [f64; 3])>,
}

/// 3MF with one mesh resource per object and a build item per placement;
/// `description` goes in the model's Description metadata.
pub fn three_mf_instanced(objects: &[Object], items: &[Placed], description: &str) -> Vec<u8> {
    use rayon::prelude::*;
    use std::fmt::Write as _;
    const CHUNK: usize = 1 << 16;
    let mut m = String::new();
    m.push_str(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<model unit=\"millimeter\" xml:lang=\"en-US\" \
         xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\" \
         xmlns:m=\"http://schemas.microsoft.com/3dmanufacturing/material/2015/02\">\n",
    );
    m.push_str("<metadata name=\"Application\">ink</metadata>\n");
    if !description.is_empty() {
        let _ = writeln!(m, "<metadata name=\"Description\">{}</metadata>", xml_escape(description));
    }
    m.push_str("<resources>\n");
    for (k, o) in objects.iter().enumerate() {
        let (cg, obj) = (2 * k + 1, 2 * k + 2);
        let _ = writeln!(m, "<m:colorgroup id=\"{cg}\">");
        let colors: Vec<String> = o.colors.par_chunks(CHUNK).map(|c| c.iter().map(|c| format!("<m:color color=\"{}\"/>\n", hex(*c))).collect()).collect();
        colors.iter().for_each(|c| m.push_str(c));
        m.push_str("</m:colorgroup>\n");
        let _ = writeln!(m, "<object id=\"{obj}\" type=\"model\" name=\"{}\" pid=\"{cg}\" pindex=\"0\"><mesh><vertices>", xml_escape(o.name));
        let verts: Vec<String> = o
            .mesh
            .verts
            .par_chunks(CHUNK)
            .map(|c| {
                let mut s = String::with_capacity(c.len() * 48);
                for v in c {
                    let _ = writeln!(s, "<vertex x=\"{:.5}\" y=\"{:.5}\" z=\"{:.5}\"/>", v[0], v[1], v[2]);
                }
                s
            })
            .collect();
        verts.iter().for_each(|c| m.push_str(c));
        m.push_str("</vertices><triangles>\n");
        let tris: Vec<String> = o
            .mesh
            .tris
            .par_chunks(CHUNK)
            .map(|c| {
                let mut s = String::with_capacity(c.len() * 80);
                for t in c {
                    let _ = writeln!(s, "<triangle v1=\"{}\" v2=\"{}\" v3=\"{}\" pid=\"{cg}\" p1=\"{}\" p2=\"{}\" p3=\"{}\"/>", t[0], t[1], t[2], t[0], t[1], t[2]);
                }
                s
            })
            .collect();
        tris.iter().for_each(|c| m.push_str(c));
        m.push_str("</triangles></mesh></object>\n");
    }
    m.push_str("</resources>\n<build>\n");
    for it in items {
        let _ = write!(m, "<item objectid=\"{}\" partnumber=\"{}\"", 2 * it.object + 2, xml_escape(it.name));
        if let Some((r, t)) = it.transform {
            // 3MF multiplies row vectors: [x y z 1] M, so M's rows are R's columns.
            let _ =
                write!(m, " transform=\"{} {} {} {} {} {} {} {} {} {} {} {}\"", r[0][0], r[1][0], r[2][0], r[0][1], r[1][1], r[2][1], r[0][2], r[1][2], r[2][2], t[0], t[1], t[2]);
        }
        m.push_str("/>\n");
    }
    m.push_str("</build>\n</model>\n");

    let types = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\
<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\
<Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/></Types>";
    let rels = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\
<Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/></Relationships>";
    zip(&[("[Content_Types].xml", types.as_bytes()), ("_rels/.rels", rels.as_bytes()), ("3D/3dmodel.model", m.as_bytes())])
}

// ---------------------------------------------------------------- containers

fn crc32(data: &[u8]) -> u32 {
    let mut h = crc32fast::Hasher::new();
    h.update(data);
    h.finalize()
}

/// Raw deflate on all cores: independent 4 MiB chunks joined at sync-flush
/// boundaries (each chunk ends byte-aligned; only the last is final).
fn deflate_parallel(data: &[u8], level: u8) -> Vec<u8> {
    use miniz_oxide::deflate::core::{CompressorOxide, TDEFLFlush, TDEFLStatus, compress_to_output, create_comp_flags_from_zip_params};
    use rayon::prelude::*;
    if data.is_empty() {
        return miniz_oxide::deflate::compress_to_vec(data, level);
    }
    let chunks: Vec<&[u8]> = data.chunks(4 << 20).collect();
    let n = chunks.len();
    let parts: Vec<Vec<u8>> = chunks
        .par_iter()
        .enumerate()
        .map(|(i, c)| {
            let mut d = CompressorOxide::new(create_comp_flags_from_zip_params(level as i32, -15, 0));
            let flush = if i + 1 == n { TDEFLFlush::Finish } else { TDEFLFlush::Sync };
            let mut out = Vec::with_capacity(c.len() / 4);
            let mut input: &[u8] = c;
            loop {
                let (st, used) = compress_to_output(&mut d, input, flush, |b| {
                    out.extend_from_slice(b);
                    true
                });
                input = &input[used..];
                match st {
                    TDEFLStatus::Done => break,
                    TDEFLStatus::Okay if input.is_empty() && flush == TDEFLFlush::Sync => break,
                    TDEFLStatus::Okay => {}
                    other => panic!("deflate failed: {other:?}"),
                }
            }
            out
        })
        .collect();
    parts.concat()
}

/// Deflate-compressed zip archive.
fn zip(files: &[(&str, &[u8])]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut central = Vec::new();
    for (name, data) in files {
        let comp = deflate_parallel(data, 3);
        let crc = crc32(data);
        let offset = out.len() as u32;
        let header = |sig: u32, central: bool| {
            let mut h = Vec::new();
            h.extend_from_slice(&sig.to_le_bytes());
            if central {
                h.extend_from_slice(&20u16.to_le_bytes());
            }
            h.extend_from_slice(&20u16.to_le_bytes()); // version needed
            h.extend_from_slice(&0u16.to_le_bytes()); // flags
            h.extend_from_slice(&8u16.to_le_bytes()); // deflate
            h.extend_from_slice(&0u16.to_le_bytes()); // time
            h.extend_from_slice(&0x21u16.to_le_bytes()); // date 1980-01-01
            h.extend_from_slice(&crc.to_le_bytes());
            h.extend_from_slice(&(comp.len() as u32).to_le_bytes());
            h.extend_from_slice(&(data.len() as u32).to_le_bytes());
            h.extend_from_slice(&(name.len() as u16).to_le_bytes());
            h.extend_from_slice(&0u16.to_le_bytes()); // extra
            if central {
                h.extend_from_slice(&[0; 6]); // comment len, disk start, internal attrs
                h.extend_from_slice(&0u32.to_le_bytes()); // external attrs
                h.extend_from_slice(&offset.to_le_bytes());
            }
            h.extend_from_slice(name.as_bytes());
            h
        };
        out.extend(header(0x0403_4b50, false));
        out.extend_from_slice(&comp);
        central.extend(header(0x0201_4b50, true));
    }
    let cd_offset = out.len() as u32;
    out.extend_from_slice(&central);
    out.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    out.extend_from_slice(&[0; 4]);
    out.extend_from_slice(&(files.len() as u16).to_le_bytes());
    out.extend_from_slice(&(files.len() as u16).to_le_bytes());
    out.extend_from_slice(&(central.len() as u32).to_le_bytes());
    out.extend_from_slice(&cd_offset.to_le_bytes());
    out.extend_from_slice(&0u16.to_le_bytes());
    out
}

pub fn png(w: usize, h: usize, rgb: &[u8]) -> Vec<u8> {
    let mut raw = Vec::with_capacity((w * 3 + 1) * h);
    for row in rgb.chunks(w * 3) {
        raw.push(0);
        raw.extend_from_slice(row);
    }
    let idat = miniz_oxide::deflate::compress_to_vec_zlib(&raw, 6);
    let mut out = vec![0x89, b'P', b'N', b'G', 0x0D, 0x0A, 0x1A, 0x0A];
    let mut chunk = |kind: &[u8], data: &[u8]| {
        out.extend_from_slice(&(data.len() as u32).to_be_bytes());
        let mut c = kind.to_vec();
        c.extend_from_slice(data);
        out.extend_from_slice(&c);
        out.extend_from_slice(&crc32(&c).to_be_bytes());
    };
    let mut ihdr = Vec::new();
    ihdr.extend_from_slice(&(w as u32).to_be_bytes());
    ihdr.extend_from_slice(&(h as u32).to_be_bytes());
    ihdr.extend_from_slice(&[8, 2, 0, 0, 0]);
    chunk(b"IHDR", &ihdr);
    chunk(b"IDAT", &idat);
    chunk(b"IEND", &[]);
    out
}

/// One outline for DXF: its layer (the target, or a mark's name), points
/// in the plane's axes, whether it closes, and whether it is a mark (a
/// line drawn, not cut: a sheet's bend).
pub struct Outline<'a> {
    pub layer: &'a str,
    pub pts: &'a [[f64; 2]],
    pub closed: bool,
    pub mark: bool,
}

/// DXF (R12, the version every cutting service and CAD reader takes), in
/// mm: one POLYLINE per outline on its layer, an open one of two points a
/// LINE; `note` is a comment at the top saying how it was made. A layer
/// of cuts is continuous, white; one of marks dashed, red: shops read a
/// bend as a dashed line at its centre on a layer of its own
/// (SendCutSend, lasercut.parts, QuickCut's BEND). Err: a cut and a mark
/// whose names are one layer as DXF reads them, naming both and it.
pub fn dxf(outlines: &[Outline], note: &str) -> Result<Vec<u8>, String> {
    // A layer name may not hold <>/\":;?*|=' (DXF reference), and is read
    // in any case: one layer per name as read, spelled as first written.
    let clean = |s: &str| s.chars().map(|c| if "<>/\\\":;?*|='".contains(c) || c.is_whitespace() { '_' } else { c }).collect::<String>();
    let mut layers: Vec<(String, bool, &str)> = Vec::new();
    for o in outlines {
        let l = clean(o.layer);
        match layers.iter().find(|(m, _, _)| m.eq_ignore_ascii_case(&l)) {
            None => layers.push((l, o.mark, o.layer)),
            Some((m, mark, first)) if *mark != o.mark => {
                let (cut, mark) = if o.mark { (first, &o.layer) } else { (&o.layer, first) };
                return Err(format!("the mark '{mark}' and the target '{cut}' would share DXF layer {}", m.to_uppercase()));
            }
            Some(_) => {}
        }
    }
    let layer = |s: &str| {
        let l = clean(s);
        layers.iter().find(|(m, _, _)| m.eq_ignore_ascii_case(&l)).map_or(l, |(m, _, _)| m.clone())
    };
    let mut s = format!("999\n{}\n0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1009\n9\n$INSUNITS\n70\n4\n0\nENDSEC\n", note.replace('\n', " "));
    // ISO's DASHED (acadiso.lin): 12.7 on, 6.35 off.
    s += "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLTYPE\n70\n2\n\
          0\nLTYPE\n2\nCONTINUOUS\n70\n0\n3\nSolid line\n72\n65\n73\n0\n40\n0.0\n\
          0\nLTYPE\n2\nDASHED\n70\n0\n3\nDashed __ __ __\n72\n65\n73\n2\n40\n19.05\n49\n12.7\n49\n-6.35\n0\nENDTAB\n";
    s += &format!("0\nTABLE\n2\nLAYER\n70\n{}\n", layers.len());
    for (l, mark, _) in &layers {
        let (color, ltype) = if *mark { (1, "DASHED") } else { (7, "CONTINUOUS") };
        s += &format!("0\nLAYER\n2\n{l}\n70\n0\n62\n{color}\n6\n{ltype}\n");
    }
    s += "0\nENDTAB\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n";
    for o in outlines {
        let l = layer(o.layer);
        if let ([a, b], false) = (o.pts, o.closed) {
            s += &format!("0\nLINE\n8\n{l}\n10\n{}\n20\n{}\n30\n0\n11\n{}\n21\n{}\n31\n0\n", a[0], a[1], b[0], b[1]);
            continue;
        }
        // 128: the linetype's pattern runs on across vertices (a mark's
        // dashes, else short chords of a curve draw solid).
        s += &format!("0\nPOLYLINE\n8\n{l}\n66\n1\n70\n{}\n10\n0\n20\n0\n30\n0\n", u8::from(o.closed) | if o.mark { 128 } else { 0 });
        for p in o.pts {
            s += &format!("0\nVERTEX\n8\n{l}\n10\n{}\n20\n{}\n30\n0\n", p[0], p[1]);
        }
        s += &format!("0\nSEQEND\n8\n{l}\n");
    }
    s += "0\nENDSEC\n0\nEOF\n";
    Ok(s.into_bytes())
}

pub fn write(path: &str, bytes: &[u8]) -> std::io::Result<()> {
    if let Some(dir) = std::path::Path::new(path).parent()
        && !dir.as_os_str().is_empty()
    {
        std::fs::create_dir_all(dir)?;
    }
    std::fs::File::create(path)?.write_all(bytes)
}
