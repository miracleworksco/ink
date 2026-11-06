//! How a body is drawn, read once when it is loaded: its look (`color`,
//! `metallic`, `roughness`, `transmission`). A malformed one is refused at
//! its line; clients (the CLI's meshes and GLB, viewers, packages) read
//! `Look`, never the attributes' layout.

use crate::Value;
use crate::graph::{Graph, NodeId};

/// A body's look: each channel a field (a number is a constant) meant in
/// 0..1 (clients clamp a field's values), or None where the body does not
/// say: then `PALETTE` (by the body's index) and the defaults below.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Look {
    pub color: Option<[NodeId; 3]>,
    pub metallic: Option<NodeId>,
    pub roughness: Option<NodeId>,
    pub transmission: Option<NodeId>,
}

/// The colours of bodies that give none, the k-th body (or target drawn)
/// taking the k-th (mod 12): twelve told apart at a glance, the first six
/// as they were.
pub const PALETTE: [[f32; 3]; 12] = [
    [0.30, 0.55, 0.85],
    [0.90, 0.55, 0.20],
    [0.35, 0.70, 0.40],
    [0.80, 0.30, 0.35],
    [0.60, 0.45, 0.80],
    [0.55, 0.55, 0.55],
    [0.20, 0.70, 0.70],
    [0.90, 0.80, 0.25],
    [0.85, 0.45, 0.70],
    [0.55, 0.40, 0.25],
    [0.65, 0.80, 0.30],
    [0.20, 0.30, 0.60],
];
/// Metallic, roughness and transmission where a body gives none: a
/// dielectric, satin, opaque.
pub const METALLIC: f32 = 0.0;
pub const ROUGHNESS: f32 = 0.4;
pub const TRANSMISSION: f32 = 0.0;

impl Look {
    /// The look among a body's (or an item's, placed) attributes, as
    /// `read` left them at load.
    pub fn of(attrs: &[(String, Value)]) -> Look {
        let get = |k: &str| attrs.iter().find(|(n, _)| n == k).map(|(_, v)| v);
        let field = |k: &str| if let Some(Value::Field(f)) = get(k) { Some(*f) } else { None };
        let color = match get("color") {
            Some(Value::Vec(c)) => match c.as_slice() {
                [Value::Field(r), Value::Field(g), Value::Field(b)] => Some([*r, *g, *b]),
                _ => None,
            },
            _ => None,
        };
        Look { color, metallic: field("metallic"), roughness: field("roughness"), transmission: field("transmission") }
    }

    /// Its channels in order (red, green, blue, metallic, roughness,
    /// transmission).
    pub fn channels(&self) -> [Option<NodeId>; 6] {
        let c = self.color.map_or([None; 3], |c| c.map(Some));
        [c[0], c[1], c[2], self.metallic, self.roughness, self.transmission]
    }
}

/// An attribute as the language keeps it: a look's `color` as (vec r g b)
/// of fields, "#RRGGBB" read into numbers; `metallic`, `roughness` and
/// `transmission` a field; a number outside 0..1 refused. Other
/// attributes are kept as written.
pub(crate) fn read(g: &mut Graph, key: &str, v: Value) -> Result<Value, String> {
    let unit = |g: &Graph, f: NodeId, what: &str| match g.as_const(f) {
        Some(c) if !(0.0..=1.0).contains(&c) => Err(format!("{what} {c} is outside 0..1")),
        _ => Ok(()),
    };
    match key {
        "color" => {
            let rgb = match &v {
                Value::Str(s) => hex(s).map(|c| c.map(|c| g.constant(c))),
                Value::Vec(c) => match c.as_slice() {
                    [Value::Field(r), Value::Field(gr), Value::Field(b)] => Some([*r, *gr, *b]),
                    _ => None,
                },
                _ => None,
            };
            let Some(rgb) = rgb else {
                return Err(format!("color {} is not a colour: give \"#RRGGBB\" (six hex digits) or (vec r g b) of numbers or fields in 0..1", shown(&v)));
            };
            for (c, name) in rgb.iter().zip(["red", "green", "blue"]) {
                unit(g, *c, &format!("color's {name}"))?;
            }
            Ok(Value::Vec(rgb.iter().map(|&c| Value::Field(c)).collect()))
        }
        "metallic" | "roughness" | "transmission" => match v {
            Value::Field(f) => unit(g, f, key).map(|_| v),
            _ => Err(format!("{key} {} is not a number or a scalar field in 0..1", shown(&v))),
        },
        _ => Ok(v),
    }
}

fn hex(s: &str) -> Option<[f64; 3]> {
    let h = s.strip_prefix('#').filter(|h| h.len() == 6 && h.is_ascii())?;
    let c = |i: usize| u8::from_str_radix(&h[i..i + 2], 16).ok().map(|b| b as f64 / 255.0);
    Some([c(0)?, c(2)?, c(4)?])
}

/// A value as an error shows it: strings quoted, a vec by its elements'
/// kinds.
fn shown(v: &Value) -> String {
    match v {
        Value::Str(s) => format!("{s:?}"),
        Value::Vec(items) => format!("(vec {})", items.iter().map(|i| i.kind()).collect::<Vec<_>>().join(" ")),
        other => format!("(a {})", other.kind()),
    }
}
