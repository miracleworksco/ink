//! Reading a study's loads, and the ranges its material's numbers may
//! take.
use super::*;

/// (load name (force fx fy fz) | (torque tx ty tz) | (pressure MPa) (region r)), or (source name (power w) (region r)).
pub(super) fn read_load(r: &mut Reader, c: &Sexp, head: &str, lname: &Sexp, rest: &[Sexp], usage: &str) -> R<StudyLoad> {
    let (mut vector, mut region, mut vector_t, mut pressure, mut torque) = (None, None, None, None, None);
    let mut from_model: Option<[NodeId; 3]> = None;
    // Components that read what the model works out are fields until the study is asked.
    let modelled = |d: &mut Reader, es: &[&Sexp], q: Quantity, what: &str| -> R<Option<[NodeId; 3]>> {
        let mut nodes = [0; 3];
        for (i, e) in es.iter().enumerate() {
            nodes[i] = d.field_of(e, q, what)?;
        }
        if !nodes.iter().any(|&n| d.reads_solved(n)) {
            return Ok(None);
        }
        if let Some(i) = nodes.iter().position(|&n| d.graph().depends_on(n, |v| matches!(v, Var::X | Var::Y | Var::Z | Var::T))) {
            return Err(
                d.err(es[i], "a load that reads what the model works out (a reaction, a tension) may use it and params, not x y z or t: it is taken at the motion's worst moment")
            );
        }
        nodes.iter().for_each(|&n| d.at_worst(n));
        Ok(Some(nodes))
    };
    // Components: numbers known at load time, or expressions of t.
    let comps = |d: &mut Reader, vector: &mut Option<[f64; 3]>, vector_t: &mut Option<[NodeId; 3]>, es: &[&Sexp], q: Quantity, what: &str| -> R<()> {
        let mut nodes = [0; 3];
        for (i, e) in es.iter().enumerate() {
            nodes[i] = d.field_of(e, q, what)?;
            if d.graph().depends_on(nodes[i], |v| matches!(v, Var::X | Var::Y | Var::Z)) {
                return Err(d.err(e, "load components may use t and params, not x y z (use the region to place it)"));
            }
        }
        if es.len() == 1 {
            nodes[1] = d.graph().constant(0.0);
            nodes[2] = nodes[1];
        }
        if nodes.iter().any(|&n| d.graph().depends_on(n, |v| v == Var::T) || d.graph().motion_reads(n).iter().any(|r| matches!(r, crate::graph::MotionRead::Drive(_)))) {
            // Its t is the transient's own time, not the model's.
            for &n in &nodes[..es.len()] {
                d.on_own_clock(n);
            }
            *vector = Some([0.0; 3]);
            *vector_t = Some(nodes);
        } else {
            let mut v = [0.0; 3];
            for (i, n) in nodes.iter().enumerate() {
                v[i] = d.load_value(*n).ok_or_else(|| d.err(es[i.min(es.len() - 1)], "expected a number"))?;
            }
            *vector = Some(v);
        }
        Ok(())
    };
    for item in rest {
        let l = item.list().unwrap_or(&[]);
        match (l.first().and_then(|h| h.sym()), &l.get(1..).unwrap_or(&[])) {
            (Some("force"), [a, b, cc]) if head == "load" => match modelled(r, &[a, b, cc], Quantity::Force, "force")? {
                Some(nodes) => (from_model, vector) = (Some(nodes), Some([0.0; 3])),
                None => comps(r, &mut vector, &mut vector_t, &[a, b, cc], Quantity::Force, "force")?,
            },
            (Some("power"), [w]) if head == "source" => comps(r, &mut vector, &mut vector_t, &[w], Quantity::Power, "power")?,
            (Some("torque"), [a, b, cc]) if head == "load" => {
                match modelled(r, &[a, b, cc], Quantity::Torque, "torque")? {
                    Some(nodes) => (from_model, torque) = (Some(nodes), Some([0.0; 3])),
                    None => {
                        let mut n = |e| r.number_of(e, Quantity::Torque, "torque");
                        torque = Some([n(a)?, n(b)?, n(cc)?]);
                    }
                }
                vector = Some([0.0; 3]);
            }
            (Some("pressure"), [pr]) if head == "load" => {
                pressure = Some(r.field_of(pr, Quantity::Pressure, "pressure")?);
                vector = Some([0.0; 3]);
            }
            (Some("region"), [e]) => region = Some(Region::read(r, e)?),
            _ => return Err(r.err_hint(item, format!("bad {head} clause"), Some(usage.into()))),
        }
    }
    let what = if head == "load" { "(force fx fy fz), (torque tx ty tz) or (pressure MPa)" } else { "(power watts)" };
    Ok(StudyLoad {
        from_model,
        vector_t,
        pressure,
        torque,
        name: r.name(lname)?,
        vector: vector.ok_or_else(|| r.err(c, format!("{head} needs {what}")))?,
        region: region.ok_or_else(|| r.err(c, format!("{head} needs (region field)")))?,
    })
}

/// Young's moduli of real materials, MPa: soft gels' (about 1e-4) to
/// diamond's (about 1.1e6).
pub const MODULUS_RANGE: (f64, f64) = (1e-4, 1.2e6);

pub(super) const MODULUS: Magnitude = Magnitude { name: "E (MPa)", range: MODULUS_RANGE, strict: false, zero: false, why: modulus_message };
pub(super) const YIELD: Magnitude = Magnitude {
    name: "yield strength (MPa)",
    range: (1e-3, 1e5),
    strict: false,
    zero: false,
    why: |v| format!("yield {v} MPa is outside the supported strength range; write its stress unit, e.g. 250MPa or 250e6Pa"),
};
pub(super) const POISSON: Magnitude = Magnitude { name: "nu", range: (-1.0, 0.5), strict: true, zero: false, why: |_| "nu must be in (-1, 0.5)".into() };

fn modulus_message(e: f64) -> String {
    let (_, hi) = MODULUS_RANGE;
    let why = if e <= 0.0 {
        "is not positive".to_string()
    } else if e > hi {
        format!("is stiffer than any material in MPa, ink's modulus unit (diamond's is about 1.1e6): if it is in Pa, write {e:e}Pa; if in psi, {e:e}psi")
    } else {
        "is softer than any solid in MPa, ink's modulus unit (a soft gel's is about 1e-3): converted twice? write it with its unit, e.g. 200GPa".to_string()
    };
    format!("E {e:e} {why}")
}
