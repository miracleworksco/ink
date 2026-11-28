//! Reading `(study ...)`.
use super::*;

/// (study name (shape f) (elastic (E v) (nu v) [(density v)] [(yield v)]) | (heat (k v))
///   (fix name region [(dofs x y z)] [(pivot ...)] [(slide ...)]) | (temperature name region value)
///   (load name (force fx fy fz) (region r)) | (source name (power w) (region r))
///   [(gravity gx gy gz)] [(res h)])
pub(crate) fn read_study(r: &mut Reader, form: &Sexp, args: &[Sexp]) -> R<Study> {
    let usage = "usage: (study name (shape f) (elastic (E v) (nu v) [(density v)] [(yield v)]) \
                 [(mirror name axis coordinate)] (fix name region [(support)] [(dofs x y z)] [(move dx dy dz)] [(pivot (at x y z) [(axis dx dy dz)])] [(slide dx dy dz)]) (load name (force fx fy fz) | (pressure MPa) (region r)) [(gravity gx gy gz)] [(res h)] [(support-margin mm)] \
); \
                 or with (heat (k v)), (temperature name region value), (source name (power w) (region r))";
    let Some((name, clauses)) = args.split_first() else { return Err(r.err_hint(form, "bad study", Some(usage.into()))) };
    let n = r.name(name)?;
    r.new_result_name(name, &n)?;
    let (mut shape, mut physics, mut gravity, mut res, mut margin, mut modes) = (None, None, None, None, None, None);
    let mut property_sites = std::collections::HashMap::new();
    // (each material property and where it came from: as written, or the material named)
    let mut material_from: Vec<(String, String)> = Vec::new();
    let mut transient = None;
    // (the first held temperature that reads t: a transient's alone)
    let mut timed_held: Option<Sexp> = None;
    let (mut fixes, mut loads) = (Vec::new(), Vec::new());
    let mut any_move = false;
    // (key value) pairs inside a clause, e.g. (E 2100) (nu 0.35).
    let pairs = |d: &mut Reader, items: &[Sexp]| -> R<Vec<(String, Sexp)>> {
        items
            .iter()
            .map(|it| match it.list() {
                Some([k, v]) => Ok((d.name(k)?, v.clone())),
                _ => Err(d.err_hint(it, "expected (key value)", Some(usage.into()))),
            })
            .collect()
    };
    for c in &r.clauses(clauses)? {
        let (head, args) = c.head_args().ok_or_else(|| r.err_hint(c, "study clauses are lists", Some(usage.into())))?;
        match (head, args) {
            ("shape", [e]) => shape = Some(r.field(e)?),
            ("elastic", kv) => {
                let (mut e, mut nu, mut density, mut ys) = (None, None, None, None);
                let kv = pairs(r, kv)?;
                if let Some((_, m)) = kv.iter().find(|(k, _)| k == "material") {
                    let m = r.material(m)?;
                    let num = |k: &str| m.get(k).and_then(|v| v.as_f64());
                    e = num("E").map(|v| r.graph().constant(v));
                    nu = num("nu").map(|v| r.graph().constant(v));
                    density = num("density").map(|v| r.graph().constant(v));
                    ys = num("yield").map(|v| r.graph().constant(v));
                    material_from = from_material(&m, ["E", "nu", "density", "yield"], &kv);
                } else {
                    material_from.clear();
                }
                property_sites.clear();
                for (k, v) in kv {
                    property_sites.insert(k.clone(), v.clone());
                    said(&mut material_from, &k, &v);
                    match k.as_str() {
                        "material" => {}
                        "E" => e = Some(r.material_field(&v, Quantity::Pressure, "E")?),
                        "nu" => match r.value(&v)? {
                            Value::Field(f, None) => nu = Some(f),
                            _ => return Err(r.err(&v, "nu must be a dimensionless number or field")),
                        },
                        "density" => density = Some(r.material_field(&v, Quantity::Density, "density")?),
                        "yield" => ys = Some(r.material_field(&v, Quantity::Pressure, "yield")?),
                        other => return Err(r.err(c, format!("unknown elastic property '{other}' (E nu density yield material)"))),
                    }
                }
                let e = e.ok_or_else(|| r.err(c, "elastic needs (E modulus-in-MPa)"))?;
                // No value stands in for it: a part held sideways answers by it (a strip clamped along its
                // edges is 1 / (1 - nu^2) stiffer), and a claim would rest on a number its file never gave.
                let nu =
                    nu.ok_or_else(|| r.err(c, "elastic needs (nu ratio), Poisson's ratio (most metals 0.27 to 0.34, plastics 0.35 to 0.42), or a (material name) that has it"))?;
                physics = Some(StudyPhysics::Elastic { e, nu, density, yield_strength: ys });
            }
            ("heat", kv) => {
                let (mut k, mut density, mut capacity) = (None, None, None);
                let kv = pairs(r, kv)?;
                if let Some((_, m)) = kv.iter().find(|(k, _)| k == "material") {
                    let m = r.material(m)?;
                    let num = |key: &str| m.get(key).and_then(|v| v.as_f64());
                    k = num("k").map(|v| r.graph().constant(v));
                    density = num("density").map(|v| r.graph().constant(v));
                    capacity = num("capacity");
                    material_from = from_material(&m, ["k", "density", "capacity"], &kv);
                    if !kv.iter().any(|(key, _)| key == "capacity")
                        && let Some(why) = capacity.and_then(|v| materials::CAPACITY.error(v))
                    {
                        return Err(r.err(c, why));
                    }
                } else {
                    material_from.clear();
                }
                property_sites.clear();
                for (key, v) in kv {
                    property_sites.insert(key.clone(), v.clone());
                    said(&mut material_from, &key, &v);
                    match key.as_str() {
                        "material" => {}
                        "k" => k = Some(r.material_field(&v, Quantity::Conductivity, "k")?),
                        "density" => density = Some(r.material_field(&v, Quantity::Density, "density")?),
                        "capacity" => capacity = Some(r.bounded_number_of(&v, Quantity::Capacity, materials::CAPACITY)?),
                        other => return Err(r.err(c, format!("unknown heat property '{other}' (k density capacity material)"))),
                    }
                }
                physics = Some(StudyPhysics::Heat { k: k.ok_or_else(|| r.err(c, "heat needs (k W-per-m-K)"))?, density, capacity });
            }
            ("transient", kv) => {
                let mut tr = crate::Transient { duration: 0.0, steps: 100, initial: 0.0, capacity: 0.0 };
                for (key, v) in pairs(r, kv)? {
                    let q = match key.as_str() {
                        "duration" => Some(Quantity::Time),
                        "initial" => Some(Quantity::Temperature),
                        _ => None,
                    };
                    let x = match q {
                        Some(q) => r.number_of(&v, q, &key)?,
                        None => r.number(&v)?,
                    };
                    match key.as_str() {
                        "duration" => tr.duration = x,
                        "steps" => tr.steps = x.round().max(1.0) as usize,
                        "initial" => tr.initial = x,
                        other => return Err(r.err(c, format!("unknown transient option '{other}' (duration steps initial)"))),
                    }
                }
                if !(tr.duration > 0.0) {
                    return Err(r.err(c, "transient needs (duration seconds) > 0"));
                }
                if tr.steps > 100_000 {
                    return Err(r.err(c, "at most 100000 steps"));
                }
                transient = Some(tr);
            }
            ("mirror", [fname, axis, coord]) => {
                let axis = match axis.sym() {
                    Some("x") => 0,
                    Some("y") => 1,
                    Some("z") => 2,
                    _ => return Err(r.err(axis, "a mirror axis must be x, y or z")),
                };
                let at = r.number_of(coord, Quantity::Length, "mirror coordinate")?;
                if fixes.iter().any(|f: &StudyFix| matches!(f.boundary, crate::Boundary::Mirror { axis: a, .. } if a == axis)) {
                    return Err(r.err(c, "parallel mirrors describe a repeating part; model its finite length and declare at most one mirror per axis"));
                }
                let zero = r.graph().constant(0.0);
                fixes.push(StudyFix {
                    name: r.name(fname)?,
                    region: Region::Field(zero),
                    dofs: std::array::from_fn(|d| d == axis),
                    value: [zero; 3],
                    film: None,
                    normal: false,
                    free: Vec::new(),
                    rigid: false,
                    boundary: crate::Boundary::Mirror { axis, at },
                });
            }
            ("fix", [fname, region, rest @ ..]) => {
                let mut boundary = crate::Boundary::Unspecified;
                let mut dofs = [true; 3];
                let mut normal = false;
                let mut free = Vec::new();
                let mut rigid = false;
                let mut seen_dofs = false;
                let zero = r.graph().constant(0.0);
                let mut moved: Option<([NodeId; 3], &Sexp)> = None;
                let fix_usage = "usage: (fix name region [(dofs x y z) | (dofs normal)] [(support)] [(move dx dy dz)] [(pivot (at x y z) [(axis dx dy dz)])] [(slide dx dy dz)] [(rigid)]). \
                                 A fix holds its region at zero displacement, or moved by (move ...) (mm, fields of x y z), along what \
                                 it holds, except the rigid motions of its region that pivot (turning about the axis, or about the point) \
                                 and slide leave free. With pivot or slide it holds its region's mean motion, the region free to deform; \
                                 (rigid) holds the region rigid but for its free motions";
                // Three numbers, a direction or a point (lengths).
                let vec3 = |r: &mut Reader, xs: &[Sexp], at: &Sexp, what: &str| -> R<V3> {
                    let n = |r: &mut Reader, e: &Sexp| if what == "at" { r.number_of(e, Quantity::Length, what) } else { r.number(e) };
                    match xs {
                        [a, b, c] => Ok([n(r, a)?, n(r, b)?, n(r, c)?]),
                        _ => Err(r.err_hint(at, format!("{what} takes three numbers"), Some(fix_usage.into()))),
                    }
                };
                let unit =
                    |r: &Reader, v: V3, at: &Sexp| -> R<V3> { ink_core::v3::direction(v).ok_or_else(|| r.err_hint(at, "a direction must be nonzero", Some(fix_usage.into()))) };
                for d in rest {
                    let head = d.list().and_then(|l| l.first()).and_then(|h| h.sym());
                    let dl = match (d.list(), head) {
                        (Some([_]), Some("support")) => {
                            boundary = crate::Boundary::Support;
                            continue;
                        }
                        (Some([_]), Some("rigid")) => {
                            rigid = true;
                            continue;
                        }
                        (Some(l), Some("pivot")) => {
                            let (mut at, mut axis) = (None, None);
                            for o in &l[1..] {
                                match o.list() {
                                    Some([k, xs @ ..]) if k.sym() == Some("at") && at.is_none() => at = Some(vec3(r, xs, o, "at")?),
                                    Some([k, xs @ ..]) if k.sym() == Some("axis") && axis.is_none() => {
                                        let v = vec3(r, xs, o, "axis")?;
                                        axis = Some(unit(r, v, o)?);
                                    }
                                    _ => return Err(r.err_hint(o, "pivot takes (at x y z) and, for a hinge, (axis dx dy dz)", Some(fix_usage.into()))),
                                }
                            }
                            let at = at.ok_or_else(|| r.err_hint(d, "pivot needs (at x y z)", Some(fix_usage.into())))?;
                            // No axis: a ball joint, free to turn about every axis through the point.
                            match axis {
                                Some(a) => free.push(crate::Motion::turn(at, a)),
                                None => free.extend((0..3).map(|k| crate::Motion::turn(at, std::array::from_fn(|j| (j == k) as u8 as f64)))),
                            }
                            continue;
                        }
                        (Some(l), Some("slide")) => {
                            let v = vec3(r, &l[1..], d, "slide")?;
                            free.push(crate::Motion::slide(unit(r, v, d)?));
                            continue;
                        }
                        (Some(l), Some("move")) if moved.is_none() => {
                            let [a, b, cc] = &l[1..] else { return Err(r.err_hint(d, "move takes three displacements (mm), x y z", Some(fix_usage.into()))) };
                            let mut v = [zero; 3];
                            for (k, e) in [a, b, cc].into_iter().enumerate() {
                                v[k] = r.field_of(e, Quantity::Length, "move")?;
                                if r.graph().depends_on(v[k], |x| x == Var::T) {
                                    return Err(r.err_hint(e, "a move may use x y z and params, not t", Some(fix_usage.into())));
                                }
                            }
                            moved = Some((v, d));
                            continue;
                        }
                        (_, Some("move")) => return Err(r.err_hint(d, "move given twice", Some(fix_usage.into()))),
                        (Some(l), Some("dofs")) if !seen_dofs => l,
                        (_, Some("dofs")) => return Err(r.err_hint(d, "dofs given twice", Some(fix_usage.into()))),
                        _ => {
                            let what = head.map_or_else(|| "a fix option must be a list".to_string(), |h| format!("unknown fix option '{h}'"));
                            return Err(r.err_hint(d, what, Some(fix_usage.into())));
                        }
                    };
                    seen_dofs = true;
                    dofs = [false; 3];
                    for a in &dl[1..] {
                        match a.sym() {
                            Some("x") => dofs[0] = true,
                            Some("y") => dofs[1] = true,
                            Some("z") => dofs[2] = true,
                            Some("normal") => normal = true,
                            _ => return Err(r.err(a, "dofs are x, y, z, or normal (only along the surface normal)")),
                        }
                    }
                    if normal && dofs.iter().any(|&b| b) {
                        return Err(r.err_hint(d, "(dofs normal) stands alone: it holds the one direction along the normal", Some(fix_usage.into())));
                    }
                    if !normal && !dofs.iter().any(|&b| b) {
                        return Err(r.err_hint(d, "(dofs) needs at least one of x, y, z, or normal", Some(fix_usage.into())));
                    }
                    if normal {
                        dofs = [true; 3];
                    }
                }
                // A move along an axis the fix leaves free would do nothing.
                if let Some((v, at)) = moved
                    && !normal
                    && let Some(k) = (0..3).find(|&k| !dofs[k] && r.graph().as_const(v[k]) != Some(0.0))
                {
                    return Err(r.err_hint(at, format!("the move is along {}, which the fix's (dofs ...) leave free", ["x", "y", "z"][k]), Some(fix_usage.into())));
                }
                any_move |= moved.is_some();
                let value = moved.map_or([zero; 3], |m| m.0);
                if rigid && free.is_empty() {
                    return Err(r.err_hint(fname, "(rigid) goes with pivot or slide: a fix without them holds its region rigidly already", Some(fix_usage.into())));
                }
                fixes.push(StudyFix { name: r.name(fname)?, region: Region::read(r, region)?, dofs, value, film: None, normal, free, rigid, boundary });
            }
            ("temperature", [fname, region, value]) => {
                // A held temperature may vary over the face (x y z and params), and in a transient in time.
                let v = r.field_of(value, Quantity::Temperature, "temperature")?;
                if timed_held.is_none() && r.graph().depends_on(v, |x| x == Var::T) {
                    timed_held = Some(value.clone());
                }
                let zero = r.graph().constant(0.0);
                fixes.push(StudyFix {
                    name: r.name(fname)?,
                    region: Region::read(r, region)?,
                    dofs: [true, false, false],
                    value: [v, zero, zero],
                    film: None,
                    boundary: crate::Boundary::Support,
                    normal: false,
                    free: Vec::new(),
                    rigid: false,
                });
            }
            ("convection", [fname, region, rest @ ..]) => {
                let (mut film, mut ambient) = (None, None);
                for (k, v) in pairs(r, rest)? {
                    match k.as_str() {
                        "coefficient" => film = Some(r.number_of(&v, Quantity::FilmCoefficient, "coefficient")?),
                        "ambient" => ambient = Some(r.number_of(&v, Quantity::Temperature, "ambient")?),
                        other => return Err(r.err(c, format!("unknown convection option '{other}' (coefficient ambient)"))),
                    }
                }
                let film = film.ok_or_else(|| r.err(c, "convection needs (coefficient W-per-m2-K)"))?;
                if !(film > 0.0) {
                    return Err(r.err(c, "the convection coefficient must be positive"));
                }
                let value = ambient.ok_or_else(|| r.err(c, "convection needs (ambient temperature)"))?;
                let (v, zero) = (r.graph().constant(value), r.graph().constant(0.0));
                fixes.push(StudyFix {
                    name: r.name(fname)?,
                    region: Region::read(r, region)?,
                    dofs: [true, false, false],
                    value: [v, zero, zero],
                    film: Some(film),
                    boundary: crate::Boundary::Support,
                    normal: false,
                    free: Vec::new(),
                    rigid: false,
                });
            }
            ("load" | "source", [lname, rest @ ..]) => loads.push(read_load(r, c, head, lname, rest, usage)?),
            ("gravity", [a, b, cc]) => {
                let mut n = |e| r.number_of(e, Quantity::Acceleration, "gravity");
                gravity = Some([n(a)?, n(b)?, n(cc)?])
            }
            ("res", [h]) => res = Some(r.cell_size(h)?),
            ("support-margin", [d]) => margin = Some(r.number_of(d, Quantity::Length, "support-margin")?),
            ("modes", [k]) => {
                let k = r.number(k)?;
                if !(1.0..=50.0).contains(&k) {
                    return Err(r.err(c, "modes must be between 1 and 50"));
                }
                modes = Some(k as usize);
            }
            _ => return Err(r.err_hint(c, format!("bad study clause '{head}'"), Some(usage.into()))),
        }
    }
    let shape = shape.ok_or_else(|| r.err(form, "study needs (shape field)"))?;
    let mut physics = physics.ok_or_else(|| r.err(form, "study needs (elastic ...) or (heat ...)"))?;
    let site = |key: &str| property_sites.get(key).or_else(|| property_sites.get("material")).unwrap_or(form);
    let density = match &mut physics {
        StudyPhysics::Elastic { e, nu, density, yield_strength, .. } => {
            r.guard_magnitude(site("E"), *e, Some(shape), MODULUS)?;
            r.guard_magnitude(site("nu"), *nu, Some(shape), POISSON)?;
            if let Some(y) = yield_strength {
                r.guard_magnitude(site("yield"), *y, Some(shape), YIELD)?;
            }
            density
        }
        StudyPhysics::Heat { k, density, .. } => {
            r.guard_magnitude(site("k"), *k, Some(shape), materials::CONDUCTIVITY)?;
            density
        }
    };
    if let Some(f) = density {
        r.guard_magnitude(site("density"), *f, Some(shape), materials::DENSITY)?;
        let kg_m3 = r.graph().constant(1e6);
        *f = r.graph().mul(*f, kg_m3);
    }
    let heat = matches!(physics, StudyPhysics::Heat { .. });
    if heat && fixes.iter().any(|f| matches!(f.boundary, crate::Boundary::Mirror { .. })) {
        return Err(r.err(form, "mirror declares an elastic symmetry plane; a thermal symmetry face already has zero normal heat flow"));
    }
    if !heat && fixes.iter().any(|f| f.film.is_some()) {
        return Err(r.err(form, "convection applies to heat studies"));
    }
    if heat && (any_move || fixes.iter().any(|f| f.normal || !f.free.is_empty())) {
        return Err(r.err(form, "(dofs normal), move, pivot and slide apply to elastic studies"));
    }
    // (A modal study's part may be free: its rigid motions are no modes.)
    if fixes.is_empty() && (heat || modes.is_none()) {
        let what = if heat { "(temperature name region value)" } else { "(fix name region)" };
        return Err(r.err(form, format!("study needs at least one {what}, or the part is free to float")));
    }
    if heat && gravity.is_some() {
        return Err(r.err(form, "gravity does not apply to heat studies"));
    }
    if gravity.is_some() && !matches!(physics, StudyPhysics::Elastic { density: Some(_), .. }) {
        return Err(r.err(form, "gravity needs (density ...) in (elastic ...), e.g. (density 1270kg/m3)"));
    }
    if modes.is_some() && !matches!(physics, StudyPhysics::Elastic { density: Some(_), .. }) {
        return Err(r.err(form, "(modes n) needs (elastic ... (density ...)), e.g. (density 1270kg/m3)"));
    }
    if let Some(tr) = &mut transient {
        match physics {
            StudyPhysics::Heat { density: Some(_), capacity: Some(c), .. } => tr.capacity = c,
            StudyPhysics::Heat { .. } => return Err(r.err(form, "transient heat needs (density ...) and (capacity J-per-kg-K) in (heat ...), e.g. (density 1270kg/m3)")),
            StudyPhysics::Elastic { .. } => return Err(r.err(form, "(transient ...) applies to heat studies")),
        }
    } else if loads.iter().any(|l| l.vector_t.is_some()) {
        return Err(r.err(form, "loads that depend on t need (transient (duration s) (steps n))"));
    } else if let Some(v) = &timed_held {
        return Err(r.err(v, "a held temperature that depends on t needs (transient (duration s) (steps n)): only a transient study steps through time"));
    }
    let moves = !heat && fixes.iter().flat_map(|f| f.value).any(|v| r.graph().as_const(v) != Some(0.0));
    if moves && modes.is_some() {
        return Err(r.err(form, "a fix's (move ...) applies to static elastic studies (a modal study takes no loads)"));
    }
    Ok(Study { name: n, shape, physics, fixes, loads, gravity, res, support_margin: margin, modes, transient, moved: moves, material_from })
}

/// The properties `keys` a `(material m)` gives (`m` its data), each said to
/// come from it as written in `kv`.
fn from_material<const N: usize>(m: &serde_json::Value, keys: [&str; N], kv: &[(String, Sexp)]) -> Vec<(String, String)> {
    let written = kv.iter().find(|(k, _)| k == "material").map_or_else(String::new, |(_, v)| format!("(material {v})"));
    keys.iter().filter(|k| m.get(**k).is_some_and(|v| v.is_number())).map(|k| (k.to_string(), written.clone())).collect()
}

/// A property written in the study: from there, as written.
fn said(from: &mut Vec<(String, String)>, key: &str, v: &Sexp) {
    if key == "material" {
        return;
    }
    from.retain(|(k, _)| k != key);
    from.push((key.to_string(), format!("({key} {v})")));
}
