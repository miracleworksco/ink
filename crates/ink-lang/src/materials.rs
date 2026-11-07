//! Engineering materials and fluids: the library std/materials.ink, plain IR
//! data (`(def petg (data (E 2100) (nu 0.38) ...))`). Each has one name:
//! another spelling (another case, spaces or underscores for hyphens, or a
//! name the library lists in `(def other-names (data (aluminium
//! aluminum-6061) ...))`) is refused with the name to write. The kernel
//! holds no table of its own; it reads the library's defs by name.

use super::*;
use ink_core::eval::eval_interval;
use ink_core::interval::I;
use std::sync::Arc;
use std::sync::OnceLock;

/// A published mass defines a uniform effective density over what has
/// none of its own: one body's (`(attr mass M)` on it), or the bodies an
/// assembly places (`(attr mass M)` on it), each placed as often as the
/// assembly places it. Bodies there with a density of their own (a steel
/// rotor) weigh as it says, and the rest share what remains. It is
/// measured at the configuration the mass is published for: the params'
/// defaults, or those `(attr mass-params ...)` gives, but the params the
/// mass reads, which keep their values (it is published at each); so a
/// size param scales the mass, as it scales the volume. A body's motion over t is a pose:
/// measured at each moment, the mass stays. A derived param, so every
/// consumer sees one density field; its volumes are measured, so claims
/// reading it are estimates, never proofs; one that does not settle still
/// sets it, and the answers say so.
#[derive(Clone)]
struct PublishedMass {
    name: String,
    /// The shapes sharing the density, each in its own coordinates, and
    /// how many times each is placed.
    shared: Vec<(Arc<Tape>, f64)>,
    /// Bodies with a density of their own: shape, density, times placed.
    own: Vec<(Arc<Tape>, Arc<Tape>, f64)>,
    mass: Tape,
    /// The published configuration: params and their values.
    at: Vec<(usize, f64)>,
    owns: [usize; 1],
    uses_t: bool,
}

impl Form for PublishedMass {
    fn name(&self) -> &str {
        &self.name
    }

    fn provider(&self, _: &Doc) -> Option<Arc<dyn ParamProvider>> {
        Some(Arc::new(self.clone()))
    }
}

impl ParamProvider for PublishedMass {
    fn name(&self) -> &str {
        &self.name
    }
    fn owns(&self) -> &[usize] {
        &self.owns
    }
    fn uses_t(&self) -> bool {
        self.uses_t
    }
    /// Measured where a mass is read, not by every command on the model.
    fn on_demand(&self) -> bool {
        true
    }
    fn warm_starts(&self) -> bool {
        false
    }

    fn provide(&self, params: &mut [f64], t: f64) -> Result<Option<String>, String> {
        let (ctx, mass) = self.at(params, t)?;
        let mut unsettled = None;
        let mut weigh = |shape: &Arc<Tape>, rho: Option<&Arc<Tape>>| -> Result<f64, String> {
            let (sum, why) = self.measure(shape, rho, &ctx)?;
            unsettled = unsettled.take().or(why);
            Ok(sum)
        };
        let own = self.own.iter().map(|(s, rho, n)| Ok(n * weigh(s, Some(rho))?)).sum::<Result<f64, String>>()?;
        let volume = self.shared.iter().map(|(s, n)| Ok(n * weigh(s, None)?)).sum::<Result<f64, String>>()?;
        if !(volume.is_finite() && volume > 0.0) {
            return Err(self.empty());
        }
        if own > mass {
            return Err(self.outweighed(own, mass));
        }
        params[self.owns[0]] = (mass - own) / volume;
        Ok(unsettled.map(|why| format!("a volume it rests on is not measured to {}%: {why}. The density is mass over the finest measurement reached", query::ACCURACY * 100.0)))
    }

    /// Without measuring the volume it shares where it can: its mass, its
    /// bodies' bounds, a volume above 0 (`holds_material`, else measured),
    /// and what its bodies with a density of their own weigh, measured only
    /// where their density's bound over their box could pass the mass.
    fn refusal(&self, params: &[f64], t: f64) -> Option<String> {
        let check = || -> Result<(), String> {
            let (ctx, mass) = self.at(params, t)?;
            // A volume above 0: a cell of one body proven inside, else
            // the bodies measured.
            let mut solid = false;
            for (s, _) in &self.shared {
                if let Some((lo, hi)) = self.bounds(s, &ctx)? {
                    solid = solid || holds_material(s, &ctx, lo, hi);
                }
            }
            if !solid {
                let volume = self.shared.iter().map(|(s, n)| Ok(n * self.measure(s, None, &ctx)?.0)).sum::<Result<f64, String>>()?;
                if !(volume.is_finite() && volume > 0.0) {
                    return Err(self.empty());
                }
            }
            let mut most = 0.0;
            for (s, rho, n) in &self.own {
                if let Some((lo, hi)) = self.bounds(s, &ctx)? {
                    let b: [I; 3] = std::array::from_fn(|k| I::new(lo[k], hi[k]));
                    let volume: f64 = (0..3).map(|k| hi[k] - lo[k]).product();
                    most += n * eval_interval(rho, &ctx, b, None).hi.max(0.0) * volume;
                }
            }
            if !(most * (1.0 + 1e-9) <= mass) {
                let own = self.own.iter().map(|(s, rho, n)| Ok(n * self.measure(s, Some(rho), &ctx)?.0)).sum::<Result<f64, String>>()?;
                if own > mass {
                    return Err(self.outweighed(own, mass));
                }
            }
            Ok(())
        };
        check().err()
    }
}

impl PublishedMass {
    /// The context it is measured in (its configuration) and its mass.
    fn at(&self, params: &[f64], t: f64) -> Result<(Ctx, f64), String> {
        let mut values = params.to_vec();
        for &(k, v) in &self.at {
            values[k] = v;
        }
        let ctx = Ctx::new(t, values);
        let mass = eval_point(&self.mass, &ctx, [0.0; 3]);
        if !mass.is_finite() || mass < 0.0 {
            return Err(format!("{}: mass must be finite and nonnegative (grams)", self.name));
        }
        Ok((ctx, mass))
    }

    /// A body's box (None: empty), or why a mass over it cannot be measured.
    fn bounds(&self, shape: &Arc<Tape>, ctx: &Ctx) -> Result<Option<(V3, V3)>, String> {
        let reach = crate::reach_of(shape, DEFAULT_DOMAIN, ctx);
        if !reach.whole {
            return Err(format!("{}: a published mass needs bounded bodies", self.name));
        }
        Ok(query::bounds(shape, ctx, reach.lo, reach.hi).map(|(lo, hi, _)| (lo, hi)))
    }

    /// A body's volume, or its mass under `rho`, as `ink mass` measures it;
    /// unsettled, the finest measurement reached, and why.
    fn measure(&self, shape: &Arc<Tape>, rho: Option<&Arc<Tape>>, ctx: &Ctx) -> Result<(f64, Option<String>), String> {
        let Some((lo, hi)) = self.bounds(shape, ctx)? else { return Ok((0.0, None)) };
        let v = query::measured(shape, &[rho.map(|r| &**r), None], true, ctx, lo, hi, None);
        Ok((v.sums[0], v.unsettled()))
    }

    fn empty(&self) -> String {
        format!("{}: a published mass needs a body with positive finite volume", self.name)
    }

    fn outweighed(&self, own: f64, mass: f64) -> String {
        format!("{}: its bodies with a density of their own weigh {own:.4} g, more than its published {mass} g", self.name)
    }
}

/// Whether a shape in `lo..hi` holds a cell of material (so a volume
/// above 0): a point the emptiness search finds inside, confirmed, and a
/// box round it a billionth of the shape's size bounded below 0. False
/// where that is not shown.
fn holds_material(shape: &Arc<Tape>, ctx: &Ctx, lo: V3, hi: V3) -> bool {
    let scale = query::max_extent(lo, hi);
    let floor = 1e-6 * scale;
    if !(floor > 0.0 && floor.is_finite()) {
        return false;
    }
    let query::Verdict::Disproven { witness, .. } = query::prove_empty(shape, ctx, lo, hi, None, floor, 0.0).verdict else { return false };
    let Some(p) = query::confirm_inside(shape, ctx, witness, scale) else { return false };
    let e = 1e-9 * scale;
    eval_interval(shape, ctx, p.map(|x| I::new(x - e, x + e)), None).all_inside()
}

/// Why an assembly's published mass is not used, where every body it
/// places has a density of its own.
pub const UNUSED_MASS: &str = "gives no body a density: every body it places has its own (a material, a density, or an inner assembly's mass), \
    so they weigh what those say and this mass is not used";

pub const LIB: &str = include_str!("../../../std/materials.ink");

/// The library's def of the names it refuses: data, not a material.
const OTHER_NAMES: &str = "other-names";

/// The library's entries, (name, properties) in library order, and the
/// names it refuses, (name, the one to write).
struct Library {
    entries: Vec<(String, serde_json::Value)>,
    instead: Vec<(String, String)>,
}

fn library() -> &'static Library {
    static T: OnceLock<Library> = OnceLock::new();
    T.get_or_init(|| {
        let forms = parse(LIB).expect("std/materials.ink parses");
        let mut lib = Library { entries: Vec::new(), instead: Vec::new() };
        for f in &forms {
            // (def name (data (key value) ...)): pairs of a symbol and a value.
            let Some([def, name, value]) = f.list() else { panic!("std/materials.ink: {f}: every form is (def name (data ...))") };
            let (Some("def"), Some(name), Some([data, pairs @ ..])) = (def.sym(), name.sym(), value.list()) else {
                panic!("std/materials.ink: {f}: every form is (def name (data ...))")
            };
            assert!(data.sym() == Some("data"), "std/materials.ink: {name}: every form is (def name (data ...))");
            let pairs = pairs.iter().filter_map(|p| match p.list() {
                Some([k, v]) => Some((k.sym()?, v)),
                _ => None,
            });
            if name == OTHER_NAMES {
                lib.instead.extend(pairs.filter_map(|(k, to)| Some((k.to_string(), to.sym()?.to_string()))));
            } else {
                let m = pairs.filter_map(|(k, v)| Some((k.to_string(), serde_json::json!(v.num()?)))).collect();
                lib.entries.push((name.to_string(), serde_json::Value::Object(m)));
            }
        }
        for (name, to) in &lib.instead {
            assert!(lib.entries.iter().any(|(n, _)| n == to), "std/materials.ink: {OTHER_NAMES}: ({name} {to}): no material {to}");
        }
        lib
    })
}

/// Solid materials' names (those with a Young's modulus), for hints.
pub fn names() -> Vec<&'static str> {
    library().entries.iter().filter(|(_, v)| v.get("E").is_some()).map(|(n, _)| n.as_str()).collect()
}

/// Fluids' names (those with a viscosity), for hints.
pub fn fluids() -> Vec<&'static str> {
    library().entries.iter().filter(|(_, v)| v.get("viscosity").is_some()).map(|(n, _)| n.as_str()).collect()
}

/// Properties of a material or fluid as a JSON object, with its `name`.
pub fn lookup(name: &str) -> Option<serde_json::Value> {
    library().entries.iter().find(|(n, _)| n == name).map(|(n, v)| {
        let mut v = v.clone();
        v["name"] = serde_json::json!(n);
        v
    })
}

/// A body's property `key` from its `(attr material ...)` (a library name,
/// or data with the keys), if it has one with that key.
pub fn material_property(attrs: &[(String, Value)], key: &str) -> Option<f64> {
    let m = match &attrs.iter().find(|(k, _)| k == "material")?.1 {
        Value::Str(n) => lookup(n)?,
        Value::Data(d) => d.value().clone(),
        _ => return None,
    };
    m.get(key)?.as_f64()
}

/// A body's attribute `key`, as every reader takes it (a claim's `(attr
/// key)`, `(by-body ...)`, `ink mass`): its own
/// `(attr key ...)`, else its material's `key`. None where it has neither;
/// an error, naming the target as `who` (Doc::named), where its own is not
/// a number or a field.
pub fn attr(g: &mut Graph, who: &str, attrs: &[(String, Value)], key: &str) -> Result<Option<NodeId>, InkError> {
    match attrs.iter().find(|(k, _)| k == key) {
        Some((_, Value::Field(f, _))) => Ok(Some(*f)),
        Some((_, other)) => Err(InkError::msg(format!("attr '{key}' of {who} is a {}, not a number or a field", other.kind()))),
        None => Ok(material_property(attrs, key).map(|v| g.constant(v))),
    }
}

/// `attr`, where the body must have it: else `missing`.
pub fn required(g: &mut Graph, who: &str, attrs: &[(String, Value)], key: &str) -> Result<NodeId, InkError> {
    attr(g, who, attrs, key)?.ok_or_else(|| missing(who, attrs, key))
}

/// The error every reader gives for a target `who` (Doc::named) whose
/// body lacks attribute `key`: the attributes it has, and for a density,
/// what gives one.
pub fn missing(who: &str, attrs: &[(String, Value)], key: &str) -> InkError {
    let names: Vec<&str> = attrs.iter().map(|(a, _)| a.as_str()).collect();
    let has = if names.is_empty() { "none".to_string() } else { names.join(", ") };
    let give = if key == "density" {
        ". Give its body (attr material NAME) or (attr density 7.85g/cm3) in the file that defines it, or give an assembly placing it a published (attr mass M)"
    } else {
        ""
    };
    InkError { hint: Some(format!("its attrs: {has}{give}")), ..InkError::msg(format!("{who} has no attr '{key}'")) }
}

/// A density's unit is g/mm^3 everywhere (a body's attribute, the
/// library, studies, fluids), so a body's mass over its mm^3 is in grams;
/// a number may say its unit, `7.85g/cm3`, `7850kg/m3` (sexp::UNITS). Every
/// real material lies between hydrogen's 9e-8 g/mm^3 and osmium's 0.0226,
/// and a density written in kg/m^3 or g/cm^3 lies outside that (but foams
/// under 23 kg/m^3 in g/cm^3), as does a kg/m^3 converted twice: so a
/// number outside it is refused, never read a thousand or a million times
/// off. 0 is massless.
pub const DENSITY_RANGE: (f64, f64) = (5e-8, 0.023);

/// Why `d` cannot be a density in g/mm^3, if it cannot.
pub fn density_error(d: f64) -> Option<String> {
    let (lo, hi) = DENSITY_RANGE;
    if d == 0.0 || (lo..=hi).contains(&d) {
        return None;
    }
    let why = if d < 0.0 {
        "is negative".to_string()
    } else if d > hi {
        format!("is denser than any material in g/mm^3, ink's density unit (osmium is 0.0226). If it is kg/m^3, write {d}kg/m3; if g/cm^3, {d}g/cm3")
    } else {
        "is lighter than any material in g/mm^3, ink's density unit (hydrogen is 9e-8). Converted twice? Write it with its unit, e.g. 7850kg/m3".to_string()
    };
    Some(format!("density {d} {why}"))
}

/// A density field's error, where it is a number (fields are not bounded).
pub fn density_field_error(g: &Graph, f: NodeId) -> Option<String> {
    g.as_const(f).and_then(density_error)
}

/// A material's density (g/mm^3), refused if it cannot be one: the
/// library's are, but `(data (density 7850))` is not.
pub fn material_density_error(m: &serde_json::Value) -> Option<String> {
    m.get("density")?.as_f64().and_then(density_error).map(|e| format!("material's {e}"))
}

/// A command's `--density`: a number (which may say its unit, `7.85g/cm3`),
/// refused as `density_error` refuses it, or None for an attribute's name.
pub fn density_option(s: &str) -> Option<Result<f64, String>> {
    if !s.starts_with(|c: char| c.is_ascii_digit() || matches!(c, '.' | '-' | '+')) {
        return None;
    }
    Some(crate::sexp::number_of(s, crate::units::Quantity::Density, "this").and_then(|d| density_error(d).map_or(Ok(d), Err)))
}

/// The density `ink mass` weighs target `who`
/// (Doc::named) by, and what it is (for the answer): the attribute asked
/// for by name (`--density NAME`), which it must have (`required`); else
/// its `density`, if it has one (None: it has no mass).
pub fn mass_density(d: &mut Doc, who: &str, attrs: &[(String, Value)], asked: Option<&str>) -> Result<Option<(NodeId, String)>, InkError> {
    let key = asked.unwrap_or("density");
    let f = match asked {
        Some(_) => required(&mut d.g, who, attrs, key)?,
        None => match attr(&mut d.g, who, attrs, key)? {
            Some(f) => f,
            None => return Ok(None),
        },
    };
    let published = d.published.iter().find(|(k, _)| d.g.var(Var::Param(*k as u32)) == f).map(|(_, what)| what.clone());
    let has = |k: &str| attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v);
    let what = match (published, has(key), has("material"), d.g.as_const(f)) {
        (Some(what), ..) if key == "density" => format!(
            "{what}, measured to {}% as ink mass measures a volume, at the configuration it is published for (the params' defaults, or its mass-params, but those the mass reads); uniform density: \
             centre of mass is the centroid, inertia assumes that uniform distribution",
            query::ACCURACY * 100.0
        ),
        (_, None, Some(m), Some(d)) => {
            let name = if let Value::Str(n) = m { n.as_str() } else { "its material" };
            if key == "density" { format!("{name}: {d} g/mm^3") } else { format!("{name}'s {key}: {d}") }
        }
        _ => format!("attr '{key}'"),
    };
    Ok(Some((f, what)))
}

/// The name to write for one the library does not take: the same name in
/// another case or with spaces or underscores for hyphens, or a name it
/// lists in `other-names`.
pub fn instead(name: &str) -> Option<&'static str> {
    let key = name.to_ascii_lowercase().replace(['_', ' '], "-");
    let lib = library();
    let own = lib.entries.iter().map(|(n, _)| (n, n));
    own.chain(lib.instead.iter().map(|(n, to)| (n, to))).find(|(n, _)| **n == key).map(|(_, to)| to.as_str())
}

/// A name that is not a material: what is wrong, and what to write.
pub fn refused(name: &str) -> (String, String) {
    match instead(name) {
        Some(to) => (format!("material '{name}' is spelled '{to}'"), format!("write {to}: each material has one name in std/materials.ink")),
        None => (format!("unknown material '{name}'"), format!("materials: {}", names().join(" "))),
    }
}

impl Doc {
    /// A target as a reader names it: its path, and the body it places,
    /// with the file that defines it where that is another (an import).
    pub fn named(&self, it: &Item) -> String {
        match it.instance.as_ref().and_then(|(b, _)| self.body(b)) {
            Some(b) if *b.file != *self.file => format!("'{}' (body '{}' in {})", it.label, b.name, b.file),
            Some(b) if b.name != it.label => format!("'{}' (body '{}')", it.label, b.name),
            _ => format!("'{}'", it.label),
        }
    }

    /// `(attr mass M)` on a body: its published mass in grams, spread
    /// uniformly over it (PublishedMass).
    /// Built with arguments (`of` its definition) and the same mass, it
    /// takes the definition's density, measured at its defaults: a size
    /// passed as an argument scales the mass as a param would.
    pub(super) fn published_mass(&mut self, at: &Sexp, name: &str, of: Option<&str>, shape: NodeId, attrs: &mut Vec<(String, Value)>) -> R<()> {
        let Some((_, mass)) = attrs.iter().find(|(k, _)| k == "mass") else { return Ok(()) };
        let mass = self.total_mass(at, mass)?;
        if attrs.iter().any(|(k, _)| k == "density") || material_property(attrs, "density").is_some() {
            return Err(self.err_hint(
                at,
                "attr mass defines density, so it cannot accompany a density or a material with density",
                Some("for a solid with a known material, keep material density and compare its measured mass with the published mass".into()),
            ));
        }
        let defined = of.and_then(|d| self.body(d)).map(|b| {
            let get = |k: &str| b.attrs.iter().find(|(a, _)| a == k).map(|(_, v)| v.clone());
            (get("mass"), get("density"))
        });
        if let Some((Some(Value::Field(m, _)), Some(density))) = defined
            && m == mass
        {
            attrs.push(("density".into(), density));
            return Ok(());
        }
        let config = self.published_config(at, attrs.iter().find(|(k, _)| k == "mass-params").map(|(_, v)| v))?;
        let k = self.publish(at, name, format!("published mass of '{name}' over its volume"), mass, vec![(shape, 1.0)], Vec::new(), config)?;
        attrs.push(("density".into(), Value::Field(self.g.var(Var::Param(k as u32)), None)));
        Ok(())
    }

    /// A body built with arguments, `of` a definition an assembly's
    /// published mass gave a density (a catalog part's root weighing its
    /// bodies), with no density, material or mass of its own: it takes that
    /// density, measured where the mass is published, so a size passed as
    /// an argument scales the mass as a param would.
    pub(super) fn weighed_as_defined(&self, mut body: Body) -> Body {
        let own = |a: &[(String, Value)]| a.iter().any(|(k, _)| k == "density" || k == "mass") || material_property(a, "density").is_some();
        let Some(def) = body.of.as_deref().and_then(|d| self.body(d)) else { return body };
        if own(&body.attrs) {
            return body;
        }
        if let Some((_, rho)) = def.attrs.iter().find(|(k, _)| k == "density") {
            let mut attrs = (*body.attrs).clone();
            attrs.push(("density".into(), rho.clone()));
            body.attrs = Attrs::new(attrs);
        }
        body
    }

    /// `(attr mass M)` on assembly `asm`, its instances read: the bodies it
    /// places, each as often as it does, weigh M together. Those with a
    /// density of their own keep it; the rest share one (PublishedMass).
    pub(super) fn assembly_mass(&mut self, at: &Sexp, asm: &str, mass: &Value, params: Option<&Value>) -> R<()> {
        let mass = self.total_mass(at, mass)?;
        let mut placed: Vec<(String, f64)> = Vec::new();
        self.count_placed(asm, &mut placed);
        let (mut shared, mut own, mut names) = (Vec::new(), Vec::new(), Vec::new());
        // A body with a density when this is read keeps it: its own, its
        // material's, or a share of an assembly inside this one (a
        // catalog part's root): this one's mass weighs the rest.
        for (b, n) in placed {
            let (shape, attrs) = self.body(&b).map(|x| (x.shape, x.attrs.clone())).unwrap();
            match attr(&mut self.g, &format!("'{b}'"), &attrs, "density").map_err(|e| self.err(at, e.message))? {
                Some(rho) => own.push((shape, rho, n)),
                None => {
                    shared.push((shape, n));
                    names.push(b);
                }
            }
        }
        // A mass that gives no body a density is not used: said, not
        // ignored (one lighter than what has its own is an error where
        // weighed: no density can make up the difference).
        if shared.is_empty() {
            let w = format!("{}:{}: the (attr mass ...) of assembly '{asm}' {UNUSED_MASS}. Remove it, or the density of the bodies it should weigh", self.file, at.line);
            if !self.load_warnings.contains(&w) {
                self.load_warnings.push(w);
            }
            self.keep_mass(asm, AssemblyMass { mass, density: None });
            return Ok(());
        }
        let config = self.published_config(at, params)?;
        let what = match own.is_empty() {
            true => format!("published mass of '{asm}' over its bodies' volume"),
            false => format!("published mass of '{asm}', less what its bodies with a density of their own weigh, over the others' volume"),
        };
        let k = self.publish(at, asm, what, mass, shared, own, config)?;
        let f = self.g.var(Var::Param(k as u32));
        for b in names {
            let body = self.bodies.iter_mut().find(|x| x.name == b).unwrap();
            let mut attrs = (*body.attrs).clone();
            attrs.push(("density".into(), Value::Field(f, None)));
            body.attrs = Attrs::new(attrs);
        }
        self.keep_mass(asm, AssemblyMass { mass, density: Some(f) });
        Ok(())
    }

    fn keep_mass(&mut self, asm: &str, mass: AssemblyMass) {
        self.assemblies.iter_mut().find(|a| a.name == asm).unwrap().mass = Some(mass);
    }

    /// The bodies assembly `asm` places, through its sub-assemblies, and
    /// how many times each.
    fn count_placed(&self, asm: &str, out: &mut Vec<(String, f64)>) {
        for i in &self.assembly(asm).unwrap().instances {
            if self.body(&i.child).is_some() {
                match out.iter_mut().find(|(b, _)| *b == i.child) {
                    Some((_, n)) => *n += 1.0,
                    None => out.push((i.child.clone(), 1.0)),
                }
            } else {
                self.count_placed(&i.child, out);
            }
        }
    }

    /// A published mass's value: grams, the same everywhere.
    fn total_mass(&self, at: &Sexp, mass: &Value) -> R<NodeId> {
        let Value::Field(mass, _) = *mass else { return Err(self.err(at, "attr mass must be a number in grams")) };
        if self.g.depends_on(mass, |v| matches!(v, Var::X | Var::Y | Var::Z)) {
            return Err(self.err(at, "attr mass is a total mass, so it cannot vary with x, y or z"));
        }
        // A number at load is checked there; what the mass is spread over
        // is measured only where a mass is read (PublishedMass::on_demand).
        if !self.g.depends_on(mass, |v| v == Var::T) {
            let m = self.load_eval(mass);
            if !m.is_finite() || m < 0.0 {
                return Err(self.err(at, format!("attr mass must be finite and nonnegative (grams), not {m}")));
            }
        }
        Ok(mass)
    }

    /// The configuration a mass is published for: every param read so far
    /// (but those a solve sets) at its default, or at the value
    /// `(attr mass-params (read-json FILE))` gives it, FILE a JSON object
    /// of param values.
    fn published_config(&self, at: &Sexp, given: Option<&Value>) -> R<Vec<(usize, f64)>> {
        let solved: HashSet<usize> = self.solves.iter().flat_map(|s| s.owns.iter().copied()).collect();
        let mut config: Vec<(usize, f64)> = (0..self.params.len()).filter(|k| !solved.contains(k)).map(|k| (k, self.params[k].default)).collect();
        let usage = "attr mass-params is the configuration the mass is published for: (read-json FILE), FILE a JSON object of param values, e.g. {\"length\": 1193.8}";
        let Some(given) = given else { return Ok(config) };
        let Some(given) = (match given {
            Value::Data(d) => d.value().as_object().cloned(),
            _ => None,
        }) else {
            return Err(self.err_hint(at, "bad attr mass-params", Some(usage.into())));
        };
        for (name, v) in given {
            let (Some(k), Some(v)) = (self.param_index(&name), v.as_f64()) else {
                return Err(self.err_hint(at, format!("attr mass-params: '{name}' is not a param with a number"), Some(usage.into())));
            };
            config.retain(|(j, _)| *j != k);
            config.push((k, v));
        }
        Ok(config)
    }

    /// The derived density `NAME.mass-density` of a published mass over
    /// `shared`, less what `own` weighs (PublishedMass): its param.
    fn publish(&mut self, at: &Sexp, name: &str, what: String, mass: NodeId, shared: Vec<(NodeId, f64)>, own: Vec<(NodeId, NodeId, f64)>, config: Vec<(usize, f64)>) -> R<usize> {
        let k = self.params.len();
        let derived = format!("{name}.mass-density");
        // Not a number until measured: a reader that misses Doc::want
        // weighs nothing it could take for a mass (0 is a density).
        self.params.push(Param { name: derived.clone(), default: f64::NAN, tol: None });
        let reads: Vec<NodeId> = std::iter::once(mass).chain(shared.iter().map(|s| s.0)).chain(own.iter().flat_map(|o| [o.0, o.1])).collect();
        let uses_t = reads.iter().any(|&r| self.g.depends_on(r, |v| v == Var::T));
        let tape = |g: &Graph, n: NodeId| Arc::new(Tape::new(g, n));
        // The params the mass reads keep their values: it is published at
        // each (a toleranced mass, a mass per length).
        let follows = self.params_read(&[mass]);
        let p = PublishedMass {
            name: derived,
            shared: shared.iter().map(|&(s, n)| (tape(&self.g, s), n)).collect(),
            own: own.iter().map(|&(s, rho, n)| (tape(&self.g, s), tape(&self.g, rho), n)).collect(),
            mass: Tape::new(&self.g, mass),
            at: config.into_iter().filter(|(k, _)| !follows.contains(k)).collect(),
            owns: [k],
            uses_t,
        };
        self.keep_form("mass", "mass", at, Arc::new(p), Reads { fields: reads, forms: vec![] })?;
        self.published.push((k, what));
        Ok(k)
    }

    /// A material: a name from the library (symbol or string), or data with
    /// the property keys.
    pub(super) fn material_of(&mut self, m: &Sexp, env: &mut Env) -> R<serde_json::Value> {
        let name = match &m.k {
            SK::Sym(s) | SK::Str(s) => Some(s.clone()),
            _ => None,
        };
        if let Some(n) = &name
            && let Some(v) = lookup(n)
        {
            return Ok(v);
        }
        let refuse = |d: &Self, n: &str| {
            let (what, hint) = refused(n);
            d.err_hint(m, what, Some(hint))
        };
        if let Some(n) = name.as_deref().filter(|n| self.def(n).is_none() && env.iter().all(|(k, _)| k != n)) {
            return Err(refuse(self, n));
        }
        match self.eval(m, env)? {
            Value::Data(d) if d.value().is_object() => match material_density_error(d.value()) {
                Some(e) => Err(self.err(m, e)),
                None => Ok(d.value().clone()),
            },
            Value::Str(s) => lookup(&s).ok_or_else(|| refuse(self, &s)),
            _ => Err(self.err_hint(m, "a material is a library name or (material name)", Some(format!("materials: {}", names().join(" "))))),
        }
    }

    /// `(material name)` and `(by-body target key [default])`.
    pub(super) fn eval_materials(&mut self, name: &str, e: &Sexp, args: &[Sexp], env: &mut Env) -> R<Value> {
        match name {
            "material" => {
                let [m] = args else { return Err(self.err_hint(e, "bad material", Some(format!("usage: (material name); names: {}", names().join(" "))))) };
                return Ok(Value::Data(crate::Data::new(self.material_of(m, env)?)));
            }
            "by-body" => {
                // A field taking each body's property where that body is:
                // (by-body target key [default]), the body's attribute `key`.
                let usage = "usage: (by-body target key [default]): a field with each body's property where it is (attr key, else its material's key)";
                let (target, key, default) = match args {
                    [t, k] => (t, k, None),
                    [t, k, d] => (t, k, Some(self.field(d, env)?)),
                    _ => return Err(self.err_hint(e, "bad by-body", Some(usage.into()))),
                };
                let target = self.path_of(target)?;
                let key = self.name_of(key)?;
                let items = self.items(&target).map_err(|m| self.err(e, m))?;
                // Each body's own value, else the default; the nearest body
                // (smallest field value) wins, inside and out: solvers
                // sample properties slightly past surfaces, where "inside
                // this body, else the rest" would give the wrong
                // neighbour's value.
                let mut acc: Option<(NodeId, NodeId)> = None;
                for it in items {
                    let who = self.named(&it);
                    let (shape, attrs) = (self.placed_field(&it), self.placed_attrs(&it));
                    let value = match default {
                        Some(d) => attr(&mut self.g, &who, &attrs, &key).map(|f| f.unwrap_or(d)),
                        None => required(&mut self.g, &who, &attrs, &key),
                    };
                    let value = value.map_err(|m| self.err_hint(e, m.message, m.hint.or(Some(usage.into()))))?;
                    acc = Some(match acc {
                        None => (value, shape),
                        Some((a, near)) => {
                            let c = self.g.sub(shape, near);
                            (self.g.select(c, value, a), self.g.min(near, shape))
                        }
                    });
                }
                let Some((f, _)) = acc else { return Err(self.err(e, format!("'{target}' has no bodies"))) };
                self.want(&[f]);
                return Ok(Value::Field(f, None));
            }
            _ => {}
        }
        unreachable!()
    }
}
