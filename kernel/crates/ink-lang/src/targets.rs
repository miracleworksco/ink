//! Targets: what a name asked of a document is (its items, each a placed
//! shape with attributes), and placing fields and attributes.

use super::*;

impl Doc {
    fn field_items(name: &str, value: &Value) -> Result<Vec<Item>, String> {
        match value {
            Value::Field(f, _) => Ok(vec![Item { label: name.into(), shape: *f, attrs: Attrs::default(), instance: None, links: Vec::new() }]),
            other => Err(format!("def '{name}' is a {}, not a field", other.kind())),
        }
    }

    /// An expression's target, before its source name is qualified. A
    /// local field def can have the same spelling as a caller's body.
    pub(super) fn scoped_items(&mut self, written: &str) -> Result<(String, Vec<Item>), String> {
        if !self.files[self.scope.file as usize].objects.contains_key(written)
            && !written.contains('/')
            && let Some(v) = self.def(written)
        {
            return Ok((written.to_string(), Self::field_items(written, v)?));
        }
        let path = self.resolve(written)?;
        let items = self.items(&path)?;
        Ok((path, items))
    }

    /// Resolve a body, assembly, def, or instance path (`asm/inst/...`) to
    /// its items: each body's own field and its placement (Item).
    pub fn items(&mut self, target: &str) -> Result<Vec<Item>, String> {
        let parts: Vec<&str> = target.split('/').collect();
        if parts.len() == 1 {
            if let Some(b) = self.body(target) {
                let (name, shape, attrs) = (b.name.clone(), b.shape, b.attrs.clone());
                let id = self.identity();
                return Ok(vec![Item { label: name.clone(), shape, attrs, instance: Some((name, id)), links: Vec::new() }]);
            }
            if self.assembly(target).is_some() {
                let id = self.identity();
                return self.flatten(target, &id, target, &[]);
            }
            if let Some(v) = self.def(target) {
                return Self::field_items(target, v);
            }
            let names = self.bodies.iter().map(|b| b.name.as_str()).chain(self.assemblies.iter().map(|a| a.name.as_str()));
            let hint = closest(target, names).unwrap_or_default();
            return Err(format!("no body, assembly, or field def named '{target}'. {hint}").trim().to_string());
        }
        let mut asm_name = parts[0].to_string();
        if self.assembly(&asm_name).is_none() {
            return Err(format!("'{asm_name}' is not an assembly"));
        }
        let mut xf = self.identity();
        // The instances from the world down to the target: an item's links
        // start with them, as its world placement does (a placement above
        // the target moves and scales what is under it too).
        let mut above: Vec<Link> = Vec::new();
        for (k, seg) in parts[1..].iter().enumerate() {
            let asm = self.assembly(&asm_name).unwrap();
            let inst = asm
                .instances
                .iter()
                .find(|i| i.name == *seg)
                .ok_or_else(|| format!("assembly '{asm_name}' has no instance '{seg}'. {}", choices(seg, "Its instances", asm.instances.iter().map(|i| i.name.as_str()))))?;
            let (child, at, joint, on) = (inst.child.clone(), inst.at, inst.joint.clone(), inst.on.clone());
            xf = self.compose(&at, &xf);
            above.push(Link { path: parts[..k + 2].join("/"), at, joint, on });
            let last = k == parts.len() - 2;
            if let Some(b) = self.body(&child) {
                if !last {
                    return Err(format!("'{seg}' is a body, with no instances"));
                }
                return Ok(vec![Item { label: target.to_string(), shape: b.shape, attrs: b.attrs.clone(), instance: Some((child, xf)), links: above }]);
            }
            if last {
                return self.flatten(&child, &xf, target, &above);
            }
            asm_name = child;
        }
        unreachable!()
    }

    /// The items of assembly `asm` placed at `parent`, labelled under
    /// `prefix`, each with the instances above it (`links`, then its own).
    pub(super) fn flatten(&mut self, asm: &str, parent: &Xf, prefix: &str, links: &[Link]) -> Result<Vec<Item>, String> {
        let a = self.assembly(asm).unwrap();
        let insts: Vec<(String, String, Xf, Option<Joint>, Vec<String>)> =
            a.instances.iter().map(|i| (i.name.clone(), i.child.clone(), i.at, i.joint.clone(), i.on.clone())).collect();
        let mut out = Vec::new();
        for (name, child, at, joint, on) in insts {
            let xf = self.compose(&at, parent);
            let label = format!("{prefix}/{name}");
            let mut path = links.to_vec();
            path.push(Link { path: label.clone(), at, joint, on });
            if let Some(b) = self.body(&child) {
                out.push(Item { label, shape: b.shape, attrs: b.attrs.clone(), instance: Some((child, xf)), links: path });
            } else {
                out.extend(self.flatten(&child, &xf, &label, &path)?);
            }
        }
        Ok(out)
    }

    /// A part's declared tolerance, its `(attr tolerance x)`: how far the
    /// made part may differ from its shape (claims::finest).
    pub fn tolerance(&self, attrs: &[(String, Value)]) -> Option<f64> {
        match attrs.iter().find(|(k, _)| k == "tolerance")? {
            (_, Value::Field(f, _)) => self.g.as_const(*f),
            _ => None,
        }
    }

    /// The parts field `f` is made of, for its accuracy (claims::finest):
    /// each body whose shape it reaches at any depth by the tag every body's
    /// shape carries, through unions, assemblies, placements, warps and the
    /// shapes of other bodies alike, with its declared tolerance (None: it
    /// declares none); and None for material of no body, a coordinate read
    /// outside every body's tag. A tag only adds a part, so naming one never
    /// hides another.
    pub fn parts_of(&self, f: NodeId) -> Vec<Option<(String, f64)>> {
        use crate::graph::Node;
        let mut bodies: HashMap<u32, Vec<&Body>> = HashMap::new();
        for b in &self.bodies {
            if let Node::Tag(t, _) = self.g.nodes[b.shape as usize] {
                bodies.entry(t).or_default().push(b);
            }
        }
        // Each node with whether a body's tag is above it.
        let (mut out, mut seen, mut tags, mut stack) = (Vec::new(), HashSet::new(), HashSet::new(), vec![(f, false)]);
        let mut loose = false;
        while let Some((n, inside)) = stack.pop() {
            if !seen.insert((n, inside)) {
                continue;
            }
            match self.g.nodes[n as usize] {
                Node::Tag(t, a) if let Some(bs) = bodies.get(&t) => {
                    if tags.insert(t) {
                        out.extend(bs.iter().map(|b| self.tolerance(&b.attrs).map(|x| (b.name.clone(), x))));
                    }
                    stack.push((a, true));
                }
                Node::Var(Var::X | Var::Y | Var::Z) => loose |= !inside,
                node => {
                    let (ch, k) = node.children();
                    stack.extend(ch[..k].iter().map(|&c| (c, inside)));
                }
            }
        }
        if loose {
            out.push(None);
        }
        out
    }

    /// An item's field where it is placed (`warped`), recorded as its
    /// body's field and placement (`Doc::placed`).
    /// What world field `w` was placed from, where it was (`Doc::placed`,
    /// or a field warped by a placement: a zone, an attribute): its field
    /// in its own frame and the placement.
    pub fn placement_of(&self, w: NodeId) -> Option<Placed> {
        self.placed
            .get(&w)
            .copied()
            .or_else(|| self.placed_fields.iter().filter(|(k, x)| **x == w && k.0 != w).map(|(k, _)| Placed { shape: k.0, at: k.1 }).min_by_key(|p| p.shape))
    }

    /// Every world field placed (`placement_of`), in node order.
    pub fn placements(&self) -> Vec<(NodeId, Placed)> {
        let mut out: Vec<(NodeId, Placed)> = self.placed.iter().map(|(w, p)| (*w, *p)).collect();
        out.extend(self.placed_fields.iter().filter(|((f, _), w)| *w != f && !self.placed.contains_key(w)).map(|(k, w)| (*w, Placed { shape: k.0, at: k.1 })));
        // A named field placed is its name on the field placed (a placement is a substitution, which
        // keeps a tag where it is), and what is made of it (max with a slab) reads the field without
        // the name: that is the unnamed shape placed the same way.
        let mut i = 0;
        while i < out.len() {
            let (w, p) = out[i];
            if let (crate::graph::Node::Tag(t, inner), crate::graph::Node::Tag(u, shape)) = (self.g.nodes[w as usize], self.g.nodes[p.shape as usize])
                && t == u
                && inner != shape
            {
                out.push((inner, Placed { shape, at: p.at }));
            }
            i += 1;
        }
        out.sort_unstable_by_key(|x| (x.0, x.1.shape));
        out.dedup_by_key(|x| x.0);
        out
    }

    pub fn placed_field(&mut self, it: &Item) -> NodeId {
        let Some((_, xf)) = &it.instance else { return it.shape };
        let f = self.warped(it.shape, xf);
        self.record_placed(f, it.shape, xf);
        self.placed_labels.entry(f).or_insert_with(|| it.label.clone());
        f
    }

    /// A field of an item's frame (an attribute) where the item is placed:
    /// one reading none of x, y and z is itself.
    pub fn place(&mut self, f: NodeId, it: &Item) -> NodeId {
        match &it.instance {
            Some((_, xf)) if self.xyz_reads(f).contains(&true) => self.warped(f, xf),
            _ => f,
        }
    }

    /// `f` warped by placement `xf`, once per field and placement (by the
    /// identity: itself).
    pub fn warped(&mut self, f: NodeId, xf: &Xf) -> NodeId {
        if let Some(&w) = self.placed_fields.get(&(f, *xf)) {
            return w;
        }
        let map = self.inverse_map(xf);
        let identity = map.iter().all(|&(v, u)| self.g.motion_reads(u).is_empty() && matches!(self.g.nodes[u as usize], crate::graph::Node::Var(w) if w == v));
        let w = if identity { f } else { self.g.subst(f, &map) };
        self.placed_fields.insert((f, *xf), w);
        w
    }

    /// A value as an error message shows it: text as written, a number (a
    /// param's at its value) in its shortest form (24, 2.5), a vec's items, a
    /// record as JSON, anything else by kind.
    pub(super) fn value_text(&self, v: &Value) -> String {
        match v {
            Value::Str(s) => s.clone(),
            Value::Field(f, _) => self.load_value(*f).map_or_else(|| "(a field of x, y, z or t)".into(), |x| format!("{x}")),
            Value::Vec(items) => format!("[{}]", items.iter().map(|x| self.value_text(x)).collect::<Vec<_>>().join(", ")),
            Value::Data(d) => d.value().to_string(),
            other => format!("(a {})", other.kind()),
        }
    }

    /// An item's attributes where it is placed: their fields `place`d.
    pub fn placed_attrs(&mut self, it: &Item) -> Attrs {
        fn placed(doc: &mut Doc, v: &Value, it: &Item) -> Value {
            match v {
                Value::Field(f, said) => {
                    let said = said.as_ref().map(|s| match s.explicit_when {
                        Some(w) => Rc::new(units::Said { explicit_when: Some(doc.place(w, it)), ..(**s).clone() }),
                        None => s.clone(),
                    });
                    Value::Field(doc.place(*f, it), said)
                }
                Value::Vec(items) => Value::list(items.iter().map(|x| placed(doc, x, it)).collect()),
                other => other.clone(),
            }
        }
        Attrs::new(it.attrs.iter().map(|(k, v)| (k.clone(), placed(self, v, it))).collect())
    }

    /// Fields reading what the items' placed fields read (params, samples,
    /// t), none of them placed: each item's own field and, for each of x, y
    /// and z it reads, what its placement puts there.
    pub fn placed_reads(&mut self, items: &[Item]) -> Vec<NodeId> {
        let mut out = Vec::new();
        for it in items {
            out.push(it.shape);
            let Some((_, xf)) = &it.instance else { continue };
            let reads = self.xyz_reads(it.shape);
            if reads.contains(&true) {
                let map = self.inverse_map(xf);
                out.extend(map.into_iter().zip(reads).filter(|(_, r)| *r).map(|((_, u), _)| u));
            }
        }
        out
    }

    /// Which of x, y and z field `f` reads (once per field).
    fn xyz_reads(&mut self, f: NodeId) -> [bool; 3] {
        let g = &self.g;
        *self.xyz_read.entry(f).or_insert_with(|| {
            let r = g.reachable(&[f]);
            [Var::X, Var::Y, Var::Z].map(|v| r.iter().any(|&n| matches!(g.nodes[n as usize], crate::graph::Node::Var(w) if w == v)))
        })
    }

    /// The default target: `root` if it exists, else every body.
    /// The assemblies written in the file loaded that none of its other
    /// assemblies places, in the order written: the models it puts
    /// together (what `ink render` draws without a target or a root).
    pub fn top_assemblies(&self) -> Vec<&str> {
        let own: Vec<&Assembly> = self.assemblies.iter().filter(|a| a.of.is_none() && self.is_own_file(&a.file)).collect();
        own.iter().filter(|a| !own.iter().any(|b| b.instances.iter().any(|i| i.child == a.name))).map(|a| a.name.as_str()).collect()
    }

    pub fn default_items(&mut self) -> Result<Vec<Item>, String> {
        if self.body("root").is_some() || self.assembly("root").is_some() {
            return self.items("root");
        }
        if self.bodies.is_empty() {
            return Err("document has no bodies; add (body name (shape expr))".into());
        }
        let names: Vec<String> = self.bodies.iter().filter(|b| b.of.is_none()).map(|b| b.name.clone()).collect();
        let mut out = Vec::new();
        for n in names {
            out.extend(self.items(&n)?);
        }
        Ok(out)
    }
}
