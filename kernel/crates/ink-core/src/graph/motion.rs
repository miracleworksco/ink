//! What has to be defined for a design expression to have a value. Kept
//! independently of the numerical formula, which is free to fold to zero.
//! A parameter can acquire a provider later, so keep its identity too.
use super::{Graph, Node, NodeId, Var};
use rustc_hash::FxHashMap as HashMap;
use std::sync::Arc;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct DriveRead {
    pub field: NodeId,
    pub at: Option<[NodeId; 4]>,
    /// A further spatial or parameter derivative needs values on both
    /// sides along its axis, even if the jump vanishes at its centre.
    pub near_space: [bool; 3],
    pub near_params: bool,
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum MotionRead {
    Param(u32),
    /// The original drive and where its derivative is read. No map means
    /// the current (x,y,z,t). Never differentiate the map to find whether
    /// the original drive jumps: even a frozen clock can read a switch.
    Drive(DriveRead),
}

#[derive(Default)]
pub(super) struct Domains {
    pub(super) of: Vec<u32>,
    has_drives: bool,
    sets: Vec<Arc<[MotionRead]>>,
    ids: HashMap<Arc<[MotionRead]>, u32>,
    pub(super) original: HashMap<NodeId, NodeId>,
    aliases: HashMap<(NodeId, u32), NodeId>,
}

impl Domains {
    fn set(&self, id: u32) -> &[MotionRead] {
        if id == 0 { &[] } else { &self.sets[id as usize - 1] }
    }

    pub(super) fn id(&mut self, mut reads: Vec<MotionRead>) -> u32 {
        reads.sort_unstable();
        reads.dedup();
        if reads.is_empty() {
            return 0;
        }
        if let Some(&id) = self.ids.get(reads.as_slice()) {
            return id;
        }
        let reads: Arc<[MotionRead]> = reads.into();
        self.sets.push(reads.clone());
        let id = self.sets.len() as u32;
        self.has_drives |= reads.iter().any(|r| matches!(r, MotionRead::Drive(_)));
        self.ids.insert(reads, id);
        id
    }

    pub(super) fn union(&mut self, nodes: &[NodeId]) -> u32 {
        let mut ids = nodes.iter().map(|&n| self.of[n as usize]).filter(|&id| id != 0);
        let first = ids.next().unwrap_or(0);
        if ids.all(|id| id == first) {
            return first;
        }
        self.id(nodes.iter().flat_map(|&n| self.set(self.of[n as usize]).iter().copied()).collect())
    }

    pub(super) fn node(&mut self, n: Node) -> u32 {
        let (ch, k) = n.children();
        let inherited = self.union(&ch[..k]);
        let read = match n {
            Node::Var(Var::Param(k)) => MotionRead::Param(k),
            Node::Deriv(_, f, 3) => MotionRead::Drive(DriveRead { field: f, at: None, near_space: [false; 3], near_params: false }),
            _ => return inherited,
        };
        let mut reads = self.set(inherited).to_vec();
        reads.push(read);
        self.id(reads)
    }
}

impl Graph {
    pub fn has_speed_reads(&self) -> bool {
        self.domains.has_drives
    }

    pub fn motion_reads(&self, n: NodeId) -> &[MotionRead] {
        self.domains.set(self.domains.of[n as usize])
    }
    /// Intern a separate design identity if folding would give two
    /// expressions with different requirements the same numerical node.
    /// Its tape and arithmetic are exactly those of the folded value.
    pub(super) fn with_domain(&mut self, value: NodeId, domain: u32) -> NodeId {
        if self.domains.of[value as usize] == domain {
            return value;
        }
        let value = self.domains.original.get(&value).copied().unwrap_or(value);
        if self.domains.of[value as usize] == domain {
            return value;
        }
        if let Some(&id) = self.domains.aliases.get(&(value, domain)) {
            return id;
        }
        let i = value as usize;
        let id = self.nodes.len() as NodeId;
        self.nodes.push(self.nodes[i]);
        self.shapes.push(self.shapes[i]);
        self.fixed.push(self.fixed[i]);
        self.partial.push(self.partial[i]);
        self.marked.push(self.marked[i]);
        self.domains.of.push(domain);
        self.domains.original.insert(id, value);
        self.domains.aliases.insert((value, domain), id);
        id
    }

    /// A piece extracted from a design field still requires the enclosing
    /// expression to be defined, including reads absent from the piece.
    pub fn with_reads_of(&mut self, value: NodeId, from: NodeId) -> NodeId {
        self.keeping_inputs(value, &[value, from])
    }

    pub(super) fn keeping_inputs(&mut self, value: NodeId, inputs: &[NodeId]) -> NodeId {
        let domain = self.domains.union(inputs);
        self.with_domain(value, domain)
    }

    pub(super) fn mapped_domain(&mut self, root: NodeId, value: NodeId, map: &[(Var, NodeId)]) -> NodeId {
        if self.motion_reads(root).is_empty() && map.iter().all(|(_, n)| self.motion_reads(*n).is_empty()) {
            return value;
        }
        let original = self.motion_reads(root).to_vec();
        let mut reads = Vec::new();
        let axes = [Var::X, Var::Y, Var::Z, Var::T];
        let replacements: HashMap<_, _> = map.iter().map(|&(v, n)| (self.var(v), n)).collect();
        for read in original {
            match read {
                MotionRead::Param(k) => match map.iter().find(|(v, _)| *v == Var::Param(k)) {
                    Some((_, n)) => reads.extend_from_slice(self.motion_reads(*n)),
                    None => reads.push(read),
                },
                MotionRead::Drive(d) => {
                    let f = d.field;
                    let identity = axes.map(|v| self.var(v));
                    let at = d.at.unwrap_or(identity).map(|n| self.rebuild(n, &replacements));
                    let params: HashMap<_, _> = map.iter().filter(|(v, _)| matches!(v, Var::Param(_))).map(|&(v, n)| (self.var(v), n)).collect();
                    let f = if params.is_empty() { f } else { self.rebuild(f, &params) };
                    for (axis, n) in at.into_iter().enumerate() {
                        if !self.motion_reads(n).is_empty() && self.depends_on(f, |v| v == axes[axis]) {
                            reads.extend_from_slice(self.motion_reads(n));
                        }
                    }
                    reads.push(MotionRead::Drive(DriveRead { field: f, at: (at != identity).then_some(at), ..d }));
                }
            }
        }
        // Coordinates specify where a field is read. Keep their domains
        // even if its numerical formula stopped reading that coordinate
        // (e.g. warp x-x to an undefined speed). Parameters are ordinary
        // substitutions; their original occurrences were recorded above.
        for &(var, n) in map {
            if !matches!(var, Var::Param(_)) {
                reads.extend_from_slice(self.motion_reads(n));
            }
        }
        let domain = self.domains.id(reads);
        self.with_domain(value, domain)
    }

    pub(super) fn derivative_domain(&mut self, value: NodeId, root: NodeId, var: Var) -> NodeId {
        if !self.motion_reads(root).iter().any(|r| matches!(r, MotionRead::Drive(_))) || var == Var::T {
            return self.keeping_inputs(value, &[root]);
        }
        let mut reads = self.motion_reads(root).to_vec();
        let axes = [Var::X, Var::Y, Var::Z];
        for &read in self.motion_reads(root) {
            let MotionRead::Drive(mut d) = read else { continue };
            for (k, axis) in axes.into_iter().enumerate() {
                d.near_space[k] |= d.at.map_or(axis == var, |at| self.depends_on(at[k], |v| v == var));
            }
            d.near_params |= matches!(var, Var::Param(_)) && self.depends_on(d.field, |v| v == var);
            reads.push(MotionRead::Drive(d));
        }
        let domain = self.domains.id(reads);
        self.with_domain(value, domain)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::{BinOp, UnOp};
    use std::collections::BTreeSet;

    #[test]
    fn every_constructor_preserves_its_input_reads() {
        let mut g = Graph::new();
        let (zero, one, x, t, p) = (g.constant(0.0), g.constant(1.0), g.var(Var::X), g.var(Var::T), g.var(Var::Param(0)));
        let drive = g.derivative(t, Var::T);
        let marked = g.mul(zero, drive);
        let mut covered = BTreeSet::new();
        macro_rules! check {
            ($method:ident($($arg:expr),*)) => {{
                let n = g.$method($($arg),*);
                for read in g.motion_reads(marked) {
                    assert!(g.motion_reads(n).contains(read), "{} lost {read:?}", stringify!($method));
                }
                covered.insert(stringify!($method).to_owned());
                n
            }};
        }
        assert!(matches!(g.motion_reads(marked), [MotionRead::Drive(_)]));
        for op in UnOp::ALL {
            check!(un(op, marked));
            check!(un_real(op, marked));
        }
        for op in BinOp::ALL {
            for (a, b) in [(marked, one), (one, marked)] {
                check!(bin(op, a, b));
                check!(bin_real(op, a, b));
            }
        }
        let two = g.constant(2.0);
        let retained_two = g.with_reads_of(two, marked);
        check!(frame_bin(BinOp::Mul, retained_two, two));
        check!(frame_bin(BinOp::Add, two, retained_two));
        for component in g.turn_components(retained_two) {
            assert!(g.motion_reads(component).contains(&g.motion_reads(marked)[0]));
        }
        covered.insert("turn_components".to_owned());
        // Numerical frame identities (a zero times a parameter) must still
        // fold identically when the zero has a drive to check. Otherwise a
        // rigid linkage can lose its proof that two placements share a pivot.
        let product = check!(bin_real(BinOp::Mul, marked, p));
        assert_eq!(g.constant_value(product), Some(0.0));
        assert!(g.motion_reads(product).contains(&MotionRead::Param(0)));
        check!(real(marked));
        check!(add(marked, zero));
        check!(sub(marked, marked));
        check!(mul(zero, marked));
        check!(min(marked, zero));
        check!(max(zero, marked));
        check!(neg(marked));
        check!(select(marked, one, zero));
        check!(select(one, zero, marked));
        check!(select(t, marked, one));
        check!(select(t, one, marked));
        let unused = g.select(one, marked, zero);
        assert!(g.motion_reads(unused).is_empty());
        let named = check!(tag("piece", marked));
        check!(tag_opaque("outer", named));
        check!(untagged(named));
        check!(lipschitz(one, marked));
        check!(lipschitz(marked, one));
        check!(corner_reach(marked, one, one, zero));
        check!(corner_reach(one, marked, one, zero));
        check!(corner_reach(one, one, marked, zero));
        check!(corner_reach(one, one, one, marked));
        let grid = Arc::new(crate::grid::GridData::new([0.0; 3], 1.0, [2; 3], vec![0.0; 8]));
        for at in [[marked, x, x], [x, marked, x], [x, x, marked]] {
            check!(grid(grid.clone(), at));
        }
        check!(replace(x, &[(x, marked)].into()));
        check!(replace(marked, &Default::default()));
        check!(replace_equal(x, &[(x, marked)].into()));
        check!(replace_equal(marked, &[(marked, zero)].into()));
        let peak = g.max(x, t);
        let named_peak = check!(with_reads_of(peak, marked));
        check!(with_emax(named_peak, &[(named_peak, vec![vec![x], vec![t]])].into_iter().collect()));
        check!(subst(p, &[(Var::Param(0), marked)]));
        check!(subst(zero, &[(Var::X, marked)]));
        for axis in [Var::X, Var::T, Var::Param(0)] {
            check!(diff(marked, axis));
            check!(derivative(marked, axis));
        }
        check!(mark_speed(marked, &[t]));
        check!(mark_speed(zero, &[marked]));
        check!(speed_read(marked, &[t]));
        check!(speed_read(t, &[marked]));
        check!(speed_read(zero, &[marked]));
        check!(with_reads_of(zero, marked));

        // A marker is read as a param is, and what folds it away keeps reading it, in its formula too.
        let m = g.marker(7);
        assert_eq!(g.motion_reads(m), [MotionRead::Param(7)]);
        covered.insert("marker".to_owned());
        let doubled = g.mul(two, m);
        let either = g.select(m, one, one);
        let corner = g.corner_reach(m, x, one, x);
        let cornered = g.mul(zero, corner);
        for kept in [g.mul(zero, m), g.mul(doubled, zero), g.sub(m, m), g.sub(doubled, doubled), either, cornered] {
            assert!(g.reachable(&[kept]).into_iter().any(|n| n == m), "{:?} drops the marker", g.nodes[kept as usize]);
        }

        let tenth = g.constant(0.1);
        let inexact = g.real_constant(tenth, crate::exact::Exact::of(crate::exact::Ratio::new(1, 10).unwrap()));
        for (method, leaf) in [("constant", g.constant(1.5)), ("own", g.own(2.5)), ("within", g.within(1.0, 2.0, 1.5, 0)), ("var", g.var(Var::X)), ("real_constant", inexact)] {
            assert!(g.motion_reads(leaf).is_empty(), "{method} has no field input");
            covered.insert(method.to_owned());
        }
        // Rust has no reflection over inherent methods. Discover public graph
        // constructors in the implementation, so adding one requires an example
        // in this test (including shortcuts which never call intern).
        let mut constructors = BTreeSet::new();
        let mut paths = vec![std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src")];
        while let Some(path) = paths.pop() {
            if path.is_dir() {
                paths.extend(std::fs::read_dir(path).unwrap().map(|e| e.unwrap().path()));
                continue;
            }
            if path.extension().is_none_or(|e| e != "rs") {
                continue;
            }
            let src = std::fs::read_to_string(path).unwrap();
            for item in syn::parse_file(&src).unwrap().items {
                let syn::Item::Impl(i) = item else { continue };
                if !matches!(*i.self_ty, syn::Type::Path(ref p) if p.path.segments.last().is_some_and(|s| s.ident == "Graph")) {
                    continue;
                }
                for m in i.items {
                    let syn::ImplItem::Fn(f) = m else { continue };
                    // Include every public builder even if it returns fields
                    // in an array, Option, or a new representation.
                    if matches!(f.vis, syn::Visibility::Public(_))
                        && f.sig.receiver().is_some_and(|r| matches!(r.kind, syn::ReceiverKind::Reference(_, _, Some(_))))
                        && matches!(f.sig.output, syn::ReturnType::Type(..))
                    {
                        constructors.insert(f.sig.ident.to_string());
                    }
                }
            }
        }
        assert_eq!(constructors, covered);
    }
}
