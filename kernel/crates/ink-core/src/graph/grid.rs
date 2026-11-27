use super::{Graph, Node, NodeId, UnOp, Var};
use rustc_hash::FxHashMap as HashMap;

impl Graph {
    /// A sampled field read at the given coordinates (Node::Read): in value,
    /// its interpolant at the coordinates clamped into the grid,
    /// plus the distance to the grid box outside it unless the grid holds
    /// its boundary value there (`GridData::holds`) (in plain ops, so the
    /// core rules give its derivatives and ranges; the distance's squares
    /// are max(x - hi, 0)^2 + max(lo - x, 0)^2, exactly 0 inside).
    pub fn grid(&mut self, data: std::sync::Arc<crate::grid::GridData>, at: [NodeId; 3]) -> NodeId {
        let id = match self.grids.iter().position(|g| std::sync::Arc::ptr_eq(g, &data)) {
            Some(i) => i,
            None => {
                self.grids.push(data.clone());
                self.grids.len() - 1
            }
        } as u32;
        let value = self.path(Vec::new());
        self.read(id, value, at)
    }

    /// The index of a derivative path (`Graph::paths`).
    pub(super) fn path(&mut self, p: Vec<u8>) -> u32 {
        if let Some(&i) = self.path_index.get(&p) {
            return i;
        }
        self.paths.push(p.clone());
        self.path_index.insert(p, (self.paths.len() - 1) as u32);
        (self.paths.len() - 1) as u32
    }

    /// Grid g read at `at` along derivative `path` (Node::Read): its value
    /// is made from the coordinates alone, so equal reads are equal nodes.
    pub(super) fn read(&mut self, g: u32, path: u32, at: [NodeId; 3]) -> NodeId {
        let form = self.read_form(g, path);
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let map: HashMap<NodeId, NodeId> = (0..3).map(|k| (xyz[k], at[k])).collect();
        let a = self.rebuild(form, &map);
        let frame = self.intern(Node::Frame(at));
        self.intern(Node::Read(g, path, a, frame))
    }

    /// The form of grid g's read along `path` (read_form), made with every
    /// read of it (Graph::read).
    pub fn read_form_of(&self, g: u32, path: u32) -> NodeId {
        self.read_forms[&(g, path)]
    }

    /// The samples' read of grid g along `path` at (x, y, z): the value
    /// (interpolant at the clamped point, plus the distance to the grid
    /// box unless the grid holds), and each further axis its derivative
    /// along it.
    fn read_form(&mut self, g: u32, path: u32) -> NodeId {
        if let Some(&f) = self.read_forms.get(&(g, path)) {
            return f;
        }
        let xyz = [Var::X, Var::Y, Var::Z].map(|v| self.var(v));
        let f = match self.paths[path as usize].clone().split_last() {
            Some((&axis, before)) => {
                let before = self.path(before.to_vec());
                let f = self.read_form(g, before);
                self.diff(f, [Var::X, Var::Y, Var::Z][axis as usize])
            }
            None => {
                let (lo, hi) = self.grids[g as usize].corners();
                let holds = self.grids[g as usize].holds();
                let zero = self.constant(0.0);
                let mut q = xyz;
                let mut out = zero;
                for k in 0..3 {
                    let (l, h) = (self.constant(lo[k]), self.constant(hi[k]));
                    let below_hi = self.min(xyz[k], h);
                    q[k] = self.max(below_hi, l);
                    if holds {
                        continue;
                    }
                    let up = self.sub(xyz[k], h);
                    let up = self.max(up, zero);
                    let dn = self.sub(l, xyz[k]);
                    let dn = self.max(dn, zero);
                    let (up2, dn2) = (self.un(UnOp::Square, up), self.un(UnOp::Square, dn));
                    let both = self.add(up2, dn2);
                    out = self.add(out, both);
                }
                let raw = self.intern(Node::Grid(g, 0, q));
                if holds {
                    raw
                } else {
                    let dist = self.un(UnOp::Sqrt, out);
                    self.add(raw, dist)
                }
            }
        };
        self.read_forms.insert((g, path), f);
        f
    }
}
