use super::{Graph, Node, NodeId, Var};
use rustc_hash::FxHashMap as HashMap;

impl Graph {
    /// The hash of `n`'s expression (`shapes`).
    pub(super) fn shape_of(&self, n: Node) -> u128 {
        let h = |a: NodeId| self.shapes[a as usize];
        let mut f = Fingerprint::default();
        match n {
            Node::Const(v) => f.word(1).f64(v),
            // Not its number: a param's number is where the model declared
            // or imported it, so an expression's form (and its tape's
            // fingerprint) is the same in any model that imports it.
            Node::Var(Var::Param(_)) => f.word(2),
            Node::Var(v) => f.word(3).word(match v {
                Var::X => 0,
                Var::Y => 1,
                Var::Z => 2,
                _ => 3,
            }),
            Node::Un(op, a) => f.word(4).word(op as u64).key(h(a)),
            Node::Bin(op, a, b) => f.word(5).word(op as u64).key(h(a)).key(h(b)),
            Node::Select(c, a, b) => f.word(6).key(h(c)).key(h(a)).key(h(b)),
            // Names are transparent: renaming a part reorders nothing.
            Node::Tag(_, a) => return h(a),
            Node::Grid(g, m, c) => {
                self.grids[g as usize].fingerprint(f.word(8).word(m as u64));
                c.iter().for_each(|n| _ = f.key(h(*n)));
                &mut f
            }
            Node::Read(g, path, a, at) => {
                self.grids[g as usize].fingerprint(f.word(11));
                self.paths[path as usize].iter().for_each(|&k| _ = f.word(k as u64));
                f.key(h(a)).key(h(at))
            }
            Node::Lipschitz(a, k, fr) => f.word(9).key(h(a)).key(h(k)).key(h(fr)),
            Node::Frame(c) => f.word(10).key(h(c[0])).key(h(c[1])).key(h(c[2])),
            Node::Deriv(a, fd, k) => f.word(12).word(k as u64).key(h(a)).key(h(fd)),
            Node::Real(a) => f.word(13).key(h(a)),
            Node::Within(l, hi, v, id) => f.word(14).word(id).key(h(l)).key(h(hi)).key(h(v)),
            Node::CornerReach(ci, r, g, at) => {
                let (a, b) = self.corners[ci as usize];
                f.word(15).key(h(a)).key(h(b)).key(h(r)).key(h(g)).key(h(at))
            }
        };
        f.value()
    }

    /// The hash of n's expression (`shapes`).
    pub(crate) fn shape(&self, n: NodeId) -> u128 {
        self.shapes[n as usize]
    }

    /// Content hash of the expressions at `roots`, with parameters replaced by
    /// their values: equal hashes mean equal functions of x, y, z, t (up to
    /// hash collisions), wherever and however they were built. Stable across
    /// runs and machines; for caching results. Construction names (tags)
    /// count only `with_names`: they matter to provenance (named faces), not
    /// to values.
    pub fn fingerprint(&self, roots: &[NodeId], params: &[f64], with_names: bool, fp: &mut Fingerprint) {
        let h = self.fingerprints(roots, params, with_names);
        for r in roots {
            fp.key(h[r]);
        }
    }

    /// `fingerprint`'s content hash of each node the `roots` reach, for a
    /// key with structure: a record writes each of its fields' hashes in
    /// the slot the field fills.
    pub fn fingerprints(&self, roots: &[NodeId], params: &[f64], with_names: bool) -> std::collections::HashMap<NodeId, u128> {
        let ids = self.reach(roots, false, true);
        let mut h: HashMap<NodeId, u128> = HashMap::with_capacity_and_hasher(ids.len(), Default::default());
        for &id in &ids {
            if let Node::Tag(_, a) = self.nodes[id as usize]
                && !with_names
            {
                h.insert(id, h[&a]); // names are transparent
                continue;
            }
            let mut f = Fingerprint::default();
            match self.nodes[id as usize] {
                Node::Const(v) => f.word(1).f64(v),
                Node::Var(Var::Param(k)) => f.word(2).f64(params.get(k as usize).copied().unwrap_or(f64::NAN)),
                Node::Var(v) => f.word(3).word(match v {
                    Var::X => 0,
                    Var::Y => 1,
                    Var::Z => 2,
                    _ => 3,
                }),
                Node::Un(op, a) => f.word(4).word(op as u64).key(h[&a]),
                Node::Bin(op, a, b) => f.word(5).word(op as u64).key(h[&a]).key(h[&b]),
                Node::Select(c, a, b) => f.word(6).key(h[&c]).key(h[&a]).key(h[&b]),
                // Whether provenance stops at the name counts too: it changes
                // which names a face reports.
                Node::Tag(t, a) => f.word(7).str(&self.tag_names[t as usize]).word(self.opaque.contains(&t) as u64).key(h[&a]),
                Node::Grid(g, m, c) => {
                    self.grids[g as usize].fingerprint(f.word(8).word(m as u64));
                    c.iter().for_each(|n| _ = f.key(h[n]));
                    &mut f
                }
                // Its value is made from these (Graph::read).
                Node::Read(g, path, _, at) => {
                    self.grids[g as usize].fingerprint(f.word(11));
                    f.word(self.paths[path as usize].len() as u64);
                    self.paths[path as usize].iter().for_each(|&k| _ = f.word(k as u64));
                    f.key(h[&at])
                }
                Node::Lipschitz(a, k, fr) => f.word(9).key(h[&a]).key(h[&k]).key(h[&fr]),
                Node::Frame(c) => f.word(10).key(h[&c[0]]).key(h[&c[1]]).key(h[&c[2]]),
                Node::Deriv(a, fd, k) => f.word(12).word(k as u64).key(h[&a]).key(h[&fd]),
                Node::Real(a) => f.word(13).key(h[&a]),
                Node::Within(l, hi, v, id) => f.word(14).word(id).key(h[&l]).key(h[&hi]).key(h[&v]),
                Node::CornerReach(ci, r, g, at) => {
                    let (a, b) = self.corners[ci as usize];
                    f.word(15).key(h[&a]).key(h[&b]).key(h[&r]).key(h[&g]).key(h[&at])
                }
            };
            h.insert(id, f.value());
        }
        h.into_iter().collect()
    }
}

/// The splitmix64 finaliser: a bijective mix of 64 bits.
pub const fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

/// A bijective mix of 128 bits, each output bit depending on every input
/// bit: spreads FNV's weak low bits.
const fn mix128(z: u128) -> u128 {
    let (lo, hi) = (z as u64, (z >> 64) as u64);
    let hi = mix64(hi ^ mix64(lo));
    let lo = mix64(lo ^ hi);
    (hi as u128) << 64 | lo as u128
}

/// A deterministic 128-bit hash (FNV-1a over bytes, then a final mix), for
/// content keys that must be stable across runs: at 128 bits, two contents
/// never share a key by chance.
#[derive(Clone, Copy, Debug)]
pub struct Fingerprint(u128);

impl Default for Fingerprint {
    fn default() -> Self {
        Fingerprint(0x6c62272e07bb014262b821756295c58d)
    }
}

impl Fingerprint {
    pub fn word(&mut self, w: u64) -> &mut Self {
        self.raw(&w.to_le_bytes())
    }
    /// Another hash (a node's, a key).
    pub fn key(&mut self, k: u128) -> &mut Self {
        self.raw(&k.to_le_bytes())
    }
    pub fn f64(&mut self, v: f64) -> &mut Self {
        self.word(v.to_bits())
    }
    pub fn str(&mut self, s: &str) -> &mut Self {
        self.word(s.len() as u64).raw(s.as_bytes())
    }
    /// Bytes as they are: the caller frames them (a length before them).
    pub fn raw(&mut self, bytes: &[u8]) -> &mut Self {
        for b in bytes {
            self.0 = (self.0 ^ *b as u128).wrapping_mul(0x0000000001000000000000000000013b);
        }
        self
    }
    pub fn value(&self) -> u128 {
        mix128(self.0)
    }
}
