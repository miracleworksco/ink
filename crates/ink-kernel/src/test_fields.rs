//! Fields for the kernel's tests that pin results against the CLI. Each is
//! written with core ops, as in the .ink file quoted with it, so the CLI and
//! the test evaluate the same graph; a shape is clipped to a window, the
//! cube from -50 to 50, as `Doc::clip` clips to its own (the CLI's proofs
//! read the clipped field).

use ink_core::V3;
use ink_core::eval::Ctx;
use ink_core::graph::{BinOp, Graph, NodeId, UnOp, Var};
use ink_core::tape::Tape;
use std::sync::Arc;

pub const LO: V3 = [-50.0; 3];
pub const HI: V3 = [50.0; 3];

pub fn at(t: f64) -> Ctx {
    Ctx::new(t, vec![])
}

pub struct B(pub Graph);

impl B {
    pub fn new() -> B {
        B(Graph::new())
    }
    pub fn c(&mut self, v: f64) -> NodeId {
        self.0.constant(v)
    }
    pub fn v(&mut self, v: Var) -> NodeId {
        self.0.var(v)
    }
    pub fn un(&mut self, op: UnOp, a: NodeId) -> NodeId {
        self.0.un(op, a)
    }
    pub fn bin(&mut self, op: BinOp, a: NodeId, b: NodeId) -> NodeId {
        self.0.bin(op, a, b)
    }
    /// (- (abs (- a c)) h): a slab of half width h about c.
    pub fn slab(&mut self, a: NodeId, c: f64, h: f64) -> NodeId {
        let c = self.c(c);
        let d = self.bin(BinOp::Sub, a, c);
        let d = self.un(UnOp::Abs, d);
        let h = self.c(h);
        self.bin(BinOp::Sub, d, h)
    }
    /// (max (max a b) c)
    pub fn max3(&mut self, a: NodeId, b: NodeId, c: NodeId) -> NodeId {
        let ab = self.bin(BinOp::Max, a, b);
        self.bin(BinOp::Max, ab, c)
    }
    /// `Doc::clip` to the window: max(f, box), the box's six planes
    /// reduced as a balanced tree, as `Doc::box_field` does.
    pub fn clip(&mut self, f: NodeId) -> Arc<Tape> {
        let mut p = Vec::new();
        for v in [Var::X, Var::Y, Var::Z] {
            let (x, l, h) = (self.v(v), self.c(LO[0]), self.c(HI[0]));
            p.push(self.bin(BinOp::Sub, l, x));
            p.push(self.bin(BinOp::Sub, x, h));
        }
        let a = self.bin(BinOp::Max, p[1], p[2]);
        let a = self.bin(BinOp::Max, p[0], a);
        let b = self.bin(BinOp::Max, p[4], p[5]);
        let b = self.bin(BinOp::Max, p[3], b);
        let bx = self.bin(BinOp::Max, a, b);
        let f = self.bin(BinOp::Max, f, bx);
        Arc::new(Tape::new(&self.0, f))
    }
    pub fn tape(&self, f: NodeId) -> Tape {
        Tape::new(&self.0, f)
    }

    /// A box of half sizes `h` about `c`: (max (max (- (abs (- x cx)) hx) ...) ...).
    pub fn cuboid(&mut self, c: V3, h: V3) -> NodeId {
        let (x, y, z) = (self.v(Var::X), self.v(Var::Y), self.v(Var::Z));
        let (a, b, d) = (self.slab(x, c[0], h[0]), self.slab(y, c[1], h[1]), self.slab(z, c[2], h[2]));
        self.max3(a, b, d)
    }
}
