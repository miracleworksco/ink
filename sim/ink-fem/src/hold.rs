//! A fix's hold by Nitsche's method. The free material's boundary with what
//! the fix holds (the held region's edge inside the part, or the face it
//! holds) carries the held value weakly and consistently: its traction is
//! the hold's force, so the hold does not slip, and only the held material
//! is held. Per cell, with the cell's free material's own stiffness.

use ink_core::{V3, v3};

/// A point of a held boundary in a cell: which of the cell's holds, where
/// (xi), its share of the boundary (mm^2), the free material's outward
/// normal there, the value held, and what is held (the projector onto the
/// held directions: the fix's axes).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Trace {
    pub hold: usize,
    pub xi: V3,
    pub w: f64,
    /// Its share of the boundary for the traction's terms: `w`, less what
    /// the cell beside it carries (solve.rs's `Lent`).
    pub wt: f64,
    pub n: V3,
    pub v: [f64; 3],
    pub p: [[f64; 3]; 3],
}

/// One hold's terms in one cell (`ne` x `ne` row-major matrices, `ne`
/// vectors): with u the cell's dofs and w the held value at its nodes, the
/// hold's energy is u^T kh u / 2 + u^T kuw w + w^T kww w / 2 (and its held
/// value's work, gu and gw for a value not at the nodes).
pub(crate) struct Terms {
    /// Added to the cell's stiffness.
    pub kh: Vec<f64>,
    /// (read by `kuw` and `kww`: empty where they are kh's, `tie`)
    kuw: Vec<f64>,
    kww: Vec<f64>,
    /// The held value's right-hand side (-kuw g), and kww g.
    pub gu: Vec<f64>,
    pub gw: Vec<f64>,
    /// The trace's mass (gamma int N_i N_j, the held directions' alike), for
    /// the solve's start and a motion's hold were it held on every axis, and
    /// its held value's (gamma int N_i P g).
    pub mass: [f64; 64],
    pub lift: Vec<f64>,
}

impl Terms {
    pub fn zero(ndof: usize) -> Terms {
        let ne = 8 * ndof;
        Terms { kh: vec![0.0; ne * ne], kuw: vec![0.0; ne * ne], kww: vec![0.0; ne * ne], gu: vec![0.0; ne], gw: vec![0.0; ne], mass: [0.0; 64], lift: vec![0.0; ne] }
    }

    /// The held value's terms alone (`kww`, `gw`, `mass`): no stiffness on u.
    pub fn held_value_alone(self) -> Terms {
        let ne = self.kh.len().isqrt();
        let (kww, gw, mass) = (self.kww().to_vec(), self.gw, self.mass);
        Terms { kh: vec![0.0; ne * ne], kuw: vec![0.0; ne * ne], kww, gu: vec![0.0; ne], gw, mass, lift: vec![0.0; ne] }
    }

    /// kuw: -kh where tied (`tie`).
    pub fn kuw(&self) -> std::borrow::Cow<'_, [f64]> {
        if self.kuw.is_empty() { self.kh.iter().map(|v| 0.0 - v).collect() } else { std::borrow::Cow::Borrowed(&self.kuw) }
    }

    /// kww: kh where tied (`tie`).
    pub fn kww(&self) -> &[f64] {
        if self.kww.is_empty() { &self.kh } else { &self.kww }
    }

    pub fn kuw_mut(&mut self) -> &mut [f64] {
        self.untie();
        &mut self.kuw
    }

    pub fn kww_mut(&mut self) -> &mut [f64] {
        self.untie();
        &mut self.kww
    }

    /// Kept once where kuw is exactly -kh and kww kh (to the bit), as a
    /// penalty's are, which every wholly held cell's are: two of its three
    /// matrices, two thirds of what a large held region's holds take.
    pub fn tie(mut self) -> Terms {
        let tied = self.kh.iter().zip(&self.kuw).zip(&self.kww).all(|((h, u), w)| u.to_bits() == (0.0 - h).to_bits() && w.to_bits() == h.to_bits());
        if tied && !self.kuw.is_empty() {
            (self.kuw, self.kww) = (Vec::new(), Vec::new());
        }
        self
    }

    fn untie(&mut self) {
        if self.kuw.is_empty() {
            self.kuw = self.kuw().into_owned();
            self.kww = self.kh.clone();
        }
    }

    pub fn add(&mut self, o: &Terms) {
        self.untie();
        let (ouw, oww) = (o.kuw(), o.kww());
        for (a, b) in
            [(&mut self.kh, &o.kh[..]), (&mut self.kuw, &ouw[..]), (&mut self.kww, oww), (&mut self.gu, &o.gu[..]), (&mut self.gw, &o.gw[..]), (&mut self.lift, &o.lift[..])]
        {
            a.iter_mut().zip(b).for_each(|(a, b)| *a += b);
        }
        self.mass.iter_mut().zip(&o.mass).for_each(|(a, b)| *a += b);
    }

    /// A hold's terms at the constant `gamma`, from those `nitsche` gave at
    /// 1 for every piece: the constant multiplies int N_i N_j P and its
    /// held value's alone (`kww`, `lift`, `mass`), which the rest hold once.
    pub fn at(mut self, gamma: f64) -> Terms {
        self.untie();
        let more = gamma - 1.0;
        for ((h, u), w) in self.kh.iter_mut().zip(self.kuw.iter_mut()).zip(self.kww.iter_mut()) {
            (*h, *u, *w) = (*h + more * *w, *u - more * *w, gamma * *w);
        }
        for ((g, w), l) in self.gu.iter_mut().zip(self.gw.iter_mut()).zip(self.lift.iter_mut()) {
            (*g, *w, *l) = (*g + more * *l, gamma * *w, gamma * *l);
        }
        self.mass.iter_mut().for_each(|m| *m *= gamma);
        self
    }
}

/// A penalty's terms: beta |P (u - w)|^2 / 2 over what it holds, M = int
/// N_i N_j there and P its axes (`dofs`), its held value's g = int N_i v.
pub(crate) fn penalty(ndof: usize, dofs: &[bool], m: &[f64; 64], g: &[f64], beta: f64) -> Terms {
    let ne = 8 * ndof;
    let mut t = Terms::zero(ndof);
    for ij in 0..64 {
        for d in (0..ndof).filter(|&d| dofs[d]) {
            let (r, c, v) = ((ij / 8) * ndof + d, (ij % 8) * ndof + d, beta * m[ij]);
            (t.kh[r * ne + c], t.kuw[r * ne + c], t.kww[r * ne + c]) = (v, -v, v);
        }
    }
    for (i, v) in g.iter().enumerate() {
        (t.gu[i], t.gw[i], t.lift[i]) = (beta * v, beta * v, beta * v);
    }
    t.mass = m.map(|v| beta * v);
    t
}

/// The polygon (offsets from the cube's centre) where the plane
/// g[0] + g[1..4] . d = 0 crosses the cube of side `s`, kept where each
/// `clips` plane c[0] + c[1..4] . d is at most 0. A plane on a face of the
/// cube is the cube's whose `owns(axis, side)` (the face at side * s / 2
/// along axis) says so: of two cubes sharing the face, exactly one, from a
/// number both compute alike (each cube's own offset of the plane, its
/// centre's, rounds differently: a plane x = 0.05 between samples of 0.05
/// was declined by both, and its hold vanished).
pub(crate) fn plane_piece(g: &[f64; 4], s: f64, clips: &[[f64; 4]], owns: &dyn Fn(usize, f64) -> bool) -> Vec<V3> {
    let n = [g[1], g[2], g[3]];
    let l2 = v3::dot(n, n);
    if !(l2 > 0.0 && l2.is_finite() && g[0].is_finite()) {
        return Vec::new();
    }
    let (half, tie) = (0.5 * s, 1e-9 * s);
    // On a face, to rounding: the plane within `tie` of it over all of it.
    // (A plane that far from a face is on one side of it in both cubes'
    // numbers alike, so either test gives one cube.)
    let len = l2.sqrt();
    let a = (0..3).max_by(|&x, &y| n[x].abs().total_cmp(&n[y].abs())).unwrap_or(0);
    let across: f64 = (0..3).filter(|&b| b != a).map(|b| half * n[b].abs()).sum();
    if across <= tie * len && (g[0].abs() - half * n[a].abs()).abs() <= tie * len {
        let side = if g[0] * n[a] < 0.0 { 1.0 } else { -1.0 };
        if !owns(a, side) {
            return Vec::new();
        }
        let (u, v) = ((a + 1) % 3, (a + 2) % 3);
        let mut poly: Vec<V3> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .iter()
            .map(|&(x, y)| {
                let mut p = [0.0; 3];
                (p[a], p[u], p[v]) = (side * half, x * half, y * half);
                p
            })
            .collect();
        for c in clips {
            let len = (c[1] * c[1] + c[2] * c[2] + c[3] * c[3]).sqrt();
            poly = clip(&poly, c, tie * len);
        }
        return poly;
    }
    // How far the cube reaches along the normal (times its length). The
    // plane through the cube, its piece in it; on a face two cubes share,
    // the one on the free side's (g0 = +span there, -span on the other):
    // each edge traced once, by whichever cube it lies in, whatever share
    // of either the region holds (a sample nearly all held kept no piece,
    // its neighbour's face took it within a reach, and a plane tilted, or
    // nearer its face than the reach's tolerance, lost its edge).
    let span = half * (n[0].abs() + n[1].abs() + n[2].abs());
    if !(g[0] > -span && g[0] <= span) {
        return Vec::new();
    }
    let g0 = g[0];
    // A square in the plane about its point nearest the centre, wider than
    // the cube, then clipped to it.
    let p0 = v3::scale(n, -g0 / l2);
    let a = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let e1 = v3::direction(v3::cross(n, a)).unwrap_or([0.0; 3]);
    let e2 = v3::direction(v3::cross(n, e1)).unwrap_or([0.0; 3]);
    let r = 2.0 * s;
    let mut poly: Vec<V3> = [(-r, -r), (r, -r), (r, r), (-r, r)].iter().map(|&(x, y)| v3::add(p0, v3::add(v3::scale(e1, x), v3::scale(e2, y)))).collect();
    for k in 0..3 {
        for way in [1.0, -1.0] {
            let mut c = [-half, 0.0, 0.0, 0.0];
            c[1 + k] = way;
            poly = clip(&poly, &c, tie);
        }
    }
    for c in clips {
        let len = (c[1] * c[1] + c[2] * c[2] + c[3] * c[3]).sqrt();
        poly = clip(&poly, c, tie * len);
    }
    poly
}

/// The part of a polygon where c[0] + c[1..4] . p <= 0, its points within
/// `tie` past it kept where they are (a plane on the cube's face).
pub(crate) fn clip(poly: &[V3], c: &[f64; 4], tie: f64) -> Vec<V3> {
    let f = |p: V3| c[0] + c[1] * p[0] + c[2] * p[1] + c[3] * p[2];
    let mut out = Vec::with_capacity(poly.len() + 1);
    for (k, &p) in poly.iter().enumerate() {
        let q = poly[(k + 1) % poly.len()];
        let (fp, fq) = (f(p), f(q));
        if fp <= tie {
            out.push(p);
        }
        if (fp <= tie) != (fq <= tie) {
            let t = (fp / (fp - fq)).clamp(0.0, 1.0);
            out.push(std::array::from_fn(|a| p[a] + t * (q[a] - p[a])));
        }
    }
    out
}

/// The part of the cube of side `s` (about its centre) where every plane
/// (c[0] + c[1..4] . d) is at most 0: its volume, and its first and second
/// moments about the centre (the second as xx, yy, zz, xy, yz, xz),
/// exactly: its faces (the cube's, clipped to the planes; each plane's
/// piece within the others) as fans of tetrahedra from the centre. (By
/// points an eighth of the cube apart, a strip of a hold 0.01 mm square
/// along a part's edge held no material.)
pub(crate) fn within(planes: &[[f64; 4]], s: f64) -> (f64, V3, [f64; 6]) {
    let half = 0.5 * s;
    let corners = |g: &[f64; 4]| -> [f64; 8] { std::array::from_fn(|i| g[0] + (0..3).map(|a| g[1 + a] * if i >> a & 1 == 1 { half } else { -half }).sum::<f64>()) };
    let tie = |g: &[f64; 4]| 1e-9 * s * (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt();
    // A plane the whole cube lies inside says nothing; one it lies outside leaves nothing.
    if planes.iter().any(|g| corners(g).iter().all(|v| *v >= -tie(g))) {
        return (0.0, [0.0; 3], [0.0; 6]);
    }
    // (and one written twice is one: each would count the other's piece)
    let unit = |g: &[f64; 4]| g.map(|v| v / (g[1] * g[1] + g[2] * g[2] + g[3] * g[3]).sqrt());
    let mut cut: Vec<[f64; 4]> = Vec::new();
    for g in planes.iter().filter(|g| corners(g).iter().any(|v| *v > tie(g))) {
        if !cut.iter().any(|o| unit(o).iter().zip(unit(g)).enumerate().all(|(k, (a, b))| (a - b).abs() <= if k == 0 { 1e-9 * s } else { 1e-9 })) {
            cut.push(*g);
        }
    }
    let (mut vol, mut first, mut second) = (0.0, [0.0; 3], [0.0; 6]);
    let mut face = |poly: &[V3], out: V3| {
        for k in 1..poly.len().saturating_sub(1) {
            let (a, mut b, mut c) = (poly[0], poly[k], poly[k + 1]);
            if v3::dot(v3::cross(v3::sub(b, a), v3::sub(c, a)), out) < 0.0 {
                (b, c) = (c, b);
            }
            let v = v3::dot(a, v3::cross(b, c)) / 6.0;
            vol += v;
            let sum = v3::add(a, v3::add(b, c));
            (0..3).for_each(|i| first[i] += v * sum[i] / 4.0);
            for (q, (i, j)) in [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)].into_iter().enumerate() {
                second[q] += v / 20.0 * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + sum[i] * sum[j]);
            }
        }
    };
    for a in 0..3 {
        for side in [1.0, -1.0] {
            let (u, v) = ((a + 1) % 3, (a + 2) % 3);
            let mut poly: Vec<V3> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
                .iter()
                .map(|&(x, y)| {
                    let mut p = [0.0; 3];
                    (p[a], p[u], p[v]) = (side * half, x * half, y * half);
                    p
                })
                .collect();
            for g in &cut {
                poly = clip(&poly, g, 0.0);
            }
            let mut out = [0.0; 3];
            out[a] = side;
            face(&poly, out);
        }
    }
    for (i, g) in cut.iter().enumerate() {
        let others: Vec<[f64; 4]> = cut.iter().enumerate().filter(|(q, _)| *q != i).map(|(_, o)| *o).collect();
        face(&plane_piece(g, s, &others, &|_, _| true), [g[1], g[2], g[3]]);
    }
    (vol.max(0.0), first, second)
}

/// Quadrature over a flat polygon: its fan of triangles, each by its edges'
/// midpoints (exact for quadratics). Points and their weights (area).
pub(crate) fn polygon_points(poly: &[V3]) -> Vec<(V3, f64)> {
    let mut out = Vec::new();
    for k in 1..poly.len().saturating_sub(1) {
        let (a, b, c) = (poly[0], poly[k], poly[k + 1]);
        let area = 0.5 * v3::norm(v3::cross(v3::sub(b, a), v3::sub(c, a)));
        if area > 0.0 {
            let mid = |p: V3, q: V3| v3::scale(v3::add(p, q), 0.5);
            for p in [mid(a, b), mid(b, c), mid(c, a)] {
                out.push((p, area / 3.0));
            }
        }
    }
    out
}

fn signs(i: usize) -> [f64; 3] {
    [((i & 1) * 2) as f64 - 1.0, (((i >> 1) & 1) * 2) as f64 - 1.0, (((i >> 2) & 1) * 2) as f64 - 1.0]
}

fn shape_n(i: usize, xi: V3) -> f64 {
    let s = signs(i);
    (1.0 + s[0] * xi[0]) * (1.0 + s[1] * xi[1]) * (1.0 + s[2] * xi[2]) / 8.0
}

/// The terms of a cell's holds (`holds` of them, each trace naming its
/// own), each trace at its own constant (`gammas`, as `traces`), `flux` a
/// hold's traction of the free material's field per dof at a point and
/// normal (`ndof` x `ne`); and the traces' tractions' Gram matrix (sum wt
/// (P t)^T (P t), `ne` x `ne`): a cell's constant is twice the largest
/// that takes of its free material's energy (`largest`), so the system is
/// positive definite however the hold's edge cuts it (each hold its own,
/// together they could take it more than once).
pub(crate) fn nitsche(ndof: usize, traces: &[Trace], flux: &dyn Fn(usize, V3, V3) -> Vec<f64>, holds: usize, gammas: &[f64]) -> (Vec<Option<Terms>>, Vec<f64>) {
    let ne = 8 * ndof;
    let mut s = vec![0.0; ne * ne];
    let mut c = vec![vec![0.0; ne * ne]; holds];
    let mut m = vec![vec![0.0; ne * ne]; holds];
    let mut mv = vec![vec![0.0; ne]; holds];
    let mut fv = vec![vec![0.0; ne]; holds];
    let mut mass = vec![[0.0; 64]; holds];
    let mut seen = vec![false; holds];
    for (t, gamma) in traces.iter().zip(gammas) {
        let k = t.hold;
        seen[k] = true;
        let nn: [f64; 8] = std::array::from_fn(|i| shape_n(i, t.xi));
        let f = flux(k, t.xi, t.n);
        // P F and P v.
        let pf: Vec<f64> = (0..ndof * ne).map(|r| (0..ndof).map(|b| t.p[r / ne][b] * f[b * ne + r % ne]).sum()).collect();
        let pv: Vec<f64> = (0..ndof).map(|a| (0..ndof).map(|b| t.p[a][b] * t.v[b]).sum()).collect();
        for i in 0..ne {
            for j in 0..ne {
                s[i * ne + j] += t.wt * (0..ndof).map(|a| pf[a * ne + i] * pf[a * ne + j]).sum::<f64>();
            }
        }
        for i in 0..8 {
            for a in 0..ndof {
                let row = i * ndof + a;
                for col in 0..ne {
                    c[k][row * ne + col] += t.wt * nn[i] * pf[a * ne + col];
                }
                for j in 0..8 {
                    for b in 0..ndof {
                        m[k][row * ne + j * ndof + b] += gamma * t.w * nn[i] * nn[j] * t.p[a][b];
                    }
                }
                mv[k][row] += gamma * t.w * nn[i] * pv[a];
            }
            for j in 0..8 {
                mass[k][i * 8 + j] += gamma * t.w * nn[i] * nn[j];
            }
        }
        for col in 0..ne {
            fv[k][col] += t.wt * (0..ndof).map(|a| f[a * ne + col] * pv[a]).sum::<f64>();
        }
    }
    let terms = (0..holds)
        .map(|k| {
            if !seen[k] {
                return None;
            }
            let (ck, mk) = (&c[k], &m[k]);
            let kh = (0..ne * ne).map(|ij| mk[ij] - ck[ij] - ck[(ij % ne) * ne + ij / ne]).collect();
            let kuw = (0..ne * ne).map(|ij| ck[(ij % ne) * ne + ij / ne] - mk[ij]).collect();
            let gu = (0..ne).map(|i| mv[k][i] - fv[k][i]).collect();
            Some(Terms { kh, kuw, kww: mk.clone(), gu, lift: mv[k].clone(), gw: mv[k].clone(), mass: mass[k] })
        })
        .collect();
    (terms, s)
}

/// The largest lambda of S v = lambda K v over v off K's null space (the
/// cell's rigid motions, which S does not see either): of S against K plus
/// a multiple of the rigid motions' projector, which leaves the rest alone.
pub(crate) fn largest(ndof: usize, k: &[f64], s: &[f64]) -> f64 {
    let ne = 8 * ndof;
    // The rigid motions at the nodes (xi coordinates).
    let modes: Vec<Vec<f64>> = if ndof == 1 {
        vec![vec![1.0; 8]]
    } else {
        (0..6)
            .map(|r| {
                (0..24)
                    .map(|q| {
                        let (i, d) = (q / 3, q % 3);
                        if r < 3 {
                            (d == r) as u8 as f64
                        } else {
                            let x = signs(i);
                            v3::cross(std::array::from_fn(|a| (a == r - 3) as u8 as f64), x)[d]
                        }
                    })
                    .collect()
            })
            .collect()
    };
    let scale = (0..ne).map(|i| k[i * ne + i]).sum::<f64>() / ne as f64;
    let mut a: Vec<f64> = k.to_vec();
    for r in &modes {
        let rr: f64 = r.iter().map(|v| v * v).sum();
        for i in 0..ne {
            for j in 0..ne {
                a[i * ne + j] += scale * r[i] * r[j] / rr;
            }
        }
    }
    let Some(l) = cholesky(&a, ne) else { return f64::INFINITY };
    // Z = L^-1 S, then L^-1 Z^T = L^-1 S L^-T.
    let solve = |b: &[f64]| -> Vec<f64> {
        let mut x = vec![0.0; ne];
        for i in 0..ne {
            x[i] = (b[i] - (0..i).map(|j| l[i * ne + j] * x[j]).sum::<f64>()) / l[i * ne + i];
        }
        x
    };
    let mut z = vec![0.0; ne * ne];
    for col in 0..ne {
        let b: Vec<f64> = (0..ne).map(|r| s[r * ne + col]).collect();
        for (r, v) in solve(&b).into_iter().enumerate() {
            z[r * ne + col] = v;
        }
    }
    let mut w = vec![0.0; ne * ne];
    for col in 0..ne {
        let b: Vec<f64> = (0..ne).map(|r| z[col * ne + r]).collect();
        for (r, v) in solve(&b).into_iter().enumerate() {
            w[r * ne + col] = v;
        }
    }
    for i in 0..ne {
        for j in 0..i {
            let v = 0.5 * (w[i * ne + j] + w[j * ne + i]);
            w[i * ne + j] = v;
            w[j * ne + i] = v;
        }
    }
    let (vals, _) = crate::eigen::dense_sym_eigen(&w, ne);
    vals.last().copied().unwrap_or(0.0).max(0.0)
}

/// x with L L^T x = b (`cholesky`'s L).
pub(crate) fn solve(l: &[f64], n: usize, b: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        y[i] = (b[i] - (0..i).map(|j| l[i * n + j] * y[j]).sum::<f64>()) / l[i * n + i];
    }
    for i in (0..n).rev() {
        y[i] = (y[i] - (i + 1..n).map(|j| l[j * n + i] * y[j]).sum::<f64>()) / l[i * n + i];
    }
    y
}

pub(crate) fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for j in 0..n {
        let v = a[j * n + j] - (0..j).map(|q| l[j * n + q] * l[j * n + q]).sum::<f64>();
        if !(v > 0.0) {
            return None;
        }
        l[j * n + j] = v.sqrt();
        for i in j + 1..n {
            l[i * n + j] = (a[i * n + j] - (0..j).map(|q| l[i * n + q] * l[j * n + q]).sum::<f64>()) / l[j * n + j];
        }
    }
    Some(l)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_plane_through_a_cube_and_on_its_face() {
        let area = |p: &[V3]| polygon_points(p).iter().map(|q| q.1).sum::<f64>();
        // x = 0.1 across a cube of side 1: a unit square.
        // (on a face: the cube on the free side owns it, as the plane's
        // field at the face is 0)
        let free_side = |g: [f64; 4]| move |_: usize, _: f64| g[0] > 0.0;
        let piece = |g: [f64; 4], clips: &[[f64; 4]]| plane_piece(&g, 1.0, clips, &free_side(g));
        let p = piece([-0.1, 1.0, 0.0, 0.0], &[]);
        assert!((area(&p) - 1.0).abs() < 1e-12, "{}", area(&p));
        // On its face x = 0.5 (held below, free above): the cube below has
        // its free side beyond the face, so not it; the cube above, its
        // face x = -0.5, holds it. Once, either way.
        assert!(piece([-0.5, 1.0, 0.0, 0.0], &[]).len() < 3);
        let p = piece([0.5, 1.0, 0.0, 0.0], &[]);
        assert!((area(&p) - 1.0).abs() < 1e-12, "{}", area(&p));
        // Past it: nothing.
        assert!(piece([-0.6, 1.0, 0.0, 0.0], &[]).len() < 3);
        // Kept where y <= 0.25: three quarters.
        let p = piece([-0.1, 1.0, 0.0, 0.0], &[[-0.25, 0.0, 1.0, 0.0]]);
        assert!((area(&p) - 0.75).abs() < 1e-12, "{}", area(&p));
        // x + y + z = 0 through the centre: a regular hexagon of side sqrt(2)/2.
        let p = piece([0.0, 1.0, 1.0, 1.0], &[]);
        let want = 1.5 * 3f64.sqrt() * 0.5;
        assert!((area(&p) - want).abs() < 1e-12, "{} {want}", area(&p));
    }

    /// A slab of samples of side `s`, `w` wide in y and 8 deep in z, held
    /// where x < x0 + ty y + tz z: the area its samples trace, against the
    /// plane's area across the slab, sqrt(1 + ty^2 + tz^2) w 8 (each cube's
    /// centre and the plane's offset from it its own float; a face owned as
    /// the solve owns it: by the region's field at the face's centre, from
    /// its lattice point).
    fn traced_on(x0: f64, ty: f64, tz: f64, w: usize, s: f64) -> f64 {
        let field = |p: V3| p[0] - x0 - ty * p[1] - tz * p[2];
        let mut a = 0.0;
        for i in -2i64..48 {
            for j in 0..w as i64 {
                for k in 0..8i64 {
                    let at = [2 * i + 1, 2 * j + 1, 2 * k + 1];
                    let c = at.map(|v| v as f64 * (0.5 * s));
                    let g = [field(c), 1.0, -ty, -tz];
                    let owns = |axis: usize, side: f64| {
                        let mut f = at;
                        f[axis] += side as i64;
                        let held_side = g[0] < 0.0;
                        (field(f.map(|v| v as f64 * (0.5 * s))) > 0.0) == held_side
                    };
                    a += polygon_points(&plane_piece(&g, s, &[], &owns)).iter().map(|q| q.1).sum::<f64>();
                }
            }
        }
        a / ((1.0 + ty * ty + tz * tz).sqrt() * w as f64 * 8.0 * s * s)
    }

    #[test]
    fn every_hold_edge_is_traced_once_at_any_tilt_or_offset() {
        // Once a sample left under 3% free was counted all held and its edge
        // traced on the free neighbour's face only where the plane missed
        // that neighbour: a tilted plane crossing a sample plane lost up to a
        // row (0.92 of its area on 16 samples, 0.72 on 4, none on 1), and at
        // 3% an aligned plane was traced twice. Every tilt about one axis and
        // two, every offset across a sample (the old band and its edges
        // included), slabs 16, 4 and 1 samples wide: its area, once.
        let offsets = [0.0, 1e-15, 1e-12, 1e-9, 1e-6, 1e-3, 0.03, 0.03 + 1e-12, 0.25, 0.5, 0.75, 0.97, 1.0 - 1e-9, 1.0 - 1e-15];
        let tilts = [0.0, 1e-12, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0];
        for w in [16, 4, 1] {
            for &o in &offsets {
                for &t in &tilts {
                    for (ty, tz) in [(t, 0.0), (t, 0.5 * t)] {
                        let r = traced_on(20.0 + o, ty, tz, w, 1.0);
                        assert!((r - 1.0).abs() < 1e-8, "slab {w}, offset {o}, tilt ({ty}, {tz}): {r}");
                    }
                }
            }
        }
        // Samples whose side binary cannot hold, the plane on their faces
        // (each cube's own offset of it rounding apart: both declined it).
        for s in [0.05, 0.1, 0.2, 0.4] {
            for x0 in [0.05, 0.1, 0.15, 0.2, 0.25, 0.4, 0.45, 1.05, 1.2, 1.65] {
                for t in [0.0, 1e-12] {
                    let r = traced_on(x0, t, 0.0, 4, s);
                    assert!((r - 1.0).abs() < 1e-8, "samples {s}, plane {x0}, tilt {t}: {r}");
                }
            }
        }
        // A cylinder's edge (radius 5 about z through samples of side 1, its
        // tangent plane at each sample's centre): within its tangent planes'
        // chordal error of the circle's length, 2 pi 5, per unit of z.
        let (r0, mut a) = (5.3, 0.0);
        for i in -8..8 {
            for j in -8..8 {
                let c = [i as f64 + 0.5, j as f64 + 0.5];
                let rho = c[0].hypot(c[1]);
                let g = [rho - r0, c[0] / rho, c[1] / rho, 0.0];
                a += polygon_points(&plane_piece(&g, 1.0, &[], &|_, _| g[0] > 0.0)).iter().map(|q| q.1).sum::<f64>();
            }
        }
        let want = 2.0 * std::f64::consts::PI * r0;
        assert!((a / want - 1.0).abs() < 0.02, "{a} against {want}");
    }

    #[test]
    fn a_cubes_part_within_planes_is_its_parts_by_points() {
        // `within` against 48^3 points: 400 sets of one to three planes at
        // random, and planes on a face, through a corner, a parallel pair,
        // one written twice (each counted the other's piece: 1.5 of the
        // part's volume), planes meeting at the centre, and two a hair apart.
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        };
        let (s, n) = (0.7, 48usize);
        let mut sets: Vec<Vec<[f64; 4]>> = (0..400).map(|k| (0..1 + k % 3).map(|_| [0.35 * s * rnd(), rnd(), rnd(), rnd()]).collect()).collect();
        sets.extend([
            vec![[-0.5 * s, 1.0, 0.0, 0.0]],
            vec![[0.5 * s, 1.0, 0.0, 0.0]],
            vec![[-0.5 * s, 1.0, 0.0, 0.0], [0.1, 0.0, 1.0, 0.0]],
            vec![[-1.5 * 0.5 * s, 1.0, 1.0, 1.0]],
            vec![[1.5 * 0.5 * s, 1.0, 1.0, 1.0]],
            vec![[-0.1, 1.0, 0.0, 0.0], [-0.05, -1.0, 0.0, 0.0]],
            vec![[-0.1, 1.0, 0.0, 0.0], [-0.1, 1.0, 0.0, 0.0]],
            vec![[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            vec![[0.0, 1.0, 0.0, 0.0], [0.0, -1.0, 0.0, 0.0]],
            vec![[1e-12, 1.0, 0.0, 0.0], [-1e-12, -1.0, 0.0, 0.0]],
            vec![[-0.3499, 1.0, 0.0, 0.0], [0.3399, -1.0, 0.0, 0.0]],
        ]);
        let (cube, w) = (s * s * s, (s / n as f64).powi(3));
        for planes in &sets {
            let (v, m1, m2) = within(planes, s);
            let (mut pv, mut p1, mut p2) = (0.0, [0.0; 3], [0.0; 6]);
            for i in 0..n * n * n {
                let d: V3 = std::array::from_fn(|a| ((([i % n, (i / n) % n, i / (n * n)][a]) as f64 + 0.5) / n as f64 - 0.5) * s);
                if planes.iter().all(|g| g[0] + g[1] * d[0] + g[2] * d[1] + g[3] * d[2] <= 0.0) {
                    pv += w;
                    (0..3).for_each(|a| p1[a] += w * d[a]);
                    for (q, (a, b)) in [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)].into_iter().enumerate() {
                        p2[q] += w * d[a] * d[b];
                    }
                }
            }
            let e1 = (0..3).map(|a| (m1[a] - p1[a]).abs()).fold(0.0, f64::max) / (cube * s);
            let e2 = (0..6).map(|a| (m2[a] - p2[a]).abs()).fold(0.0, f64::max) / (cube * s * s);
            assert!(
                (v - pv).abs() / cube < 0.02 && e1 < 0.01 && e2 < 0.005,
                "{planes:?}: volume {v} against {pv} by points; first moments {m1:?} against {p1:?}; second {m2:?} against {p2:?}"
            );
        }
    }

    #[test]
    fn a_held_face_is_positive_definite_and_reacts_with_its_traction() {
        // Heat in a unit cell held at 1 on its face x = 1: the terms with the
        // cell's own stiffness are positive definite, and a uniform field of
        // the held value costs nothing.
        let ne = 8;
        let mut k = vec![0.0; 64];
        let g = |i: usize, xi: V3| -> V3 {
            let s = signs(i);
            let f = [1.0 + s[0] * xi[0], 1.0 + s[1] * xi[1], 1.0 + s[2] * xi[2]];
            [s[0] * f[1] * f[2] / 8.0, f[0] * s[1] * f[2] / 8.0, f[0] * f[1] * s[2] / 8.0]
        };
        let q = 1.0 / 3f64.sqrt();
        for gp in 0..8 {
            let xi = signs(gp).map(|v| v * q);
            for i in 0..8 {
                for j in 0..8 {
                    k[i * 8 + j] += v3::dot(g(i, xi), g(j, xi)) * 0.5;
                }
            }
        }
        let flux = |_: usize, xi: V3, n: V3| -> Vec<f64> { (0..8).map(|i| v3::dot(g(i, xi), n) * 2.0).collect() };
        let mut traces = Vec::new();
        for (p, w) in polygon_points(&[[1.0, -1.0, -1.0], [1.0, 1.0, -1.0], [1.0, 1.0, 1.0], [1.0, -1.0, 1.0]]) {
            traces.push(Trace { hold: 0, xi: p, w: w / 4.0, wt: w / 4.0, n: [1.0, 0.0, 0.0], v: [1.0, 0.0, 0.0], p: [[1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]] });
        }
        let gamma = 2.0 * largest(1, &k, &nitsche(1, &traces, &flux, 1, &vec![0.0; traces.len()]).1);
        let t = nitsche(1, &traces, &flux, 1, &vec![gamma; traces.len()]).0.pop().flatten().unwrap();
        let total: Vec<f64> = (0..64).map(|ij| k[ij] + t.kh[ij]).collect();
        assert!(cholesky(&total, ne).is_some());
        // T = 1 everywhere: K 1 + kh 1 = gu.
        for i in 0..8 {
            let row: f64 = (0..8).map(|j| total[i * 8 + j]).sum();
            assert!((row - t.gu[i]).abs() < 1e-9 * t.gu[i].abs().max(1.0), "{row} {}", t.gu[i]);
        }
    }

    #[test]
    fn tied_terms_read_and_add_as_the_full_ones_to_the_bit() {
        // Penalties (a wholly held cell's) keep kh alone; read back, and
        // added to, they are the three matrices bit for bit.
        let m: [f64; 64] = std::array::from_fn(|ij| 0.1 + ((ij * 37) % 11) as f64 / 7.0);
        let g: Vec<f64> = (0..24).map(|i| (i as f64 * 0.37).sin()).collect();
        let full = |beta: f64| {
            let mut t = Terms::zero(3);
            t.add(&penalty(3, &[true, false, true], &m, &g, beta));
            t.add(&penalty(3, &[true, true, true], &m, &g, 1.0 / 3.0));
            t
        };
        let (a, b) = (full(1e4), full(0.7).tie());
        assert!(b.kuw.is_empty() && b.kww.is_empty(), "a sum of penalties ties");
        let mut sum = full(1e4);
        sum.add(&full(0.7));
        let mut tied_sum = full(1e4).tie();
        tied_sum.add(&b);
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        for (x, y) in [(&sum, &tied_sum), (&a, &full(1e4).tie())] {
            assert_eq!(bits(&x.kh), bits(&y.kh));
            assert_eq!(bits(&x.kuw()), bits(&y.kuw()));
            assert_eq!(bits(x.kww()), bits(y.kww()));
        }
        // Nitsche's are not tied.
        let mut n = full(1.0);
        n.kuw_mut()[1] += 0.5;
        let n = n.tie();
        assert!(!n.kuw.is_empty() && n.kuw()[1] != 0.0 - n.kh[1]);
    }
}
