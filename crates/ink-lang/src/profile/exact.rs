//! A flat `polygon` of straight edges, arcs and rounds, as drawn, built
//! exactly: every decision of its construction is an exact predicate on
//! its points (each float an exact rational), or an interval test on the
//! bounds of a number it constructs, and everything it emits is exact in
//! real numbers or, a circle's centre and radius, known within bounds on
//! their exact values (Graph::within), so its sign is the drawn profile's
//! wherever interval arithmetic decides it, at any magnitude, with no
//! tolerance. Drafts, and what these cannot bound, are still polygon's
//! construction (profile.rs, walls.rs), within its allowance.
//!
//! Its chord polygon is profile.rs's construction made exact: its convex
//! hull less the pockets between the two, each pocket read the same way, a
//! pocket leaving out its lid. A corner is dropped only where it lies
//! exactly on the segment between its neighbours (both edges straight); a
//! point lies on a hull line only if exactly on it (it then splits the
//! pocket there); a loop that doubles back, crosses or touches itself is
//! refused. Predicates are filtered: a float test with its error bound
//! decides most, expansions (exact sums of floats, Shewchuk's) the rest.
//!
//! A line through u and v is emitted as ((x - ux) dy - (y - uy) dx) s, its
//! differences dx = vx - ux and dy = vy - uy each the exact sum of two
//! floats (a two-sum, written as two terms so nothing folds them), s any
//! positive float near 1 / |v - u|: zero exactly on the line through u and
//! v, whatever s rounds to.
//!
//! Each arc then adds its segment (between arc and chord) to the chord
//! polygon, or takes it away, as its disc within its kite (`Kite`: the
//! chord's ends, where the end tangents meet, and a point across the
//! chord), all rational in the data. Each kite is proven clear of every
//! other edge, so the disc within it is the segment on its arc's side and
//! on the chord polygon's own side across; the chord lies inside the kite
//! and the disc, so no zero sheet runs along it.

use super::*;

/// An expansion: floats whose exact sum is the value, nonoverlapping, in
/// increasing magnitude, no zeros.
pub(super) type E = Vec<f64>;

fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    (p, a.mul_add(b, -p))
}

/// e + b, exactly (Shewchuk's grow-expansion, zeros eliminated).
fn grow(e: &[f64], b: f64) -> E {
    let mut q = b;
    let mut h = Vec::with_capacity(e.len() + 1);
    for &x in e {
        let (s, err) = two_sum(q, x);
        if err != 0.0 {
            h.push(err);
        }
        q = s;
    }
    if q != 0.0 {
        h.push(q);
    }
    h
}

/// The same value in as few terms as it takes (Shewchuk's compress).
fn compress(e: &[f64]) -> E {
    if e.len() < 2 {
        return e.to_vec();
    }
    let mut g = e.to_vec();
    let m = g.len();
    let mut q = g[m - 1];
    let mut bottom = m - 1;
    for i in (0..m - 1).rev() {
        let (s, err) = two_sum(q, g[i]);
        if err != 0.0 {
            g[bottom] = s;
            bottom -= 1;
            q = err;
        } else {
            q = s;
        }
    }
    g[bottom] = q;
    let mut h = Vec::with_capacity(m);
    let mut q = g[bottom];
    for &gi in &g[bottom + 1..] {
        let (s, err) = two_sum(gi, q);
        if err != 0.0 {
            h.push(err);
        }
        q = s;
    }
    if q != 0.0 {
        h.push(q);
    }
    h
}

pub(super) fn sum(e: &[f64], f: &[f64]) -> E {
    f.iter().fold(e.to_vec(), |acc, &b| grow(&acc, b))
}

pub(super) fn mul(e: &[f64], f: &[f64]) -> E {
    let mut out = E::new();
    for &a in e {
        for &b in f {
            let (p, err) = two_prod(a, b);
            out = grow(&grow(&out, err), p);
        }
    }
    out
}

fn neg(e: &[f64]) -> E {
    e.iter().map(|x| -x).collect()
}

/// a - b, exactly.
pub(super) fn diff(a: f64, b: f64) -> E {
    grow(&[a], -b)
}

pub(super) fn sign(e: &[f64]) -> i8 {
    match e.last() {
        None => 0,
        Some(&x) => x.signum() as i8,
    }
}

/// Whether expansions of products of these differences are exact: every
/// difference 0 or within 2^-480 to 2^480 in size, so no product or its
/// error term leaves the normal floats.
pub(super) fn safe(ds: &[E]) -> bool {
    ds.iter().flatten().all(|&v| v == 0.0 || (v.abs() >= 2f64.powi(-480) && v.abs() <= 2f64.powi(480)))
}

/// The sign of (b - a) x (c - a): positive where a, b, c turn
/// counter-clockwise. A float test with its error bound decides most (where
/// its products are far from underflow); expansions the rest, or integers
/// (`big`) where an expansion's products could leave the normal floats.
pub(super) fn orient(a: P, b: P, c: P) -> i8 {
    let l = (b[0] - a[0]) * (c[1] - a[1]);
    let r = (b[1] - a[1]) * (c[0] - a[0]);
    let det = l - r;
    // (3 + 16 eps) eps of the terms bounds the float determinant's error
    // (Shewchuk's ccwerrboundA), doubled: differences of floats far apart
    // round, and this test is only a filter.
    let bound = 7e-16 * (l.abs() + r.abs());
    if bound >= 2f64.powi(-900) {
        if det > bound {
            return 1;
        }
        if -det > bound {
            return -1;
        }
    }
    let ds = [diff(b[0], a[0]), diff(c[1], a[1]), diff(b[1], a[1]), diff(c[0], a[0])];
    if !safe(&ds) {
        return big::orient(a, b, c);
    }
    let lhs = mul(&ds[0], &ds[1]);
    let rhs = mul(&ds[2], &ds[3]);
    sign(&sum(&lhs, &neg(&rhs)))
}

/// The sign of (b - a) . (c - b): positive where a -> b -> c runs on.
pub(super) fn onward(a: P, b: P, c: P) -> i8 {
    let ds = [diff(b[0], a[0]), diff(c[0], b[0]), diff(b[1], a[1]), diff(c[1], b[1])];
    if !safe(&ds) {
        return big::onward(a, b, c);
    }
    let x = mul(&ds[0], &ds[1]);
    let y = mul(&ds[2], &ds[3]);
    sign(&sum(&x, &y))
}

/// The predicates in integers: each float an integer times a power of two,
/// all of one call's floats read at their least power, so sums and products
/// are exact at any magnitudes.
mod big {
    use super::P;
    use std::cmp::Ordering;

    /// A signed integer: its sign and its magnitude's 32-bit limbs, least
    /// first.
    #[derive(Clone)]
    pub(super) struct Z {
        neg: bool,
        mag: Vec<u32>,
    }

    fn trim(mut m: Vec<u32>) -> Vec<u32> {
        while m.last() == Some(&0) {
            m.pop();
        }
        m
    }

    fn cmp_mag(a: &[u32], b: &[u32]) -> Ordering {
        a.len().cmp(&b.len()).then_with(|| a.iter().rev().cmp(b.iter().rev()))
    }

    fn add_mag(a: &[u32], b: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(a.len().max(b.len()) + 1);
        let mut carry = 0u64;
        for i in 0..a.len().max(b.len()) {
            let s = *a.get(i).unwrap_or(&0) as u64 + *b.get(i).unwrap_or(&0) as u64 + carry;
            out.push(s as u32);
            carry = s >> 32;
        }
        out.push(carry as u32);
        trim(out)
    }

    /// a - b, |a| >= |b|.
    fn sub_mag(a: &[u32], b: &[u32]) -> Vec<u32> {
        let mut out = Vec::with_capacity(a.len());
        let mut borrow = 0i64;
        for (i, &x) in a.iter().enumerate() {
            let mut d = x as i64 - *b.get(i).unwrap_or(&0) as i64 - borrow;
            borrow = 0;
            if d < 0 {
                d += 1 << 32;
                borrow = 1;
            }
            out.push(d as u32);
        }
        trim(out)
    }

    impl Z {
        /// m 2^shift for a float's integer mantissa m.
        pub(super) fn of(m: u64, neg: bool, shift: u32) -> Z {
            let mut mag = vec![0u32; (shift / 32) as usize];
            let (q, r) = (m as u128).checked_shl(shift % 32).map(|v| (v as u64 as u32, (v >> 32) as u64)).unwrap_or((0, 0));
            mag.push(q);
            mag.push(r as u32);
            mag.push((r >> 32) as u32);
            Z { neg, mag: trim(mag) }
        }
        pub(super) fn sub(&self, o: &Z) -> Z {
            self.add(&Z { neg: !o.neg, mag: o.mag.clone() })
        }
        pub(super) fn add(&self, o: &Z) -> Z {
            if self.neg == o.neg {
                return Z { neg: self.neg, mag: add_mag(&self.mag, &o.mag) };
            }
            match cmp_mag(&self.mag, &o.mag) {
                Ordering::Less => Z { neg: o.neg, mag: sub_mag(&o.mag, &self.mag) },
                _ => Z { neg: self.neg, mag: sub_mag(&self.mag, &o.mag) },
            }
        }
        pub(super) fn mul(&self, o: &Z) -> Z {
            let mut out = vec![0u64; self.mag.len() + o.mag.len() + 1];
            for (i, &x) in self.mag.iter().enumerate() {
                let mut carry = 0u64;
                for (j, &y) in o.mag.iter().enumerate() {
                    let t = out[i + j] + x as u64 * y as u64 + carry;
                    out[i + j] = t & 0xffff_ffff;
                    carry = t >> 32;
                }
                let mut k = i + o.mag.len();
                while carry != 0 {
                    let t = out[k] + carry;
                    out[k] = t & 0xffff_ffff;
                    carry = t >> 32;
                    k += 1;
                }
            }
            Z { neg: self.neg != o.neg, mag: trim(out.into_iter().map(|v| v as u32).collect()) }
        }
        pub(super) fn sign(&self) -> i8 {
            if self.mag.is_empty() {
                0
            } else if self.neg {
                -1
            } else {
                1
            }
        }
    }

    /// A finite float as (integer mantissa, exponent, sign): v = m 2^e.
    pub(super) fn parts(v: f64) -> (u64, i32, bool) {
        let bits = v.to_bits();
        let (neg, exp, frac) = (bits >> 63 == 1, ((bits >> 52) & 0x7ff) as i32, bits & ((1 << 52) - 1));
        if exp == 0 { (frac, -1074, neg) } else { (frac | (1 << 52), exp - 1075, neg) }
    }

    /// Each coordinate as an integer at the least power of two among them.
    pub(super) fn ints<const N: usize>(vs: [f64; N]) -> [Z; N] {
        let ps = vs.map(parts);
        let least = ps.iter().filter(|p| p.0 != 0).map(|p| p.1).min().unwrap_or(0);
        ps.map(|(m, e, neg)| if m == 0 { Z { neg: false, mag: vec![] } } else { Z::of(m, neg, (e - least) as u32) })
    }

    pub(super) fn orient(a: P, b: P, c: P) -> i8 {
        let [ax, ay, bx, by, cx, cy] = ints([a[0], a[1], b[0], b[1], c[0], c[1]]);
        bx.sub(&ax).mul(&cy.sub(&ay)).sub(&by.sub(&ay).mul(&cx.sub(&ax))).sign()
    }

    pub(super) fn onward(a: P, b: P, c: P) -> i8 {
        let [ax, ay, bx, by, cx, cy] = ints([a[0], a[1], b[0], b[1], c[0], c[1]]);
        bx.sub(&ax).mul(&cx.sub(&bx)).add(&by.sub(&ay).mul(&cy.sub(&by))).sign()
    }
}

fn lex_less(a: P, b: P) -> bool {
    a[0] < b[0] || (a[0] == b[0] && a[1] < b[1])
}

/// The loop's points and its edges' bulges (edge i from point i to the
/// next), repeats dropped (the edge on from a repeated point is its last
/// copy's) and straight corners exactly on the segment between their
/// neighbours dropped, counter-clockwise; refused where its chords double
/// back, cross or touch, or it has no area.
fn points(input: &[[f64; 3]]) -> Result<(Vec<P>, Vec<f64>), String> {
    let mut pts: Vec<P> = Vec::with_capacity(input.len());
    let mut bs: Vec<f64> = Vec::with_capacity(input.len());
    for p in input {
        let q = [p[0], p[1]];
        if pts.last() == Some(&q) {
            *bs.last_mut().expect("a point") = p[2];
        } else {
            pts.push(q);
            bs.push(p[2]);
        }
    }
    while pts.len() > 1 && pts[0] == pts[pts.len() - 1] {
        pts.pop();
        bs.pop();
    }
    loop {
        let n = pts.len();
        if n < 3 && bs.iter().all(|&b| b == 0.0) || n < 2 {
            return Err(format!("a polygon needs 3 corners not in line; it has {n}"));
        }
        if n < 3 {
            return Err("a polygon of two corners and arcs is read by polygon's construction".into());
        }
        let (mut keep, mut kb): (Vec<P>, Vec<f64>) = (Vec::with_capacity(n), Vec::with_capacity(n));
        let mut dropped = false;
        for i in 0..n {
            let (a, ba) = match keep.last() {
                Some(&a) => (a, *kb.last().expect("a bulge")),
                None => (pts[n - 1], bs[n - 1]),
            };
            let (b, c) = (pts[i], pts[(i + 1) % n]);
            if ba == 0.0 && bs[i] == 0.0 && orient(a, b, c) == 0 {
                if onward(a, b, c) <= 0 {
                    return Err(format!("the polygon doubles back on itself at ({}, {})", b[0], b[1]));
                }
                dropped = true;
                continue;
            }
            keep.push(b);
            kb.push(bs[i]);
        }
        (pts, bs) = (keep, kb);
        if !dropped {
            break;
        }
    }
    let n = pts.len();
    if let Some((i, j)) = meet(&pts) {
        let (a, b) = (pts[i], pts[j]);
        return Err(format!("the polygon crosses or touches itself: its edges from ({}, {}) and from ({}, {}) meet", a[0], a[1], b[0], b[1]));
    }
    let lo = (0..n).fold(0, |m, i| if lex_less(pts[i], pts[m]) { i } else { m });
    if orient(pts[(lo + n - 1) % n], pts[lo], pts[(lo + 1) % n]) < 0 {
        // Edge k of the reversed loop is edge n - 2 - k (n - 1 for the
        // last) run backwards: its bulge negated.
        pts.reverse();
        bs = (0..n).map(|k| -bs[if k + 1 < n { n - 2 - k } else { n - 1 }]).collect();
    }
    Ok((pts, bs))
}

/// Two edges of the loop, not neighbours, that meet (cross or touch),
/// exactly: by the edges' boxes swept along x, each pair the boxes leave
/// tested by orientation.
fn meet(pts: &[P]) -> Option<(usize, usize)> {
    let n = pts.len();
    let bx = |i: usize| {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        ([a[0].min(b[0]), a[1].min(b[1])], [a[0].max(b[0]), a[1].max(b[1])])
    };
    let mut ix: Vec<usize> = (0..n).collect();
    ix.sort_by(|&i, &j| bx(i).0[0].total_cmp(&bx(j).0[0]));
    let within = |p: P, a: P, b: P| a[0].min(b[0]) <= p[0] && p[0] <= a[0].max(b[0]) && a[1].min(b[1]) <= p[1] && p[1] <= a[1].max(b[1]);
    for (k, &i) in ix.iter().enumerate() {
        let (lo, hi) = bx(i);
        for &j in &ix[k + 1..] {
            let (lo2, hi2) = bx(j);
            if lo2[0] > hi[0] {
                break;
            }
            if lo2[1] > hi[1] || hi2[1] < lo[1] || (i + 1) % n == j || (j + 1) % n == i {
                continue;
            }
            let (a, b, c, d) = (pts[i], pts[(i + 1) % n], pts[j], pts[(j + 1) % n]);
            let (o1, o2, o3, o4) = (orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b));
            let touch = (o1 == 0 && within(c, a, b)) || (o2 == 0 && within(d, a, b)) || (o3 == 0 && within(a, c, d)) || (o4 == 0 && within(b, c, d));
            if touch || (o1 * o2 < 0 && o3 * o4 < 0) {
                return Some((i.min(j), i.max(j)));
            }
        }
    }
    None
}

/// The indices of the strict convex hull of `pts` (no point on an edge),
/// in loop order.
fn hull(pts: &[P]) -> Vec<usize> {
    let mut ix: Vec<usize> = (0..pts.len()).collect();
    ix.sort_by(|&a, &b| pts[a][0].total_cmp(&pts[b][0]).then(pts[a][1].total_cmp(&pts[b][1])));
    let mut h: Vec<usize> = Vec::new();
    for pass in 0..2 {
        let start = h.len();
        let it: Box<dyn Iterator<Item = &usize>> = if pass == 0 { Box::new(ix.iter()) } else { Box::new(ix.iter().rev()) };
        for &i in it {
            while h.len() >= start + 2 && orient(pts[h[h.len() - 2]], pts[h[h.len() - 1]], pts[i]) <= 0 {
                h.pop();
            }
            h.push(i);
        }
        h.pop();
    }
    h.sort_unstable();
    h.dedup();
    h
}

/// The strict convex hull of `pts`, counter-clockwise (`apart` reads its
/// sides in order).
fn ccw_hull(pts: &[P]) -> Vec<P> {
    let mut ps = pts.to_vec();
    ps.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    ps.dedup();
    let mut h: Vec<P> = Vec::new();
    for pass in 0..2 {
        let start = h.len();
        let it: Box<dyn Iterator<Item = &P>> = if pass == 0 { Box::new(ps.iter()) } else { Box::new(ps.iter().rev()) };
        for &q in it {
            while h.len() >= start + 2 && orient(h[h.len() - 2], h[h.len() - 1], q) <= 0 {
                h.pop();
            }
            h.push(q);
        }
        h.pop();
    }
    h
}

/// The tree of a loop (counter-clockwise, simple): its hull's half-planes
/// (but its lid, the edge from its last point to its first, if `lid`),
/// less its pockets. `eid` names each loop edge's edge of the profile (none
/// for a lid). An arc's chord that is a hull edge takes the arc's disc as
/// its leaf (`Leaf::arc`) where the triangle of its chord and end tangents
/// lies within its ancestors' hull lines (`up`): bulging in, the half-plane
/// less the disc (the disc on the arc's side is its segment; any of it
/// past this piece's other hull lines lies where the parent reads it as
/// the arc's side already, across an edge of the piece); bulging out, the
/// half-plane or the disc, each of this
/// piece's lines it pokes past giving way to the segment (`Leaf::relax`:
/// the line or the segment, so the piece gains exactly the segment, as the
/// intersection of each line or the segment is the piece or the segment).
/// Arcs left out are joined to the whole profile within their kites
/// (`Doc::flat_polygon`).
/// A hull line of a piece, for its descendants: its ends, its piece's depth
/// (0 the profile's own hull, material; odd depths air), its leaf (none for
/// a lid), and whether it is a drawn edge's (no pocket under it).
#[derive(Clone)]
struct HullLine {
    u: P,
    v: P,
    depth: usize,
    leaf: Option<usize>,
    /// The pockets under it (`Build::excl`).
    pockets: Vec<usize>,
}

/// What a piece's construction shares across its recursion: the arcs'
/// kites, which arcs sit at a leaf, the leaves, how many pockets there are,
/// and the arcs' segments each pocket must leave out (`Tree::Excl`).
struct Build<'a> {
    ks: &'a [Option<(Kite, f64)>],
    placed: &'a mut [bool],
    out: &'a mut Vec<Leaf>,
    pockets: usize,
    excl: Vec<Vec<(usize, P, P, bool)>>,
}

fn piece(pts: &[P], eid: &[Option<usize>], lid: Option<usize>, depth: usize, up: &[HullLine], bd: &mut Build) -> Tree {
    let n = pts.len();
    let h = hull(pts);
    let nh = h.len();
    let is_lid = |k: usize| lid.is_some() && h[k] == n - 1 && h[(k + 1) % nh] == 0;
    let mut lines: Vec<HullLine> = Vec::with_capacity(nh);
    let mut terms = Vec::new();
    // The pockets under each hull line (the chain from one hull point to
    // the next, cut where it touches the line), each given its number.
    let mut under: Vec<Vec<(Vec<P>, Vec<Option<usize>>, usize)>> = Vec::with_capacity(nh);
    for k in 0..nh {
        let (i, j) = (h[k], h[(k + 1) % nh]);
        let (u, v) = (pts[i], pts[j]);
        let leaf = (!is_lid(k)).then(|| {
            bd.out.push(Leaf { u, v, arc: None, relax: Vec::new() });
            terms.push(Tree::Leaf(bd.out.len() - 1));
            bd.out.len() - 1
        });
        let len = (j + n - i) % n;
        let (mut seg, mut se, mut ps) = (vec![u], Vec::new(), Vec::new());
        for s in 1..=len {
            let p = pts[(i + s) % n];
            seg.push(p);
            se.push(eid[(i + s - 1) % n]);
            if s == len || orient(u, v, p) == 0 {
                if seg.len() > 2 {
                    // Run backwards; the lid last.
                    seg.reverse();
                    let mut re: Vec<Option<usize>> = se.iter().rev().copied().collect();
                    re.push(None);
                    bd.excl.push(Vec::new());
                    bd.pockets += 1;
                    ps.push((seg.clone(), re, bd.pockets - 1));
                }
                (seg, se) = (vec![p], Vec::new());
            }
        }
        lines.push(HullLine { u, v, depth, leaf, pockets: ps.iter().map(|p| p.2).collect() });
        under.push(ps);
    }
    // Each hull edge's arc: its segment is the region of depth `depth` (out
    // of the piece) or `depth + 1` (into it). Past a hull line (this
    // piece's or an ancestor's) of a piece of the same parity it must be
    // added there: the line gives way to it, and each pocket under the line
    // (which reads past its lid, its mouth) leaves it out; past one of the
    // other parity it lies where its side is read already.
    for k in (0..nh).filter(|&k| !is_lid(k)) {
        let (i, j) = (h[k], h[(k + 1) % nh]);
        let Some(e) = eid[i].filter(|_| j == (i + 1) % n) else { continue };
        let Some((kt, _)) = &bd.ks[e] else { continue };
        let (a, b) = (lines[k].u, lines[k].v);
        let out_ = orient_h(&ph(a), &ph(b), &kt.t) < 0;
        let parity = (depth + usize::from(!out_)) % 2;
        // Its triangle's corners are floats: the filtered predicate.
        let tri = [a, kt.tf, b];
        let past = |l: &HullLine| !tri.iter().all(|&p| orient(l.u, l.v, p) >= 0);
        let poked: Vec<&HullLine> = lines.iter().enumerate().filter(|&(l, _)| l != k).map(|(_, l)| l).chain(up.iter()).filter(|l| past(l)).collect();
        let give: Vec<&HullLine> = poked.into_iter().filter(|l| l.depth % 2 == parity).collect();
        // A lid (a line its own piece leaves out) cannot give way; its
        // parent reads the same line and gives way there.
        if give.iter().any(|l| l.leaf.is_none()) {
            continue;
        }
        for l in &give {
            bd.out[l.leaf.expect("a leaf")].relax.push((e, a, b, out_));
            for &q in &l.pockets {
                bd.excl[q].push((e, a, b, out_));
            }
        }
        bd.out[lines[k].leaf.expect("not a lid")].arc = Some((e, out_));
        bd.placed[e] = true;
    }
    let mut anc = up.to_vec();
    anc.extend(lines.iter().cloned());
    for ps in under {
        for (seg, re, q) in ps {
            terms.push(Tree::Neg(Box::new(piece(&seg, &re, Some(q), depth + 1, &anc, bd))));
        }
    }
    if let Some(q) = lid {
        terms.push(Tree::Excl(q));
    }
    Tree::Max(terms)
}

struct Leaf {
    u: P,
    v: P,
    /// The profile's arc on its chord, and whether its segment is outside
    /// the piece.
    arc: Option<(usize, bool)>,
    /// Arcs whose segments poke past this line where it must give way:
    /// their edge, chord (in their piece's direction) and whether the
    /// segment is right of it.
    relax: Vec<(usize, P, P, bool)>,
}

enum Tree {
    Leaf(usize),
    Max(Vec<Tree>),
    Neg(Box<Tree>),
    /// The arcs' segments pocket q leaves out (none: no term).
    Excl(usize),
}

/// A point in homogeneous coordinates (X, Y, W), each exact, W > 0:
/// (X / W, Y / W).
type H = [E; 3];

fn ph(p: P) -> H {
    let e = |v: f64| if v == 0.0 { vec![] } else { vec![v] };
    [e(p[0]), e(p[1]), vec![1.0]]
}

fn scale(e: &[f64], k: f64) -> E {
    mul(e, &[k])
}

/// The sign of the orientation of three homogeneous points, exactly.
fn orient_h(a: &H, b: &H, c: &H) -> i8 {
    let m2 = |p: &E, q: &E, r: &E, s: &E| sum(&mul(p, q), &neg(&mul(r, s)));
    let d = sum(&sum(&mul(&a[0], &m2(&b[1], &c[2], &b[2], &c[1])), &neg(&mul(&a[1], &m2(&b[0], &c[2], &b[2], &c[0])))), &mul(&a[2], &m2(&b[0], &c[1], &b[1], &c[0])));
    sign(&d) * sign(&a[2]) * sign(&b[2]) * sign(&c[2])
}

/// Whether convex polygons A and B (counter-clockwise; a segment is two
/// points) have no interior point in common: some edge's line has the
/// other wholly on its closed outer side.
fn apart(a: &[P], b: &[P]) -> bool {
    let sep = |x: &[P], y: &[P]| {
        let n = x.len();
        (0..n).any(|i| {
            let (p, q) = (&x[i], &x[(i + 1) % n]);
            y.iter().all(|&r| orient(*p, *q, r) <= 0) || (n == 2 && y.iter().all(|&r| orient(*q, *p, r) <= 0))
        })
    };
    sep(a, b) || sep(b, a)
}

/// An arc's kite: its chord's ends a and q, t where its end tangents meet
/// (beyond the arc), and t' = m - k (t - m) across the chord from it (m
/// the chord's middle), counter-clockwise. The arc's segment (between arc
/// and chord) is the disc within the triangle a q t; the kite holds the
/// chord inside it.
pub(super) struct Kite {
    /// Where the end tangents meet (exact), a float corner beyond it, the
    /// chord's middle (a float near it).
    t: H,
    pub(super) tf: P,
    pub(super) m: P,
    out: bool,
}

impl Kite {
    /// The arc from a to q of bulge b: t is the chord's middle plus its
    /// right normal (ey, -ex) times b / (1 - b^2) (the sweep is 4 atan b, each
    /// end tangent turned half of it from the chord); none from half a turn.
    /// Its float corner tf is beyond t from the chord, the triangle a q tf
    /// proven to hold t, so the arc's segment, exactly.
    pub(super) fn of(a: P, q: P, b: f64) -> Option<Kite> {
        if b.abs() >= 1.0 {
            return None;
        }
        let (ex, ey) = (diff(q[0], a[0]), diff(q[1], a[1]));
        let w = grow(&neg(&mul(&[b], &[b])), 1.0);
        let mx = scale(&grow(&[a[0]], q[0]), 0.5);
        let my = scale(&grow(&[a[1]], q[1]), 0.5);
        let t = [compress(&sum(&mul(&mx, &w), &scale(&ey, b))), compress(&sum(&mul(&my, &w), &neg(&scale(&ex, b)))), compress(&w)];
        let m = [0.5 * (a[0] + q[0]), 0.5 * (a[1] + q[1])];
        let out = b > 0.0;
        // Along the chord's right normal (left, bulging in), just past the
        // apex's height and on, until the triangle holds t: at any size,
        // however far below a float's step at its coordinates the arc's sag.
        let (fx, fy) = (q[0] - a[0], q[1] - a[1]);
        let l = fx.hypot(fy);
        let sg = if out { 1.0 } else { -1.0 };
        let nrm = [sg * fy / l, -sg * fx / l];
        let step = f64::EPSILON * a[0].abs().max(a[1].abs()).max(q[0].abs()).max(q[1].abs());
        let apex = l * b.abs() / (1.0 - b * b);
        let mut more = step.max(apex * 1e-12);
        for _ in 0..40 {
            let h = apex + more;
            let tf = [m[0] + nrm[0] * h, m[1] + nrm[1] * h];
            let tri = if out { [ph(a), ph(tf), ph(q)] } else { [ph(a), ph(q), ph(tf)] };
            if (0..3).all(|i| orient_h(&tri[i], &tri[(i + 1) % 3], &t) >= 0) {
                return Some(Kite { t, tf, m, out });
            }
            more *= 4.0;
        }
        None
    }

    /// The float corner across the chord for k: m - k (tf - m).
    pub(super) fn across(&self, k: f64) -> P {
        [self.m[0] - k * (self.tf[0] - self.m[0]), self.m[1] - k * (self.tf[1] - self.m[1])]
    }

    /// The triangle a q tf, counter-clockwise (bulging out of the loop,
    /// right of a -> q, tf is right of it).
    fn tri(&self, a: P, q: P) -> Vec<P> {
        if self.out { vec![a, self.tf, q] } else { vec![a, q, self.tf] }
    }

    /// The triangle across: a q t'; none where t' is not strictly across.
    fn below(&self, a: P, q: P, k: f64) -> Option<Vec<P>> {
        let t2 = self.across(k);
        let side = orient(a, q, t2);
        if (self.out && side <= 0) || (!self.out && side >= 0) {
            return None;
        }
        Some(if self.out { vec![a, q, t2] } else { vec![a, t2, q] })
    }

    fn kite(&self, a: P, q: P, k: f64) -> Vec<P> {
        let t2 = self.across(k);
        if self.out { vec![a, self.tf, q, t2] } else { vec![a, t2, q, self.tf] }
    }
}

/// Each arc's kite, and its k: the triangle a q t clear of every other
/// edge (a straight one, or another arc's kite), and the triangle across,
/// k halved until it is too. Then no other edge enters an arc's segment or
/// the kite's part across its chord, the profile is simple, and the kite
/// across the chord lies wholly on the chord polygon's side there. Boxes
/// first, then exactly.
fn kites(pts: &[P], bs: &[f64]) -> Result<Vec<Option<(Kite, f64)>>, String> {
    let n = pts.len();
    let mut ks: Vec<Option<(Kite, f64)>> = Vec::with_capacity(n);
    let mut boxes: Vec<(P, P)> = Vec::with_capacity(n);
    for i in 0..n {
        let (a, q, b) = (pts[i], pts[(i + 1) % n], bs[i]);
        let (mut lo, mut hi) = ([a[0].min(q[0]), a[1].min(q[1])], [a[0].max(q[0]), a[1].max(q[1])]);
        if b == 0.0 {
            ks.push(None);
        } else {
            if !disc_in_range(a, q, b) {
                return Err("an arc's exact circle needs a number outside ink's range".into());
            }
            let kt = Kite::of(a, q, b).ok_or("an arc of half a turn or more")?;
            // tf, and t' for k = 1 (k only shrinks it toward the chord).
            for p in [kt.tf, kt.across(1.0)] {
                for c in 0..2 {
                    lo[c] = lo[c].min(p[c]);
                    hi[c] = hi[c].max(p[c]);
                }
            }
            // A kite's corner across, halved toward m (in the box of a and
            // q), rounds within a few ulps of the largest coordinate of the
            // segment from m to t': its box grown by 4 of them.
            let big = [a, q, kt.tf, kt.across(1.0)].iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()));
            let e = 4.0 * f64::EPSILON * big;
            for c in 0..2 {
                lo[c] -= e;
                hi[c] += e;
            }
            let fits = |p: P| p.iter().all(|&c| c == 0.0 || (1e-100..=1e100).contains(&c.abs()));
            if !fits(kt.tf) || !fits(kt.across(1.0)) || !fits(kt.m) {
                return Err("an arc's kite needs a number outside ink's range".into());
            }
            ks.push(Some((kt, 1.0)));
        }
        boxes.push((lo, hi));
    }
    let near = |i: usize, j: usize| {
        let ((l1, h1), (l2, h2)) = (boxes[i], boxes[j]);
        !(l2[0] > h1[0] || l1[0] > h2[0] || l2[1] > h1[1] || l1[1] > h2[1])
    };
    let shape = |ks: &[Option<(Kite, f64)>], j: usize| match &ks[j] {
        None => vec![pts[j], pts[(j + 1) % n]],
        Some((kt, k)) => kt.kite(pts[j], pts[(j + 1) % n], *k),
    };
    for i in (0..n).filter(|&i| bs[i] != 0.0) {
        let (a, q) = (pts[i], pts[(i + 1) % n]);
        let tri = ks[i].as_ref().expect("an arc").0.tri(a, q);
        for j in (0..n).filter(|&j| j != i && near(i, j)) {
            let other = match &ks[j] {
                None => vec![pts[j], pts[(j + 1) % n]],
                Some((kt, _)) => kt.tri(pts[j], pts[(j + 1) % n]),
            };
            if !apart(&tri, &other) {
                return Err(format!("the arc from ({}, {}) comes near another edge", a[0], a[1]));
            }
        }
    }
    for _ in 0..64 {
        let mut halved = false;
        for i in (0..n).filter(|&i| bs[i] != 0.0) {
            let (a, q) = (pts[i], pts[(i + 1) % n]);
            for j in (0..n).filter(|&j| j != i && near(i, j)) {
                let (kt, k) = ks[i].as_ref().expect("an arc");
                let Some(below) = kt.below(a, q, *k) else {
                    return Err("an arc's kite across its chord rounds onto the chord".into());
                };
                if !apart(&below, &shape(&ks, j)) {
                    let k = ks[i].as_mut().expect("an arc");
                    k.1 *= 0.5;
                    halved = true;
                }
            }
        }
        if !halved {
            return Ok(ks);
        }
    }
    Err("an arc's chord meets its neighbours at a vanishing angle".into())
}

/// An arc's circle in exact parts: W = 4 |b|, W times its centre (C), and
/// G = W u - C (u its start), each an expansion.
struct DiscParts {
    w: f64,
    cx: E,
    cy: E,
    gx: E,
    gy: E,
}

fn disc_parts(u: P, v: P, b: f64) -> DiscParts {
    let (ex, ey) = (diff(v[0], u[0]), diff(v[1], u[1]));
    let b2m1 = grow(&mul(&[b], &[b]), -1.0);
    let sg = b.signum();
    let cx = compress(&scale(&sum(&scale(&grow(&[u[0]], v[0]), 2.0 * b), &mul(&ey, &b2m1)), sg));
    let cy = compress(&scale(&sum(&scale(&grow(&[u[1]], v[1]), 2.0 * b), &neg(&mul(&ex, &b2m1))), sg));
    let w = 4.0 * b.abs();
    let gx = compress(&sum(&mul(&[w], &[u[0]]), &neg(&cx)));
    let gy = compress(&sum(&mul(&[w], &[u[1]]), &neg(&cy)));
    // The field is the same for any W, C, G scaled alike: by a power of two
    // (exact) putting |G| (W times the radius) near 1, so its numbers stay
    // in range however small or large the arc.
    let big = gx.iter().chain(gy.iter()).fold(0.0f64, |m, v| m.max(v.abs()));
    let k = if big > 0.0 { (-big.log2().round()).clamp(-900.0, 900.0) as i32 } else { 0 };
    let s = 2f64.powi(k);
    let sc = |e: &E| e.iter().map(|v| v * s).collect::<E>();
    DiscParts { w: w * s, cx: sc(&cx), cy: sc(&cy), gx: sc(&gx), gy: sc(&gy) }
}

/// A number's identity (Graph::within) from what defines it: what it is
/// (`what`), the exact inputs of its construction, and which of its
/// numbers it is (`role`).
fn ident(what: &str, inputs: &[f64], role: u64) -> u64 {
    let mut f = ink_core::graph::Fingerprint::default();
    f.str(what).word(inputs.len() as u64);
    inputs.iter().for_each(|&v| _ = f.word(v.to_bits()));
    f.word(role);
    f.value() as u64
}

/// The circle of the arc from u to v of bulge b: its centre less u (its
/// size the circle's, however far out the arc lies) and its radius, each
/// as bounds on its exact value (-G / W and |G| / W, bounded outward) and a
/// float within them.
pub(super) fn circle(u: P, v: P, b: f64) -> [[f64; 3]; 3] {
    use ink_core::interval::{self, I, add_o};
    let DiscParts { w, gx, gy, .. } = disc_parts(u, v, b);
    let total = |e: &E| e.iter().fold(I::point(0.0), |s, &t| add_o(s, I::point(t)));
    let over_w = |x: I| interval::bin(BinOp::Div, x, I::point(w)).0;
    let q = compress(&sum(&mul(&gx, &gx), &mul(&gy, &gy)));
    let r = over_w(interval::un(UnOp::Sqrt, total(&q)));
    let off = |g: &E| over_w(total(&neg(g)));
    [off(&gx), off(&gy), r].map(|i| [i.lo, i.hi, i.mid().clamp(i.lo, i.hi)])
}

/// Whether every number an arc's exact disc writes is 0 or within ink's
/// range (1e-100 to 1e100 in size; past it every claim on the model is
/// unknown): an expansion's low terms can fall below it. Else the profile
/// is polygon's construction.
fn disc_in_range(u: P, v: P, b: f64) -> bool {
    let DiscParts { w, cx, cy, gx, gy } = disc_parts(u, v, b);
    let mut all: Vec<E> = Vec::new();
    if b.abs() >= 1.0 / 1024.0 {
        all.extend(circle(u, v, b).map(|c| c.to_vec()));
    } else {
        all.extend([vec![w], gx.clone(), gy.clone()]);
        let fl = super::arc_of(u, v, b);
        all.extend([mul(&[w], &[w]), compress(&scale(&gx, 2.0 * w)), compress(&scale(&gy, 2.0 * w)), cx, cy, vec![fl.c[0], fl.c[1], fl.r]]);
    }
    all.iter().flatten().all(|&c| c == 0.0 || (1e-100..=1e100).contains(&c.abs()))
}

impl Doc {
    /// A flat profile (x, y, bulge at each point), built exactly: of
    /// straight edges always; with drawn arcs where each arc's tangent
    /// triangle is clear of the other edges and within its pieces' hull
    /// lines (else, and with rounds' solved arcs, by polygon's construction,
    /// within its allowance).
    pub(crate) fn flat_polygon(&mut self, input: &[[f64; 3]], drawn: bool) -> Result<NodeId, String> {
        let arcs = input.iter().any(|p| p[2] != 0.0);
        if (drawn || !arcs)
            && let Some(f) = self.exact_flat(input)?
        {
            return Ok(f);
        }
        let t = polygon(input)?;
        self.lower_profile(&t, None)
    }

    /// `flat_polygon` built exactly; none where it is polygon's
    /// construction (an arc it cannot build exactly).
    fn exact_flat(&mut self, input: &[[f64; 3]]) -> Result<Option<NodeId>, String> {
        let arcs = input.iter().any(|p| p[2] != 0.0);
        // Far from 1 in size (past 2^+-400), built at a power-of-two scale
        // (exact) putting its largest coordinate near 1, so no product its
        // lines emit leaves the normal floats; read at that scale. (Its
        // predicates are exact at any size: `big`.) An arc whose exact
        // circle needs a number outside ink's range there: built again at
        // the scale putting the profile's size near 1.
        let big = input.iter().fold(0.0f64, |m, p| m.max(p[0].abs()).max(p[1].abs()));
        let far = big > 0.0 && !(2f64.powi(-400)..=2f64.powi(400)).contains(&big);
        let k = if far { (-big.log2().ceil()).clamp(-1000.0, 1000.0) as i32 } else { 0 };
        let size = [0, 1].map(|c| {
            let (lo, hi) = input.iter().fold((f64::MAX, f64::MIN), |(l, h), p| (l.min(p[c]), h.max(p[c])));
            hi - lo
        });
        let size = size[0].max(size[1]);
        // (No larger than puts its largest coordinate near 2^400.)
        let fit = (arcs && size > 0.0 && size.is_finite()).then(|| (-size.log2().ceil()).min(400.0 - big.log2().ceil()).clamp(-1000.0, 1000.0) as i32).filter(|&fit| fit != k);
        let build = |sc: f64| {
            let scaled: Vec<[f64; 3]> = input.iter().map(|p| [p[0] * sc, p[1] * sc, p[2]]).collect();
            points(&scaled).and_then(|(pts, bs)| {
                let ks = if arcs { kites(&pts, &bs)? } else { (0..pts.len()).map(|_| None).collect() };
                let mut leaves = Vec::new();
                let mut placed = vec![false; pts.len()];
                let eid: Vec<Option<usize>> = (0..pts.len()).map(Some).collect();
                let mut bd = Build { ks: &ks, placed: &mut placed, out: &mut leaves, pockets: 0, excl: Vec::new() };
                let t = piece(&pts, &eid, None, 0, &[], &mut bd);
                let excl = bd.excl;
                Ok((pts, bs, ks, t, leaves, placed, excl))
            })
        };
        let mut sc = 2f64.powi(k);
        let mut built = build(sc);
        if let (Err(why), Some(fit)) = (&built, fit)
            && why.contains("outside ink's range")
        {
            sc = 2f64.powi(fit);
            built = build(sc);
        }
        let (pts, bs, ks, t, leaves, placed, excl) = match built {
            Ok(b) => b,
            Err(_) if arcs => return Ok(None),
            Err(why) => return Err(why),
        };
        let (x0, y0) = (self.g.var(Var::X), self.g.var(Var::Y));
        let c = self.g.constant(sc);
        let (x, y) = (self.g.mul(c, x0), self.g.mul(c, y0));
        let n = pts.len();
        let mut discs: Vec<Option<NodeId>> = vec![None; n];
        for (i, kt) in ks.iter().enumerate() {
            if kt.is_some() && placed[i] {
                discs[i] = Some(self.exact_disc(pts[i], pts[(i + 1) % n], bs[i], [x, y]));
            }
        }
        let mut f = self.lower_exact(&t, &leaves, &discs, &excl, [x, y]);
        // Each arc not taken by its chord's leaf: its segment added to the
        // chord polygon (bulging out) or taken from it (in), as its disc
        // within its kite. The kite holds the chord inside it, where the
        // disc is negative too: no zero sheet along the chord.
        let (mut adds, mut takes) = (Vec::new(), Vec::new());
        for (i, kt) in ks.iter().enumerate() {
            let Some((kt, k)) = kt else { continue };
            if placed[i] {
                continue;
            }
            let piece = self.kite_disc(pts[i], pts[(i + 1) % n], bs[i], kt, *k, [x, y]);
            if kt.out {
                adds.push(piece);
            } else {
                takes.push(self.g.neg(piece));
            }
        }
        // The kites are apart, so the order of these does not matter: each
        // set as one balanced min or max.
        if !adds.is_empty() {
            adds.push(f);
            f = self.reduce(BinOp::Min, &adds);
        }
        if !takes.is_empty() {
            takes.push(f);
            f = self.reduce(BinOp::Max, &takes);
        }
        // Its box: its corners', and its arcs' tangent triangles'.
        let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
        for i in 0..n {
            let (mut l, mut h) = (pts[i], pts[i]);
            if let Some((kt, _)) = &ks[i] {
                (l, h) = ([kt.tf[0].min(kt.m[0]), kt.tf[1].min(kt.m[1])], [kt.tf[0].max(kt.m[0]), kt.tf[1].max(kt.m[1])]);
            }
            for k in 0..2 {
                lo[k] = lo[k].min(pts[i][k]).min(l[k]);
                hi[k] = hi[k].max(pts[i][k]).max(h[k]);
            }
        }
        let mut terms = vec![f];
        for (v, k) in [(x, 0), (y, 1)] {
            terms.push(self.axis_face(v, hi[k], true));
            terms.push(self.axis_face(v, lo[k], false));
        }
        let f = self.reduce(BinOp::Max, &terms);
        let back = self.g.constant(1.0 / sc);
        Ok(Some(self.g.mul(back, f)))
    }

    fn lower_exact(&mut self, t: &Tree, leaves: &[Leaf], kd: &[Option<NodeId>], excl: &[Vec<(usize, P, P, bool)>], xy: [NodeId; 2]) -> NodeId {
        match t {
            Tree::Leaf(k) => {
                let Leaf { u, v, arc, ref relax } = leaves[*k];
                let mut h = self.exact_line(u, v, xy);
                if let Some((e, out)) = arc {
                    let d = kd[e].expect("an arc's disc");
                    h = if out {
                        self.g.min(h, d)
                    } else {
                        let nd = self.g.neg(d);
                        self.g.max(h, nd)
                    };
                }
                // The segment of each arc poking past: its disc beyond its
                // chord.
                for &(e, a, b, right) in relax {
                    let d = kd[e].expect("an arc's disc");
                    let c = self.exact_line(a, b, xy);
                    let side = if right { self.g.neg(c) } else { c };
                    let seg = self.g.max(d, side);
                    h = self.g.min(h, seg);
                }
                h
            }
            Tree::Neg(a) => {
                let a = self.lower_exact(a, leaves, kd, excl, xy);
                self.g.neg(a)
            }
            // Each segment left out: its field negated (positive on it).
            Tree::Excl(q) => {
                let mut ids = Vec::new();
                for &(e, a, b, right) in &excl[*q] {
                    let d = kd[e].expect("an arc's disc");
                    let c = self.exact_line(a, b, xy);
                    let side = if right { self.g.neg(c) } else { c };
                    let seg = self.g.max(d, side);
                    ids.push(self.g.neg(seg));
                }
                self.reduce(BinOp::Max, &ids)
            }
            Tree::Max(v) => {
                let mut ids = Vec::with_capacity(v.len());
                for a in v {
                    if matches!(a, Tree::Excl(q) if excl[*q].is_empty()) {
                        continue;
                    }
                    let id = self.lower_exact(a, leaves, kd, excl, xy);
                    if !ids.contains(&id) {
                        ids.push(id);
                    }
                }
                self.reduce(BinOp::Max, &ids)
            }
        }
    }

    /// A line square to an axis, at coordinate c of that axis (v): v - c,
    /// positive past c (`above`), or its negation; exact, at most
    /// 1-Lipschitz, one value for every edge, box side or slab written on
    /// that line.
    pub(super) fn axis_face(&mut self, v: NodeId, c: f64, above: bool) -> NodeId {
        let at = self.g.constant(c);
        let f = self.g.sub(v, at);
        if above { f } else { self.g.neg(f) }
    }

    /// An arc's disc within its kite: negative on its segment, on its kite's
    /// part across its chord, and on the chord between.
    fn kite_disc(&mut self, a: P, q: P, b: f64, kt: &Kite, k: f64, [x, y]: [NodeId; 2]) -> NodeId {
        let d = self.exact_disc(a, q, b, [x, y]);
        // Its four sides, counter-clockwise, each exact through its float
        // corners.
        let (t, t2) = (kt.tf, kt.across(k));
        let poly = if kt.out { [a, t, q, t2] } else { [a, t2, q, t] };
        let [l1, l2, l3, l4] = [0, 1, 2, 3].map(|e| self.exact_line(poly[e], poly[(e + 1) % 4], [x, y]));
        self.reduce(BinOp::Max, &[d, l1, l2, l3, l4])
    }

    /// The disc of the arc from u to v of bulge b, its drawn circle's.
    /// Its centre c = mid + (ey, -ex) (b^2 - 1) / 4b, scaled by W = 4 |b|,
    /// is W c = C, a polynomial in the data (an expansion, its terms written
    /// apart), and G = W u - C = W (u - c), so c - u and R are -G / W and
    /// |G| / W, each bounded outward. A flat arc's centre is far out, so its
    /// disc is read about u in real numbers: u is on the circle, so W^2
    /// (|p - c|^2 - R^2) = W^2 |p - u|^2 + 2 W (p - u) . G (no term of the
    /// circle's size cancels), over (|W p - C| + R') W, R' a float at least
    /// W R: zero exactly on the circle, |p - c| - R near it, its slope at
    /// most 1 (1 / W rounded down).
    fn exact_disc(&mut self, u: P, v: P, b: f64, [x, y]: [NodeId; 2]) -> NodeId {
        let DiscParts { w, cx, cy, gx, gy } = disc_parts(u, v, b);
        // An arc of radius up to 1024 times its chord: |(p - u) - (c - u)|
        // - R, c - u and R each a number known within its exact value's
        // bounds (Graph::within), so ranges read the drawn circle, a point
        // within those bounds of it is undecided, and nothing of the
        // coordinates' size rounds (p - u is exact near u).
        if b.abs() >= 1.0 / 1024.0 {
            let mut role = 0;
            let [cx, cy, r] = circle(u, v, b).map(|[lo, hi, m]| {
                role += 1;
                self.g.within(lo, hi, m, ident("arc circle", &[u[0], u[1], v[0], v[1], b], role))
            });
            let g = &mut self.g;
            let (ux, uy) = (g.constant(u[0]), g.constant(u[1]));
            let (px, py) = (g.sub(x, ux), g.sub(y, uy));
            let (dx, dy) = (g.sub(px, cx), g.sub(py, cy));
            let (sx, sy) = (g.un(UnOp::Square, dx), g.un(UnOp::Square, dy));
            let q = g.add(sx, sy);
            let d = g.un(UnOp::Sqrt, q);
            let f = g.sub(d, r);
            let one = g.constant(1.0);
            return g.lipschitz(f, one);
        }
        // Flatter, about its end u (no term of the circle's size cancels),
        // held to the float circle below.
        let g = &mut self.g;
        // Terms k_i * a, each written apart.
        let terms = |g: &mut Graph, a: NodeId, k: &[f64], acc: &mut Option<NodeId>| {
            for &c in k.iter().rev() {
                let c = g.constant(c);
                let t = g.mul(a, c);
                *acc = Some(match *acc {
                    None => t,
                    Some(s) => g.add(s, t),
                });
            }
        };
        let (ux, uy) = (g.constant(u[0]), g.constant(u[1]));
        let (px, py) = (g.sub(x, ux), g.sub(y, uy));
        let (sx, sy) = (g.un(UnOp::Square, px), g.un(UnOp::Square, py));
        let s2 = g.add(sx, sy);
        let mut num = None;
        terms(g, s2, &mul(&[w], &[w]), &mut num);
        terms(g, px, &compress(&scale(&gx, 2.0 * w)), &mut num);
        terms(g, py, &compress(&scale(&gy, 2.0 * w)), &mut num);
        let num = num.expect("terms");
        // |W p - C|, and R' >= W R = |W u - C| = |G|.
        let wc = g.constant(w);
        let mut axis = |var: NodeId, c: &E| {
            let mut d = g.mul(wc, var);
            for &t in c.iter().rev() {
                let t = g.constant(t);
                d = g.sub(d, t);
            }
            g.un(UnOp::Square, d)
        };
        let (dx2, dy2) = (axis(x, &cx), axis(y, &cy));
        let d2 = g.add(dx2, dy2);
        let dist = g.un(UnOp::Sqrt, d2);
        // |G|^2, bounded outward (its expansion's sum rounds).
        let gg = sum(&mul(&gx, &gx), &mul(&gy, &gy)).iter().fold(ink_core::interval::I::point(0.0), |a, &t| ink_core::interval::add_o(a, ink_core::interval::I::point(t)));
        let r_up = ink_core::interval::un(UnOp::Sqrt, gg).hi;
        let r = g.constant(r_up);
        let den = g.add(dist, r);
        let f = g.bin(BinOp::Div, num, den);
        let inv = g.constant(ink_core::interval::inv_norm_down(&[w]));
        let exact = g.mul(inv, f);
        // Held to the float circle's field |p - cf| - rf within delta, a
        // bound on how far the two are apart within the circle's size of it
        // (the centres', the radii's, the denominator's radius and the 1 / W
        // scale's): the identity there in real numbers, and past it (where
        // both are far from 0, of one sign) a value of that sign; but its
        // bounds over a box are the float circle's, as tight as a distance
        // to a point is (the exact form reads x in a square and a product
        // apart, and its bounds over a box are loose).
        use ink_core::interval::{add_up, mul_up, sqrt_up, sum_up};
        let fl = super::arc_of(u, v, b);
        let up_size = |e: &E| sum_up(e.iter().map(|v| v.abs()));
        let dcx = sum(&cx, &neg(&mul(&[w], &[fl.c[0]])));
        let dcy = sum(&cy, &neg(&mul(&[w], &[fl.c[1]])));
        let off_c = (sqrt_up(add_up(mul_up(up_size(&dcx), up_size(&dcx)), mul_up(up_size(&dcy), up_size(&dcy)))) / w).next_up();
        let over_w = |v: f64| ink_core::interval::bin(BinOp::Div, ink_core::interval::I::point(v), ink_core::interval::I::point(w)).0;
        let r_lo = over_w(ink_core::interval::un(UnOp::Sqrt, ink_core::interval::I::new(gg.lo.max(0.0), gg.hi)).lo).lo;
        let r_hi = over_w(r_up).hi;
        let off_r = (r_hi - fl.r).abs().max((fl.r - r_lo).abs()).next_up();
        let delta = sum_up([off_c, off_r, (r_hi - r_lo).next_up(), mul_up(r_hi, 8.0 * f64::EPSILON)].into_iter());
        let (fcx, fcy, fr) = (g.constant(fl.c[0]), g.constant(fl.c[1]), g.constant(fl.r));
        let (ddx, ddy) = (g.sub(x, fcx), g.sub(y, fcy));
        let (qx, qy) = (g.un(UnOp::Square, ddx), g.un(UnOp::Square, ddy));
        let q = g.add(qx, qy);
        let fd = g.un(UnOp::Sqrt, q);
        let float = g.sub(fd, fr);
        let dl = g.constant(delta);
        let (lo, hi) = (g.sub(float, dl), g.add(float, dl));
        let held = g.min(hi, exact);
        let held = g.max(lo, held);
        let one = g.constant(1.0);
        g.lipschitz(held, one)
    }

    /// ((x - ux) dy - (y - uy) dx) s: positive right of u -> v, zero on
    /// the line through them in real numbers; square to an axis, that
    /// axis's face (`axis_face`).
    fn exact_line(&mut self, u: P, v: P, [x, y]: [NodeId; 2]) -> NodeId {
        if u[0] == v[0] {
            return self.axis_face(x, u[0], v[1] > u[1]);
        }
        if u[1] == v[1] {
            return self.axis_face(y, u[1], v[0] < u[0]);
        }
        let (dx, dy) = (diff(v[0], u[0]), diff(v[1], u[1]));
        let (ux, uy) = (self.g.constant(u[0]), self.g.constant(u[1]));
        let (xu, yu) = (self.g.sub(x, ux), self.g.sub(y, uy));
        let mut acc: Option<NodeId> = None;
        for (w, d, sg) in [(xu, &dy, 1.0), (yu, &dx, -1.0)] {
            for &c in d.iter().rev() {
                let c = self.g.constant(sg * c);
                let term = self.g.mul(w, c);
                acc = Some(match acc {
                    None => term,
                    Some(a) => self.g.add(a, term),
                });
            }
        }
        let size = |e: &E| ink_core::interval::sum_up(e.iter().map(|v| v.abs()));
        let s = self.g.constant(ink_core::interval::inv_norm_down(&[size(&dx), size(&dy)]));
        let a = acc.expect("u and v differ");
        self.g.mul(s, a)
    }
}

/// A round at corner p between edges a -> p and p -> b (the loop
/// counter-clockwise), radius r, convex or not: c - p (`ck`) bounded
/// outward, the line m beyond p, and its triangle's bound (`hold`: a float
/// convex polygon holding it); a straight side's normal is square to its
/// edge through c.
struct Round {
    p: P,
    a: P,
    b: P,
    r: f64,
    convex: bool,
    ck: [ink_core::interval::I; 2],
    m: [P; 2],
    hold: Vec<P>,
    /// A side along an arc (edge a -> p, then p -> b): its normal is the
    /// line through c along the unit direction from the arc's centre to c
    /// (bounds given), positive toward p times the sign given.
    arcs: [Option<([ink_core::interval::I; 2], f64)>; 2],
    /// Its exact inputs: a, b1, p, b2, b, r.
    inputs: [f64; 10],
}

/// A round at corner p: edge a -> p of bulge b1 (0: a line), then
/// p -> b of bulge b2 (the loop counter-clockwise), radius r. Each edge's
/// carrier is its line or circle; c is where both, offset by r to the
/// round's side, meet nearest p, its bounds certified by Krawczyk's test
/// (an interval Newton step mapping a box around the float root into
/// itself, so one root lies in it and in the step's image). Its tangent
/// points inside their edges, its triangle T (past the normals through
/// them: a line edge's square to it, an arc's through its centre, cut by a
/// line m beyond p), each arc's piece from its tangent point to p within T
/// (by the triangle of its end tangents), and T within `hold` (a float
/// convex polygon). None where an interval cannot decide it.
fn round_at(a: P, b1: f64, p: P, b2: f64, b: P, r: f64) -> Option<Round> {
    use ink_core::interval::{self as iv, I, add_o, mul_o, sub_o};
    if !(r > 0.0 && r.is_finite()) || b1.abs() >= 1.0 || b2.abs() >= 1.0 {
        return None;
    }
    let pt = I::point;
    let neg = |x: I| I::new(-x.hi, -x.lo);
    let dot = |u: [I; 2], v: [I; 2]| add_o(mul_o(u[0], v[0]), mul_o(u[1], v[1]));
    let cross = |u: [I; 2], v: [I; 2]| sub_o(mul_o(u[0], v[1]), mul_o(u[1], v[0]));
    let div = |u: I, v: I| iv::bin(BinOp::Div, u, v).0;
    let len = |u: [I; 2]| iv::un(UnOp::Sqrt, dot(u, u));
    let vsub = |u: [I; 2], v: [I; 2]| [sub_o(u[0], v[0]), sub_o(u[1], v[1])];
    let rel = |q: P| [sub_o(pt(q[0]), pt(p[0])), sub_o(pt(q[1]), pt(p[1]))];
    let total = |e: &E| e.iter().fold(pt(0.0), |s, &t| add_o(s, pt(t)));
    // A carrier: a line through p (its inward unit normal), or a circle
    // (its centre less p, radius, and 1 where the profile is inside it).
    // A circle also carries its arc's start less p and that start less
    // its centre (G / W), to read a distance to it about the start: no term
    // of the circle's size cancels, however flat the arc.
    #[derive(Clone, Copy)]
    enum Car {
        Line([I; 2]),
        Circle([I; 2], I, f64, [I; 2], [I; 2]),
    }
    let circle = |u: P, v: P, bb: f64, start: bool| {
        let DiscParts { w, gx, gy, .. } = disc_parts(u, v, bb);
        let gu = [div(total(&gx), pt(w)), div(total(&gy), pt(w))];
        // c - u; less p: plus u - p (none where u is p).
        let ur = if start { rel(u) } else { ptv0() };
        let o = vsub(ur, gu);
        let rad = div(iv::un(UnOp::Sqrt, total(&compress(&sum(&mul(&gx, &gx), &mul(&gy, &gy))))), pt(w));
        Car::Circle(o, rad, bb.signum(), ur, gu)
    };
    fn ptv0() -> [ink_core::interval::I; 2] {
        [ink_core::interval::I::point(0.0); 2]
    }
    fn add2(u: [ink_core::interval::I; 2], v: [ink_core::interval::I; 2]) -> [ink_core::interval::I; 2] {
        [ink_core::interval::add_o(u[0], v[0]), ink_core::interval::add_o(u[1], v[1])]
    }
    let unit = |u: [I; 2]| {
        let l = len(u);
        (l.lo > 0.0).then(|| [div(u[0], l), div(u[1], l)])
    };
    let left = |e: [I; 2]| [neg(e[1]), e[0]];
    let c1 = if b1 == 0.0 { Car::Line(left(unit(rel(a).map(neg))?)) } else { circle(a, p, b1, true) };
    let c2 = if b2 == 0.0 { Car::Line(left(unit(rel(b))?)) } else { circle(p, b, b2, false) };
    // Inward distance and its gradient, q less p.
    let g = |c: Car, q: [I; 2]| match c {
        Car::Line(m) => dot(q, m),
        Car::Circle(o, rad, s, u, gu) => {
            // rad - |q - o| = -(|q - u|^2 + 2 (q - u) . (u - o)) / (|q - o| + rad).
            let qu = vsub(q, u);
            let num = add_o(dot(qu, qu), mul_o(pt(2.0), dot(qu, gu)));
            let d = neg(div(num, add_o(len(vsub(q, o)), rad)));
            if s > 0.0 { d } else { neg(d) }
        }
    };
    let grad = |c: Car, q: [I; 2]| match c {
        Car::Line(m) => Some(m),
        Car::Circle(o, _, s, _, _) => {
            let u = unit(vsub(q, o))?;
            Some(if s > 0.0 { u.map(neg) } else { u })
        }
    };
    // Tangents at p, along the loop: convex where it turns left.
    let tangent = |c: Car, toward: [I; 2]| match c {
        Car::Line(_) => toward,
        Car::Circle(o, _, s, _, _) => {
            let po = o.map(neg);
            let t = [neg(po[1]), po[0]];
            if s > 0.0 { t } else { t.map(neg) }
        }
    };
    let turn = cross(tangent(c1, rel(a).map(neg)), tangent(c2, rel(b)));
    let convex = if turn.lo > 0.0 {
        true
    } else if turn.hi < 0.0 {
        false
    } else {
        return None;
    };
    let k = if convex { r } else { -r };
    // The float root: Newton's steps from the lines' answer.
    let mid = |x: [I; 2]| [x[0].mid(), x[1].mid()];
    let ptv = |q: P| [pt(q[0]), pt(q[1])];
    let (m1, m2) = (mid(grad(c1, ptv([0.0, 0.0]))?), mid(grad(c2, ptv([0.0, 0.0]))?));
    let den = 1.0 + m1[0] * m2[0] + m1[1] * m2[1];
    if den <= 0.0 {
        return None;
    }
    let mut c = [k * (m1[0] + m2[0]) / den, k * (m1[1] + m2[1]) / den];
    let step = |c: P| -> Option<(P, [P; 2])> {
        let q = ptv(c);
        let f = [g(c1, q).mid() - k, g(c2, q).mid() - k];
        let j = [mid(grad(c1, q)?), mid(grad(c2, q)?)];
        let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if det == 0.0 || !det.is_finite() {
            return None;
        }
        let y = [[j[1][1] / det, -j[0][1] / det], [-j[1][0] / det, j[0][0] / det]];
        Some(([c[0] - (y[0][0] * f[0] + y[0][1] * f[1]), c[1] - (y[1][0] * f[0] + y[1][1] * f[1])], y))
    };
    for _ in 0..30 {
        c = step(c)?.0;
    }
    let y = step(c)?.1;
    // Krawczyk: K = c - Y F(c) + (I - Y J(X)) (X - c) inside X.
    let mut ck = None;
    for e in [1e-12, 1e-9, 1e-6] {
        let h = e * r + 4.0 * f64::EPSILON * (c[0].abs() + c[1].abs() + r);
        let x = [I::new(c[0] - h, c[0] + h), I::new(c[1] - h, c[1] + h)];
        let fc = [sub_o(g(c1, ptv(c)), pt(k)), sub_o(g(c2, ptv(c)), pt(k))];
        let (j1, j2) = (grad(c1, x)?, grad(c2, x)?);
        let jx = [j1, j2];
        let yi = y.map(|row| row.map(pt));
        let dx = [sub_o(x[0], pt(c[0])), sub_o(x[1], pt(c[1]))];
        let mut kk = [pt(0.0); 2];
        for i in 0..2 {
            let yf = add_o(mul_o(yi[i][0], fc[0]), mul_o(yi[i][1], fc[1]));
            let mut acc = sub_o(pt(c[i]), yf);
            for jcol in 0..2 {
                // (I - Y J)_{i, jcol}
                let yj = add_o(mul_o(yi[i][0], jx[0][jcol]), mul_o(yi[i][1], jx[1][jcol]));
                let m = sub_o(pt(if i == jcol { 1.0 } else { 0.0 }), yj);
                acc = add_o(acc, mul_o(m, dx[jcol]));
            }
            kk[i] = acc;
        }
        if (0..2).all(|i| kk[i].lo > x[i].lo && kk[i].hi < x[i].hi) {
            ck = Some(kk);
            break;
        }
    }
    let ck = ck?;
    let foot = |car: Car| match car {
        Car::Line(m) => {
            let h = dot(ck, m);
            [sub_o(ck[0], mul_o(h, m[0])), sub_o(ck[1], mul_o(h, m[1]))]
        }
        Car::Circle(o, _, s, _, _) => {
            // c moved along its normal by its distance to the circle.
            let d = unit(vsub(ck, o)).unwrap_or([I::ENTIRE; 2]);
            let t = if s > 0.0 { g(car, ck) } else { neg(g(car, ck)) };
            add2(ck, [mul_o(d[0], t), mul_o(d[1], t)])
        }
    };
    let (f1, f2) = (foot(c1), foot(c2));
    // Each side's normal, positive toward p: a line edge's (q - c) . d (d
    // from its far end to p); an arc's s cross(c - o, q - o).
    let side = |car: Car, far: P| -> Option<(Box<dyn Fn([I; 2]) -> I>, Option<([I; 2], f64)>)> {
        match car {
            Car::Line(_) => {
                let d = rel(far).map(neg);
                Some((Box::new(move |q: [I; 2]| dot(vsub(q, ck), d)), None))
            }
            Car::Circle(o, _, _, _, _) => {
                let co = unit(vsub(ck, o))?;
                let at0 = cross(co, ck.map(neg));
                let sg = if at0.lo > 0.0 {
                    1.0
                } else if at0.hi < 0.0 {
                    -1.0
                } else {
                    return None;
                };
                let f = move |q: [I; 2]| {
                    let v = cross(co, vsub(q, ck));
                    if sg > 0.0 { v } else { neg(v) }
                };
                Some((Box::new(f), Some((co, sg))))
            }
        }
    };
    let (s1, arc1) = side(c1, a)?;
    let (s2, arc2) = side(c2, b)?;
    let zero = ptv([0.0, 0.0]);
    let ok = s1(rel(a)).hi < 0.0 && s2(rel(b)).hi < 0.0 && s1(zero).lo > 0.0 && s2(zero).lo > 0.0 && s1(f2).lo > 0.0 && s2(f1).lo > 0.0;
    // Each tangent point before p along its edge (on its arc, not where
    // the circle runs on past p).
    let (t1, t2) = (tangent(c1, rel(a).map(neg)), tangent(c2, rel(b)));
    if !ok || dot(f1, t1).hi >= 0.0 || dot(f2, t2).lo <= 0.0 {
        return None;
    }
    // m: through p - (c - p), square to c - p, p on its left.
    let cm = mid(ck);
    let q1 = [p[0] - cm[0], p[1] - cm[1]];
    let mut ends = [q1, [q1[0] + cm[1], q1[1] - cm[0]]];
    match orient(ends[0], ends[1], p) {
        0 => return None,
        o if o < 0 => ends.swap(0, 1),
        _ => {}
    }
    let [q1, q2] = ends;
    let qd = vsub(ptv(q2), ptv(q1));
    let r1 = rel(q1);
    let left_of = |q: [I; 2]| cross(qd, vsub(q, r1));
    if !(left_of(ck).lo > 0.0 && left_of(f1).lo > 0.0 && left_of(f2).lo > 0.0) {
        return None;
    }
    // Each arc's piece from its tangent point to p lies in the triangle of
    // its end tangents: where they meet, in T.
    for (car, f) in [(c1, f1), (c2, f2)] {
        if let Car::Circle(o, ..) = car {
            // p + l tp, tp and tf the tangents' unit directions there.
            let (np, nf) = (unit(o.map(neg))?, unit(vsub(f, o))?);
            let (tp, tf) = ([neg(np[1]), np[0]], [neg(nf[1]), nf[0]]);
            let det = cross(tf, tp);
            if det.lo <= 0.0 && det.hi >= 0.0 {
                return None;
            }
            let l = div(cross(tf, f), det);
            let t = [mul_o(l, tp[0]), mul_o(l, tp[1])];
            if !(s1(t).lo >= 0.0 && s2(t).lo >= 0.0 && left_of(t).lo > 0.0) {
                return None;
            }
        }
    }
    // T is the wedge at c between the rays along each normal into the
    // other side's half-plane, cut by m: bounded, its corners c and where
    // each ray meets m, only where m crosses both rays beyond c (else the
    // wedge runs on past m, and no hull of corners holds it).
    let meet = |dir: [I; 2], other: &dyn Fn([I; 2]) -> I| {
        let rise = sub_o(other(add2(ck, dir)), other(ck));
        let dir = if rise.lo > 0.0 {
            dir
        } else if rise.hi < 0.0 {
            dir.map(neg)
        } else {
            return None;
        };
        let den = cross(qd, dir);
        if den.lo <= 0.0 && den.hi >= 0.0 {
            return None;
        }
        let t = div(neg(left_of(ck)), den);
        (t.lo > 0.0).then(|| add2(ck, [mul_o(t, dir[0]), mul_o(t, dir[1])]))
    };
    let dir = |car: Car, far: P| match car {
        Car::Line(_) => {
            let d = rel(far).map(neg);
            [neg(d[1]), d[0]]
        }
        Car::Circle(o, ..) => vsub(ck, o),
    };
    let (x1, x2) = (meet(dir(c1, a), &*s2)?, meet(dir(c2, b), &*s1)?);
    let mut corners: Vec<P> = Vec::new();
    for v in [ck, x1, x2] {
        // Each corner of its box, rounded out on its own side (a point's
        // two ends are equal: each copy goes its own way).
        for (hx, hy) in [(false, false), (true, false), (true, true), (false, true)] {
            let out = |c: usize, high: bool| if high { (p[c] + v[c].hi).next_up() } else { (p[c] + v[c].lo).next_down() };
            let q = [out(0, hx), out(1, hy)];
            if !(q[0].is_finite() && q[1].is_finite()) {
                return None;
            }
            corners.push(q);
        }
    }
    let hold = ccw_hull(&corners);
    if hold.len() < 3 {
        return None;
    }
    let inputs = [a[0], a[1], b1, p[0], p[1], b2, b[0], b[1], r, 0.0];
    Some(Round { p, a, b, r, convex, ck, m: [q1, q2], hold, arcs: [arc1, arc2], inputs })
}

impl Doc {
    /// A flat profile with rounds (`round_at`): its
    /// corners' polygon built exactly (`exact_flat`), and each round
    /// certified apart from every other edge and round. In a
    /// round's triangle the profile is its disc (convex) or the disc's
    /// outside, and elsewhere the corners' polygon: max(f, min(D, -T)) or
    /// min(f, max(-D, T)), T the triangle's field (its normals' lines exact
    /// in real numbers through c, known within bounds), D the disc. None
    /// where it is polygon's construction.
    pub(super) fn rounded_polygon(&mut self, input: &[[f64; 4]]) -> Option<NodeId> {
        let n = input.len();
        if n < 3 || (0..n).any(|i| input[i][..2] == input[(i + 1) % n][..2]) {
            return None;
        }
        let pts: Vec<P> = input.iter().map(|p| [p[0], p[1]]).collect();
        let lo = (0..n).fold(0, |m, i| if lex_less(pts[i], pts[m]) { i } else { m });
        let ccw = match orient(pts[(lo + n - 1) % n], pts[lo], pts[(lo + 1) % n]) {
            0 => return None,
            o => o > 0,
        };
        let mut rounds = Vec::new();
        for (i, q) in input.iter().enumerate().filter(|(_, q)| q[3] != 0.0) {
            let (h, j) = ((i + n - 1) % n, (i + 1) % n);
            // Counter-clockwise: a -> p of bulge b1, p -> b of bulge b2.
            let (a, b1, b, b2) = if ccw { (pts[h], input[h][2], pts[j], q[2]) } else { (pts[j], -q[2], pts[h], -input[h][2]) };
            let rd = round_at(a, b1, pts[i], b2, b, q[3]);
            rounds.push((i, rd?));
        }
        // Each round's triangle apart from every edge but its own two, and
        // from every other round's.
        for (k, (i, rd)) in rounds.iter().enumerate() {
            for e in 0..n {
                let (u, v) = (pts[e], pts[(e + 1) % n]);
                if e == *i || (e + 1) % n == *i {
                    continue;
                }
                // Another edge, or the triangle holding its arc
                // (counter-clockwise).
                let (s, t, bb) = if ccw { (u, v, input[e][2]) } else { (v, u, -input[e][2]) };
                let other = if bb == 0.0 { vec![s, t] } else { Kite::of(s, t, bb)?.tri(s, t) };
                if !apart(&rd.hold, &other) {
                    return None;
                }
            }
            if rounds[k + 1..].iter().any(|(_, o)| !apart(&rd.hold, &o.hold)) {
                return None;
            }
        }
        let sharp: Vec<[f64; 3]> = input.iter().map(|p| [p[0], p[1], p[2]]).collect();
        let mut f = self.exact_flat(&sharp).ok()??;
        let (x, y) = (self.g.var(Var::X), self.g.var(Var::Y));
        let (mut cut, mut fill) = (Vec::new(), Vec::new());
        for (_, rd) in &rounds {
            // Each of its numbers by its role: c - p (1, 2), each arc
            // side's direction (3 to 6).
            let w = |g: &mut Graph, i: ink_core::interval::I, role: u64| g.within(i.lo, i.hi, i.mid().clamp(i.lo, i.hi), ident("round", &rd.inputs, role));
            let (px, py) = (self.g.constant(rd.p[0]), self.g.constant(rd.p[1]));
            let (xp, yp) = (self.g.sub(x, px), self.g.sub(y, py));
            // The disc, about p.
            let (cx, cy) = (w(&mut self.g, rd.ck[0], 1), w(&mut self.g, rd.ck[1], 2));
            let (dx, dy) = (self.g.sub(xp, cx), self.g.sub(yp, cy));
            let (sx, sy) = (self.g.un(UnOp::Square, dx), self.g.un(UnOp::Square, dy));
            let q = self.g.add(sx, sy);
            let dist = self.g.un(UnOp::Sqrt, q);
            let rr = self.g.constant(rd.r);
            let disc = self.g.sub(dist, rr);
            let one = self.g.constant(1.0);
            let disc = self.g.lipschitz(disc, one);
            // Each normal: ((p - c) . (p_corner - q)) s, positive toward
            // the corner, its differences exact (two terms each).
            let mut sides = Vec::new();
            for (side, (q, arc)) in [(rd.a, rd.arcs[0]), (rd.b, rd.arcs[1])].into_iter().enumerate() {
                let side = side as u64;
                // Along an arc: sg cross(c - o, (p - corner) - o) k, k at
                // most 1 / |c - o|.
                if let Some((co, sg)) = arc {
                    use ink_core::interval::{self as iv, I};
                    let big = |i: I| i.lo.abs().max(i.hi.abs());
                    let k = iv::inv_norm_down(&[big(co[0]), big(co[1])]) * (1.0 - 4.0 * f64::EPSILON);
                    let (ox, oy) = (w(&mut self.g, rd.ck[0], 1), w(&mut self.g, rd.ck[1], 2));
                    let (cx, cy) = (w(&mut self.g, co[0], 3 + 2 * side), w(&mut self.g, co[1], 4 + 2 * side));
                    let (qx, qy) = (self.g.sub(xp, ox), self.g.sub(yp, oy));
                    let (t1, t2) = (self.g.mul(cx, qy), self.g.mul(cy, qx));
                    let cr = self.g.sub(t1, t2);
                    let kk = self.g.constant(sg * k);
                    let side = self.g.mul(kk, cr);
                    sides.push(self.g.neg(side));
                    continue;
                }
                let (ex, ey) = (diff(rd.p[0], q[0]), diff(rd.p[1], q[1]));
                let size = |e: &E| ink_core::interval::sum_up(e.iter().map(|v| v.abs()));
                let s = ink_core::interval::inv_norm_down(&[size(&ex), size(&ey)]);
                let tot = |e: &E| e.iter().fold(ink_core::interval::I::point(0.0), |a, &t| ink_core::interval::add_o(a, ink_core::interval::I::point(t)));
                let at = ink_core::interval::add_o(ink_core::interval::mul_o(rd.ck[0], tot(&ex)), ink_core::interval::mul_o(rd.ck[1], tot(&ey)));
                let mut acc = None;
                for (v, e) in [(xp, &ex), (yp, &ey)] {
                    for &t in e.iter().rev() {
                        let t = self.g.constant(t);
                        let term = self.g.mul(v, t);
                        acc = Some(match acc {
                            None => term,
                            Some(a) => self.g.add(a, term),
                        });
                    }
                }
                let at = w(&mut self.g, at, 7 + side);
                let lin = self.g.sub(acc?, at);
                let s = self.g.constant(s);
                let side = self.g.mul(s, lin);
                sides.push(self.g.neg(side));
            }
            sides.push(self.exact_line(rd.m[0], rd.m[1], [x, y]));
            let tri = self.reduce(BinOp::Max, &sides);
            if rd.convex {
                let out = self.g.neg(tri);
                cut.push(self.g.min(disc, out));
            } else {
                let nd = self.g.neg(disc);
                fill.push(self.g.max(nd, tri));
            }
        }
        if !cut.is_empty() {
            cut.push(f);
            f = self.reduce(BinOp::Max, &cut);
        }
        if !fill.is_empty() {
            fill.push(f);
            f = self.reduce(BinOp::Min, &fill);
            // Its box outermost (a fill's terms reach past it): its
            // corners', and each arc's kite corner (the triangle of its
            // chord and that float corner proven to hold the arc).
            let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
            for (k, q) in input.iter().enumerate() {
                let mut ps = vec![[q[0], q[1]]];
                if q[2] != 0.0 {
                    let v = input[(k + 1) % n];
                    ps.push(Kite::of([q[0], q[1]], [v[0], v[1]], q[2])?.tf);
                }
                for p in ps {
                    for c in 0..2 {
                        lo[c] = lo[c].min(p[c]);
                        hi[c] = hi[c].max(p[c]);
                    }
                }
            }
            let mut terms = vec![f];
            for (v, c) in [(x, 0), (y, 1)] {
                terms.push(self.axis_face(v, hi[c], true));
                terms.push(self.axis_face(v, lo[c], false));
            }
            f = self.reduce(BinOp::Max, &terms);
        }
        Some(f)
    }
}

/// A field's exact value at a point, for tests: a dyadic rational (an
/// integer times a power of two), the graph's operations taken exactly
/// (sums, products, min, max, negation of its constants and coordinates).
#[cfg(test)]
fn dyadic_at(g: &Graph, root: NodeId, p: [f64; 2]) -> (big::Z, i32) {
    use ink_core::graph::Node;
    type D = (big::Z, i32);
    fn of(v: f64) -> D {
        let (m, e, neg) = big::parts(v);
        (big::Z::of(m, neg, 0), e)
    }
    fn align(a: &D, b: &D) -> (big::Z, big::Z, i32) {
        let e = a.1.min(b.1);
        let up = |d: &D| d.0.mul(&big::Z::of(1, false, (d.1 - e) as u32));
        (up(a), up(b), e)
    }
    let mut vals: Vec<Option<D>> = vec![None; root as usize + 1];
    for id in 0..=root as usize {
        let get = |vals: &Vec<Option<D>>, a: NodeId| vals[a as usize].clone().expect("an exact value");
        let v = match &g.nodes[id] {
            Node::Const(c) => Some(of(*c)),
            Node::Var(Var::X) => Some(of(p[0])),
            Node::Var(Var::Y) => Some(of(p[1])),
            Node::Var(_) => Some(of(0.0)),
            Node::Tag(_, a) | Node::Lipschitz(a, _, _) => vals[*a as usize].clone(),
            Node::Un(UnOp::Neg, a) => vals[*a as usize].clone().map(|(z, e)| (big::Z::of(0, false, 0).sub(&z), e)),
            Node::Bin(op @ (BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Min | BinOp::Max), a, b) if vals[*a as usize].is_some() && vals[*b as usize].is_some() => {
                let (x, y) = (get(&vals, *a), get(&vals, *b));
                Some(match op {
                    BinOp::Mul => (x.0.mul(&y.0), x.1 + y.1),
                    _ => {
                        let (zx, zy, e) = align(&x, &y);
                        match op {
                            BinOp::Add => (zx.add(&zy), e),
                            BinOp::Sub => (zx.sub(&zy), e),
                            BinOp::Min => (if zx.sub(&zy).sign() <= 0 { zx } else { zy }, e),
                            _ => (if zx.sub(&zy).sign() >= 0 { zx } else { zy }, e),
                        }
                    }
                })
            }
            _ => None,
        };
        vals[id] = v;
    }
    vals[root as usize].clone().expect("the field's exact value")
}

/// Whether |f(q) - f(p)| <= d exactly (f the field at root).
#[cfg(test)]
pub(super) fn differs_by_at_most(g: &Graph, root: NodeId, p: [f64; 2], q: [f64; 2], d: f64) -> bool {
    let (p, q) = (&dyadic_at(g, root, p), &dyadic_at(g, root, q));
    let (m, e, neg) = big::parts(d);
    let dd = (big::Z::of(m, neg, 0), e);
    let least = p.1.min(q.1).min(dd.1);
    let up = |z: &big::Z, e: i32| z.mul(&big::Z::of(1, false, (e - least) as u32));
    let diff = up(&q.0, q.1).sub(&up(&p.0, p.1));
    let dz = up(&dd.0, dd.1);
    diff.sub(&dz).sign() <= 0 && big::Z::of(0, false, 0).sub(&diff).sub(&dz).sign() <= 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orient_is_exact_where_floats_round() {
        // Three points on a line a float determinant reads as turning.
        let (a, b, c) = ([0.5, 0.5], [12.0, 12.0], [24.0, 24.0]);
        assert_eq!(orient(a, b, c), 0);
        let e = f64::EPSILON;
        assert_eq!(orient([0.5 + e, 0.5], b, c), -1);
        assert_eq!(orient([0.5, 0.5 + e], b, c), 1);
        // Far out, where a float determinant cancels to 0.
        let k = 6e16;
        assert_eq!(orient([k, k], [k + 8.0, k], [k, k + 8.0]), 1);
    }

    #[test]
    fn a_corner_off_its_neighbours_line_is_kept() {
        let pts = [[0.0, 0.0, 0.0], [500.0, -5e-7, 0.0], [1000.0, 0.0, 0.0], [1000.0, 1000.0, 0.0], [0.0, 1000.0, 0.0]];
        assert_eq!(points(&pts).unwrap().0.len(), 5);
        let on = [[0.0, 0.0, 0.0], [500.0, 0.0, 0.0], [1000.0, 0.0, 0.0], [1000.0, 1000.0, 0.0], [0.0, 1000.0, 0.0]];
        assert_eq!(points(&on).unwrap().0.len(), 4);
    }

    #[test]
    fn an_arc_keeps_its_corners_and_turns_with_its_loop() {
        // A corner in line with its neighbours is kept where an edge on
        // from it is an arc.
        let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.1], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]];
        let (p, b) = points(&pts).unwrap();
        assert_eq!((p.len(), b), (5, vec![0.0, 0.1, 0.0, 0.0, 0.0]));
        // Clockwise, reversed: the arc runs the other way, its bulge negated.
        let cw = [[0.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 2.0, 0.0], [2.0, 0.0, 0.3]];
        let (p, b) = points(&cw).unwrap();
        let k = b.iter().position(|&b| b != 0.0).unwrap();
        assert_eq!((p[k], p[(k + 1) % 4], b[k]), ([0.0, 0.0], [2.0, 0.0], -0.3));
    }

    #[test]
    fn an_arc_s_triangle_must_be_clear_of_other_edges() {
        // A deep arc on the bottom of a thin rectangle reaches its top.
        let (p, b) = points(&[[0.0, 0.0, -0.5], [10.0, 0.0, 0.0], [10.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).unwrap();
        assert!(kites(&p, &b).is_err());
        let (p, b) = points(&[[0.0, 0.0, -0.05], [10.0, 0.0, 0.0], [10.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).unwrap();
        assert!(kites(&p, &b).is_ok());
    }

    #[test]
    fn orientation_is_exact_across_any_magnitudes() {
        // A convex loop from 1e-100 to 1e70: one scale to its largest
        // coordinate leaves products of its small differences below the
        // floats, where expansions round; it read as doubling back.
        let pts = [[0.0, 0.0, 0.0], [1e-100, 0.0, 0.0], [2e-100, 1e-100, 0.0], [1e70, 1e70, 0.0], [0.0, 1e70, 0.0]];
        let big = 1e70f64;
        let sc = 2f64.powi(-(big.log2().ceil() as i32));
        let scaled: Vec<[f64; 3]> = pts.iter().map(|p| [p[0] * sc, p[1] * sc, 0.0]).collect();
        let (p, _) = points(&scaled).unwrap();
        assert_eq!(p.len(), 5);
        // The integers agree with the expansions where both are exact.
        let mut r = 0x9e3779b97f4a7c15u64;
        let mut f = || {
            r ^= r << 13;
            r ^= r >> 7;
            r ^= r << 17;
            (r >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for _ in 0..2000 {
            let (a, b) = ([f(), f()], [f(), f()]);
            let c = [a[0] + (b[0] - a[0]) * 0.5 + f() * 1e-17, a[1] + (b[1] - a[1]) * 0.5];
            assert_eq!(big::orient(a, b, c), orient(a, b, c), "{a:?} {b:?} {c:?}");
            assert_eq!(big::onward(a, b, c), onward(a, b, c));
        }
    }

    #[test]
    fn a_unit_normal_s_scale_is_proven_at_most_one() {
        // s^2 (x^2 + y^2) <= 1 in integers, for s = inv_norm_down.
        let check = |x: f64, y: f64| {
            let s = ink_core::interval::inv_norm_down(&[x, y]);
            let [zs, zx, zy] = big::ints([s, x, y]);
            let lhs = zs.mul(&zs).mul(&zx.mul(&zx).add(&zy.mul(&zy)));
            // 1 at the same scale: s, x and y read at their least power e,
            // so s^2 (x^2 + y^2) carries 2^(4e); 1 is 2^(-4e) of that unit.
            let least = [s, x, y].iter().filter(|v| **v != 0.0).map(|v| big::parts(*v).1).min().unwrap();
            let one = big::Z::of(1, false, (-4 * least) as u32);
            assert!(lhs.sub(&one).sign() <= 0, "{x} {y}: s {s}");
        };
        check(3.0, 4.0);
        check(1.0, 0.0);
        check(0.1, 0.7);
        for k in 1..500 {
            check(k as f64 * 0.37, 1.0 / k as f64);
        }
        assert_eq!(ink_core::interval::inv_norm_down(&[0.0, 2.0]), 0.5);
    }

    #[test]
    fn a_loop_that_crosses_or_touches_itself_is_refused() {
        let bow = [[0.0, 0.0, 0.0], [1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(points(&bow).is_err());
        // A corner touching the opposite edge exactly.
        let touch = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        assert!(points(&touch).is_err());
        // The same corner a float above it: simple.
        let near = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 2.0, 0.0], [1.0, 1e-100, 0.0], [0.0, 2.0, 0.0]];
        assert!(points(&near).is_ok());
    }
}
