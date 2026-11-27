//! A cell holding two materials across a face: its element as the laminate
//! they make.
//!
//! A cell's element is one trilinear field, so it cannot hold the jump in
//! strain a face between two materials makes: one stiffness for the cell
//! places the face on a plane of the cells (its centre's material), and
//! one weighed over where the materials lie (each sample's, a mean) holds a
//! face across the load as stiff as one along it, the soft side's strain
//! unseen (a 0.1 mm insulating layer read first-order hot). Two materials
//! side by side across a flat face, each strained uniformly, are exactly a
//! laminate: the strain along the face is one in both, the traction across
//! it one, and the jump between their strains follows from those (a weak
//! discontinuity, condensed into the cell). So the cell's element takes the
//! laminate's stiffness, which carries a pull across the face in series and
//! one along it in parallel exactly, and each material's stress is its own
//! strain, recovered from the cell's mean strain, times its own stiffness.
//! A face oblique or curved across the cell takes its fitted plane's normal.

use ink_core::V3;

/// The two materials of a cell and how they lie (`measure`).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Laminate {
    /// The face's unit normal, out of material `a` (cell coordinates, a
    /// cube's, so a direction is the same in space).
    pub n: V3,
    /// The face's plane, `n . xi = d`: `a` below it, `b` above. None where
    /// the cell holds more than one face (a layer within it), which reads
    /// take either material at any point (`Laminate::fields`).
    pub d: Option<f64>,
    /// Material `a`'s share of the cell's material.
    pub fa: f64,
    /// Each material's coefficient (E, or k) and Poisson's ratio.
    pub a: (f64, f64),
    pub b: (f64, f64),
    /// Each material's yield strength where the study's is a field and it
    /// is one value over each (`Solution::ratio_between_reads`).
    #[serde(default)]
    pub yields: Option<[f64; 2]>,
}

/// A dense row-major matrix of `r` rows.
type Mat = Vec<f64>;

fn mul(a: &[f64], b: &[f64], n: usize, k: usize, m: usize) -> Mat {
    let mut out = vec![0.0; n * m];
    for i in 0..n {
        for q in 0..k {
            let v = a[i * k + q];
            if v != 0.0 {
                for j in 0..m {
                    out[i * m + j] += v * b[q * m + j];
                }
            }
        }
    }
    out
}

fn transpose(a: &[f64], n: usize, m: usize) -> Mat {
    let mut out = vec![0.0; n * m];
    for i in 0..n {
        for j in 0..m {
            out[j * n + i] = a[i * m + j];
        }
    }
    out
}

/// `a`^-1 `b` for `a` n x n (symmetric positive definite here), `b` n x m.
fn solve(a: &[f64], b: &[f64], n: usize, m: usize) -> Mat {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    for c in 0..n {
        let p = (c..n).max_by(|&x, &y| a[x * n + c].abs().total_cmp(&a[y * n + c].abs())).unwrap_or(c);
        for j in 0..n {
            a.swap(c * n + j, p * n + j);
        }
        for j in 0..m {
            b.swap(c * m + j, p * m + j);
        }
        let piv = a[c * n + c];
        if piv == 0.0 {
            continue;
        }
        for r in 0..n {
            if r != c {
                let f = a[r * n + c] / piv;
                if f != 0.0 {
                    for j in 0..n {
                        a[r * n + j] -= f * a[c * n + j];
                    }
                    for j in 0..m {
                        b[r * m + j] -= f * b[c * m + j];
                    }
                }
            }
        }
    }
    for r in 0..n {
        let piv = a[r * n + r];
        for j in 0..m {
            b[r * m + j] = if piv != 0.0 { b[r * m + j] / piv } else { 0.0 };
        }
    }
    b
}

impl Laminate {
    /// The size of the strain (Voigt, engineering shear) or gradient.
    fn nv(elastic: bool) -> usize {
        if elastic { 6 } else { 3 }
    }

    /// A material's stiffness (E D(nu), Voigt) or conductivity (k I).
    fn stiffness(elastic: bool, (c, nu): (f64, f64)) -> Mat {
        if elastic { crate::element::d_matrix(nu).iter().flatten().map(|v| c * v).collect() } else { (0..9usize).map(|i| if i.is_multiple_of(4) { c } else { 0.0 }).collect() }
    }

    /// The strain a jump `c` across the face makes, sym(c n) (nv x nc), and
    /// whose transpose takes a stress to its traction on the face.
    fn jump(&self, elastic: bool) -> (Mat, usize) {
        let n = self.n;
        if elastic {
            #[rustfmt::skip]
            let s = vec![
                n[0], 0.0, 0.0,
                0.0, n[1], 0.0,
                0.0, 0.0, n[2],
                n[1], n[0], 0.0,
                0.0, n[2], n[1],
                n[2], 0.0, n[0],
            ];
            (s, 3)
        } else {
            (n.to_vec(), 1)
        }
    }

    /// Each material's strain from the cell's mean, as matrices (A_a, A_b):
    /// eps_a = eps + fb S c, eps_b = eps - fa S c, c the jump that makes the
    /// traction across the face one: S^T (fb D_a + fa D_b) S c = S^T (D_b - D_a) eps.
    fn localization(&self, elastic: bool) -> (Mat, Mat) {
        let nv = Self::nv(elastic);
        let (s, nc) = self.jump(elastic);
        let (da, db) = (Self::stiffness(elastic, self.a), Self::stiffness(elastic, self.b));
        let (fa, fb) = (self.fa, 1.0 - self.fa);
        let st = transpose(&s, nv, nc);
        let mix: Mat = da.iter().zip(&db).map(|(a, b)| fb * a + fa * b).collect();
        let m = mul(&mul(&st, &mix, nc, nv, nv), &s, nc, nv, nc);
        let diff: Mat = db.iter().zip(&da).map(|(b, a)| b - a).collect();
        let c = solve(&m, &mul(&st, &diff, nc, nv, nv), nc, nv);
        let sc = mul(&s, &c, nv, nc, nv);
        let eye = |i: usize| if i.is_multiple_of(nv + 1) { 1.0 } else { 0.0 };
        let aa = (0..nv * nv).map(|i| eye(i) + fb * sc[i]).collect();
        let ab = (0..nv * nv).map(|i| eye(i) - fa * sc[i]).collect();
        (aa, ab)
    }

    /// The laminate's stiffness (nv x nv, symmetric): fa D_a A_a + fb D_b A_b.
    pub fn effective(&self, elastic: bool) -> Mat {
        let nv = Self::nv(elastic);
        let (aa, ab) = self.localization(elastic);
        let (da, db) = (Self::stiffness(elastic, self.a), Self::stiffness(elastic, self.b));
        let (ka, kb) = (mul(&da, &aa, nv, nv, nv), mul(&db, &ab, nv, nv, nv));
        let k: Mat = ka.iter().zip(&kb).map(|(a, b)| self.fa * a + (1.0 - self.fa) * b).collect();
        // Symmetric exactly (it is, but for rounding).
        (0..nv * nv).map(|i| 0.5 * (k[i] + k[(i % nv) * nv + i / nv])).collect()
    }

    /// The elastic laminate's stiffness as an array.
    pub fn effective6(&self) -> [[f64; 6]; 6] {
        let k = self.effective(true);
        std::array::from_fn(|i| std::array::from_fn(|j| k[i * 6 + j]))
    }

    /// Each material's stress (or flux, `-k grad`: the first 3) for the
    /// cell's mean strain (or temperature gradient) `e`, a then b.
    pub fn stresses(&self, elastic: bool, e: &[f64]) -> [[f64; 6]; 2] {
        let nv = Self::nv(elastic);
        let (aa, ab) = self.localization(elastic);
        let one = |m: (f64, f64), a: &Mat| -> [f64; 6] {
            let d = Self::stiffness(elastic, m);
            let da = mul(&d, a, nv, nv, nv);
            let sign = if elastic { 1.0 } else { -1.0 };
            std::array::from_fn(|r| if r < nv { sign * (0..nv).map(|q| da[r * nv + q] * e[q]).sum::<f64>() } else { 0.0 })
        };
        [one(self.a, &aa), one(self.b, &ab)]
    }

    /// The materials a point `xi` (cell coordinates) may be read as: its
    /// side of the face, or both where the cell holds more than one face.
    /// (a, b): whether each.
    pub fn at(&self, xi: V3) -> (bool, bool) {
        match self.d {
            Some(d) => {
                let below = self.n[0] * xi[0] + self.n[1] * xi[1] + self.n[2] * xi[2] < d;
                (below, !below)
            }
            None => (true, true),
        }
    }

    /// Which material a key (coefficient, nu) is: a's or not.
    pub fn is_a(&self, k: (f64, f64), same: &dyn Fn(f64, f64) -> bool) -> bool {
        same(k.0, self.a.0) && same(k.1, self.a.1)
    }

    /// The coefficient a cell of this laminate counts as (its materials'
    /// mean, by share): what scales its holds' penalties and the like.
    pub fn mean(&self) -> (f64, f64) {
        (self.fa * self.a.0 + (1.0 - self.fa) * self.b.0, self.fa * self.a.1 + (1.0 - self.fa) * self.b.1)
    }
}

/// The lines a cell's materials are measured along (`measure`): per axis,
/// the 16 through its samples' columns, each read at these places along
/// it: its two faces (a hair in) and 15 between, an eighth of the cell
/// apart, the samples among them (a layer thinner than that may lie
/// between them).
pub const ALONG: [f64; 17] = {
    let mut t = [0.0; 17];
    let mut i = 0;
    while i < 17 {
        t[i] = -1.0 + i as f64 / 8.0;
        i += 1;
    }
    t[0] += 1e-6;
    t[16] -= 1e-6;
    t
};

/// Bisection steps placing a face between two places along a line.
const BISECT: usize = 40;

/// A cell's two materials as a laminate, from `key`, the material at points
/// of the cell (its coordinates, -1 to 1): a pair of values (coefficient,
/// nu), alike within `same`. `weight`: each of its 64 samples' share of
/// its material (0 outside its piece; `sample_xi(4)`'s order). None where
/// the lines through it find one material, or more than two.
///
/// Along the axis whose lines most often cross a face, each line's faces
/// are placed by bisection between the places it is read at; each
/// sample's share of material a is the part of its span along the line
/// that holds it, and a's share is theirs, weighed. The face's normal is
/// the plane fitted through where the lines cross it (each one face): its
/// mean place over the lines is exact for a plane however it leans, and
/// the fit is its normal.
pub fn measure(key: &dyn Fn(&[V3]) -> Vec<(f64, f64)>, same: &dyn Fn(f64, f64) -> bool, weight: &[f64]) -> Option<Laminate> {
    let lateral = |j: usize| -0.75 + 0.5 * j as f64;
    // Point k along line (axis a, column j).
    let point = |a: usize, j: usize, t: f64| -> V3 {
        let (u, v) = ((a + 1) % 3, (a + 2) % 3);
        let mut p = [0.0; 3];
        p[a] = t;
        p[u] = lateral(j % 4);
        p[v] = lateral(j / 4);
        p
    };
    let pts: Vec<V3> = (0..3).flat_map(|a| (0..16).flat_map(move |j| ALONG.iter().map(move |&t| point(a, j, t)))).collect();
    let keys = key(&pts);
    let eq = |p: (f64, f64), q: (f64, f64)| same(p.0, q.0) && same(p.1, q.1);
    let mut materials: Vec<(f64, f64)> = Vec::new();
    for &k in &keys {
        if !materials.iter().any(|&m| eq(m, k)) {
            materials.push(k);
        }
        if materials.len() > 2 {
            return None;
        }
    }
    let [ma, mb] = materials[..] else { return None };
    let is_a = |k: (f64, f64)| eq(k, ma);
    let at = |a: usize, j: usize, i: usize| keys[(a * 16 + j) * ALONG.len() + i];
    // The axis whose lines cross a face most often.
    let crossings = |a: usize| (0..16).filter(|&j| (1..ALONG.len()).any(|i| is_a(at(a, j, i)) != is_a(at(a, j, i - 1)))).count();
    let axis = (0..3).max_by_key(|&a| (crossings(a), std::cmp::Reverse(a)))?;
    if crossings(axis) == 0 {
        return None;
    }
    // Each line's spans between changes, bisected together.
    let mut spans: Vec<(usize, usize, f64, f64)> = Vec::new();
    for j in 0..16 {
        for i in 1..ALONG.len() {
            if is_a(at(axis, j, i)) != is_a(at(axis, j, i - 1)) {
                spans.push((j, i, ALONG[i - 1], ALONG[i]));
            }
        }
    }
    for _ in 0..BISECT {
        let mids: Vec<V3> = spans.iter().map(|&(j, _, lo, hi)| point(axis, j, 0.5 * (lo + hi))).collect();
        let k = key(&mids);
        for (s, k) in spans.iter_mut().zip(k) {
            let mid = 0.5 * (s.2 + s.3);
            // Material at the span's low end: a's or not.
            if is_a(k) == is_a(at(axis, s.0, s.1 - 1)) { s.2 = mid } else { s.3 = mid }
        }
    }
    // Per line: material a's intervals along it, from -1 to 1.
    let mut share = [0.0f64; 64];
    let mut faces: Vec<(f64, f64, f64)> = Vec::new();
    let mut layered = false;
    for j in 0..16 {
        let cuts: Vec<f64> = spans.iter().filter(|s| s.0 == j).map(|s| 0.5 * (s.2 + s.3)).collect();
        if cuts.len() > 1 {
            layered = true;
        }
        if cuts.len() == 1 {
            faces.push((lateral(j % 4), lateral(j / 4), cuts[0]));
        }
        // Material a's pieces of [-1, 1].
        let mut a_parts = Vec::new();
        let mut start = -1.0;
        let mut in_a = is_a(at(axis, j, 0));
        for &c in cuts.iter().chain(std::iter::once(&1.0)) {
            if in_a {
                a_parts.push((start, c));
            }
            start = c;
            in_a = !in_a;
        }
        for k in 0..4 {
            let (lo, hi) = (-1.0 + 0.5 * k as f64, -0.5 + 0.5 * k as f64);
            let len: f64 = a_parts.iter().map(|&(p, q): &(f64, f64)| (q.min(hi) - p.max(lo)).max(0.0)).sum();
            // Sample index: along `axis` k, laterally j.
            let mut idx = [0usize; 3];
            idx[axis] = k;
            idx[(axis + 1) % 3] = j % 4;
            idx[(axis + 2) % 3] = j / 4;
            share[idx[0] + 4 * idx[1] + 16 * idx[2]] = len / 0.5;
        }
    }
    let total: f64 = weight.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let fa = (weight.iter().zip(&share).map(|(w, s)| w * s).sum::<f64>() / total).clamp(0.0, 1.0);
    // The face's plane: t = c0 + c1 u + c2 v by least squares, where three
    // or more lines not in a row cross it; else across the axis.
    let mut n = [0.0; 3];
    n[axis] = 1.0;
    let mut d = faces.iter().map(|f| f.2).sum::<f64>() / faces.len().max(1) as f64;
    if faces.len() >= 3 {
        let (mut ata, mut atb) = ([0.0; 9], [0.0; 3]);
        for &(u, v, t) in &faces {
            let r = [1.0, u, v];
            for p in 0..3 {
                for q in 0..3 {
                    ata[p * 3 + q] += r[p] * r[q];
                }
                atb[p] += r[p] * t;
            }
        }
        let det = crate::element::m3_det(&[[ata[0], ata[1], ata[2]], [ata[3], ata[4], ata[5]], [ata[6], ata[7], ata[8]]]);
        if det.abs() > 1e-9 {
            let c = solve(&ata, &atb, 3, 1);
            let mut m = [0.0; 3];
            m[axis] = 1.0;
            m[(axis + 1) % 3] = -c[1];
            m[(axis + 2) % 3] = -c[2];
            let len = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
            n = m.map(|v| v / len);
            d = c[0] / len;
        }
    }
    // Oriented out of a: a lies below where the first line starts in it.
    let below_is_a = {
        let j = (0..16).find(|&j| spans.iter().any(|s| s.0 == j)).unwrap_or(0);
        let p = point(axis, j, ALONG[0]);
        let below = n[0] * p[0] + n[1] * p[1] + n[2] * p[2] < d;
        below == is_a(at(axis, j, 0))
    };
    if !below_is_a {
        n = n.map(|v| -v);
        d = -d;
    }
    Some(Laminate { n, d: (!layered && !faces.is_empty()).then_some(d), fa, a: ma, b: mb, yields: None })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lam(n: V3, fa: f64, a: (f64, f64), b: (f64, f64)) -> Laminate {
        Laminate { n, d: Some(0.0), fa, a, b, yields: None }
    }

    #[test]
    fn a_laminate_is_its_layers_in_series_across_and_in_parallel_along() {
        // nu 0: across the face, 1 / (fa / Ea + fb / Eb); along it, fa Ea + fb Eb.
        let l = lam([1.0, 0.0, 0.0], 0.3, (2000.0, 0.0), (200000.0, 0.0));
        let k = l.effective(true);
        let series = 1.0 / (0.3 / 2000.0 + 0.7 / 200000.0);
        let parallel = 0.3 * 2000.0 + 0.7 * 200000.0;
        assert!((k[0] / series - 1.0).abs() < 1e-12, "{}", k[0]);
        assert!((k[7] / parallel - 1.0).abs() < 1e-12, "{}", k[7]);
        // Shear across the face in series (G = E / 2 at nu 0), along it in parallel.
        assert!((k[3 * 6 + 3] / (0.5 * series) - 1.0).abs() < 1e-12 && (k[4 * 6 + 4] / (0.5 * parallel) - 1.0).abs() < 1e-12);
        // Each material's stress: one across the face; along it, its own E times the strain.
        let s = l.stresses(true, &[1e-3, 1e-3, 0.0, 0.0, 0.0, 0.0]);
        assert!((s[0][0] - s[1][0]).abs() < 1e-9 * s[0][0].abs() && (s[0][1] - 2.0).abs() < 1e-9 && (s[1][1] - 200.0).abs() < 1e-9, "{s:?}");
        // Heat: harmonic across, arithmetic along; the flux across one.
        let h = l.effective(false);
        assert!((h[0] / series - 1.0).abs() < 1e-12 && (h[4] / parallel - 1.0).abs() < 1e-12);
        // Turned, the same tensor turned: the normal's own stiffness the series one.
        let r = 0.5f64.sqrt();
        let t = lam([r, r, 0.0], 0.3, (2000.0, 0.25), (200000.0, 0.35)).effective(false);
        let nn = (0..3).map(|i| (0..3).map(|j| [r, r, 0.0][i] * t[i * 3 + j] * [r, r, 0.0][j]).sum::<f64>()).sum::<f64>();
        assert!((nn / series - 1.0).abs() < 1e-12, "{nn}");
    }

    #[test]
    fn with_poisson_s_ratio_the_layers_share_the_traction_and_the_strain_along_the_face() {
        // nu 0.25 and 0.35, the face across x, a general mean strain (shears
        // engineering, xy yz zx): what makes it a laminate, each condition
        // by itself. The traction on the face (xx, xy, zx) one in both; the
        // strain along it (yy, zz, yz) one in both; the shares' mean strain
        // the cell's, their mean stress its stiffness times it.
        let (fa, a, b) = (0.3, (2000.0, 0.25), (200000.0, 0.35));
        let l = lam([1.0, 0.0, 0.0], fa, a, b);
        let e = [1e-3, -2e-4, 3e-4, 5e-4, -1e-4, 2e-4];
        let [sa, sb] = l.stresses(true, &e);
        let strain = |m: (f64, f64), s: &[f64; 6]| solve(&Laminate::stiffness(true, m), s, 6, 1);
        let (ea, eb) = (strain(a, &sa), strain(b, &sb));
        let scale = sa.iter().chain(&sb).fold(0.0f64, |m, v| m.max(v.abs()));
        for i in [0, 3, 5] {
            assert!((sa[i] - sb[i]).abs() < 1e-12 * scale, "traction {i}: {sa:?} {sb:?}");
        }
        for i in [1, 2, 4] {
            assert!((ea[i] - e[i]).abs() < 1e-12 && (eb[i] - e[i]).abs() < 1e-12, "along {i}: {ea:?} {eb:?}");
        }
        let k = l.effective(true);
        for i in 0..6 {
            assert!((fa * ea[i] + (1.0 - fa) * eb[i] - e[i]).abs() < 1e-12, "mean strain {i}");
            let ke: f64 = (0..6).map(|j| k[i * 6 + j] * e[j]).sum();
            assert!((fa * sa[i] + (1.0 - fa) * sb[i] - ke).abs() < 1e-9 * scale, "mean stress {i}: {ke}");
        }
        // Poisson's coupling is in it: a pull across the face alone narrows
        // each layer along it by its own ratio, held to one strain.
        assert!(k[1].abs() > 1.0 && sa[1] != sb[1]);
    }

    #[test]
    fn a_face_s_share_and_normal_are_measured_off_the_samples() {
        let same = |a: f64, b: f64| (a - b).abs() <= 1e-9 * a.abs().max(b.abs());
        // A plane through the cell leaning in two directions: a below.
        let (n, d) = ([0.2f64, -0.1, 1.0], 0.13);
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let key = |p: &[V3]| p.iter().map(|p| if n[0] * p[0] + n[1] * p[1] + n[2] * p[2] < d { (5.0, 0.3) } else { (7.0, 0.3) }).collect();
        let l = measure(&key, &same, &[1.0; 64]).unwrap();
        // The share of a cube below a plane leaving through its top and bottom
        // only: its height at the centre over 2.
        let want = (1.0 + d / n[2]) / 2.0;
        assert!((l.fa - want).abs() < 1e-9, "{} {want}", l.fa);
        assert!(l.a.0 == 5.0 && (0..3).all(|i| (l.n[i] - n[i] / len).abs() < 1e-9), "{l:?}");
        assert!(l.at([0.0, 0.0, -0.5]) == (true, false) && l.at([0.0, 0.0, 0.5]) == (false, true));
        // A layer within the cell: two faces, either material anywhere.
        let key = |p: &[V3]| p.iter().map(|p| if p[0].abs() < 0.3 { (5.0, 0.3) } else { (7.0, 0.3) }).collect();
        let l = measure(&key, &same, &[1.0; 64]).unwrap();
        let fa = if l.a.0 == 5.0 { l.fa } else { 1.0 - l.fa };
        assert!((fa - 0.3).abs() < 1e-9 && l.d.is_none() && l.n[0].abs() == 1.0, "{l:?}");
        // Three materials: none.
        let key = |p: &[V3]| p.iter().map(|p| ((p[0] * 2.0).round(), 0.3)).collect();
        assert!(measure(&key, &same, &[1.0; 64]).is_none());
    }
}
