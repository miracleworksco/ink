//! Arithmetic on points and directions (`V3`), in floats rounded to
//! nearest: for sampling, meshing and pictures, not for proofs (those use
//! `interval`).

use crate::V3;

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
/// Twice the area vector of a triangle, along its normal (counter-clockwise
/// seen from where it points).
pub fn normal(t: [V3; 3]) -> V3 {
    cross(sub(t[1], t[0]), sub(t[2], t[0]))
}
/// The box (lo, hi) holding nothing: box_union's identity.
pub const EMPTY_BOX: (V3, V3) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
/// The smallest box (lo, hi) holding both boxes (exact: min and max).
pub fn box_union(a: (V3, V3), b: (V3, V3)) -> (V3, V3) {
    (std::array::from_fn(|k| a.0[k].min(b.0[k])), std::array::from_fn(|k| a.1[k].max(b.1[k])))
}
/// `a` scaled to length 1, if it has a length (nonzero and finite): a
/// direction an author gave.
pub fn direction(a: V3) -> Option<V3> {
    let l = norm(a);
    (l > 0.0 && l.is_finite()).then(|| a.map(|c| c / l))
}
/// `a` scaled to length 1 (the zero vector, or one not finite, gives zero).
/// Divided by its largest component first, so its squares neither
/// underflow nor overflow at any length.
pub fn unit(a: V3) -> V3 {
    let m = a.iter().fold(0.0f64, |m, c| m.max(c.abs()));
    if !(m > 0.0 && m.is_finite()) {
        return [0.0; 3];
    }
    let b = scale(a, 1.0 / m);
    scale(b, 1.0 / norm(b))
}

#[cfg(test)]
mod tests {
    #[test]
    fn unit_has_length_one_at_any_length() {
        // A 1e-300 floor on the norm made a tiny vector huge (its squares
        // underflow: [1e130, 0, 0]) and a huge one zero (they overflow).
        assert_eq!(super::unit([1e-170, 0.0, 0.0]), [1.0, 0.0, 0.0]);
        let u = super::unit([1e200, -1e200, 0.0]);
        assert!((u[0] - 0.5f64.sqrt()).abs() < 1e-15 && (u[1] + 0.5f64.sqrt()).abs() < 1e-15 && u[2] == 0.0, "{u:?}");
        assert_eq!(super::unit([0.0; 3]), [0.0; 3]);
    }
}
