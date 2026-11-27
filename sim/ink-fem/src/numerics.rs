//! Ordered parallel floating-point sums for the solvers.

use rayon::prelude::*;

/// Sum in fixed 8192-value chunks, then sum their totals in index order.
/// This is the reduction used by the simulation CG and conjugate-heat
/// BiCGSTAB solvers: the grouping is independent of the Rayon pool size
/// and worker scheduling. Each input value must itself be deterministic.
/// An indexed iterator prevents an unordered producer from entering the sum.
/// This is ordinary f64 arithmetic, not an interval bound or an exact sum.
pub fn ordered_sum(values: impl IndexedParallelIterator<Item = f64>) -> f64 {
    let parts: Vec<f64> = values.fold_chunks(8192, || -0.0, |sum, x| sum + x).collect();
    parts.iter().sum()
}

/// Ordered dot product. As with a serial zip, use the shorter input.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    ordered_sum(a.par_iter().zip(b).map(|(p, q)| p * q))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sums_match_the_serial_chunk_order_across_pools() {
        // Cancellation, multiple chunks and a partial final chunk: compare
        // bits with a serial sum of chunk sums, including empty input and
        // signed zero.
        let values: Vec<f64> = (0..16_387)
            .map(|i| match i % 4 {
                0 => 1e16,
                1 => 0.1 * i as f64,
                2 => -1e16,
                _ => -0.3,
            })
            .collect();
        for v in [&values[..], &values[..8192], &[-0.0][..], &[][..]] {
            let want: f64 = v.chunks(8192).map(|c| c.iter().sum::<f64>()).sum();
            for threads in [1, 2, 4] {
                let got = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| ordered_sum(v.par_iter().copied()));
                assert_eq!(got.to_bits(), want.to_bits(), "{threads} threads, {} values", v.len());
            }
        }
    }
}
