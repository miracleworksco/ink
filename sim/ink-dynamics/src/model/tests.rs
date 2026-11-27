use super::pinv;

#[test]
fn a_shared_direction_does_not_hide_a_unique_drive() {
    for unit in [1e-100, 1.0, 1e100] {
        let drive = [vec![unit, 1.0], vec![0.0, 1.0], vec![0.0, 1.0]];
        // The shared drive alone carries the first load. The others
        // supply the remaining 5, split 1:4 by their squared ratings.
        for rating in [1e-100, 4e4, 1e100] {
            let pm = super::shares(&drive, &[rating, 1.0, 2.0], 2).unwrap();
            let f = [7.0 * unit, 12.0];
            for (row, expected) in pm.iter().zip([7.0, 1.0, 4.0]) {
                assert!((super::dot_n(row, &f) - expected).abs() < 1e-12, "{unit}, {rating}: {pm:?}");
            }
        }
    }
}

#[test]
fn a_fixed_drive_supplies_its_only_effort_when_the_other_can_balance_it() {
    let drives = [vec![1.0], vec![1.0]];
    for (spans, fixed, efforts) in [([0.0, 10.0], [3.0, 0.0], [3.0, 7.0]), ([0.0, 0.0], [3.0, 5.0], [4.0, 6.0])] {
        let (pm, offset) = super::drives::effort_split(&drives, &spans, &fixed, 1);
        for i in 0..2 {
            assert!((pm[i][0] * 10.0 + offset[i] - efforts[i]).abs() < 1e-12);
        }
    }
}

#[test]
fn pinv_inverts_a_regular_matrix_and_drops_a_null_direction() {
    let a = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
    let p = pinv(&a);
    for i in 0..2 {
        for j in 0..2 {
            let x: f64 = (0..2).map(|k| a[i][k] * p[k][j]).sum();
            assert!((x - if i == j { 1.0 } else { 0.0 }).abs() < 1e-14, "{x}");
        }
    }
    // [1 1; 1 1] = 2 u u^T, u = (1, 1)/sqrt 2: its pseudo-inverse is u u^T / 2.
    let p = pinv(&[vec![1.0, 1.0], vec![1.0, 1.0]]);
    assert!(p.iter().flatten().all(|x| (x - 0.25).abs() < 1e-15), "{p:?}");
}
