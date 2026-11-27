use super::sig;

#[test]
fn three_grids_moving_apart_give_no_grid_error_but_within_their_noise() {
    use super::grid_error;
    // A 100 N load on a 0.1 mm patch: its displacement read 1.022e-4,
    // 3.262e-4 and 7.736e-4 mm on cells of 4, 2 and 1 mm, and 3.556e-3
    // and 6.336e-3 on 0.25 and 0.125: 3 |x - x2| allowed 2.12e-3.
    assert_eq!(grid_error(7.736e-4, 3.262e-4, Some(1.022e-4), 1e-11), (f64::INFINITY, "three grids moving apart"));
    // Within what the solves resolve, alike.
    assert_eq!(grid_error(1.0 + 3e-9, 1.0 + 2e-9, Some(1.0 + 1e-9), 1e-8).1, "three grids alike within what the solves resolve");
}

/// Two coarse grids that agree and a finest that leaves them: a post's displacement read 1.73608e-3, 1.73650e-3 and 1.73218e-3 mm on
/// cells of 10, 5 and 2.5, against a limit of 1 mm. Whether that is no error at all is the open question.
#[test]
#[ignore]
fn swinging_wider_refuses_a_limit_far_past_it() {
    let e = super::grid_error(1.7321785392548188e-3, 1.7364988796356353e-3, Some(1.7360779735187686e-3), 5.6146797108218886e-8).0;
    assert!(e < 1.0 - 1.7321785392548188e-3, "{e}");
}

#[test]
fn three_grids_swinging_wider_give_no_grid_error() {
    // 100, 99, 140 MPa on cells of 4, 2 and 1 mm: swings of 1 and 41, the newest far the larger. As
    // 100, 101, 142 they move apart; this was 140 +- 123, and `< 270` held.
    assert_eq!(super::grid_error(140.0, 99.0, Some(100.0), 1e-6), (f64::INFINITY, "three grids swinging wider"));
    // A swing growing by no more than three stays inside the error it is given on the next grid: a
    // truss's 0.082629, 0.082596, 0.082630 (swings of 3.3e-5 and 3.4e-5).
    assert_eq!(super::grid_error(0.082630, 0.082596, Some(0.082629), 1e-10).1, "three grids oscillating");
}

/// Three grids of `values` finest first on cells of 1, 2 and 4 mm, each with `read` and `band`.
fn grids(values: [f64; 3], read: Option<crate::ext::ReadError>, band: Option<(f64, f64)>) -> Vec<crate::ext::Measured> {
    resolved(values, read, band, 0.0)
}

/// Its grids, each resolving `share` of the largest value.
fn resolved(values: [f64; 3], read: Option<crate::ext::ReadError>, band: Option<(f64, f64)>, share: f64) -> Vec<crate::ext::Measured> {
    let most = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    values
        .iter()
        .zip([1.0, 2.0, 4.0])
        .map(|(&v, h)| crate::ext::Measured {
            value: Ok(v),
            at: None,
            res: Some(h),
            finest: None,
            warnings: vec![],
            excluded: vec![],
            unseen: None,
            unseen_finer: false,
            resolves: share * most,
            unbounded: Some(crate::ext::NOT_CLASSIFIED.to_string()),
            extent: None,
            read_error: read,
            proven: false,
            held: None,
            band: band.map(|(low, high)| crate::ext::Band { low, high, parts: vec![], basis: String::new(), narrows: vec![] }),
        })
        .collect()
}

#[test]
fn growth_under_a_reads_error_or_a_band_is_not_noise() {
    // A peak stress 100, 103, 107 MPa, its read known to 4% over and 10% under: decided on 95 to 123.3,
    // `< 125` held, though it grows by a third more each halving (a corner the crease rule does not see).
    let read = crate::ext::ReadError { below: Some(0.1), above: Some(0.04), proven: false, how: "" };
    let g = grids([107.0, 103.0, 100.0], Some(read), None);
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 1.0).unwrap().0, f64::INFINITY);
    // A transient peak 300, 304, 312, its time steps' band 307 to 317 on every grid: floor 10, `< 345` held.
    let g = grids([312.0, 304.0, 300.0], None, Some((307.0, 317.0)));
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap().0, f64::INFINITY);
}

#[test]
fn three_grids_converging_keep_roaches_factor_on_their_last_change() {
    // Converging at an order over 1.17, Roache's 1.25 |x - x2| / (2^p - 1) is under the last change: the
    // error is 1.25 times it, not the bare change (a bar's hold, 0.1860121 +- 1.11e-5 for 0.186).
    let g = resolved([0.1860121, 0.1860232, 0.1860565], None, None, 0.0);
    let (e, how) = super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap();
    assert!(how == "three grids converging" && (e - 1.25 * 1.11e-5).abs() < 1e-10, "{e} {how}");
}

#[test]
fn grids_close_but_not_in_order_are_given_three_times_the_larger_change_where_nothing_can_grow_without_bound() {
    // A column's buckling factor on cells off its faces: 7.3492, 7.3536, 7.3532 (truth 7.3524), swinging
    // wider. Its solve found nothing that could grow without bound: three times the larger change, 0.0124.
    let mut g = grids([7.349483341, 7.353622606, 7.353158818], None, None);
    g.iter_mut().for_each(|m| m.unbounded = None);
    let (e, how) = super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap();
    assert!((e - 3.0 * 0.004139265).abs() < 1e-6 && how == "three grids close, not in order", "{e} {how}");
    // The patch load: its load lies inside less than a cell, so its growth is refused.
    let mut g = grids([7.736e-4, 3.262e-4, 1.022e-4], None, None);
    g.iter_mut().for_each(|m| m.unbounded = Some("the load 'tip' is applied over less than a cell".into()));
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap(), (f64::INFINITY, "three grids moving apart"));
    // One grid that says so is enough.
    let mut g = grids([7.349483341, 7.353622606, 7.353158818], None, None);
    g[1].unbounded = None;
    g[2].unbounded = None;
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap().0, f64::INFINITY);
}

#[test]
fn differences_within_the_noise_are_not_an_order() {
    // hs-r1: 3.6e-15 then 1.6e-14 apart, under a floor of 3.2e-4 of the reading: no order, faster or not.
    let g = resolved([10.0 + 3.6e-15 + 1.6e-14, 10.0 + 1.6e-14, 10.0], None, None, 3.2e-5);
    let (e, how) = super::grids_error(&g.iter().collect::<Vec<_>>(), 1.0).unwrap();
    assert_eq!(how, "three grids alike within what the solves resolve", "{e}");
    // One difference in the noise: no order either, judged by the real one.
    let g = resolved([10.0 + 0.02 + 1e-6, 10.0 + 0.02, 10.0], None, None, 3.2e-5);
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 1.0).unwrap().1, "three grids, one difference within what the solves resolve");
}

#[test]
fn grids_alike_within_what_their_solve_resolves_are_not_moving_apart() {
    // A bar pulled to 7.1875 mm read 7.18749875, 7.18749935 and 7.1875 on cells of 2, 1 and 0.5 mm: 9e-8 of
    // it apart, under the 1e-6 stiffness kept in void, and a claim with a margin of seven million was not
    // settled. Within the share a reading's solve does not resolve, the grids are alike.
    let g = grids([7.1875, 7.18749935, 7.18749875], None, None);
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap().1, "three grids moving apart");
    let g = resolved([7.1875, 7.18749935, 7.18749875], None, None, 3.2e-5);
    let (e, how) = super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap();
    assert!((e - 3.2e-5 * 7.1875).abs() < 1e-12 && how != "three grids moving apart", "{e} {how}");
    // A divergence is not: the patch load's 1.022e-4, 3.262e-4, 7.736e-4.
    let g = resolved([7.736e-4, 3.262e-4, 1.022e-4], None, None, 3.2e-5);
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap().0, f64::INFINITY);
    // Nor carried 1000 mm whole: what the solves resolve is of the deflection, not the offset.
    let mut g = grids([1000.0007736, 1000.0003262, 1000.0001022], None, None);
    g.iter_mut().for_each(|m| m.resolves = 3.2e-5 * 7.736e-4);
    assert_eq!(super::grids_error(&g.iter().collect::<Vec<_>>(), 2.0).unwrap().0, f64::INFINITY);
}

#[test]
fn a_difference_within_what_the_solves_resolve_is_no_evidence_of_a_trend() {
    use super::grid_error;
    // A motion's peak 0.06278 and 0.06278 (8e-7 apart) on cells of 4 and 2 mm, 0.06285 on 1: one real
    // difference, judged as two grids are, not a ratio against noise.
    let floor = 3.2e-5 * 0.06285;
    assert_eq!(grid_error(0.06285, 0.0627808, Some(0.0627800), floor), (3.0 * (0.06285 - 0.0627808_f64).abs(), "three grids, one difference within what the solves resolve"));
    // Two real differences growing are still refused: 100, 103, 107.
    assert_eq!(grid_error(107.0, 103.0, Some(100.0), 3.2e-5 * 107.0).0, f64::INFINITY);
}

#[test]
fn an_infinite_value_has_no_grid_error_where_grids_agree_and_any_where_not() {
    // A buckling factor where nothing is compressed is infinite: its
    // grids' differences are no numbers (inf - inf), and NaN failed
    // every claim on it.
    let inf = f64::INFINITY;
    assert_eq!(grid_error(inf, inf, Some(inf), 0.0).0, 0.0);
    assert_eq!(grid_error(inf, inf, None, 0.0).0, 0.0);
    assert_eq!(grid_error(inf, 2000.0, Some(inf), 0.0).0, inf);
    assert_eq!(grid_error(1500.0, inf, Some(inf), 0.0).0, inf);
    assert_eq!(widened(inf, 0.0), (inf, inf));
    assert_eq!(widened(inf, inf), (-inf, inf));
    assert_eq!(widened(5.0, inf), (-inf, inf));
    use super::{Cmp, Comparison, Status, grid_error, point_status, widened};
    let gt = Comparison::exact(Cmp::Gt, 3.0);
    let lt = Comparison::exact(Cmp::Lt, 10.0);
    assert_eq!((point_status(gt, widened(inf, 0.0).0), point_status(lt, widened(inf, 0.0).1)), (Status::Holds, Status::Fails));
    assert_eq!((point_status(gt, widened(inf, inf).0), point_status(gt, widened(inf, inf).1)), (Status::Fails, Status::Holds));
}

#[test]
fn a_read_error_widens_a_finite_value_and_leaves_an_infinite_one() {
    // A safety factor's band with no top (cells not of one material):
    // widened by a fraction of an infinite end, it read NaN, and a
    // claim that holds failed.
    use super::{ReadError, spread};
    let inf = f64::INFINITY;
    let read = Some(ReadError { below: Some(0.5), above: Some(0.1), proven: false, how: "" });
    assert_eq!(spread(inf, 0.4, read), (inf, inf));
    assert_eq!(spread(-inf, 0.4, read), (-inf, -inf));
    assert_eq!(spread(20.0, 1.0, read), (9.0, 23.0));
    let proven = Some(ReadError { below: None, above: Some(0.0), proven: true, how: "" });
    assert_eq!(spread(20.0, 1.0, proven), (19.0, 21.0));
}

#[test]
fn claim_numbers_keep_four_significant_figures() {
    // A displacement of 1.05e-5 mm read "0" to 1e-4 absolute, and so its
    // grid error and its grids' values: "holds at 0, grid error about 0".
    let said: Vec<String> = [1.0549064e-5, 1.7266e-4, 3.7413581, 4211.3, 123456.7, -0.25, 0.0, f64::INFINITY].map(sig).to_vec();
    assert_eq!(said, ["1.055e-5", "1.727e-4", "3.741", "4211", "123500", "-0.25", "0", "unbounded"]);
}
