//! The tensions of couplings' lines (`(couple ...)` in a relation: a belt's
//! spans, a gear's lines of action, a screw's axis either way), by statics:
//! each passive joint a coupling reaches passes nothing along its motion,
//! so what its bodies need along it is the lines' pull there. A tension
//! that balance leaves free, a mode (a belt's spans tightened together, a
//! gear's two flanks, a screw's push against its pull), is set as the
//! coupling was installed: locked, the belt keeps its length (its spans'
//! tensions less their pretension, weighted by their lengths, sum to 0:
//! Gates' T0 +- Te / 2 on two equal spans), and a line that would push goes
//! slack at 0 (a gear's idle flank, a screw's unloaded direction); held by
//! a constant-force tensioner, its least line carries the pretension.
//! Where the joints' balances cannot all hold (a part driven through
//! something no coupling declares), or leave a tension undetermined beyond
//! that, the tensions are refused, saying which.

use crate::model::eigen;

/// A passive joint's balance along its motion: the coefficient of each
/// line's tension (the torque about its axis, or the force along it, a line
/// of unit tension exerts on what the joint carries), what those bodies need
/// along it (their reaction's part), and that need's rounding.
pub struct Balance {
    pub joint: String,
    pub row: Vec<f64>,
    pub need: f64,
    pub noise: f64,
}

/// A line as the solve reads it: its pretension (N), its length (mm: the
/// weight of its stretch where its mode keeps its length), and whether a
/// constant-force tensioner holds it.
#[derive(Clone, Copy)]
pub struct Line {
    pub t0: f64,
    pub length: f64,
    pub held: bool,
}

/// Each line's tension (N), and the mode each is in (lines that tighten
/// together; None: fixed by the balances alone); Err why none is found.
pub fn tensions(eqs: &[Balance], lines: &[Line], line: &dyn Fn(usize) -> String) -> Result<(Vec<f64>, Vec<Option<usize>>), String> {
    let n = lines.len();
    // A line no balance reads: nothing fixes its coupling's force.
    if let Some(l) = (0..n).find(|&l| eqs.iter().all(|e| e.row[l] == 0.0)) {
        return Err(format!(
            "nothing fixes the force of {}: no part it acts on rides a passive joint (one with no effort-limits), whose balance along its motion \
             would give it",
            line(l)
        ));
    }
    // Rows made unit, so that torque and force balances weigh alike.
    let scale: Vec<f64> = eqs.iter().map(|e| e.row.iter().map(|x| x * x).sum::<f64>().sqrt()).collect();
    let a: Vec<Vec<f64>> = eqs.iter().zip(&scale).map(|(e, s)| e.row.iter().map(|x| if *s > 0.0 { x / s } else { 0.0 }).collect()).collect();
    let b: Vec<f64> = eqs.iter().zip(&scale).map(|(e, s)| if *s > 0.0 { e.need / s } else { e.need }).collect();
    let ata: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| a.iter().map(|r| r[i] * r[j]).sum()).collect()).collect();
    let (vals, vecs) = eigen(&ata);
    let small = 1e-10 * vals.iter().fold(1.0f64, |m, v| m.max(v.abs()));
    let atb: Vec<f64> = (0..n).map(|i| a.iter().zip(&b).map(|(r, y)| r[i] * y).sum()).collect();
    // The least-norm solution, and what the balances leave free.
    let mut x = vec![0.0; n];
    let mut magnitude = vec![0.0; n];
    let mut free: Vec<usize> = Vec::new();
    for k in 0..n {
        if vals[k] > small {
            let c: f64 = (0..n).map(|i| vecs[i][k] * atb[i]).sum::<f64>() / vals[k];
            for i in 0..n {
                let term = c * vecs[i][k];
                x[i] += term;
                magnitude[i] += term.abs();
            }
        } else {
            free.push(k);
        }
    }
    // Every balance met, to rounding: else something no coupling declares
    // drives a joint.
    let worst = eqs
        .iter()
        .enumerate()
        .map(|(e, q)| {
            let s = if scale[e] > 0.0 { scale[e] } else { 1.0 };
            let got: f64 = (0..n).map(|l| a[e][l] * x[l]).sum();
            let size: f64 = b[e].abs() + (0..n).map(|l| (a[e][l] * x[l]).abs()).sum::<f64>();
            ((got - b[e]).abs() - (q.noise / s + 1e-8 * size), e, (got - b[e]) * s)
        })
        .max_by(|p, q| p.0.total_cmp(&q.0));
    if let Some((over, e, by)) = worst
        && over > 0.0
    {
        return Err(format!(
            "the balance of {} along its motion is not met by the couplings acting on what it carries (off by {by:.6}): something no coupling \
             declares drives it (a part rigidly driving it: place it on its driver's frame, or relate the two)",
            eqs[e].joint
        ));
    }
    // The free tensions: lines whose rows in the free directions are
    // parallel move together (a belt's spans), each such set one mode.
    let rows: Vec<Vec<f64>> = (0..n).map(|l| free.iter().map(|&k| vecs[l][k]).collect()).collect();
    let len = |r: &[f64]| r.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut modes: Vec<(Vec<f64>, Vec<(usize, f64)>)> = Vec::new();
    for l in 0..n {
        let size = len(&rows[l]);
        if size < 1e-6 {
            continue;
        }
        let dot = |d: &[f64]| d.iter().zip(&rows[l]).map(|(p, q)| p * q).sum::<f64>();
        match modes.iter_mut().find(|(d, _)| dot(d).abs() >= (1.0 - 1e-6) * size) {
            Some((d, ls)) => {
                let c = dot(d);
                ls.push((l, c));
            }
            None => {
                let d: Vec<f64> = rows[l].iter().map(|x| x / size).collect();
                modes.push((d, vec![(l, size)]));
            }
        }
    }
    let names = |ls: &[(usize, f64)]| ls.iter().map(|&(l, _)| line(l)).collect::<Vec<_>>().join(", ");
    if modes.len() != free.len() {
        let all: Vec<String> = modes.iter().map(|(_, ls)| names(ls)).collect();
        return Err(format!(
            "the parts' balances leave the couplings' tensions undetermined beyond each one's pretension ({} free, in {} sets: {}): couplings in \
             parallel (two belts between one pair of shafts) share by their stiffness, and two actuators driving one coupled train as their \
             control does, neither of which the model gives",
            free.len(),
            modes.len(),
            all.join("; ")
        ));
    }
    let mut mode = vec![None; n];
    for (k, (_, ls)) in modes.iter().enumerate() {
        // A self-balanced tension pulls every line of it, or it is no
        // pretension: one line would push.
        let sign = ls[0].1.signum();
        if ls.iter().any(|&(_, c)| c.signum() != sign) {
            return Err(format!("a tension the balances leave free would pull some of {} and push the others, so no pretension sets it", names(ls)));
        }
        let held = ls.iter().filter(|&&(l, _)| lines[l].held).count();
        if held != 0 && held != ls.len() {
            return Err(format!("{} tighten together, and some are held by a tensioner and some locked: state one for the belt", names(ls)));
        }
        // Raised as installed: held, until its least line carries its
        // pretension; locked, until its lines' stretch, their tensions less
        // their pretension weighted by their lengths, sums to 0, and at
        // least until its least line is slack (none pushes).
        let least = |floor: &dyn Fn(usize) -> f64| ls.iter().map(|&(l, c)| (floor(l) - x[l]) / (c * sign)).fold(f64::NEG_INFINITY, f64::max);
        let lift = if held > 0 {
            least(&|l| lines[l].t0)
        } else {
            let (num, den) = ls.iter().fold((0.0, 0.0), |(p, q), &(l, c)| (p + lines[l].length * (lines[l].t0 - x[l]), q + lines[l].length * c * sign));
            (num / den).max(least(&|_| 0.0))
        };
        for &(l, c) in ls {
            x[l] += lift * c * sign;
            magnitude[l] += (lift * c).abs();
            mode[l] = Some(k);
        }
    }
    // A line no mode moves pulls or holds: it cannot push.
    if let Some(l) = (0..n).find(|&l| x[l] < -1e-8 * magnitude[l]) {
        return Err(format!("the motion needs {} to push, by {:.6} N: a belt, chain or tooth flank only pulls or presses one way", line(l), -x[l]));
    }
    Ok((x, mode))
}

/// How far changes `dneeds` of what the balances need (one vector per
/// independent change: a body's measured moment within its bounds) can move
/// the tensions `tensions` finds: per change, each line's least-norm part
/// exactly (it is linear in the need), and the pretension lift of each set
/// of lines that tighten together, a largest over its lines (held) or over
/// their weighted stretch and their slack (locked), by no more than its
/// largest line's move. Per change, each line's linear move and, for each
/// set, the bound on its lift with the set's lines and their coefficients.
pub fn tension_moves(eqs: &[Balance], lines: &[Line], dneeds: &[Vec<f64>]) -> Vec<(Vec<f64>, Vec<(f64, Vec<(usize, f64)>)>)> {
    let n = lines.len();
    let scale: Vec<f64> = eqs.iter().map(|e| e.row.iter().map(|x| x * x).sum::<f64>().sqrt()).collect();
    let a: Vec<Vec<f64>> = eqs.iter().zip(&scale).map(|(e, s)| e.row.iter().map(|x| if *s > 0.0 { x / s } else { 0.0 }).collect()).collect();
    let ata: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| a.iter().map(|r| r[i] * r[j]).sum()).collect()).collect();
    let (vals, vecs) = eigen(&ata);
    let small = 1e-10 * vals.iter().fold(1.0f64, |m, v| m.max(v.abs()));
    let free: Vec<usize> = (0..n).filter(|&k| vals[k] <= small).collect();
    // The sets, as `tensions` finds them: lines whose rows in the free directions are parallel.
    let rows: Vec<Vec<f64>> = (0..n).map(|l| free.iter().map(|&k| vecs[l][k]).collect()).collect();
    let len = |r: &[f64]| r.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut sets: Vec<(Vec<f64>, Vec<(usize, f64)>)> = Vec::new();
    for l in 0..n {
        let size = len(&rows[l]);
        if size < 1e-6 {
            continue;
        }
        let dot = |d: &[f64]| d.iter().zip(&rows[l]).map(|(p, q)| p * q).sum::<f64>();
        match sets.iter_mut().find(|(d, _)| dot(d).abs() >= (1.0 - 1e-6) * size) {
            Some((d, ls)) => {
                let c = dot(d);
                ls.push((l, c));
            }
            None => sets.push((rows[l].iter().map(|x| x / size).collect(), vec![(l, size)])),
        }
    }
    dneeds
        .iter()
        .map(|dneed| {
            let b: Vec<f64> = dneed.iter().zip(&scale).map(|(y, s)| if *s > 0.0 { y / s } else { *y }).collect();
            let atb: Vec<f64> = (0..n).map(|i| a.iter().zip(&b).map(|(r, y)| r[i] * y).sum()).collect();
            let mut dx = vec![0.0; n];
            for k in (0..n).filter(|&k| vals[k] > small) {
                let c: f64 = (0..n).map(|i| vecs[i][k] * atb[i]).sum::<f64>() / vals[k];
                (0..n).for_each(|i| dx[i] += c * vecs[i][k]);
            }
            let lifts = sets
                .iter()
                .map(|(_, ls)| {
                    let sign = ls[0].1.signum();
                    let line_move = ls.iter().map(|&(l, c)| dx[l].abs() / c.abs()).fold(0.0, f64::max);
                    let held = ls.iter().any(|&(l, _)| lines[l].held);
                    let (num, den) = ls.iter().fold((0.0, 0.0), |(p, q), &(l, c)| (p + lines[l].length * dx[l], q + lines[l].length * c * sign));
                    let lift = if held || den == 0.0 { line_move } else { line_move.max(num.abs() / den.abs()) };
                    (lift, ls.iter().map(|&(l, c)| (l, c * sign)).collect())
                })
                .collect();
            (dx, lifts)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{Balance, Line, tensions};

    #[test]
    fn tension_in_another_line_cannot_hide_a_push() {
        let eqs = [eq(vec![1.0, 0.0], 1e16), eq(vec![0.0, 1.0], -1.0)];
        let err = tensions(&eqs, &belt(0.0, false), &|l| format!("line {l}")).unwrap_err();
        assert!(err.contains("line 1 to push, by 1.000000 N"), "{err}");
    }

    fn eq(row: Vec<f64>, need: f64) -> Balance {
        Balance { joint: "j".into(), row, need, noise: 0.0 }
    }

    fn belt(t0: f64, held: bool) -> [Line; 2] {
        [Line { t0, length: 100.0, held }; 2]
    }

    /// A driven pulley of radius 10 needing 500 N mm, its spans at +-10:
    /// installed at 40 N and locked, Gates' 40 +- 25 (Formulas 5 and 6 at
    /// their Tst); needing 1000, more than locked 40 N carries, the slack
    /// span goes to 0 and the tight one takes 100; held by a tensioner, the
    /// slack span keeps 40.
    #[test]
    fn a_belts_spans_share_its_torque_as_installed() {
        let at = |need: f64, held: bool| tensions(&[eq(vec![10.0, -10.0], need)], &belt(40.0, held), &|l| format!("line {l}")).unwrap();
        let (t, mode) = at(500.0, false);
        assert!((t[0] - 65.0).abs() < 1e-9 && (t[1] - 15.0).abs() < 1e-9 && mode == [Some(0), Some(0)], "{t:?} {mode:?}");
        let (t, _) = at(-500.0, false);
        assert!((t[0] - 15.0).abs() < 1e-9 && (t[1] - 65.0).abs() < 1e-9, "{t:?}");
        let (t, _) = at(1000.0, false);
        assert!((t[0] - 100.0).abs() < 1e-9 && t[1].abs() < 1e-9, "{t:?}");
        let (t, _) = at(500.0, true);
        assert!((t[0] - 90.0).abs() < 1e-9 && (t[1] - 40.0).abs() < 1e-9, "{t:?}");
    }

    #[test]
    fn a_line_nothing_balances_and_a_balance_nothing_meets_are_refused() {
        let free = [Line { t0: 0.0, length: 1.0, held: false }; 2];
        assert!(tensions(&[eq(vec![10.0, 0.0], 5.0)], &free, &|l| format!("line {l}")).unwrap_err().contains("line 1"));
        let e = tensions(&[eq(vec![10.0, -10.0], 500.0), eq(vec![10.0, -10.0], 700.0)], &free, &|l| format!("line {l}")).unwrap_err();
        assert!(e.contains("not met"), "{e}");
    }
}
