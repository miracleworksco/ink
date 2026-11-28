//! A mechanism's natural frequencies at a pose, and how it moves in each
//! (its modes): small oscillations about the pose, at rest, undamped, as a
//! multibody model is linearised (Adams/Linear, Simscape's linearize). The
//! rigid mechanism moves by its coordinates (Model); a coupling line that
//! states a stiffness (a belt span's EA / L) is a spring, so each joint one
//! of its ends rides gets a coordinate of its own, a small turn or slide
//! beside what the coordinates move it by (and so does each joint a rigid
//! line reaches from one of those: it stays a constraint), and a joint that
//! states one, `(stiffness K)`, holds its value (a stepper's holding
//! stiffness). Mass: each rigid body's inertia along all of them (its
//! partial velocities, as Kane's sums give the efforts'). Stiffness: each
//! spring's stretch per unit of each coordinate, squared, times its
//! stiffness. What moves no body (a coordinate whose every body rides
//! released joints) is dropped, and what has stiffness but no mass (an
//! idler given no inertia) is condensed out statically (Guyan); then K phi
//! = omega^2 M phi, its shapes by Jacobi's method on the mass-normalised
//! stiffness. Whether a spring holds a direction, and how much, is read
//! from each spring alone (its stretch along it), never against the
//! stiffest: a soft mode beside a stiff one keeps its own frequency.

use super::*;
use ink_lang::Reload;

/// A placement's partial along coordinate k (dR, dT), and a body's (its
/// centre's velocity and its turning per unit of k).
type Partial = (usize, M3, V3);
type Velocity = (usize, V3, V3);

/// What a mode's stiffness comes from: a coupling's line or a joint, its
/// stiffness, in N/mm (a line, a slide) or N mm/rad (a turn), and where it
/// is stated (file:line: the relation, the joint's clause).
pub struct Spring {
    pub what: String,
    pub stiffness: f64,
    pub turn: bool,
    pub site: String,
}

/// A joint's part in a mode: its instance path, its motion (deg a turn, mm
/// a slide; the mode scaled so the joint whose bodies carry most of its
/// kinetic energy moves 1), and the share of that energy its bodies carry
/// (those it is the nearest joint of). A body no joint moves (placed by a
/// coordinate) is its own, its path in `joint` and no motion.
pub struct Moves {
    pub joint: String,
    pub turn: bool,
    pub motion: Option<f64>,
    pub energy: f64,
}

/// One mode: its frequency (Hz; 0 where nothing holds it) and what moves.
pub struct Mode {
    pub frequency: f64,
    pub moves: Vec<Moves>,
}

/// A mechanism's modes at moment t, lowest first, and what they rest on.
pub struct Modes {
    pub t: f64,
    /// Each coordinate's value at the pose.
    pub pose: Vec<(String, f64)>,
    pub modes: Vec<Mode>,
    pub springs: Vec<Spring>,
    /// The joints given a coordinate of their own, and the lines with no
    /// stiffness that tie some of them (rigid).
    pub released: Vec<String>,
    pub rigid: Vec<String>,
    /// The model the masses were measured in (its carried and deforming
    /// bodies are left out of the modes).
    pub model: Model,
}

/// The modes of `target` in the file `reload` loads, at moment t (its
/// coordinates' values there; their speeds are not read).
pub fn modes(reload: &Reload, target: &str, t: f64) -> Result<Modes, InkError> {
    let doc = reload.with_coordinates(target)?;
    if doc.assembly(target).is_none() {
        return Err(InkError::msg(format!("'{target}' is not a top-level assembly of the file: modes are of the mechanism an assembly's instances pass arguments to")));
    }
    let (ctx, unsolved) = doc.ctx_at(&reload.overrides, t).map_err(InkError::msg)?;
    let mut v = ModelView::new(doc, ctx);
    v.overrides = reload.overrides.clone();
    v.unsolved = unsolved;
    let items = v.doc.items(target).map_err(InkError::msg)?;
    let relations: Vec<crate::coupling::Related> = crate::coupling::related(&mut v.doc, target).into_iter().filter(|r| r.couple.is_some()).collect();
    let stiff = |r: &crate::coupling::Related, k: bool| r.couple.as_ref().unwrap().lines.iter().any(|l| l.stiffness.is_some() == k);
    if !relations.iter().any(|r| stiff(r, true)) && !items.iter().flat_map(|it| &it.links).any(|l| l.joint.as_ref().is_some_and(|j| crate::joint::stiffness(j).is_some())) {
        return Err(InkError {
            hint: Some(
                "a belt's spans are springs given its tensile stiffness, (belt-span SPAN T0 E (tensile-stiffness EA)); a joint is held by (stiffness K) \
                 (a catalog stepper's drive states its holding stiffness)"
                    .into(),
            ),
            ..InkError::msg(format!("nothing under {target} states a stiffness (no coupling line, no joint), so it has no natural frequency"))
        });
    }
    // The joints a part's motion rests on: the one it rides and every one
    // supporting it, up to the target (a pulley's turn and the carriage
    // carrying it).
    let jointed: HashSet<&str> = items.iter().flat_map(|it| &it.links).filter(|l| l.joint.is_some()).map(|l| l.path.as_str()).collect();
    let resting = |path: &str| super::springs::resting(&items, &jointed, path, target);
    // The joints released from `seed`: each one, then those a rigid line
    // reaches from one of them, until none is added.
    let close = |seed: Vec<String>| -> Vec<String> {
        let mut released = seed;
        let release = |released: &mut Vec<String>, r: &crate::coupling::Related| {
            for e in resting(&r.a).into_iter().chain(resting(&r.b)) {
                if !released.contains(&e) {
                    released.push(e);
                }
            }
        };
        loop {
            let (sets, _) = carried(&items, &released);
            let before = released.len();
            for r in relations.iter().filter(|r| stiff(r, false) && sets.iter().any(|s| s.contains(&r.a) || s.contains(&r.b))) {
                release(&mut released, r);
            }
            if released.len() == before {
                return released;
            }
        }
    };
    // Released: each joint a stiff line's end rests on, and each one its
    // own stiffness holds that no coordinate drives (a constant drive, or
    // one reading only sizes). One a coordinate drives keeps it as its
    // motion, with whatever its drive's expression ties to it (a stepper
    // geared to an arm by their drives). Which arguments are coordinates
    // the model finds, so every joint stating a stiffness is among those
    // whose bodies it measures (as their reactions' are), moving along a
    // coordinate or not.
    let mut ends: Vec<String> = Vec::new();
    for r in relations.iter().filter(|r| stiff(r, true)) {
        for e in resting(&r.a).into_iter().chain(resting(&r.b)) {
            if !ends.contains(&e) {
                ends.push(e);
            }
        }
    }
    let held: Vec<String> =
        items.iter().flat_map(|it| &it.links).filter(|l| l.joint.as_ref().is_some_and(|j| crate::joint::stiffness(j).is_some())).map(|l| l.path.clone()).collect();
    let mut seed = ends.clone();
    for h in &held {
        if !seed.contains(h) {
            seed.push(h.clone());
        }
    }
    let measured = close(seed);
    let (model, drives) = Model::compile_in(&mut v, target, [0.0; 3], &measured, false, false, false, |v, eps, joints| {
        joints
            .into_iter()
            .map(|(path, j)| {
                let along: Vec<_> = eps.iter().map(|&k| v.doc.derivative(j.drive, Var::Param(k))).collect();
                (path, j, along)
            })
            .collect::<Vec<_>>()
    })?;
    let ctx = model.ctx(t).map_err(InkError::msg)?;
    let nq = model.coords.len();
    let at = |d: &Doc, n: NodeId| Num::of(&d.g, n).at(&ctx);
    // Every joint under the target: its path, kind, frame in the world, its
    // drive's rate along each coordinate (radians a turn) and stiffness.
    struct J {
        path: String,
        turn: bool,
        axis: V3,
        origin: V3,
        along: Vec<f64>,
        stiffness: Option<(f64, String)>,
    }
    let mut joints: Vec<J> = Vec::new();
    for (path, j, along) in drives {
        let frame = path_frame(&mut v, &items, target, &path).expect("a joint placed under the target");
        let f: Vec<f64> = frame.iter().map(|&n| at(&v.doc, n)).collect();
        let unit = if j.turn { std::f64::consts::PI / 180.0 } else { 1.0 };
        let along = along.into_iter().map(|n| unit * at(&v.doc, n)).collect();
        joints.push(J {
            path,
            turn: j.turn,
            axis: [f[j.axis], f[3 + j.axis], f[6 + j.axis]],
            origin: [f[9], f[10], f[11]],
            along,
            stiffness: crate::joint::stiffness(&j).map(|s| (at(&v.doc, s.k), s.site.clone())),
        });
    }
    // Each coupling line, in the model's order: its point, direction, pitch
    // and stiffness.
    struct L {
        name: String,
        p: V3,
        u: V3,
        h: f64,
        k: Option<f64>,
        tension: f64,
        site: String,
    }
    let mut lines: Vec<L> = Vec::new();
    for (c, r) in relations.iter().enumerate() {
        let name = &model.couplings[c].name;
        let couple = r.couple.as_ref().unwrap();
        for (i, l) in couple.lines.iter().enumerate() {
            let u: V3 = l.dir.map(|n| at(&v.doc, n));
            let len = norm(u);
            let site = couple.site.clone();
            lines.push(L {
                name: format!("line {} of {name}", i + 1),
                p: l.at.map(|n| at(&v.doc, n)),
                u: u.map(|x| x / len),
                h: at(&v.doc, l.pitch),
                k: l.stiffness.map(|n| at(&v.doc, n)),
                tension: at(&v.doc, couple.pretension),
                site,
            });
        }
    }
    if let Some(l) = lines.iter().find(|l| l.k.is_some_and(|k| !(k >= 0.0 && k.is_finite()))) {
        return Err(InkError::msg(format!("the stiffness of {} is {}: a stiffness is at least 0 and finite", l.name, l.k.unwrap())));
    }
    if let Some(j) = joints.iter().find(|j| j.stiffness.as_ref().is_some_and(|k| !(k.0 >= 0.0 && k.0.is_finite()))) {
        return Err(InkError::msg(format!("the stiffness of {} is {}: a stiffness is at least 0 and finite", j.path, j.stiffness.as_ref().unwrap().0)));
    }
    // A stiff joint keeps its coordinate where it is the only one on it;
    // two stiffnesses on one coordinate are two actuators (a gripper's two
    // jaws on one width), each its own spring, so each is released.
    let sprung: Vec<&J> = joints.iter().filter(|j| j.stiffness.is_some()).collect();
    let shares = |j: &J| sprung.iter().any(|o| o.path != j.path && (0..nq).any(|k| j.along[k] != 0.0 && o.along[k] != 0.0));
    let free: Vec<String> = sprung.iter().filter(|j| (j.along.iter().all(|x| *x == 0.0) || shares(j)) && !ends.contains(&j.path)).map(|j| j.path.clone()).collect();
    let mut seed = ends.clone();
    seed.extend(free);
    let released = close(seed);
    let (sets, _) = carried(&items, &released);
    let nr = released.len();
    let n = nq + nr;
    let rj: Vec<&J> = released.iter().map(|p| joints.iter().find(|j| j.path == *p).expect("a released joint")).collect();
    // A rigid motion's partial along a released joint (a unit turn about its
    // axis, a unit slide along it), at a body placed R p + T: dR, dT.
    let unit_motion = |j: &J, r: &M3, x: V3| -> (M3, V3) {
        if j.turn {
            let w: M3 = [[0.0, -j.axis[2], j.axis[1]], [j.axis[2], 0.0, -j.axis[0]], [-j.axis[1], j.axis[0], 0.0]];
            let arm: V3 = std::array::from_fn(|k| x[k] - j.origin[k]);
            (mul(&w, r), mv(&w, &arm))
        } else {
            ([[0.0; 3]; 3], j.axis)
        }
    };
    // M: each rigid body's, along the coordinates and the released joints
    // carrying it; each body's partials kept for its share of a mode's
    // energy (its velocity per unit of each coordinate, and its turning).
    let mut m = vec![vec![0.0; n]; n];
    // Each one's heft, the size of what its inertia sums: each body's mass
    // times its speed squared, and its turning squared times its inertia's
    // diagonal (as a simulation weighs a released joint, simulate::simulation).
    let mut heft = vec![0.0; n];
    let mut partials: Vec<(f64, M3, Vec<Velocity>)> = Vec::new();
    // A body a coordinate moves otherwise than through the released joints
    // carrying it, where a released joint's drive reads that coordinate too:
    // its tie to them is its drive's expression alone (an arm turned a fifth
    // of a belt-driven shaft by their drives), which releasing them cuts.
    let mut cut: Vec<(String, String)> = Vec::new();
    for b in &model.bodies {
        let (r, t) = (rot(&b.at.r.iter().map(|e| e.at(&ctx)).collect::<Vec<_>>()), b.at.p.each_ref().map(|e| e.at(&ctx)));
        let mut dq: Vec<Partial> = b.at.dq.iter().map(|(k, d)| (*k, split(&d.iter().map(|e| e.at(&ctx)).collect::<Vec<_>>()))).map(|(k, (dr, dp))| (k, dr, dp)).collect();
        for (i, j) in rj.iter().enumerate() {
            if sets[i].contains(&b.label) {
                let (dr, dp) = unit_motion(j, &r, t);
                dq.push((nq + i, dr, dp));
            }
        }
        for &(k, dr, dp) in dq.iter().filter(|x| x.0 < nq) {
            let (mut er, mut ep, mut size) = (dr, dp, norm(dp) + dr.iter().flatten().map(|x| x.abs()).sum::<f64>());
            for &(i, ur, up) in dq.iter().filter(|x| x.0 >= nq) {
                let a = rj[i - nq].along[k];
                (0..3).for_each(|x| (0..3).for_each(|y| er[x][y] -= a * ur[x][y]));
                (0..3).for_each(|x| ep[x] -= a * up[x]);
                size += a.abs() * (norm(up) + ur.iter().flatten().map(|x| x.abs()).sum::<f64>());
            }
            let off = norm(ep) + er.iter().flatten().map(|x| x.abs()).sum::<f64>();
            if off > 1e-9 * size && rj.iter().any(|j| j.along[k] != 0.0) {
                cut.push((b.label.clone(), model.coords[k].name.clone()));
            }
        }
        let rt = tr(&r);
        let iw = mul(&mul(&r, &b.inertia), &rt);
        let vel: Vec<Velocity> = dq.iter().map(|(k, dr, dp)| (*k, std::array::from_fn(|i| mv(dr, &b.center)[i] + dp[i]), vee(&mul(dr, &rt)))).collect();
        let diagonal = (0..3).map(|i| iw[i][i].abs()).sum::<f64>();
        for (k, jv, jw) in &vel {
            heft[*k] += b.mass * dot(*jv, *jv) + dot(*jw, *jw) * diagonal;
            for (l, kv, kw) in &vel {
                m[*k][*l] += b.mass * dot(*jv, *kv) + dot(*jw, mv(&iw, kw));
            }
        }
        partials.push((b.mass, iw, vel));
    }
    // One whose inertia is within rounding of its heft carries none, as a
    // simulation reads it: an idler given none, its measured centre a hair
    // off its axis, has m r^2 there and m r^2 less in its own inertia (its
    // spin-inertia, about the axis), which sum to rounding, of either sign.
    // Left, a crumb above 0 was scaled to a unit inertia (below) and rang at
    // 4e27 Hz; now it is condensed out as one at 0 is.
    let none: Vec<usize> = (0..n).filter(|&i| m[i][i] <= ROUNDING * heft[i]).collect();
    for i in none {
        (0..n).for_each(|j| (m[i][j], m[j][i]) = (0.0, 0.0));
    }
    if !cut.is_empty() {
        let names: Vec<String> = cut.iter().take(3).map(|(b, q)| format!("{b} (with {q})")).collect();
        let more = if cut.len() > 3 { format!(" and {} more", cut.len() - 3) } else { String::new() };
        return Err(InkError {
            hint: Some("relate it to what drives it (gear-mesh, belt-drive, lead-screw: its line then ties it, rigid or a spring), or place it on its driver's frame".into()),
            ..InkError::msg(format!(
                "{}{more} move with a coordinate that also drives joints released for the springs, but not through them: tied to them only by the drives' \
                 expressions, which a mode cannot follow (released: {})",
                names.join(", "),
                released.join(", ")
            ))
        });
    }
    // A coupling's numbers as an effort checks them (Coupled::unsound): a
    // second number above 1, where its back efficiency goes, may be a
    // tensile stiffness meant (EA), which is tagged.
    if let Some((c, why)) = model.couplings.iter().find_map(|c| Some((c, c.unsound(&ctx)?))) {
        let (forward, back) = c.efficiency.as_ref().map_or((0.0, 0.0), |(e, b)| (e.at(&ctx), b.at(&ctx)));
        let hint = (back > 1.0 && back != forward)
            .then(|| "if a tensile stiffness was meant, it is (tensile-stiffness EA) among the optional arguments, not a number after the efficiency".to_string());
        return Err(InkError { hint, ..InkError::msg(why) });
    }
    // Each line's stretch per unit of each coordinate: its ends' points on
    // it drawing apart along it (and turning apart about it, times its
    // pitch), along the coordinates by their placements' exact partials,
    // along a released joint by its unit motion. A rigid line's first is 0
    // (the effort's coupling check: refused past 1e-3 of the terms); a
    // coordinate whose every body rides released joints moves nothing and
    // stretches nothing, to rounding (its ends' partials are their joints'
    // motions, by the chain rule).
    let vals = |ns: &[Num]| -> Vec<f64> { ns.iter().map(|e| e.at(&ctx)).collect() };
    let mut stretch: Vec<Vec<f64>> = Vec::new();
    let mut li = 0;
    for c in &model.couplings {
        let ends: Vec<(M3, V3, Vec<Partial>)> = c
            .ends
            .iter()
            .map(|e| {
                let (r, o) = split(&vals(&e.on));
                (r, o, e.dq.iter().map(|(k, d)| (*k, split(&vals(d)))).map(|(k, (dr, dp))| (k, dr, dp)).collect())
            })
            .collect();
        for _ in &c.lines {
            let l = &lines[li];
            li += 1;
            let mut g = vec![0.0; n];
            let mut size = vec![0.0; n];
            // End e's velocity at P and turning per unit of a motion (dR, dT)
            // of its placement (R, o).
            let moved = |r: &M3, o: &V3, dr: &M3, dp: &V3| -> (V3, V3, f64) {
                let w = mul(dr, &tr(r));
                let arm: V3 = std::array::from_fn(|j| l.p[j] - o[j]);
                let turn = mv(&w, &arm);
                (std::array::from_fn(|j| turn[j] + dp[j]), vee(&w), norm(turn) + norm(*dp))
            };
            for (e, sign) in [(0, 1.0), (1, -1.0)] {
                let (r, o, dq) = &ends[e];
                let mut each: Vec<Partial> = dq.clone();
                for (i, j) in rj.iter().enumerate() {
                    if sets[i].contains(&c.ends[e].path) {
                        let (dr, dp) = unit_motion(j, r, *o);
                        each.push((nq + i, dr, dp));
                    }
                }
                for (k, dr, dp) in &each {
                    let (vel, w, s) = moved(r, o, dr, dp);
                    g[*k] += sign * (dot(l.u, vel) + l.h * dot(l.u, w));
                    size[*k] += s + l.h.abs() * norm(w);
                }
            }
            // A spring may stretch as the coordinates move (a carriage on a
            // spring to ground); a rigid line may not (else it does not
            // match the joints' motion).
            if let Some(k) = (0..nq).find(|&k| l.k.is_none() && g[k].abs() > 1e-3 * size[k]) {
                return Err(InkError::msg(format!(
                    "{} separates its parts at {:.6} mm per unit of {}: the coupling does not match the joints' motion, so its stretch is not the \
                     joints' (ink effort names the same)",
                    l.name, g[k], model.coords[k].name
                )));
            }
            stretch.push(g);
        }
    }
    // K: each stiff line's stretch squared times its stiffness; each joint
    // stating one, its value's change (what the coordinates move it by,
    // plus its own where released) squared times it. A rigid line touching
    // a released joint is a constraint: its stretch is 0.
    // Each spring's stretch per unit of each coordinate and its stiffness,
    // kept: whether it holds a mode, and how much, is its own.
    let mut terms: Vec<(Vec<f64>, f64)> = Vec::new();
    let mut springs = Vec::new();
    let mut constraints: Vec<Vec<f64>> = Vec::new();
    let mut rigid = Vec::new();
    let (tangent, unpaired) = springs::line_tangents(&stretch, lines.iter().map(|l| (l.k.unwrap_or(0.0), l.tension)));
    for (i, (l, g)) in lines.iter().zip(&stretch).enumerate() {
        match l.k {
            Some(s) => {
                terms.push((g.clone(), tangent[i]));
                springs.push(Spring { what: l.name.clone(), stiffness: s, turn: false, site: l.site.clone() });
            }
            None if g[nq..].iter().any(|x| *x != 0.0) => {
                constraints.push(g.clone());
                rigid.push(l.name.clone());
            }
            None => {}
        }
    }
    if let Some(j) = unpaired {
        return Err(InkError::msg(format!(
            "{} has zero pretension and no matching opposed spring: a common stiffness for both directions of displacement could not be established",
            lines[j].name
        )));
    }
    let value = |j: &J| -> Vec<f64> {
        let mut h = vec![0.0; n];
        h[..nq].copy_from_slice(&j.along);
        if let Some(i) = released.iter().position(|p| *p == j.path) {
            h[nq + i] = 1.0;
        }
        h
    };
    for j in &joints {
        if let Some((s, site)) = &j.stiffness {
            terms.push((value(j), *s));
            springs.push(Spring { what: format!("joint {}", j.path), stiffness: *s, turn: j.turn, site: site.clone() });
        }
    }
    // Scaled so each coordinate's own inertia is 1 (mixed units: deg, mm,
    // radians): what rounds to 0 is then 0 alike in each.
    let d: Vec<f64> = (0..n).map(|i| if m[i][i] > 0.0 { 1.0 / m[i][i].sqrt() } else { 1.0 }).collect();
    let scale = |a: &[Vec<f64>]| -> Vec<Vec<f64>> { (0..n).map(|i| (0..n).map(|j| d[i] * a[i][j] * d[j]).collect()).collect() };
    // The springs in the scaled coordinates (z = D y: a stretch g z is (g D) y).
    let terms: Vec<(Vec<f64>, f64)> = terms.into_iter().map(|(g, s)| ((0..n).map(|i| g[i] * d[i]).collect(), s)).collect();
    let mut ks = vec![vec![0.0; n]; n];
    for (g, s) in &terms {
        (0..n).for_each(|a| (0..n).for_each(|b| ks[a][b] += s * g[a] * g[b]));
    }
    let ms = scale(&m);
    // A motion y's stiffness, each spring's stretch along it squared times
    // it, and whether any spring holds it (stretches): each spring judged
    // on its own scale, so a stiff one elsewhere hides no soft one.
    let stiffness =
        |y: &[f64]| -> (f64, bool) { terms.iter().fold((0.0, false), |(k, held), (g, s)| (k + s * dot_n(g, y).powi(2), held || (*s > 0.0 && springs::stretches(g, y)))) };
    // What the constraints leave free (their null space).
    let basis: Vec<Vec<f64>> = if constraints.is_empty() {
        (0..n).map(|i| (0..n).map(|j| f64::from(u8::from(i == j))).collect()).collect()
    } else {
        // z = D y: a constraint g z = 0 is (g D) y = 0; each made unit, so
        // a line between heavy parts is as much a constraint as one between
        // light ones.
        let c: Vec<Vec<f64>> = constraints
            .iter()
            .map(|g| {
                let r: Vec<f64> = (0..n).map(|i| g[i] * d[i]).collect();
                let size = r.iter().map(|x| x * x).sum::<f64>().sqrt();
                r.iter().map(|x| x / size).collect()
            })
            .collect();
        let ctc: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| c.iter().map(|r| r[i] * r[j]).sum()).collect()).collect();
        let (e, vecs) = eigen(&ctc);
        let top = e.iter().fold(0.0f64, |a, x| a.max(x.abs()));
        (0..n).filter(|&i| e[i] <= 1e-12 * top).map(|i| (0..n).map(|r| vecs[r][i]).collect()).collect()
    };
    let project = |a: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let nb = basis.len();
        (0..nb).map(|p| (0..nb).map(|q| (0..n).map(|i| basis[p][i] * (0..n).map(|j| a[i][j] * basis[q][j]).sum::<f64>()).sum()).collect()).collect()
    };
    let (m1, k1) = (project(&ms), project(&ks));
    let nb = basis.len();
    // Massive directions, mass-normalised, and massless ones.
    let (e, vecs) = eigen(&m1);
    let top = e.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    let heavy: Vec<usize> = (0..nb).filter(|&i| e[i] > 1e-10 * top).collect();
    let light: Vec<usize> = (0..nb).filter(|&i| e[i] <= 1e-10 * top).collect();
    let col = |i: usize, s: f64| -> Vec<f64> { (0..nb).map(|r| s * vecs[r][i]).collect() };
    let w: Vec<Vec<f64>> = heavy.iter().map(|&i| col(i, 1.0 / e[i].sqrt())).collect();
    let y: Vec<Vec<f64>> = light.iter().map(|&i| col(i, 1.0)).collect();
    let form = |a: &[Vec<f64>], b: &[Vec<f64>]| -> Vec<Vec<f64>> {
        a.iter().map(|p| b.iter().map(|q| (0..nb).map(|i| p[i] * (0..nb).map(|j| k1[i][j] * q[j]).sum::<f64>()).sum()).collect()).collect()
    };
    let kww = form(&w, &w);
    let kwy = form(&w, &y);
    let kyy = form(&y, &y);
    // A direction of the constraints' basis in the scaled coordinates.
    let full = |zb: &[f64]| -> Vec<f64> { (0..n).map(|r| (0..nb).map(|p| basis[p][r] * zb[p]).sum()).collect() };
    // Condensed: the massless directions where the springs leave them,
    // K_yy^+ dropping each that no spring holds, its stiffness its springs'.
    let kyy_inv = {
        let (_, v) = eigen(&kyy);
        let ny = kyy.len();
        let keep: Vec<f64> = (0..ny)
            .map(|q| {
                let along: Vec<f64> = (0..nb).map(|r| (0..ny).map(|c| v[c][q] * y[c][r]).sum()).collect();
                match stiffness(&full(&along)) {
                    (k, true) if k > 0.0 => 1.0 / k,
                    _ => 0.0,
                }
            })
            .collect();
        (0..ny).map(|i| (0..ny).map(|j| (0..ny).map(|q| v[i][q] * keep[q] * v[j][q]).sum()).collect()).collect::<Vec<Vec<f64>>>()
    };
    let nw = w.len();
    let ny = y.len();
    // y = -K_yy^+ K_yw w.
    let follow: Vec<Vec<f64>> = (0..ny).map(|a| (0..nw).map(|b| -(0..ny).map(|c| kyy_inv[a][c] * kwy[b][c]).sum::<f64>()).collect()).collect();
    let kc: Vec<Vec<f64>> = (0..nw).map(|a| (0..nw).map(|b| kww[a][b] + (0..ny).map(|c| kwy[a][c] * follow[c][b]).sum::<f64>()).collect()).collect();
    let (_, phi) = eigen(&kc);
    // Each mode in the scaled coordinates: w's direction, the massless ones
    // following, back through the constraints' basis; its stiffness over
    // its inertia (Rayleigh's quotient: each spring's own stretch, so a
    // soft mode is not read to a stiff one's rounding), and whether any
    // spring holds it. One no spring holds is free: 0, whatever its
    // rounding.
    let shapes: Vec<Vec<f64>> = (0..nw)
        .map(|i| {
            let wi: Vec<f64> = (0..nw).map(|a| phi[a][i]).collect();
            let mut zb = vec![0.0; nb];
            for (a, wv) in w.iter().enumerate() {
                (0..nb).for_each(|r| zb[r] += wi[a] * wv[r]);
            }
            for (c, yv) in y.iter().enumerate() {
                let s: f64 = (0..nw).map(|b| follow[c][b] * wi[b]).sum();
                (0..nb).for_each(|r| zb[r] += s * yv[r]);
            }
            full(&zb)
        })
        .collect();
    // Each read again among those below it, stiffest first (Rayleigh-Ritz in their span, its stiffness
    // each spring's own): the eigensolve above reads each shape to its stiffest spring's rounding, so
    // a mode far softer than that leans toward the others below it by that much (a soft span beside
    // one 1e15 times stiffer, toward the free mode: 2 parts in 10,000 low). Among those below the
    // stiffest, the stiff spring barely stretches, and the next eigensolve reads them to their own
    // scale; and so on down.
    let quad = |a: &[f64], b: &[f64]| -> f64 { (0..n).map(|i| a[i] * (0..n).map(|j| ms[i][j] * b[j]).sum::<f64>()).sum() };
    let ritz = |vs: &[Vec<f64>]| -> Vec<Vec<f64>> {
        let k = vs.len();
        let (mu, u) = eigen(&(0..k).map(|a| (0..k).map(|b| quad(&vs[a], &vs[b])).collect()).collect::<Vec<Vec<f64>>>());
        let unit: Vec<Vec<f64>> = (0..k).map(|c| (0..n).map(|i| (0..k).map(|a| u[a][c] * vs[a][i]).sum::<f64>() / mu[c].sqrt()).collect()).collect();
        let along = |a: &[f64], b: &[f64]| -> f64 { terms.iter().map(|(g, s)| s * dot_n(g, a) * dot_n(g, b)).sum() };
        let (lam, q) = eigen(&(0..k).map(|a| (0..k).map(|b| along(&unit[a], &unit[b])).collect()).collect::<Vec<Vec<f64>>>());
        let mut by: Vec<usize> = (0..k).collect();
        by.sort_by(|&a, &b| lam[a].total_cmp(&lam[b]));
        by.iter().map(|&c| (0..n).map(|i| (0..k).map(|a| q[a][c] * unit[a][i]).sum()).collect()).collect()
    };
    let mut below = shapes;
    let mut shapes: Vec<Vec<f64>> = Vec::new();
    while !below.is_empty() {
        below = ritz(&below);
        shapes.extend(below.pop());
    }
    // Each one's stiffness over its inertia (Rayleigh's quotient, each spring's own stretch).
    let shapes: Vec<(Vec<f64>, f64)> = shapes
        .into_iter()
        .map(|zs| {
            let inertia = quad(&zs, &zs);
            let lam = match stiffness(&zs) {
                (k, _) if inertia > 0.0 => k / inertia,
                _ => 0.0,
            };
            (zs, lam)
        })
        .collect();
    let mut order: Vec<usize> = (0..nw).collect();
    order.sort_by(|&a, &b| shapes[a].1.total_cmp(&shapes[b].1));
    // How many of them no spring holds: as many as the springs' stretches leave free among the
    // shapes, which is a count of which springs reach which motions, whatever their stiffnesses
    // (each spring's row and each shape made unit). The shapes themselves cannot say: a free one
    // reads its springs' rounding, not 0. The free ones are the lowest.
    let unit = |v: &[f64]| -> Vec<f64> {
        let size = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        v.iter().map(|x| if size > 0.0 { x / size } else { 0.0 }).collect()
    };
    let units: Vec<Vec<f64>> = shapes.iter().map(|(zs, _)| unit(zs)).collect();
    let mut rows: Vec<Vec<f64>> = terms.iter().filter(|(_, s)| *s > 0.0).map(|(g, _)| unit(g)).map(|g| units.iter().map(|z| dot_n(&g, z)).collect()).collect();
    let mut held = 0;
    for col in 0..nw {
        let Some(p) = (held..rows.len()).max_by(|&a, &b| rows[a][col].abs().total_cmp(&rows[b][col].abs())).filter(|&p| rows[p][col].abs() > springs::HELD) else { continue };
        rows.swap(held, p);
        let pivot = rows[held].clone();
        for r in rows.iter_mut().skip(held + 1) {
            let f = r[col] / pivot[col];
            r.iter_mut().zip(&pivot).for_each(|(x, y)| *x -= f * y);
        }
        held += 1;
    }
    let shapes: Vec<(Vec<f64>, f64)> = {
        let free: Vec<usize> = order.iter().copied().take(nw - held.min(nw)).collect();
        shapes.into_iter().enumerate().map(|(i, (zs, lam))| (zs, if free.contains(&i) { 0.0 } else { lam })).collect()
    };
    // Each body's nearest joint, whose share of a mode's energy it carries.
    let group: Vec<String> = model.bodies.iter().map(|b| rides(&items, &b.label, target).0.map_or_else(|| b.label.clone(), |(p, _)| p)).collect();
    let mut out = Vec::new();
    for &i in &order {
        // The mode in the coordinates, unscaled.
        let (zs, lam) = &shapes[i];
        let z: Vec<f64> = (0..n).map(|r| d[r] * zs[r]).collect();
        // Energy by body, summed by nearest joint.
        let mut shares: Vec<(String, f64)> = Vec::new();
        for ((mass, iw, vel), g) in partials.iter().zip(&group) {
            let v: V3 = std::array::from_fn(|a| vel.iter().map(|(k, jv, _)| jv[a] * z[*k]).sum());
            let om: V3 = std::array::from_fn(|a| vel.iter().map(|(k, _, jw)| jw[a] * z[*k]).sum());
            let e = mass * dot(v, v) + dot(om, mv(iw, &om));
            match shares.iter_mut().find(|s| s.0 == *g) {
                Some(s) => s.1 += e,
                None => shares.push((g.clone(), e)),
            }
        }
        let total: f64 = shares.iter().map(|s| s.1).sum();
        shares.iter_mut().for_each(|s| s.1 /= total);
        shares.sort_by(|a, b| b.1.total_cmp(&a.1));
        let motion = |path: &str| -> Option<(bool, f64)> {
            let j = joints.iter().find(|j| j.path == path)?;
            let dv: f64 = value(j).iter().zip(&z).map(|(h, x)| h * x).sum();
            Some((j.turn, if j.turn { dv.to_degrees() } else { dv }))
        };
        let lead = shares.iter().find_map(|s| motion(&s.0).filter(|m| m.1 != 0.0)).map_or(1.0, |m| m.1);
        let moves = shares
            .iter()
            .filter(|s| s.1 >= 1e-3)
            .map(|(p, e)| {
                let m = motion(p);
                Moves { joint: p.clone(), turn: m.is_some_and(|m| m.0), motion: m.map(|m| m.1 / lead), energy: *e }
            })
            .collect();
        out.push(Mode { frequency: (1e6 * lam).sqrt() / (2.0 * std::f64::consts::PI), moves });
    }
    let pose = model.coords.iter().map(|c| (c.name.clone(), c.q.at(&ctx))).collect();
    Ok(Modes { t, pose, modes: out, springs, released, rigid, model })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    static REGISTRY: ink_lang::Registry = ink_lang::Registry::new(&[&crate::Dynamics]);

    fn of(src: &str) -> Result<Modes, InkError> {
        let d = Doc::load_str_in(&REGISTRY, src, "modes.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{} at {}", e.message, e.line));
        modes(&d.reloader().unwrap(), "root", 0.0)
    }

    fn mass(m: &Modes, label: &str) -> f64 {
        m.model.bodies.iter().find(|b| b.label == label).unwrap_or_else(|| panic!("no body {label}")).mass
    }

    /// The roots, in omega^2, of m1 m2 w^4 - (m1 k2 + m2 (k1 + k2)) w^2 +
    /// k1 k2 = 0: mass m1 on a spring k1 to ground and k2 to mass m2.
    fn two_mass(m1: f64, k1: f64, m2: f64, k2: f64) -> (f64, f64) {
        let (a, b, c) = (m1 * m2, -(m1 * k2 + m2 * (k1 + k2)), k1 * k2);
        let s = (b * b - 4.0 * a * c).sqrt();
        ((-b - s) / (2.0 * a), (-b + s) / (2.0 * a))
    }

    fn hz(w2: f64) -> f64 {
        w2.sqrt() / (2.0 * std::f64::consts::PI)
    }

    /// A printer axis: a stepper's 20-tooth GT2 pulley (pitch radius 20 /
    /// pi) at x = 0 and an idler of the same radius at x = 600, a carriage
    /// clamped into the belt's lower run at x = s = 300 (its clamps 5 mm
    /// either side: 295 mm to each wheel), the upper run 600 mm; its rotor the
    /// LDO-42STH40-1684AC's 53 g cm2 and its pulley's 178.7 g mm2 (the
    /// triage's), held by 50 x 0.45 N m (its holding stiffness), the belt
    /// 6 mm fibreglass GT2, EA 16903 N (Gates).
    fn axis(rotor: &str, idler: &str, hold: &str) -> String {
        format!(
            r#"(import "std")
            (def R (/ 20 pi))
            (def EA 16903N)
            (body cart (shape (box 40 30 20)) (attr density 1.27g/cm3))
            (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia {rotor}g*mm2))
            (body idler (shape (cylinder-z 6 10)) (attr density 2.7g/cm3) (attr spin-inertia (* 1g*mm2 {idler})))
            (assembly (axis (s 300))
              (instance cart cart (at) (slide x s))
              (instance rotor rotor (at) (turn z (rad->deg (/ s R))) (effort-limits -450N*mm 450N*mm){hold})
              (instance idler idler (at (move 600 0 0)) (turn z (rad->deg (/ s R))))
              (relate rotor idler (belt-span (vec (vec 0 R 0) (vec 600 R 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate idler cart (belt-span (vec (vec 600 (neg R) 0) (vec (+ s 5) (neg R) 0)) 50N 0.95 (tensile-stiffness EA)))
              (relate cart rotor (belt-span (vec (vec (- s 5) (neg R) 0) (vec 0 (neg R) 0)) 50N 0.95 (tensile-stiffness EA))))
            (assembly root (instance a (axis 300)))"#
        )
    }

    /// The two-mass belt axis (the triage's hand check): the idler without
    /// inertia, so its spans (295 mm under the carriage, 600 over) act in
    /// series, k2 = EA (1 / 295 + 1 / 895), and the rotor and pulley, J =
    /// 5478.7 g mm2, are m1 = J / R^2 on k1 = 22500 N mm/rad / R^2. Both
    /// modes match the closed form's to 1e-9. Without the belt the rotor
    /// alone on its holding stiffness rings at 322.5 Hz; with the rotor held
    /// rigidly, the carriage on its two spans at (1/2 pi) sqrt(k2 / m).
    #[test]
    fn a_belt_axis_rings_at_its_two_mass_frequencies() {
        let r = 20.0 / std::f64::consts::PI;
        let ea = 16903.0;
        let k2 = ea * (1.0 / 295.0 + 1.0 / 895.0);
        let hold = " (stiffness (* 1N*mm (* 50 450)))";
        let m = of(&axis("5478.7", "0", hold)).unwrap_or_else(|e| panic!("{}", e.message));
        let cart = mass(&m, "root/a/cart");
        let (lo, hi) = two_mass(5478.7 / (r * r), 22500.0 / (r * r), cart, k2);
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        assert_eq!(f.len(), 2, "{f:?}");
        let near = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-9;
        assert!(near(f[0], hz(1e6 * lo)) && near(f[1], hz(1e6 * hi)), "{f:?} against {} {}", hz(1e6 * lo), hz(1e6 * hi));
        assert_eq!(m.released, ["root/a/rotor", "root/a/idler", "root/a/cart"]);
        // The lowest mode is the carriage's, the higher the rotor's.
        assert!(m.modes[0].moves[0].joint == "root/a/cart" && m.modes[1].moves[0].joint == "root/a/rotor");
        // The rotor alone on its holding stiffness: 322.5 Hz.
        let alone = hz(1e6 * 22500.0 / 5478.7);
        assert!((alone - 322.5).abs() < 0.05, "{alone}");
        // A rotor held rigidly: the carriage on its spans.
        let fixed = of(&axis("5478.7", "0", " (stiffness 1e12N*mm)")).unwrap();
        let f0 = fixed.modes[0].frequency;
        assert!((f0 / hz(1e6 * k2 / cart) - 1.0).abs() < 1e-8, "{f0} against {}", hz(1e6 * k2 / cart));
    }

    /// The axis at its pose with a strap whose shape reads s: s is then a
    /// size, no coordinate (Dynamics, Coordinates), so the carriage, rotor
    /// and idler move only by their own released joints, and each is
    /// measured as their reactions carry it: the same two modes.
    #[test]
    fn a_pose_a_shape_reads_rings_alike() {
        let plain = of(&axis("5478.7", "0", " (stiffness 22500N*mm)")).unwrap();
        let src = axis("5478.7", "0", " (stiffness 22500N*mm)").replace(
            "(assembly (axis (s 300))",
            "(body (strap (s 0)) (shape (box (+ 1 (* 0.01 s)) 1 1)) (attr density 1g/cm3))\n(assembly (axis (s 300)) (instance strap (strap s) (at (move 0 50 0)))",
        );
        let sized = of(&src).unwrap();
        assert!(sized.pose.is_empty() && plain.pose.len() == 1, "{:?} {:?}", sized.pose, plain.pose);
        let f = |m: &Modes| m.modes.iter().map(|x| x.frequency).collect::<Vec<_>>();
        assert!(f(&sized).len() == 2 && f(&sized).iter().zip(f(&plain)).all(|(a, b)| (a / b - 1.0).abs() < 1e-9), "{:?} {:?}", f(&sized), f(&plain));
    }

    /// An idler whose inertia along its turn is within rounding of its heft
    /// (1e-12 g mm2, beside the 106 g mm2 its own diagonal holds) carries
    /// none, as a simulation reads it: condensed out, the axis rings at the
    /// same two modes as with none. Before the fix a third rang at the
    /// springs over that crumb, as one did in the pose above where the
    /// idler's centre measured a hair off its axis left m r^2 less its own
    /// as a crumb above 0 (4e27 Hz).
    #[test]
    fn an_idler_with_inertia_within_rounding_of_its_heft_has_none() {
        let none = of(&axis("5478.7", "0", " (stiffness 22500N*mm)")).unwrap();
        let crumb = of(&axis("5478.7", "1e-12", " (stiffness 22500N*mm)")).unwrap();
        let f = |m: &Modes| m.modes.iter().map(|x| x.frequency).collect::<Vec<_>>();
        assert!(f(&crumb).len() == 2 && f(&crumb).iter().zip(f(&none)).all(|(a, b)| (a / b - 1.0).abs() < 1e-9), "{:?} {:?}", f(&crumb), f(&none));
    }

    /// Two slides of 1000 g and 9000 g tied by a rigid line, the first held
    /// to ground by a taut 100 N/mm line: they move together (1 : 1) on the
    /// spring, (1/2 pi) sqrt(k / (m_a + m_b)), 15.92 Hz (before the fix a
    /// constraint read in the wrong scaling gave 47.75 Hz, 1 : 0.111).
    #[test]
    fn a_rigid_line_ties_unequal_masses_together() {
        let m = of(r#"(import "std")
            (body anchor (shape (box 1 1 1)))
            (body a (shape (box 100 100 100)) (attr density 1g/cm3))
            (body b (shape (box 100 100 100)) (attr density 9g/cm3))
            (fn elastic (fa fb sa sb) (relation elastic (couple (vec (vec (vec 0 0 0) (vec 1 0 0) 0 100 100)) 1N)))
            (fn rigid (fa fb sa sb) (relation rigid (couple (vec (vec (vec 0 0 0) (vec 1 0 0))) 0)))
            (assembly root
              (instance anchor anchor)
              (instance a a (at) (slide x 0))
              (instance b b (at (move 0 200 0)) (slide x 0))
              (relate anchor a elastic)
              (relate a b rigid))"#)
        .unwrap_or_else(|e| panic!("{}", e.message));
        let want = hz(1e6 * 100.0 / (mass(&m, "root/a") + mass(&m, "root/b")));
        assert!(m.modes.len() == 1 && (m.modes[0].frequency / want - 1.0).abs() < 1e-9, "{} against {want}", m.modes[0].frequency);
        assert!((want - 15.9155).abs() < 1e-3, "{want}");
        let motion = |j: &str| m.modes[0].moves.iter().find(|x| x.joint == j).unwrap().motion.unwrap();
        assert!((motion("root/b") / motion("root/a") - 1.0).abs() < 1e-9, "{} {}", motion("root/a"), motion("root/b"));
        assert_eq!(m.rigid, ["line 1 of root/a~root/b (rigid)"]);
    }

    /// A carriage held by a taut 100 N/mm line, tied rigidly at radius 10 to a
    /// pinion of 40000 g mm2 (a rack's tooth line along x at the pinion's
    /// top): one mode, the pinion turning -1 / r rad per mm of the carriage,
    /// on the carriage's mass plus J / r^2 = 400 g: (1/2 pi) sqrt(k / (m + J
    /// / r^2)), by hand.
    #[test]
    fn a_rigid_ratio_reflects_its_inertia() {
        let m = of(r#"(import "std")
            (body anchor (shape (box 1 1 1)))
            (body cart (shape (box 100 100 100)) (attr density 1g/cm3))
            (body pinion (shape (cylinder-z 10 10)) (attr density 7.85g/cm3) (attr spin-inertia 40000g*mm2))
            (fn elastic (fa fb sa sb) (relation elastic (couple (vec (vec (vec -300 0 0) (vec 1 0 0) 0 300 100)) 1N)))
            (fn rack (fa fb sa sb) (relation rack (couple (vec (vec (vec 0 10 0) (vec 1 0 0))) 0)))
            (assembly (axis (s 0))
              (instance anchor anchor (at (move -300 0 0)))
              (instance cart cart (at (move 0 200 0)) (slide x s))
              (instance pinion pinion (at) (turn z (rad->deg (/ (neg s) 10))))
              (relate anchor cart elastic)
              (relate cart pinion rack))
            (assembly root (instance a (axis 0)))"#)
        .unwrap_or_else(|e| panic!("{}", e.message));
        let want = hz(1e6 * 100.0 / (mass(&m, "root/a/cart") + 40000.0 / 100.0));
        assert!(m.modes.len() == 1 && (m.modes[0].frequency / want - 1.0).abs() < 1e-9, "{} against {want}", m.modes[0].frequency);
        let motion = |j: &str| m.modes[0].moves.iter().find(|x| x.joint == j).unwrap().motion.unwrap();
        assert!((motion("root/a/pinion") / motion("root/a/cart") + (0.1f64).to_degrees()).abs() < 1e-9, "{} {}", motion("root/a/cart"), motion("root/a/pinion"));
    }

    /// A joint its own stiffness holds that no coordinate drives (a constant
    /// drive) is released and its body measured: a 1 kg slide on 100 N/mm
    /// at 50.33 Hz; beside one a coordinate drives, both modes (before
    /// the fix: no modes, then one).
    #[test]
    fn a_joint_its_stiffness_holds_rings_without_a_coordinate() {
        let alone = of(r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly root (instance b b (at) (slide x 0) (stiffness 100N/mm)))"#)
        .unwrap();
        let k = |k: f64, m: &Modes, b: &str| hz(1e6 * k / mass(m, b));
        assert!(
            alone.modes.len() == 1 && (alone.modes[0].frequency / k(100.0, &alone, "root/b") - 1.0).abs() < 1e-12,
            "{:?}",
            alone.modes.iter().map(|x| x.frequency).collect::<Vec<_>>()
        );
        let mixed = of(r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (s 0))
              (instance a b (at) (slide x s) (stiffness 100N/mm))
              (instance b b (at (move 0 200 0)) (slide x 0) (stiffness 25N/mm)))
            (assembly root (instance a (axis 0)))"#)
        .unwrap();
        let f: Vec<f64> = mixed.modes.iter().map(|x| x.frequency).collect();
        let want = [k(25.0, &mixed, "root/a/b"), k(100.0, &mixed, "root/a/a")];
        assert!(f.len() == 2 && f.iter().zip(want).all(|(a, b)| (a / b - 1.0).abs() < 1e-12), "{f:?} against {want:?}");
        assert_eq!(mixed.released, ["root/a/b"]);
    }

    /// A belt's optional arguments keep their meaning: (belt-drive 10 10 50
    /// 0.95 0.8) back-drives at EB 0.8 (100 N mm on the output passes 80 to
    /// the motor), and its tensile stiffness is a tagged argument, its lines
    /// each EA / L; given twice, or on lines with no length, it is refused.
    #[test]
    fn a_belts_efficiencies_keep_their_meaning_beside_its_stiffness() {
        let src = |more: &str| {
            format!(
                r#"(import "std") (gravity 0 0 0)
                (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 1000g*mm2))
                (assembly (axis (q 0))
                  (instance a rotor (at) (turn z q) (effort-limits -100000N*mm 100000N*mm))
                  (instance b rotor (at (move 100 0 0)) (turn z q))
                  (relate a b (belt-drive 10 10 50N {more})))
                (assembly root (instance a (axis (* 10 (/ t 1s)))))
                (load backdrive root/a/b (torque 0 0 100N*mm))"#
            )
        };
        let effort = |more: &str| {
            let d = Doc::load_str_in(&REGISTRY, &src(more), "modes.ink", Path::new("."), &[]).unwrap_or_else(|e| panic!("{}", e.message));
            crate::build(&d.reloader().unwrap(), "root", [0.0; 3], 0.5, &[], false).unwrap().at(0.5).effort[0]
        };
        assert!((effort("0.95 0.8") + 80.0).abs() < 1e-6, "{}", effort("0.95 0.8"));
        assert!((effort("0.95 (tensile-stiffness 16903N)") + 95.0).abs() < 1e-6, "{}", effort("0.95 (tensile-stiffness 16903N)"));
        assert!((effort("0.95 0.8 (tensile-stiffness 16903N)") + 80.0).abs() < 1e-6);
        let m = of(&src("0.95 (tensile-stiffness 16903N)")).unwrap();
        assert!(
            m.springs.len() == 2 && m.springs.iter().all(|s| (s.stiffness - 169.03).abs() < 1e-9 && s.site == "modes.ink:6"),
            "{:?}",
            m.springs.iter().map(|s| (s.stiffness, &s.site)).collect::<Vec<_>>()
        );
        let err = |src: &str| Doc::load_str_in(&REGISTRY, src, "modes.ink", Path::new("."), &[]).err().unwrap_or_else(|| panic!("loaded: {src}")).message;
        assert!(err(&src("(tensile-stiffness 1N) (tensile-stiffness 2N)")).contains("given twice"));
        let screw = src("0.95").replace("(belt-drive 10 10 50N 0.95)", "(lead-screw 8 0.9 (tensile-stiffness 1000N))");
        assert!(err(&screw).contains("its length"), "{}", err(&screw));
    }

    /// A belt's tensile stiffness is a force, as a joint's stiffness says its
    /// unit: 10000N, 10kN and 2248.09lbf are each EA / L = 100 N/mm on its
    /// 100 mm spans; a length, a torque, a mass or a bare number is refused
    /// (it was read as N whatever it said), as is a negative or infinite one.
    #[test]
    fn a_belts_tensile_stiffness_is_a_force() {
        let src = |ea: &str| {
            format!(
                r#"(import "std") (gravity 0 0 0)
                (body rotor (shape (cylinder-z 6 10)) (attr density 7.85g/cm3) (attr spin-inertia 1000g*mm2))
                (assembly (axis (q 0))
                  (instance a rotor (at) (turn z q) (effort-limits -100000N*mm 100000N*mm))
                  (instance b rotor (at (move 100 0 0)) (turn z q))
                  (relate a b (belt-drive 10 10 50N 0.95 (tensile-stiffness {ea}))))
                (assembly root (instance a (axis (* 10 (/ t 1s)))))"#
            )
        };
        let lbf = 10000.0 / 4.4482216152605;
        for (ea, k) in [("10000N", 100.0), ("10kN", 100.0), (&format!("{lbf}lbf") as &str, 100.0), ("(* 1000 3.8lbf)", 38.0 * 4.4482216152605)] {
            let m = of(&src(ea)).unwrap_or_else(|e| panic!("{ea}: {}", e.message));
            assert!(m.springs.len() == 2 && m.springs.iter().all(|s| (s.stiffness - k).abs() < 1e-9 * k), "{ea}: {:?}", m.springs.iter().map(|s| s.stiffness).collect::<Vec<_>>());
        }
        for ea in ["10mm", "10N*m", "10kg", "10000"] {
            let e = Doc::load_str_in(&REGISTRY, &src(ea), "modes.ink", Path::new("."), &[]).err().unwrap_or_else(|| panic!("loaded: {ea}"));
            assert!(e.message.contains("a tensile stiffness (EA)") && e.message.contains("force"), "{ea}: {}", e.message);
        }
        // A negative or no number is refused: at load where it is known
        // there, else where the modes read it (a param's).
        let e = Doc::load_str_in(&REGISTRY, &src("-1000N"), "modes.ink", Path::new("."), &[]).err().expect("loaded: -1000N");
        assert!(e.message.contains("a tensile stiffness (EA) is -1000 N: a stiffness is at least 0 and finite"), "{}", e.message);
        // (Ink has no NaN: x / 0, inf - inf are 0, so (/ 1N 0) is no
        // stiffness, as 0N is; an overflow is infinite.)
        for ea in ["(* 1N (exp 1000))", "(neg (* 1N (exp 1000)))"] {
            let e = Doc::load_str_in(&REGISTRY, &src(ea), "modes.ink", Path::new("."), &[]).and_then(|_| of(&src(ea)).map(|_| ())).err().unwrap_or_else(|| panic!("modes of {ea}"));
            assert!(e.message.contains("a stiffness is at least 0 and finite"), "{ea}: {}", e.message);
        }
        let param = src("ea").replace("(gravity 0 0 0)", "(gravity 0 0 0) (param ea -1000N)");
        let e = of(&param).err().expect("a negative EA refused");
        assert!(e.message.contains("a stiffness is at least 0 and finite"), "{}", e.message);
    }

    /// A slide's stiffness is a force per length, a turn's a torque per
    /// radian: 100N*m on a slide is refused, as is 5N/mm on a turn; 2N*m on
    /// a turn is 2000 N mm per radian.
    #[test]
    fn a_joints_stiffness_has_its_motions_units() {
        let src = |m: &str, k: &str| {
            format!(
                r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (s 0)) (instance b b (at) ({m} s) (stiffness {k})))
            (assembly root (instance a (axis 0)))"#
            )
        };
        let err = |src: &str| Doc::load_str_in(&REGISTRY, src, "modes.ink", Path::new("."), &[]).err().map(|e| e.message);
        assert!(err(&src("slide x", "100N*m")).is_some_and(|e| e.contains("a slide's stiffness")), "{:?}", err(&src("slide x", "100N*m")));
        assert!(err(&src("turn z", "5N/mm")).is_some_and(|e| e.contains("a turn's stiffness")), "{:?}", err(&src("turn z", "5N/mm")));
        assert!(err(&src("slide x", "100N/mm")).is_none());
        let m = of(&src("turn z", "2N*m")).unwrap();
        assert!((m.springs[0].stiffness - 2000.0).abs() < 1e-9 && m.springs[0].turn, "{}", m.springs[0].stiffness);
    }

    /// A rotor of 1000 g mm2 on 20000 N mm/rad, a rigid belt-drive 10:20 to a
    /// 3000 g mm2 shaft, and a stiff belt-drive 10:10 (EA 10000, spans 200 mm)
    /// on to an 8000 g mm2 flywheel. By hand, in the
    /// rotor's and the flywheel's angles (the shaft turning half the rotor):
    /// M = diag(J_r + J_s / 4, J_f), K = [[k_h + k_b / 4, -k_b / 2], [-k_b / 2,
    /// k_b]], k_b = 2 (EA / 200) 10^2: 166.7819772 and 574.0399595 Hz (the
    /// rigid line read in the wrong scaling held nothing: 0, 340.73, 711.76).
    #[test]
    fn a_rigid_belt_on_a_released_shaft_reflects_its_ratio() {
        let m = of(r#"(import "std")
            (body rotor (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 1000g*mm2))
            (body shaft (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 3000g*mm2))
            (body fly (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 8000g*mm2))
            (assembly (m (s 0))
              (instance rotor rotor (at) (turn z s) (stiffness 20000N*mm))
              (instance shaft shaft (at (move 100 0 0)) (turn z (/ s 2)))
              (instance fly fly (at (move 100 200 0)) (turn z (/ s 2)))
              (relate rotor shaft (belt-drive 10 20 50N 0.95))
              (relate shaft fly (belt-drive 10 10 50N 0.95 (tensile-stiffness 10000N))))
            (assembly root (instance a (m 0)))"#)
        .unwrap_or_else(|e| panic!("{}", e.message));
        let (kh, kb): (f64, f64) = (20000.0, 2.0 * (10000.0 / 200.0) * 100.0);
        let (m1, m2) = (1000.0 + 3000.0 / 4.0, 8000.0);
        let (k11, k12, k22) = (kh + kb / 4.0, -kb / 2.0, kb);
        let (a, b, c) = (m1 * m2, -(k11 * m2 + k22 * m1), k11 * k22 - k12 * k12);
        let sq = (b * b - 4.0 * a * c).sqrt();
        let want = [hz(1e6 * (-b - sq) / (2.0 * a)), hz(1e6 * (-b + sq) / (2.0 * a))];
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        assert!(f.len() == 2 && f.iter().zip(want).all(|(x, y)| (x / y - 1.0).abs() < 1e-9), "{f:?} against {want:?}");
        assert!((want[0] - 166.7819772).abs() < 1e-6 && (want[1] - 574.0399595).abs() < 1e-6, "{want:?}");
        assert_eq!(m.rigid, ["line 1 of root/a/rotor~root/a/shaft (belt-drive)", "line 2 of root/a/rotor~root/a/shaft (belt-drive)"]);
    }

    /// A carriage bolted at x = 0 (a constant slide) carrying a pulley of
    /// 500 g mm2 on its own turn, held to an anchor by two spans of EA / 200
    /// = 50 N/mm at +-10: the slide supports the spans'
    /// end, so it is released too, and both modes are there: the carriage
    /// and pulley along x on 2k, (1/2 pi) sqrt(2k / (m_c + m_p)), and the
    /// pulley's turn on 2k r^2, 711.76 Hz (only the turn's, before).
    #[test]
    fn a_joint_supporting_a_springs_end_is_released() {
        let m = of(r#"(import "std")
            (body anchor (shape (box 2 2 2)) (attr density 1g/cm3))
            (body carriage (shape (box 40 40 10)) (attr density 1g/cm3))
            (body pulley (shape (cylinder-z 10 6)) (attr density 1g/cm3) (attr spin-inertia 500g*mm2))
            (assembly (m (c 0))
              (instance a anchor (at (move -300 0 0)))
              (instance cart carriage (at (move 0 0 -20)) (slide x 0))
              (instance p pulley (at (move 0 0 20) (frame-of cart)) (turn z 0))
              (relate a p (belt-span (vec (vec -200 10 0) (vec c 10 0)) 50N 0.95 (tensile-stiffness 10000N)))
              (relate p a (belt-span (vec (vec c -10 0) (vec -200 -10 0)) 50N 0.95 (tensile-stiffness 10000N))))
            (assembly root (instance m (m 0)))"#)
        .unwrap_or_else(|e| panic!("{}", e.message));
        let along = hz(1e6 * 100.0 / (mass(&m, "root/m/cart") + mass(&m, "root/m/p")));
        let turn = hz(1e6 * 100.0 * 100.0 / 500.0);
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        assert!(f.len() == 2 && (f[0] / along - 1.0).abs() < 1e-9 && (f[1] / turn - 1.0).abs() < 1e-9, "{f:?} against {along} {turn}");
        assert!((turn - 711.76).abs() < 0.01 && (along - 376.3).abs() < 1.0, "{along} {turn}");
        assert_eq!(m.released, ["root/m/p", "root/m/cart"]);
    }

    /// An arm turned a fifth of a belt-driven shaft by their drives alone,
    /// no relation tying it: releasing the shaft for its
    /// belt cuts that tie, so it is refused, naming the arm (it was a 0 Hz
    /// mode, its inertia gone from the others).
    #[test]
    fn a_body_tied_only_by_drives_to_a_released_joint_is_refused() {
        let e = of(r#"(import "std")
            (body rotor (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 1000g*mm2))
            (body shaft (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 3000g*mm2))
            (body arm (shape (cylinder-z 5 10)) (attr density 1g/cm3) (attr spin-inertia 50000g*mm2))
            (assembly (m (s 0))
              (instance rotor rotor (at) (turn z (* 2 s)) (stiffness 20000N*mm))
              (instance shaft shaft (at (move 200 0 0)) (turn z s))
              (instance arm arm (at (move 200 0 50)) (turn z (/ s 5)))
              (relate rotor shaft (belt-drive 10 20 50N 0.95 (tensile-stiffness 10000N))))
            (assembly root (instance a (m 0)))"#)
        .err()
        .expect("refused");
        assert!(e.message.contains("root/a/arm (with a.s)") && e.hint.as_ref().is_some_and(|h| h.contains("gear-mesh")), "{} | {:?}", e.message, e.hint);
    }

    /// A body a coordinate places with no joint, beside a slide held by
    /// 100 N/mm on that coordinate: one mode, both bodies' mass on the
    /// spring, 50.33 / sqrt 2 Hz; the body gives its energy share alone.
    #[test]
    fn a_body_no_joint_moves_shares_its_coordinates_mode() {
        let m = of(r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (axis (s 0))
              (instance b b (at) (slide x s) (stiffness 100N/mm))
              (instance c b (at (move s 300 0))))
            (assembly root (instance a (axis 0)))"#)
        .unwrap();
        let want = hz(1e6 * 100.0 / (mass(&m, "root/a/b") + mass(&m, "root/a/c")));
        assert!(m.modes.len() == 1 && (m.modes[0].frequency / want - 1.0).abs() < 1e-12, "{}", m.modes[0].frequency);
        let c = m.modes[0].moves.iter().find(|x| x.joint == "root/a/c").expect("the body's share");
        assert!(c.motion.is_none() && (c.energy - 0.5).abs() < 1e-6, "{:?} {}", c.motion, c.energy);
    }

    /// A joint's stiffness clause: once per instance, on a motion, of params
    /// and not of x y z or t; below 0 its modes are refused.
    #[test]
    fn a_joints_stiffness_clause_is_checked() {
        let src = |clauses: &str| {
            format!(
                r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
                (assembly (axis (s 0)) (instance b b (at) {clauses}))
                (assembly root (instance a (axis 0)))"#
            )
        };
        let err = |src: &str| Doc::load_str_in(&REGISTRY, src, "modes.ink", Path::new("."), &[]).err().map(|e| e.message).unwrap_or_default();
        assert!(err(&src("(slide x s) (stiffness 1N/mm) (stiffness 2N/mm)")).contains("stiffness"));
        assert!(err(&src("(stiffness 100N/mm)")).contains("stiffness"));
        assert!(err(&src("(slide x s) (stiffness (* 1N/mm (* 100 (/ t 1s))))")).contains("not of x y z or t"));
        let e = of(&src("(slide x s) (stiffness -5N/mm)")).err().expect("refused");
        assert!(e.message.contains("at least 0"), "{}", e.message);
    }

    /// A gripper's two 1 kg jaws, each on its own actuator, held by 100 and
    /// 25 N/mm and written on one width w: two springs,
    /// two modes, 50.33 and 25.16 Hz (one was welded into one at 39.79 Hz).
    #[test]
    fn two_stiff_joints_on_one_coordinate_ring_apart() {
        let m = of(r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
            (assembly (grip (w 0))
              (instance l b (at (move -100 0 0)) (slide x (neg w)) (stiffness 100N/mm))
              (instance r b (at (move 100 0 0)) (slide x w) (stiffness 25N/mm)))
            (assembly root (instance g (grip 0)))"#)
        .unwrap();
        let k = |k: f64, b: &str| hz(1e6 * k / mass(&m, b));
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        let want = [k(25.0, "root/g/r"), k(100.0, "root/g/l")];
        assert!(f.len() == 2 && f.iter().zip(want).all(|(a, b)| (a / b - 1.0).abs() < 1e-12), "{f:?} against {want:?}");
    }

    /// A belt's tensile stiffness written as a number after its efficiency
    /// (the slot is EB) is refused as effort refuses it, not read as a rigid
    /// belt with efficiency 16903 back (291 Hz before the fix).
    /// So is a pretension below 0, and an efficiency forward above 1, with
    /// no hint about a stiffness (effort's own check, Coupled::unsound).
    #[test]
    fn a_belts_numbers_are_checked_as_effort_checks_them() {
        let src = || axis("5478.7", "0", " (stiffness 22500N*mm)");
        let e = of(&src().replace("(tensile-stiffness EA)", "EA")).err().expect("refused");
        assert!(e.message.contains("16903 back") && e.hint.as_ref().is_some_and(|h| h.contains("(tensile-stiffness EA)")), "{} | {:?}", e.message, e.hint);
        let e = of(&src().replacen("50N 0.95 (tensile-stiffness EA)", "-5N 0.95 (tensile-stiffness EA)", 1)).err().expect("refused");
        assert!(e.message.contains("is -5: it is a tension"), "{}", e.message);
        let e = of(&src().replacen("50N 0.95 (tensile-stiffness EA)", "50N 1.5 (tensile-stiffness EA)", 1)).err().expect("refused");
        assert!(e.message.contains("1.5 forward") && e.hint.is_none(), "{} | {:?}", e.message, e.hint);
    }

    /// The same axis with an idler of 2000 g mm2 and the rotor on 5000 N
    /// mm/rad: three masses (the carriage x, the rotor's and the idler's
    /// rims a and b), k_L (x - a)^2 + k_R (b - x)^2 + k_T (a - b)^2 + (k_h /
    /// R^2) a^2, whose generalized eigenproblem, written by hand, has the
    /// three roots of its characteristic cubic: ink's within 1e-9.
    #[test]
    fn a_belt_axis_with_an_idler_rings_at_its_three_mass_frequencies() {
        let r = 20.0 / std::f64::consts::PI;
        let ea = 16903.0;
        let m = of(&axis("5478.7", "2000", " (stiffness 5000N*mm)")).unwrap_or_else(|e| panic!("{}", e.message));
        let cart = mass(&m, "root/a/cart");
        let (kl, kr, kt, kh) = (ea / 295.0, ea / 295.0, ea / 600.0, 5000.0 / (r * r));
        let mm = [cart, 5478.7 / (r * r), 2000.0 / (r * r)];
        let kk = [[kl + kr, -kl, -kr], [-kl, kl + kt + kh, -kt], [-kr, -kt, kr + kt]];
        let det = |w2: f64| {
            let a: Vec<Vec<f64>> = (0..3).map(|i| (0..3).map(|j| kk[i][j] - if i == j { w2 * mm[i] } else { 0.0 }).collect()).collect();
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
        };
        // The cubic's roots, by bisection between its sign changes on a
        // fine logarithmic grid.
        let mut roots = Vec::new();
        let grid: Vec<f64> = (0..4000).map(|i| 1e-4 * 10f64.powf(i as f64 / 400.0)).collect();
        for w in grid.windows(2) {
            let (mut a, mut b) = (w[0], w[1]);
            if det(a).signum() == det(b).signum() {
                continue;
            }
            for _ in 0..200 {
                let c = 0.5 * (a + b);
                if det(a).signum() == det(c).signum() { a = c } else { b = c }
            }
            roots.push(hz(1e6 * 0.5 * (a + b)));
        }
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        assert_eq!(roots.len(), 3, "{roots:?}");
        assert!(f.len() == 3 && f.iter().zip(&roots).all(|(a, b)| (a / b - 1.0).abs() < 1e-9), "{f:?} against {roots:?}");
    }

    /// A slide held by a spring of 100 N/mm carrying 1 kg: 50.33 Hz; the
    /// same without the spring has nothing stating a stiffness, refused.
    #[test]
    fn a_joint_held_by_its_stiffness_rings_and_nothing_held_is_refused() {
        let src = |k: &str| {
            format!(
                r#"(import "std") (body b (shape (box 100 100 100)) (attr density 1g/cm3))
                (assembly (axis (s 0)) (instance b b (at) (slide x s){k}))
                (assembly root (instance a (axis 0)))"#
            )
        };
        let m = of(&src(" (stiffness 100N/mm)")).unwrap();
        let want = hz(1e6 * 100.0 / mass(&m, "root/a/b"));
        assert!(m.modes.len() == 1 && (m.modes[0].frequency / want - 1.0).abs() < 1e-12, "{}", m.modes[0].frequency);
        assert!((want - 50.33).abs() < 0.01, "{want}");
        let e = of(&src("")).err().expect("refused");
        assert!(e.message.contains("states a stiffness"), "{}", e.message);
    }

    /// Two 1 kg sliders on their own springs, 1 and 1e15 N/mm: the soft
    /// one's mode is its own, sqrt(1e6 x 1 / 1000) / 2 pi = 5.0329 Hz,
    /// however stiff the other (once read as 0, under 1e-13 of the stiff
    /// one's).
    #[test]
    fn a_soft_modes_frequency_is_independent_of_a_stiff_neighbor() {
        let m = of(r#"(import "std") (body b (shape (box 50 50 50)) (attr mass 1000g))
(assembly root (instance soft b (at) (slide x 0) (stiffness 1N/mm))
 (instance stiff b (at (move 0 100 0)) (slide x 0) (stiffness 1e15N/mm)))"#)
        .unwrap();
        assert!((m.modes[0].frequency - hz(1000.0)).abs() < 1e-7, "{}", m.modes[0].frequency);
        // Beside them two 1 kg carts on one taut span, 1e15 N/mm (EA 1e17 over 100 mm): together they are free,
        // 0, the span's own stretch being its rounding; apart, k (1/m + 1/m), 1e6 x 1e15 x 2 / 1000.
        let m = of(r#"(import "std") (body b (shape (box 50 50 50)) (attr mass 1000g))
(assembly root (instance soft b (at (move 0 300 0)) (slide x 0) (stiffness 1N/mm))
 (instance a b (at) (slide x 0)) (instance c b (at (move 100 0 0)) (slide x 0))
 (relate a c (belt-span (vec (vec 0 0 0) (vec 100 0 0)) 1N 1 (tensile-stiffness 1e17N))))"#)
        .unwrap();
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        assert!(f.len() == 3 && f[0] == 0.0 && (f[1] - hz(1000.0)).abs() < 1e-7 && (f[2] / hz(2e18) - 1.0).abs() < 1e-9, "{f:?}");
    }

    /// Three 1 kg carts in a row, nothing to ground: a to c on a taut span of 1e15 N/mm, c to e on one of 1 N/mm
    /// (EA 100 over 100 mm). Together they are free, 0 Hz, however far apart the two stiffnesses lie: the
    /// soft span's mode is a and c as one 2 kg mass against e, sqrt(1e6 x 1 x (1/2000 + 1/1000)) / 2 pi.
    #[test]
    fn a_free_chain_of_a_stiff_and_a_soft_span_is_free() {
        let m = of(r#"(import "std") (body b (shape (box 50 50 50)) (attr mass 1000g))
(assembly root (instance a b (at) (slide x 0)) (instance c b (at (move 100 0 0)) (slide x 0)) (instance e b (at (move 200 0 0)) (slide x 0))
 (relate a c (belt-span (vec (vec 0 0 0) (vec 100 0 0)) 1N 1 (tensile-stiffness 1e17N)))
 (relate c e (belt-span (vec (vec 100 0 0) (vec 200 0 0)) 1N 1 (tensile-stiffness 100N))))"#)
        .unwrap();
        let f: Vec<f64> = m.modes.iter().map(|x| x.frequency).collect();
        // The soft mode to rounding: read among the modes below the stiff one, where the stiff span barely stretches.
        assert!(f.len() == 3 && f[0] == 0.0 && (f[1] / hz(1500.0) - 1.0).abs() < 1e-9, "{f:?}");
    }

    #[test]
    fn opposed_unstretched_belts_have_one_spans_stiffness() {
        let src = r#"(import "std")
        (body block (shape (box 50 50 50)) (attr mass 1000g))
        (assembly root
          (instance left block (at (move -100 0 0)))
          (instance cart block (at) (slide x 0))
          (instance right block (at (move 100 0 0)))
          (relate left cart (belt-span (vec (vec -100 0 0) (vec 0 0 0)) 0 1 (tensile-stiffness 1000N)))
          (relate cart right (belt-span (vec (vec 0 0 0) (vec 100 0 0)) 0 1 (tensile-stiffness 1000N))))"#;
        let m = of(src).unwrap();
        assert!((m.modes[0].frequency - hz(10000.0)).abs() < 1e-6, "{}", m.modes[0].frequency);
        let rounded = src.replace("(vec 100 0 0)", "(vec 100 0.000005 0)");
        let m = of(&rounded).unwrap();
        assert!((m.modes[0].frequency - hz(10000.0)).abs() < 1e-6);
        let mismatched = src.replacen("tensile-stiffness 1000", "tensile-stiffness 500", 1);
        assert!(of(&mismatched).err().unwrap().message.contains("both directions of displacement"));
        let taut = of(&src.replace("0 1 (tensile-stiffness", "1N 1 (tensile-stiffness")).unwrap();
        assert!((taut.modes[0].frequency - hz(20000.0)).abs() < 1e-6);
    }
}
