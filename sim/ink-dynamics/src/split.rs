//! A reaction divided among the parts that carry it (`(reaction PATH (by S
//! ...))`), by statics: each support carries, at the frame its part
//! declares its load at (`load`), the components its record's
//! `load-ratings` say it rates (Carries), less what PATH's joint leaves free
//! (no part carries a force along a slide's motion or a moment about a
//! turn's axis: a roller rolls, a carriage slides); the joint's drive passes
//! that, along a slide's line where an axial part (a screw's nut) lies on
//! it; their shares sum to the reaction. Guide carriages of one record (two
//! on a rail, four under a table) share by the rigid-table rule HIWIN's and
//! THK's catalogs size them by (HIWIN G99TE24-2410, Table 1-3): a rigid
//! table on identical blocks of equal stiffness, each block's moments
//! carrying only what the blocks' forces cannot (a single rail's roll), then
//! the forces of least sum of squares: the equations' solution as the
//! blocks' moment stiffness is small beside their forces over their
//! spacing, the assumption the answer states. Where they do not fix the
//! shares otherwise (more unknowns than the arrangement's equations: three
//! bearings on a shaft, carriages of different records) or leave some of
//! the reaction to none of them (an axial force on two radial bearings), the
//! split is refused, saying which: a share is never guessed.

use crate::model::{FLOAT, Wrench};
use ink_core::v3::{cross, dot, norm};

type V3 = [f64; 3];

/// What a support carries, in its load frame, at its origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Carries {
    /// A radial bearing or bushing (direction radial): force across its z.
    Radial,
    /// A thrust bearing, a screw's nut (direction axial): force along its z.
    Axial,
    /// A guide carriage (moment ratings MR, MP, MY): force across its x (its
    /// travel) and moment about each of its axes.
    Guide,
}

impl Carries {
    /// What a record's `load-ratings` say its part carries; None where they
    /// say neither its direction nor a guide's moments.
    pub fn of(ratings: &serde_json::Value) -> Option<Carries> {
        let moments = ratings.get("moment-rating-Nm").and_then(|m| m.as_object());
        if moments.is_some_and(|m| ["MR", "MP", "MY"].iter().any(|k| m.contains_key(*k))) {
            return Some(Carries::Guide);
        }
        match ratings.get("direction").and_then(|d| d.as_str()) {
            Some("radial") => Some(Carries::Radial),
            Some("axial") => Some(Carries::Axial),
            _ => None,
        }
    }

    /// The axes of its load frame it carries a force along, and a moment
    /// about.
    fn axes(self) -> (&'static [usize], &'static [usize]) {
        match self {
            Carries::Radial => (&[0, 1], &[]),
            Carries::Axial => (&[2], &[]),
            Carries::Guide => (&[1, 2], &[0, 1, 2]),
        }
    }

    /// Whether an axial part carries a pull as well as a push along its
    /// load frame's z: where its record says its thrust is double-direction
    /// (a screw's nut); a single-direction one (a thrust ball bearing,
    /// pushed from its housing washer toward its shaft washer, +z), or one
    /// whose record does not say, pushes only.
    pub fn pulls(ratings: &serde_json::Value) -> bool {
        ratings.get("thrust").and_then(|t| t.as_str()) == Some("double-direction")
    }

    pub fn says(self) -> &'static str {
        match self {
            Carries::Radial => "force across its axis (z)",
            Carries::Axial => "force along its axis (z)",
            Carries::Guide => "force across its travel (x) and moment about each axis",
        }
    }
}

/// The path's joint: a turn passes a torque about its frame's axis, a slide
/// a force along it, through its origin.
#[derive(Clone, Copy)]
pub struct Drive {
    pub turn: bool,
    pub axis: usize,
}

/// A part a reaction is split among, at a moment: its path, what it
/// carries, its load frame and the frame its share is given in (each R row
/// major, then T, in the world).
pub struct Part<'a> {
    pub path: &'a str,
    pub carries: Carries,
    /// An axial part that carries a pull too (Carries::pulls).
    pub pulls: bool,
    pub load: [f64; 12],
    pub given: [f64; 12],
    /// The force and moment (about the world's origin) a coupling it is a
    /// part of exerts on it: an axial part (a screw's nut, its thread the
    /// coupling's) carries its force along its axis besides its share.
    pub coupled: Option<(V3, V3)>,
    /// Its record's load-ratings, as written: carriages alike in them are
    /// identical blocks (the rigid-table rule).
    pub rating: &'a str,
}

/// Whether parts (what each carries, its record's load-ratings) are guide
/// carriages of one record, two or more, beside only axial parts (a screw's
/// nut): a table whose carriages' shares the rigid-table rule gives (split
/// refuses it where their rails are not parallel, or where statics leaves
/// another part's share open).
pub fn table(parts: &[(Carries, &str)]) -> bool {
    let guides: Vec<&str> = parts.iter().filter(|p| p.0 == Carries::Guide).map(|p| p.1).collect();
    guides.len() >= 2 && guides.iter().all(|r| *r == guides[0]) && parts.iter().all(|p| matches!(p.0, Carries::Guide | Carries::Axial))
}

/// What a split by the rigid-table rule says of itself.
pub const TABLE: &str = "rest on a rigid table on identical blocks of equal stiffness (HIWIN G99TE24-2410, Table 1-3; THK's load calculation): each \
    block's moments carry only what the blocks' forces cannot (a single rail's roll, alike), and the forces are those of least sum of squares, \
    the catalogs' formulas; that is the catalogs' assumption, not the model's stiffness";

/// The least-norm solution map of the columns `cols` (A^+, one row per
/// column): (A^T A)^+ A^T, directions whose squared size is under 1e-12 of
/// `scale`^2 (or of the largest) taken as none. `scale` is the size of the
/// columns they were made from: the forces' for the forces, the moments'
/// for the moments, so a moment column is never judged against a force's.
fn least_norm(cols: &[[f64; 6]], scale: f64) -> Vec<[f64; 6]> {
    let n = cols.len();
    if n == 0 {
        return Vec::new();
    }
    let ata: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| (0..6).map(|k| cols[i][k] * cols[j][k]).sum()).collect()).collect();
    let (e, v) = crate::model::eigen(&ata);
    let largest = e.iter().map(|x| x.abs()).fold(scale * scale, f64::max);
    let keep: Vec<f64> = e.iter().map(|&x| if x > 1e-12 * largest { 1.0 / x } else { 0.0 }).collect();
    let inv: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j| (0..n).map(|k| v[i][k] * keep[k] * v[j][k]).sum()).collect()).collect();
    (0..n).map(|i| std::array::from_fn(|k| (0..n).map(|j| inv[i][j] * cols[j][k]).sum())).collect()
}

/// The rigid-table rule's map from the reaction to the unknowns (one row
/// per column of `a`, `moment[j]` whether column j is a carriage's moment):
/// the carriages' moments of least norm that leave the rest to what the
/// other columns can carry, then those of least norm.
fn rigid_table(a: &[[f64; 6]], moment: &[bool]) -> Vec<[f64; 6]> {
    let (ys, ms): (Vec<usize>, Vec<usize>) = (0..a.len()).partition(|&j| !moment[j]);
    let scale = |js: &[usize]| js.iter().map(|&j| a[j].iter().map(|x| x * x).sum::<f64>().sqrt()).fold(0.0, f64::max);
    let ay: Vec<[f64; 6]> = ys.iter().map(|&j| a[j]).collect();
    let xy = least_norm(&ay, scale(&ys));
    // What the other columns cannot carry: the projection off their range.
    let off: [[f64; 6]; 6] = std::array::from_fn(|r| std::array::from_fn(|c| f64::from(u8::from(r == c)) - (0..ay.len()).map(|i| ay[i][r] * xy[i][c]).sum::<f64>()));
    let bm: Vec<[f64; 6]> = ms.iter().map(|&j| std::array::from_fn(|r| (0..6).map(|c| off[r][c] * a[j][c]).sum())).collect();
    let xm = least_norm(&bm, scale(&ms));
    let gm: Vec<[f64; 6]> = xm.iter().map(|row| std::array::from_fn(|c| (0..6).map(|r| row[r] * off[r][c]).sum())).collect();
    // The rest: what the moments leave, b - A_m G_m b.
    let rest: [[f64; 6]; 6] = std::array::from_fn(|r| std::array::from_fn(|c| f64::from(u8::from(r == c)) - ms.iter().zip(&gm).map(|(&j, g)| a[j][r] * g[c]).sum::<f64>()));
    let gy: Vec<[f64; 6]> = xy.iter().map(|row| std::array::from_fn(|c| (0..6).map(|r| row[r] * rest[r][c]).sum())).collect();
    let mut g = vec![[0.0; 6]; a.len()];
    for (&j, row) in ys.iter().zip(gy) {
        g[j] = row;
    }
    for (&j, row) in ms.iter().zip(gm) {
        g[j] = row;
    }
    g
}

/// A frame's 12 numbers (R row major, then T): its axis k and origin.
fn axis(f: &[f64; 12], k: usize) -> V3 {
    [f[k], f[3 + k], f[6 + k]]
}

fn origin(f: &[f64; 12]) -> V3 {
    [f[9], f[10], f[11]]
}

/// What of the directions `vs` is square to unit d: an orthonormal basis of
/// their span with d taken out.
fn without(vs: &[V3], d: V3) -> Vec<V3> {
    let mut out: Vec<V3> = Vec::new();
    for v in vs {
        let mut u: V3 = std::array::from_fn(|k| v[k] - dot(*v, d) * d[k]);
        for b in &out {
            let s = dot(u, *b);
            u = std::array::from_fn(|k| u[k] - s * b[k]);
        }
        let n = norm(u);
        if n > 1e-9 {
            out.push(u.map(|x| x / n));
        }
    }
    out
}

/// The shares of reaction `w` its parts carry, and the path's drive's along
/// its joint: per part its (fx fy fz mx my mz) in its `given` frame, N and N
/// mm about that frame's origin, and the linear map from the reaction
/// (force, then moment about `w`'s frame's origin, world axes) to those
/// components (for a band on the reaction); Err where statics leaves them
/// undetermined, or a part of the reaction to none of them.
#[cfg(test)]
pub fn split(w: &Wrench, parts: &[Part], drive: Option<Drive>) -> Result<(Vec<[f64; 6]>, Vec<[[f64; 6]; 6]>), String> {
    split_located(w, parts, drive, false).map(|(shares, map, _)| (shares, map))
}

/// `split`, and with `located` the thrust radial parts that locate a shaft carry among them: where the parts
/// are radial bearings about a turning joint's axis, none an axial part, and some of them carry thrust (their
/// record's `thrust`, double-direction: a deep-groove ball bearing), the force along the axis is theirs, its
/// size the third value, and their radial shares are what they are without it. Which of them carries it,
/// statics does not say (it is the one the design locates the shaft with): it is in no share.
pub fn split_located(w: &Wrench, parts: &[Part], drive: Option<Drive>, located: bool) -> Result<(Vec<[f64; 6]>, Vec<[[f64; 6]; 6]>, Option<f64>), String> {
    let o = origin(&w.frame);
    let free = drive.map(|d| (d.turn, axis(&w.frame, d.axis)));
    // On a slide, each part is set to its motion: a carriage's rail (its x)
    // along it; a bushing's or a nut's bore (its z) along it, or a
    // roller's or a thrust bearing's across it. One set at a slant would
    // move against what it carries.
    if let Some((false, d)) = free {
        for p in parts {
            let (k, across) = if p.carries == Carries::Guide { (0, false) } else { (2, true) };
            let c = dot(axis(&p.load, k), d).abs();
            if !(c > 1.0 - 1e-6 || (across && c < 1e-6)) {
                return Err(format!(
                    "{} is not set to the slide's motion: its {} is at {:.3} deg to it, where it must be along it{}; a part at a slant would move \
                     against what it carries",
                    p.path,
                    if k == 0 { "x (its rail)" } else { "z (its bore)" },
                    c.min(1.0).acos().to_degrees(),
                    if across { " or across it" } else { "" }
                ));
            }
        }
    }
    // The axial parts along a slide (screws' nuts) are what pass its drive,
    // each along its own line: two of them share it as their screws'
    // stiffness says, which statics does not.
    let nuts: Vec<usize> = match free {
        Some((false, d)) => (0..parts.len()).filter(|&i| parts[i].carries == Carries::Axial && dot(axis(&parts[i].load, 2), d).abs() > 1.0 - 1e-6).collect(),
        _ => Vec::new(),
    };
    // Unknowns: each part's force along each direction it carries at its
    // load frame's origin and moment about each axis, less the joint's
    // freedom; then the drive's.
    let mut cols: Vec<[f64; 6]> = Vec::new();
    let mut whose: Vec<(usize, V3, bool)> = Vec::new();
    for (s, part) in parts.iter().enumerate() {
        let p: V3 = std::array::from_fn(|k| origin(&part.load)[k] - o[k]);
        let (fs, ms) = part.carries.axes();
        let mut forces: Vec<V3> = fs.iter().map(|&k| axis(&part.load, k)).collect();
        let mut moments: Vec<V3> = ms.iter().map(|&k| axis(&part.load, k)).collect();
        match free {
            Some((false, d)) if !nuts.contains(&s) => forces = without(&forces, d),
            Some((true, d)) => moments = without(&moments, d),
            _ => {}
        }
        for u in forces {
            cols.push(force(u, p));
            whose.push((s, u, true));
        }
        for u in moments {
            cols.push([0.0, 0.0, 0.0, u[0], u[1], u[2]]);
            whose.push((s, u, false));
        }
    }
    // The thrust of radial bearings that locate the shaft: one unknown along their axis, through the first's centre.
    let locating = located
        && drive.is_some_and(|d| d.turn)
        && !parts.iter().any(|p| p.carries == Carries::Axial)
        && parts.iter().all(|p| p.carries == Carries::Radial && dot(axis(&p.load, 2), axis(&parts[0].load, 2)).abs() > 1.0 - 1e-6)
        && parts.iter().any(|p| p.pulls);
    if let Some(first) = parts.first().filter(|_| locating) {
        let u = axis(&first.load, 2);
        cols.push(force(u, std::array::from_fn(|k| origin(&first.load)[k] - o[k])));
        whose.push((parts.len() + 1, u, true));
    }
    let drove = drive.filter(|_| nuts.is_empty());
    if let Some(d) = drove {
        let u = axis(&w.frame, d.axis);
        cols.push(if d.turn { [0.0, 0.0, 0.0, u[0], u[1], u[2]] } else { force(u, [0.0; 3]) });
        whose.push((parts.len(), u, !d.turn));
    }
    let names = || parts.iter().map(|p| p.path).collect::<Vec<_>>().join(", ");
    let kinds = || {
        let mut k: Vec<String> = parts.iter().map(|p| format!("{}: {}", p.path, p.carries.says())).collect();
        k.extend(match drive {
            Some(_) if !nuts.is_empty() => {
                Some(format!("its joint's drive: along its motion, on the line of {}", nuts.iter().map(|&n| parts[n].path).collect::<Vec<_>>().join(", ")))
            }
            Some(d) => Some(format!("its joint's drive: along its motion, which none of them carries ({})", if d.turn { "a torque about its axis" } else { "a force along it" })),
            None => None,
        });
        k.join("; ")
    };
    let undetermined = || {
        format!(
            "statics leaves the shares of {} undetermined ({} unknowns, {} of them independent, against the 6 equations of what the path carries: {}): \
             they follow the parts' stiffness, which the model does not give, so none is guessed (guide carriages of one record share as a rigid \
             table on identical blocks)",
            names(),
            cols.len(),
            rank(&cols),
            kinds()
        )
    };
    let table = table(&parts.iter().map(|p| (p.carries, p.rating)).collect::<Vec<_>>());
    if cols.len() > 6 && !table {
        return Err(undetermined());
    }
    // Moments in N mm over a length scale, so that force and moment rows
    // weigh alike in the solve.
    let l = parts.iter().map(|p| norm(std::array::from_fn(|k| origin(&p.load)[k] - o[k]))).fold(1.0, f64::max);
    let scaled = |c: &[f64; 6]| -> [f64; 6] { std::array::from_fn(|k| if k < 3 { c[k] } else { c[k] / l }) };
    let a: Vec<[f64; 6]> = cols.iter().map(scaled).collect();
    let n = a.len();
    // G: the unknowns from the reaction (scaled). Carriages of one record
    // under a table: the rigid-table rule; else statics alone, by QR
    // (modified Gram-Schmidt, twice over each column), G = R^-1 Q^T.
    let g = if table {
        // The rule fixes how carriages on parallel rails share what they
        // carry, and nothing else: every other part's share (a screw's nut,
        // the drive) must be fixed by statics, whatever the carriages carry.
        let guide = |&(s, _, _): &(usize, V3, bool)| s < parts.len() && parts[s].carries == Carries::Guide;
        let rails = parts.iter().filter(|p| p.carries == Carries::Guide).map(|p| axis(&p.load, 0)).collect::<Vec<_>>();
        if rails.iter().any(|x| dot(*x, rails[0]).abs() < 1.0 - 1e-6) {
            return Err(undetermined());
        }
        let ours: Vec<[f64; 6]> = a.iter().zip(&whose).filter(|(_, w)| guide(w)).map(|(c, _)| *c).collect();
        if rank(&a) - rank(&ours) < a.len() - ours.len() {
            return Err(undetermined());
        }
        let moment: Vec<bool> = whose.iter().map(|w| guide(w) && !w.2).collect();
        rigid_table(&a, &moment)
    } else {
        let mut q: Vec<[f64; 6]> = Vec::with_capacity(n);
        let mut r = vec![vec![0.0; n]; n];
        for j in 0..n {
            let mut v = a[j];
            for _ in 0..2 {
                for (i, qi) in q.iter().enumerate() {
                    let d: f64 = (0..6).map(|k| qi[k] * v[k]).sum();
                    r[i][j] += d;
                    for k in 0..6 {
                        v[k] -= d * qi[k];
                    }
                }
            }
            let size = (0..6).map(|k| v[k] * v[k]).sum::<f64>().sqrt();
            if size <= 1e-9 * (0..6).map(|k| a[j][k] * a[j][k]).sum::<f64>().sqrt() {
                return Err(undetermined());
            }
            r[j][j] = size;
            q.push(v.map(|x| x / size));
        }
        let mut g = vec![[0.0; 6]; n];
        for j in 0..6 {
            let qt: Vec<f64> = q.iter().map(|qi| qi[j]).collect();
            for i in (0..n).rev() {
                let s: f64 = ((i + 1)..n).map(|m| r[i][m] * g[m][j]).sum();
                g[i][j] = (qt[i] - s) / r[i][i];
            }
        }
        g
    };
    let mo: V3 = std::array::from_fn(|k| w.moment[k] - cross(o, w.force)[k]);
    let b = scaled(&[w.force[0], w.force[1], w.force[2], mo[0], mo[1], mo[2]]);
    let lambda: Vec<f64> = g.iter().map(|gi| (0..6).map(|k| gi[k] * b[k]).sum()).collect();
    // What none of them carries: the reaction less their shares, beyond
    // what the arithmetic can leave (FLOAT of the sums the reaction is made
    // of, and of the solve). Statics makes it zero; a share of the
    // reaction's own size would let a large load across the parts hide a
    // small one along what none of them carries.
    let mut left = b;
    for (c, x) in a.iter().zip(&lambda) {
        for k in 0..6 {
            left[k] -= c[k] * x;
        }
    }
    let carried: f64 = a.iter().zip(&lambda).map(|(c, x)| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt() * x.abs()).sum();
    let (ff, fm) = w.float();
    let tol_f = ff + FLOAT * (norm(w.force) + carried);
    let tol_m = fm + norm(o) * ff + FLOAT * (norm(mo) + norm(o) * norm(w.force) + l * carried);
    let (lf, lm): (V3, V3) = ([left[0], left[1], left[2]], [left[3] * l, left[4] * l, left[5] * l]);
    if norm(lf) > tol_f || norm(lm) > tol_m {
        let local = |v: V3| -> V3 { std::array::from_fn(|k| dot(axis(&w.frame, k), v)) };
        let (f, m) = (local(w.force), local(mo));
        // What is left, by what the arrangement is: an axial force on radial
        // bearings alone; a path no joint moves; else unsaid.
        let turn_axis = free.filter(|f| f.0).map(|f| f.1);
        let axial = !parts.iter().any(|p| p.carries == Carries::Axial)
            && parts.iter().any(|p| {
                let z = axis(&p.load, 2);
                p.carries == Carries::Radial && turn_axis.is_none_or(|a| dot(a, z).abs() > 1.0 - 1e-6) && drive.is_none_or(|d| d.turn) && dot(w.force, z).abs() > tol_f
            });
        let hint = if axial {
            "the force along the bearings' axis needs a part that carries it: name the thrust bearing. The split does not assign thrust to a radial bearing"
        } else if drive.is_some_and(|d| !d.turn) {
            "on a slide a bushing carries across its bore, a roller across its axle, a carriage across its rail, and a nut along its screw: name the \
             parts that carry the rest"
        } else if drive.is_none() {
            "the path has no joint, so nothing passes a force along a motion or a torque about it, and a carriage carries nothing along its rail: ask the \
             reaction of the path whose joint moves it (its slide), or name the part that carries the rest"
        } else {
            "name the part that carries the rest"
        };
        return Err(format!(
            "{} cannot carry all of the reaction, a force ({:.6}, {:.6}, {:.6}) N and a moment ({:.6}, {:.6}, {:.6}) N mm in the path's frame about its \
             origin: {hint}",
            kinds(),
            f[0],
            f[1],
            f[2],
            m[0],
            m[1],
            m[2],
        ));
    }
    // A part that only pushes along its z, pulled.
    for (s, part) in parts.iter().enumerate().filter(|(_, p)| p.carries == Carries::Axial && !p.pulls) {
        let z = axis(&part.load, 2);
        let along: f64 =
            whose.iter().zip(&lambda).filter(|((w, _, f), _)| *w == s && *f).map(|((_, u, _), x)| dot(*u, z) * x).sum::<f64>() + part.coupled.map_or(0.0, |c| dot(c.0, z));
        if along < -tol_f {
            return Err(format!(
                "{} would pull what it carries by {:.6} N along its axis, and it pushes only, along its load frame's +z (a single-direction thrust bearing, \
                 from its housing washer toward its shaft washer; its record does not say it carries a pull, thrust double-direction): turn it over, or \
                 name the part that carries the pull",
                part.path, -along
            ));
        }
    }
    // Each part's share, a force at its load frame's origin and a moment,
    // in its given frame about that frame's origin; its rows of G, unscaled
    // to the reaction's units, so.
    let mut out = vec![[0.0; 6]; parts.len()];
    let mut map = vec![[[0.0; 6]; 6]; parts.len()];
    let thrust = whose.iter().zip(&lambda).find(|((s, ..), _)| *s == parts.len() + 1).map(|(_, x)| x.abs());
    for (i, &(s, u, is_force)) in whose.iter().enumerate() {
        if s >= parts.len() {
            continue;
        }
        let given = &parts[s].given;
        let arm: V3 = std::array::from_fn(|k| origin(&parts[s].load)[k] - origin(given)[k]);
        let m = if is_force { cross(arm, u) } else { u };
        let row: [f64; 6] = std::array::from_fn(|j| if j < 3 { g[i][j] } else { g[i][j] / l });
        for k in 0..3 {
            let (cf, cm) = (if is_force { dot(axis(given, k), u) } else { 0.0 }, dot(axis(given, k), m));
            out[s][k] += cf * lambda[i];
            out[s][3 + k] += cm * lambda[i];
            for j in 0..6 {
                map[s][k][j] += cf * row[j];
                map[s][3 + k][j] += cm * row[j];
            }
        }
    }
    // An axial part's coupling (a screw's nut: the thread is what it rates)
    // passes its force along the part's axis, at its load frame: in its
    // given frame about that frame's origin. Any other part a coupling acts
    // on (a carriage a belt is clamped to) carries what its rated contact
    // does, which the split gives: the coupling's lines are not on it.
    for (s, part) in parts.iter().enumerate().filter(|(_, p)| p.carries == Carries::Axial) {
        let Some((f, _)) = part.coupled else { continue };
        let z = axis(&part.load, 2);
        let f: V3 = z.map(|x| x * dot(f, z));
        let arm: V3 = std::array::from_fn(|k| origin(&part.load)[k] - origin(&part.given)[k]);
        let m = cross(arm, f);
        for k in 0..3 {
            out[s][k] += dot(axis(&part.given, k), f);
            out[s][3 + k] += dot(axis(&part.given, k), m);
        }
    }
    Ok((out, map, thrust))
}

/// A unit force along u at p (from the moment's origin): its six numbers.
fn force(u: V3, p: V3) -> [f64; 6] {
    let m = cross(p, u);
    [u[0], u[1], u[2], m[0], m[1], m[2]]
}

/// How many of the columns are independent (scaled as the solve's are not:
/// for a message).
fn rank(cols: &[[f64; 6]]) -> usize {
    let mut q: Vec<[f64; 6]> = Vec::new();
    for c in cols {
        let mut v = *c;
        for _ in 0..2 {
            for qi in &q {
                let d: f64 = (0..6).map(|k| qi[k] * v[k]).sum();
                for k in 0..6 {
                    v[k] -= d * qi[k];
                }
            }
        }
        let size = (0..6).map(|k| v[k] * v[k]).sum::<f64>().sqrt();
        if size > 1e-9 * (0..6).map(|k| c[k] * c[k]).sum::<f64>().sqrt().max(1e-300) {
            q.push(v.map(|x| x / size));
        }
    }
    q.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    const I: [f64; 12] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];

    fn at(x: f64, y: f64) -> [f64; 12] {
        let mut f = I;
        (f[9], f[10]) = (x, y);
        f
    }

    /// What the table's weight W (at a, b, h from its centre) and a push
    /// Fy across the rails at the same point make of the reaction on what
    /// the table carries, about its centre.
    fn load(w: f64, fy: f64, (a, b, h): (f64, f64, f64)) -> Wrench {
        let (f, r) = ([0.0, fy, w], [a, b, h]);
        Wrench::of(f, cross(r, f), I)
    }

    fn shares(at: &[[f64; 12]], w: &Wrench, rating: &str) -> Vec<[f64; 6]> {
        let parts: Vec<Part> = at.iter().map(|&f| Part { path: "c", carries: Carries::Guide, pulls: false, load: f, given: f, coupled: None, rating }).collect();
        split(w, &parts, Some(Drive { turn: false, axis: 0 })).unwrap().0
    }

    /// Four identical blocks under a rigid table, rails along x at y = +-d/2,
    /// blocks at x = +-c/2 (c = 80, d = 100): HIWIN's Table 1-3 (case 1)
    /// and THK's formulas, P = W/4 +- W a / 2c +- W b / 2d down each block
    /// and Pt = F/4 +- F a / 2c across, the blocks' moments none; two on one
    /// rail, P = W/2 +- W a / c and Pt = F/2 +- F a / c, the roll W b - F h
    /// alike on both as MR. Blocks of different records are refused.
    #[test]
    fn carriages_share_as_a_rigid_table() {
        let (c, d, (w, fy), p) = (80.0, 100.0, (200.0, 30.0), (25.0, -15.0, 40.0));
        let (a, b, h) = p;
        let four = shares(&[at(40.0, 50.0), at(-40.0, 50.0), at(-40.0, -50.0), at(40.0, -50.0)], &load(w, fy, p), "MGN12H");
        let roll = w * b - fy * h;
        for (s, (x, y)) in four.iter().zip([(1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0), (1.0, -1.0)]) {
            let down = w / 4.0 + x * w * a / (2.0 * c) + y * roll / (2.0 * d);
            let across = fy / 4.0 + x * fy * a / (2.0 * c);
            let want = [0.0, across, down, 0.0, 0.0, 0.0];
            for k in 0..6 {
                assert!((s[k] - want[k]).abs() < 1e-9 * w, "four blocks: {s:?} against {want:?}");
            }
        }
        let two = shares(&[at(40.0, 0.0), at(-40.0, 0.0)], &load(w, fy, p), "MGN12H");
        for (s, x) in two.iter().zip([1.0, -1.0]) {
            let want = [0.0, fy / 2.0 + x * fy * a / c, w / 2.0 + x * w * a / c, roll / 2.0, 0.0, 0.0];
            for k in 0..6 {
                assert!((s[k] - want[k]).abs() < 1e-9 * w * c, "two blocks on a rail: {s:?} against {want:?}");
            }
        }
        let parts: Vec<Part> = [("MGN12H", 40.0), ("MGN12C", -40.0)]
            .iter()
            .map(|&(r, x)| Part { path: "c", carries: Carries::Guide, pulls: false, load: at(x, 0.0), given: at(x, 0.0), coupled: None, rating: r })
            .collect();
        let why = split(&load(w, fy, p), &parts, Some(Drive { turn: false, axis: 0 })).err().unwrap();
        assert!(why.contains("undetermined") && why.contains("one record"), "{why}");
    }

    /// Two radial bearings on a shaft along z, 80 mm apart, under 100 N across at their middle and 20 N along
    /// the shaft: each carries 50 N across whatever the thrust. Where one of them carries thrust (a deep-groove
    /// ball bearing) the 20 N is theirs between them, in no share, and the split says its size; where neither
    /// does (bushings), or nobody asks for it, the split is refused: no thrust is assigned to a radial bearing.
    #[test]
    fn a_located_shafts_thrust_is_its_bearings_and_their_radial_loads_stand() {
        let on_axis = |z: f64| {
            let mut f = I;
            f[11] = z;
            f
        };
        let bearings = |locates: bool| -> Vec<Part<'static>> {
            [0.0, 80.0].iter().map(|&z| Part { path: "b", carries: Carries::Radial, pulls: locates, load: on_axis(z), given: on_axis(z), coupled: None, rating: "608" }).collect()
        };
        let (f, r) = ([0.0, -100.0, 20.0], [0.0, 0.0, 40.0]);
        let w = Wrench::of(f, cross(r, f), I);
        let turn = Some(Drive { turn: true, axis: 2 });
        let (shares, _, thrust) = split_located(&w, &bearings(true), turn, true).unwrap();
        for share in &shares {
            assert!((share[1] + 50.0).abs() < 1e-9 && share[0].abs() < 1e-9 && share[2].abs() < 1e-9, "{shares:?}");
        }
        assert!(thrust.is_some_and(|t| (t - 20.0).abs() < 1e-9), "{thrust:?}");
        for refused in [split_located(&w, &bearings(false), turn, true).err(), split(&w, &bearings(true), turn).err()] {
            assert!(refused.is_some_and(|why| why.contains("name the thrust bearing")));
        }
        // With no thrust there is none to share, and none is reported beyond rounding.
        let across = Wrench::of([0.0, -100.0, 0.0], cross(r, [0.0, -100.0, 0.0]), I);
        assert!(split_located(&across, &bearings(true), turn, true).unwrap().2.is_some_and(|t| t.abs() < 1e-9));
    }
}
