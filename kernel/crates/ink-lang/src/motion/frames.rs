use super::*;

/// The frame a shape of a pair claim rests in: the placement of a shape
/// placed in `a` (an instance's, a construction's in it, a zone's: the
/// first in node order), else in `b`, that moves over `span`. A claim about two shapes is the same claim of them
/// both moved rigidly at each moment (disjoint, a clearance, a subset,
/// overlaps, touches: distances and depths are kept), so it is asked of
/// them in that frame: `a` rests there and does not move with t, `b` moves
/// only as the two move apart (a pin turning in its hole moves nothing),
/// and the search covers `a`'s own box, not the space the pair sweeps.
/// In world bounds over a piece of the motion each shape blurs by how far
/// its own points move, and their max does not see that they move
/// together; and the space the pair sweeps is searched cell by cell, the
/// contact passing through each (benchmark 32's pin-b-seated: 4.2 million
/// boxes and left to samples; in the crank's frame, proven in seconds).
///
/// Each shape placed (Doc::placements: a field f at R_x^T (p - t_x)) is f
/// placed relative to the frame, at R_a^T R_x, R_a^T (t_x - t_a), the
/// frame's own f itself; anything else read at x, y, z is read at the
/// frame's point R_a p + t_a. Exact in real numbers where R_a is a
/// rotation (its inverse is R_a^T), so a placement is taken only where it
/// is one over the span (rigid::turn_bounds: det exactly 1, r^T r exactly
/// I). A shape turning about z relative to the frame whose move is not one
/// point is placed about its pivot (`pivot`), shown to drift less than the
/// claim's accuracy `tol`, where the claim can be charged its drift
/// (`charge`). Gives the node map, the frame's point (x, y, z of the world
/// as fields of the frame's p and t: `to_world`) and that charge (0: the
/// frame is exact). None: no shape placed in `a` or `b` moves (a shape
/// placed at a moving frame of another, a belt's tooth at its shift, is no
/// instance's: it is searched where it is).
pub(super) fn moving_frame(v: &mut ModelView, claim: &Claim, span: (f64, f64), tol: Option<f64>) -> Option<(std::collections::HashMap<u32, u32>, [u32; 3], f64)> {
    let (Claim::Disjoint(a, b) | Claim::Clearance(a, b, _) | Claim::Touches(a, b) | Claim::Subset(a, b) | Claim::Overlaps(a, b)) = *claim else { return None };
    let over = Ctx { t_span: Some(span), ..v.ctx.clone() };
    // The placements of the shapes placed in `a` (Doc::placed), in node
    // order, then in `b` (a still shape against a moving one): the first
    // that moves and is an isometry over the scope (rigid::Isometry: a
    // rotation, its translation finite). Either is exact: every shape is
    // placed relative to it.
    let placements = v.doc.placements();
    let moves = |xf: &Xf| xf.r.iter().flatten().chain(&xf.t).any(|&n| v.doc.g.depends_on(n, |x| matches!(x, Var::T)));
    let first = |s: u32| {
        let within: std::collections::HashSet<u32> = v.doc.g.reachable(&[s]).into_iter().collect();
        placements.iter().filter(|(w, _)| within.contains(w)).map(|(_, p)| p.at).filter(moves).find_map(|xf| Isometry::over(&v.doc.g, &xf.r, &xf.t, &over))
    };
    let mut frames: Vec<Xf> = [a, b].into_iter().filter_map(first).map(|iso| Xf { r: *iso.r(), t: *iso.t() }).collect();
    frames.dedup();
    if frames.is_empty() {
        return None;
    }
    let read: std::collections::HashSet<u32> = v.doc.g.reachable(&[a, b]).into_iter().collect();
    // In node order, so the nodes made are the same each run.
    let mut placed: Vec<((u32, Xf), u32)> = placements.iter().filter(|(w, _)| read.contains(w)).map(|(w, p)| ((p.shape, p.at), *w)).collect();
    // Each placed relative to the frame, its entries finite over the scope:
    // one infinite at a moment places nothing there, which no relative
    // placement need say alike.
    if !placed.iter().all(|((_, xf), _)| ink_kernel::rigid::finite_over(&v.doc.g, xf.r.iter().flatten().chain(&xf.t), &over)) {
        return None;
    }
    // The world's frame too, where a shape rests in it: a shaft turning
    // about its own axis in a still block's bore moves nothing of the
    // shaft there, where in the shaft's frame the block's screws sweep
    // round it. Asked in the first frame where nothing of the claim moves
    // (reads t), else in the first.
    let g = &mut v.doc.g;
    let (z, o) = (g.constant(0.0), g.constant(1.0));
    let identity = Xf { r: [[o, z, z], [z, o, z], [z, z, o]], t: [z, z, z] };
    // A shape no instance places (a still part's own field) is placed
    // where it is, so that in another's frame it is one placed shape too:
    // a still bush round about a turning shaft is unchanged by its turn.
    for s in [a, b] {
        let within: std::collections::HashSet<u32> = g.reachable(&[s]).into_iter().collect();
        if !placed.iter().any(|(_, w)| within.contains(w)) && g.depends_on(s, |x| matches!(x, Var::X | Var::Y | Var::Z)) {
            placed.push(((s, identity), s));
        }
    }
    frames.push(identity);
    // Nothing of the claim rests in any frame where a shape's own field
    // reads t (a belt's teeth), nor in one where a shape turns about no
    // axis of its own fixed on the frame's (`Axes`): those frames, but the
    // first, are not built.
    let own = placed.iter().any(|((f, _), _)| v.doc.g.depends_on(*f, |x| matches!(x, Var::T)));
    let mut first = None;
    for (i, fa) in frames.into_iter().enumerate() {
        let rels = relative(v, fa, &placed, &over);
        let may = !own && rels.iter().all(|(_, _, rel, _)| rel.as_ref().is_none_or(|rel| !moving(&v.doc.g, rel) || Axes::of(&v.doc.g, rel, &over).is_some()));
        if i > 0 && !may {
            continue;
        }
        let m = in_frame(v, claim, fa, rels, span, tol, &over);
        let g = &mut v.doc.g;
        let still = may
            && [a, b].into_iter().all(|s| {
                let s = g.replace(s, &m.0);
                !g.depends_on(s, |x| matches!(x, Var::T))
            });
        if still {
            return Some(m);
        }
        first.get_or_insert(m);
    }
    first
}

/// A placed shape seen from a frame: the node placing it, the shape, its placement relative to the
/// frame (None: the frame's own shape, where it rests), and the frame's origin in the shape's own
/// frame where that is one point over the scope (`relative`).
type Rel = (u32, u32, Option<Xf>, Option<[u32; 3]>);

/// Each placed shape's placement relative to the frame `fa`: R_a^T R_x,
/// R_a^T (t_x - t_a), as (node placing it, shape, placement); the frame's
/// own shape where it rests (None).
fn relative(v: &mut ModelView, fa: Xf, placed: &[((u32, Xf), u32)], over: &Ctx) -> Vec<Rel> {
    let over = over.clone();
    let mut rels: Vec<Rel> = Vec::with_capacity(placed.len());
    for &((f, xf), w) in placed {
        if xf == fa {
            rels.push((w, f, None, None));
            continue;
        }
        let g = &mut v.doc.g;
        let dot = |g: &mut ink_core::graph::Graph, u: [u32; 3], w: [u32; 3]| {
            let p: Vec<u32> = (0..3).map(|k| g.bin_real(BinOp::Mul, u[k], w[k])).collect();
            let s = g.bin_real(BinOp::Add, p[0], p[1]);
            g.bin_real(BinOp::Add, s, p[2])
        };
        let col = |m: &[[u32; 3]; 3], j: usize| [m[0][j], m[1][j], m[2][j]];
        let d: [u32; 3] = std::array::from_fn(|k| g.bin_real(BinOp::Sub, xf.t[k], fa.t[k]));
        let r: [[u32; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| dot(g, col(&fa.r, i), col(&xf.r, j))));
        let t: [u32; 3] = std::array::from_fn(|i| dot(g, col(&fa.r, i), d));
        // What the frame's turn cancels, cancelled (R_a^T R_a = I).
        let entries: Vec<u32> = r.iter().flatten().copied().chain(t).collect();
        let rel = match ink_kernel::rigid::simplified(g, &entries, &over) {
            Some(e) => Xf { r: std::array::from_fn(|i| std::array::from_fn(|j| e[3 * i + j])), t: [e[9], e[10], e[11]] },
            None => Xf { r, t },
        };
        // Where the shape does not turn in the world and stands a fixed way from the frame's origin,
        // that way in its own frame: -R_x^T (t_x - t_a), every number of it one number over the scope.
        let fixed = |g: &ink_core::graph::Graph, n: u32| g.enclosed_constant(n).is_some() || ink_kernel::rigid::constant_over(g, n, &over).is_some();
        let own = xf.r.iter().flatten().chain(&d).all(|&n| fixed(g, n)).then(|| {
            let zero = g.constant(0.0);
            std::array::from_fn(|i| {
                let s = dot(g, col(&xf.r, i), d);
                g.bin_real(BinOp::Sub, zero, s)
            })
        });
        rels.push((w, f, Some(rel), own));
    }
    rels
}

/// Whether a placement reads t.
fn moving(g: &ink_core::graph::Graph, xf: &Xf) -> bool {
    xf.r.iter().flatten().chain(&xf.t).any(|&n| g.depends_on(n, |x| matches!(x, Var::T)))
}

/// `moving_frame` in frame `fa` (an isometry over the scope `over`), of
/// the shapes placed in the claim (`placed`: each shape and its placement,
/// with the node placing it).
fn in_frame(v: &mut ModelView, claim: &Claim, fa: Xf, rels: Vec<Rel>, span: (f64, f64), tol: Option<f64>, over: &Ctx) -> (std::collections::HashMap<u32, u32>, [u32; 3], f64) {
    let over = over.clone();
    let moves = |v: &ModelView, xf: &Xf| moving(&v.doc.g, xf);
    // A proper turn about the shape's own z (an isometry, r_zz exactly 1):
    // its parts every such turn leaves unchanged (a pin, a bore about the
    // axis) placed by the move alone, equal in real numbers
    // (`split_about`). Where its move is not one point over the scope,
    // the turn about a pivot it holds still to within a certified drift
    // (`pivot`), where the claim can be charged it (`charge`).
    // A turn about the shape's own x or y fixed on one of the frame's axes
    // (a shaft along y turned by a rot-y) is the same turn about z seen
    // through relabelled axes (`Axes`): each shape so turning is asked as
    // its own (`seen`: the shape and its placement relabelled), and what is
    // made of it relabelled back.
    let about_z = |v: &ModelView, rel: &Xf| Isometry::over(&v.doc.g, &rel.r, &rel.t, &over).filter(|_| ink_kernel::rigid::constant_over(&v.doc.g, rel.r[2][2], &over) == Some(1.0));
    let mut seen: std::collections::HashMap<u32, (Axes, u32, Xf, Xf)> = std::collections::HashMap::new();
    for (w, f, rel, _) in &rels {
        let Some(rel) = rel.filter(|rel| moves(v, rel)) else { continue };
        let Some(axes) = Axes::of(&v.doc.g, &rel, &over) else { continue };
        let (f, rel) = axes.seen(&mut v.doc.g, *f, &rel);
        if let Some(turn) = about_z(v, &rel) {
            seen.insert(*w, (axes, f, rel, Xf { r: *turn.r(), t: *turn.t() }));
        }
    }
    // A shape that does not turn in the world, a fixed way from the frame's origin, is the frame's turn
    // seen backwards about one point of its own (`about_own`): exact, so no pivot is fitted to it.
    let exact: std::collections::HashMap<u32, [u32; 3]> = rels
        .iter()
        .filter_map(|(w, _, rel, own)| {
            let rel = rel.as_ref().filter(|rel| moving(&v.doc.g, rel))?;
            let unfixed = !rel.t.iter().all(|&n| v.doc.g.enclosed_constant(n).is_some() || ink_kernel::rigid::constant_over(&v.doc.g, n, &over).is_some());
            own.filter(|_| unfixed).map(|k| (*w, k))
        })
        .collect();
    let mut pivots: std::collections::HashMap<u32, Pivot> = std::collections::HashMap::new();
    for (w, (_, f, rel, _)) in &seen {
        if !exact.contains_key(w)
            && !rel.t.iter().all(|&n| v.doc.g.enclosed_constant(n).is_some() || ink_kernel::rigid::constant_over(&v.doc.g, n, &over).is_some())
            && let Some(p) = tol.and_then(|tol| pivot(v, *f, rel, span, tol))
        {
            pivots.insert(*w, p);
        }
    }
    let charge = match pivots.is_empty() {
        true => 0.0,
        false => {
            let moved: std::collections::HashMap<u32, f64> = pivots.iter().map(|(w, p)| (*w, p.delta)).collect();
            match self::charge(&v.doc.g, claim, &moved) {
                Some(c) => c,
                None => {
                    pivots.clear();
                    0.0
                }
            }
        }
    };
    let mut map = std::collections::HashMap::new();
    for (w, f, rel, _) in rels {
        let to = match rel {
            None => f,
            Some(rel) if exact.contains_key(&w) => about_own(v, f, &rel, exact[&w], seen.get(&w).is_some_and(|s| s.0.identity())),
            Some(rel) => match (seen.get(&w).cloned(), pivots.get(&w)) {
                (Some((axes, f, rel, _)), Some(p)) => {
                    let g = &mut v.doc.g;
                    // The turn about the pivot: p -> R (q - l) + c.
                    let l = [g.constant(p.l[0]), g.constant(p.l[1])];
                    let t: [u32; 3] = std::array::from_fn(|k| {
                        let (a, b) = (g.bin_real(BinOp::Mul, rel.r[k][0], l[0]), g.bin_real(BinOp::Mul, rel.r[k][1], l[1]));
                        let rl = g.bin_real(BinOp::Add, a, b);
                        let c = g.constant(p.c[k]);
                        g.bin_real(BinOp::Sub, c, rl)
                    });
                    let turn = Xf { r: rel.r, t };
                    // c - l as a real number (Graph::real), not folded to
                    // its float.
                    let moved: [u32; 3] = std::array::from_fn(|k| {
                        let c = g.constant(p.c[k]);
                        if k < 2 {
                            let c = g.real(c);
                            g.sub(c, l[k])
                        } else {
                            c
                        }
                    });
                    let placed = split_about(v, f, &turn, p.l, moved).unwrap_or_else(|| v.doc.warped(f, &turn));
                    axes.back(&mut v.doc.g, placed)
                }
                (Some((axes, seen_f, _, turn)), None) => match split_about(v, seen_f, &turn, [0.0; 2], turn.t) {
                    Some(placed) => axes.back(&mut v.doc.g, placed),
                    None => v.doc.warped(f, &rel),
                },
                (None, _) => v.doc.warped(f, &rel),
            },
        };
        map.insert(w, to);
    }
    // Anything else in the world at the frame's point: R_a p + t_a.
    let g = &mut v.doc.g;
    let p = [Var::X, Var::Y, Var::Z].map(|x| g.var(x));
    let world: [u32; 3] = std::array::from_fn(|i| {
        let terms: Vec<u32> = (0..3).map(|k| g.mul(fa.r[i][k], p[k])).collect();
        let s = g.add(terms[0], terms[1]);
        let s = g.add(s, terms[2]);
        g.add(s, fa.t[i])
    });
    for i in 0..3 {
        map.entry(p[i]).or_insert(world[i]);
    }
    (map, world, charge)
}

/// A placement relative to a frame (R, T) whose turn keeps an axis of the
/// shape's own frame, k, on an axis of the frame's, j (R e_k = s e_j, s
/// +-1, exactly over the scope), seen as a turn about z: the shape's axes
/// relabelled by P (a cycle, P e_z = e_k) and the frame's by M (a signed
/// cycle, M e_z = s e_j, det 1). The shape f placed, f(R^T (p - T)), is f'
/// = f o P placed by (M^T R P, M^T T), read at M^T p: a turn about z,
/// exactly (its r_zz is s s), and relabelling is a substitution of x, y
/// and z by signed x, y and z, so every value is the same in real numbers.
#[derive(Clone, Copy)]
struct Axes {
    own: [[f64; 3]; 3],
    frame: [[f64; 3]; 3],
}

impl Axes {
    /// The relabelling of `rel`, its own z first: None where no column of
    /// its turn is a signed axis of the frame's over the scope.
    fn of(g: &ink_core::graph::Graph, rel: &Xf, over: &Ctx) -> Option<Axes> {
        // The cycle taking z to axis k: (P q)_i = q_((i + 2 - k) mod 3).
        let cycle = |k: usize| -> [[f64; 3]; 3] { std::array::from_fn(|i| std::array::from_fn(|j| if j == (i + 2 - k) % 3 { 1.0 } else { 0.0 })) };
        [2, 0, 1].into_iter().find_map(|k| {
            let col: Vec<f64> = (0..3).map(|i| ink_kernel::rigid::constant_over(g, rel.r[i][k], over)).collect::<Option<_>>()?;
            let j = (0..3).find(|&j| col[j].abs() == 1.0 && (0..3).all(|i| i == j || col[i] == 0.0))?;
            let s = col[j];
            // M = P_j diag(1, s, s): M e_z = s e_j, det s^2.
            let frame = cycle(j).map(|row| std::array::from_fn(|c| row[c] * if c == 0 { 1.0 } else { s }));
            Some(Axes { own: cycle(k), frame })
        })
    }

    fn identity(&self) -> bool {
        let i: [[f64; 3]; 3] = std::array::from_fn(|a| std::array::from_fn(|b| if a == b { 1.0 } else { 0.0 }));
        self.own == i && self.frame == i
    }

    /// The shape and its placement seen with z the axis: f o P and
    /// (M^T R P, M^T T).
    fn seen(&self, g: &mut ink_core::graph::Graph, f: u32, rel: &Xf) -> (u32, Xf) {
        if self.identity() {
            return (f, *rel);
        }
        let (p, m) = (self.own, self.frame);
        // The one term of a signed permutation's product: sign times node.
        let term = |g: &mut ink_core::graph::Graph, terms: Vec<(f64, u32)>| {
            let (s, n) = terms.into_iter().find(|(s, _)| *s != 0.0).expect("a signed permutation has one term a row");
            if s > 0.0 { n } else { g.neg(n) }
        };
        // (M^T R P)_ab = sum_i,c M_ia R_ic P_cb.
        let r = std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                let terms = (0..3).flat_map(|i| (0..3).map(move |c| (i, c))).map(|(i, c)| (m[i][a] * p[c][b], rel.r[i][c])).collect();
                term(g, terms)
            })
        });
        let t = std::array::from_fn(|a| term(g, (0..3).map(|i| (m[i][a], rel.t[i])).collect()));
        (relabel(g, f, &p), Xf { r, t })
    }

    /// A field made at the relabelled point read at the frame's: q = M^T p.
    fn back(&self, g: &mut ink_core::graph::Graph, f: u32) -> u32 {
        if self.identity() {
            return f;
        }
        let mt = std::array::from_fn(|a| std::array::from_fn(|b| self.frame[b][a]));
        relabel(g, f, &mt)
    }
}
