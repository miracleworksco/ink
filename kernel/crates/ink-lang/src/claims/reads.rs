use super::*;

/// The bodies at a failure's witness, for each shape the claim compares (its first, then its
/// second): each placed body (by the path of the instance placing it) that shape is made of whose
/// own shape holds the witness there, at its moment and params; for a clearance, the nearest
/// body to each point of its pair. Which bodies meet, where the point itself moves whenever the
/// proof does. None where it does not fail at a point.
pub(super) fn bodies_at(v: &ModelView, operands: &[u32], clearance: bool, check: &Check) -> Vec<Vec<String>> {
    if check.status != Status::Fails || operands.is_empty() {
        return Vec::new();
    }
    let (pts, t, params): ([V3; 2], f64, &[(u32, f64)]) = match &check.evidence {
        Evidence::Empty(Proof { verdict: Verdict::Disproven { witness, t, params, .. }, .. }) => ([*witness; 2], *t, params),
        Evidence::ClearanceSearch { verdict: Verdict::Disproven { t, params, .. }, pair: Some((a, b)), .. } => ([*a, *b], *t, params),
        Evidence::Point { at, .. } => ([*at; 2], v.ctx.t, &[]),
        _ => return Vec::new(),
    };
    let mut c = Ctx { t, ..v.ctx.clone() };
    for &(k, x) in params {
        c.params[k as usize] = x;
    }
    let c = c.point();
    let g = &v.doc.g;
    operands
        .iter()
        .enumerate()
        .map(|(i, &f)| {
            let read: Vec<(f64, &String)> = g
                .reachable(&[f])
                .into_iter()
                .filter_map(|n| v.doc.placed_labels.get(&n).map(|l| (n, l)))
                .map(|(n, l)| {
                    let mut out = Vec::new();
                    eval_points(&Tape::new(g, n), &c, &[pts[i]], &mut out);
                    (out[0], l)
                })
                .collect();
            // A clearance's pair lies on its shapes' surfaces: the nearest body each side.
            let near = read.iter().map(|x| x.0).fold(f64::INFINITY, f64::min);
            let at = if clearance { near.max(0.0) } else { 0.0 };
            let mut names: Vec<String> = read.into_iter().filter(|x| x.0 <= at).map(|x| x.1.clone()).collect();
            names.sort();
            names.dedup();
            names
        })
        .collect()
}

/// A claim's terms at `c`: what a solve on t sets taken as it sets it
/// at c's moment (an actuator's effort at a witness), the rest at c.
pub(super) fn terms_at(v: &ModelView, terms: &[(String, u32)], mut c: Ctx) -> Vec<(String, f64)> {
    let nodes: Vec<u32> = terms.iter().map(|t| t.1).collect();
    let moving: Vec<usize> = v.doc.solves_read(&nodes).into_iter().filter(|&i| v.doc.solves[i].on_t).collect();
    if !moving.is_empty()
        && c.t != v.ctx.t
        && let Ok((at, _)) = v.doc.ctx_at(&v.overrides, c.t)
    {
        for k in moving.iter().flat_map(|&i| v.doc.solves[i].owns.clone()) {
            c.params[k] = at.params[k];
        }
    }
    terms.iter().map(|(n, f)| (n.clone(), v.doc.scalar(*f, &c))).collect()
}

/// The ranged params (in the view's scope) these fields read, with their spans.
pub(super) fn ranged_reads(v: &ModelView, nodes: &[u32]) -> Vec<(u32, (f64, f64))> {
    v.ranged_reads(nodes)
}

/// Doc::at_one_value for the ranged params these fields read.
pub(super) fn at_one_value(v: &ModelView, nodes: &[u32]) -> String {
    v.doc.at_one_value(&ranged_reads(v, nodes).into_iter().map(|(k, _)| k as usize).collect::<Vec<_>>(), &v.ctx.p_spans)
}

/// Whether these fields read a param the view ranges over. A method
/// evaluated at one point cannot cover such a scope.
pub(super) fn scope_reads(v: &ModelView, nodes: &[u32]) -> bool {
    !ranged_reads(v, nodes).is_empty()
}

/// The one moment a question reading `roots` is asked at, though it reads
/// t (`reads_t`): no span covers it and the run asked for none
/// (ModelView::unspanned_moment). The one rule a claim's and an interfere
/// pair's answer record (Check::one_moment, Pair::one_moment), and `ink
/// check` warns by.
pub fn one_moment_of(v: &ModelView, roots: &[u32]) -> Option<f64> {
    v.unspanned_moment().filter(|_| model_reads_t(v) && reads_t(v, roots)).map(|_| v.ctx.t)
}

/// Whether the model reads t anywhere (a t in its graph, a solve on t):
/// where it does not, no question does.
pub fn model_reads_t(v: &ModelView) -> bool {
    v.doc.g.nodes.iter().any(|n| matches!(n, ink_core::graph::Node::Var(Var::T))) || v.doc.solves.iter().any(|s| s.on_t)
}

/// The claims `ink require` asks (every one, imported or not) that it would
/// ask at one moment though they read t, in a view scoped as assert scopes
/// it: each over its file's span (ModelView::motion_of), at the view's
/// params, its bound included, by `one_moment_of`. Named, for `ink check`.
pub fn claims_at_one_moment(v: &mut ModelView) -> Result<Vec<String>, InkError> {
    let list: Vec<(String, Claim, Option<u32>, crate::Scope)> = v.doc.asserts.iter().map(|a| (a.name.clone(), a.claim.clone(), a.bound, a.scope)).collect();
    let mut out = Vec::new();
    for (name, claim, bound, file) in list {
        v.motion_of(file)?;
        let c = v.doc.claim_at(&claim, bound, &v.ctx);
        let roots = v.doc.claim_motion_inputs(&c, bound);
        if one_moment_of(v, &roots).is_some() {
            out.push(format!("claim {name}"));
        }
    }
    Ok(out)
}

/// Time may enter through a provider's param as well as through the graph.
pub fn reads_t(v: &ModelView, nodes: &[u32]) -> bool {
    nodes.iter().any(|&f| v.doc.g.depends_on(f, |x| x == Var::T)) || v.doc.solves_read(nodes).iter().any(|&i| v.doc.solves[i].on_t)
}

/// Whether fields `nodes` change with t over `span`: through a solve that
/// reads t, or by a rate in t not 0 somewhere (bounded over all of space,
/// the span and every value of the ranged params). One reading t only
/// through a term that is 0 at these params (a motion times a param set
/// to 0) does not: a claim on it is the claim at any moment. One run of
/// each field (motion::still): a field switched by t moves.
pub fn moves(v: &ModelView, nodes: &[u32], span: (f64, f64)) -> bool {
    v.doc.solves_read(nodes).iter().any(|&i| v.doc.solves[i].on_t)
        || nodes.iter().filter(|&&f| v.doc.g.depends_on(f, |x| x == Var::T)).any(|&f| !ink_kernel::motion::still(&Tape::new(&v.doc.g, f), &v.ctx, [I::ENTIRE; 3], span))
}

/// A proof's answer about a shape that may not be proven bounded: a proof
/// in the window is unknown for all of it (a witness there still counts).
pub(crate) fn in_window(status: Status, whole: bool, evidence: Evidence) -> Check {
    let mut c = checked(status, evidence);
    if !whole && status != Status::Fails {
        if status == Status::Holds {
            (c.status, c.in_window) = (Status::Unknown, Some(status));
        }
        c.unbounded = true;
    }
    c
}
