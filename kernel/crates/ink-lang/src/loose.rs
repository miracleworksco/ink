//! Parts that touch nothing: a target's items fall into groups, two items
//! in one where they come within a gap of each other (a `clearance` of
//! the gap that does not hold) or a relation of the target joins them
//! (a contact, a mesh, a belt). Each group but the largest is loose: no
//! part of it is within the gap of a part outside it, proven, and nothing
//! the file declares joins them. What holds a part is not seen: a part
//! resting where its screw was still touches.

use crate::claims::{self, ClaimOptions, Status};
use crate::view::ModelView;
use crate::{Claim, InkError, Item};
use std::collections::HashMap;

/// The groups of `items` (labels) that touch nothing of the rest, largest
/// first, each its labels in the target's order: every group but the one
/// with most items (the first of those, on a tie), where there are two or
/// more. Two items touch where a clearance of `gap` between them does not
/// hold (fails, is unknown or an estimate, or cannot be asked), or a
/// relation joins them (`relations`: their ends, by instance path). Two
/// items neither of which is proven bounded touch, and an item resting on
/// inexact inputs or a refused solve is treated as touching everything so
/// it is never called loose without proof. The view is scoped as a claim's
/// question (ModelView::scope_question), at the moment it is at.
pub fn loose(v: &mut ModelView, items: &[Item], relations: &[(String, String)], gap: f64) -> Result<Vec<Vec<String>>, InkError> {
    v.scope_question()?;
    let n = items.len();
    let mut groups = Groups((0..n).collect());
    // A relation's end covers an item if it names it or one of its parents.
    let covers = |end: &str, label: &str| label == end || label.starts_with(&format!("{end}/"));
    for (a, b) in relations {
        let ends = |end: &str| (0..n).filter(|&i| covers(end, &items[i].label)).collect::<Vec<_>>();
        for i in ends(a).into_iter().chain(ends(b)).collect::<Vec<_>>().windows(2) {
            groups.join(i[0], i[1]);
        }
    }
    let certified = v.certified();
    let mut reads_inexact: HashMap<u32, bool> = HashMap::new();
    let unproven: Vec<bool> = items
        .iter()
        .map(|it| {
            let reads = v.doc.placed_reads(std::slice::from_ref(it));
            v.refusal(&reads).is_some() || reads.iter().any(|&f| *reads_inexact.entry(f).or_insert_with(|| v.doc.inexact(&[f], &certified).any()))
        })
        .collect();
    let fields: Vec<u32> = items.iter().map(|it| v.doc.placed_field(it)).collect();
    let whole: Vec<bool> = fields.iter().enumerate().map(|(i, &f)| !unproven[i] && v.whole(f)).collect();
    // Each item's box grown by the gap: a pair whose grown boxes do not meet is apart by more than it.
    let bounds: Vec<Option<_>> =
        fields.iter().enumerate().map(|(i, &f)| (!unproven[i]).then(|| v.bounds(f).map(|(lo, hi, e)| (lo.map(|x| x - gap), hi.map(|x| x + gap), e))).flatten()).collect();
    let (mut near, _, _) = crate::pairs::candidate_pairs(&bounds, &whole, &[], &[]);
    // Most overlapping first: those likely touch, and a pair already in one group is not asked.
    let overlap = |(i, j): (usize, usize)| match (bounds[i], bounds[j]) {
        (Some(a), Some(b)) => (0..3).map(|k| a.1[k].min(b.1[k]) - a.0[k].max(b.0[k])).fold(f64::INFINITY, f64::min),
        _ => f64::INFINITY,
    };
    near.sort_by(|&p, &q| overlap(q).total_cmp(&overlap(p)));
    let opts = ClaimOptions { accuracy: Some(gap / 4.0), ..Default::default() };
    for (i, j) in near {
        if groups.find(i) == groups.find(j) {
            continue;
        }
        // Apart only where proven so over all of space (a clearance covers the one of them proven bounded),
        // at every value of what they read; an unproven item or refused question is no proof.
        let apart = !unproven[i]
            && !unproven[j]
            && (whole[i] || whole[j])
            && claims::check(v, Claim::Clearance(fields[i], fields[j], gap), None, opts).is_ok_and(|c| c.status == Status::Holds && c.accuracy.is_some() && !c.unbounded);
        if !apart {
            groups.join(i, j);
        }
    }
    // An item with no solid is no part: in no group.
    let mut by: Vec<Vec<usize>> = Vec::new();
    for i in (0..n).filter(|&i| !whole[i] || bounds[i].is_some()) {
        let root = groups.find(i);
        match by.iter_mut().find(|g| groups.find(g[0]) == root) {
            Some(g) => g.push(i),
            None => by.push(vec![i]),
        }
    }
    by.sort_by_key(|g| std::cmp::Reverse(g.len()));
    Ok(by.into_iter().skip(1).map(|g| g.into_iter().map(|i| items[i].label.clone()).collect()).collect())
}

/// Union-find over item indices.
struct Groups(Vec<usize>);

impl Groups {
    fn find(&self, mut r: usize) -> usize {
        while self.0[r] != r {
            r = self.0[r];
        }
        r
    }

    fn join(&mut self, i: usize, j: usize) {
        let (a, b) = (self.find(i), self.find(j));
        self.0[a.max(b)] = a.min(b);
    }
}
