//! Enclosing fields for queries that need whole targets.

use crate::{InkError, Item, Session};
use ink::view::Local;

/// The usual cause of a body that is unexpectedly empty, for its warnings.
pub(super) const EMPTY_CAUSE: &str = "(max a b) is the intersection of a and b, (min a b) their union. Shapes that only touch don't intersect";

/// Which field can supply a whole mesh: the body with a proven rigid
/// placement, otherwise the placed field, whose own reach must hold it.
pub(super) fn bounded_local(s: &mut Session, it: &Item) -> Result<Option<Local>, InkError> {
    let local = s.local(it);
    if local.is_none() {
        let f = s.doc.placed_field(it);
        s.reach(f).require_whole(&it.label)?;
    }
    Ok(local)
}

/// A query in world coordinates can use the body's proof too. Intersect
/// with the placed enclosure only where the world field has no proof of
/// its own. Every material point is inside it, even far from the origin.
pub(super) fn whole_placed_field(s: &mut Session, it: &Item) -> Result<u32, InkError> {
    let f = s.doc.placed_field(it);
    let reach = s.reach(f);
    if !reach.whole {
        if let Some(local) = s.local(it) {
            let Some((lo, hi, _)) = s.bounds(local.shape) else { return Ok(s.doc.g.constant(1.0)) };
            let (lo, hi) = local.world_box(lo, hi);
            // A clip is zero on its faces. Keep even a material point on
            // a face of the proven box strictly inside the clip.
            let (lo, hi) = (lo.map(f64::next_down), hi.map(f64::next_up));
            if lo.into_iter().chain(hi).all(f64::is_finite) {
                let enclosed = s.doc.clip(f, lo, hi);
                s.reach(enclosed).require_whole(&it.label)?;
                return Ok(enclosed);
            }
        }
        reach.require_whole(&it.label)?;
    }
    Ok(f)
}
