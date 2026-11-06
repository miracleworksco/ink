//! ink kernel services: general queries over the core's fields. The one
//! search, the octree, queries (bounds, integrals, emptiness, clearance, rays, slices,
//! moments), proofs over ranged params (time among them), rigid
//! placements and mass properties, points sampled on the surface and wall
//! thickness, meshing, the disk cache of results, and the call's time
//! limit.
//! Depends only on ink-core (docs/CORE.md).

pub mod cache;
pub mod centred;
pub mod feature;
pub mod mesh;
pub mod motion;
pub mod octree;
pub mod progress;
pub mod query;
pub mod ranged;
pub mod rigid;
pub mod runners;
pub mod scope;
pub mod search;
pub mod stop;
#[cfg(test)]
mod test_fields;
