//! Sparse manifold dual contouring.
//!
//! The octree finds bricks of `B^3` cells that straddle the surface. Each brick
//! evaluates its corner grid in SIMD batches. Within each sign-changing cell,
//! edge crossings are grouped into surface patches (connected through the
//! cell's faces; ambiguous faces are resolved by the field at the face centre,
//! which both neighbouring cells compute identically). Each patch gets one
//! vertex placed by a quadratic error function (QEF) over surface normals, which
//! keeps sharp edges and corners. One quad is emitted per sign-changing grid
//! edge, using the patch in each of its four cells that owns that edge. One
//! vertex per patch (not per cell) keeps thin walls and sheets manifold.
//! Where sheets come within a cell of each other, a vertex the QEF placed
//! can make triangles pass through others; each vertex of those goes to its
//! crossings' mean inside its cell instead (`untangle`).
//!
//! What sampling at `h` cannot see is reported, not hidden: `mesh_checked`
//! takes every cell whose corners agree in sign and proves, by the field's
//! range over it (interval arithmetic), that it holds that sign throughout;
//! the cells it does not prove are searched, all at once level by level:
//! split in eight, the parts' corners sampled and their ranges taken, again
//! in the parts not proven, while each level costs within what the search
//! is allowed. A cell where a sample has the other sign holds material or a
//! gap the mesh leaves out (and how thin it is, is measured); one neither
//! found nor proven is undecided, down to the finest samples taken
//! (`Unresolved::searched`). Every such cell is reported. A mesh whose grid
//! would take more bricks' samples than it is allowed is not made (`Limit`).
//!
//! Made to an accuracy (`mesh_whole`), a mesh is not sampled where bounds
//! show the surface a height field flat to within half of it: slabs of
//! blocks drawn whole (`Flat`), joined into regions, meshed only round
//! their rims and triangulated inside (`disc`). The cells inside are shown
//! to hold nothing but that surface, so none is open.

mod analysis;
mod brick;
mod build;
mod contour;
mod flat;
mod proof;
mod report;
mod search;
mod solid;

#[cfg(test)]
mod tests;

pub use analysis::{Analysis, Mesh, Shell, analyze, find};
pub use build::{mesh, mesh_checked, mesh_in_time, mesh_whole};
pub use flat::Even;
pub use report::{BRICK_TRIANGLES, Limit, MAX_REGIONS, Region, Unresolved};
pub use solid::crossing;
