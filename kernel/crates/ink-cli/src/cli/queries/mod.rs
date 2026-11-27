//! Geometry queries: check, eval, bounds, mass, slice, mesh, components,
//! ray, distance, render.

mod budget;
mod check;
mod mass;
mod mesh;
mod mesh_output;
mod render;
mod slice;
mod spatial;
mod target;
#[cfg(test)]
mod tests;

pub use self::budget::{Budget, SEEN};
pub(crate) use self::check::cmd_check;
pub(crate) use self::mass::cmd_mass;
pub use self::mesh::{MeshHow, Meshed, Plan, Quality, REPORTED, drawn, looks, make_meshes, motion, paint_tapes, plan_meshes, reported, span_times, tape_of};
pub(crate) use self::mesh_output::cmd_mesh;
pub use self::mesh_output::{glb_attrs, left_out, left_out_json, mesh_answer};
pub(crate) use self::render::cmd_render;
pub(crate) use self::slice::cmd_slice;
pub(crate) use self::spatial::{cmd_bounds, cmd_components, cmd_distance, cmd_eval, cmd_ray};
