//! ink input and output: meshes simplified to a tolerance (and how far one
//! lies from its field), mesh export (STL, 3MF, GLB) and meshes for seeing,
//! PNG, the software renderer, and STL import as an exact distance field.
//! Depends on ink-core and ink-kernel (docs/CORE.md).

pub mod export;
pub mod gltf;
pub mod render;
pub mod seen;
pub mod simplify;
pub mod stl;
