//! GEOS 3.13.1 `TopologyPreservingSimplifier` for every geometry type.
//!
//! Only the Python bindings use it, so it builds with the `python` feature.
mod segment;
mod shape;
mod simplifier;
pub(crate) mod wkb;

pub(crate) use shape::{Kind, Shape};
pub(crate) use simplifier::simplify;
