//! Pure Rust GEOS 3.13.1 coverage and topology-preserving simplification.
//!
//! No native geometry libraries. Results match GEOS 3.13.1 byte for byte, so the
//! ports preserve GEOS encounter order and floating-point operation order.
//!
//! ```
//! use coarsen::{Boundary, coverage_simplify, wkb};
//! # fn main() -> Result<(), coarsen::wkb::WkbError> {
//! # // A unit square at `x0` with a slight bump on its top edge, as WKB.
//! # let square = |x0: f64| -> Vec<u8> {
//! #     let ring: [(f64, f64); 6] = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.01), (0.0, 1.0), (0.0, 0.0)];
//! #     let mut wkb = vec![1, 3, 0, 0, 0, 1, 0, 0, 0, 6, 0, 0, 0];
//! #     for (x, y) in ring {
//! #         wkb.extend((x0 + x).to_le_bytes());
//! #         wkb.extend(y.to_le_bytes());
//! #     }
//! #     wkb
//! # };
//! let mut coverage = vec![wkb::read(&square(0.0))?, wkb::read(&square(1.0))?];
//! coverage_simplify(&mut coverage, 0.1, Boundary::Simplify);
//! assert_eq!(coverage[0].polygons()[0].rings()[0].len(), 5);
//! # Ok(())
//! # }
//! ```
mod edges;
pub mod geom;
mod index;
mod predicates;
#[cfg(feature = "python")]
mod python;
#[cfg(feature = "python")]
mod topology;
mod tpvw;
pub mod validate;
pub mod wkb;

pub use geom::Geometry;

/// Which coverage edges [`coverage_simplify`] may move.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Boundary {
    /// Simplify every edge, including the coverage's outer boundary.
    #[default]
    Simplify,
    /// Freeze edges used by one ring; simplify only edges shared by exactly two.
    Preserve,
}

/// Simplifies the geometries as one coverage, like GEOS `CoverageSimplifier`.
///
/// Shared edges are simplified once, so neighbors stay matched. Invalid coverages
/// are simplified as GEOS would, not repaired; check with
/// [`validate::invalid_mask`] first if that matters. Runs on the current Rayon
/// pool; use a pool's `install` to control parallelism.
///
/// # Panics
///
/// Panics if the coverage has `u32::MAX` or more distinct vertices.
pub fn coverage_simplify(geometries: &mut [Geometry], tolerance: f64, boundary: Boundary) {
    let coverage = edges::Coverage::extract(geometries);
    let rings = tpvw::simplify(&coverage.edges, tolerance, boundary);
    coverage.rebuild(geometries, &rings);
}
