use std::{error::Error, fmt};

/// Why a WKB buffer was rejected.
///
/// The Python bindings surface these messages verbatim, so they are stable.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum WkbError {
    /// The buffer ended inside a value.
    Truncated,
    /// The byte-order marker was neither 0 nor 1.
    InvalidByteOrder,
    /// The geometry was not a `Polygon`, or a `MultiPolygon` of them.
    NotPolygonal,
    /// The geometry type code is not one of the seven OGC types.
    UnsupportedType,
    /// An element count exceeds what the remaining bytes could hold.
    InvalidCount,
    /// A ring's point count exceeds what the remaining bytes could hold.
    InvalidRingSize,
    /// An X or Y ordinate was NaN or infinite.
    NonFiniteXy,
    /// A ring was not closed or had fewer than three points.
    UnclosedRing,
    /// A line had exactly one point.
    ShortLine,
    /// A multi-geometry contained a member of the wrong type.
    InvalidChildType,
    /// Collections were nested too deeply.
    NestingTooDeep,
    /// Bytes remained after the geometry.
    TrailingBytes,
}

impl fmt::Display for WkbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Truncated => "truncated WKB",
            Self::InvalidByteOrder => "invalid WKB endian",
            Self::NotPolygonal => "expected Polygon/MultiPolygon WKB",
            Self::UnsupportedType => "unsupported WKB type",
            Self::InvalidCount => "invalid WKB count",
            Self::InvalidRingSize => "invalid ring size",
            Self::NonFiniteXy => "nonfinite XY coordinate",
            Self::UnclosedRing => "unclosed or short ring",
            Self::ShortLine => "line must be empty or contain at least two points",
            Self::InvalidChildType => "invalid child type in WKB multi-geometry",
            Self::NestingTooDeep => "WKB nesting exceeds 256",
            Self::TrailingBytes => "trailing WKB",
        })
    }
}

impl Error for WkbError {}
