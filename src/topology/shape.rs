//! General geometry storage for the GEOS simplifiers (rings retain their identity).
use crate::geom::{Dims, Point};
use std::{error::Error, fmt};

/// OGC geometry types, plus `LinearRing`, which WKB encodes as a `LineString`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Kind {
    Point,
    LineString,
    LinearRing,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
}
impl Kind {
    pub(crate) const fn wkb_code(self) -> u32 {
        match self {
            Self::Point => 1,
            Self::LineString | Self::LinearRing => 2,
            Self::Polygon => 3,
            Self::MultiPoint => 4,
            Self::MultiLineString => 5,
            Self::MultiPolygon => 6,
            Self::GeometryCollection => 7,
        }
    }
    pub(crate) const fn from_wkb_code(code: u32) -> Option<Self> {
        Some(match code {
            1 => Self::Point,
            2 => Self::LineString,
            3 => Self::Polygon,
            4 => Self::MultiPoint,
            5 => Self::MultiLineString,
            6 => Self::MultiPolygon,
            7 => Self::GeometryCollection,
            _ => return None,
        })
    }
    pub(crate) const fn is_linear(self) -> bool {
        matches!(self, Self::LineString | Self::LinearRing)
    }
    pub(crate) const fn is_collection(self) -> bool {
        matches!(
            self,
            Self::MultiPoint
                | Self::MultiLineString
                | Self::MultiPolygon
                | Self::GeometryCollection
        )
    }
    /// The collection type that holds members of this kind.
    const fn multi(self) -> Self {
        match self {
            Self::Point => Self::MultiPoint,
            Self::LineString | Self::LinearRing => Self::MultiLineString,
            Self::Polygon => Self::MultiPolygon,
            _ => Self::GeometryCollection,
        }
    }
    /// The required member type of a homogeneous multi-geometry.
    pub(crate) const fn member(self) -> Option<Self> {
        match self {
            Self::MultiPoint => Some(Self::Point),
            Self::MultiLineString => Some(Self::LineString),
            Self::MultiPolygon => Some(Self::Polygon),
            _ => None,
        }
    }
}

/// GEOS rejects nonfinite line coordinates when inserting them into its index.
///
/// The message is GEOS's own, which the Python bindings re-raise as `GEOSException`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) struct NonFiniteBoundsError;
impl fmt::Display for NonFiniteBoundsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("IllegalArgumentException: Non-finite envelope bounds passed to index insert")
    }
}
impl Error for NonFiniteBoundsError {}

/// A geometry of any type: points for points and lines, parts for everything else.
#[derive(Clone, Debug)]
pub(crate) struct Shape {
    pub(crate) kind: Kind,
    pub(crate) dims: Dims,
    pub(crate) points: Vec<Point>,
    pub(crate) parts: Vec<Shape>,
}

impl Shape {
    pub(crate) fn empty(kind: Kind, dims: Dims) -> Self {
        Self {
            kind,
            dims,
            points: Vec::new(),
            parts: Vec::new(),
        }
    }
    pub(crate) fn is_empty(&self) -> bool {
        self.points.is_empty() && self.parts.iter().all(Self::is_empty)
    }
    /// Own dimensions for leaves, the union of member dimensions otherwise.
    pub(crate) fn dimensions(&self) -> Dims {
        if self.parts.is_empty() {
            self.dims
        } else {
            self.parts.iter().fold(Dims::XY, |d, p| d | p.dimensions())
        }
    }
    /// Checks that every line coordinate is finite, as GEOS's index insert does.
    pub(crate) fn check_indexable(&self) -> Result<(), NonFiniteBoundsError> {
        if self.kind.is_linear()
            && !self
                .points
                .iter()
                .all(|p| p.x.is_finite() && p.y.is_finite())
        {
            return Err(NonFiniteBoundsError);
        }
        self.parts.iter().try_for_each(Self::check_indexable)
    }
    /// GEOS's newly populated `CoordinateSequence` infers Z from its first point
    /// (Z when empty) and drops M.
    pub(crate) fn infer_dims(&mut self) {
        self.dims = if self.points.first().is_none_or(|p| !p.z.is_nan()) {
            Dims::Z
        } else {
            Dims::XY
        };
        for p in &mut self.points {
            p.m = f64::NAN;
        }
    }
    /// GEOS `GeometryFactory::buildGeometry`: a single part as itself, otherwise
    /// the narrowest collection, or an empty `empty_kind` when there are no parts.
    pub(crate) fn build(mut parts: Vec<Self>, empty_kind: Kind) -> Self {
        if parts.len() == 1 {
            return parts.pop().expect("length is one");
        }
        let kind = parts.first().map_or(empty_kind, |first| {
            if parts.iter().all(|p| p.kind == first.kind) {
                first.kind.multi()
            } else {
                Kind::GeometryCollection
            }
        });
        Self {
            parts,
            ..Self::empty(kind, Dims::XY)
        }
    }
}
