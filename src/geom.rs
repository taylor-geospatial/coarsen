//! Polygonal geometry storage and bounds. The hot path never constructs GEOS objects.
use std::{
    cmp::Ordering,
    ops::{BitOr, BitOrAssign},
    slice,
};

/// A coordinate. Equality compares only X and Y, as GEOS does.
///
/// Missing Z and M ordinates are NaN.
#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct Point {
    /// X ordinate.
    pub x: f64,
    /// Y ordinate.
    pub y: f64,
    /// Z ordinate, or NaN when absent.
    pub z: f64,
    /// M ordinate, or NaN when absent.
    pub m: f64,
}
impl Default for Point {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: f64::NAN,
            m: f64::NAN,
        }
    }
}
impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}
impl Point {
    /// A hashable XY key; `-0.0` and `0.0` share a key, matching `==`.
    pub(crate) fn key(self) -> (u64, u64) {
        (
            if self.x == 0.0 { 0 } else { self.x.to_bits() },
            if self.y == 0.0 { 0 } else { self.y.to_bits() },
        )
    }
    /// Lexicographic XY order, as GEOS `Coordinate::compareTo`.
    ///
    /// # Panics
    ///
    /// Panics if an ordinate is NaN; WKB decoding rejects nonfinite XY.
    pub(crate) fn cmp_xy(self, other: Self) -> Ordering {
        let order = |a: f64, b: f64| a.partial_cmp(&b).expect("XY ordinates are finite");
        order(self.x, other.x).then(order(self.y, other.y))
    }
}

/// A closed coordinate sequence.
pub type Ring = Vec<Point>;

/// EWKB coordinate-dimension flags, stored in the high bits of the type code.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub(crate) struct Dims(u32);
impl Dims {
    pub(crate) const XY: Self = Self(0);
    pub(crate) const Z: Self = Self(0x8000_0000);
    pub(crate) const M: Self = Self(0x4000_0000);

    /// Keeps only the Z and M flags of an EWKB type code.
    pub(crate) const fn from_ewkb(code: u32) -> Self {
        Self(code & (Self::Z.0 | Self::M.0))
    }
    pub(crate) const fn bits(self) -> u32 {
        self.0
    }
    pub(crate) const fn has_z(self) -> bool {
        self.0 & Self::Z.0 != 0
    }
    pub(crate) const fn has_m(self) -> bool {
        self.0 & Self::M.0 != 0
    }
}
impl BitOr for Dims {
    type Output = Self;
    fn bitor(self, other: Self) -> Self {
        Self(self.0 | other.0)
    }
}
impl BitOrAssign for Dims {
    fn bitor_assign(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

/// A polygon: a shell followed by holes. Empty polygons have no rings.
#[derive(Clone, Debug)]
pub struct Polygon {
    pub(crate) rings: Vec<Ring>,
    /// WKB dimension flags, retained even when the polygon is empty.
    pub(crate) dims: Dims,
}
impl Polygon {
    /// Returns the shell followed by the holes.
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }
}

/// A `Polygon` or `MultiPolygon`, as decoded by [`crate::wkb::read`].
#[derive(Clone, Debug)]
pub enum Geometry {
    /// A single polygon.
    Polygon(Polygon),
    /// Zero or more polygons.
    MultiPolygon(Vec<Polygon>),
}
impl Geometry {
    /// Returns the member polygons; a [`Geometry::Polygon`] has exactly one.
    pub fn polygons(&self) -> &[Polygon] {
        match self {
            Self::Polygon(p) => slice::from_ref(p),
            Self::MultiPolygon(ps) => ps,
        }
    }
    pub(crate) fn polygons_mut(&mut self) -> &mut [Polygon] {
        match self {
            Self::Polygon(p) => slice::from_mut(p),
            Self::MultiPolygon(ps) => ps,
        }
    }
    pub(crate) fn rings(&self) -> impl Iterator<Item = &Ring> {
        self.polygons().iter().flat_map(|p| &p.rings)
    }
    pub(crate) fn bounds(&self) -> Bounds {
        Bounds::from_points(self.rings().flatten().copied())
    }
}

/// An axis-aligned XY envelope. The default is empty.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Bounds {
    pub(crate) min: [f64; 2],
    pub(crate) max: [f64; 2],
}
impl Default for Bounds {
    fn default() -> Self {
        Self {
            min: [f64::INFINITY; 2],
            max: [f64::NEG_INFINITY; 2],
        }
    }
}
impl Bounds {
    pub(crate) fn from_points(pts: impl IntoIterator<Item = Point>) -> Self {
        let mut b = Self::default();
        for p in pts {
            b.add(p);
        }
        b
    }
    pub(crate) fn add(&mut self, p: Point) {
        self.min[0] = self.min[0].min(p.x);
        self.min[1] = self.min[1].min(p.y);
        self.max[0] = self.max[0].max(p.x);
        self.max[1] = self.max[1].max(p.y);
    }
    pub(crate) fn union(&mut self, b: Self) {
        if !b.is_empty() {
            self.add(Point {
                x: b.min[0],
                y: b.min[1],
                ..Point::default()
            });
            self.add(Point {
                x: b.max[0],
                y: b.max[1],
                ..Point::default()
            });
        }
    }
    pub(crate) fn is_empty(self) -> bool {
        self.min[0] > self.max[0]
    }
    pub(crate) fn intersects(self, b: Self) -> bool {
        self.min[0] <= b.max[0]
            && self.max[0] >= b.min[0]
            && self.min[1] <= b.max[1]
            && self.max[1] >= b.min[1]
    }
    pub(crate) fn contains(self, p: Point) -> bool {
        p.x >= self.min[0] && p.x <= self.max[0] && p.y >= self.min[1] && p.y <= self.max[1]
    }
}
