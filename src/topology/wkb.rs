//! WKB codec for every geometry type, decoding into [`Shape`].
use super::shape::{Kind, Shape};
use crate::{
    geom::Point,
    wkb::{
        WkbError,
        codec::{Reader, point_bytes, put_header, put_len, put_point},
    },
};

/// Collections nested deeper than this are rejected to bound recursion.
const MAX_DEPTH: usize = 256;
/// Smallest encoded size of a ring count.
const MIN_RING_BYTES: usize = 4;
/// Smallest encoded size of a collection member (byte order plus type code).
const MIN_MEMBER_BYTES: usize = 5;

/// Decodes any OGC geometry. Rings decode as [`Kind::LinearRing`] polygon parts.
pub(crate) fn read(bytes: &[u8]) -> Result<Shape, WkbError> {
    let mut reader = Reader::new(bytes);
    let shape = read_shape(&mut reader, 0)?;
    reader.finish()?;
    Ok(shape)
}

fn read_shape(reader: &mut Reader<'_>, depth: usize) -> Result<Shape, WkbError> {
    if depth >= MAX_DEPTH {
        return Err(WkbError::NestingTooDeep);
    }
    let header = reader.header()?;
    let kind = Kind::from_wkb_code(header.type_id).ok_or(WkbError::UnsupportedType)?;
    reader.skip_srid(header)?;
    let dims = header.dims;
    let points = |reader: &mut Reader<'_>, n: usize| -> Result<Vec<Point>, WkbError> {
        (0..n).map(|_| reader.point(dims)).collect()
    };
    let mut g = Shape::empty(kind, dims);
    match kind {
        Kind::Point => {
            let p = reader.point(dims)?;
            if !(p.x.is_nan() && p.y.is_nan()) {
                g.points.push(p);
            }
        }
        Kind::LineString | Kind::LinearRing => {
            let n = reader.count(point_bytes(dims), WkbError::InvalidCount)?;
            if n == 1 {
                return Err(WkbError::ShortLine);
            }
            g.points = points(reader, n)?;
        }
        Kind::Polygon => {
            for _ in 0..reader.count(MIN_RING_BYTES, WkbError::InvalidCount)? {
                let n = reader.count(point_bytes(dims), WkbError::InvalidCount)?;
                let mut ring = Shape::empty(Kind::LinearRing, dims);
                ring.points = points(reader, n)?;
                if !(n == 0 || (n >= 3 && ring.points.first() == ring.points.last())) {
                    return Err(WkbError::UnclosedRing);
                }
                g.parts.push(ring);
            }
        }
        Kind::MultiPoint
        | Kind::MultiLineString
        | Kind::MultiPolygon
        | Kind::GeometryCollection => {
            for _ in 0..reader.count(MIN_MEMBER_BYTES, WkbError::InvalidCount)? {
                let child = read_shape(reader, depth + 1)?;
                if kind.member().is_some_and(|member| member != child.kind) {
                    return Err(WkbError::InvalidChildType);
                }
                g.parts.push(child);
            }
        }
    }
    Ok(g)
}

/// Encodes as little-endian WKB; rings are written as `LineString`s.
///
/// # Panics
///
/// Panics if a count exceeds `u32::MAX`, which WKB cannot represent.
pub(crate) fn write(g: &Shape) -> Vec<u8> {
    let mut out = Vec::new();
    write_shape(&mut out, g);
    out
}

fn write_shape(out: &mut Vec<u8>, g: &Shape) {
    let dims = g.dimensions();
    put_header(out, g.kind.wkb_code(), dims);
    match g.kind {
        Kind::Point => {
            let empty = Point {
                x: f64::NAN,
                y: f64::NAN,
                ..Point::default()
            };
            put_point(out, g.points.first().copied().unwrap_or(empty), dims);
        }
        Kind::LineString | Kind::LinearRing => {
            put_len(out, g.points.len());
            for &p in &g.points {
                put_point(out, p, dims);
            }
        }
        Kind::Polygon => {
            put_len(out, g.parts.len());
            for ring in &g.parts {
                put_len(out, ring.points.len());
                for &p in &ring.points {
                    put_point(out, p, dims);
                }
            }
        }
        Kind::MultiPoint
        | Kind::MultiLineString
        | Kind::MultiPolygon
        | Kind::GeometryCollection => {
            put_len(out, g.parts.len());
            for p in &g.parts {
                write_shape(out, p);
            }
        }
    }
}
