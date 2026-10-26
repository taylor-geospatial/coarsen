//! Polygon/MultiPolygon codec for [`Geometry`], including Z/M and SRID.
use super::{
    WkbError,
    codec::{Header, Reader, point_bytes, put_header, put_len, put_point},
};
use crate::geom::{Dims, Geometry, Polygon};

/// OGC WKB type code of a `Polygon`.
const POLYGON: u32 = 3;
/// OGC WKB type code of a `MultiPolygon`.
const MULTI_POLYGON: u32 = 6;

/// Decodes a `Polygon` or `MultiPolygon`. Accepts either byte order and EWKB or ISO
/// dimension flags; an SRID is skipped.
///
/// # Errors
///
/// Returns [`WkbError`] if the buffer is malformed, holds another geometry type,
/// has a nonfinite X or Y, has an unclosed ring, or has trailing bytes.
pub fn read(bytes: &[u8]) -> Result<Geometry, WkbError> {
    let mut reader = Reader::new(bytes);
    let header = reader.header()?;
    let geometry = if header.type_id == MULTI_POLYGON {
        reader.skip_srid(header)?;
        let n = reader.count(4, WkbError::InvalidCount)?;
        let polygons = (0..n)
            .map(|_| {
                let header = reader.header()?;
                read_polygon(&mut reader, header)
            })
            .collect::<Result<_, _>>()?;
        Geometry::MultiPolygon(polygons)
    } else {
        Geometry::Polygon(read_polygon(&mut reader, header)?)
    };
    reader.finish()?;
    Ok(geometry)
}

fn read_polygon(reader: &mut Reader<'_>, header: Header) -> Result<Polygon, WkbError> {
    if header.type_id != POLYGON {
        return Err(WkbError::NotPolygonal);
    }
    reader.skip_srid(header)?;
    let dims = header.dims;
    let n = reader.count(4, WkbError::InvalidCount)?;
    let mut rings = Vec::with_capacity(n);
    for _ in 0..n {
        let count = reader.count(point_bytes(dims), WkbError::InvalidRingSize)?;
        let mut ring = Vec::with_capacity(count);
        for _ in 0..count {
            let p = reader.point(dims)?;
            if !(p.x.is_finite() && p.y.is_finite()) {
                return Err(WkbError::NonFiniteXy);
            }
            ring.push(p);
        }
        if !(ring.is_empty() || (ring.len() >= 3 && ring.first() == ring.last())) {
            return Err(WkbError::UnclosedRing);
        }
        rings.push(ring);
    }
    Ok(Polygon { rings, dims })
}

/// Encodes a geometry as little-endian WKB. A `MultiPolygon` header carries the
/// union of its members' dimensions.
///
/// # Panics
///
/// Panics if a count exceeds `u32::MAX`, which WKB cannot represent.
pub fn write(geometry: &Geometry) -> Vec<u8> {
    let mut bytes = Vec::new();
    match geometry {
        Geometry::Polygon(p) => write_polygon(&mut bytes, p),
        Geometry::MultiPolygon(polygons) => {
            let dims = polygons.iter().fold(Dims::XY, |d, p| d | p.dims);
            put_header(&mut bytes, MULTI_POLYGON, dims);
            put_len(&mut bytes, polygons.len());
            for p in polygons {
                write_polygon(&mut bytes, p);
            }
        }
    }
    bytes
}

fn write_polygon(bytes: &mut Vec<u8>, polygon: &Polygon) {
    put_header(bytes, POLYGON, polygon.dims);
    put_len(bytes, polygon.rings.len());
    for ring in &polygon.rings {
        put_len(bytes, ring.len());
        for &p in ring {
            put_point(bytes, p, polygon.dims);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malformed_counts_and_recursive_collections() {
        assert_eq!(
            read(&[1, 3, 0, 0, 0, 255, 255, 255, 255]).unwrap_err(),
            WkbError::InvalidCount
        );
        let mut nested = vec![1, 6, 0, 0, 0, 1, 0, 0, 0];
        nested.extend_from_slice(&[1, 6, 0, 0, 0, 0, 0, 0, 0]);
        assert_eq!(read(&nested).unwrap_err(), WkbError::NotPolygonal);
    }

    #[test]
    fn round_trips_big_endian_iso_z() {
        let mut wkb = vec![0];
        wkb.extend_from_slice(&1003u32.to_be_bytes());
        wkb.extend_from_slice(&1u32.to_be_bytes());
        wkb.extend_from_slice(&4u32.to_be_bytes());
        for (x, y) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)] {
            for v in [x, y, 7.0f64] {
                wkb.extend_from_slice(&v.to_be_bytes());
            }
        }
        let geometry = read(&wkb).unwrap();
        let [polygon] = geometry.polygons() else {
            panic!("expected one polygon, got {geometry:?}");
        };
        assert!(polygon.dims.has_z() && !polygon.dims.has_m());
        assert_eq!(
            read(&write(&geometry)).unwrap().polygons()[0].rings(),
            polygon.rings()
        );
    }
}
