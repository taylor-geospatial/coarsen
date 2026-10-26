//! Byte-level WKB primitives shared by the polygonal and general codecs.
use super::WkbError;
use crate::geom::{Dims, Point};

type Result<T> = std::result::Result<T, WkbError>;

/// EWKB flag: an SRID follows the type code.
const SRID_FLAG: u32 = 0x2000_0000;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ByteOrder {
    Big,
    Little,
}

/// A decoded geometry header. The type id is the base OGC code (1-7) with ISO
/// dimension offsets removed; validating it is the caller's job.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Header {
    pub(crate) type_id: u32,
    pub(crate) dims: Dims,
    has_srid: bool,
}

/// Bounds-checked cursor over a WKB buffer.
#[derive(Debug)]
pub(crate) struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
    order: ByteOrder,
}
impl<'a> Reader<'a> {
    pub(crate) fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            pos: 0,
            order: ByteOrder::Little,
        }
    }
    #[inline]
    fn remaining(&self) -> usize {
        self.bytes.len() - self.pos
    }
    #[inline]
    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        if N > self.remaining() {
            return Err(WkbError::Truncated);
        }
        let bytes = self.bytes[self.pos..self.pos + N]
            .try_into()
            .expect("slice length equals N");
        self.pos += N;
        Ok(bytes)
    }
    #[inline]
    fn u32(&mut self) -> Result<u32> {
        let bytes = self.take()?;
        Ok(match self.order {
            ByteOrder::Little => u32::from_le_bytes(bytes),
            ByteOrder::Big => u32::from_be_bytes(bytes),
        })
    }
    #[inline]
    fn f64(&mut self) -> Result<f64> {
        let bytes = self.take()?;
        Ok(match self.order {
            ByteOrder::Little => f64::from_le_bytes(bytes),
            ByteOrder::Big => f64::from_be_bytes(bytes),
        })
    }
    /// Reads a byte-order marker and type code. Call [`Reader::skip_srid`] after
    /// validating the type, so type errors take precedence over truncation.
    pub(crate) fn header(&mut self) -> Result<Header> {
        self.order = match self.take::<1>()?[0] {
            0 => ByteOrder::Big,
            1 => ByteOrder::Little,
            _ => return Err(WkbError::InvalidByteOrder),
        };
        let code = self.u32()?;
        let iso = (code & 0xffff) / 1000;
        let mut dims = Dims::from_ewkb(code);
        if iso == 1 || iso == 3 {
            dims |= Dims::Z;
        }
        if iso == 2 || iso == 3 {
            dims |= Dims::M;
        }
        Ok(Header {
            type_id: (code & 0xffff) % 1000,
            dims,
            has_srid: code & SRID_FLAG != 0,
        })
    }
    pub(crate) fn skip_srid(&mut self, header: Header) -> Result<()> {
        if header.has_srid {
            self.u32()?;
        }
        Ok(())
    }
    /// Reads an element count, rejecting counts of `item_bytes`-sized elements
    /// that the remaining bytes could not hold with `too_many`.
    pub(crate) fn count(&mut self, item_bytes: usize, too_many: WkbError) -> Result<usize> {
        let n = usize::try_from(self.u32()?).map_err(|_| too_many)?;
        if n > self.remaining() / item_bytes {
            return Err(too_many);
        }
        Ok(n)
    }
    #[inline]
    pub(crate) fn point(&mut self, dims: Dims) -> Result<Point> {
        Ok(Point {
            x: self.f64()?,
            y: self.f64()?,
            z: if dims.has_z() { self.f64()? } else { f64::NAN },
            m: if dims.has_m() { self.f64()? } else { f64::NAN },
        })
    }
    pub(crate) fn finish(self) -> Result<()> {
        if self.remaining() == 0 {
            Ok(())
        } else {
            Err(WkbError::TrailingBytes)
        }
    }
}

/// Bytes per point with the given dimensions.
pub(crate) fn point_bytes(dims: Dims) -> usize {
    16 + 8 * usize::from(dims.has_z()) + 8 * usize::from(dims.has_m())
}

/// Writes a little-endian geometry header.
pub(crate) fn put_header(out: &mut Vec<u8>, type_id: u32, dims: Dims) {
    out.push(1);
    put_u32(out, type_id | dims.bits());
}

/// Writes an element count.
///
/// # Panics
///
/// Panics if `len` exceeds `u32::MAX`, which WKB cannot represent.
pub(crate) fn put_len(out: &mut Vec<u8>, len: usize) {
    put_u32(
        out,
        u32::try_from(len).unwrap_or_else(|_| panic!("WKB count {len} exceeds u32")),
    );
}

#[inline]
pub(crate) fn put_point(out: &mut Vec<u8>, p: Point, dims: Dims) {
    out.extend_from_slice(&p.x.to_le_bytes());
    out.extend_from_slice(&p.y.to_le_bytes());
    if dims.has_z() {
        out.extend_from_slice(&p.z.to_le_bytes());
    }
    if dims.has_m() {
        out.extend_from_slice(&p.m.to_le_bytes());
    }
}

#[inline]
fn put_u32(out: &mut Vec<u8>, n: u32) {
    out.extend_from_slice(&n.to_le_bytes());
}
