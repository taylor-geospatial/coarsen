//! Bounded, endian-aware WKB codec for polygonal coverages: [`read`] and [`write()`].
pub(crate) mod codec;
mod error;
mod polygonal;

pub use error::WkbError;
pub use polygonal::{read, write};
