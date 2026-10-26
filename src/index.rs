//! Spatial indexes: an R-tree entry type and a packed, order-preserving envelope tree.
use crate::geom::Bounds;
use rstar::{AABB, RTreeObject};

/// Entries per leaf and per internal node of [`PackedTree`]. Larger nodes mean
/// fewer envelopes to store but more candidates tested per intersecting node.
const NODE_CAPACITY: usize = 16;

impl Bounds {
    pub(crate) fn aabb(self) -> AABB<[f64; 2]> {
        AABB::from_corners(self.min, self.max)
    }
}

/// An R-tree entry: a caller-defined id and its envelope.
#[derive(Clone, Debug)]
pub(crate) struct BoundsEntry {
    pub(crate) id: usize,
    pub(crate) bounds: Bounds,
}
impl RTreeObject for BoundsEntry {
    type Envelope = AABB<[f64; 2]>;
    fn envelope(&self) -> Self::Envelope {
        self.bounds.aabb()
    }
}

/// Packed envelope tree over a sequence, retaining sequence order.
///
/// Stores only one envelope per node, not a heap allocation or R-tree entry per item.
#[derive(Debug)]
pub(crate) struct PackedTree {
    levels: Vec<Vec<Bounds>>,
    len: usize,
}
impl PackedTree {
    pub(crate) fn new(bounds: impl IntoIterator<Item = Bounds>) -> Self {
        let mut leaves = Vec::new();
        let mut b = Bounds::default();
        let mut n = 0;
        for v in bounds {
            b.union(v);
            n += 1;
            if n % NODE_CAPACITY == 0 {
                leaves.push(b);
                b = Bounds::default();
            }
        }
        if n % NODE_CAPACITY != 0 {
            leaves.push(b);
        }
        let mut levels = vec![leaves];
        while let Some(top) = levels.last().filter(|top| top.len() > 1) {
            let upper = top
                .chunks(NODE_CAPACITY)
                .map(|c| {
                    let mut b = Bounds::default();
                    for &v in c {
                        b.union(v);
                    }
                    b
                })
                .collect();
            levels.push(upper);
        }
        Self { levels, len: n }
    }
    /// Whether `test` holds for any item whose envelope intersects `b`, visiting
    /// candidates in sequence order and stopping at the first match.
    pub(crate) fn any(&self, b: Bounds, mut test: impl FnMut(usize) -> bool) -> bool {
        if self.len == 0 {
            return false;
        }
        self.visit(self.levels.len() - 1, 0, b, &mut test)
    }
    fn visit(
        &self,
        level: usize,
        node: usize,
        b: Bounds,
        test: &mut impl FnMut(usize) -> bool,
    ) -> bool {
        if !self.levels[level][node].intersects(b) {
            return false;
        }
        let children = if level == 0 {
            self.len
        } else {
            self.levels[level - 1].len()
        };
        let mut range = node * NODE_CAPACITY..((node + 1) * NODE_CAPACITY).min(children);
        if level == 0 {
            return range.any(test);
        }
        range.any(|j| self.visit(level - 1, j, b, test))
    }
}
