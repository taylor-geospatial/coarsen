//! GEOS `CoveragePolygonValidator` at `gapWidth=0`: matching, interactions, interiors.
use crate::{
    geom::{Bounds, Geometry, Point, Ring},
    index::{BoundsEntry, PackedTree},
    predicates::orient,
};
use rayon::prelude::*;
use rstar::RTree;
use rustc_hash::FxHashMap;
use std::mem;

type Key = ((u64, u64), (u64, u64));
fn key(a: Point, b: Point) -> Key {
    if a.cmp_xy(b).is_lt() {
        (a.key(), b.key())
    } else {
        (b.key(), a.key())
    }
}
/// A deduplicated ring with a segment index; `right` when its interior lies to
/// the right of its direction.
#[derive(Debug)]
struct IndexedRing {
    pts: Ring,
    right: bool,
    index: PackedTree,
    bounds: Bounds,
}
#[derive(Debug)]
struct IndexedPolygon {
    rings: Vec<IndexedRing>,
    bounds: Bounds,
}
#[derive(Debug)]
struct IndexedGeom {
    polys: Vec<IndexedPolygon>,
    bounds: Bounds,
}
impl IndexedRing {
    /// Indexes the ring at `position` in its polygon; position 0 is the shell.
    fn new(ring: &Ring, position: usize) -> Self {
        let mut pts = ring.clone();
        pts.dedup();
        let index = PackedTree::new(
            pts.windows(2)
                .map(|w| Bounds::from_points(w.iter().copied())),
        );
        let bounds = Bounds::from_points(pts.iter().copied());
        Self {
            right: is_ccw(&pts) != (position == 0),
            pts,
            index,
            bounds,
        }
    }
}
impl IndexedGeom {
    fn new(g: &Geometry) -> Self {
        let polys = g
            .polygons()
            .iter()
            .map(|p| {
                let rings: Vec<_> = p
                    .rings()
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| !r.is_empty())
                    .map(|(i, r)| IndexedRing::new(r, i))
                    .collect();
                let bounds = rings.first().map_or(Bounds::default(), |r| r.bounds);
                IndexedPolygon { rings, bounds }
            })
            .collect();
        Self {
            polys,
            bounds: g.bounds(),
        }
    }
}

/// Where a point lies relative to a ring.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Location {
    Exterior,
    Boundary,
    Interior,
}

/// GEOS `Orientation::isCCW` cap algorithm (also defined for degenerate rings).
fn is_ccw(r: &[Point]) -> bool {
    if r.len() < 4 {
        return false;
    }
    let n = r.len() - 1;
    let mut hi = 0;
    for i in 1..=n {
        if r[i].y > r[i - 1].y && r[i].y >= r[hi].y {
            hi = i;
        }
    }
    if hi == 0 {
        return false;
    }
    let mut low = (hi + 1) % n;
    while r[low].y == r[hi].y {
        low = (low + 1) % n;
    }
    let down = (low + n - 1) % n;
    if r[hi] == r[down] {
        orient(r[hi - 1], r[hi], r[low]) > 0
    } else {
        r[down].x - r[hi].x < 0.0
    }
}
fn greater(o: Point, p: Point, q: Point) -> bool {
    let quadrant = |v: Point| {
        let x = v.x - o.x;
        let y = v.y - o.y;
        if x >= 0. {
            if y >= 0. { 0 } else { 3 }
        } else if y >= 0. {
            1
        } else {
            2
        }
    };
    let a = quadrant(p);
    let b = quadrant(q);
    a > b || (a == b && orient(o, q, p) > 0)
}
fn interior(o: Point, mut a: Point, mut b: Point, p: Point) -> bool {
    let between = !greater(o, a, b);
    if !between {
        mem::swap(&mut a, &mut b);
    }
    (greater(o, p, a) && !greater(o, p, b)) == between
}
fn on(a: Point, b: Point, p: Point) -> bool {
    Bounds::from_points([a, b]).contains(p) && orient(a, b, p) == 0
}
fn interacting(a: Point, b: Point, r: &IndexedRing, j: usize) -> bool {
    let c = r.pts[j];
    let d = r.pts[j + 1];
    if a == b || c == d || key(a, b) == key(c, d) {
        return false;
    }
    if !Bounds::from_points([a, b]).intersects(Bounds::from_points([c, d])) {
        return false;
    }
    let ac = orient(a, b, c);
    let ad = orient(a, b, d);
    let ca = orient(c, d, a);
    let cb = orient(c, d, b);
    if ac * ad < 0 && ca * cb < 0 {
        return true;
    }
    for (p, x, y) in [(a, c, d), (b, c, d), (c, a, b), (d, a, b)] {
        if p != x && p != y && on(x, y, p) {
            return true;
        }
    }
    let shared = if a == c || a == d {
        Some((a, b))
    } else if b == c || b == d {
        Some((b, a))
    } else {
        None
    };
    if let Some((o, p)) = shared {
        let n = r.pts.len() - 1;
        let at = if o == c { j } else { (j + 1) % n };
        let mut prev = r.pts[(at + n - 1) % n];
        let mut next = r.pts[(at + 1) % n];
        if p == prev || p == next {
            return false;
        }
        if !r.right {
            mem::swap(&mut prev, &mut next);
        }
        return interior(o, prev, next, p);
    }
    false
}
/// Indexed GEOS-style ray crossing.
fn locate(p: Point, r: &IndexedRing) -> Location {
    if !r.bounds.contains(p) {
        return Location::Exterior;
    }
    let ray = Bounds {
        min: [p.x, p.y],
        max: [r.bounds.max[0], p.y],
    };
    let mut crossings = 0;
    let boundary = r.index.any(ray, |i| {
        let a = r.pts[i];
        let b = r.pts[i + 1];
        if on(a, b, p) {
            return true;
        }
        if (a.y > p.y && b.y <= p.y) || (b.y > p.y && a.y <= p.y) {
            let mut s = orient(a, b, p);
            if b.y < a.y {
                s = -s;
            }
            if s > 0 {
                crossings += 1;
            }
        }
        false
    });
    if boundary {
        Location::Boundary
    } else if crossings % 2 == 1 {
        Location::Interior
    } else {
        Location::Exterior
    }
}
fn contains(p: Point, poly: &IndexedPolygon) -> bool {
    let Some((shell, holes)) = poly.rings.split_first() else {
        return false;
    };
    poly.bounds.contains(p)
        && locate(p, shell) == Location::Interior
        && holes.iter().all(|r| locate(p, r) == Location::Exterior)
}

/// Flags each polygonal geometry that is invalid within the coverage: one with a
/// duplicated or same-direction shared segment, or whose unmatched segments
/// cross, touch the interior of, or lie inside a neighbor.
pub fn invalid_mask(geoms: &[Geometry]) -> Vec<bool> {
    let indexed: Vec<_> = geoms.par_iter().map(IndexedGeom::new).collect();
    let tree = RTree::bulk_load(
        indexed
            .iter()
            .enumerate()
            .filter(|(_, g)| !g.bounds.is_empty())
            .map(|(id, g)| BoundsEntry {
                id,
                bounds: g.bounds,
            })
            .collect(),
    );
    (0..geoms.len())
        .into_par_iter()
        .map(|id| {
            let g = &indexed[id];
            if g.bounds.is_empty() {
                return false;
            }
            let near: Vec<_> = tree
                .locate_in_envelope_intersecting(g.bounds.aabb())
                .filter(|v| v.id != id)
                .map(|v| &indexed[v.id])
                .collect();
            is_invalid(g, &near)
        })
        .collect()
}

/// A segment's direction relative to its ring interior, and its position in the
/// target list.
type SegmentSide = (bool, usize);

fn is_invalid(g: &IndexedGeom, near: &[&IndexedGeom]) -> bool {
    let rings: Vec<_> = g.polys.iter().flat_map(|p| &p.rings).collect();
    let mut targets = Vec::new();
    let mut matches = Vec::new();
    let mut map: FxHashMap<Key, SegmentSide> = FxHashMap::default();
    for (ri, r) in rings.iter().enumerate() {
        for i in 0..r.pts.len().saturating_sub(1) {
            let a = r.pts[i];
            let b = r.pts[i + 1];
            let dir = a.cmp_xy(b).is_lt() == r.right;
            let k = key(a, b);
            let pos = targets.len();
            targets.push((ri, i));
            matches.push(false);
            if let Some(&(other, first)) = map.get(&k) {
                if dir == other {
                    return true;
                }
                matches[pos] = true;
                matches[first] = true;
            } else {
                map.insert(k, (dir, pos));
            }
        }
    }
    for adj in near {
        for p in &adj.polys {
            for r in &p.rings {
                let bad = r.index.any(g.bounds, |i| {
                    let a = r.pts[i];
                    let b = r.pts[i + 1];
                    if let Some(&(dir, pos)) = map.get(&key(a, b)) {
                        if (a.cmp_xy(b).is_lt() == r.right) == dir {
                            return true;
                        }
                        matches[pos] = true;
                    }
                    false
                });
                if bad {
                    return true;
                }
            }
        }
    }
    targets
        .iter()
        .zip(&matches)
        .filter(|(_, matched)| !**matched)
        .any(|(&(ri, i), _)| {
            let r = rings[ri];
            unmatched_is_invalid(r.pts[i], r.pts[i + 1], near)
        })
}

/// Whether an unmatched segment `a`→`b` interacts with or lies inside a neighbor.
fn unmatched_is_invalid(a: Point, b: Point, near: &[&IndexedGeom]) -> bool {
    let env = Bounds::from_points([a, b]);
    near.iter()
        .flat_map(|adj| &adj.polys)
        .filter(|poly| poly.bounds.intersects(env))
        .any(|poly| {
            poly.rings
                .iter()
                .any(|ring| ring.index.any(env, |j| interacting(a, b, ring, j)))
                || contains(a, poly)
        })
}
