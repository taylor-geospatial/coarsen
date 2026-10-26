//! Ordered edge-dependency rounds. Intersecting envelopes retain GEOS edge order.
//! A removed corner and all future corners remain inside the original envelope.
//! Therefore edges in the same round cannot observe each other's mutations.
use crate::{
    Boundary,
    edges::Edge,
    geom::{Bounds, Ring},
    index::{BoundsEntry, PackedTree},
    predicates::{triangle_area, triangle_contains},
};
use rayon::prelude::*;
use rstar::RTree;
use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed},
};

/// Edges with more vertices than this get a packed vertex index.
const INDEX_THRESHOLD: usize = 16;
/// Linked-list neighbor of a removed vertex or of an open edge's endpoint.
const DETACHED: usize = usize::MAX;

/// An edge plus the vertex state other edges observe while it simplifies.
#[derive(Debug)]
struct Indexed<'a> {
    edge: &'a Edge,
    index: Option<PackedTree>,
    alive: Vec<AtomicBool>,
    size: AtomicUsize,
}
/// A removable vertex `i` with the neighbors its area was computed from.
#[derive(Clone, Copy, Debug)]
struct Corner {
    area: f64,
    i: usize,
    prev: usize,
    next: usize,
}
impl PartialEq for Corner {
    fn eq(&self, b: &Self) -> bool {
        self.area == b.area && self.i == b.i
    }
}
impl Eq for Corner {}
impl PartialOrd for Corner {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        Some(self.cmp(b))
    }
}
/// Reversed so [`BinaryHeap`] pops the smallest area, then the lowest index.
impl Ord for Corner {
    fn cmp(&self, b: &Self) -> Ordering {
        b.area
            .partial_cmp(&self.area)
            .expect("only areas <= tolerance, never NaN, are queued")
            .then(b.i.cmp(&self.i))
    }
}
fn is_ring(e: &Edge) -> bool {
    e.points.len() >= 4 && e.points.first() == e.points.last()
}
/// Removes consecutive duplicates and, if `e` is a ring, restores the closing point.
fn normalize(mut points: Ring, e: &Edge) -> Ring {
    points.dedup();
    if is_ring(e) && !points.is_empty() {
        points.push(points[0]);
    }
    points
}

/// An edge without an initial removable-size corner cannot ever change: no
/// vertex removal can create a new corner until its own first removal occurs.
fn can_change(e: &Edge, area_tolerance: f64) -> bool {
    let ring = is_ring(e);
    let n = e.points.len() - usize::from(ring);
    if n <= if ring { 3 } else { 2 } {
        return false;
    }
    let start = usize::from(!(ring && e.free_ring));
    (start..e.points.len() - 1).any(|i| {
        triangle_area(
            e.points[(i + n - 1) % n],
            e.points[i],
            e.points[(i + 1) % n],
        ) <= area_tolerance
    })
}

/// GEOS `TPVWSimplifier` over coverage edges, returning one ring per edge.
///
/// With [`Boundary::Preserve`], only edges shared by exactly two rings change.
pub(crate) fn simplify(edges: &[Edge], tolerance: f64, boundary: Boundary) -> Vec<Ring> {
    let all = boundary == Boundary::Simplify;
    let area_tolerance = tolerance * tolerance;
    let active: Vec<bool> = edges
        .par_iter()
        .map(|e| (all || e.ring_count == 2) && can_change(e, area_tolerance))
        .collect();
    // No removals means no spatial constraints or dependency rounds are needed.
    // Retain the same duplicate and ring-closure handling as the indexed path.
    if !active.contains(&true) {
        return edges
            .par_iter()
            .map(|e| {
                let n = e.points.len() - usize::from(is_ring(e));
                normalize(e.points[..n].to_vec(), e)
            })
            .collect();
    }
    let index: Vec<_> = edges
        .par_iter()
        .map(|e| {
            let n = e.points.len() - usize::from(is_ring(e));
            Indexed {
                edge: e,
                index: (n > INDEX_THRESHOLD).then(|| {
                    PackedTree::new(e.points[..n].iter().map(|&p| Bounds::from_points([p])))
                }),
                alive: (0..n).map(|_| AtomicBool::new(true)).collect(),
                size: AtomicUsize::new(n),
            }
        })
        .collect();
    let neighbors = neighbors(edges, &active, boundary);
    for round in dependency_rounds(&active, &neighbors) {
        round
            .par_iter()
            .for_each(|&id| simplify_edge(id, &index, &neighbors[id], area_tolerance));
    }
    index
        .par_iter()
        .map(|e| {
            let pts = e
                .edge
                .points
                .iter()
                .zip(&e.alive)
                .filter(|(_, alive)| alive.load(Relaxed))
                .map(|(&p, _)| p)
                .collect();
            normalize(pts, e.edge)
        })
        .collect()
}

/// Edges whose envelopes intersect each active edge's envelope. Immutable edges
/// still constrain others, but need no list of their own.
fn neighbors(edges: &[Edge], active: &[bool], boundary: Boundary) -> Vec<Vec<usize>> {
    let tree = RTree::bulk_load(
        edges
            .iter()
            .enumerate()
            .map(|(id, e)| BoundsEntry {
                id,
                bounds: Bounds::from_points(e.points.iter().copied()),
            })
            .collect(),
    );
    edges
        .par_iter()
        .enumerate()
        .map(|(i, e)| {
            if !active[i] {
                return Vec::new();
            }
            tree.locate_in_envelope_intersecting(
                Bounds::from_points(e.points.iter().copied()).aabb(),
            )
            // simplifyInner excludes edges used by >2 rings from constraints.
            .filter(|v| boundary == Boundary::Simplify || edges[v.id].ring_count <= 2)
            .map(|v| v.id)
            .collect()
        })
        .collect()
}

/// Groups active edges so each runs one round after every earlier intersecting
/// edge, preserving GEOS's sequential edge order where it matters.
fn dependency_rounds(active: &[bool], neighbors: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut levels = vec![0; active.len()];
    let mut rounds: Vec<Vec<usize>> = vec![];
    for i in 0..active.len() {
        if !active[i] {
            continue;
        }
        levels[i] = neighbors[i]
            .iter()
            .filter(|&&j| j < i && active[j])
            .map(|&j| levels[j] + 1)
            .max()
            .unwrap_or(0);
        if rounds.len() <= levels[i] {
            rounds.resize_with(levels[i] + 1, Vec::new);
        }
        rounds[levels[i]].push(i);
    }
    rounds
}
fn simplify_edge(id: usize, edges: &[Indexed<'_>], neighbors: &[usize], area_tolerance: f64) {
    let this = &edges[id];
    let e = this.edge;
    let pts = &e.points;
    let ring = is_ring(e);
    let n = this.alive.len();
    let min = if ring { 3 } else { 2 };
    if n <= min {
        return;
    }
    let mut prev: Vec<usize> = (0..n)
        .map(|i| {
            if i == 0 {
                if ring { n - 1 } else { DETACHED }
            } else {
                i - 1
            }
        })
        .collect();
    let mut next: Vec<usize> = (0..n)
        .map(|i| {
            if i + 1 == n {
                if ring { 0 } else { DETACHED }
            } else {
                i + 1
            }
        })
        .collect();
    let mut q = BinaryHeap::new();
    let add = |i: usize, prev: &[usize], next: &[usize], q: &mut BinaryHeap<Corner>| {
        if e.free_ring || (i != 0 && i != pts.len() - 1) {
            let a = triangle_area(pts[prev[i]], pts[i], pts[next[i]]);
            if a <= area_tolerance {
                q.push(Corner {
                    area: a,
                    i,
                    prev: prev[i],
                    next: next[i],
                });
            }
        }
    };
    for i in usize::from(!(ring && e.free_ring))..pts.len() - 1 {
        add(i, &prev, &next, &mut q);
    }
    let mut size = n;
    while size > min {
        let Some(c) = q.pop() else {
            break;
        };
        if prev[c.i] != c.prev || next[c.i] != c.next {
            continue;
        }
        if c.area > area_tolerance {
            break;
        }
        let a = pts[c.prev];
        let b = pts[c.i];
        let d = pts[c.next];
        let env = Bounds::from_points([a, b, d]);
        let blocked = neighbors.iter().any(|&j| {
            let other = &edges[j];
            let intersects = |v: usize| {
                if !other.alive[v].load(Relaxed) {
                    return false;
                }
                let p = other.edge.points[v];
                env.contains(p) && p != a && p != b && p != d && triangle_contains(a, b, d, p)
            };
            let hit = match &other.index {
                Some(index) => index.any(env, intersects),
                None => (0..other.alive.len()).any(intersects),
            };
            if hit {
                return true;
            }
            if j != id && other.size.load(Relaxed) == 2 {
                let p = other.edge.points[0];
                let q = *other.edge.points.last().expect("edges have points");
                return (a == p && d == q) || (a == q && d == p);
            }
            false
        });
        if blocked {
            continue;
        }
        next[c.prev] = c.next;
        prev[c.next] = c.prev;
        prev[c.i] = DETACHED;
        next[c.i] = DETACHED;
        this.alive[c.i].store(false, Relaxed);
        size -= 1;
        this.size.store(size, Relaxed);
        add(c.prev, &prev, &next, &mut q);
        add(c.next, &prev, &next, &mut q);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::Point;

    fn pt(x: f64, y: f64) -> Point {
        Point {
            x,
            y,
            ..Point::default()
        }
    }

    #[test]
    fn parallel_order_and_all_boundaries() {
        let edges: Vec<_> = (0..40u32)
            .map(|i| {
                let x = f64::from(i % 10) * 10.;
                let y = f64::from(i / 10) * 10.;
                Edge {
                    points: vec![
                        pt(x, y),
                        pt(x + 2., y + 0.1),
                        pt(x + 4., y),
                        pt(x + 4., y + 4.),
                        pt(x, y + 4.),
                        pt(x, y),
                    ],
                    free_ring: true,
                    ring_count: 1,
                }
            })
            .collect();
        let run = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| simplify(&edges, 1., Boundary::Simplify))
        };
        let a = run(1);
        assert_eq!(a, run(8));
        assert!(a.iter().all(|r| r.len() == 5));
        let preserved = simplify(&edges, 1., Boundary::Preserve);
        assert!(
            preserved.iter().all(|r| r.len() == 6),
            "boundary-only edges must not change"
        );
    }
}
