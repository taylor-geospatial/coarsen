//! `CoverageRingEdges` / `CoverageEdge`, including GEOS 3.13.1 key conventions.
use crate::geom::{Dims, Geometry, Point, Ring};
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// A maximal run of ring boundary between nodes, shared by every ring that uses it.
#[derive(Debug)]
pub(crate) struct Edge {
    pub(crate) points: Ring,
    /// A whole ring without nodes, so its closing vertex may also move.
    pub(crate) free_ring: bool,
    /// Rings using this edge: one on the coverage boundary, two inside.
    pub(crate) ring_count: usize,
}

/// Unique edges plus, per input ring in encounter order, the edges forming it.
#[derive(Debug)]
pub(crate) struct Coverage {
    pub(crate) edges: Vec<Edge>,
    rings: Vec<Vec<usize>>,
}

/// GEOS `CoverageEdge::key`: an endpoint vertex id and the id of its neighbor
/// along the edge, which identifies an edge regardless of direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct EdgeKey(u32, u32);

/// One ring's section between consecutive nodes, before deduplication.
#[derive(Debug)]
struct RingEdge {
    points: Ring,
    free_ring: bool,
    key: EdgeKey,
}

/// A ring segment by its sorted vertex ids, and where it occurs.
#[derive(Clone, Copy, Debug)]
struct Segment {
    a: u32,
    b: u32,
    ring: u32,
    pos: u32,
}

fn to_id(i: usize) -> u32 {
    u32::try_from(i).unwrap_or_else(|_| panic!("id {i} exceeds u32"))
}
fn index(id: u32) -> usize {
    usize::try_from(id).expect("u32 fits in usize")
}

impl Coverage {
    /// Splits every ring at its nodes and deduplicates the resulting edges.
    ///
    /// # Panics
    ///
    /// Panics if the coverage has `u32::MAX` or more distinct vertices.
    pub(crate) fn extract(geoms: &[Geometry]) -> Self {
        let rings: Vec<&Ring> = geoms.iter().flat_map(Geometry::rings).collect();
        let (ring_ids, vertex_count) = vertex_ids(&rings);
        let nodes = find_nodes(&ring_ids, vertex_count);
        let split: Vec<Vec<RingEdge>> = ring_ids
            .into_par_iter()
            // Indexing beats `.zip(&rings)`, measured ~25% faster at 8 threads.
            .enumerate()
            .map(|(ri, ids)| split_ring(ids, rings[ri], &nodes))
            .collect();
        let mut edges = Vec::new();
        let mut ring_edges = Vec::with_capacity(split.len());
        let mut unique = FxHashMap::default();
        for ring in split {
            let mut refs = Vec::with_capacity(ring.len());
            for edge in ring {
                let id = *unique.entry(edge.key).or_insert_with(|| {
                    edges.push(Edge {
                        points: edge.points,
                        free_ring: edge.free_ring,
                        ring_count: 0,
                    });
                    edges.len() - 1
                });
                edges[id].ring_count += 1;
                refs.push(id);
            }
            ring_edges.push(refs);
        }
        Self {
            edges,
            rings: ring_edges,
        }
    }

    /// Replaces each ring of `geoms`, which must be the geometries passed to
    /// [`Coverage::extract`], with its edges from `simplified`, indexed like `edges`.
    pub(crate) fn rebuild(&self, geoms: &mut [Geometry], simplified: &[Ring]) {
        let mut rings = self.rings.iter();
        for polygon in geoms.iter_mut().flat_map(Geometry::polygons_mut) {
            let mut dims = Dims::XY;
            for ring in &mut polygon.rings {
                let refs = rings.next().expect("one edge list per extracted ring");
                if refs.is_empty() {
                    dims |= polygon.dims;
                    continue;
                }
                let pts = assemble(refs, simplified);
                // GEOS's dimension-unknown CoordinateSequence inspects only
                // the first coordinate when inferring Z, even for mixed rings.
                if pts.first().is_some_and(|p| !p.z.is_nan()) {
                    dims |= Dims::Z;
                }
                *ring = pts;
            }
            if !polygon.rings.is_empty() {
                polygon.dims = dims;
            }
        }
    }
}

/// Maps every ring vertex to a dense id. A parallel sort gives deterministic ids
/// in XY order. Returns the per-ring ids and the number of distinct vertices.
fn vertex_ids(rings: &[&Ring]) -> (Vec<Vec<u32>>, usize) {
    let mut vertices: Vec<Point> = rings
        .par_iter()
        .flat_map_iter(|r| r.iter().copied())
        .collect();
    vertices.par_sort_unstable_by(|a, b| a.cmp_xy(*b));
    vertices.dedup();
    assert!(
        vertices.len() < index(u32::MAX),
        "{} distinct vertices exceed u32 ids",
        vertices.len()
    );
    let ids: FxHashMap<_, _> = vertices
        .iter()
        .enumerate()
        .map(|(i, p)| (p.key(), to_id(i)))
        .collect();
    let ring_ids = rings
        .par_iter()
        .map(|r| r.iter().map(|p| ids[&p.key()]).collect())
        .collect();
    (ring_ids, vertices.len())
}

/// GEOS `findNodes`: vertices used by more than two rings, with boundary degree
/// above three, or where a ring switches between boundary and shared segments.
fn find_nodes(ring_ids: &[Vec<u32>], vertex_count: usize) -> Vec<bool> {
    let mut counts = vec![0u32; vertex_count];
    for r in ring_ids {
        for &v in r.iter().skip(usize::from(r.len() >= 4)) {
            counts[index(v)] += 1;
        }
    }
    let (boundary, degree) = boundary_segments(ring_ids, vertex_count);
    // findBoundaryNodes initializes each counter at zero, not one, in 3.13.1.
    let mut nodes: Vec<bool> = counts
        .iter()
        .zip(degree)
        .map(|(&n, d)| n > 2 || d > 3)
        .collect();
    drop(counts);
    for (r, b) in ring_ids.iter().zip(&boundary) {
        let Some(&last) = b.last() else {
            continue;
        };
        let mut prev = last;
        for (i, &v) in b.iter().enumerate() {
            if prev != v {
                nodes[index(r[i])] = true;
            }
            prev = v;
        }
    }
    nodes
}

/// Flags each ring segment used an odd number of times (a coverage boundary
/// segment) and counts boundary segments per vertex.
fn boundary_segments(ring_ids: &[Vec<u32>], vertex_count: usize) -> (Vec<Vec<bool>>, Vec<u32>) {
    let mut segments: Vec<Segment> = ring_ids
        .par_iter()
        .enumerate()
        .flat_map_iter(|(ri, r)| {
            let ring = to_id(ri);
            r.windows(2).enumerate().map(move |(i, w)| Segment {
                a: w[0].min(w[1]),
                b: w[0].max(w[1]),
                ring,
                pos: to_id(i),
            })
        })
        .collect();
    segments.par_sort_unstable_by_key(|s| (s.a, s.b));
    let mut boundary: Vec<Vec<bool>> = ring_ids
        .iter()
        .map(|r| vec![false; r.len().saturating_sub(1)])
        .collect();
    let mut degree = vec![0u32; vertex_count];
    for group in segments.chunk_by(|s, t| (s.a, s.b) == (t.a, t.b)) {
        if group.len() % 2 == 1 {
            degree[index(group[0].a)] += 1;
            degree[index(group[0].b)] += 1;
            for s in group {
                boundary[index(s.ring)][index(s.pos)] = true;
            }
        }
    }
    (boundary, degree)
}

/// Splits one ring, given as vertex ids `r` and coordinates, at its nodes. A ring
/// without nodes is one free edge.
fn split_ring(mut r: Vec<u32>, coords: &Ring, nodes: &[bool]) -> Vec<RingEdge> {
    r.dedup();
    if r.len() < 3 {
        return Vec::new();
    }
    let mut coords = coords.clone();
    coords.dedup();
    // GEOS CoverageEdge copies Coordinate (XYZ), dropping M.
    for p in &mut coords {
        p.m = f64::NAN;
    }
    let n = r.len() - 1;
    let cuts: Vec<usize> = (0..n).filter(|&i| nodes[index(r[i])]).collect();
    let ring_key = || {
        // GEOS comment says lowest; implementation selects highest.
        let high = (0..n).max_by_key(|&i| r[i]).expect("ring has vertices");
        EdgeKey(r[high], r[(high + 1) % n].min(r[(high + n - 1) % n]))
    };
    if cuts.is_empty() {
        return vec![RingEdge {
            points: coords,
            free_ring: true,
            key: ring_key(),
        }];
    }
    cuts.iter()
        .enumerate()
        .map(|(ci, &a)| {
            let b = cuts[(ci + 1) % cuts.len()];
            let key = if a == b {
                ring_key()
            } else if r[a] < r[b] {
                EdgeKey(r[a], r[(a + 1) % n])
            } else {
                EdgeKey(r[b], r[(b + n - 1) % n])
            };
            let len = if b > a { b - a + 1 } else { n - a + b + 1 };
            RingEdge {
                points: (0..len).map(|i| coords[(a + i) % n]).collect(),
                free_ring: false,
                key,
            }
        })
        .collect()
}

/// Concatenates a ring's simplified edges, orienting each to continue the last.
fn assemble(refs: &[usize], simplified: &[Ring]) -> Ring {
    let mut pts: Ring = Vec::new();
    for (i, &e) in refs.iter().enumerate() {
        let edge = &simplified[e];
        let forward = if refs.len() <= 2 && i == 0 {
            true
        } else if i == 0 {
            let next = &simplified[refs[1]];
            edge.last() == next.first() || edge.last() == next.last()
        } else {
            pts.last() == edge.first()
        };
        let mut append = |p: Point| {
            if pts.last() != Some(&p) {
                pts.push(p);
            }
        };
        if forward {
            for &p in edge {
                append(p);
            }
        } else {
            for &p in edge.iter().rev() {
                append(p);
            }
        }
    }
    pts
}
