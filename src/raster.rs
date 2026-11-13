//! Pair complementary one-pixel corners by their exact diagonal endpoints.
use crate::geom::{Dims, Geometry, Point};
use rustc_hash::FxHashMap;

type XY = (u64, u64);
type Diagonal = (XY, XY);

#[derive(Clone, Copy)]
struct Corner {
    shell: usize,
    index: usize,
    tip: Point,
}

// Keep only diagonals traversed by exactly two different shells.
#[derive(Default)]
enum Hits {
    #[default]
    Empty,
    One(Corner),
    Two(Corner, Corner),
    Ambiguous,
}

fn diagonal(a: Point, b: Point) -> Diagonal {
    let (u, v) = (a.key(), b.key());
    if u < v { (u, v) } else { (v, u) }
}

fn midpoint(a: Point, b: Point) -> Point {
    Point {
        x: a.x.midpoint(b.x),
        y: a.y.midpoint(b.y),
        ..Point::default()
    }
}

/// Round repeated, complementary exterior staircases on an axis-aligned XY raster.
/// Supply valid, non-overlapping cell polygons on one grid, in grid coordinate units.
/// Holes are never edited. Coordinates at multi-parcel junctions are protected.
pub(crate) fn smooth(geoms: &mut [Geometry], pixel: f64) {
    let eps = pixel * 1e-5;
    let mut hits = FxHashMap::<Diagonal, Hits>::default();
    let mut owners = FxHashMap::<XY, (usize, usize)>::default();
    let mut pairs = Vec::new();
    let mut marks = Vec::new();
    for (shell, poly) in geoms.iter().flat_map(Geometry::polygons).enumerate() {
        let n = poly.rings.first().map_or(0, |r| r.len().saturating_sub(1));
        pairs.push(vec![None; n]);
        marks.push(vec![false; n]);
        for (k, ring) in poly.rings.iter().enumerate() {
            for p in ring.iter().skip(1) {
                let entry = owners.entry(p.key()).or_insert((shell, shell));
                if k != 0 || (entry.0 != shell && entry.1 != entry.0 && entry.1 != shell) {
                    *entry = (usize::MAX, usize::MAX);
                } else if entry.0 != shell && entry.0 != usize::MAX {
                    entry.1 = shell;
                }
            }
        }
        if n < 5 || poly.dims != Dims::XY {
            continue;
        }
        let r = &poly.rings[0];
        for i in 0..n {
            let (a, tip, b) = (r[(i + n - 1) % n], r[i], r[(i + 1) % n]);
            let (dx, dy) = ((tip.x - a.x).abs(), (tip.y - a.y).abs());
            let (ex, ey) = ((b.x - tip.x).abs(), (b.y - tip.y).abs());
            if !((dx <= eps && (dy - pixel).abs() <= eps && ey <= eps && (ex - pixel).abs() <= eps)
                || (dy <= eps
                    && (dx - pixel).abs() <= eps
                    && ex <= eps
                    && (ey - pixel).abs() <= eps))
            {
                continue;
            }
            let hit = Corner {
                shell,
                index: i,
                tip,
            };
            let count = hits.entry(diagonal(a, b)).or_default();
            *count = match *count {
                Hits::Empty => Hits::One(hit),
                Hits::One(prev)
                    if prev.shell != shell
                        && (prev.tip.x - tip.x).abs() > pixel / 2.
                        && (prev.tip.y - tip.y).abs() > pixel / 2. =>
                {
                    Hits::Two(prev, hit)
                }
                _ => Hits::Ambiguous,
            };
        }
    }
    for ((a, b), hit) in hits {
        if let Hits::Two(u, v) = hit {
            let endpoints_clear = [a, b].iter().all(|xy| owners[xy].0 != usize::MAX);
            let tips_clear = [u, v]
                .iter()
                .all(|k| owners[&k.tip.key()] == (k.shell, k.shell));
            if endpoints_clear && tips_clear {
                pairs[u.shell][u.index] = Some(v.shell);
                pairs[v.shell][v.index] = Some(u.shell);
            }
        }
    }
    for (shell, poly) in geoms
        .iter_mut()
        .flat_map(Geometry::polygons_mut)
        .enumerate()
    {
        let Some(r) = poly.rings.first_mut() else {
            continue;
        };
        let n = marks[shell].len();
        for i in 0..n {
            let j = (i + 2) % n;
            if let Some(partner) = pairs[shell][i] {
                if pairs[shell][j] != Some(partner) {
                    continue;
                }
                let (a, m, b) = (r[(i + n - 1) % n], r[(i + 1) % n], r[(j + 1) % n]);
                if (b.x - m.x - (m.x - a.x)).abs() <= eps && (b.y - m.y - (m.y - a.y)).abs() <= eps
                {
                    for k in [i, (i + 1) % n, j] {
                        marks[shell][k] = true;
                    }
                }
            }
        }
        if !marks[shell].contains(&true) {
            continue;
        }
        round_marked(r, &marks[shell]);
    }
}

fn round_marked(r: &mut Vec<Point>, marked: &[bool]) {
    let n = marked.len();
    let mut out = Vec::with_capacity(n + 1);
    for i in 0..n {
        if marked[i] {
            for p in [
                midpoint(r[(i + n - 1) % n], r[i]),
                midpoint(r[i], r[(i + 1) % n]),
            ] {
                if out.last() != Some(&p) {
                    out.push(p);
                }
            }
        } else {
            out.push(r[i]);
        }
    }
    if out.last() != out.first() {
        out.push(out[0]);
    }
    *r = out;
}
