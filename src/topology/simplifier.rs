//! `TaggedLinesSimplifier` and `ComponentJumpChecker`, GEOS 3.13.1 encounter order.
use super::{
    segment::{crossing, distance, furthest, invalid_intersection},
    shape::{Kind, Shape},
};
use crate::{
    geom::{Bounds, Dims, Point},
    index::{BoundsEntry, PackedTree},
    predicates::orient,
};
use rstar::RTree;

#[derive(Clone, Copy, Debug)]
struct Segment {
    a: Point,
    b: Point,
    output: Option<usize>,
}
impl Segment {
    fn new(a: Point, b: Point) -> Self {
        Self { a, b, output: None }
    }
    fn bounds(self) -> Bounds {
        Bounds::from_points([self.a, self.b])
    }
    fn intersects(self, other: Self) -> bool {
        invalid_intersection(self.a, self.b, other.a, other.b)
    }
}
/// A tagged line: input points, the minimum output size, the index of its first
/// segment in [`Simplifier::input`], and the simplified segments so far.
#[derive(Debug)]
struct Line {
    points: Vec<Point>,
    ring: bool,
    minimum: usize,
    offset: usize,
    result: Vec<Segment>,
}
impl Line {
    fn size(&self) -> usize {
        if self.result.is_empty() {
            0
        } else {
            self.result.len() + 1
        }
    }
    fn component_point(&self) -> Point {
        self.result.first().map_or(self.points[1], |s| s.a)
    }
}
/// An input segment: index `index` of line `line`.
#[derive(Clone, Copy, Debug)]
struct InputSegment {
    segment: Segment,
    line: usize,
    index: usize,
}
/// Segments `start..end` of line `line`, which a candidate replaces.
#[derive(Clone, Copy, Debug)]
struct Span {
    line: usize,
    start: usize,
    end: usize,
}
#[derive(Debug)]
struct Simplifier {
    lines: Vec<Line>,
    input: Vec<InputSegment>,
    alive: Vec<bool>,
    packed: PackedTree,
    output: Vec<Segment>,
    output_alive: Vec<bool>,
    tree: RTree<BoundsEntry>,
}
fn collect_lines(g: &Shape, lines: &mut Vec<Line>) {
    if g.kind.is_linear() && !g.points.is_empty() {
        lines.push(Line {
            points: g.points.clone(),
            ring: g.kind == Kind::LinearRing,
            minimum: if g.points.first() == g.points.last() {
                4
            } else {
                2
            },
            offset: 0,
            result: Vec::new(),
        });
    }
    for part in &g.parts {
        collect_lines(part, lines);
    }
}
impl Simplifier {
    fn new(g: &Shape) -> Self {
        let mut lines = Vec::new();
        collect_lines(g, &mut lines);
        let mut input = Vec::new();
        for (id, line) in lines.iter_mut().enumerate() {
            line.offset = input.len();
            for (i, p) in line.points.windows(2).enumerate() {
                input.push(InputSegment {
                    segment: Segment::new(p[0], p[1]),
                    line: id,
                    index: i,
                });
            }
        }
        let packed = PackedTree::new(input.iter().map(|s| s.segment.bounds()));
        Self {
            lines,
            alive: vec![true; input.len()],
            packed,
            input,
            output: Vec::new(),
            output_alive: Vec::new(),
            tree: RTree::new(),
        }
    }
    fn add_output(&mut self, mut s: Segment) -> Segment {
        let id = self.output.len();
        s.output = Some(id);
        self.output.push(s);
        self.output_alive.push(true);
        self.tree.insert(BoundsEntry {
            id,
            bounds: s.bounds(),
        });
        s
    }
    /// Whether `s` crosses no live output segment and no live input segment
    /// outside `exclude`.
    fn valid(&self, s: Segment, exclude: Option<Span>) -> bool {
        if self
            .tree
            .locate_in_envelope_intersecting(s.bounds().aabb())
            .any(|item| self.output_alive[item.id] && s.intersects(self.output[item.id]))
        {
            return false;
        }
        !self.packed.any(s.bounds(), |id| {
            let InputSegment {
                segment,
                line,
                index,
            } = self.input[id];
            if !self.alive[id] {
                return false;
            }
            if exclude.is_some_and(|x| x.line == line && x.start <= index && index < x.end) {
                return false;
            }
            s.intersects(segment)
        })
    }
    /// `ComponentJumpChecker::hasJump`: whether replacing `points` with `flat`
    /// would move another line's component point to the other side.
    fn jumps(&self, line: usize, points: &[Point], flat: Segment) -> bool {
        let bounds = Bounds::from_points(points.iter().copied());
        self.lines.iter().enumerate().any(|(id, other)| {
            if id == line {
                return false;
            }
            let p = other.component_point();
            bounds.contains(p)
                && points
                    .windows(2)
                    .map(|s| crossing(p, s[0], s[1]))
                    .sum::<usize>()
                    % 2
                    != crossing(p, flat.a, flat.b) % 2
        })
    }
    fn simplify_line(&mut self, id: usize, tolerance: f64) {
        // Explicit stack has GEOS's recursive left-before-right visit order.
        let mut stack = vec![(0, self.lines[id].points.len() - 1, 0)];
        while let Some((i, j, depth)) = stack.pop() {
            let depth = depth + 1;
            let line = &self.lines[id];
            if i + 1 >= j {
                if i < j {
                    let s = Segment::new(line.points[i], line.points[j]);
                    self.lines[id].result.push(s);
                }
                continue;
            }
            let (k, max) = furthest(&line.points, i, j);
            if max < 0.0 {
                let segments: Vec<_> = line.points[i..=j]
                    .windows(2)
                    .map(|p| Segment::new(p[0], p[1]))
                    .collect();
                self.lines[id].result.extend(segments);
                continue;
            }
            let flat = Segment::new(line.points[i], line.points[j]);
            let enough = line.size() >= line.minimum || depth + 1 >= line.minimum;
            if enough
                && max <= tolerance
                && self.valid(
                    flat,
                    Some(Span {
                        line: id,
                        start: i,
                        end: j,
                    }),
                )
                && !self.jumps(id, &line.points[i..=j], flat)
            {
                let offset = self.lines[id].offset;
                self.alive[offset + i..offset + j].fill(false);
                let flat = self.add_output(flat);
                self.lines[id].result.push(flat);
            } else {
                stack.push((k, j, depth));
                stack.push((i, k, depth));
            }
        }
        let line = &self.lines[id];
        if !line.ring
            || line.points.len() < 4
            || line.points.first() != line.points.last()
            || line.size() <= line.minimum
        {
            return;
        }
        let (Some(&first), Some(&last)) = (line.result.first(), line.result.last()) else {
            return;
        };
        let flat = Segment::new(last.a, first.b);
        if distance(first.a, flat.a, flat.b) > tolerance {
            return;
        }
        if orient(flat.a, flat.b, first.a) != 0
            && (!self.valid(flat, None) || self.endpoint_jumps(id, first, last, flat))
        {
            return;
        }
        for s in [first, last] {
            if let Some(index) = s.output {
                self.output_alive[index] = false;
            }
        }
        let flat = self.add_output(flat);
        self.lines[id].result.pop();
        self.lines[id].result[0] = flat;
    }
    fn endpoint_jumps(&self, line: usize, a: Segment, b: Segment, flat: Segment) -> bool {
        let bounds = Bounds::from_points([a.a, a.b, b.a, b.b]);
        self.lines.iter().enumerate().any(|(id, other)| {
            if id == line {
                return false;
            }
            let p = other.component_point();
            bounds.contains(p)
                && (crossing(p, a.a, a.b) + crossing(p, b.a, b.b)) % 2
                    != crossing(p, flat.a, flat.b) % 2
        })
    }
}
/// GEOS `TopologyPreservingSimplifier`: simplifies every line and ring in `g`
/// without introducing intersections or moving components across each other.
pub(crate) fn simplify(g: Shape, tolerance: f64) -> Shape {
    if g.is_empty() {
        return g;
    }
    let mut simplifier = Simplifier::new(&g);
    for id in 0..simplifier.lines.len() {
        simplifier.simplify_line(id, tolerance);
    }
    transform(g, &mut simplifier.lines.into_iter())
}

/// GEOS `GeometryTransformer`: rebuilds `g` from the simplified lines, in order,
/// demoting short rings and invalid polygons as GEOS does.
fn transform(mut g: Shape, lines: &mut impl Iterator<Item = Line>) -> Shape {
    match g.kind {
        Kind::LineString | Kind::LinearRing => {
            if g.points.is_empty() {
                g.dims = Dims::XY;
            } else {
                let line = lines.next().expect("one simplified line per nonempty line");
                g.points = line.result.iter().map(|s| s.a).collect();
                if let Some(s) = line.result.last() {
                    g.points.push(s.b);
                }
                g.infer_dims();
            }
            if g.kind == Kind::LinearRing && !g.points.is_empty() && g.points.len() < 4 {
                g.kind = Kind::LineString;
            }
        }
        Kind::Polygon => {
            let mut parts = g.parts.into_iter().map(|p| transform(p, lines));
            let shell = parts
                .next()
                .unwrap_or_else(|| Shape::empty(Kind::LinearRing, Dims::XY));
            let mut valid = shell.kind == Kind::LinearRing && !shell.is_empty();
            let mut components = vec![shell];
            for hole in parts {
                if hole.is_empty() {
                    continue;
                }
                if hole.kind == Kind::LinearRing {
                    components.push(hole);
                } else {
                    valid = false;
                }
            }
            if valid {
                g.parts = components;
                g.dims = Dims::XY;
            } else {
                g = Shape::build(components, Kind::GeometryCollection);
            }
        }
        Kind::MultiPoint
        | Kind::MultiLineString
        | Kind::MultiPolygon
        | Kind::GeometryCollection => {
            let parts: Vec<_> = g
                .parts
                .into_iter()
                .map(|p| transform(p, lines))
                .filter(|p| !p.is_empty())
                .collect();
            if g.kind == Kind::GeometryCollection {
                g.parts = parts;
                g.dims = Dims::XY;
            } else {
                g = Shape::build(parts, g.kind);
            }
        }
        Kind::Point => {}
    }
    g
}
