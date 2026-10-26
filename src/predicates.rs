//! GEOS 3.13.1 `CGAlgorithmsDD` orientation with DD arithmetic, operation order
//! preserved, and triangle predicates built on it.
use crate::geom::Point;

/// Dekker's splitter, `2^27 + 1`, which splits a double into two 26-bit halves.
const SPLITTER: f64 = 134_217_729.0;

/// A double-double: high and low parts of an unevaluated sum.
#[derive(Clone, Copy, Debug)]
struct Dd(f64, f64);
impl Dd {
    fn add(self, b: Self) -> Self {
        let s = self.0 + b.0;
        let t = self.1 + b.1;
        let e = s - self.0;
        let f = t - self.1;
        let low_s = (b.0 - e) + (self.0 - (s - e));
        let low_t = (b.1 - f) + (self.1 - (t - f));
        let e = low_s + t;
        let h = s + e;
        let low_h = e + (s - h);
        let e = low_t + low_h;
        let hi = h + e;
        Self(hi, e + (h - hi))
    }
    fn sub(self, b: Self) -> Self {
        self.add(Self(-b.0, -b.1))
    }
    fn mul(self, b: Self) -> Self {
        let c = SPLITTER * self.0;
        let hx = c - (c - self.0);
        let tx = self.0 - hx;
        let c = SPLITTER * b.0;
        let hy = c - (c - b.0);
        let ty = b.0 - hy;
        let product = self.0 * b.0;
        let c =
            ((((hx * hy - product) + hx * ty) + tx * hy) + tx * ty) + (self.0 * b.1 + self.1 * b.0);
        let hi = product + c;
        Self(hi, c + (product - hi))
    }
}
fn sign(x: f64) -> i32 {
    if x < 0.0 { -1 } else { i32::from(x > 0.0) }
}
/// Orientation of `c` relative to the directed line `a`→`b`: `1` counterclockwise,
/// `-1` clockwise, `0` collinear. Returned as a sign so callers can multiply results.
pub(crate) fn orient(a: Point, b: Point, c: Point) -> i32 {
    let left = (a.x - c.x) * (b.y - c.y);
    let right = (a.y - c.y) * (b.x - c.x);
    let det = left - right;
    let sum = if left > 0.0 {
        if right <= 0.0 {
            return sign(det);
        }
        left + right
    } else if left < 0.0 {
        if right >= 0.0 {
            return sign(det);
        }
        -left - right
    } else {
        return sign(det);
    };
    let bound = 1e-15 * sum;
    if det >= bound || -det >= bound {
        return sign(det);
    }
    let dx1 = Dd(b.x, 0.0).add(Dd(-a.x, 0.0));
    let dy1 = Dd(b.y, 0.0).add(Dd(-a.y, 0.0));
    let dx2 = Dd(c.x, 0.0).add(Dd(-b.x, 0.0));
    let dy2 = Dd(c.y, 0.0).add(Dd(-b.y, 0.0));
    let d = dx1.mul(dy2).sub(dy1.mul(dx2));
    if d.0 < 0.0 || (d.0 == 0.0 && d.1 < 0.0) {
        -1
    } else {
        i32::from(d.0 > 0.0 || (d.0 == 0.0 && d.1 > 0.0))
    }
}
/// Whether `p` lies inside or on triangle `a`, `b`, `c`.
pub(crate) fn triangle_contains(a: Point, b: Point, c: Point, p: Point) -> bool {
    let exterior = if orient(a, b, c) > 0 { -1 } else { 1 };
    orient(a, b, p) != exterior && orient(b, c, p) != exterior && orient(c, a, p) != exterior
}
/// Area of triangle `a`, `b`, `c`.
pub(crate) fn triangle_area(a: Point, b: Point, c: Point) -> f64 {
    (((c.x - a.x) * (b.y - a.y) - (b.x - a.x) * (c.y - a.y)) / 2.0).abs()
}
