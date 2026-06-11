//! Core 2D point type used throughout this crate.

/// A point in 2D Euclidean space.
///
/// `Point2D` is `Copy` — pass it by value everywhere; no references needed.
///
/// # Examples
///
/// ```
/// use rusty_convex_hull::Point2D;
///
/// let p = Point2D { x: 1.0, y: 2.0 };
/// let q = Point2D::from((3.0_f64, 4.0_f64));
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl From<(f64, f64)> for Point2D {
    fn from((x, y): (f64, f64)) -> Self {
        Self { x, y }
    }
}

impl From<Point2D> for (f64, f64) {
    fn from(p: Point2D) -> Self {
        (p.x, p.y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_tuple() {
        let p = Point2D::from((1.0, 2.0));
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
    }

    #[test]
    fn test_into_tuple() {
        let p = Point2D { x: 3.0, y: 4.0 };
        let t: (f64, f64) = p.into();
        assert_eq!(t, (3.0, 4.0));
    }

    #[test]
    fn test_copy_semantics() {
        let p = Point2D { x: 1.0, y: 2.0 };
        let q = p; // copy, not move
        assert_eq!(p, q);
    }
}
