//! Graham Scan convex hull algorithm for 2D point sets.
//!
//! Provides an implementation of the Graham Scan algorithm to compute the
//! convex hull of a 2D point set represented as [`Point2D`] values.

use std::cmp::Ordering::Equal;

use crate::Point2D;

// helper functions
#[allow(dead_code)]
/// Returns `true` if two values differ by no more than `tol`.
fn within_tolerance(p: f64, q: f64, tol: f64) -> bool {
    (p - q).abs() <= tol
}

#[allow(dead_code)]
/// Returns `true` if three 2D points are collinear within floating-point tolerance.
///
/// Compares slopes of segments pq and qr using a tolerance of 1e-12 to handle
/// floating-point rounding error.
fn are_collinear(p: Point2D, q: Point2D, r: Point2D) -> bool {
    let slope12 = (q.y - p.y) / (q.x - p.x);
    let slope23 = (r.y - q.y) / (r.x - q.x);
    within_tolerance(slope12, slope23, 1e-12)
}

#[allow(dead_code)]
/// Removes collinear interior points from a sorted point sequence in place.
///
/// For each consecutive triple a → b → c, removes b if the three points are
/// collinear. Assumes `points` is sorted so that the collinear point is always
/// the middle element of the triple.
fn remove_collinear_points(points: &mut Vec<Point2D>) {
    let mut i = 0;
    while i + 2 < points.len() {
        if are_collinear(points[i], points[i + 1], points[i + 2]) {
            points.remove(i + 1);
        } else {
            i += 1;
        }
    }
}

/// Sorts points by polar angle relative to `min`, breaking ties by distance.
///
/// Returns a new vector ordered counter-clockwise from `min`, with closer
/// points appearing first when angles are equal.
fn sort_by_min_angle(pts: &[Point2D], min: Point2D) -> Vec<Point2D> {
    let mut points: Vec<(f64, f64, Point2D)> = pts
        .iter()
        .map(|&p| {
            (
                (p.y - min.y).atan2(p.x - min.x), // angle
                (p.y - min.y).hypot(p.x - min.x), // distance (closest first)
                p,
            )
        })
        .collect();
    points.sort_by(|a, b| (a.0, a.1).partial_cmp(&(b.0, b.1)).unwrap_or(Equal));
    points.into_iter().map(|t| t.2).collect()
}

/// Computes the z-component of the cross product of vectors pq and pr.
///
/// A positive value indicates a counter-clockwise turn at q, zero indicates
/// collinearity, and a negative value indicates a clockwise turn.
fn cross_z(p: Point2D, q: Point2D, r: Point2D) -> f64 {
    (q.x - p.x) * (r.y - p.y) - (r.x - p.x) * (q.y - p.y)
}

/// Computes the convex hull of a 2D point set using the Graham Scan algorithm.
///
/// Returns hull vertices in counter-clockwise order, starting from the point
/// with the lowest y-coordinate (x-coordinate breaks ties). Returns an empty
/// vector when `pts` is empty.
///
/// # Examples
///
/// ```
/// use rusty_convex_hull::{Point2D, convex_hull_graham::convex_hull_graham};
///
/// let points = [
///     Point2D { x: -1.0, y: -10.0 },
///     Point2D { x:  0.0, y:   0.0 },
///     Point2D { x: -3.0, y:   5.0 },
///     Point2D { x:  6.0, y:   1.0 },
/// ];
/// let hull = convex_hull_graham(&points);
/// assert_eq!(hull, vec![
///     Point2D { x: -1.0, y: -10.0 },
///     Point2D { x:  6.0, y:   1.0 },
///     Point2D { x: -3.0, y:   5.0 },
/// ]);
/// ```
pub fn convex_hull_graham(pts: &[Point2D]) -> Vec<Point2D> {
    if pts.is_empty() {
        return vec![];
    }

    let mut stack: Vec<Point2D> = vec![];
    let &min = pts
        .iter()
        .min_by(|a, b| {
            let ord = a.y.partial_cmp(&b.y).unwrap_or(Equal);
            match ord {
                Equal => a.x.partial_cmp(&b.x).unwrap_or(Equal),
                o => o,
            }
        })
        .unwrap();
    let points = sort_by_min_angle(pts, min);

    if points.len() <= 3 {
        return points;
    }

    for point in points {
        while stack.len() > 1
            && cross_z(stack[stack.len() - 2], stack[stack.len() - 1], point) < 0.0
        {
            stack.pop();
        }
        stack.push(point);
    }

    stack
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- within_tolerance ---

    #[test]
    fn test_within_tolerance_equal() {
        assert!(within_tolerance(1.0, 1.0, 1e-12));
    }

    #[test]
    fn test_within_tolerance_inside() {
        assert!(within_tolerance(1.0, 1.0 + 1e-13, 1e-12));
    }

    #[test]
    fn test_within_tolerance_near_boundary() {
        // 0.5e-12 is well inside the 1e-12 tolerance without hitting f64 rounding
        assert!(within_tolerance(1.0, 1.0 + 0.5e-12, 1e-12));
    }

    #[test]
    fn test_within_tolerance_outside() {
        assert!(!within_tolerance(1.0, 2.0, 1e-12));
    }

    // --- are_collinear ---

    #[test]
    fn test_are_collinear_diagonal() {
        assert!(are_collinear(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 2.0 },
        ));
    }

    #[test]
    fn test_are_collinear_horizontal() {
        assert!(are_collinear(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 2.0, y: 0.0 },
        ));
    }

    #[test]
    fn test_are_collinear_false() {
        assert!(!are_collinear(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 0.0, y: 1.0 },
        ));
    }

    // --- remove_collinear_points ---

    #[test]
    fn test_remove_collinear_points_removes_middle() {
        let mut pts = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 2.0 },
        ];
        remove_collinear_points(&mut pts);
        assert_eq!(
            pts,
            vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 2.0, y: 2.0 }]
        );
    }

    #[test]
    fn test_remove_collinear_points_chain() {
        let mut pts = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 2.0 },
            Point2D { x: 3.0, y: 3.0 },
        ];
        remove_collinear_points(&mut pts);
        assert_eq!(
            pts,
            vec![Point2D { x: 0.0, y: 0.0 }, Point2D { x: 3.0, y: 3.0 }]
        );
    }

    #[test]
    fn test_remove_collinear_points_no_change() {
        let mut pts = vec![
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 0.0, y: 1.0 },
        ];
        remove_collinear_points(&mut pts);
        assert_eq!(
            pts,
            vec![
                Point2D { x: 0.0, y: 0.0 },
                Point2D { x: 1.0, y: 0.0 },
                Point2D { x: 0.0, y: 1.0 },
            ]
        );
    }

    // --- sort_by_min_angle ---

    #[test]
    fn test_sort_by_min_angle_orders_by_angle() {
        let min = Point2D { x: 0.0, y: 0.0 };
        // angles: (1,0)→0, (1,1)→π/4, (0,1)→π/2
        let pts = [
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 0.0, y: 1.0 },
            Point2D { x: 1.0, y: 0.0 },
        ];
        let sorted = sort_by_min_angle(&pts, min);
        assert_eq!(
            sorted,
            vec![
                Point2D { x: 1.0, y: 0.0 },
                Point2D { x: 1.0, y: 1.0 },
                Point2D { x: 0.0, y: 1.0 },
            ]
        );
    }

    #[test]
    fn test_sort_by_min_angle_tie_broken_by_distance() {
        let min = Point2D { x: 0.0, y: 0.0 };
        let pts = [Point2D { x: 2.0, y: 0.0 }, Point2D { x: 1.0, y: 0.0 }];
        let sorted = sort_by_min_angle(&pts, min);
        assert_eq!(
            sorted,
            vec![Point2D { x: 1.0, y: 0.0 }, Point2D { x: 2.0, y: 0.0 }]
        );
    }

    // --- cross_z ---

    #[test]
    fn test_cross_z_ccw_is_positive() {
        let z = cross_z(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
        );
        assert!(z > 0.0);
    }

    #[test]
    fn test_cross_z_cw_is_negative() {
        let z = cross_z(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 1.0, y: 0.0 },
        );
        assert!(z < 0.0);
    }

    #[test]
    fn test_cross_z_collinear_is_zero() {
        let z = cross_z(
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 2.0, y: 2.0 },
        );
        assert_eq!(z, 0.0);
    }

    // --- convex_hull_graham ---

    #[test]
    fn test_hull_empty_input() {
        assert!(convex_hull_graham(&[]).is_empty());
    }

    #[test]
    fn test_hull_single_point() {
        let p = Point2D { x: 3.0, y: 4.0 };
        assert_eq!(convex_hull_graham(&[p]), vec![p]);
    }

    #[test]
    fn test_hull_three_points_returned_as_is() {
        let pts = [
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 0.0, y: 1.0 },
        ];
        assert_eq!(convex_hull_graham(&pts).len(), 3);
    }

    #[test]
    fn test_hull_excludes_interior_point() {
        let pts = [
            Point2D { x: -1.0, y: -10.0 },
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: -3.0, y: 5.0 },
            Point2D { x: 6.0, y: 1.0 },
        ];
        let hull = convex_hull_graham(&pts);
        assert_eq!(
            hull,
            vec![
                Point2D { x: -1.0, y: -10.0 },
                Point2D { x: 6.0, y: 1.0 },
                Point2D { x: -3.0, y: 5.0 },
            ]
        );
        assert!(!hull.contains(&Point2D { x: 0.0, y: 0.0 }));
    }

    #[test]
    fn test_hull_square_with_interior_point() {
        let pts = [
            Point2D { x: 0.0, y: 0.0 },
            Point2D { x: 1.0, y: 0.0 },
            Point2D { x: 1.0, y: 1.0 },
            Point2D { x: 0.0, y: 1.0 },
            Point2D { x: 0.5, y: 0.5 },
        ];
        let hull = convex_hull_graham(&pts);
        assert_eq!(hull.len(), 4);
        assert!(!hull.contains(&Point2D { x: 0.5, y: 0.5 }));
    }
}
