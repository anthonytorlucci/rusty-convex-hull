//! Graham Scan convex hull algorithm for 2D point sets.
//!
//! Provides an implementation of the Graham Scan algorithm to compute the
//! convex hull of a 2D point set represented as `(f64, f64)` tuples.

use std::cmp::Ordering::Equal;

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
fn are_collinear(p: &(f64, f64), q: &(f64, f64), r: &(f64, f64)) -> bool {
    let slope12 = (q.1 - p.1) / (q.0 - p.0);
    let slope23 = (r.1 - q.1) / (r.0 - q.0);
    within_tolerance(slope12, slope23, 1e-12)
}

#[allow(dead_code)]
/// Removes collinear interior points from a sorted point sequence in place.
///
/// For each consecutive triple a → b → c, removes b if the three points are
/// collinear. Assumes `points` is sorted so that the collinear point is always
/// the middle element of the triple.
fn remove_collinear_points(points: &mut Vec<(f64, f64)>) {
    let mut i = 0;
    while i + 2 < points.len() {
        if are_collinear(&points[i], &points[i + 1], &points[i + 2]) {
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
fn sort_by_min_angle(pts: &[(f64, f64)], min: &(f64, f64)) -> Vec<(f64, f64)> {
    let mut points: Vec<(f64, f64, (f64, f64))> = pts
        .iter()
        .map(|x| {
            (
                (x.1 - min.1).atan2(x.0 - min.0),
                // angle
                (x.1 - min.1).hypot(x.0 - min.0),
                // distance (we want the closest to be first)
                *x,
            )
        })
        .collect();
    points.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Equal));
    points.into_iter().map(|x| x.2).collect()
}

/// Computes the z-component of the cross product of vectors pq and pr.
///
/// A positive value indicates a counter-clockwise turn at q, zero indicates
/// collinearity, and a negative value indicates a clockwise turn.
fn calc_z_coord_vector_product(p: &(f64, f64), q: &(f64, f64), r: &(f64, f64)) -> f64 {
    (q.0 - p.0) * (r.1 - p.1) - (r.0 - p.0) * (q.1 - p.1)
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
/// use rusty_convex_hull::convex_hull_graham::convex_hull_graham;
///
/// let points = [(-1.0_f64, -10.0_f64), (0.0, 0.0), (-3.0, 5.0), (6.0, 1.0)];
/// let hull = convex_hull_graham(&points);
/// assert_eq!(hull, vec![(-1.0, -10.0), (6.0, 1.0), (-3.0, 5.0)]);
/// ```
pub fn convex_hull_graham(pts: &[(f64, f64)]) -> Vec<(f64, f64)> {
    if pts.is_empty() {
        return vec![];
    }

    let mut stack: Vec<(f64, f64)> = vec![];
    let min = pts
        .iter()
        .min_by(|a, b| {
            let ord = a.1.partial_cmp(&b.1).unwrap_or(Equal);
            match ord {
                Equal => a.0.partial_cmp(&b.0).unwrap_or(Equal),
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
            && calc_z_coord_vector_product(&stack[stack.len() - 2], &stack[stack.len() - 1], &point)
                < 0.
        {
            stack.pop();
        }
        stack.push(point);
    }

    stack
}

// TODO: is there a better way to structure this? pure functions, structs, or even traits? For example Point2D {x,y}, Segment {Point2D, Point2D}
// elaborate on the purpose of the cross product

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
        assert!(are_collinear(&(0.0, 0.0), &(1.0, 1.0), &(2.0, 2.0)));
    }

    #[test]
    fn test_are_collinear_horizontal() {
        assert!(are_collinear(&(0.0, 0.0), &(1.0, 0.0), &(2.0, 0.0)));
    }

    #[test]
    fn test_are_collinear_false() {
        assert!(!are_collinear(&(0.0, 0.0), &(1.0, 0.0), &(0.0, 1.0)));
    }

    // --- remove_collinear_points ---

    #[test]
    fn test_remove_collinear_points_removes_middle() {
        let mut pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)];
        remove_collinear_points(&mut pts);
        assert_eq!(pts, vec![(0.0, 0.0), (2.0, 2.0)]);
    }

    #[test]
    fn test_remove_collinear_points_chain() {
        let mut pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)];
        remove_collinear_points(&mut pts);
        assert_eq!(pts, vec![(0.0, 0.0), (3.0, 3.0)]);
    }

    #[test]
    fn test_remove_collinear_points_no_change() {
        let mut pts = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        remove_collinear_points(&mut pts);
        assert_eq!(pts, vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]);
    }

    // --- sort_by_min_angle ---

    #[test]
    fn test_sort_by_min_angle_orders_by_angle() {
        let min = (0.0, 0.0);
        // angles: (1,0)→0, (1,1)→π/4, (0,1)→π/2
        let pts = [(1.0, 1.0), (0.0, 1.0), (1.0, 0.0)];
        let sorted = sort_by_min_angle(&pts, &min);
        assert_eq!(sorted, vec![(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]);
    }

    #[test]
    fn test_sort_by_min_angle_tie_broken_by_distance() {
        let min = (0.0, 0.0);
        let pts = [(2.0, 0.0), (1.0, 0.0)];
        let sorted = sort_by_min_angle(&pts, &min);
        assert_eq!(sorted, vec![(1.0, 0.0), (2.0, 0.0)]);
    }

    // --- calc_z_coord_vector_product ---

    #[test]
    fn test_calc_z_ccw_is_positive() {
        // p→q→r is a left (CCW) turn
        let z = calc_z_coord_vector_product(&(0.0, 0.0), &(1.0, 0.0), &(1.0, 1.0));
        assert!(z > 0.0);
    }

    #[test]
    fn test_calc_z_cw_is_negative() {
        // p→q→r is a right (CW) turn
        let z = calc_z_coord_vector_product(&(0.0, 0.0), &(1.0, 1.0), &(1.0, 0.0));
        assert!(z < 0.0);
    }

    #[test]
    fn test_calc_z_collinear_is_zero() {
        let z = calc_z_coord_vector_product(&(0.0, 0.0), &(1.0, 1.0), &(2.0, 2.0));
        assert_eq!(z, 0.0);
    }

    // --- convex_hull_graham ---

    #[test]
    fn test_hull_empty_input() {
        let hull = convex_hull_graham(&[]);
        assert!(hull.is_empty());
    }

    #[test]
    fn test_hull_single_point() {
        let hull = convex_hull_graham(&[(3.0, 4.0)]);
        assert_eq!(hull, vec![(3.0, 4.0)]);
    }

    #[test]
    fn test_hull_three_points_returned_as_is() {
        let pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
        let hull = convex_hull_graham(&pts);
        assert_eq!(hull.len(), 3);
    }

    #[test]
    fn test_hull_excludes_interior_point() {
        let pts = [(-1.0, -10.0), (0.0, 0.0), (-3.0, 5.0), (6.0, 1.0)];
        let hull = convex_hull_graham(&pts);
        assert_eq!(hull, vec![(-1.0, -10.0), (6.0, 1.0), (-3.0, 5.0)]);
        assert!(!hull.contains(&(0.0, 0.0)));
    }

    #[test]
    fn test_hull_square_with_interior_point() {
        let pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5)];
        let hull = convex_hull_graham(&pts);
        assert_eq!(hull.len(), 4);
        assert!(!hull.contains(&(0.5, 0.5)));
    }
}

