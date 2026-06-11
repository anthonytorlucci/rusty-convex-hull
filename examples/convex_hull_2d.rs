//! Convex hull of a randomly scattered 2D point cloud.
//!
//! This example generates a cloud of points scattered inside a disk, shuffles
//! them into a random order (the Graham Scan does not care about input order),
//! computes the convex hull, reports a few geometric statistics, and renders an
//! ASCII plot with the hull boundary drawn around the interior points.
//!
//! It is fully self-contained — a tiny `xorshift` PRNG stands in for the `rand`
//! crate so the example pulls in no dependencies and is reproducible run to run.

use rusty_convex_hull::{Point2D, convex_hull_graham::convex_hull_graham};

/// Minimal `xorshift64*` pseudo-random generator.
///
/// Deterministic given a seed, which keeps the example reproducible while still
/// producing a convincingly scattered cloud.
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        // Avoid the all-zero state, which xorshift cannot escape.
        Self {
            state: seed | 1,
        }
    }

    /// Returns the next raw 64-bit value.
    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Returns a float in the half-open range `[0.0, 1.0)`.
    fn next_f64(&mut self) -> f64 {
        // Use the top 53 bits — the full mantissa of an f64.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Returns an integer in `[0, n)`.
    fn next_usize(&mut self, n: usize) -> usize {
        (self.next_f64() * n as f64) as usize
    }
}

/// Generates `count` points scattered uniformly inside a disk.
///
/// Sampling the radius as `sqrt(u)` keeps the points uniform over the area
/// rather than clustering them near the center.
fn random_disk_cloud(rng: &mut Rng, count: usize, center: Point2D, radius: f64) -> Vec<Point2D> {
    (0..count)
        .map(|_| {
            let angle = rng.next_f64() * std::f64::consts::TAU;
            let r = radius * rng.next_f64().sqrt();
            Point2D {
                x: center.x + r * angle.cos(),
                y: center.y + r * angle.sin(),
            }
        })
        .collect()
}

/// Shuffles a slice in place using the Fisher–Yates algorithm.
fn shuffle(rng: &mut Rng, points: &mut [Point2D]) {
    for i in (1..points.len()).rev() {
        let j = rng.next_usize(i + 1);
        points.swap(i, j);
    }
}

/// Returns the signed area of a polygon via the shoelace formula.
///
/// Positive for counter-clockwise vertex order, as produced by the hull.
fn polygon_area(poly: &[Point2D]) -> f64 {
    if poly.len() < 3 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        sum += a.x * b.y - b.x * a.y;
    }
    sum / 2.0
}

/// Returns the total edge length around a closed polygon.
fn polygon_perimeter(poly: &[Point2D]) -> f64 {
    if poly.len() < 2 {
        return 0.0;
    }
    let mut total = 0.0;
    for i in 0..poly.len() {
        let a = poly[i];
        let b = poly[(i + 1) % poly.len()];
        total += (b.x - a.x).hypot(b.y - a.y);
    }
    total
}

/// Renders the cloud and its hull as ASCII art.
///
/// Interior points are drawn as `.`, hull edges as `o`, and hull vertices as
/// `#`. The plot is flipped vertically so that +y points up, as on paper.
fn render(cloud: &[Point2D], hull: &[Point2D], width: usize, height: usize) -> String {
    // Determine the bounding box of the whole cloud.
    let (mut min_x, mut min_y) = (f64::INFINITY, f64::INFINITY);
    let (mut max_x, mut max_y) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in cloud {
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
    }
    let span_x = (max_x - min_x).max(1e-9);
    let span_y = (max_y - min_y).max(1e-9);

    // Map a point into integer grid coordinates.
    let to_cell = |p: Point2D| -> (i32, i32) {
        let col = ((p.x - min_x) / span_x * (width as f64 - 1.0)).round() as i32;
        // Flip y so larger y is higher on screen.
        let row = ((max_y - p.y) / span_y * (height as f64 - 1.0)).round() as i32;
        (col, row)
    };

    let mut grid = vec![vec![' '; width]; height];
    let plot = |grid: &mut Vec<Vec<char>>, col: i32, row: i32, ch: char| {
        if (0..width as i32).contains(&col) && (0..height as i32).contains(&row) {
            grid[row as usize][col as usize] = ch;
        }
    };

    // Interior cloud points first, so hull marks draw on top.
    for &p in cloud {
        let (c, r) = to_cell(p);
        plot(&mut grid, c, r, '.');
    }

    // Hull edges via a simple Bresenham line between consecutive vertices.
    for i in 0..hull.len() {
        let (x0, y0) = to_cell(hull[i]);
        let (x1, y1) = to_cell(hull[(i + 1) % hull.len()]);
        let (dx, dy) = ((x1 - x0).abs(), -(y1 - y0).abs());
        let (sx, sy) = (if x0 < x1 { 1 } else { -1 }, if y0 < y1 { 1 } else { -1 });
        let (mut x, mut y, mut err) = (x0, y0, dx + dy);
        loop {
            plot(&mut grid, x, y, 'o');
            if x == x1 && y == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x += sx;
            }
            if e2 <= dx {
                err += dx;
                y += sy;
            }
        }
    }

    // Hull vertices last so they are never overwritten.
    for &p in hull {
        let (c, r) = to_cell(p);
        plot(&mut grid, c, r, '#');
    }

    grid.into_iter()
        .map(|row| row.into_iter().collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

pub fn main() {
    println!("Starting rusty convex hull 2d example ...\n");

    // Fixed seed → reproducible cloud. Change it to explore other clouds.
    let mut rng = Rng::new(42);

    let center = Point2D { x: 0.0, y: 0.0 };
    let radius = 20.0;
    let count = 120;

    let mut cloud = random_disk_cloud(&mut rng, count, center, radius);

    // The algorithm is order-independent — shuffle to prove the point.
    shuffle(&mut rng, &mut cloud);

    let hull = convex_hull_graham(&cloud);

    let area = polygon_area(&hull);
    let perimeter = polygon_perimeter(&hull);

    println!("Generated {} points in a disk of radius {radius:.1}.", cloud.len());
    println!("Hull has {} vertices.\n", hull.len());

    println!("Hull vertices (counter-clockwise from lowest point):");
    for (i, p) in hull.iter().enumerate() {
        println!("  {i:>2}: ({:>7.3}, {:>7.3})", p.x, p.y);
    }

    println!("\nHull area:      {area:.3}");
    println!("Hull perimeter: {perimeter:.3}");
    println!(
        "Disk coverage:  {:.1}% of the bounding disk area\n",
        area / (std::f64::consts::PI * radius * radius) * 100.0
    );

    println!("Plot ('.' = point, 'o' = hull edge, '#' = hull vertex):\n");
    println!("{}", render(&cloud, &hull, 72, 28));
}
