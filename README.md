# A Rusty Convex Hull in 2D

Imagine you have a scattered set of data points, like cities on a map, stars in the sky, or even measurements from an experiment. You have been given the task to find the smallest convex shape that completely encloses all of these points. This encompassing boundary is known as the *convex hull*.

Implementing a robust and efficient solution for calculating the convex hull is a cornerstone skill for any programmer dealing with geometric algorithms. This article will delve into the details of implementing the Graham Scan algorithm, exploring the core steps and computational considerations necessary to transform a cloud of points into its defining geometric boundary.

> [!NOTE]
> We only consider a collection of 2D points on a Cartesian grid for this example. Higher dimensional convex hull solutions exist, but this implementation is strictly 2D.

## Problem Statement
Given a list of 2D points, compute the **convex hull** of these points. The task is to return the vertices of the bounding polygon as a subset of the input points, ordered counter-clockwise (CCW) to form a closed loop.

Points are represented by the [`Point2D`] struct rather than bare tuples, which keeps the geometry explicit at the call site:

```rust
use rusty_convex_hull::Point2D;

let input_set = [
    Point2D { x: -1.0, y: -10.0 },
    Point2D { x:  0.0, y:   0.0 },
    Point2D { x: -3.0, y:   5.0 },
    Point2D { x:  6.0, y:   1.0 },
];
```

`Point2D` is `Copy` and implements `From<(f64, f64)>` / `Into<(f64, f64)>`, so you can still build it from tuples (`Point2D::from((1.0, 2.0))`) when that is more convenient.

## Solution

### Convex Set and Convex Hull

A set of points $C$ is called **convex** if for any two points $p$ and $q$ in $C$, the entire segment $[pq]$ is contained in $C$. The **convex hull** of a set $S$ is the smallest convex set containing $S$ — equivalently, the *intersection of all convex sets that contain $S$*. The two definitions are provably equivalent.^[1](#ref-item-1) For more background on convex sets, see [[2](#ref-item-2)].

We can think of a convex hull like a rubber band wrapped around all the points without any inward bends (every interior angle is less than 180 degrees). For a finite point set that is not entirely collinear, the hull boundary is also the polygon of **minimum perimeter** enclosing every point.^[3](#ref-item-3)

> [!NOTE]
> By the canonical convention, a point that lies *on* a hull edge strictly between two vertices is a collinear interior point and is **not** reported as a hull vertex — only the segment endpoints are. This implementation follows that convention (see the implementation note below).

### Graham Scan Algorithm

The Graham Scan is arguably the most famous textbook convex hull algorithm, primarily because it clearly demonstrates the power of sorting combined with geometric checks.^[4](#ref-item-4)

* **Description:**
  1. Find the point with the minimum y-coordinate (the **anchor point** $P\_0$), breaking ties by the minimum x-coordinate. This bottom-most point is guaranteed to lie on the hull.
  2. **Sort** all other points by the **polar angle** they make with $P\_0$ and the positive x-axis. This is the key step and can be done efficiently using the cross-product to avoid explicit trigonometric angle calculations. When several points share the same angle, ties must be broken consistently by distance from $P\_0$.
  3. Iterate through the sorted points using a stack. For each new point, check whether the last three points (second-to-last, last, new) make a **left turn**, tested via the sign of the cross product.
     * On a **right turn** (a non-left turn), the last point on the stack lies inside the hull and is popped.
     * On a **left turn**, the new point is pushed onto the stack.
* **Ease of implementation:** The stack-based scan is straightforward; the **angular sort** is the subtle part. It requires a correct choice of anchor point $P\_0$, measuring angles implicitly via cross-products, and consistent handling of collinear ties. If done incorrectly, the sorted order is wrong and the resulting hull is wrong.

#### Why the algorithm is correct and efficient

The anchor $P\_0$ is the bottom-most (then left-most) point, so it must be a hull vertex. Sorting the remaining points by polar angle around $P\_0$ lays them out in counter-clockwise order. The single stack pass then keeps only the points that preserve a left turn at every step; any point that would create a right turn is an interior point and is discarded. When the scan finishes, the stack holds exactly the hull vertices in counter-clockwise order.

* **Time complexity — $O(n \log n)$.** The cost is dominated by the angular sort. The stack scan is $O(n)$ because each point is pushed and popped at most once.
* **Space complexity — $O(n)$** for the sorted points and the output stack.
* **Optimality.** $O(n \log n)$ is asymptotically optimal: any algorithm for the 2D convex hull problem in the *algebraic decision tree* model requires $\Omega(n \log n)$ time in the worst case. This follows from a reduction from sorting — mapping $n$ numbers onto the parabola $y = x^2$ shows that computing their hull also sorts them — first proved by Yao (1981).^[5](#ref-item-5)

#### Implementation note

This crate's orientation test pops on a **non-strict** comparison (`cross_z(...) <= 0.0`), so a point that turns clockwise *or* lies exactly collinear with the current edge is discarded. As a result, only the endpoints of each edge are reported and collinear interior-of-edge points are excluded — the canonical convex hull. (The alternative, a strict `< 0.0`, would instead retain collinear boundary points; both are valid Graham Scans and differ only in this respect.) Collinearity is decided by the sign of the cross product, which avoids the division and slope special-cases that a slope-comparison approach would require for vertical edges.

## Usage

Run the bundled example, which generates a random 2D point cloud inside a disk, shuffles it into arbitrary order (the algorithm is order-independent), computes the hull, and renders it as ASCII art:

```sh
cargo run --example convex_hull_2d
```

The example also reports geometric statistics — hull vertex count, area (via the shoelace formula), perimeter, and how much of the bounding disk the hull covers:

```text
Generated 120 points in a disk of radius 20.0.
Hull has 17 vertices.

Hull area:      1098.732
Hull perimeter: 119.764
Disk coverage:  87.4% of the bounding disk area

Plot ('.' = point, 'o' = hull edge, '#' = hull vertex):

                                  oooooo#ooo
                      oooooooooooo          ooo#oo
               #oooooo                            oooo
             oo            .                ..        oooo
           oo      .              .                       o#o#o
          o                                          .       . o#o
        oo       .                .       .         .             o#
      oo                .         .                 .              o
     o  .                          .                  .     .       o
   oo                                                     .       ..o
 ...
```

Using the library directly is a single call:

```rust
use rusty_convex_hull::{Point2D, convex_hull_graham::convex_hull_graham};

let points = [
    Point2D { x: -1.0, y: -10.0 },
    Point2D { x:  0.0, y:   0.0 },
    Point2D { x: -3.0, y:   5.0 },
    Point2D { x:  6.0, y:   1.0 },
];
let hull = convex_hull_graham(&points);
// Returns the hull vertices in counter-clockwise order, starting from the
// lowest point; the interior point (0.0, 0.0) is excluded.
```

## References {#references}
1. M. de Berg, O. Cheong, M. van Kreveld, M. Overmars. *Computational Geometry: Algorithms and Applications*, 3rd ed. Springer, 2008. ISBN 978-3-540-77973-5. (Chapter 1 covers the convex hull definitions and their equivalence.) <span id="ref-item-1"></span>
2. Wikipedia — [Convex set](https://en.wikipedia.org/wiki/Convex_set) <span id="ref-item-2"></span>
3. Wikipedia — [Convex hull](https://en.wikipedia.org/wiki/Convex_hull) (minimum-perimeter property) <span id="ref-item-3"></span>
4. R. L. Graham. "An Efficient Algorithm for Determining the Convex Hull of a Finite Planar Set." *Information Processing Letters*, 1(4):132–133, 1972. DOI: [10.1016/0020-0190(72)90045-2](https://doi.org/10.1016/0020-0190(72)90045-2) <span id="ref-item-4"></span>
5. A. C.-C. Yao. "A Lower Bound to Finding Convex Hulls." *Journal of the ACM*, 28(4):780–787, 1981. DOI: [10.1145/322276.322289](https://doi.org/10.1145/322276.322289) <span id="ref-item-5"></span>

---
### Additional Sources of information
* T. H. Cormen, C. E. Leiserson, R. L. Rivest, C. Stein. *Introduction to Algorithms*, 3rd ed., §33.3 "Finding the Convex Hull." MIT Press, 2009. ISBN 978-0-262-03384-8.
* F. P. Preparata, M. I. Shamos. *Computational Geometry: An Introduction*. Springer, 1985. ISBN 978-0-387-96131-6. (Theorem 3.3 proves the $\Omega(n \log n)$ lower bound and details collinear-degeneracy handling.)
* Algorithmic Solutions LEDA Guide — [Convex Hulls](https://www.algorithmic-solutions.info/leda_guide/geo_algs/convex_hull.html)
* [The Algorithms — convex_hull.rs](https://github.com/TheAlgorithms/Rust/blob/master/src/general/convex_hull.rs)
* [YouTube — Convex Hull Algorithms](https://www.youtube.com/watch?v=fTqPVjy0rzU) (part of a course on computational geometry)

[`Point2D`]: src/point.rs
