pub fn main() {
    println!("Starting rusty convex hull 2d example ...");

    let input_set = [
        rusty_convex_hull::Point2D { x: -1.0, y: -10.0 },
        rusty_convex_hull::Point2D { x:  0.0, y:   0.0 },
        rusty_convex_hull::Point2D { x: -3.0, y:   5.0 },
        rusty_convex_hull::Point2D { x:  6.0, y:   1.0 },
    ];
    let hull = rusty_convex_hull::convex_hull_graham::convex_hull_graham(&input_set);
    println!("{:?}", hull);
}
