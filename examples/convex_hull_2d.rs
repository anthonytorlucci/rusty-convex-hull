pub fn main() {
    println!("Starting rusty convex hull 2d example ...");

    let input_set = [(-1.0, -10.0), (0.0, 0.0), (-3.0, 5.0), (6.0, 1.0)];
    let hull: Vec<(f64, f64)> =
        rusty_convex_hull::convex_hull_graham::convex_hull_graham(&input_set);
    println!("{:?}", hull);
}
