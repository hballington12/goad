use geo::algorithm::simplify::Simplify;
use geo::{polygon, CoordsIter, MultiPolygon};

fn main() {
    // Define a MultiPolygon
    let multipolygon = MultiPolygon(vec![polygon![
        (x: 0.0, y: 0.0),
        (x: 5.0, y: 0.0),
        (x: 5.0, y: 5.0),
        (x: 5.0, y: 4.99),
        (x: 0.0, y: 5.0),
        (x: 0.0, y: 0.0),
    ]]);

    println!("Original MultiPolygon: {:?}", multipolygon);

    let cleaned = multipolygon.simplify(&0.01);

    // Print the cleaned polygon
    println!("Cleaned MultiPolygon: {:?}", cleaned);

    // Assert that the number of vertices in the cleaned exterior is 5
    let cleaned_exterior = &cleaned.0[0].exterior();
    assert_eq!(cleaned_exterior.coords_count(), 5);
}
