use geo::{polygon, CoordsIter, MultiLineString, MultiPolygon, Polygon};
use geo_clipper::Clipper;
use pbt::config;

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

    let cleaned_multipolygon = simplify(&multipolygon);

    // Print the cleaned polygon
    println!("Original MultiPolygon: {:?}", multipolygon);
    println!("Cleaned Polygon: {:?}", cleaned_multipolygon);

    // Assert that the number of vertices in the cleaned exterior is 5
    let cleaned_exterior = &cleaned_multipolygon.0[0].exterior();
    assert_eq!(cleaned_exterior.coords_count(), 5);
}

fn simplify(multipolygon: &MultiPolygon) -> MultiPolygon {
    // Clean the MultiPolygon, which returns a MultiLineString
    let cleaned_multipolygon: MultiLineString<f64> = multipolygon.simplify(
        geo_clipper::PolyFillType::EvenOdd,
        config::CLIP_TOLERANCE.into(),
    );

    // Convert the MultiLineString back into a MultiPolygon
    let cleaned_multipolygon: MultiPolygon<f64> = MultiPolygon(
        cleaned_multipolygon
            .0
            .into_iter()
            .map(|line_string| Polygon::new(line_string, vec![]))
            .collect(),
    );
    cleaned_multipolygon
}
