use crate::geom::Face;
use geo_types::{Coord, Polygon};
use macroquad::prelude::*;

/// Draws the polygons from a MultiPolygon result onto the screen.
///
/// # Arguments
/// * `multi_polygon` - A reference to a `MultiPolygon` containing the polygons to draw.
pub fn draw_multipolygon(polygon: &Polygon<f32>, color: Color) {
    let scale = 25.0; // modify this depending on widget size
    let offset_x = 400.0;
    let offset_y = 300.0;
    // Extract the exterior LineString
    let points = &polygon.exterior().0;

    // Convert the points into macroquad-compatible coordinates
    let mut screen_points: Vec<(f32, f32)> = Vec::new();
    for coord in points {
        let screen_x = -coord.x as f32 * scale + offset_x; // Scale and center
        let screen_y = coord.y as f32 * scale + offset_y; // Scale and center
        screen_points.push((screen_x, screen_y));
    }

    // Draw the polygon by connecting the points
    for i in 0..screen_points.len() {
        let (x1, y1) = screen_points[i];
        let (x2, y2) = screen_points[(i + 1) % screen_points.len()]; // Wrap around
        draw_line(x1, y1, x2, y2, 2.0, color);
    }
}

/// Draws a polygon from a Face onto the screen.
///
/// # Arguments
/// * `face` - A reference to a `Face` containing the polygon to draw.
pub fn draw_face(face: &Face, color: Color) {
    let scale = 25.0; // modify this depending on widget size
    let offset_x = 400.0;
    let offset_y = 300.0;

    // Extract the exterior LineString
    let mut points = Vec::new();
    match face {
        Face::Simple(data) => {
            for vertex in &data.exterior {
                points.push(Coord {
                    x: vertex.x,
                    y: vertex.y,
                });
            }
        }
        Face::Complex { data, interiors } => {
            for vertex in &data.exterior {
                points.push(Coord {
                    x: vertex.x,
                    y: vertex.y,
                });
            }

            for interior in interiors {
                for vertex in interior {
                    points.push(Coord {
                        x: vertex.x,
                        y: vertex.y,
                    });
                }
            }
        }
    }

    // Convert the points into macroquad-compatible coordinates
    let mut screen_points: Vec<(f32, f32)> = Vec::new();
    for coord in points {
        let screen_x = -coord.x as f32 * scale + offset_x; // Scale and center
        let screen_y = coord.y as f32 * scale + offset_y; // Scale and center
        screen_points.push((screen_x, screen_y));
    }

    // Draw the polygon by connecting the points
    for i in 0..screen_points.len() {
        let (x1, y1) = screen_points[i];
        let (x2, y2) = screen_points[(i + 1) % screen_points.len()]; // Wrap around
        draw_line(x1, y1, x2, y2, 2.0, color);
    }
}
