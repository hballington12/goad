use crate::geom::Face;
use geo_types::{Coord, Polygon};
#[cfg(feature = "macroquad")]
use macroquad::prelude::*;

const SCALE: f32 = 25.0; // modify this depending on widget size
const OFFSET_X: f32 = 400.0;
const OFFSET_Y: f32 = 300.0;
/// Draws the polygons from a MultiPolygon result onto the screen.
///
/// # Arguments
/// * `multi_polygon` - A reference to a `MultiPolygon` containing the polygons to draw.
#[cfg(feature = "macroquad")]
pub fn draw_multipolygon(polygon: &Polygon<f32>, color: Color) {
    // Extract the exterior LineString
    let points = &polygon.exterior().0;

    // Convert the points into macroquad-compatible coordinates
    let mut screen_points = Vec::new();
    for coord in points {
        let screen_x = -coord.x as f32 * SCALE + OFFSET_X; // Scale and center
        let screen_y = coord.y as f32 * SCALE + OFFSET_Y; // Scale and center
        screen_points.push((screen_x as f32, screen_y as f32));
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
#[cfg(feature = "macroquad")]
pub fn draw_face(face: &Face, color: Color, thickness: f32) {
    // Extract the exterior LineString
    let mut line_strings = Vec::new();
    match face {
        Face::Simple(data) => {
            let mut line_string = Vec::new();
            for vertex in &data.exterior {
                line_string.push(Coord {
                    x: vertex.x,
                    y: vertex.y,
                });
            }
            line_strings.push(line_string);
        }
        Face::Complex { data, interiors } => {
            let mut line_string = Vec::new();
            for vertex in &data.exterior {
                line_string.push(Coord {
                    x: vertex.x,
                    y: vertex.y,
                });
            }
            line_strings.push(line_string);

            for interior in interiors {
                let mut line_string = Vec::new();
                for vertex in interior {
                    line_string.push(Coord {
                        x: vertex.x,
                        y: vertex.y,
                    });
                }
                line_strings.push(line_string);
            }
        }
    }

    lines_to_screen(line_strings, color, thickness);
}

#[cfg(feature = "macroquad")]
pub fn lines_to_screen(line_strings: Vec<Vec<Coord<f32>>>, color: Color, thickness: f32) {
    // Convert the points into macroquad-compatible coordinates
    for points in line_strings {
        let mut screen_points = Vec::new();
        for coord in points {
            let screen_x = -coord.x as f32 * SCALE + OFFSET_X; // Scale and center
            let screen_y = coord.y as f32 * SCALE + OFFSET_Y; // Scale and center
            screen_points.push((screen_x as f32, screen_y as f32));
        }

        // Draw the polygon by connecting the points
        for i in 0..screen_points.len() {
            let (x1, y1) = screen_points[i];
            let (x2, y2) = screen_points[(i + 1) % screen_points.len()]; // Wrap around
            draw_line(x1, y1, x2, y2, thickness as f32, color);
        }
    }
}

pub fn add_one(x: i32) -> i32 {
    x + 1
}
