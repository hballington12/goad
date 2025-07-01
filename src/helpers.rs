//! Utility functions for visualization and debugging support.
//!
//! This module provides helper functions primarily focused on interactive
//! visualization and debugging of geometric operations during development.
//! It includes coordinate transformations, rendering utilities, and simple
//! testing functions to support the development and validation process.
//!
//! The helper system provides:
//! - Interactive geometric visualization
//! - Coordinate system transformations for display
//! - Face and polygon rendering utilities
//! - Screen coordinate mapping
//! - Development and testing support functions
//!
//! # Visualization Features
//!
//! - Real-time rendering of faces and polygons
//! - Coordinate transformation for screen display
//! - Color-coded visualization for different geometric elements
//! - Interactive debugging support for geometric operations
//! - Simple utility functions for framework validation

use crate::geom::Face;
use geo_types::{Coord, Polygon};
use macroquad::prelude::*;

const SCALE: f32 = 25.0; // modify this depending on widget size
const OFFSET_X: f32 = 400.0;
const OFFSET_Y: f32 = 300.0;
/// Renders a geometric polygon for interactive visualization.
/// 
/// **Context**: Interactive debugging and educational visualization require
/// rendering geometric shapes to help understand clipping operations,
/// beam propagation, and geometric transformations.
/// 
/// **How it Works**: Converts geometric coordinates to screen coordinates
/// with scaling and offset, then draws line segments between consecutive
/// vertices using the graphics library.
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

/// Renders a Face structure for geometric visualization.
/// 
/// **Context**: Faces represent beam cross-sections and particle surface
/// elements that need visualization for debugging clipping operations,
/// geometric intersections, and beam propagation paths.
/// 
/// **How it Works**: Handles both simple and complex faces with holes,
/// extracting vertex coordinates and converting them to line strings
/// for rendering with specified color and thickness.
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

/// Converts geometric line strings to screen rendering.
/// 
/// **Context**: Geometric visualization requires coordinate transformation
/// from world space to screen space with appropriate scaling and offset
/// for readable display.
/// 
/// **How it Works**: Applies scale factor and offset transformation to
/// convert geometric coordinates to screen pixels, then renders line
/// segments between consecutive points to visualize the geometry.
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

/// Simple utility function for testing framework validation.
/// 
/// **Context**: Rust testing frameworks require simple functions to
/// validate that the testing infrastructure is working correctly.
/// 
/// **How it Works**: Performs trivial arithmetic to test basic
/// function calling and return value handling.
pub fn add_one(x: i32) -> i32 {
    x + 1
}
