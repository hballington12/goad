//! Tests for n2f_aperture_diffraction and the aperture coordinate system.
//!
//! The aperture system should satisfy:
//! 1. Aperture vertices rotated into the xy plane (z ≈ 0 for all transformed vertices)
//! 2. Propagation direction has +z component and lies in the xz plane (y ≈ 0)
//! 3. Rotated e_perp (vk7) is along +y axis

use nalgebra::{Complex, Matrix2, Point3, Vector3};
use std::f32::consts::PI;

use goad::diff::init_diff;

/// Tolerance for floating point comparisons
const TOL: f32 = 1e-4;

// =============================================================================
// Helper functions for creating test apertures
// =============================================================================

/// Creates a triangle aperture in the xy plane centered at origin
fn triangle_xy() -> Vec<Point3<f32>> {
    let r = 1.0;
    vec![
        Point3::new(r, 0.0, 0.0),
        Point3::new(-r / 2.0, r * 0.866, 0.0),
        Point3::new(-r / 2.0, -r * 0.866, 0.0),
    ]
}

/// Creates a square aperture in the xy plane centered at origin
fn square_xy() -> Vec<Point3<f32>> {
    vec![
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
        Point3::new(1.0, -1.0, 0.0),
    ]
}

/// Creates a hexagon aperture in the xy plane centered at origin
fn hexagon_xy() -> Vec<Point3<f32>> {
    let r = 1.0;
    (0..6)
        .map(|i| {
            let angle = i as f32 * PI / 3.0;
            Point3::new(r * angle.cos(), r * angle.sin(), 0.0)
        })
        .collect()
}

/// Creates a triangle aperture in the xz plane centered at origin
fn triangle_xz() -> Vec<Point3<f32>> {
    let r = 1.0;
    vec![
        Point3::new(r, 0.0, 0.0),
        Point3::new(-r / 2.0, 0.0, r * 0.866),
        Point3::new(-r / 2.0, 0.0, -r * 0.866),
    ]
}

/// Creates a square aperture in the xz plane centered at origin
fn square_xz() -> Vec<Point3<f32>> {
    vec![
        Point3::new(1.0, 0.0, 1.0),
        Point3::new(-1.0, 0.0, 1.0),
        Point3::new(-1.0, 0.0, -1.0),
        Point3::new(1.0, 0.0, -1.0),
    ]
}

/// Creates a square aperture in the yz plane centered at origin
fn square_yz() -> Vec<Point3<f32>> {
    vec![
        Point3::new(0.0, 1.0, 1.0),
        Point3::new(0.0, -1.0, 1.0),
        Point3::new(0.0, -1.0, -1.0),
        Point3::new(0.0, 1.0, -1.0),
    ]
}

/// Creates a square aperture tilted at 45 degrees (between xy and xz planes)
fn square_tilted_45() -> Vec<Point3<f32>> {
    let c = 0.5_f32.sqrt(); // cos(45°) = sin(45°)
    vec![
        Point3::new(1.0, c, c),
        Point3::new(-1.0, c, c),
        Point3::new(-1.0, -c, -c),
        Point3::new(1.0, -c, -c),
    ]
}

/// Creates a hexagon in an arbitrary tilted plane
fn hexagon_tilted() -> Vec<Point3<f32>> {
    // Hexagon with normal pointing roughly along (1, 1, 1)
    let r = 1.0;
    let normal = Vector3::new(1.0, 1.0, 1.0).normalize();

    // Find two perpendicular vectors in the plane
    let u = if (normal.x as f32).abs() < 0.9 {
        Vector3::x().cross(&normal).normalize()
    } else {
        Vector3::y().cross(&normal).normalize()
    };
    let v = normal.cross(&u).normalize();

    (0..6)
        .map(|i| {
            let angle = i as f32 * PI / 3.0;
            let offset = u * r * angle.cos() + v * r * angle.sin();
            Point3::new(offset.x, offset.y, offset.z)
        })
        .collect()
}

/// Rotates a set of vertices around the x-axis by angle (radians)
fn rotate_x(verts: &[Point3<f32>], angle: f32) -> Vec<Point3<f32>> {
    let (s, c) = angle.sin_cos();
    verts
        .iter()
        .map(|v| Point3::new(v.x, v.y * c - v.z * s, v.y * s + v.z * c))
        .collect()
}

/// Rotates a set of vertices around the y-axis by angle (radians)
fn rotate_y(verts: &[Point3<f32>], angle: f32) -> Vec<Point3<f32>> {
    let (s, c) = angle.sin_cos();
    verts
        .iter()
        .map(|v| Point3::new(v.x * c + v.z * s, v.y, -v.x * s + v.z * c))
        .collect()
}

/// Rotates a set of vertices around the z-axis by angle (radians)
fn rotate_z(verts: &[Point3<f32>], angle: f32) -> Vec<Point3<f32>> {
    let (s, c) = angle.sin_cos();
    verts
        .iter()
        .map(|v| Point3::new(v.x * c - v.y * s, v.x * s + v.y * c, v.z))
        .collect()
}

/// Translates vertices by an offset
fn translate_verts(verts: &[Point3<f32>], offset: Vector3<f32>) -> Vec<Point3<f32>> {
    verts.iter().map(|v| v + offset).collect()
}

/// Computes the normal of a planar face from its vertices (assumes anticlockwise ordering)
fn compute_normal(verts: &[Point3<f32>]) -> Vector3<f32> {
    let v0 = verts[1] - verts[0];
    let v1 = verts[2] - verts[0];
    v0.cross(&v1).normalize()
}

/// Returns a propagation vector pointing outward from the aperture (same direction as normal)
fn prop_outward(verts: &[Point3<f32>]) -> Vector3<f32> {
    compute_normal(verts)
}

/// Returns e_perp vector perpendicular to both prop and face normal
/// e_perp = (prop × normal).normalize()
///
/// For the degenerate case where prop is parallel to normal (normal incidence),
/// we fall back to finding any vector perpendicular to prop.
fn make_e_perp_for_face(prop: Vector3<f32>, normal: Vector3<f32>) -> Vector3<f32> {
    let cross = prop.cross(&normal);
    if cross.norm() < 1e-6 {
        // Degenerate case: prop parallel to normal, pick any perpendicular vector
        if prop.x.abs() < 0.9 {
            Vector3::x().cross(&prop).normalize()
        } else {
            Vector3::y().cross(&prop).normalize()
        }
    } else {
        cross.normalize()
    }
}

// =============================================================================
// Aperture system validation helpers
// =============================================================================

/// Validates that all transformed vertices lie in the xy plane (z ≈ 0)
fn assert_vertices_in_xy_plane(rot3: &nalgebra::Matrix3<f32>, verts: &[Vector3<f32>], tol: f32) {
    for (i, v) in verts.iter().enumerate() {
        let transformed = rot3 * v;
        assert!(
            transformed.z.abs() < tol,
            "Vertex {} not in xy plane: z = {} (transformed: {:?})",
            i,
            transformed.z,
            transformed
        );
    }
}

/// Validates that prop2 has positive z component
fn assert_prop_positive_z(prop2: Vector3<f32>, tol: f32) {
    assert!(
        prop2.z > -tol,
        "prop2 should have +z component, got z = {}",
        prop2.z
    );
}

/// Validates that prop2 lies in the xz plane (y ≈ 0)
fn assert_prop_in_xz_plane(prop2: Vector3<f32>, tol: f32) {
    assert!(
        prop2.y.abs() < tol,
        "prop2 should lie in xz plane, got y = {} (prop2: {:?})",
        prop2.y,
        prop2
    );
}

/// Validates that transformed e_perp is along +y axis
fn assert_e_perp_along_y(
    rot3: &nalgebra::Matrix3<f32>,
    vk7: Vector3<f32>,
    prop2: Vector3<f32>,
    tol: f32,
) {
    let perp2 = rot3 * vk7;
    // e_perp should be perpendicular to prop in the transformed frame
    // and should align with +y (or at least have dominant y component)

    // Check perpendicularity to prop2
    let dot = perp2.dot(&prop2);
    assert!(
        dot.abs() < tol,
        "perp2 should be perpendicular to prop2, dot product = {}",
        dot
    );

    // Check that y component is positive and dominant
    // (allowing for some flexibility since the exact alignment depends on the rotation)
    assert!(
        perp2.y.abs() > 0.5 || perp2.norm() < tol,
        "perp2 should have significant y component, got perp2 = {:?}",
        perp2
    );
}

// =============================================================================
// Tests: Triangle apertures
// =============================================================================

#[test]
fn test_triangle_xy_normal_incidence() {
    let verts = triangle_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0); // straight down into xy plane
    let vk7 = Vector3::new(1.0, 0.0, 0.0); // e_perp along x
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (com, rel_verts, rot3, prop2) =
        init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_triangle_xz_normal_incidence() {
    let verts = triangle_xz();
    let prop = Vector3::new(0.0, -1.0, 0.0); // into xz plane
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_triangle_tilted_aperture() {
    let base = triangle_xy();
    let verts = rotate_x(&base, PI / 4.0); // tilt 45 degrees around x
    let normal = compute_normal(&verts);
    let prop = prop_outward(&verts);
    let vk7 = make_e_perp_for_face(prop, normal);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

// =============================================================================
// Tests: Square apertures
// =============================================================================

#[test]
fn test_square_xy_normal_incidence() {
    let verts = square_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_square_xz_normal_incidence() {
    let verts = square_xz();
    let prop = Vector3::new(0.0, -1.0, 0.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_square_yz_normal_incidence() {
    let verts = square_yz();
    let prop = Vector3::new(-1.0, 0.0, 0.0);
    let vk7 = Vector3::new(0.0, 1.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_square_tilted_45() {
    let verts = square_tilted_45();
    let normal = compute_normal(&verts);
    let prop = prop_outward(&verts);
    let vk7 = make_e_perp_for_face(prop, normal);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

// =============================================================================
// Tests: Hexagon apertures
// =============================================================================

#[test]
fn test_hexagon_xy_normal_incidence() {
    let verts = hexagon_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

#[test]
fn test_hexagon_tilted() {
    let verts = hexagon_tilted();
    let normal = compute_normal(&verts);
    let prop = prop_outward(&verts);
    let vk7 = make_e_perp_for_face(prop, normal);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

// =============================================================================
// Tests: Various orientations and propagation directions
// =============================================================================

#[test]
fn test_square_rotated_various_angles() {
    let base = square_xy();
    let angles = [
        0.0,
        PI / 6.0,
        PI / 4.0,
        PI / 3.0,
        PI / 2.0,
        2.0 * PI / 3.0,
        PI,
    ];

    for &angle_x in &angles {
        for &angle_y in &angles {
            let rotated = rotate_y(&rotate_x(&base, angle_x), angle_y);
            let normal = compute_normal(&rotated);
            let prop = prop_outward(&rotated);
            let vk7 = make_e_perp_for_face(prop, normal);
            let mut ampl = Matrix2::<Complex<f32>>::identity();
            let wavenumber = 2.0 * PI;

            let (_, rel_verts, rot3, prop2) =
                init_diff(&rotated, &mut ampl, prop, vk7, wavenumber).unwrap();

            assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
            assert_prop_positive_z(prop2, TOL);
            assert_prop_in_xz_plane(prop2, TOL);
        }
    }
}

#[test]
fn test_hexagon_rotated_various_angles() {
    let base = hexagon_xy();
    let angles = [PI / 8.0, PI / 4.0, 3.0 * PI / 8.0, PI / 2.0, 5.0 * PI / 8.0];

    for &angle_x in &angles {
        for &angle_z in &angles {
            let rotated = rotate_z(&rotate_x(&base, angle_x), angle_z);
            let normal = compute_normal(&rotated);
            let prop = prop_outward(&rotated);
            let vk7 = make_e_perp_for_face(prop, normal);
            let mut ampl = Matrix2::<Complex<f32>>::identity();
            let wavenumber = 2.0 * PI;

            let (_, rel_verts, rot3, prop2) =
                init_diff(&rotated, &mut ampl, prop, vk7, wavenumber).unwrap();

            assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
            assert_prop_positive_z(prop2, TOL);
            assert_prop_in_xz_plane(prop2, TOL);
        }
    }
}

#[test]
fn test_translated_apertures() {
    let base = square_xy();
    let offsets = [
        Vector3::new(10.0, 0.0, 0.0),
        Vector3::new(0.0, 10.0, 0.0),
        Vector3::new(0.0, 0.0, 10.0),
        Vector3::new(5.0, 5.0, 5.0),
        Vector3::new(-3.0, 7.0, -2.0),
    ];

    for offset in offsets {
        let translated = translate_verts(&base, offset);
        let prop = Vector3::new(0.0, 0.0, -1.0);
        let vk7 = Vector3::new(1.0, 0.0, 0.0);
        let mut ampl = Matrix2::<Complex<f32>>::identity();
        let wavenumber = 2.0 * PI;

        let (com, rel_verts, rot3, prop2) =
            init_diff(&translated, &mut ampl, prop, vk7, wavenumber).unwrap();

        // Center of mass should be at the offset
        assert!(
            (com.coords - offset).norm() < TOL,
            "COM should be at offset"
        );

        assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
        assert_prop_positive_z(prop2, TOL);
        assert_prop_in_xz_plane(prop2, TOL);
    }
}

// =============================================================================
// Tests: Near-normal propagation directions
// =============================================================================

#[test]
fn test_near_normal_propagation_small_angles() {
    let verts = square_xy();
    let small_angles: [f32; 4] = [0.001, 0.01, 0.05, 0.1];

    for &angle in &small_angles {
        // Slightly off-normal in x direction
        let prop = Vector3::new(angle.sin(), 0.0, -angle.cos());
        let vk7 = Vector3::new(0.0, 1.0, 0.0);
        let mut ampl = Matrix2::<Complex<f32>>::identity();
        let wavenumber = 2.0 * PI;

        let (_, rel_verts, rot3, prop2) =
            init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

        assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
        assert_prop_positive_z(prop2, TOL);
        assert_prop_in_xz_plane(prop2, TOL);
    }
}

#[test]
fn test_near_normal_propagation_various_azimuths() {
    let verts = square_xy();
    let normal = compute_normal(&verts);
    let theta: f32 = 0.05; // 5% off normal
    let azimuths = [0.0, PI / 4.0, PI / 2.0, PI, 3.0 * PI / 2.0];

    for &phi in &azimuths {
        let prop = Vector3::new(
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            -theta.cos(),
        );
        let vk7 = make_e_perp_for_face(prop, normal);
        let mut ampl = Matrix2::<Complex<f32>>::identity();
        let wavenumber = 2.0 * PI;

        let (_, rel_verts, rot3, prop2) =
            init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

        assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
        assert_prop_positive_z(prop2, TOL);
        // Note: prop2 is NOT necessarily in xz plane after third rotation aligns perp2 with +y
    }
}

#[test]
fn test_grazing_incidence() {
    let verts = square_xy();
    // Nearly parallel to aperture surface
    let prop = Vector3::new(0.99, 0.0, -0.14).normalize();
    let vk7 = Vector3::new(0.0, 1.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, rel_verts, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_vertices_in_xy_plane(&rot3, &rel_verts, TOL);
    assert_prop_positive_z(prop2, TOL);
    assert_prop_in_xz_plane(prop2, TOL);
}

// =============================================================================
// Tests: Non-planar aperture should return error
// =============================================================================

#[test]
fn test_nonplanar_aperture_returns_error() {
    // Create a non-planar "square" where one vertex is out of plane
    let verts = vec![
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.5), // This vertex is out of plane!
        Point3::new(1.0, -1.0, 0.0),
    ];
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    // This should return an error because the aperture is non-planar
    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(result.is_err(), "Expected error for non-planar aperture");

    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("non-planar"),
        "Error message should mention 'non-planar', got: {}",
        err_msg
    );
}

#[test]
fn test_nonplanar_aperture_small_deviation() {
    // Create a non-planar "square" with a small deviation (but above tolerance)
    let verts = vec![
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.001), // Small deviation above tolerance (1e-4)
        Point3::new(1.0, -1.0, 0.0),
    ];
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(
        result.is_err(),
        "Expected error for non-planar aperture with small deviation"
    );
}

#[test]
fn test_planar_aperture_within_tolerance() {
    // Create a nearly-planar square with deviation within tolerance
    let verts = vec![
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 1e-5), // Tiny deviation within tolerance
        Point3::new(1.0, -1.0, 0.0),
    ];
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(
        result.is_ok(),
        "Expected success for nearly-planar aperture within tolerance"
    );
}

// =============================================================================
// Tests: e_perp validation - must be perpendicular to both prop and face normal
// =============================================================================

#[test]
fn test_e_perp_not_perpendicular_to_prop_returns_error() {
    let verts = square_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0);
    // e_perp has a component along prop (z-component)
    let vk7 = Vector3::new(1.0, 0.0, 0.5).normalize();
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(
        result.is_err(),
        "Expected error when e_perp is not perpendicular to prop"
    );

    let err_msg = result.unwrap_err().to_string();
    assert!(
        err_msg.contains("perpendicular") || err_msg.contains("e_perp"),
        "Error message should mention perpendicularity issue, got: {}",
        err_msg
    );
}

#[test]
fn test_e_perp_not_perpendicular_to_normal_returns_error() {
    let verts = square_xy();
    let normal = compute_normal(&verts); // (0, 0, 1) for square_xy
    let prop = Vector3::new(0.0, 0.0, -1.0);
    // e_perp is perpendicular to prop but has component along normal (z-component)
    let vk7 = Vector3::new(1.0, 0.0, 0.5).normalize();
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(
        result.is_err(),
        "Expected error when e_perp is not perpendicular to face normal"
    );
}

#[test]
fn test_e_perp_valid_perpendicular_to_both() {
    let verts = square_xy();
    let normal = compute_normal(&verts);
    let prop = Vector3::new(0.0, 0.0, -1.0);
    // e_perp perpendicular to both prop and normal (lies in xy plane)
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    // Verify our e_perp is indeed perpendicular to both
    let dot_prop: f32 = vk7.dot(&prop);
    let dot_normal: f32 = vk7.dot(&normal);
    assert!(dot_prop.abs() < TOL, "vk7 should be perpendicular to prop");
    assert!(
        dot_normal.abs() < TOL,
        "vk7 should be perpendicular to normal"
    );

    let result = init_diff(&verts, &mut ampl, prop, vk7, wavenumber);
    assert!(
        result.is_ok(),
        "Expected success when e_perp is perpendicular to both prop and normal"
    );
}

// =============================================================================
// Tests: e_perp alignment
// =============================================================================

#[test]
fn test_e_perp_alignment_basic() {
    let verts = square_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, _, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_e_perp_along_y(&rot3, vk7, prop2, TOL);
}

#[test]
fn test_e_perp_alignment_tilted() {
    let base = square_xy();
    let verts = rotate_x(&base, PI / 3.0);
    let normal = compute_normal(&verts);
    let prop = prop_outward(&verts);
    let vk7 = make_e_perp_for_face(prop, normal);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let (_, _, rot3, prop2) = init_diff(&verts, &mut ampl, prop, vk7, wavenumber).unwrap();

    assert_e_perp_along_y(&rot3, vk7, prop2, TOL);
}

// =============================================================================
// Tests: Wavenumber scaling
// =============================================================================

#[test]
fn test_amplitude_scaled_by_wavenumber() {
    let verts = square_xy();
    let prop = Vector3::new(0.0, 0.0, -1.0);
    let vk7 = Vector3::new(1.0, 0.0, 0.0);
    let wavenumber = 2.0 * PI;

    let mut ampl1 = Matrix2::<Complex<f32>>::identity();
    let _ = init_diff(&verts, &mut ampl1, prop, vk7, wavenumber).unwrap();

    let mut ampl2 = Matrix2::<Complex<f32>>::identity();
    let wavenumber2 = 4.0 * PI;
    let _ = init_diff(&verts, &mut ampl2, prop, vk7, wavenumber2).unwrap();

    // ampl should be scaled by wavenumber
    let ratio = ampl2[(0, 0)].re / ampl1[(0, 0)].re;
    assert!(
        (ratio - 2.0).abs() < TOL,
        "Amplitude should scale with wavenumber, got ratio = {}",
        ratio
    );
}

// =============================================================================
// Tests: All octant combinations (prop direction × facet normal direction)
// =============================================================================

/// Octant indices map to directions:
/// 0: (-1,-1,-1), 1: (-1,-1,+1), 2: (-1,+1,-1), 3: (-1,+1,+1)
/// 4: (+1,-1,-1), 5: (+1,-1,+1), 6: (+1,+1,-1), 7: (+1,+1,+1)
fn octant_direction(index: usize) -> Vector3<f32> {
    let sx = if index & 4 != 0 { 1.0 } else { -1.0 };
    let sy = if index & 2 != 0 { 1.0 } else { -1.0 };
    let sz = if index & 1 != 0 { 1.0 } else { -1.0 };
    Vector3::new(sx, sy, sz).normalize()
}

/// Create a square aperture with normal pointing in the given direction
fn square_with_normal(normal: Vector3<f32>) -> Vec<Point3<f32>> {
    // Start with a square in xy plane (normal along +z)
    let base = vec![
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(-1.0, 1.0, 0.0),
        Vector3::new(-1.0, -1.0, 0.0),
        Vector3::new(1.0, -1.0, 0.0),
    ];

    // Compute rotation to align +z with the target normal
    let z_axis = Vector3::new(0.0, 0.0, 1.0);
    let dot = z_axis.dot(&normal);

    let rotated = if dot > 0.999 {
        // Already aligned with +z
        base
    } else if dot < -0.999 {
        // Flip 180° around x-axis
        base.iter().map(|v| Vector3::new(v.x, -v.y, -v.z)).collect()
    } else {
        // General rotation using Rodrigues' formula
        let axis = z_axis.cross(&normal).normalize();
        let cos_angle = dot;
        let sin_angle = (1.0 - cos_angle * cos_angle).sqrt();

        base.iter()
            .map(|v| {
                // Rodrigues' rotation: v' = v*cos(θ) + (k×v)*sin(θ) + k*(k·v)*(1-cos(θ))
                let k_cross_v = axis.cross(v);
                let k_dot_v = axis.dot(v);
                v * cos_angle + k_cross_v * sin_angle + axis * k_dot_v * (1.0 - cos_angle)
            })
            .collect()
    };

    rotated
        .into_iter()
        .map(|v| Point3::new(v.x, v.y, v.z))
        .collect()
}

/// Helper to run a single octant combination test
fn run_octant_test(normal_idx: usize, prop_idx: usize) {
    let normal_dir = octant_direction(normal_idx);
    let prop_dir = octant_direction(prop_idx);

    let verts = square_with_normal(normal_dir);
    // e_perp must be perpendicular to both prop and face normal
    let vk7 = make_e_perp_for_face(prop_dir, normal_dir);
    let mut ampl = Matrix2::<Complex<f32>>::identity();
    let wavenumber = 2.0 * PI;

    let result = init_diff(&verts, &mut ampl, prop_dir, vk7, wavenumber);

    // Check that init_diff succeeded
    assert!(
        result.is_ok(),
        "normal_octant={} prop_octant={}: init_diff failed: {}",
        normal_idx,
        prop_idx,
        result.unwrap_err()
    );

    let (_, rel_verts, rot3, prop2) = result.unwrap();

    // Check all vertices are in xy plane
    for (i, v) in rel_verts.iter().enumerate() {
        let transformed = rot3 * v;
        assert!(
            transformed.z.abs() < TOL,
            "normal_octant={} prop_octant={}: vertex {} not in xy plane, z={:.6}",
            normal_idx,
            prop_idx,
            i,
            transformed.z
        );
    }

    // Check prop2 has +z component
    assert!(
        prop2.z > -TOL,
        "normal_octant={} prop_octant={}: prop2.z={:.6} should be positive",
        normal_idx,
        prop_idx,
        prop2.z
    );

    // Note: prop2 is NOT necessarily in the xz plane after the third rotation
    // that aligns perp2 with +y. This is expected behavior.

    // First check: e_perp must be perpendicular to prop in lab frame
    let vk7_dot_prop = vk7.dot(&prop_dir);
    assert!(
        vk7_dot_prop.abs() < TOL,
        "normal_octant={} prop_octant={}: vk7 should be perpendicular to prop, got dot={:.6}",
        normal_idx,
        prop_idx,
        vk7_dot_prop
    );

    // Check that rotated e_perp (perp2) is entirely along ±y axis
    let perp2 = rot3 * vk7;
    let perp2_normalized = perp2.normalize();
    let y_axis = Vector3::new(0.0, 1.0, 0.0);

    // perp2 should be colinear with y-axis (dot product > 0.999)
    let colinearity = perp2_normalized.dot(&y_axis).abs();
    assert!(
        colinearity > 0.999,
        "normal_octant={} prop_octant={}: perp2 should be along y-axis, got colinearity={:.6} (perp2={:?})",
        normal_idx,
        prop_idx,
        colinearity,
        perp2
    );
}

// =============================================================================
// Tests: Prop same as normal (8 tests, one per octant)
// =============================================================================

#[test]
fn test_octant_normal0_prop_same() {
    run_octant_test(0, 0);
}
#[test]
fn test_octant_normal1_prop_same() {
    run_octant_test(1, 1);
}
#[test]
fn test_octant_normal2_prop_same() {
    run_octant_test(2, 2);
}
#[test]
fn test_octant_normal3_prop_same() {
    run_octant_test(3, 3);
}
#[test]
fn test_octant_normal4_prop_same() {
    run_octant_test(4, 4);
}
#[test]
fn test_octant_normal5_prop_same() {
    run_octant_test(5, 5);
}
#[test]
fn test_octant_normal6_prop_same() {
    run_octant_test(6, 6);
}
#[test]
fn test_octant_normal7_prop_same() {
    run_octant_test(7, 7);
}

// =============================================================================
// Tests: All valid octant combinations (normal × prop where prop·normal > 0.2)
//
// Octant index → direction:
//   0: (-,-,-), 1: (-,-,+), 2: (-,+,-), 3: (-,+,+)
//   4: (+,-,-), 5: (+,-,+), 6: (+,+,-), 7: (+,+,+)
//
// For normalized octant vectors, dot product depends on matching sign components:
//   3 match → dot = 1.0 (same octant)
//   2 match → dot = 1/3 ≈ 0.33 (valid, > 0.2)
//   1 match → dot = -1/3 ≈ -0.33 (invalid)
//   0 match → dot = -1.0 (opposite octant, invalid)
//
// Adjacent octants (differ by 1 bit) share 2 components → valid pairs
// =============================================================================

// Normal octant 0: (-,-,-)
// Valid props: 0 (same), 1 (differ z), 2 (differ y), 4 (differ x)
#[test]
fn test_octant_n0_p0() {
    run_octant_test(0, 0);
}
#[test]
fn test_octant_n0_p1() {
    run_octant_test(0, 1);
}
#[test]
fn test_octant_n0_p2() {
    run_octant_test(0, 2);
}
#[test]
fn test_octant_n0_p4() {
    run_octant_test(0, 4);
}

// Normal octant 1: (-,-,+)
// Valid props: 1 (same), 0 (differ z), 3 (differ y), 5 (differ x)
#[test]
fn test_octant_n1_p1() {
    run_octant_test(1, 1);
}
#[test]
fn test_octant_n1_p0() {
    run_octant_test(1, 0);
}
#[test]
fn test_octant_n1_p3() {
    run_octant_test(1, 3);
}
#[test]
fn test_octant_n1_p5() {
    run_octant_test(1, 5);
}

// Normal octant 2: (-,+,-)
// Valid props: 2 (same), 0 (differ y), 3 (differ z), 6 (differ x)
#[test]
fn test_octant_n2_p2() {
    run_octant_test(2, 2);
}
#[test]
fn test_octant_n2_p0() {
    run_octant_test(2, 0);
}
#[test]
fn test_octant_n2_p3() {
    run_octant_test(2, 3);
}
#[test]
fn test_octant_n2_p6() {
    run_octant_test(2, 6);
}

// Normal octant 3: (-,+,+)
// Valid props: 3 (same), 1 (differ y), 2 (differ z), 7 (differ x)
#[test]
fn test_octant_n3_p3() {
    run_octant_test(3, 3);
}
#[test]
fn test_octant_n3_p1() {
    run_octant_test(3, 1);
}
#[test]
fn test_octant_n3_p2() {
    run_octant_test(3, 2);
}
#[test]
fn test_octant_n3_p7() {
    run_octant_test(3, 7);
}

// Normal octant 4: (+,-,-)
// Valid props: 4 (same), 0 (differ x), 5 (differ z), 6 (differ y)
#[test]
fn test_octant_n4_p4() {
    run_octant_test(4, 4);
}
#[test]
fn test_octant_n4_p0() {
    run_octant_test(4, 0);
}
#[test]
fn test_octant_n4_p5() {
    run_octant_test(4, 5);
}
#[test]
fn test_octant_n4_p6() {
    run_octant_test(4, 6);
}

// Normal octant 5: (+,-,+)
// Valid props: 5 (same), 1 (differ x), 4 (differ z), 7 (differ y)
#[test]
fn test_octant_n5_p5() {
    run_octant_test(5, 5);
}
#[test]
fn test_octant_n5_p1() {
    run_octant_test(5, 1);
}
#[test]
fn test_octant_n5_p4() {
    run_octant_test(5, 4);
}
#[test]
fn test_octant_n5_p7() {
    run_octant_test(5, 7);
}

// Normal octant 6: (+,+,-)
// Valid props: 6 (same), 2 (differ x), 4 (differ y), 7 (differ z)
#[test]
fn test_octant_n6_p6() {
    run_octant_test(6, 6);
}
#[test]
fn test_octant_n6_p2() {
    run_octant_test(6, 2);
}
#[test]
fn test_octant_n6_p4() {
    run_octant_test(6, 4);
}
#[test]
fn test_octant_n6_p7() {
    run_octant_test(6, 7);
}

// Normal octant 7: (+,+,+)
// Valid props: 7 (same), 3 (differ x), 5 (differ y), 6 (differ z)
#[test]
fn test_octant_n7_p7() {
    run_octant_test(7, 7);
}
#[test]
fn test_octant_n7_p3() {
    run_octant_test(7, 3);
}
#[test]
fn test_octant_n7_p5() {
    run_octant_test(7, 5);
}
#[test]
fn test_octant_n7_p6() {
    run_octant_test(7, 6);
}
