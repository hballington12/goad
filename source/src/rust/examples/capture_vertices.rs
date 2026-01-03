use goad::clip::Clipping;
use goad::geom::{Geom, Face};
use nalgebra::{Point3, Vector3};

fn main() {
    println!("=== projection-debug (hex.obj, face 4, projection (0,0,-1)) ===");
    test_projection_debug();
    
    println!("\n=== projection1 (concave1.obj, face 4, projection (-0.3,0,-1)) ===");
    test_projection1();
    
    println!("\n=== projection2 (cube_inside_ico.obj) ===");
    test_projection2();
    
    println!("\n=== projection_multi (multiple.obj, shape[0].face[5], projection (-1,0,0)) ===");
    test_projection_multi();
    
    println!("\n=== clip_test (clip_test.obj) ===");
    test_clip_test();
    
    println!("\n=== remainder (multiple.obj, custom rect, projection (0,0,-1)) ===");
    test_remainder();
}

fn test_projection_debug() {
    let geoms = Geom::load("./examples/data/hex.obj").unwrap();
    let mut geom = geoms[0].clone();
    let mut clip = geom.shapes[0].faces.remove(4);
    let projection = Vector3::new(0.0, 0.0, -1.0);
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn test_projection1() {
    let geoms = Geom::load("./examples/data/concave1.obj").unwrap();
    let mut geom = geoms[0].clone();
    let mut clip = geom.shapes[0].faces.remove(4);
    let projection = Vector3::new(-0.3, 0.0, -1.0);
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn test_projection2() {
    let geoms = Geom::load("./examples/data/cube_inside_ico.obj").unwrap();
    let mut geom = geoms[0].clone();
    let mut clip = geom.shapes[0].faces.remove(5);
    let projection = Vector3::new(-0.2, 0.0, -1.0);
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn test_projection_multi() {
    let geoms = Geom::load("./examples/data/multiple.obj").unwrap();
    let mut geom = geoms[0].clone();
    let mut clip = geom.shapes[0].faces.remove(5);
    let projection = Vector3::new(-1.0, 0.0, 0.0);
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn test_clip_test() {
    let geoms = Geom::load("./examples/data/clip_test.obj").unwrap();
    let mut geom = geoms[0].clone();
    let mut clip = geom.shapes[1].faces.remove(1);
    let projection = Vector3::new(1.0, 1.0, 0.0);
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn test_remainder() {
    let geoms = Geom::load("./examples/data/multiple.obj").unwrap();
    let mut geom = geoms[0].clone();
    let projection = Vector3::new(0.0, 0.0, -1.0);
    
    let mut clip_vertices = vec![
        Point3::new(-19.0, 3.0, 10.0),
        Point3::new(-19.0, -3.0, 10.0),
        Point3::new(10.0, -3.0, 10.0),
        Point3::new(10.0, 3.0, 10.0),
    ];
    clip_vertices.reverse();
    let mut clip = Face::new_simple(clip_vertices, None, None).unwrap();
    
    let mut clipping = Clipping::new(&mut geom, &mut clip, &projection);
    let _ = clipping.clip(0.01);
    
    print_vertices(&clipping);
}

fn print_vertices(clipping: &Clipping) {
    for (i, intsn) in clipping.intersections.iter().enumerate() {
        let data = intsn.data();
        println!("intersection[{}]:", i);
        println!("  midpoint: ({:.6}, {:.6}, {:.6})", data.midpoint.x, data.midpoint.y, data.midpoint.z);
        println!("  normal: ({:.6}, {:.6}, {:.6})", data.normal.x, data.normal.y, data.normal.z);
        println!("  num_vertices: {}", data.num_vertices);
        println!("  vertices:");
        for (j, v) in data.exterior.iter().enumerate() {
            println!("    [{}]: ({:.6}, {:.6}, {:.6})", j, v.x, v.y, v.z);
        }
    }
    for (i, rem) in clipping.remaining.iter().enumerate() {
        let data = rem.data();
        println!("remaining[{}]:", i);
        println!("  midpoint: ({:.6}, {:.6}, {:.6})", data.midpoint.x, data.midpoint.y, data.midpoint.z);
        println!("  num_vertices: {}", data.num_vertices);
        println!("  vertices:");
        for (j, v) in data.exterior.iter().enumerate() {
            println!("    [{}]: ({:.6}, {:.6}, {:.6})", j, v.x, v.y, v.z);
        }
    }
}
