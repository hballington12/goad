//! Consistent, outward face orientation for loaded meshes.
//!
//! Face normals follow the vertex winding, and the clipping uses them to decide
//! from which side a beam can meet a face. A mesh with mixed winding, or wound
//! inward as a whole, would give wrong normals. `orient_outward` makes the
//! winding consistent across shared edges and then turns each connected
//! component outward, so results do not depend on the winding of the input.

use std::collections::{HashMap, VecDeque};

use nalgebra::Point3;

/// Reorder the vertex indices of `faces` so that the winding is consistent across
/// shared edges and every connected component encloses a positive volume.
///
/// Returns the number of faces whose orientation was reversed.
pub fn orient_outward(vertices: &[Point3<f32>], faces: &mut [Vec<usize>]) -> usize {
    let mut edge_faces: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (fi, face) in faces.iter().enumerate() {
        for (a, b) in edges(face) {
            edge_faces.entry((a.min(b), a.max(b))).or_default().push(fi);
        }
    }

    let mut flipped = vec![false; faces.len()];
    let mut visited = vec![false; faces.len()];
    for start in 0..faces.len() {
        if visited[start] {
            continue;
        }
        let component = orient_component(start, faces, &edge_faces, &mut visited, &mut flipped);
        if signed_volume(vertices, faces, &component) < 0.0 {
            for &fi in &component {
                faces[fi].reverse();
                flipped[fi] = !flipped[fi];
            }
        }
    }
    flipped.iter().filter(|&&f| f).count()
}

/// Breadth-first walk over shared edges from `start`, reversing any neighbour that
/// traverses a shared edge in the same direction. Returns the faces of the component.
fn orient_component(
    start: usize,
    faces: &mut [Vec<usize>],
    edge_faces: &HashMap<(usize, usize), Vec<usize>>,
    visited: &mut [bool],
    flipped: &mut [bool],
) -> Vec<usize> {
    let mut component = vec![start];
    let mut queue = VecDeque::from([start]);
    visited[start] = true;
    while let Some(fi) = queue.pop_front() {
        let face_edges: Vec<_> = edges(&faces[fi]).collect();
        for (a, b) in face_edges {
            for &nj in &edge_faces[&(a.min(b), a.max(b))] {
                if visited[nj] {
                    continue;
                }
                // consistent winding: the neighbour traverses the shared edge as b -> a
                if edges(&faces[nj]).any(|e| e == (a, b)) {
                    faces[nj].reverse();
                    flipped[nj] = !flipped[nj];
                }
                visited[nj] = true;
                component.push(nj);
                queue.push_back(nj);
            }
        }
    }
    component
}

fn edges(face: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    face.iter()
        .zip(face.iter().cycle().skip(1))
        .map(|(&a, &b)| (a, b))
}

/// Signed volume enclosed by the faces of one component (fan triangulation).
fn signed_volume(vertices: &[Point3<f32>], faces: &[Vec<usize>], component: &[usize]) -> f64 {
    let mut volume = 0.0f64;
    for &fi in component {
        let f = &faces[fi];
        let p0 = vertices[f[0]].coords.cast::<f64>();
        for k in 1..f.len().saturating_sub(1) {
            let p1 = vertices[f[k]].coords.cast::<f64>();
            let p2 = vertices[f[k + 1]].coords.cast::<f64>();
            volume += p0.dot(&p1.cross(&p2));
        }
    }
    volume / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube() -> (Vec<Point3<f32>>, Vec<Vec<usize>>) {
        let v = (0..8)
            .map(|i| Point3::new((i & 1) as f32, ((i >> 1) & 1) as f32, ((i >> 2) & 1) as f32))
            .collect();
        // outward (counter-clockwise seen from outside)
        let f = vec![
            vec![0, 2, 3, 1],
            vec![4, 5, 7, 6],
            vec![0, 1, 5, 4],
            vec![2, 6, 7, 3],
            vec![0, 4, 6, 2],
            vec![1, 3, 7, 5],
        ];
        (v, f)
    }

    #[test]
    fn outward_cube_is_unchanged() {
        let (v, mut f) = cube();
        let before = f.clone();
        assert_eq!(orient_outward(&v, &mut f), 0);
        assert_eq!(f, before);
        assert!(signed_volume(&v, &f, &[0, 1, 2, 3, 4, 5]) > 0.0);
    }

    #[test]
    fn inward_cube_is_turned_outward() {
        let (v, mut f) = cube();
        let outward = f.clone();
        f.iter_mut().for_each(|face| face.reverse());
        assert_eq!(orient_outward(&v, &mut f), 6);
        assert!(signed_volume(&v, &f, &[0, 1, 2, 3, 4, 5]) > 0.0);
        for (a, b) in f.iter().zip(&outward) {
            assert!(edges(a).all(|e| edges(b).any(|g| g == e)));
        }
    }

    #[test]
    fn mixed_winding_is_made_consistent() {
        let (v, mut f) = cube();
        f[1].reverse();
        f[4].reverse();
        assert_eq!(orient_outward(&v, &mut f), 2);
        assert!(signed_volume(&v, &f, &[0, 1, 2, 3, 4, 5]) > 0.0);
    }
}
