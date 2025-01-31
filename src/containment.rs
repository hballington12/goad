use nalgebra::Point3;

#[derive(Debug, Clone, PartialEq)]
pub struct ContainmentGraph {
    parent: Vec<Option<usize>>, // Maps each shape index to its containing shape
}

impl ContainmentGraph {
    /// Creates a new containment graph with a given number of shapes.
    pub fn new(num_shapes: usize) -> Self {
        Self {
            parent: vec![None; num_shapes], // Initially, no shapes are contained in others
        }
    }

    /// Sets the parent of a given shape.
    pub fn set_parent(&mut self, child: usize, parent: usize) {
        assert!(
            child < self.parent.len(),
            "child id is {}, but the containment graph only has space for {} shapes",
            child,
            self.parent.len()
        );
        assert!(
            parent < self.parent.len(),
            "parent id is {}, but the containment graph only has space for {} shapes",
            parent,
            self.parent.len()
        );
        self.parent[child] = Some(parent);
    }

    /// Returns the index of the shape that contains the given shape, if any.
    pub fn get_parent(&self, shape: usize) -> Option<usize> {
        self.parent[shape]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AABB {
    pub min: Point3<f32>,
    pub max: Point3<f32>,
}
