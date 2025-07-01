//! Geometric containment relationships and spatial acceleration structures.
//!
//! This module provides tools for managing hierarchical containment relationships
//! between geometric shapes in complex particle structures. It supports nested
//! geometries like coated particles, inclusions, and multi-material systems
//! where understanding the containment hierarchy is essential for proper
//! electromagnetic field propagation.
//!
//! The containment system provides:
//! - Parent-child relationship tracking
//! - Efficient containment queries
//! - Axis-aligned bounding boxes for spatial acceleration
//! - Support for arbitrary nesting levels
//! - Integration with material property assignment
//!
//! # Key Components
//!
//! - [`ContainmentGraph`]: Hierarchical shape relationships
//! - [`AABB`]: Axis-aligned bounding boxes for spatial queries
//! - Containment validation and consistency checking
//! - Material index resolution for nested structures

use nalgebra::Point3;

/// Hierarchical containment relationships between geometric shapes.
/// 
/// **Context**: Complex particle geometries often involve nested structures where
/// shapes are contained within other shapes (e.g., inclusions within a host medium).
/// These containment relationships affect refractive index transitions and beam
/// propagation paths, requiring explicit tracking.
/// 
/// **How it Works**: Maintains a parent-child relationship graph where each shape
/// can have at most one parent shape that contains it. This enables efficient
/// queries for determining the surrounding medium when beams exit a shape.
#[derive(Debug, Clone, PartialEq)]
pub struct ContainmentGraph {
    parent: Vec<Option<usize>>, // Maps each shape index to its containing shape
}

impl ContainmentGraph {
    /// Creates a new containment graph for tracking shape relationships.
    /// 
    /// **Context**: Containment graphs must be initialized with the total number
    /// of shapes in the geometry to allocate proper storage for relationships.
    /// 
    /// **How it Works**: Allocates a vector of optional parent indices, initially
    /// all None to indicate no containment relationships.
    pub fn new(num_shapes: usize) -> Self {
        Self {
            parent: vec![None; num_shapes], // Initially, no shapes are contained in others
        }
    }

    /// Establishes a containment relationship between shapes.
    /// 
    /// **Context**: Containment relationships must be explicitly specified based
    /// on geometric analysis. This allows modeling of complex nested structures
    /// like coated particles or embedded inclusions.
    /// 
    /// **How it Works**: Records that the child shape is contained within the
    /// parent shape, with bounds checking to prevent invalid relationships.
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

    /// Queries the containing shape for a given shape.
    /// 
    /// **Context**: When beams exit a shape, the simulation needs to determine
    /// the refractive index of the surrounding medium, which depends on the
    /// containing shape if any.
    /// 
    /// **How it Works**: Returns the parent shape index if one exists, or None
    /// if the shape is at the top level of the hierarchy.
    pub fn get_parent(&self, shape: usize) -> Option<usize> {
        self.parent[shape]
    }
}

/// Axis-aligned bounding box for spatial acceleration structures.
/// 
/// **Context**: Testing containment relationships between arbitrary 3D shapes
/// is computationally expensive. Axis-aligned bounding boxes provide a fast
/// preliminary test - if bounding boxes don't overlap, the shapes cannot
/// have a containment relationship.
/// 
/// **How it Works**: Stores minimum and maximum coordinates along each axis,
/// defining the smallest box that completely contains a shape. Used for
/// rapid spatial queries and containment pre-filtering.
#[derive(Debug, Clone, PartialEq)]
pub struct AABB {
    pub min: Point3<f32>,
    pub max: Point3<f32>,
}
