//! Recording layer for inspecting the near-field beam tree.
//!
//! A `Recording` is the captured output of a single-threaded
//! `Problem::solve_near_with_recording` (or `Problem::run_with_recording`)
//! pass. Each propagation event — one input beam consumed, zero or more
//! output beams produced — is stored as a `BeamEvent`. Truncated beams
//! (below thresholds, over recursion/TIR limits, or that produced no
//! output) are NOT recorded; their fate is implicit in the producing
//! parent's `RecordedOutput` having no follow-up event.
//!
//! `Recording::view()` returns a builder-style `BeamView` that filters
//! events by recursion depth, TIR count, variant, power, etc., and
//! exposes iterators over the filtered set. Materialise into owned data
//! via the `collect_*` methods when handing across an FFI boundary.

use nalgebra::{Point3, Vector3};

use crate::beam::{Beam, BeamId, BeamVariant, OutputKind};
use crate::field::Ampl;
use crate::geom::{Face, Plane};

pub type EventId = usize;

/// One output beam emitted by a propagation event, tagged with the queue it
/// would be dispatched to in the production path.
#[derive(Debug, Clone)]
pub struct RecordedOutput {
    pub beam: Beam,
    pub kind: OutputKind,
}

/// A single propagation event: an input beam that was consumed and the
/// outputs it produced. Outputs are in the order emitted by
/// `Beam::propagate` (intersections first, then remainders).
#[derive(Debug, Clone)]
pub struct BeamEvent {
    pub id: EventId,
    pub input: Beam,
    pub outputs: Vec<RecordedOutput>,
}

impl BeamEvent {
    /// Tree depth from the initial illumination event. Derived from the
    /// input beam's counters — every propagation event increments exactly
    /// one of `rec_count` / `tir_count`, so their sum is the path length
    /// from the initial event.
    pub fn depth(&self) -> i32 {
        self.input.rec_count + self.input.tir_count
    }

    /// Geometric decomposition of the input beam's column at this event.
    /// One entry per output. Each entry holds both:
    ///
    /// - `back`: a polygon on the input face's plane — for OutGoing this
    ///   is the output verbatim (already coplanar with the input); for
    ///   Internal / ExternalDiff outputs the output face is back-projected
    ///   along the input beam's propagation direction onto the input
    ///   plane.
    /// - `front`: the output's actual face vertices (the polygon on the
    ///   downstream plane).
    ///
    /// Use `back` alone for the 2D-on-plane subdivision view. Combine
    /// `back` + `front` + `DecomposedRegion::side_quads` for the full 3D
    /// prism.
    pub fn decompose(&self) -> Vec<DecomposedRegion> {
        let input_plane = self.input.face.plane();
        let dir = self.input.field.prop();

        self.outputs
            .iter()
            .enumerate()
            .map(|(i, ro)| {
                let front = polygon_rings_of(&ro.beam.face);
                let back = match ro.kind {
                    OutputKind::OutGoing => front.clone(),
                    OutputKind::Internal | OutputKind::ExternalDiff => PolygonRings {
                        exterior: back_project_ring(&front.exterior, &input_plane, dir),
                        interiors: front
                            .interiors
                            .iter()
                            .map(|ring| back_project_ring(ring, &input_plane, dir))
                            .collect(),
                    },
                };
                DecomposedRegion {
                    output_index: i,
                    kind: ro.kind,
                    back,
                    front,
                }
            })
            .collect()
    }
}

/// Exterior ring plus any interior (hole) rings of a polygon in 3D.
#[derive(Debug, Clone)]
pub struct PolygonRings {
    pub exterior: Vec<Point3<f32>>,
    pub interiors: Vec<Vec<Point3<f32>>>,
}

/// One output's contribution to the spatial decomposition of a parent
/// beam's column. See `BeamEvent::decompose`.
#[derive(Debug, Clone)]
pub struct DecomposedRegion {
    /// Index into the parent event's `outputs`.
    pub output_index: usize,
    pub kind: OutputKind,
    /// Polygon on the input face's plane.
    pub back: PolygonRings,
    /// Polygon on the output's actual plane (equal to `back` for OutGoing).
    pub front: PolygonRings,
}

impl DecomposedRegion {
    /// Side-wall quads of the prism connecting `back` to `front`. Each
    /// quad is `[back[i], back[i+1], front[i+1], front[i]]` — i.e. one
    /// quad per polygon edge, wound consistently. Exterior walls first,
    /// then any interior (hole) walls in declaration order. Returns an
    /// empty vec when `back == front` (OutGoing — the prism has zero
    /// length and no side faces).
    pub fn side_quads(&self) -> Vec<[Point3<f32>; 4]> {
        if matches!(self.kind, OutputKind::OutGoing) {
            return Vec::new();
        }
        let mut quads = Vec::new();
        push_ring_quads(&self.back.exterior, &self.front.exterior, &mut quads);
        for (b, f) in self.back.interiors.iter().zip(self.front.interiors.iter()) {
            push_ring_quads(b, f, &mut quads);
        }
        quads
    }
}

fn polygon_rings_of(face: &Face) -> PolygonRings {
    match face {
        Face::Simple(data) => PolygonRings {
            exterior: data.exterior.clone(),
            interiors: Vec::new(),
        },
        Face::Complex { data, interiors } => PolygonRings {
            exterior: data.exterior.clone(),
            interiors: interiors.clone(),
        },
    }
}

/// Back-project a polygon ring onto `plane` along `dir`. Each input
/// point `p` becomes `p + t*dir` where `t` solves
/// `(p + t*dir) · plane.normal + plane.offset = 0`.
fn back_project_ring(ring: &[Point3<f32>], plane: &Plane, dir: Vector3<f32>) -> Vec<Point3<f32>> {
    let n_dot_d = plane.normal.dot(&dir);
    // Ray parallel to plane should never happen in practice — the input
    // beam's prop direction is incident on its own face plane at some
    // nonzero angle (otherwise the beam has zero csa).
    debug_assert!(
        n_dot_d.abs() > 1e-6,
        "back_project_ring: ray direction parallel to target plane"
    );
    ring.iter()
        .map(|p| {
            let t = -(plane.normal.dot(&p.coords) + plane.offset) / n_dot_d;
            *p + dir * t
        })
        .collect()
}

fn push_ring_quads(back: &[Point3<f32>], front: &[Point3<f32>], out: &mut Vec<[Point3<f32>; 4]>) {
    debug_assert_eq!(
        back.len(),
        front.len(),
        "back and front rings must have matching vertex counts"
    );
    let n = back.len();
    if n < 2 {
        return;
    }
    for i in 0..n {
        let j = (i + 1) % n;
        out.push([back[i], back[j], front[j], front[i]]);
    }
}

/// Owned log of every propagation event from a recorded near-field solve.
#[derive(Debug, Default)]
pub struct Recording {
    pub events: Vec<BeamEvent>,
}

impl Recording {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    pub fn get(&self, id: EventId) -> Option<&BeamEvent> {
        self.events.get(id)
    }

    /// Event that emitted `beam_id` as one of its outputs.
    /// Linear scan — fine for typical recording sizes; introduce a
    /// `HashMap<BeamId, EventId>` index if this becomes hot.
    pub fn producer_event(&self, beam_id: BeamId) -> Option<EventId> {
        self.events
            .iter()
            .find(|e| e.outputs.iter().any(|o| o.beam.id == beam_id))
            .map(|e| e.id)
    }

    /// Event whose input is `beam_id` (i.e. the event that subsequently
    /// propagated this beam).
    pub fn consumer_event(&self, beam_id: BeamId) -> Option<EventId> {
        self.events
            .iter()
            .find(|e| e.input.id == beam_id)
            .map(|e| e.id)
    }

    /// Parent event of `id` — the event that produced this event's input
    /// beam as one of its outputs. None for the initial illumination event.
    pub fn parent(&self, id: EventId) -> Option<EventId> {
        let beam_id = self.events.get(id)?.input.id;
        self.producer_event(beam_id)
    }

    /// Direct child events of `id`, in event order.
    pub fn children(&self, id: EventId) -> Vec<EventId> {
        self.events
            .iter()
            .filter(|e| self.parent(e.id) == Some(id))
            .map(|e| e.id)
            .collect()
    }

    /// Begin a builder-style filtered view.
    pub fn view(&self) -> BeamView<'_> {
        BeamView::new(self)
    }
}

// ============================================================================
// BeamView — builder-style filtered view
// ============================================================================

#[derive(Default, Clone)]
struct Filter {
    rec_range: Option<(i32, i32)>,
    tir_range: Option<(i32, i32)>,
    depth_range: Option<(i32, i32)>,
    min_power: Option<f32>,
    variants: Option<Vec<BeamVariant>>,
    kinds: Option<Vec<OutputKind>>,
    descendants_of_beam: Option<BeamId>,
    ancestors_of_beam: Option<BeamId>,
}

/// Filtered view over a `Recording`. Chain builder methods, then iterate
/// or `collect_*` to materialise.
pub struct BeamView<'a> {
    recording: &'a Recording,
    filter: Filter,
}

impl<'a> BeamView<'a> {
    fn new(recording: &'a Recording) -> Self {
        Self {
            recording,
            filter: Filter::default(),
        }
    }

    // --- builders ---

    pub fn recursion(mut self, lo: i32, hi: i32) -> Self {
        self.filter.rec_range = Some((lo, hi));
        self
    }

    pub fn tir(mut self, lo: i32, hi: i32) -> Self {
        self.filter.tir_range = Some((lo, hi));
        self
    }

    pub fn depth(mut self, lo: i32, hi: i32) -> Self {
        self.filter.depth_range = Some((lo, hi));
        self
    }

    pub fn min_power(mut self, p: f32) -> Self {
        self.filter.min_power = Some(p);
        self
    }

    pub fn variant(mut self, v: BeamVariant) -> Self {
        self.filter.variants.get_or_insert_with(Vec::new).push(v);
        self
    }

    pub fn kind(mut self, k: OutputKind) -> Self {
        self.filter.kinds.get_or_insert_with(Vec::new).push(k);
        self
    }

    pub fn descendants_of_beam(mut self, b: BeamId) -> Self {
        self.filter.descendants_of_beam = Some(b);
        self
    }

    pub fn ancestors_of_beam(mut self, b: BeamId) -> Self {
        self.filter.ancestors_of_beam = Some(b);
        self
    }

    // --- iteration (borrowed) ---

    /// Events whose input beam passes the filter.
    pub fn events(&self) -> impl Iterator<Item = &BeamEvent> + '_ {
        let ancestry = self.ancestry_set();
        self.recording
            .events
            .iter()
            .filter(move |e| self.event_passes(e, ancestry.as_ref()))
    }

    /// Input beams of matching events.
    pub fn inputs(&self) -> impl Iterator<Item = &Beam> + '_ {
        self.events().map(|e| &e.input)
    }

    /// Output beams across matching events, paired with their parent event.
    /// The `kind` and `min_power` filters apply here (in addition to the
    /// event-level filters).
    pub fn outputs(&self) -> impl Iterator<Item = (&BeamEvent, &RecordedOutput)> + '_ {
        let ancestry = self.ancestry_set();
        self.recording
            .events
            .iter()
            .filter(move |e| self.event_passes(e, ancestry.as_ref()))
            .flat_map(|e| e.outputs.iter().map(move |o| (e, o)))
            .filter(move |(_, o)| self.output_passes(o))
    }

    // --- materialisation ---

    pub fn collect_events(&self) -> Vec<BeamEvent> {
        self.events().cloned().collect()
    }

    pub fn collect_inputs(&self) -> Vec<Beam> {
        self.inputs().cloned().collect()
    }

    pub fn collect_outputs(&self) -> Vec<RecordedOutput> {
        self.outputs().map(|(_, o)| o.clone()).collect()
    }

    // --- aggregation ---

    /// Sum the field at `point` over the filtered beams.
    /// TODO: physics — pick a convention (naive GO plane-wave evaluation
    /// vs Kirchhoff/Green-theorem). Stub for now.
    pub fn field_at(&self, _point: Point3<f32>) -> Option<Ampl> {
        todo!("field-at-point evaluation not yet implemented")
    }

    // --- internals ---

    /// If a descendants_of / ancestors_of filter is set, materialise the
    /// set of admissible event ids once so per-event predicate is O(1).
    fn ancestry_set(&self) -> Option<std::collections::HashSet<EventId>> {
        if self.filter.descendants_of_beam.is_none() && self.filter.ancestors_of_beam.is_none() {
            return None;
        }
        let mut set = std::collections::HashSet::new();

        if let Some(beam_id) = self.filter.descendants_of_beam {
            // Root = event that consumed beam_id. If no such event, no
            // descendants exist (beam was terminal).
            if let Some(root) = self.recording.consumer_event(beam_id) {
                // BFS over children.
                let mut stack = vec![root];
                while let Some(id) = stack.pop() {
                    if set.insert(id) {
                        stack.extend(self.recording.children(id));
                    }
                }
            }
        }

        if let Some(beam_id) = self.filter.ancestors_of_beam {
            // Walk from the event that consumed beam_id (or produced it,
            // if it never re-propagated) up through parent links.
            let start = self
                .recording
                .consumer_event(beam_id)
                .or_else(|| self.recording.producer_event(beam_id));
            let mut cur = start;
            while let Some(id) = cur {
                // If we're intersecting two ancestry filters, only keep events
                // already present in `set` (filled by descendants pass above).
                let intersect_mode = self.filter.descendants_of_beam.is_some();
                if intersect_mode {
                    if !set.contains(&id) {
                        // Stop walking — outside descendants subtree.
                        break;
                    }
                } else {
                    set.insert(id);
                }
                cur = self.recording.parent(id);
            }
        }

        Some(set)
    }

    fn event_passes(
        &self,
        e: &BeamEvent,
        ancestry: Option<&std::collections::HashSet<EventId>>,
    ) -> bool {
        if let Some(set) = ancestry {
            if !set.contains(&e.id) {
                return false;
            }
        }
        let b = &e.input;
        if let Some((lo, hi)) = self.filter.rec_range {
            if b.rec_count < lo || b.rec_count > hi {
                return false;
            }
        }
        if let Some((lo, hi)) = self.filter.tir_range {
            if b.tir_count < lo || b.tir_count > hi {
                return false;
            }
        }
        if let Some((lo, hi)) = self.filter.depth_range {
            let d = e.depth();
            if d < lo || d > hi {
                return false;
            }
        }
        if let Some(p) = self.filter.min_power {
            if b.power() < p {
                return false;
            }
        }
        if let Some(ref vs) = self.filter.variants {
            if !vs.iter().any(|v| v == &b.variant) {
                return false;
            }
        }
        true
    }

    fn output_passes(&self, o: &RecordedOutput) -> bool {
        if let Some(ref ks) = self.filter.kinds {
            if !ks.iter().any(|k| *k == o.kind) {
                return false;
            }
        }
        if let Some(p) = self.filter.min_power {
            if o.beam.power() < p {
                return false;
            }
        }
        true
    }
}
