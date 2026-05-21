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

use geo::Contains;
use nalgebra::{Matrix4, Point3, Vector3};

use crate::beam::{Beam, BeamId, BeamVariant, OutputKind};
use crate::field::Field;
use crate::geom::{look_along, Face, Plane};

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

/// One beam's contribution to the field at a query point. Produced by
/// `BeamView::field_at` / `fields_at`. The `field` is the parent beam's
/// field propagated forward to the query point through the beam's medium.
#[derive(Debug, Clone)]
pub struct FieldContribution {
    /// Event that produced this contribution (the parent beam being
    /// propagated through the query point's location).
    pub event_id: EventId,
    /// Index into `event.outputs` of the region whose prism the query
    /// point fell inside — useful for attributing the contribution to a
    /// specific downstream branch (e.g. refracted vs reflected).
    pub output_index: usize,
    /// Kind of the output region (mirrors `event.outputs[output_index].kind`).
    pub kind: OutputKind,
    /// Parent beam's field, propagated forward by `distance` through its
    /// medium.
    pub field: Field,
    /// Signed propagation distance from the input-face midpoint to the
    /// query point, along the parent beam's prop direction. This is the
    /// argument passed to `Field::propagate` when computing the
    /// contribution.
    pub distance: f32,
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

/// Build a 2D polygon in the clip frame from a 3D `PolygonRings` by
/// transforming each vertex through `transform` and dropping the z
/// coordinate. Used for point-in-polygon tests against the back cap of a
/// decomposed region.
fn ring_to_polygon_2d(rings: &PolygonRings, transform: &Matrix4<f32>) -> geo::Polygon<f32> {
    let to_coord = |p: &Point3<f32>| {
        let pc = transform.transform_point(p);
        geo::Coord { x: pc.x, y: pc.y }
    };
    let exterior: Vec<_> = rings.exterior.iter().map(to_coord).collect();
    let interiors: Vec<_> = rings
        .interiors
        .iter()
        .map(|ring| geo::LineString(ring.iter().map(to_coord).collect()))
        .collect();
    geo::Polygon::new(geo::LineString(exterior), interiors)
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

    /// Evaluate the field at a single query point, summing contributions
    /// from every beam in the filtered view whose column contains `point`.
    /// Returns one entry per contributing beam — callers can inspect them
    /// individually or sum after rotating to a common reference frame.
    /// See `fields_at` for the batch entry point.
    pub fn field_at(&self, point: Point3<f32>) -> Vec<FieldContribution> {
        self.fields_at(&[point]).pop().unwrap_or_default()
    }

    /// Batch field evaluation. Returns one `Vec<FieldContribution>` per
    /// query point, in input order. Iterates events once and tests every
    /// point against each event's prism regions, which is `O(events ×
    /// points)`. No spatial index yet; revisit if it becomes hot.
    ///
    /// For each event, a point lies inside an output's prism when:
    /// 1. It is downstream of the back (input) face *plane* and (for
    ///    non-OutGoing outputs) upstream of the front (output) face
    ///    *plane*, both measured along prop. The plane equation is used
    ///    rather than a midpoint-along-prop projection so the test holds
    ///    at oblique incidence — at the back plane the boundary cuts
    ///    obliquely across the prism's columnar interior.
    ///    OutGoing prisms are semi-infinite forward (geometric-optics
    ///    column, no diffraction smearing).
    /// 2. It falls inside the back polygon's xy footprint in the event's
    ///    clip frame.
    ///
    /// When both hold, the parent beam's field is propagated forward by
    /// `(X - input_mid)·prop` — the natural plane-wave phase variation
    /// from the reference midpoint where `event.input.field` is wound to.
    /// This is distinct from the slab-test distance above; at oblique
    /// incidence the back plane and prop are not aligned, so "distance
    /// along prop from the plane" and "distance along prop from a point
    /// on the plane" disagree.
    pub fn fields_at(&self, points: &[Point3<f32>]) -> Vec<Vec<FieldContribution>> {
        let mut out = vec![Vec::<FieldContribution>::new(); points.len()];

        for event in self.events() {
            let prop = event.input.field.prop();
            let n = event.input.refr_index;
            let k = event.input.wavenumber();
            let input_mid = event.input.face.data().midpoint;
            let transform = look_along(prop);

            // Back plane = input face's plane. True for every output kind:
            // OutGoing remainders lie on the input plane by construction;
            // Internal / ExternalDiff regions are back-projected onto it.
            let back_plane = event.input.face.plane();
            let n_b_dot_prop = back_plane.normal.dot(&prop);

            // Pre-transform query points into the clip frame once per
            // event so the per-region inner loop only does the xy test.
            let points_clip: Vec<Point3<f32>> = points
                .iter()
                .map(|p| transform.transform_point(p))
                .collect();

            for region in event.decompose() {
                let semi_infinite = matches!(region.kind, OutputKind::OutGoing);
                // Front plane info — None for OutGoing (no downstream cap).
                let front_plane_info = if semi_infinite {
                    None
                } else {
                    let plane = event.outputs[region.output_index].beam.face.plane();
                    let n_f_dot_prop = plane.normal.dot(&prop);
                    Some((plane, n_f_dot_prop))
                };
                let back_poly_2d = ring_to_polygon_2d(&region.back, &transform);

                for (i, &x) in points.iter().enumerate() {
                    // Slab test using plane equations. `t` is the signed
                    // distance from the plane to X along prop — positive
                    // means X is downstream of the plane.
                    let t_back = (x.coords.dot(&back_plane.normal) + back_plane.offset)
                        / n_b_dot_prop;
                    if t_back < 0.0 {
                        continue;
                    }
                    if let Some((ref fp, n_f_dot_prop)) = front_plane_info {
                        let t_front = (x.coords.dot(&fp.normal) + fp.offset) / n_f_dot_prop;
                        if t_front > 0.0 {
                            continue;
                        }
                    }

                    let x_clip = points_clip[i];
                    let p2d = geo::Point::new(x_clip.x, x_clip.y);
                    if !back_poly_2d.contains(&p2d) {
                        continue;
                    }

                    // Phase distance = signed projection of (X - input_mid)
                    // onto prop. Always relative to input_mid (where
                    // event.input.field is wound to) regardless of kind.
                    let dist = (x - input_mid).dot(&prop);
                    let mut field = event.input.field.clone();
                    field.propagate(dist, k, n);

                    out[i].push(FieldContribution {
                        event_id: event.id,
                        output_index: region.output_index,
                        kind: region.kind,
                        field,
                        distance: dist,
                    });
                }
            }
        }

        out
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
