//! Dump Re(E_x) on a 2D grid at y=0 for the near-far-mapping schematic in
//! the goad-bench manuscript.
//!
//! Mirrors the schematic's geometry:
//!   * `examples/data/cube.obj` rescaled to half-side 1.5
//!   * rotated 15deg about y (XYZ Euler convention)
//!   * illuminated along -z with wavelength 0.532
//!   * particle refractive index 1.31 (real), medium 1.0
//!
//! Output JSON (`tmp/schematic_field.json`) carries:
//!   * the rotated cube's vertices/faces, so the python plotter can draw
//!     the silhouette without redoing the rotation;
//!   * the initial-beam clip face (i.e. F_0), so the TikZ schematic can be
//!     aligned to whatever z-plane goad's `basic_initial_beam` chose;
//!   * a flat row-major `nx*nz` array of Re(E_x) summed over every
//!     contributing beam at each (x, 0, z) sample.

use goad::cancel::CancelToken;
use goad::geom::Geom;
use goad::orientation::{Euler, EulerConvention, Orientation, Scheme};
use goad::problem::Problem;
use goad::settings::load_default_config;
use nalgebra::Point3;
use num_complex::Complex;
use serde_json::json;
use std::fs::File;
use std::io::Write;

const CUBE_HALF_NATIVE: f32 = 6.025435; // half-side of examples/data/cube.obj
const HALF_SIDE_TARGET: f32 = 1.5; // matches the TikZ schematic
const WAVELENGTH: f32 = 0.532;
const N_REFR_RE: f32 = 1.31;
const EULER_BETA_DEG: f32 = 15.0;

const NX: usize = 480;
const NZ: usize = 520;
const X_MIN: f32 = -3.0;
const X_MAX: f32 = 3.0;
const Z_MIN: f32 = -3.0;
const Z_MAX: f32 = 3.5;

// Wider field grid for the diffraction figure -- same y-slice, but a
// 20x20 extent so OutGoing beams are visible well into the surrounding
// air. Coarser angular resolution per wavelength than the close-up
// grid, but still plenty for the diffraction visualisation.
const WIDE_NX: usize = 800;
const WIDE_NZ: usize = 800;
const WIDE_X_MIN: f32 = -10.0;
const WIDE_X_MAX: f32 = 10.0;
const WIDE_Z_MIN: f32 = -10.0;
const WIDE_Z_MAX: f32 = 10.0;

fn main() {
    // 1. Settings.
    let mut settings = load_default_config().expect("load default config");
    settings.scale = Some(1.0); // bypass auto-rescale-to-unit-cube
    settings.wavelength = WAVELENGTH;
    settings.medium_refr_index = Complex::new(1.0, 0.0);

    let s = HALF_SIDE_TARGET / CUBE_HALF_NATIVE;
    settings.geom_scale = Some(vec![s, s, s]);
    settings.orientation = Orientation {
        scheme: Scheme::Discrete {
            eulers: vec![Euler::new(0.0, EULER_BETA_DEG, 0.0)],
        },
        euler_convention: EulerConvention::XYZ,
    };

    // 2. Build problem and run the recorded pipeline.
    let geom_path = "./examples/data/cube.obj";
    let geoms = Geom::load(geom_path, vec![Complex::new(N_REFR_RE, 0.0)]).expect("load cube.obj");
    let geom = geoms[0].clone();
    let mut problem = Problem::new(geom, Some(settings)).expect("build problem");

    let euler = Euler::new(0.0, EULER_BETA_DEG, 0.0);
    let recording = problem
        .run_with_recording(Some(&euler), &CancelToken::noop())
        .expect("recorded solve");

    // 3. Build sampling grid in (x, z) at y=0.
    let dx = (X_MAX - X_MIN) / (NX as f32 - 1.0);
    let dz = (Z_MAX - Z_MIN) / (NZ as f32 - 1.0);
    let mut points = Vec::with_capacity(NX * NZ);
    for j in 0..NZ {
        let z = Z_MIN + dz * j as f32;
        for i in 0..NX {
            let x = X_MIN + dx * i as f32;
            points.push(Point3::new(x, 0.0, z));
        }
    }

    // 4. Identify the sub-beam pathway we want to visualise: the
    //    side-facet (T--L) refraction from the initial event, and its
    //    immediate child event (the beam propagating inside the cube
    //    until its first exit). Other branches -- the top-facet (T--R)
    //    refraction, external-diffraction remainders, deeper internal
    //    bounces -- are deliberately excluded so the figure shows one
    //    clean pathway.
    let initial_id: usize = 0;
    let initial = &recording.events[initial_id];

    let tl_output_idx = initial
        .outputs
        .iter()
        .enumerate()
        .filter(|(_, o)| matches!(o.kind, goad::beam::OutputKind::NearField))
        .min_by(|(_, a), (_, b)| {
            let ax = a.beam.face.data().midpoint.x;
            let bx = b.beam.face.data().midpoint.x;
            ax.partial_cmp(&bx).unwrap()
        })
        .map(|(i, _)| i)
        .expect("no Internal output on initial event");
    let tl_beam_id = initial.outputs[tl_output_idx].beam.id;
    let child_id = recording.consumer_event(tl_beam_id);

    // External reflection sibling: same facet as beam #4, opposite
    // handedness -- the Internal output on the same T-L face whose
    // prop has a positive component along the face's outward direction
    // (i.e. heading back into air).
    let tl_face_mid = initial.outputs[tl_output_idx].beam.face.data().midpoint;
    let refl_output_idx = initial
        .outputs
        .iter()
        .enumerate()
        .find(|(i, o)| {
            *i != tl_output_idx && matches!(o.kind, goad::beam::OutputKind::NearField) && {
                let face_mid = o.beam.face.data().midpoint;
                let same_face = (face_mid - tl_face_mid).norm() < 1e-3;
                let prop = o.beam.field.prop();
                // outward from cube centre (recentred at origin) to face mid
                let outward = face_mid;
                let outgoing = prop.x * outward.x + prop.y * outward.y + prop.z * outward.z > 0.0;
                same_face && outgoing
            }
        })
        .map(|(i, _)| i);
    let refl_child_id = refl_output_idx.and_then(|i| {
        let bid = initial.outputs[i].beam.id;
        recording.consumer_event(bid)
    });

    println!(
        "focus pathway: initial event {} output {} (T--L refraction), child event {:?}",
        initial_id, tl_output_idx, child_id
    );
    println!(
        "reflection sibling: output {:?} on initial event, child event {:?}",
        refl_output_idx, refl_child_id
    );

    // Event 3's grandchildren: the transmission out of L--B (beam #4
    // refracted into air below the cube) and the internal reflection
    // (beam #4 bouncing back inside). Distinguish by prop projection
    // onto L--B's outward direction (= face midpoint since the cube is
    // recentred at origin).
    let (trans_grandchild_id, refl_grandchild_id) =
        match child_id.and_then(|cid| recording.events.get(cid)) {
            Some(e3) => {
                let outward = e3.outputs[0].beam.face.data().midpoint;
                let mut t_idx: Option<usize> = None;
                let mut r_idx: Option<usize> = None;
                for (i, o) in e3.outputs.iter().enumerate() {
                    let p = o.beam.field.prop();
                    let dot = p.x * outward.x + p.y * outward.y + p.z * outward.z;
                    if dot > 0.0 {
                        t_idx = t_idx.or(Some(i));
                    } else {
                        r_idx = r_idx.or(Some(i));
                    }
                }
                let t = t_idx.and_then(|i| recording.consumer_event(e3.outputs[i].beam.id));
                let r = r_idx.and_then(|i| recording.consumer_event(e3.outputs[i].beam.id));
                (t, r)
            }
            None => (None, None),
        };
    println!(
        "depth-2 children of event 3: transmission {:?}, internal-reflection {:?}",
        trans_grandchild_id, refl_grandchild_id
    );

    // Depth-3 events: every direct child of either depth-2 event. These
    // are the next refractions / internal reflections after beam #4 has
    // either exited the cube or bounced internally once.
    let depth3_event_ids: Vec<usize> = [trans_grandchild_id, refl_grandchild_id]
        .iter()
        .copied()
        .flatten()
        .flat_map(|eid| recording.children(eid))
        .collect();
    println!("depth-3 events: {:?}", depth3_event_ids);

    // 5. Sample fields at every grid point (unfiltered view -- we filter
    //    contribution-by-contribution next).
    let view = recording.view();
    let results = view.fields_at(&points);

    // 6. Reduce each grid point to a scalar Re(E_x), keeping only the
    //    pathway selected above. With the default initial polarization
    //    (e_perp = x̂, prop = -ẑ), the wave is TM-polarized w.r.t. the
    //    facets whose normals lie in (x, z), so E lives in the (x, z)
    //    plane and E_x is the natural scalar to plot.
    let depth3_set: std::collections::HashSet<usize> = depth3_event_ids.iter().copied().collect();
    let in_pathway = |event_id: usize, output_index: usize| -> bool {
        (event_id == initial_id && output_index == tl_output_idx)
            || Some(event_id) == child_id
            || Some(event_id) == refl_child_id
            || Some(event_id) == trans_grandchild_id
            || Some(event_id) == refl_grandchild_id
            || depth3_set.contains(&event_id)
    };
    // Scalar field proxy: |E| * cos(phase). Projecting onto a fixed
    // global direction (x̂ etc.) introduces a sign flip at facets where
    // goad rotates the (e_perp, e_par) basis -- a basis artefact, not
    // physical. Plotting the magnitude of E times the oscillating phase
    // factor sidesteps that: |E|² = |ampl[0,0]|² + |ampl[1,0]|² for the
    // unit-perp input, which is gauge-invariant.
    //
    // Where two beams overlap, summing |E_i| * cos(phase_i) is *not* the
    // correct coherent sum -- you'd need to add the complex amplitudes
    // and take Re of the result. We're not doing that; instead, we pick
    // the lowest-depth contribution (initial illumination beats any of
    // its descendants), so each grid point shows exactly one beam's
    // wave. Same priority is used for the dumped phase so contour lines
    // (if revived later) and field bands stay consistent.
    let reduce = |contribs_grid: &[Vec<goad::recording::FieldContribution>]| -> (
        Vec<f32>,
        Vec<f32>,
        Vec<f32>,
    ) {
        let n = contribs_grid.len();
        let mut field_re = Vec::with_capacity(n);
        let mut phase = Vec::with_capacity(n);
        let mut field_mag = Vec::with_capacity(n);
        for contribs in contribs_grid {
            let best = contribs
                .iter()
                .filter(|c| in_pathway(c.event_id, c.output_index))
                .min_by_key(|c| recording.events[c.event_id].depth());
            match best {
                Some(c) => {
                    let a_perp = c.field.ampl_wo_phase()[(0, 0)];
                    let a_par = c.field.ampl_wo_phase()[(1, 0)];
                    let mag = (a_perp.norm_sqr() + a_par.norm_sqr()).sqrt();
                    let phi = c.field.phase();
                    field_re.push(mag * phi.cos());
                    phase.push(phi);
                    field_mag.push(mag);
                }
                None => {
                    field_re.push(0.0);
                    phase.push(f32::NAN);
                    field_mag.push(f32::NAN);
                }
            }
        }
        (field_re, phase, field_mag)
    };
    let (field_re, phase, field_mag) = reduce(&results);

    // Wide-grid sample for the diffraction figure.
    let wide_dx = (WIDE_X_MAX - WIDE_X_MIN) / (WIDE_NX as f32 - 1.0);
    let wide_dz = (WIDE_Z_MAX - WIDE_Z_MIN) / (WIDE_NZ as f32 - 1.0);
    let mut wide_points = Vec::with_capacity(WIDE_NX * WIDE_NZ);
    for j in 0..WIDE_NZ {
        let z = WIDE_Z_MIN + wide_dz * j as f32;
        for i in 0..WIDE_NX {
            let x = WIDE_X_MIN + wide_dx * i as f32;
            wide_points.push(Point3::new(x, 0.0, z));
        }
    }
    let wide_results = view.fields_at(&wide_points);
    let (wide_field_re, wide_phase, wide_field_mag) = reduce(&wide_results);

    // 6. Geometry dump (post-rotation vertices and per-face exteriors).
    let shapes: Vec<_> = problem
        .geom
        .shapes
        .iter()
        .map(|s| {
            let vertices: Vec<[f32; 3]> = s.vertices.iter().map(|v| [v.x, v.y, v.z]).collect();
            let faces: Vec<Vec<[f32; 3]>> = s
                .faces
                .iter()
                .map(|f| f.data().exterior.iter().map(|p| [p.x, p.y, p.z]).collect())
                .collect();
            json!({ "vertices": vertices, "faces": faces })
        })
        .collect();

    // 7. Initial-beam face (F_0) and its decomposition into sub-beams.
    let initial_face = recording.events.first().map(|e| {
        let exterior: Vec<[f32; 3]> = e
            .input
            .face
            .data()
            .exterior
            .iter()
            .map(|p| [p.x, p.y, p.z])
            .collect();
        let prop = e.input.field.prop();
        let decomp: Vec<_> = e
            .decompose()
            .into_iter()
            .map(|r| {
                let back: Vec<[f32; 3]> = r.back.exterior.iter().map(|p| [p.x, p.y, p.z]).collect();
                let front: Vec<[f32; 3]> =
                    r.front.exterior.iter().map(|p| [p.x, p.y, p.z]).collect();
                let out_prop = e.outputs[r.output_index].beam.field.prop();
                json!({
                    "output_index": r.output_index,
                    "kind": format!("{:?}", r.kind),
                    "back_exterior": back,
                    "front_exterior": front,
                    "prop": [out_prop.x, out_prop.y, out_prop.z],
                })
            })
            .collect();
        json!({
            "exterior": exterior,
            "prop": [prop.x, prop.y, prop.z],
            "decomposition": decomp,
        })
    });

    // Per-event decomposition dumper -- reusable for every event we
    // surface in the JSON.
    let dump_event = |eid: usize| -> serde_json::Value {
        let e = &recording.events[eid];
        let regions: Vec<_> = e
            .decompose()
            .into_iter()
            .map(|r| {
                let back: Vec<[f32; 3]> = r.back.exterior.iter().map(|p| [p.x, p.y, p.z]).collect();
                let front: Vec<[f32; 3]> =
                    r.front.exterior.iter().map(|p| [p.x, p.y, p.z]).collect();
                let out_prop = e.outputs[r.output_index].beam.field.prop();
                json!({
                    "output_index": r.output_index,
                    "kind": format!("{:?}", r.kind),
                    "back_exterior": back,
                    "front_exterior": front,
                    "prop": [out_prop.x, out_prop.y, out_prop.z],
                })
            })
            .collect();
        let input_prop = e.input.field.prop();
        json!({
            "event_id": e.id,
            "input_prop": [input_prop.x, input_prop.y, input_prop.z],
            "decomposition": regions,
        })
    };

    // 7b. Decomposition of the child event that consumed beam #4 (the
    //     T--L refraction). This gives us the prism of beam #4 propagating
    //     inside the cube until it hits its next facet.
    let child_decomp = child_id.map(dump_event);

    // 7c. Decomposition of the consumer of the T--L reflection (beam #5).
    //     This is the next event for the externally-reflected beam, whose
    //     prism extends from the T--L facet outward into air in the
    //     reflection direction.
    let refl_child_decomp = refl_child_id.map(dump_event);

    // 7d. Decompositions of event 3's two children: transmission out of
    //     the cube through L--B (refracted into air below) and the
    //     internal reflection (back inside the cube, heading to the next
    //     facet).
    let trans_grandchild_decomp = trans_grandchild_id.map(dump_event);
    let refl_grandchild_decomp = refl_grandchild_id.map(dump_event);

    // 7e. Depth-3 events -- one decomposition entry per descendant of a
    //     depth-2 event.
    let depth3_events_decomp: Vec<_> = depth3_event_ids.iter().copied().map(dump_event).collect();

    let dump = json!({
        "geometry": { "shapes": shapes },
        "initial_beam": initial_face,
        "child_event": child_decomp,
        "reflection_child_event": refl_child_decomp,
        "transmission_grandchild_event": trans_grandchild_decomp,
        "reflection_grandchild_event": refl_grandchild_decomp,
        "depth3_events": depth3_events_decomp,
        "reflection_initial_output_index": refl_output_idx,
        "grid": {
            "x_min": X_MIN, "x_max": X_MAX, "nx": NX,
            "z_min": Z_MIN, "z_max": Z_MAX, "nz": NZ,
            "field_re": field_re,
            "phase": phase,
            "field_mag": field_mag,
        },
        "wide_grid": {
            "x_min": WIDE_X_MIN, "x_max": WIDE_X_MAX, "nx": WIDE_NX,
            "z_min": WIDE_Z_MIN, "z_max": WIDE_Z_MAX, "nz": WIDE_NZ,
            "field_re": wide_field_re,
            "phase": wide_phase,
            "field_mag": wide_field_mag,
        },
        "params": {
            "wavelength": WAVELENGTH,
            "particle_refr_index_re": N_REFR_RE,
            "particle_refr_index_im": 0.0,
            "geom_scale": s,
            "half_side_target": HALF_SIDE_TARGET,
            "euler_alpha": 0.0,
            "euler_beta": EULER_BETA_DEG,
            "euler_gamma": 0.0,
            "euler_convention": "XYZ",
        }
    });

    std::fs::create_dir_all("tmp").ok();
    let path = "tmp/schematic_field.json";
    let mut f = File::create(path).expect("create json");
    f.write_all(serde_json::to_string(&dump).unwrap().as_bytes())
        .expect("write json");
    println!(
        "wrote {} ({} events, {}x{} grid)",
        path,
        recording.len(),
        NX,
        NZ
    );
}
