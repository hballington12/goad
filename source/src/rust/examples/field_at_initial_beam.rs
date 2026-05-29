//! Sample the near-field at a handful of query points along the initial
//! incident beam.
//!
//! Runs a recorded solve on a hex geometry, filters the recording down to
//! just the initial illumination event, then asks `BeamView::fields_at`
//! for the field contribution at each query point.
//!
//! `settings.scale = Some(1.0)` opts out of the automatic rescale-to-unit-
//! cube that `Problem::init` would otherwise apply, so query points and
//! reported distances are in the same units as the geometry file.

use goad::beam::BeamVariant;
use goad::cancel::CancelToken;
use goad::geom::Geom;
use goad::problem::Problem;
use goad::settings::load_default_config;
use nalgebra::Point3;

fn main() {
    // 1. Load default settings and pin the scale so the geometry isn't
    //    rescaled. Coordinates below are then in the hex.obj's native
    //    units (vertices at z = ±5, hex radius ≈ 5).
    let mut settings = load_default_config().expect("load default config");
    settings.scale = Some(1.0);

    let geoms = Geom::load("./examples/data/hex.obj").expect("load hex");
    let geom = geoms[0].clone();
    let mut problem = Problem::new(Some(geom), Some(settings)).expect("build problem");

    // 2. Run the full recorded pipeline (init → orient → illuminate →
    //    near-field with recording → far-field → mueller → params).
    //    `None` for euler → no rotation; geometry sits as loaded.
    let recording = problem
        .run_with_recording(None, &CancelToken::noop())
        .expect("recorded solve");

    println!("recorded {} propagation events", recording.len());

    // 3. Filter the view to just the initial illumination event.
    let view = recording.view().variant(BeamVariant::Initial);
    let initial_events: Vec<_> = view.events().collect();
    println!("initial events in view: {}", initial_events.len());

    // 4. Pick a few query points. The default illumination travels along
    //    -z from a quad above the geometry. With FAC=1.1 in
    //    basic_initial_beam and hex top at z = 5, the illumination face
    //    sits at z ≈ 5.5, so the initial beam's column spans
    //    z ∈ [5.0, 5.5]. We sample inside that slab, plus two negative
    //    controls (upstream of the start, off-footprint in xy).
    let query_points = vec![
        Point3::new(0.0, 0.0, 5.4),  // inside the initial slab, near start
        Point3::new(0.0, 0.0, 5.2),  // inside, mid-slab
        Point3::new(0.0, 0.0, 5.05), // inside, near front cap
        Point3::new(0.0, 0.0, 10.0), // upstream of illumination — outside
        Point3::new(0.0, 0.0, 0.0),  // inside the geometry — past front cap
        Point3::new(20.0, 0.0, 5.2), // off the illumination footprint in xy
    ];

    let results = view.fields_at(&query_points);

    // 5. Report per-point contributions.
    for (i, contribs) in results.iter().enumerate() {
        let p = &query_points[i];
        println!(
            "\npoint ({:.3}, {:.3}, {:.3}): {} contribution(s)",
            p.x,
            p.y,
            p.z,
            contribs.len()
        );
        for c in contribs {
            println!(
                "  event {} output {} kind={:?} distance={:.4} intensity={:.4e} phase={:.4}",
                c.event_id,
                c.output_index,
                c.kind,
                c.distance,
                c.field.intensity(),
                c.field.phase(),
            );
        }
    }
}
