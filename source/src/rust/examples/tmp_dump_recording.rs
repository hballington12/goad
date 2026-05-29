//! Temporary: dump geometry + recording to JSON for the python visualiser.
//! Not intended to be committed.

use goad::cancel::CancelToken;
use goad::geom::{Face, Geom};
use goad::orientation::Euler;
use goad::problem::Problem;
use goad::recording::Recording;
use goad::settings::load_default_config;
use serde_json::json;
use std::fs::File;
use std::io::Write;

fn face_exterior(face: &Face) -> Vec<[f32; 3]> {
    face.data()
        .exterior
        .iter()
        .map(|p| [p.x, p.y, p.z])
        .collect()
}

fn dump_geometry(geom: &Geom) -> serde_json::Value {
    let shapes: Vec<_> = geom
        .shapes
        .iter()
        .map(|s| {
            let faces: Vec<Vec<[f32; 3]>> = s.faces.iter().map(face_exterior).collect();
            json!({ "faces": faces })
        })
        .collect();
    json!({ "shapes": shapes })
}

fn ring(points: &[nalgebra::Point3<f32>]) -> Vec<[f32; 3]> {
    points.iter().map(|p| [p.x, p.y, p.z]).collect()
}

fn dump_recording(rec: &Recording) -> serde_json::Value {
    let events: Vec<_> = rec
        .events
        .iter()
        .map(|e| {
            let prop = e.input.field.prop();
            let regions: Vec<_> = e
                .decompose()
                .into_iter()
                .map(|r| {
                    let back_interiors: Vec<Vec<[f32; 3]>> =
                        r.back.interiors.iter().map(|ring_pts| ring(ring_pts)).collect();
                    let front_interiors: Vec<Vec<[f32; 3]>> =
                        r.front.interiors.iter().map(|ring_pts| ring(ring_pts)).collect();
                    json!({
                        "output_index": r.output_index,
                        "kind": format!("{:?}", r.kind),
                        "back_exterior": ring(&r.back.exterior),
                        "back_interiors": back_interiors,
                        "front_exterior": ring(&r.front.exterior),
                        "front_interiors": front_interiors,
                    })
                })
                .collect();
            json!({
                "id": e.id,
                "input_variant": format!("{:?}", e.input.variant),
                "input_face_exterior": face_exterior(&e.input.face),
                "input_prop": [prop.x, prop.y, prop.z],
                "input_rec_count": e.input.rec_count,
                "input_tir_count": e.input.tir_count,
                "regions": regions,
            })
        })
        .collect();
    json!({ "events": events })
}

fn main() {
    let mut settings = load_default_config().expect("load default config");
    settings.scale = Some(1.0);

    let geoms = Geom::load("./examples/data/hex.obj").expect("load hex");
    let geom = geoms[0].clone();
    let mut problem = Problem::new(Some(geom), Some(settings)).expect("build problem");

    let euler = Euler::new(0.0, 30.0, 30.0);
    let recording = problem
        .run_with_recording(Some(&euler), &CancelToken::noop())
        .expect("recorded solve");

    let dump = json!({
        "geometry": dump_geometry(&problem.geom),
        "recording": dump_recording(&recording),
    });

    let path = "tmp/recording.json";
    let mut f = File::create(path).expect("create json");
    f.write_all(serde_json::to_string_pretty(&dump).unwrap().as_bytes())
        .expect("write json");
    println!(
        "wrote {} ({} events, {} shapes)",
        path,
        recording.len(),
        problem.geom.shapes.len()
    );
}
