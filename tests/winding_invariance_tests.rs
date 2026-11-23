//! The result must not depend on the winding of the input mesh: an outward,
//! an inward and a mixed-winding copy of the same geometry give the same
//! cross sections and the same optical-theorem extinction.

use goad::{
    geom::Geom,
    multiproblem::MultiProblem,
    orientation::{Euler, EulerConvention, Orientation, Scheme as OrientationScheme},
    params::Param,
    result::GOComponent,
    settings,
    zones::ZoneType,
};
use num_complex::Complex32;
use std::path::PathBuf;

/// Copy of an OBJ with the vertex order of selected faces reversed.
fn rewound_copy(src: &str, tag: &str, reverse: impl Fn(usize) -> bool) -> PathBuf {
    let text = std::fs::read_to_string(src).expect("read obj");
    let mut face_index = 0;
    let lines: Vec<String> = text
        .lines()
        .map(|line| {
            if let Some(rest) = line.strip_prefix("f ") {
                let mut verts: Vec<&str> = rest.split_whitespace().collect();
                if reverse(face_index) {
                    verts.reverse();
                }
                face_index += 1;
                format!("f {}", verts.join(" "))
            } else {
                line.to_string()
            }
        })
        .collect();
    let out = std::env::temp_dir().join(format!("goad_winding_{tag}_{}.obj", std::process::id()));
    std::fs::write(&out, lines.join("\n") + "\n").expect("write obj");
    out
}

/// (C_sca integrated, C_ext optical theorem, input power) for one orientation.
fn run(path: &str) -> (f32, f32, f32) {
    let mut settings = settings::load_default_config().unwrap();
    settings.wavelength = 0.532;
    settings.medium_refr_index = Complex32::new(1.0, 0.0);
    settings.orientation = Orientation {
        scheme: OrientationScheme::Discrete {
            eulers: vec![Euler::new(20.0, 50.0, 10.0)],
        },
        euler_convention: EulerConvention::ZYZ,
    };
    let geoms = Geom::load(path, vec![Complex32::new(1.31, 0.0)]).expect("load geom");
    let mut mp = MultiProblem::new(geoms, Some(settings)).expect("multiproblem");
    mp.solve();
    let zone = |t: ZoneType| mp.result.zones.iter().find(|z| z.zone_type == t).unwrap();
    let sca = zone(ZoneType::Full).params.get(&Param::ScatCross, &GOComponent::Total).unwrap();
    let ot = zone(ZoneType::Forward)
        .params
        .get(&Param::ExtCrossOpticalTheorem, &GOComponent::Total)
        .unwrap();
    (sca, ot, mp.result.powers.input)
}

fn assert_close(a: (f32, f32, f32), b: (f32, f32, f32), what: &str) {
    for (x, y, name) in [(a.0, b.0, "C_sca"), (a.1, b.1, "C_ext OT"), (a.2, b.2, "input power")] {
        assert!((x - y).abs() <= 1e-4 * x.abs().max(1.0), "{what}: {name} {x} vs {y}");
    }
}

#[test]
fn test_results_independent_of_mesh_winding() {
    for src in ["examples/data/hex.obj", "examples/data/concave1.obj"] {
        let reference = run(src);
        let inward = rewound_copy(src, "inward", |_| true);
        let mixed = rewound_copy(src, "mixed", |i| i % 3 == 1);
        assert_close(reference, run(inward.to_str().unwrap()), &format!("{src} inward"));
        assert_close(reference, run(mixed.to_str().unwrap()), &format!("{src} mixed"));
        let _ = std::fs::remove_file(inward);
        let _ = std::fs::remove_file(mixed);
    }
}
