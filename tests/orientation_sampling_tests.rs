//! Tests for orientation sampling methods (Uniform, Sobol, Halton)

use goad::orientation::{OrientationSampler, Orientations, Scheme};

/// Test that Sobol sampler produces deterministic results with same seed
#[test]
fn test_sobol_sampler_deterministic() {
    let mut sampler1 = OrientationSampler::sobol(Some(42));
    let mut sampler2 = OrientationSampler::sobol(Some(42));

    for _ in 0..10 {
        let e1 = sampler1.next().unwrap();
        let e2 = sampler2.next().unwrap();
        assert_eq!(e1.alpha, e2.alpha);
        assert_eq!(e1.beta, e2.beta);
        assert_eq!(e1.gamma, e2.gamma);
    }
}

/// Test that Halton sampler produces deterministic results
#[test]
fn test_halton_sampler_deterministic() {
    let mut sampler1 = OrientationSampler::halton();
    let mut sampler2 = OrientationSampler::halton();

    for _ in 0..10 {
        let e1 = sampler1.next().unwrap();
        let e2 = sampler2.next().unwrap();
        assert_eq!(e1.alpha, e2.alpha);
        assert_eq!(e1.beta, e2.beta);
        assert_eq!(e1.gamma, e2.gamma);
    }
}

/// Test that different seeds produce different sequences
#[test]
fn test_sobol_different_seeds() {
    let mut sampler1 = OrientationSampler::sobol(Some(1));
    let mut sampler2 = OrientationSampler::sobol(Some(2));

    let e1 = sampler1.next().unwrap();
    let e2 = sampler2.next().unwrap();

    // Different seeds should produce different results
    assert!(e1.alpha != e2.alpha || e1.beta != e2.beta || e1.gamma != e2.gamma);
}

/// Test that sampler reset works correctly for Sobol
#[test]
fn test_sobol_sampler_reset() {
    let mut sampler = OrientationSampler::sobol(Some(42));

    let first_batch: Vec<_> = (0..5).map(|_| sampler.next().unwrap()).collect();

    sampler.reset();

    let second_batch: Vec<_> = (0..5).map(|_| sampler.next().unwrap()).collect();

    for (e1, e2) in first_batch.iter().zip(second_batch.iter()) {
        assert_eq!(e1.alpha, e2.alpha);
        assert_eq!(e1.beta, e2.beta);
        assert_eq!(e1.gamma, e2.gamma);
    }
}

/// Test that sampler reset works correctly for Halton
#[test]
fn test_halton_sampler_reset() {
    let mut sampler = OrientationSampler::halton();

    let first_batch: Vec<_> = (0..5).map(|_| sampler.next().unwrap()).collect();

    sampler.reset();

    let second_batch: Vec<_> = (0..5).map(|_| sampler.next().unwrap()).collect();

    for (e1, e2) in first_batch.iter().zip(second_batch.iter()) {
        assert_eq!(e1.alpha, e2.alpha);
        assert_eq!(e1.beta, e2.beta);
        assert_eq!(e1.gamma, e2.gamma);
    }
}

/// Test that Sobol produces angles in valid ranges
#[test]
fn test_sobol_angle_ranges() {
    let mut sampler = OrientationSampler::sobol(Some(123));

    for _ in 0..100 {
        let e = sampler.next().unwrap();
        assert!(
            e.alpha >= 0.0 && e.alpha <= 360.0,
            "alpha out of range: {}",
            e.alpha
        );
        assert!(
            e.beta >= 0.0 && e.beta <= 180.0,
            "beta out of range: {}",
            e.beta
        );
        assert!(
            e.gamma >= 0.0 && e.gamma <= 360.0,
            "gamma out of range: {}",
            e.gamma
        );
    }
}

/// Test that Halton produces angles in valid ranges
#[test]
fn test_halton_angle_ranges() {
    let mut sampler = OrientationSampler::halton();

    for _ in 0..100 {
        let e = sampler.next().unwrap();
        assert!(
            e.alpha >= 0.0 && e.alpha <= 360.0,
            "alpha out of range: {}",
            e.alpha
        );
        assert!(
            e.beta >= 0.0 && e.beta <= 180.0,
            "beta out of range: {}",
            e.beta
        );
        assert!(
            e.gamma >= 0.0 && e.gamma <= 360.0,
            "gamma out of range: {}",
            e.gamma
        );
    }
}

/// Test batch generation with Sobol scheme
#[test]
fn test_orientations_sobol_batch() {
    let scheme = Scheme::Sobol { num_orients: 100 };
    let orientations = Orientations::generate(&scheme, Some(42));

    assert_eq!(orientations.num_orientations, 100);
    assert_eq!(orientations.eulers.len(), 100);

    // Check all angles are in valid ranges
    for (alpha, beta, gamma) in &orientations.eulers {
        assert!(*alpha >= 0.0 && *alpha <= 360.0);
        assert!(*beta >= 0.0 && *beta <= 180.0);
        assert!(*gamma >= 0.0 && *gamma <= 360.0);
    }
}

/// Test batch generation with Halton scheme
#[test]
fn test_orientations_halton_batch() {
    let scheme = Scheme::Halton { num_orients: 100 };
    let orientations = Orientations::generate(&scheme, None);

    assert_eq!(orientations.num_orientations, 100);
    assert_eq!(orientations.eulers.len(), 100);

    // Check all angles are in valid ranges
    for (alpha, beta, gamma) in &orientations.eulers {
        assert!(*alpha >= 0.0 && *alpha <= 360.0);
        assert!(*beta >= 0.0 && *beta <= 180.0);
        assert!(*gamma >= 0.0 && *gamma <= 360.0);
    }
}

/// Test that Sobol batch is deterministic with same seed
#[test]
fn test_orientations_sobol_deterministic() {
    let scheme = Scheme::Sobol { num_orients: 50 };
    let o1 = Orientations::generate(&scheme, Some(999));
    let o2 = Orientations::generate(&scheme, Some(999));

    assert_eq!(o1.eulers, o2.eulers);
}

/// Test that Halton batch is deterministic (no seed needed)
#[test]
fn test_orientations_halton_deterministic() {
    let scheme = Scheme::Halton { num_orients: 50 };
    let o1 = Orientations::generate(&scheme, None);
    let o2 = Orientations::generate(&scheme, None);

    assert_eq!(o1.eulers, o2.eulers);
}

/// Test that Sobol and Halton produce different sequences
#[test]
fn test_sobol_halton_differ() {
    let mut sobol = OrientationSampler::sobol(Some(0));
    let mut halton = OrientationSampler::halton();

    let e_sobol = sobol.next().unwrap();
    let e_halton = halton.next().unwrap();

    // They should produce different results (extremely unlikely to be equal)
    assert!(
        e_sobol.alpha != e_halton.alpha
            || e_sobol.beta != e_halton.beta
            || e_sobol.gamma != e_halton.gamma
    );
}

/// Test that uniform and quasi-random samplers produce valid orientations
#[test]
fn test_all_samplers_valid_orientations() {
    let samplers: Vec<(&str, OrientationSampler)> = vec![
        ("uniform", OrientationSampler::uniform(Some(42))),
        ("sobol", OrientationSampler::sobol(Some(42))),
        ("halton", OrientationSampler::halton()),
    ];

    for (name, mut sampler) in samplers {
        for i in 0..50 {
            let e = sampler
                .next()
                .expect(&format!("{} sampler exhausted at {}", name, i));

            assert!(
                e.alpha >= 0.0 && e.alpha <= 360.0,
                "{} alpha out of range at {}: {}",
                name,
                i,
                e.alpha
            );
            assert!(
                e.beta >= 0.0 && e.beta <= 180.0,
                "{} beta out of range at {}: {}",
                name,
                i,
                e.beta
            );
            assert!(
                e.gamma >= 0.0 && e.gamma <= 360.0,
                "{} gamma out of range at {}: {}",
                name,
                i,
                e.gamma
            );
        }
    }
}
