//! Tests for the 32-bit Owen-scrambled Sobol sequence (`goad::sobol`).
//!
//! Regression target: `sobol_burley` 0.5 supported only 2^16 points, so long
//! convergence runs re-weighted the same 65,536 orientations.

use goad::orientation::{OrientationSampler, Orientations};
use goad::sobol::{euler_angles, sample, sobol_u32};
use std::collections::HashSet;

const N16: u64 = 1 << 16;

fn point_bits(index: u64, seed: u64) -> (u32, u32, u32) {
    let (a, b, g) = euler_angles(index, seed);
    (a.to_bits(), b.to_bits(), g.to_bits())
}

#[test]
fn direction_numbers_match_sobol_burley_top_16_bits() {
    // sobol_burley uses the same Joe-Kuo numbers truncated to 16 bits.
    for d in 0..3u32 {
        for i in 0..N16 as u32 {
            let ours = sobol_u32(i, d as usize) >> 16;
            let theirs = sobol_burley::parts::sobol_rev(i.reverse_bits(), d).reverse_bits() >> 16;
            assert_eq!(ours, theirs, "dimension {d} index {i}");
        }
    }
}

#[test]
fn distinct_points_beyond_2_pow_16() {
    // The old sampler gave exactly 65,536 distinct points here.
    for seed in [0u64, 1] {
        let n = 4 * N16;
        let set: HashSet<_> = (0..n).map(|i| point_bits(i, seed)).collect();
        assert_eq!(set.len() as u64, n, "seed {seed}");
    }
}

#[test]
fn sampler_draws_new_orientations_after_2_pow_16() {
    let mut sampler = OrientationSampler::sobol(Some(0));
    let first: HashSet<_> = (0..N16)
        .map(|_| {
            let e = sampler.next().unwrap();
            (e.alpha.to_bits(), e.beta.to_bits(), e.gamma.to_bits())
        })
        .collect();
    let repeats = (0..N16)
        .filter(|_| {
            let e = sampler.next().unwrap();
            first.contains(&(e.alpha.to_bits(), e.beta.to_bits(), e.gamma.to_bits()))
        })
        .count();
    assert_eq!(repeats, 0);
}

/// Every elementary interval of area 2^-m in dims (0, 1) holds exactly one point.
fn assert_02_net(start: u64, m: u32, seed: u64) {
    let n = 1u64 << m;
    for a in 0..=m {
        let b = m - a;
        let mut counts = vec![0u32; n as usize];
        for i in start..start + n {
            let x = ((sample(i, 0, seed) as f64) * (1u64 << a) as f64) as u64;
            let y = ((sample(i, 1, seed) as f64) * (1u64 << b) as f64) as u64;
            counts[(x * (1u64 << b) + y) as usize] += 1;
        }
        assert!(counts.iter().all(|&c| c == 1), "start {start} m {m} split ({a},{b})");
    }
}

#[test]
fn prefixes_and_blocks_are_02_nets() {
    for seed in [0u64, 7] {
        assert_02_net(0, 12, seed);
        assert_02_net(N16, 12, seed); // aligned block past the old limit
        assert_02_net(5 * N16 + (3 << 12), 12, seed);
    }
}

#[test]
fn seeds_are_independent() {
    let a: HashSet<_> = (0..N16).map(|i| point_bits(i, 1)).collect();
    let shared = (0..N16).filter(|&i| a.contains(&point_bits(i, 2))).count();
    assert_eq!(shared, 0);
}

#[test]
fn orientations_are_uniform() {
    let n = 1u64 << 20;
    let (mut s_alpha, mut s_cosb, mut s_gamma) = (0f64, 0f64, 0f64);
    for i in 0..n {
        let (a, b, g) = euler_angles(i, 3);
        assert!((0.0..360.0).contains(&a) && (0.0..=180.0).contains(&b) && (0.0..360.0).contains(&g));
        s_alpha += a as f64;
        s_cosb += (b as f64).to_radians().cos();
        s_gamma += g as f64;
    }
    let n = n as f64;
    assert!((s_alpha / n - 180.0).abs() < 0.01);
    assert!((s_cosb / n).abs() < 1e-4);
    assert!((s_gamma / n - 180.0).abs() < 0.01);
}

#[test]
fn block_past_2_pow_32_uses_fresh_points() {
    let base = 1u64 << 32;
    let early: HashSet<_> = (0..4096).map(|i| point_bits(i, 0)).collect();
    let shared = (base..base + 4096).filter(|&i| early.contains(&point_bits(i, 0))).count();
    assert_eq!(shared, 0);
}

#[test]
fn batch_and_streaming_samplers_agree() {
    let batch = Orientations::sobol(1000, Some(11));
    let mut stream = OrientationSampler::sobol(Some(11));
    for &(a, b, g) in &batch.eulers {
        let e = stream.next().unwrap();
        assert_eq!((a, b, g), (e.alpha, e.beta, e.gamma));
    }
}
