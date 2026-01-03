//! Tests for ConvergenceTracker merge functionality.
//!
//! These tests verify that Chan's parallel algorithm produces identical
//! results to sequential Welford updates.
//!
//! Run with: cargo test --release -p goad --test convergence_merge_tests -- --ignored

use goad::{
    convergence::ConvergenceTracker,
    orientation::Euler,
    problem::Problem,
    result::{GOComponent, Results},
    settings,
};

/// Generate deterministic results using fixed orientations.
fn generate_results(count: usize) -> Vec<Results> {
    let settings = settings::load_default_config().unwrap();

    // Fixed orientations for reproducibility
    let eulers: Vec<Euler> = (0..count)
        .map(|i| {
            let alpha = (i as f32 * 37.0) % 360.0;
            let beta = (i as f32 * 23.0) % 180.0;
            let gamma = (i as f32 * 41.0) % 360.0;
            Euler::new(alpha, beta, gamma)
        })
        .collect();

    eulers
        .iter()
        .map(|euler| {
            let mut problem = Problem::new(None, Some(settings.clone())).unwrap();
            problem.run(Some(euler)).unwrap();
            problem.result
        })
        .collect()
}

/// Helper to compare two Results for approximate equality
fn results_approx_eq(a: &Results, b: &Results, epsilon: f32) -> bool {
    // Compare key scalar parameters
    let params_to_check = [
        (
            a.params.asymmetry(&GOComponent::Total),
            b.params.asymmetry(&GOComponent::Total),
        ),
        (
            a.params.scatt_cross(&GOComponent::Total),
            b.params.scatt_cross(&GOComponent::Total),
        ),
        (
            a.params.ext_cross(&GOComponent::Total),
            b.params.ext_cross(&GOComponent::Total),
        ),
        (
            a.params.albedo(&GOComponent::Total),
            b.params.albedo(&GOComponent::Total),
        ),
    ];

    for (val_a, val_b) in params_to_check {
        match (val_a, val_b) {
            (Some(va), Some(vb)) => {
                let diff = (va - vb).abs();
                let max_val = va.abs().max(vb.abs()).max(1e-10);
                if diff / max_val > epsilon {
                    eprintln!(
                        "Mismatch: {} vs {}, relative diff: {}",
                        va,
                        vb,
                        diff / max_val
                    );
                    return false;
                }
            }
            (None, None) => {}
            _ => return false,
        }
    }

    // Compare powers
    let power_diff = (a.powers.input - b.powers.input).abs()
        + (a.powers.output - b.powers.output).abs()
        + (a.powers.absorbed - b.powers.absorbed).abs();
    if power_diff > epsilon {
        eprintln!("Powers mismatch: {:?} vs {:?}", a.powers, b.powers);
        return false;
    }

    true
}

/// Test that merging two trackers gives same result as sequential updates.
#[test]
#[ignore] // Run with: cargo test --release -- --ignored
fn test_merge_equals_sequential() {
    let results = generate_results(10);

    // Method 1: Sequential update
    let mut sequential = ConvergenceTracker::new(&results[0]);
    for r in &results {
        sequential.update(r);
    }

    // Method 2: Split 5+5 and merge
    let mut batch_a = ConvergenceTracker::new(&results[0]);
    let mut batch_b = ConvergenceTracker::new(&results[0]);
    for r in &results[0..5] {
        batch_a.update(r);
    }
    for r in &results[5..10] {
        batch_b.update(r);
    }
    batch_a.merge(&batch_b);

    // Compare means
    let seq_mean = sequential.mean();
    let merged_mean = batch_a.mean();

    assert!(
        results_approx_eq(&seq_mean, &merged_mean, 1e-5),
        "Mean mismatch between sequential and merged (5+5)"
    );

    // Compare SEMs
    let seq_sem = sequential.sem();
    let merged_sem = batch_a.sem();

    assert!(
        results_approx_eq(&seq_sem, &merged_sem, 1e-5),
        "SEM mismatch between sequential and merged (5+5)"
    );

    // Verify counts match
    assert_eq!(sequential.count(), batch_a.count());

    println!("Sequential count: {}", sequential.count());
    println!(
        "Sequential asymmetry mean: {:?}",
        seq_mean.params.asymmetry(&GOComponent::Total)
    );
    println!(
        "Merged asymmetry mean: {:?}",
        merged_mean.params.asymmetry(&GOComponent::Total)
    );
}

/// Test merging with different batch sizes (3+4+3).
#[test]
#[ignore]
fn test_merge_uneven_batches() {
    let results = generate_results(10);

    // Sequential baseline
    let mut sequential = ConvergenceTracker::new(&results[0]);
    for r in &results {
        sequential.update(r);
    }

    // Three uneven batches: 3 + 4 + 3
    let mut batch_a = ConvergenceTracker::new(&results[0]);
    let mut batch_b = ConvergenceTracker::new(&results[0]);
    let mut batch_c = ConvergenceTracker::new(&results[0]);

    for r in &results[0..3] {
        batch_a.update(r);
    }
    for r in &results[3..7] {
        batch_b.update(r);
    }
    for r in &results[7..10] {
        batch_c.update(r);
    }

    // Merge in sequence
    batch_a.merge(&batch_b);
    batch_a.merge(&batch_c);

    let seq_mean = sequential.mean();
    let merged_mean = batch_a.mean();

    assert!(
        results_approx_eq(&seq_mean, &merged_mean, 1e-5),
        "Mean mismatch between sequential and merged (3+4+3)"
    );

    let seq_sem = sequential.sem();
    let merged_sem = batch_a.sem();

    assert!(
        results_approx_eq(&seq_sem, &merged_sem, 1e-5),
        "SEM mismatch between sequential and merged (3+4+3)"
    );

    assert_eq!(sequential.count(), batch_a.count());
}

/// Test merging single-element batches (degenerate case).
#[test]
#[ignore]
fn test_merge_single_elements() {
    let results = generate_results(5);

    // Sequential baseline
    let mut sequential = ConvergenceTracker::new(&results[0]);
    for r in &results {
        sequential.update(r);
    }

    // Merge one at a time (simulates batch_size=1)
    let mut merged = ConvergenceTracker::new(&results[0]);
    for r in &results {
        let mut single = ConvergenceTracker::new(&results[0]);
        single.update(r);
        merged.merge(&single);
    }

    let seq_mean = sequential.mean();
    let merged_mean = merged.mean();

    assert!(
        results_approx_eq(&seq_mean, &merged_mean, 1e-5),
        "Mean mismatch for single-element merges"
    );

    let seq_sem = sequential.sem();
    let merged_sem = merged.sem();

    assert!(
        results_approx_eq(&seq_sem, &merged_sem, 1e-5),
        "SEM mismatch for single-element merges"
    );
}

/// Test merging empty tracker (edge case).
#[test]
#[ignore]
fn test_merge_empty_tracker() {
    let results = generate_results(5);

    let mut tracker = ConvergenceTracker::new(&results[0]);
    for r in &results {
        tracker.update(r);
    }

    let mean_before = tracker.mean();
    let count_before = tracker.count();

    // Merge empty tracker - should be no-op
    let empty = ConvergenceTracker::new(&results[0]);
    tracker.merge(&empty);

    assert_eq!(tracker.count(), count_before);
    assert!(results_approx_eq(&tracker.mean(), &mean_before, 1e-10));
}

/// Test merging into empty tracker.
#[test]
#[ignore]
fn test_merge_into_empty() {
    let results = generate_results(5);

    let mut source = ConvergenceTracker::new(&results[0]);
    for r in &results {
        source.update(r);
    }

    let mut empty = ConvergenceTracker::new(&results[0]);
    empty.merge(&source);

    assert_eq!(empty.count(), source.count());
    assert!(results_approx_eq(&empty.mean(), &source.mean(), 1e-10));
    assert!(results_approx_eq(&empty.sem(), &source.sem(), 1e-10));
}

/// Test with larger sample size (100 orientations) and various batch sizes.
#[test]
#[ignore]
fn test_merge_large_sample() {
    let results = generate_results(100);

    // Sequential baseline
    let mut sequential = ConvergenceTracker::new(&results[0]);
    for r in &results {
        sequential.update(r);
    }

    // Simulate what workers would do: batches of 10
    let batch_size = 10;
    let mut batched = ConvergenceTracker::new(&results[0]);
    for chunk in results.chunks(batch_size) {
        let mut batch_tracker = ConvergenceTracker::new(&results[0]);
        for r in chunk {
            batch_tracker.update(r);
        }
        batched.merge(&batch_tracker);
    }

    let seq_mean = sequential.mean();
    let batched_mean = batched.mean();
    let seq_sem = sequential.sem();
    let batched_sem = batched.sem();

    println!("\n=== Large Sample Test (100 orientations, batch_size=10) ===");
    println!(
        "Count: seq={}, batched={}",
        sequential.count(),
        batched.count()
    );
    println!(
        "Asymmetry Mean: seq={:?}, batched={:?}",
        seq_mean.params.asymmetry(&GOComponent::Total),
        batched_mean.params.asymmetry(&GOComponent::Total)
    );
    println!(
        "Asymmetry SEM:  seq={:?}, batched={:?}",
        seq_sem.params.asymmetry(&GOComponent::Total),
        batched_sem.params.asymmetry(&GOComponent::Total)
    );

    assert!(
        results_approx_eq(&seq_mean, &batched_mean, 1e-5),
        "Mean mismatch in large sample test"
    );

    assert!(
        results_approx_eq(&seq_sem, &batched_sem, 1e-5),
        "SEM mismatch in large sample test"
    );

    assert_eq!(sequential.count(), batched.count());

    // Also test with uneven batch sizes (simulating partial batches)
    let mut uneven_batched = ConvergenceTracker::new(&results[0]);
    let batch_sizes = [7, 13, 23, 17, 11, 9, 8, 6, 4, 2]; // sums to 100
    let mut offset = 0;
    for &size in &batch_sizes {
        let mut batch_tracker = ConvergenceTracker::new(&results[0]);
        for r in &results[offset..offset + size] {
            batch_tracker.update(r);
        }
        uneven_batched.merge(&batch_tracker);
        offset += size;
    }

    let uneven_mean = uneven_batched.mean();
    let uneven_sem = uneven_batched.sem();

    println!("\n=== Uneven batches (7+13+23+17+11+9+8+6+4+2) ===");
    println!(
        "Asymmetry Mean: seq={:?}, uneven={:?}",
        seq_mean.params.asymmetry(&GOComponent::Total),
        uneven_mean.params.asymmetry(&GOComponent::Total)
    );
    println!(
        "Asymmetry SEM:  seq={:?}, uneven={:?}",
        seq_sem.params.asymmetry(&GOComponent::Total),
        uneven_sem.params.asymmetry(&GOComponent::Total)
    );

    assert!(
        results_approx_eq(&seq_mean, &uneven_mean, 1e-5),
        "Mean mismatch with uneven batches"
    );

    assert!(
        results_approx_eq(&seq_sem, &uneven_sem, 1e-5),
        "SEM mismatch with uneven batches"
    );
}

/// Test the batch size calculation formula.
/// batch_size = ceil(num_workers * merge_time / compute_time * 1.2)
#[test]
fn test_batch_size_formula() {
    // Helper to compute batch size using the same formula as run_prognosis
    fn compute_batch_size(num_workers: usize, compute_time_ns: u64, merge_time_ns: u64) -> usize {
        const MIN_BATCH_SIZE: usize = 2;
        const MAX_BATCH_SIZE: usize = 64;

        if compute_time_ns == 0 {
            return 10; // fallback
        }

        let ratio = (num_workers as f64 * merge_time_ns as f64) / compute_time_ns as f64;
        let batch_size = (ratio * 1.2).ceil() as usize;
        batch_size.clamp(MIN_BATCH_SIZE, MAX_BATCH_SIZE)
    }

    // Case 1: Master easily keeps up (slow compute, fast merge)
    // compute=50ms, merge=0.5ms, 8 workers
    // ratio = 8 * 0.5 / 50 = 0.08 → batch_size = ceil(0.08 * 1.2) = 1 → clamped to 2
    let batch = compute_batch_size(8, 50_000_000, 500_000);
    assert_eq!(batch, 2, "Slow compute should yield MIN_BATCH_SIZE");

    // Case 2: Master is bottleneck (fast compute, slow merge)
    // compute=5ms, merge=0.5ms, 32 workers
    // ratio = 32 * 0.5 / 5 = 3.2 → batch_size = ceil(3.2 * 1.2) = 4
    let batch = compute_batch_size(32, 5_000_000, 500_000);
    assert_eq!(
        batch, 4,
        "Fast compute with many workers needs larger batch"
    );

    // Case 3: Very fast compute (high throughput)
    // compute=1ms, merge=0.5ms, 64 workers
    // ratio = 64 * 0.5 / 1 = 32 → batch_size = ceil(32 * 1.2) = 39
    let batch = compute_batch_size(64, 1_000_000, 500_000);
    assert_eq!(batch, 39, "Very fast compute needs big batches");

    // Case 4: Extreme case - hits MAX_BATCH_SIZE
    // compute=0.1ms, merge=1ms, 64 workers
    // ratio = 64 * 1 / 0.1 = 640 → clamped to 64
    let batch = compute_batch_size(64, 100_000, 1_000_000);
    assert_eq!(batch, 64, "Extreme case should clamp to MAX_BATCH_SIZE");

    // Case 5: Single worker
    // compute=10ms, merge=1ms, 1 worker
    // ratio = 1 * 1 / 10 = 0.1 → batch_size = ceil(0.1 * 1.2) = 1 → clamped to 2
    let batch = compute_batch_size(1, 10_000_000, 1_000_000);
    assert_eq!(batch, 2, "Single worker should yield MIN_BATCH_SIZE");

    // Case 6: Realistic scenario
    // compute=20ms, merge=0.2ms, 16 workers
    // ratio = 16 * 0.2 / 20 = 0.16 → batch_size = ceil(0.16 * 1.2) = 1 → clamped to 2
    let batch = compute_batch_size(16, 20_000_000, 200_000);
    assert_eq!(batch, 2, "Realistic GOAD scenario with slow compute");

    println!("All batch size formula tests passed!");
}

/// Print detailed comparison for debugging.
#[test]
#[ignore]
fn test_merge_detailed_comparison() {
    let results = generate_results(10);

    let mut sequential = ConvergenceTracker::new(&results[0]);
    for r in &results {
        sequential.update(r);
    }

    let mut batch_a = ConvergenceTracker::new(&results[0]);
    let mut batch_b = ConvergenceTracker::new(&results[0]);
    for r in &results[0..5] {
        batch_a.update(r);
    }
    for r in &results[5..10] {
        batch_b.update(r);
    }
    batch_a.merge(&batch_b);

    let seq_mean = sequential.mean();
    let merged_mean = batch_a.mean();
    let seq_sem = sequential.sem();
    let merged_sem = batch_a.sem();

    println!("\n=== Detailed Comparison ===");
    println!(
        "Count: seq={}, merged={}",
        sequential.count(),
        batch_a.count()
    );

    println!("\nAsymmetry:");
    println!(
        "  Mean: seq={:?}, merged={:?}",
        seq_mean.params.asymmetry(&GOComponent::Total),
        merged_mean.params.asymmetry(&GOComponent::Total)
    );
    println!(
        "  SEM:  seq={:?}, merged={:?}",
        seq_sem.params.asymmetry(&GOComponent::Total),
        merged_sem.params.asymmetry(&GOComponent::Total)
    );

    println!("\nScatCross:");
    println!(
        "  Mean: seq={:?}, merged={:?}",
        seq_mean.params.scatt_cross(&GOComponent::Total),
        merged_mean.params.scatt_cross(&GOComponent::Total)
    );
    println!(
        "  SEM:  seq={:?}, merged={:?}",
        seq_sem.params.scatt_cross(&GOComponent::Total),
        merged_sem.params.scatt_cross(&GOComponent::Total)
    );

    println!("\nPowers:");
    println!(
        "  Input:  seq={}, merged={}",
        seq_mean.powers.input, merged_mean.powers.input
    );
    println!(
        "  Output: seq={}, merged={}",
        seq_mean.powers.output, merged_mean.powers.output
    );
}
