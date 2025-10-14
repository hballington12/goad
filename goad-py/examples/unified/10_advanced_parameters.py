#!/usr/bin/env python3
"""
Example 10: Advanced Parameters - Comprehensive Demonstration

Demonstrates ALL modifiable parameters in run_convergence():
- Optical parameters (wavelength, refractive indices)
- Convergence parameters (batch size, max orientations, min batches)
- Beam tracing parameters (thresholds, recursion limits, cutoff)
- Geometry transformations (problem scale, per-axis scaling, distortion)
- Advanced configuration (mapping, coherence, field of view)
- Reproducibility (random seed)

This example showcases the full flexibility of the unified convergence API.
"""

from pathlib import Path
import goad_py as goad

print("=" * 80)
print("Example 10: Advanced Parameters - ALL Modifiable Parameters")
print("=" * 80)

# Get paths relative to this script's location
script_dir = Path(__file__).parent
geom_file = script_dir / "../../../examples/data/hex.obj"
phips_bins_file = script_dir / "../../../phips_bins_edges.toml"

# Resolve to absolute paths
geom_file = str(geom_file.resolve())
phips_bins_file = str(phips_bins_file.resolve())

# Configuration summary
print("\n" + "=" * 80)
print("Configuration Summary:")
print("=" * 80)
print("OPTICAL PARAMETERS:")
print("  Wavelength: 0.633 μm (HeNe laser)")
print("  Particle refractive index: 1.5 + 0.01i (absorbing particle)")
print("  Medium refractive index: 1.33 + 0.0i (water)")
print()
print("CONVERGENCE PARAMETERS:")
print("  Batch size: 12 orientations per batch")
print("  Max orientations: 500")
print("  Min batches: 3")
print()
print("BEAM TRACING PARAMETERS:")
print("  Beam power threshold: 0.01 (stricter than default 0.05)")
print("  Beam area threshold factor: 2.0 (half of default 4.0)")
print("  Cutoff: 0.0005 (stricter than default 0.001)")
print("  Max recursion depth: 6 (limited for demonstration)")
print("  Max TIR bounces: 50 (half of default 100)")
print()
print("GEOMETRY TRANSFORMATIONS:")
print("  Problem scale: 2.0 (scales entire problem by 2x)")
print("  Per-axis geometry scale: [1.0, 1.0, 1.5] (stretch z-axis by 1.5x)")
print("  Distortion: 0.1 (adds 10% geometric distortion)")
print()
print("ADVANCED CONFIGURATION:")
print("  Mapping: GeometricOptics (instead of default ApertureDiffraction)")
print("  Coherence: True (enable coherent scattering calculations)")
print("  Field of view factor: 2.0 (doubled FOV)")
print()
print("REPRODUCIBILITY:")
print("  Random seed: 42 (for reproducible results)")
print()
print("CONVERGENCE TARGET:")
print("  Mode: PHIPS detector convergence")
print("  Detectors: Indices 1-19 (26° to 170°, skipping 18° forward scatter)")
print("  Tolerance: 50% relative SEM (relaxed for demo)")
print("=" * 80 + "\n")

# Converge on all PHIPS detectors except the first one (indices 1-19)
detector_indices = list(range(1, 20))  # [1, 2, 3, ..., 19]

print(f"Starting convergence on {len(detector_indices)} PHIPS detectors...\n")

# Run PHIPS convergence with ALL MODIFIABLE PARAMETERS
results = goad.run_convergence(
    # Required parameters
    geometry=geom_file,
    targets=[
        {
            "tolerance": 0.25,  # 50% relative SEM (relaxed for demo)
            "tolerance_type": "relative",
            "detector_indices": detector_indices,
        }
    ],
    mode=goad.PHIPSMode(bins_file=phips_bins_file),
    # ========================================================================
    # OPTICAL PARAMETERS
    # ========================================================================
    wavelength=0.633,  # HeNe laser wavelength in microns
    particle_refr_index_re=1.5,  # Real part of particle refractive index
    particle_refr_index_im=0.01,  # Imaginary part (absorption)
    medium_refr_index_re=1.33,  # Real part of medium refractive index (water)
    medium_refr_index_im=0.0,  # Imaginary part (no absorption in medium)
    # ========================================================================
    # CONVERGENCE PARAMETERS
    # ========================================================================
    batch_size=12,  # Orientations per batch (smaller batches for finer control)
    max_orientations=500,  # Maximum total orientations
    min_batches=10,  # Minimum batches before checking convergence
    mueller_1d=False,  # Mueller matrix output (N/A for PHIPS mode)
    # ========================================================================
    # BEAM TRACING PARAMETERS
    # ========================================================================
    beam_power_threshold=0.01,  # Beam power threshold (stricter than default)
    beam_area_threshold_fac=0.001,  # Beam area threshold factor (tighter threshold)
    cutoff=0.0005,  # Ray power cutoff (more accurate but slower)
    max_rec=6,  # Maximum recursion depth (limited for demo)
    max_tir=20,  # Maximum total internal reflection bounces
    # ========================================================================
    # GEOMETRY TRANSFORMATIONS
    # ========================================================================
    scale=1.0,  # Problem scaling factor (numerical changes only)
    geom_scale=[3.0, 3.0, 6.0],  # Per-axis geometry scaling
    distortion=0.0,  # Geometry distortion factor
    # ========================================================================
    # ADVANCED CONFIGURATION
    # ========================================================================
    mapping=goad.Mapping.ApertureDiffraction,  # or GeometricOptics
    coherence=True,  # Enable coherent scattering calculations
    fov_factor=2.0,  # Field of view factor (doubled)
)

# Print summary
print(results.summary())

# Save results
results.save("10_advanced_parameters_results.json")

# Show detailed PHIPS output
print("\n" + "=" * 80)
print("PHIPS Detector DSCS Values (all 20 detectors)")
print("=" * 80)

import json

with open("10_advanced_parameters_results.json", "r") as f:
    data = json.load(f)

phips_dscs = data.get("phips_dscs")
detector_angles = data.get("detector_angles")

if phips_dscs and detector_angles:
    print(
        f"\n{'Detector':<10} {'Angle':<10} {'DSCS Mean':<15} {'DSCS SEM':<15} {'Rel. SEM':<10}"
    )
    print("-" * 70)

    for i, (dscs_data, angle) in enumerate(zip(phips_dscs, detector_angles)):
        mean, sem = dscs_data[0], dscs_data[1]
        rel_sem = (sem / abs(mean) * 100) if mean != 0 else float("inf")

        # Mark which detectors were converged on
        marker = "✓" if i in detector_indices else "○"

        print(
            f"{marker} {i:<8} {angle:<10.1f} {mean:<15.4e} {sem:<15.4e} {rel_sem:<10.2f}%"
        )

print("\n" + "=" * 80)
print("Legend:")
print("  ✓ = Converged on this detector")
print("  ○ = Not a convergence target (but still computed)")
print("=" * 80)

# Show configuration that was used
print("\n" + "=" * 80)
print("Configuration Used (from saved results):")
print("=" * 80)

config = data.get("config", {})
if config:
    print(f"\nOptical parameters:")
    print(f"  Wavelength: {config.get('wavelength')} μm")
    print(
        f"  Particle refractive index: {config.get('particle_refr_index_re')} + {config.get('particle_refr_index_im')}i"
    )
    print(
        f"  Medium refractive index: {config.get('medium_refr_index_re')} + {config.get('medium_refr_index_im')}i"
    )

    beam_tracing = config.get("beam_tracing", {})
    if beam_tracing:
        print(f"\nBeam tracing parameters:")
        print(f"  Max recursion depth: {beam_tracing.get('max_rec')}")
        print(f"  Max TIR bounces: {beam_tracing.get('max_tir')}")
        print(f"  Cutoff: {beam_tracing.get('cutoff')}")
        print(f"  Beam power threshold: {beam_tracing.get('beam_power_threshold')}")
        print(
            f"  Beam area threshold factor: {beam_tracing.get('beam_area_threshold_fac')}"
        )

    geom_transform = config.get("geometry_transform", {})
    if geom_transform:
        print(f"\nGeometry transformations:")
        print(f"  Problem scale: {geom_transform.get('scale')}")
        if geom_transform.get("distortion"):
            print(f"  Distortion: {geom_transform.get('distortion')}")
        if geom_transform.get("geom_scale"):
            print(f"  Per-axis geometry scale: {geom_transform.get('geom_scale')}")

    advanced = config.get("advanced_config", {})
    if advanced:
        print(f"\nAdvanced configuration:")
        print(f"  Mapping: {advanced.get('mapping')}")
        print(f"  Coherence: {advanced.get('coherence')}")
        if advanced.get("fov_factor"):
            print(f"  Field of view factor: {advanced.get('fov_factor')}")

    if config.get("seed"):
        print(f"\nReproducibility:")
        print(f"  Random seed: {config.get('seed')}")

print("\n✓ Example completed successfully!")
