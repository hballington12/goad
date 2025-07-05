#!/usr/bin/env python3
"""
Simple GOAD example demonstrating the clean new API.

This shows how easy it is to get started with GOAD light scattering simulations.
"""

import goad_py as goad
import os

# Get path to example geometry
current_dir = os.path.dirname(os.path.abspath(__file__))
geom_path = os.path.join(current_dir, "..", "examples", "data", "hex.obj")
geom_path = os.path.abspath(geom_path)

print("🔬 Simple GOAD Light Scattering Example")
print("=" * 50)

# Ultra-simple usage - just provide geometry!
print("\n1. Creating simulation with minimal code:")
print(f"   settings = goad.Settings('{os.path.basename(geom_path)}')")

settings = goad.Settings(geom_path)
mp = goad.MultiProblem(settings)

print(f"   ✅ Created simulation with {mp.num_orientations} random orientation")
print(f"   ✅ Using {len(mp.results.bins)} angular bins (interval binning)")

# Run the simulation
print("\n2. Running simulation:")
print("   mp.py_solve()")

mp.py_solve()

# Get results
results = mp.results
print("   ✅ Simulation completed!")

# Display some results
print("\n3. Results:")
print(f"   • Scattering cross-section: {results.scat_cross:.6f}")
print(f"   • Extinction cross-section: {results.ext_cross:.6f}")
absorption_cross = results.ext_cross - results.scat_cross
print(f"   • Absorption cross-section: {absorption_cross:.6f}")
print(f"   • Single scattering albedo: {results.albedo:.6f}")
print(f"   • Asymmetry parameter: {results.asymmetry:.6f}")

print("\n4. Power balance:")
powers = results.powers
print(f"   • Input power: {powers['input']:.6f}")
print(f"   • Output power: {powers['output']:.6f}")
print(f"   • Absorbed power: {powers['absorbed']:.6f}")
print(f"   • Missing power: {powers['missing']:.6f}")

print("\n🎉 That's it! Just 3 lines of code:")
print("   settings = goad.Settings('particle.obj')")
print("   mp = goad.MultiProblem(settings)")
print("   mp.py_solve()")
print("\n   All the complexity is handled by sensible defaults!")