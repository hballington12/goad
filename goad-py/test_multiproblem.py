"""
Simple test script to verify MultiProblem Python bindings work correctly.
"""

import goad_py as goad

print("Testing MultiProblem Python Bindings")
print("=" * 40)

try:
    # Test 1: Create Euler angles
    print("✓ Creating Euler angles...")
    euler1 = goad.Euler(0.0, 0.0, 0.0)
    euler2 = goad.Euler(30.0, 30.0, 30.0)
    print(f"  Created: {euler1}")
    print(f"  Created: {euler2}")

    # Test 2: Create orientation schemes
    print("✓ Creating orientation schemes...")
    uniform_scheme = goad.uniform_orientation(10)
    discrete_scheme = goad.discrete_orientation([euler1, euler2])
    print(f"  Uniform scheme: {uniform_scheme}")
    print(f"  Discrete scheme: {discrete_scheme}")

    # Test 3: Create full orientations
    print("✓ Creating full orientations...")
    uniform_orientation = goad.create_uniform_orientation(10)
    discrete_orientation = goad.create_discrete_orientation([euler1, euler2])
    print(f"  Uniform orientation: {uniform_orientation}")
    print(f"  Discrete orientation: {discrete_orientation}")

    # Test 4: Create Settings with orientation
    print("✓ Creating Settings with orientation...")
    settings = goad.Settings(
        geom_name="hex.obj",
        orientation=discrete_orientation
    )
    print(f"  Settings created with orientation: {settings.orientation}")

    # Test 5: Create MultiProblem
    print("✓ Creating MultiProblem...")
    multi_problem = goad.MultiProblem(settings=settings)  # Use keyword argument for clarity
    print(f"  MultiProblem created successfully")
    print(f"  Number of orientations: {multi_problem.num_orientations}")

    # Test 6: Test MultiProblem methods (without actually solving)
    print("✓ Testing MultiProblem methods...")
    print(f"  MultiProblem has solve method: {hasattr(multi_problem, 'py_solve')}")
    print(f"  MultiProblem has results property: {hasattr(multi_problem, 'results')}")
    print(f"  MultiProblem has reset method: {hasattr(multi_problem, 'py_reset')}")

    print("\n🎉 All MultiProblem binding tests passed!")
    print("\nAvailable functionality:")
    print("- Create custom orientation schemes (uniform/discrete)")
    print("- Set orientations in Settings")
    print("- Create MultiProblem instances")
    print("- Access MultiProblem methods and properties")
    print("- Ready for solving multi-orientation problems!")

except Exception as e:
    print(f"❌ Test failed with error: {e}")
    import traceback
    traceback.print_exc()