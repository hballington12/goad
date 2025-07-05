#!/bin/bash

# build_wheels_local.sh - Build wheels locally using cibuildwheel
# This mimics the CI build process for local testing

set -e

echo "🎯 Building GOAD-PY wheels locally with cibuildwheel"
echo "===================================================="

# Check if cibuildwheel is installed
if ! command -v cibuildwheel &> /dev/null; then
    echo "📦 Installing cibuildwheel..."
    pip install cibuildwheel
fi

# Detect platform
PLATFORM=""
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="macos"
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "win32" ]]; then
    PLATFORM="windows"
else
    echo "❌ Unsupported platform: $OSTYPE"
    exit 1
fi

echo "🖥️  Detected platform: $PLATFORM"

# Clean previous builds
echo "🧹 Cleaning previous builds..."
rm -rf dist/ wheelhouse/ build/

# Build the parent Rust project first
echo "🦀 Building parent Rust project..."
cd ..
cargo build --release
cd goad-py

# Set environment variables for cibuildwheel
export CIBW_BUILD="cp38-* cp39-* cp310-* cp311-* cp312-*"
export CIBW_SKIP="*-musllinux* *-win32 *-manylinux_i686"
export CIBW_ENVIRONMENT="MATURIN_PEP517_ARGS='--release'"

# Platform-specific settings
if [ "$PLATFORM" == "linux" ]; then
    export CIBW_MANYLINUX_X86_64_IMAGE="manylinux2014"
    export CIBW_MANYLINUX_AARCH64_IMAGE="manylinux2014"
fi

# Run cibuildwheel
echo "🔨 Building wheels with cibuildwheel..."
cibuildwheel --platform $PLATFORM

# Move wheels to dist/
mkdir -p dist/
mv wheelhouse/* dist/ 2>/dev/null || true

# List built wheels
echo ""
echo "✅ Built wheels:"
ls -la dist/*.whl 2>/dev/null || echo "No wheels found in dist/"

# Run the test script if it exists
if [ -f "test_wheels.sh" ]; then
    echo ""
    echo "🧪 Running wheel tests..."
    ./test_wheels.sh
fi

echo ""
echo "🎉 Build complete! Wheels are in the dist/ directory."
echo "To test a specific wheel: pip install dist/goad_py-*.whl"