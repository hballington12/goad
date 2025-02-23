#!/bin/bash
# This script builds the Rust python module and creates the Blender addon zip file.

# Set paths
ROOT_DIR=$(dirname "$(dirname "$(realpath "$0")")")
ADDON_DIR="$ROOT_DIR/blender/addon"
WHEEL_DIR="$ROOT_DIR/target/wheels"

# Ensure clean start
echo "Cleaning up..."
rm -rf "$ADDON_DIR/goad_py"
rm -f "$ROOT_DIR/blender/addon.zip"

# Build Rust module
echo "Building Rust module..."
cd "$ROOT_DIR/goad-py"
./maturin.sh

# Get the most recent wheel file
echo "Extracting wheel file..."
WHEEL_FILE=$(ls -t "$WHEEL_DIR"/goad_py*.whl | head -n1)
TMP_DIR=$(mktemp -d)
mkdir -p "$TMP_DIR"
unzip "$WHEEL_FILE" -d "$TMP_DIR"

# Copy the module to addon directory
cp -r "$TMP_DIR/goad_py" "$ADDON_DIR/"

# Clean up temp directory
rm -rf "$TMP_DIR"

# Create the addon zip
cd "$ROOT_DIR/blender"
zip -r addon.zip addon/

echo "Build complete. Addon zip created at $ROOT_DIR/blender/addon.zip"