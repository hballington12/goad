#!/bin/bash

echo "🔧 Building and testing GOAD Python bindings..."

# Detect Python version
PYTHON_VERSION=$(python3 --version 2>&1 | grep -oP '3\.\d+')
echo "📍 Detected Python version: $PYTHON_VERSION"

# Check if virtual environment exists
if [ ! -d ".venv" ]; then
    echo "📦 Creating virtual environment..."
    python3 -m venv .venv
fi

# Activate virtual environment
echo "🔌 Activating virtual environment..."
source .venv/bin/activate

# Upgrade pip and install maturin
echo "⬆️ Upgrading pip and installing maturin..."
pip install --upgrade pip
pip install --upgrade maturin

# Build the Rust project first
echo "🦀 Building Rust project..."
cd ..
cargo build --release
cd goad-py

# Build Python wheels
echo "🎯 Building Python wheels..."
# Use detected Python version for wheel building
maturin build --release -i python$PYTHON_VERSION

# Install the wheel
echo "📦 Installing wheel..."
WHEEL_FILE=$(ls -t ../target/wheels/goad_py*.whl | head -n1)
echo "Installing: $WHEEL_FILE"
pip install "$WHEEL_FILE" --force-reinstall

# Test the installation
echo "🧪 Testing installation..."
python3 -c "import goad_py as goad; print('✅ GOAD Python bindings imported successfully!')"

echo "🎉 Build and installation complete!"
echo "🐍 To run tests: source .venv/bin/activate && python test_multiproblem.py"