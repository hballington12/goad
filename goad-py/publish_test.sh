#!/bin/bash

# publish_test.sh - Publish to TestPyPI for testing
# Run this script to test publishing before doing a real release

set -e

echo "🧪 Publishing GOAD-PY to TestPyPI"
echo "=================================="

# Check if we're in the right directory
if [ ! -f "pyproject.toml" ]; then
    echo "❌ Error: pyproject.toml not found. Run this from the goad-py directory."
    exit 1
fi

# Check if maturin is installed
if ! command -v maturin &> /dev/null; then
    echo "📦 Installing maturin..."
    pip install maturin
fi

# Clean any existing builds
echo "🧹 Cleaning previous builds..."
rm -rf dist/ target/wheels/

# Build wheels
echo "🔨 Building wheels..."
maturin build --release

# Create dist directory and copy wheels
mkdir -p dist/
cp ../target/wheels/*.whl dist/

# Show what we're about to upload
echo ""
echo "📦 Built packages:"
ls -la dist/

echo ""
echo "🚀 Uploading to TestPyPI..."
echo "Note: You'll need TestPyPI credentials configured:"
echo "  - Create account at https://test.pypi.org"
echo "  - Generate API token"
echo "  - Run: pip install twine"
echo "  - Configure: twine configure (or use --username __token__ --password <token>)"
echo ""

read -p "Continue with upload? (y/N): " confirm
if [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]; then
    # Upload to TestPyPI
    if command -v twine &> /dev/null; then
        twine upload --repository testpypi dist/*
    else
        echo "Installing twine..."
        pip install twine
        twine upload --repository testpypi dist/*
    fi
    
    echo ""
    echo "✅ Upload complete!"
    echo "🔗 View at: https://test.pypi.org/project/goad-py/"
    echo ""
    echo "🧪 Test installation with:"
    echo "pip install --index-url https://test.pypi.org/simple/ goad-py"
else
    echo "Upload cancelled."
fi

echo ""
echo "📝 Next steps for real release:"
echo "1. Test the TestPyPI package works correctly"
echo "2. Create and push a git tag: git tag v0.2.0 && git push origin v0.2.0"
echo "3. GitHub Actions will automatically publish to PyPI"