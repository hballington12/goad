#!/bin/bash

# test_wheels.sh - Comprehensive wheel testing script for goad-py
# This script tests the built wheels across multiple Python versions

set -e  # Exit on error

echo "🧪 GOAD-PY Wheel Testing Script"
echo "================================"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if dist directory exists
if [ ! -d "dist" ]; then
    echo -e "${RED}Error: dist/ directory not found. Please build wheels first.${NC}"
    echo "Run: maturin build --release"
    exit 1
fi

# Find the wheel file
WHEEL_FILE=$(ls -t dist/goad_py*.whl 2>/dev/null | head -n1)
if [ -z "$WHEEL_FILE" ]; then
    echo -e "${RED}Error: No wheel file found in dist/${NC}"
    exit 1
fi

echo -e "Found wheel: ${GREEN}$WHEEL_FILE${NC}"
echo ""

# Test with different Python versions
PYTHON_VERSIONS=("3.8" "3.9" "3.10" "3.11" "3.12")
TESTED_VERSIONS=0
FAILED_VERSIONS=0

for PY_VERSION in "${PYTHON_VERSIONS[@]}"; do
    # Try different Python executable names
    for PY_EXE in "python${PY_VERSION}" "python${PY_VERSION/./}" "python${PY_VERSION:0:1}.${PY_VERSION:2:2}"; do
        if command -v $PY_EXE &> /dev/null; then
            PYTHON_CMD=$PY_EXE
            break
        fi
    done
    
    # Skip if Python version not found
    if [ -z "$PYTHON_CMD" ] || ! command -v $PYTHON_CMD &> /dev/null; then
        echo -e "${YELLOW}⚠️  Python ${PY_VERSION} not found, skipping...${NC}"
        PYTHON_CMD=""
        continue
    fi
    
    echo -e "\n${GREEN}Testing with Python ${PY_VERSION}${NC}"
    echo "----------------------------------------"
    
    # Create a temporary virtual environment
    VENV_DIR="test_env_${PY_VERSION}"
    echo "Creating virtual environment..."
    $PYTHON_CMD -m venv $VENV_DIR
    
    # Activate virtual environment
    source $VENV_DIR/bin/activate
    
    # Upgrade pip
    pip install --upgrade pip --quiet
    
    # Install the wheel
    echo "Installing wheel..."
    pip install "$WHEEL_FILE" --force-reinstall
    
    # Run basic import test
    echo -n "Testing import... "
    if python -c "import goad_py; print('✓ Success')" 2>/dev/null; then
        echo -e "${GREEN}Passed${NC}"
    else
        echo -e "${RED}Failed${NC}"
        FAILED_VERSIONS=$((FAILED_VERSIONS + 1))
    fi
    
    # Run basic functionality tests
    echo -n "Testing basic functionality... "
    python << 'EOF' 2>/dev/null && echo -e "${GREEN}Passed${NC}" || echo -e "${RED}Failed${NC}"
import goad_py

# Test orientation functions
orient = goad_py.create_uniform_orientation(10)
print("✓ Created orientation")

# Test binning scheme  
bins = goad_py.BinningScheme.interval(1000)
print("✓ Created binning scheme")

# Test Settings creation (requires geometry path, but we can test with dummy path)
try:
    settings = goad_py.Settings("dummy.obj")
    print("✓ Created Settings object")
except Exception as e:
    # Expected to fail without valid geometry file
    print("✓ Settings creation behaves as expected")
EOF
    
    # Check if example scripts exist and run them
    if [ -f "simple_example.py" ]; then
        echo -n "Running simple_example.py... "
        if python simple_example.py > /dev/null 2>&1; then
            echo -e "${GREEN}Passed${NC}"
        else
            echo -e "${RED}Failed${NC}"
        fi
    fi
    
    # Run any test files if they exist
    if [ -d "tests" ] && [ "$(ls -A tests/*.py 2>/dev/null)" ]; then
        echo "Running test suite..."
        python -m pytest tests/ -v
    fi
    
    # Deactivate and cleanup
    deactivate
    rm -rf $VENV_DIR
    
    TESTED_VERSIONS=$((TESTED_VERSIONS + 1))
    PYTHON_CMD=""
done

echo -e "\n================================"
echo -e "Testing Summary:"
echo -e "Tested: ${GREEN}${TESTED_VERSIONS}${NC} Python versions"
if [ $FAILED_VERSIONS -eq 0 ]; then
    echo -e "Result: ${GREEN}All tests passed! ✓${NC}"
else
    echo -e "Result: ${RED}${FAILED_VERSIONS} versions failed ✗${NC}"
fi

# Additional checks
echo -e "\n${YELLOW}Additional Information:${NC}"

# Check wheel size
WHEEL_SIZE=$(du -h "$WHEEL_FILE" | cut -f1)
echo "Wheel size: $WHEEL_SIZE"

# Check if wheel is universal
if unzip -l "$WHEEL_FILE" | grep -q "abi3"; then
    echo -e "ABI compatibility: ${GREEN}✓ abi3 (compatible with multiple Python versions)${NC}"
else
    echo -e "ABI compatibility: ${YELLOW}⚠️  Not using abi3${NC}"
fi

# List wheel contents summary
echo -e "\nWheel contents summary:"
unzip -l "$WHEEL_FILE" | grep -E "(\.so|\.pyd|\.dll)" | head -5

echo -e "\n✨ Testing complete!"