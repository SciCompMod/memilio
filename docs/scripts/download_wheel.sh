#!/bin/bash
# Download pre-built memilio wheels from GitHub CI artifacts
# Uses nightly.link service for public access to GitHub Actions artifacts

set -e

REPO="SciCompMod/memilio"
WORKFLOW="main.yml"
# Use main branch by default, can be overridden via WHEEL_BRANCH env var
BRANCH="${WHEEL_BRANCH:-main}"
PYTHON_VERSION="${1:-312}"  # Default to cp312 (Python 3.12)

# Destination directory
WHEELHOUSE_DIR="${WHEELHOUSE_DIR:-pycode/wheelhouse}"

echo "Downloading memilio wheels from GitHub CI..."
echo "  Repository: $REPO"
echo "  Branch: $BRANCH"
echo "  Python version: cp$PYTHON_VERSION"

# Create wheelhouse directory
mkdir -p "$WHEELHOUSE_DIR"

# Download function
download_artifact() {
    local artifact_name=$1
    local download_url="https://nightly.link/${REPO}/workflows/${WORKFLOW}/${BRANCH}/${artifact_name}.zip"

    echo ""
    echo "Downloading $artifact_name..."
    echo "  URL: $download_url"

    local temp_zip=$(mktemp)
    if curl -L -f -o "$temp_zip" "$download_url"; then
        echo "  Extracting..."
        unzip -o "$temp_zip" -d "$WHEELHOUSE_DIR"
        rm "$temp_zip"
        echo "  Done."
    else
        echo "  Warning: Failed to download $artifact_name (may not exist for this branch)"
        rm -f "$temp_zip"
        return 1
    fi
}

# Download simulation wheel (native extension, version-specific)
download_artifact "python-wheels-simulation"

# Download generation wheel (pure Python, works with any Python 3)
download_artifact "python-wheels-generation" || true

echo ""
echo "Downloaded wheels:"
ls -la "$WHEELHOUSE_DIR"/*.whl 2>/dev/null || echo "No wheels found!"

# Check for simulation wheel
SIMULATION_WHEEL=$(ls "$WHEELHOUSE_DIR"/*simulation*cp${PYTHON_VERSION}*.whl 2>/dev/null | head -1)
if [ -n "$SIMULATION_WHEEL" ]; then
    echo ""
    echo "Found simulation wheel for Python $PYTHON_VERSION: $SIMULATION_WHEEL"
else
    echo ""
    echo "Warning: No simulation wheel found for cp$PYTHON_VERSION"
    echo "Available wheels:"
    ls "$WHEELHOUSE_DIR"/*.whl 2>/dev/null || echo "None"
fi

# Check for generation wheel (may be pure Python py3-*.whl or version-specific cp*-*.whl)
GENERATION_WHEEL=$(ls "$WHEELHOUSE_DIR"/*generation*cp${PYTHON_VERSION}*.whl 2>/dev/null | head -1)
if [ -z "$GENERATION_WHEEL" ]; then
    # Fallback to pure Python wheel
    GENERATION_WHEEL=$(ls "$WHEELHOUSE_DIR"/*generation*py3*.whl 2>/dev/null | head -1)
fi
if [ -n "$GENERATION_WHEEL" ]; then
    echo "Found generation wheel: $GENERATION_WHEEL"
fi

echo ""
echo "Wheel download complete!"
