#!/bin/bash

#!/bin/bash

# compile.sh - Script to build the project in Debug or Release mode
#
# Usage:
#   1. Navigate to the script directory:
#        cd cpp/scripts/
#   2. Make the script executable (if not already):
#        chmod +x compile.sh
#   3. Compile in Release mode:
#        ./compile.sh Release
#   4. Compile in Debug mode:
#        ./compile.sh Debug
#   5. Compile using the current build mode:
#        ./compile.sh

# Function to display usage information
usage() {
    echo "Usage: $0 [Debug|Release]"
    exit 1
}

# Define CMake configuration options
CMAKE_OPTIONS="
    -DMEMILIO_BUILD_TESTS=OFF
    -DMEMILIO_BUILD_EXAMPLES=ON
    -DMEMILIO_BUILD_MODELS=ON
    -DMEMILIO_BUILD_SIMULATIONS=OFF
    -DMEMILIO_BUILD_BENCHMARKS=OFF
    -DMEMILIO_USE_BUNDLED_SPDLOG=ON
    -DMEMILIO_USE_BUNDLED_EIGEN=ON
    -DMEMILIO_USE_BUNDLED_BOOST=ON
    -DMEMILIO_USE_BUNDLED_JSONCPP=ON
    -DMEMILIO_SANITIZE_ADDRESS=OFF
    -DMEMILIO_SANITIZE_UNDEFINED=OFF
    -DMEMILIO_BUILD_SHARED_LIBS=ON
    -DMEMILIO_BUILD_STATIC_LIBS=ON
    -DMEMILIO_ENABLE_MPI=ON
    -DMEMILIO_ENABLE_OPENMP=ON
    -DMEMILIO_ENABLE_WARNINGS=ON
    -DMEMILIO_ENABLE_WARNINGS_AS_ERRORS=OFF
    -DMEMILIO_ENABLE_IPOPT=ON
    -DMEMILIO_ENABLE_PROFILING=OFF
    -DMEMILIO_ENABLE_LIKWID_MARKER=OFF
"

# Check if build directory exists in the parent directory
if [ -d "../build" ]; then
    build_exists=true
else
    build_exists=false
fi

# Check if build directory exists and delete if it does (only if an argument is provided)
if [ -n "$1" ] && [ -d "../build" ]; then
    echo "Removing existing build directory..."
    rm -rf ../build/CMake* || { echo "Failed to remove build directory"; exit 1; }
    build_exists=false
fi

# Create build directory in the parent directory if it doesn't exist and if a build folder didn't exist before
if ! $build_exists; then
    echo "Creating build directory..."
    mkdir -p ../build || { echo "Failed to create build directory"; exit 1; }
    build_exists=true
fi

# Determine build type
if [ -n "$1" ]; then
    case "$1" in
        Debug)
            build_type="Debug"
            ;;
        Release)
            build_type="Release"
            ;;
        *)
            echo "Invalid build type. Please specify Debug or Release."
            usage
            ;;
    esac
else
    if [ "$build_exists" != true ]; then
        echo "No build type specified. Please provide either Debug or Release."
        usage
    fi
fi

# Run cmake configuration
if [ -n "$build_type" ]; then
    echo "Configuring with $build_type build type..."

    cmake -S "${PWD}/.." -B "${PWD}/../build" \
        -DCMAKE_BUILD_TYPE="$build_type" \
        $CMAKE_OPTIONS \
        || { echo "CMake configuration failed"; exit 1; }
fi

# Build the project
cmake --build ../build -j 16
