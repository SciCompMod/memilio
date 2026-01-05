#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
EXECUTABLE="${SCRIPT_DIR}/../../../build/bin/seair_pc_global_max"
# Run executable
"$EXECUTABLE"
