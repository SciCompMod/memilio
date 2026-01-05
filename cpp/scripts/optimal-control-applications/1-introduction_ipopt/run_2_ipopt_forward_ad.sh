#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
EXECUTABLE="${SCRIPT_DIR}/../../../build/bin/ipopt_forward_ad"
# Run executable
"$EXECUTABLE"
