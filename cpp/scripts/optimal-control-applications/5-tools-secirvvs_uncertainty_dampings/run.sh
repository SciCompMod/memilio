
#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
EXECUTABLE="${SCRIPT_DIR}/../../../build/bin/optimal_control_uncertainty"
# Run executable
"$EXECUTABLE" \
    "${SCRIPT_DIR}/../../../../data" \
    "${SCRIPT_DIR}/../../../../cpp/tools/optimal_control/secirvvs_uncertainty_dampings/optimal_control_settings"
