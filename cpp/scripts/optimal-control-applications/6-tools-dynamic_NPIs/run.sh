
#!/bin/bash

# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Set path to the GMGPolar executable
EXECUTABLE="${SCRIPT_DIR}/../../../build/bin/optimal_control_dynamic_NPIs"
# Run executable
"$EXECUTABLE" \
    "${SCRIPT_DIR}/../../../../data" \
    "${SCRIPT_DIR}/../../../../cpp/tools/optimal_control/secir_dynamic_NPIs/optimization_model"
