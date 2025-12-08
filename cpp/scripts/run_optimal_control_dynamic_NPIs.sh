#!/bin/bash
cd ../../ 

CURRENT_PATH=$(pwd)

"$CURRENT_PATH/cpp/build/bin/optimal_control_dynamic_NPIs" \
    "$CURRENT_PATH/data" \
    "$CURRENT_PATH/cpp/tools/optimal_control/secir_dynamic_NPIs"
