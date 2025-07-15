#!/bin/zsh

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
SIM_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/sim_viz.py"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
PYTHON3_DIR="$MAIN_PATH/m_v/bin/python3"

# Simulation parameters
EVENT_TYPE="restaurant_table_equals_household"
SIM_TYPE="memilio"
VIZ_OPTIONS="--s90percentile"

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local results_dir="$RESULTS_BASE_DIR/run_$timestamp"
    echo "$results_dir"
}

run_simulation() {
    local results_dir=$1
    echo "Running simulation..."
    $MAIN_EXECUTABLE --event "$EVENT_TYPE" --sim "$SIM_TYPE" --output_dir "$results_dir"
    return $?
}

run_visualization() {
    local results_dir=$1
    echo "Visualizing results..."
    $PYTHON3_DIR $SIM_VIZ_SCRIPT \
        --path-to-infection-states "$results_dir/infection_state_per_age_group" \
        --path-to-loc-types "$results_dir/infection_per_location_type_per_age_group" \
        $VIZ_OPTIONS
    return $?
}

#copy results to last_result folder
copy_results_to_last_result() { 
    local results_dir=$1
    local parent_dir=$(dirname "$results_dir")
    local last_result_dir="$parent_dir/last_result"
    echo "Copying results to last_result directory: $last_result_dir"
    if [ -d "$last_result_dir" ]; then
        rm -rf "$last_result_dir"
    fi
    mkdir -p "$last_result_dir"
    cp -r "$results_dir"/* "$last_result_dir"
    echo "Results copied to $last_result_dir"
}

# Main execution
main() {
    local results_dir=$(create_results_dir)
    
    if ! run_simulation "$results_dir"; then
        echo "Error: Simulation failed."
        exit 1
    fi
    echo "Simulation completed. Results saved in $results_dir."
    
    if ! run_visualization "$results_dir"; then
        echo "Error: Visualization failed."
        exit 1
    fi

    if [ -d "$results_dir/infection_state_per_age_group" ]; then
        copy_results_to_last_result "$results_dir"
    else
        echo "No infection state results found to copy."
    fi

    echo "Visualization completed successfully."
}

# Run main function
main "$@"