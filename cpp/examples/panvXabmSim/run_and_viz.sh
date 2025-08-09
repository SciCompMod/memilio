#!/bin/zsh

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
SIM_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/sim_viz.py"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"

# Simulation parameters
# EVENT_TYPE="restaurant_table_equals_random_household"
# EVENT_TYPE="restaurant_table_equals_half_household"  # Options: restaurant_table_equals_half_household, restaurant_table_equals_full_household
EVENT_TYPE="restaurant_table_equals_household"  # Options: restaurant_table_equals_half
# EVENT_TYPE="work_meeting_many"
# EVENT_TYPE="work_meeting_few"  # Options: work_meeting_many, work_meeting_few, work_meeting_low
# EVENT_TYPE="work_meeting_low"  # Options: work_meeting_many, work
SIM_TYPE="memilio"  # Options: panvadere, memilio
# SIM_TYPE="panvadere"  # Options: panvadere, memilio
VIZ_OPTIONS="--s90percentile"
NUM_DAYS=10
NUM_PERSONS=3000
# Set environment variables
export OMP_NUM_THREADS=$OMP_NUM_THREADS

# BOOL for visualization
VISUALIZE=true

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local results_dir="$RESULTS_BASE_DIR/run_$timestamp"
    echo "$results_dir"
}

run_simulation() {
    local results_dir=$1
    echo "Running simulation..."
    $MAIN_EXECUTABLE --event "$EVENT_TYPE" --sim "$SIM_TYPE" --output_dir "$results_dir"    \
        --days "$NUM_DAYS" --n_persons "$NUM_PERSONS"  
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
    
    if [ "$VISUALIZE" = true ]; then
        echo "Visualization is enabled."
        if ! run_visualization "$results_dir"; then
            echo "Error: Visualization failed."
            exit 1
        fi
        echo "Visualization completed successfully."
    fi
   

    if [ -d "$results_dir/infection_state_per_age_group" ]; then
        copy_results_to_last_result "$results_dir"
    else
        echo "No infection state results found to copy."
    fi

    echo "All operations completed successfully."
  
}

# Run main function
main "$@"