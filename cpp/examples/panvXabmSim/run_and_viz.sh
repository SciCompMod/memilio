#!/bin/zsh

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
SIM_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/sim_viz.py"
COMPARISON_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/infection_cumulated_comp.py"
CONTACT_NETWORK_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/contact_network.py"
INFECTION_TIMELINE_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/infection_timeline.py"
STAT_TABLE_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/stat_table.py"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
VIZ_OUTPUT_DIR="$MAIN_PATH/examples/panvXabmSim/results/results_viz"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"

# Simulation parameters
EVENT_TYPE="restaurant_table_equals_half_household"  # Options: restaurant_table_equals_household restaurant_table_equals_half_household restaurant_table_equals_random_household work_meeting_many work_meeting_baseline
VIZ_OPTIONS="--s90percentile"
NUM_DAYS=10
NUM_PERSONS=1000
RUNS=100

# BOOL for visualization
VISUALIZE=true

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local sim_type=$1
    local results_dir="$RESULTS_BASE_DIR/run_${timestamp}_${sim_type}"
    echo "$results_dir"
}

create_viz_output_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local viz_dir="$VIZ_OUTPUT_DIR/comparison_$timestamp"
    mkdir -p "$viz_dir"
    echo "$viz_dir"
}

run_simulation() {
    local results_dir=$1
    local sim_type=$2
    echo "Running $sim_type simulation..."
    $MAIN_EXECUTABLE --event "$EVENT_TYPE" --sim "$sim_type" --output_dir "$results_dir" \
        --days "$NUM_DAYS" --n_persons "$NUM_PERSONS" --runs "$RUNS"
    return $?
}

run_individual_visualization() {
    local results_dir=$1
    local sim_type=$2
    echo "Visualizing individual $sim_type results..."
    $PYTHON3_DIR $SIM_VIZ_SCRIPT \
        --path-to-infection-states "$results_dir/infection_state_per_age_group" \
        --path-to-loc-types "$results_dir/infection_per_location_type_per_age_group" \
        --label1 "$sim_type" \
        $VIZ_OPTIONS
    return $?
}

run_comparison_visualizations() {
    local memilio_dir=$1
    local panvadere_dir=$2
    local viz_output_dir=$3
    
    echo "Creating comparison visualizations in $viz_output_dir..."
    
    # 1. Main cumulative infections comparison
    echo "  -> Running cumulative infections comparison..."
    $PYTHON3_DIR $COMPARISON_VIZ_SCRIPT \
        --path-to-memilio-sim "$memilio_dir/amount_of_infections" \
        --path-to-panvXabmSim "$panvadere_dir/amount_of_infections" \
        --output-path "$viz_output_dir" \
        $VIZ_OPTIONS

    # 2. Side-by-side infection states comparison
    echo "  -> Running infection states comparison..."
    $PYTHON3_DIR $SIM_VIZ_SCRIPT \
        --path-to-infection-states "$memilio_dir/infection_state_per_age_group" \
        --path-to-infection-states-2 "$panvadere_dir/infection_state_per_age_group" \
        --label1 "Memilio" \
        --label2 "Panvadere" \
        --output-path "$viz_output_dir" \
        $VIZ_OPTIONS
    
    # 3. Contact network analysis (if CSV files exist)
    # if [ -f "$memilio_dir/contact_intensiveness.csv" ] && [ -f "$memilio_dir/infection_count.csv" ]; then
    #     echo "  -> Running contact network analysis..."
    #     $PYTHON3_DIR $CONTACT_NETWORK_SCRIPT \
    #         --data-dir "$memilio_dir" \
    #         --output-path "$viz_output_dir/contact_network_memilio.png" \
    #         --scenario-name "Memilio"
    # fi
    
    # if [ -f "$panvadere_dir/contact_intensiveness.csv" ] && [ -f "$panvadere_dir/infection_count.csv" ]; then
    #     echo "  -> Running contact network analysis for Panvadere..."
    #     $PYTHON3_DIR $CONTACT_NETWORK_SCRIPT \
    #         --data-dir "$panvadere_dir" \
    #         --output-path "$viz_output_dir/contact_network_panvadere.png" \
    #         --scenario-name "Panvadere"
    # fi
    
    # # 4. Infection timeline analysis (if detailed infection files exist)
    # if [ -f "$memilio_dir/best_run_detailed_infection.csv" ] && [ -f "$memilio_dir/best_run_contact_data.csv" ]; then
    #     echo "  -> Running infection timeline analysis for Memilio..."
    #     $PYTHON3_DIR $INFECTION_TIMELINE_SCRIPT \
    #         --contact-file "$memilio_dir/best_run_contact_data.csv" \
    #         --infection-file "$memilio_dir/best_run_detailed_infection.csv" \
    #         --output-path "$viz_output_dir/infection_timeline_memilio.png" \
    #         --scenario-name "Memilio"
    # fi
    
    # if [ -f "$panvadere_dir/best_run_detailed_infection.csv" ] && [ -f "$panvadere_dir/best_run_contact_data.csv" ]; then
    #     echo "  -> Running infection timeline analysis for Panvadere..."
    #     $PYTHON3_DIR $INFECTION_TIMELINE_SCRIPT \
    #         --contact-file "$panvadere_dir/best_run_contact_data.csv" \
    #         --infection-file "$panvadere_dir/best_run_detailed_infection.csv" \
    #         --output-path "$viz_output_dir/infection_timeline_panvadere.png" \
    #         --scenario-name "Panvadere"
    # fi
    
    # # 5. Statistical analysis tables
    # echo "  -> Running statistical analysis..."
    # if [ -d "$memilio_dir" ]; then
    #     $PYTHON3_DIR $STAT_TABLE_SCRIPT \
    #         --input_dir "$memilio_dir" \
    #         --scenario_name "Memilio"
    # fi
    
    # if [ -d "$panvadere_dir" ]; then
    #     $PYTHON3_DIR $STAT_TABLE_SCRIPT \
    #         --input_dir "$panvadere_dir" \
    #         --scenario_name "Panvadere"
    # fi
    
    # # Copy statistical reports to viz output directory
    # if [ -f "$memilio_dir/epidemic_metrics.txt" ]; then
    #     cp "$memilio_dir/epidemic_metrics.txt" "$viz_output_dir/epidemic_metrics_memilio.txt"
    # fi
    # if [ -f "$panvadere_dir/epidemic_metrics.txt" ]; then
    #     cp "$panvadere_dir/epidemic_metrics.txt" "$viz_output_dir/epidemic_metrics_panvadere.txt"
    # fi
    
    return $?
}

copy_results_to_last_result() { 
    local results_dir=$1
    local sim_type=$2
    local parent_dir=$(dirname "$results_dir")
    local last_result_dir="$parent_dir/last_result_${sim_type}"
    echo "Copying $sim_type results to last_result directory: $last_result_dir"
    if [ -d "$last_result_dir" ]; then
        rm -rf "$last_result_dir"
    fi
    mkdir -p "$last_result_dir"
    cp -r "$results_dir"/* "$last_result_dir"
    echo "Results copied to $last_result_dir"
}

# Main execution
main() {
    # Get timestamp for this run
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    # Create results directories for both simulations
    local memilio_dir="$RESULTS_BASE_DIR/run_${timestamp}_memilio"
    local panvadere_dir="$RESULTS_BASE_DIR/run_${timestamp}_panvadere"
    
    echo "=== Running Dual Simulation Comparison ==="
    echo "Event Type: $EVENT_TYPE"
    echo "Days: $NUM_DAYS, Persons: $NUM_PERSONS"
    echo "Memilio results: $memilio_dir"
    echo "Panvadere results: $panvadere_dir"
    echo ""
    
    # Run Memilio simulation
    echo "=== STEP 1: Running Memilio Simulation ==="
    if ! run_simulation "$memilio_dir" "memilio"; then
        echo "Error: Memilio simulation failed."
        exit 1
    fi
    echo "Memilio simulation completed. Results saved in $memilio_dir"
    echo ""
    
    # Run Panvadere simulation  
    echo "=== STEP 2: Running Panvadere Simulation ==="
    if ! run_simulation "$panvadere_dir" "panvadere"; then
        echo "Error: Panvadere simulation failed."
        exit 1
    fi
    echo "Panvadere simulation completed. Results saved in $panvadere_dir"
    echo ""
    
    if [ "$VISUALIZE" = true ]; then
        echo "=== STEP 3: Creating Visualizations ==="
        
        # Create visualization output directory
        local viz_output_dir=$(create_viz_output_dir)
        echo "Visualization output directory: $viz_output_dir"
        
        # Run individual visualizations
        echo "Creating individual visualizations..."
        run_individual_visualization "$memilio_dir" "Memilio"
        run_individual_visualization "$panvadere_dir" "Panvadere"
        
        # Run comparison visualizations
        echo "Creating comparison visualizations..."
        if ! run_comparison_visualizations "$memilio_dir" "$panvadere_dir" "$viz_output_dir"; then
            echo "Warning: Some comparison visualizations may have failed."
        fi
        
        echo "All visualizations completed. Results saved in $viz_output_dir"
        # echo ""
    fi
    
    # Copy results to last_result folders
    echo "=== STEP 4: Copying Results ==="
    if [ -d "$memilio_dir/infection_state_per_age_group" ]; then
        copy_results_to_last_result "$memilio_dir" "memilio"
    else
        echo "No Memilio infection state results found to copy."
    fi
    
    if [ -d "$panvadere_dir/infection_state_per_age_group" ]; then
        copy_results_to_last_result "$panvadere_dir" "panvadere"
    else
        echo "No Panvadere infection state results found to copy."
    fi
    
    echo ""
    echo "=== SIMULATION COMPARISON COMPLETED SUCCESSFULLY ==="
    echo "Memilio results: $memilio_dir"
    echo "Panvadere results: $panvadere_dir"
    if [ "$VISUALIZE" = true ]; then
        echo "Comparison visualizations: $viz_output_dir"
    fi
}

# Run main function
main "$@"