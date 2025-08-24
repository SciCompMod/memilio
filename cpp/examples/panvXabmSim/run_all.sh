# Helper function to find result path by key
find_result_path() {
    local search_key=$1
    >&2 echo "    DEBUG: Searching for key '$search_key' in array of ${#RESULTS_KEYS[@]} elements"
    
    for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
        local current_key="${RESULTS_KEYS[i]}"
        >&2 echo "    DEBUG: Checking [$i] '$current_key'"
        
        if [[ "$current_key" == "$search_key" ]]; then
            local found_path="${RESULTS_PATHS[i]}"
            >&2 echo "    DEBUG: FOUND at index $i: '$found_path'"
            echo "$found_path"  # This goes to stdout (return value)
            return 0
        fi
    done
    
    >&2 echo "    DEBUG: NOT FOUND - '$search_key' not in array"
    echo ""  # Empty result to stdout
    return 1
}

# =============================================================================
# COMPREHENSIVE SIMULATION AND VISUALIZATION SCRIPT
# =============================================================================
# This script runs all combinations of event types and simulation types,
# then creates comprehensive comparison visualizations.
#
# VISUALIZATION CONTROL:
# Set any of these to 'false' to disable specific visualization types:
# - ENABLE_INDIVIDUAL_VIZ: Per-simulation visualizations during Step 1
# - ENABLE_PAIRWISE_VIZ: Memilio vs Panvadere comparisons in Step 2  
# - ENABLE_COMPREHENSIVE_VIZ: Multi-panel overview in Step 3
# - ENABLE_CONTACT_NETWORK_VIZ: Contact network analysis (part of pairwise)
# - ENABLE_INFECTION_TIMELINE_VIZ: Infection timeline analysis (part of pairwise)
# - ENABLE_STATISTICAL_ANALYSIS: Statistical tables and reports (part of pairwise)
# =============================================================================

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
COMPARISON_PYTHON_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/compare_all_scenarios.py"

# All available event types (based on your actual system)
EVENT_TYPES=(
    "restaurant_table_equals_household"
    "restaurant_table_equals_half_household" 
    "work_meeting_many"
    "work_meeting_baseline"
)

# Both simulation types
SIM_TYPES=("memilio" "panvadere")

# Simulation parameters
VIZ_OPTIONS="--s90percentile"
NUM_DAYS=10
NUM_PERSONS=1000
RUNS=10
INFECTION_K=22.6

# Visualization control flags
ENABLE_INDIVIDUAL_VIZ=true          # Individual visualization for each simulation
ENABLE_PAIRWISE_VIZ=true            # Pairwise comparisons (memilio vs panvadere)
ENABLE_COMPREHENSIVE_VIZ=true       # Comprehensive multi-panel comparison
ENABLE_CONTACT_NETWORK_VIZ=false     # Contact network analysis
ENABLE_INFECTION_TIMELINE_VIZ=true  # Infection timeline analysis  
ENABLE_STATISTICAL_ANALYSIS=false    # Statistical tables and reports

# Results storage - using simple arrays with consistent indexing
# Initialize as truly empty arrays (no pre-filled elements)
declare -a RESULTS_KEYS=()
declare -a RESULTS_PATHS=()

# Validate arrays are clean at startup
if [ ${#RESULTS_KEYS[@]} -ne 0 ] || [ ${#RESULTS_PATHS[@]} -ne 0 ]; then
    echo "WARNING: Arrays not empty at startup. Cleaning..."
    declare -a RESULTS_KEYS=()
    declare -a RESULTS_PATHS=()
fi

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local sim_type=$1
    local event_type=$2
    local results_dir="$RESULTS_BASE_DIR/run_${timestamp}_${event_type}_${sim_type}"
    echo "$results_dir"
}

create_viz_output_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local viz_dir="$VIZ_OUTPUT_DIR/comprehensive_analysis_$timestamp"
    mkdir -p "$viz_dir"
    echo "$viz_dir"
}

run_simulation() {
    local results_dir=$1
    local sim_type=$2
    local event_type=$3
    echo "  -> Running $sim_type simulation with $event_type event..."
    $MAIN_EXECUTABLE --event "$event_type" --sim "$sim_type" --output_dir "$results_dir" \
        --days "$NUM_DAYS" --n_persons "$NUM_PERSONS" --runs "$RUNS" --infection_k "$INFECTION_K"
    return $?
}

run_individual_visualization() {
    local results_dir=$1
    local sim_type=$2
    local event_type=$3
    
    if [ "$ENABLE_INDIVIDUAL_VIZ" = false ]; then
        echo "    -> Skipping individual visualization (disabled)"
        return 0
    fi
    
    echo "    -> Creating individual visualization..."
    $PYTHON3_DIR $SIM_VIZ_SCRIPT \
        --path-to-infection-states "$results_dir/infection_state_per_age_group" \
        --path-to-loc-types "$results_dir/infection_per_location_type_per_age_group" \
        --label1 "${sim_type}_${event_type}" \
        $VIZ_OPTIONS
    return $?
}

run_pairwise_comparison() {
    local memilio_dir=$1
    local panvadere_dir=$2
    local event_type=$3
    local viz_output_dir=$4
    
    if [ "$ENABLE_PAIRWISE_VIZ" = false ]; then
        echo "  -> Skipping pairwise comparison visualizations (disabled)"
        return 0
    fi
    
    echo "  -> Creating pairwise comparison visualizations for $event_type..."
    
    local output_subdir="$viz_output_dir/${event_type}"
    mkdir -p "$output_subdir"
    
    # 1. Main cumulative infections comparison
    echo "    -> Running cumulative infections comparison for $event_type..."
    $PYTHON3_DIR $COMPARISON_VIZ_SCRIPT \
        --path-to-memilio-sim "$memilio_dir/amount_of_infections" \
        --path-to-panvXabmSim "$panvadere_dir/amount_of_infections" \
        --output-path "$output_subdir" \
        $VIZ_OPTIONS

    # 2. Side-by-side infection states comparison
    echo "    -> Running infection states comparison for $event_type..."
    $PYTHON3_DIR $SIM_VIZ_SCRIPT \
        --path-to-infection-states "$memilio_dir/infection_state_per_age_group" \
        --path-to-infection-states-2 "$panvadere_dir/infection_state_per_age_group" \
        --label1 "Memilio" \
        --label2 "Panvadere" \
        --output-path "$output_subdir" \
        $VIZ_OPTIONS
    
    # 3. Contact network analysis (if enabled and CSV files exist)
    if [ "$ENABLE_CONTACT_NETWORK_VIZ" = true ]; then
        if [ -f "$memilio_dir/contact_intensiveness.csv" ] && [ -f "$memilio_dir/infection_count.csv" ]; then
            echo "    -> Running contact network analysis for Memilio ($event_type)..."
            $PYTHON3_DIR $CONTACT_NETWORK_SCRIPT \
                --data-dir "$memilio_dir" \
                --output-path "$output_subdir/contact_network_memilio.png" \
                --scenario-name "Memilio_${event_type}" \
                --total-runs "$RUNS" \
                --max-persons 1000
        fi
        
        if [ -f "$panvadere_dir/contact_intensiveness.csv" ] && [ -f "$panvadere_dir/infection_count.csv" ]; then
            echo "    -> Running contact network analysis for Panvadere ($event_type)..."
            $PYTHON3_DIR $CONTACT_NETWORK_SCRIPT \
                --data-dir "$panvadere_dir" \
                --output-path "$output_subdir/contact_network_panvadere.png" \
                --scenario-name "Panvadere_${event_type}" \
                --total-runs "$RUNS" \
                --max-persons 1000
        fi
    else
        echo "    -> Skipping contact network analysis (disabled)"
    fi
    
    # 4. Infection timeline analysis (if enabled and detailed infection files exist)
    if [ "$ENABLE_INFECTION_TIMELINE_VIZ" = true ]; then
        if [ -f "$memilio_dir/best_run_detailed_infection.csv" ] && [ -f "$memilio_dir/best_run_contact_data.csv" ]; then
            echo "    -> Running infection timeline analysis for Memilio ($event_type)..."
            $PYTHON3_DIR $INFECTION_TIMELINE_SCRIPT \
                --contact-file "$memilio_dir/best_run_contact_data.csv" \
                --infection-file "$memilio_dir/best_run_detailed_infection.csv" \
                --output-path "$output_subdir/infection_timeline_memilio.png" \
                --scenario-name "Memilio_${event_type}" \
                --show-all-potential-infectors \
                --max-display 30
        fi
        
        if [ -f "$panvadere_dir/best_run_detailed_infection.csv" ] && [ -f "$panvadere_dir/best_run_contact_data.csv" ]; then
            echo "    -> Running infection timeline analysis for Panvadere ($event_type)..."
            $PYTHON3_DIR $INFECTION_TIMELINE_SCRIPT \
                --contact-file "$panvadere_dir/best_run_contact_data.csv" \
                --infection-file "$panvadere_dir/best_run_detailed_infection.csv" \
                --output-path "$output_subdir/infection_timeline_panvadere.png" \
                --scenario-name "Panvadere_${event_type}" \
                --show-all-potential-infectors \
                --max-display 30
        fi
    else
        echo "    -> Skipping infection timeline analysis (disabled)"
    fi
    
    # 5. Statistical analysis tables (if enabled)
    if [ "$ENABLE_STATISTICAL_ANALYSIS" = true ]; then
        echo "    -> Running statistical analysis for $event_type..."
        if [ -d "$memilio_dir" ]; then
            $PYTHON3_DIR $STAT_TABLE_SCRIPT \
                --input_dir "$memilio_dir" \
                --scenario_name "Memilio_${event_type}" \
                --output_dir "$output_subdir"
        fi
        
        if [ -d "$panvadere_dir" ]; then
            $PYTHON3_DIR $STAT_TABLE_SCRIPT \
                --input_dir "$panvadere_dir" \
                --scenario_name "Panvadere_${event_type}" \
                --output_dir "$output_subdir"
        fi
        
        # Copy statistical reports to viz output directory
        if [ -f "$memilio_dir/epidemic_metrics.txt" ]; then
            cp "$memilio_dir/epidemic_metrics.txt" "$output_subdir/epidemic_metrics_memilio.txt"
        fi
        if [ -f "$panvadere_dir/epidemic_metrics.txt" ]; then
            cp "$panvadere_dir/epidemic_metrics.txt" "$output_subdir/epidemic_metrics_panvadere.txt"
        fi
    else
        echo "    -> Skipping statistical analysis (disabled)"
    fi
    
    return $?
}

copy_results_to_last_result() { 
    local results_dir=$1
    local sim_type=$2
    local event_type=$3
    local parent_dir=$(dirname "$results_dir")
    local last_result_dir="$parent_dir/last_result_${event_type}_${sim_type}"
    echo "  -> Copying $sim_type results for $event_type to: $last_result_dir"
    if [ -d "$last_result_dir" ]; then
        rm -rf "$last_result_dir"
    fi
    mkdir -p "$last_result_dir"
    cp -r "$results_dir"/* "$last_result_dir"
}

# Main execution
main() {
    # Get timestamp for this run
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    # Create main visualization output directory
    local viz_output_dir=$(create_viz_output_dir)
    echo "=== Running Comprehensive Simulation Analysis ==="
    echo "Timestamp: $timestamp"
    echo "Event Types: ${EVENT_TYPES[*]}"
    echo "Simulation Types: ${SIM_TYPES[*]}"
    echo "Days: $NUM_DAYS, Persons: $NUM_PERSONS, Runs: $RUNS"
    echo "Main visualization output: $viz_output_dir"
    echo ""
    
    # Run all combinations of event types and simulation types
    echo "=== STEP 1: Running All Simulations ==="
    for event_type in "${EVENT_TYPES[@]}"; do
        echo ""
        echo "Processing event type: $event_type"
        
        for sim_type in "${SIM_TYPES[@]}"; do            
            # Create results directory for this combination
            local results_dir="$RESULTS_BASE_DIR/run_${timestamp}_${event_type}_${sim_type}"
            
            # Run simulation
            if ! run_simulation "$results_dir" "$sim_type" "$event_type"; then
                echo "    ERROR: $sim_type simulation failed for $event_type"
                echo "    -> Simulation exit code: $?"
                echo "    -> This result will NOT be stored in arrays"
                continue
            fi
            
            echo "    SUCCESS: Results saved in $results_dir"
            
            # Verify the result directory actually contains expected output
            if [ ! -d "$results_dir" ]; then
                echo "    ERROR: Results directory not found: $results_dir"
                continue
            fi
            
            # Store the results using parallel arrays
            local key="${sim_type}_${event_type}"
            
            # Double-check we're not adding empty values
            if [[ -z "$key" || -z "$results_dir" ]]; then
                echo "    ERROR: Attempted to store empty key or path"
                echo "    -> key: '$key'"
                echo "    -> results_dir: '$results_dir'"
                continue
            fi
            
            # Debug: Show what we're about to add
            echo "    -> About to add to arrays:"
            echo "       key: '$key'"
            echo "       path: '$results_dir'"
            echo "       current array length: ${#RESULTS_KEYS[@]}"
            
            # Add to arrays - use explicit indexing to ensure proper assignment
            RESULTS_KEYS[${#RESULTS_KEYS[@]}]="$key"
            RESULTS_PATHS[${#RESULTS_PATHS[@]}]="$results_dir"
            
            # Verify what was actually added
            local new_length=${#RESULTS_KEYS[@]}
            local last_index=$((new_length - 1))
            echo "    -> After adding:"
            echo "       new array length: $new_length"
            echo "       last key: '${RESULTS_KEYS[$last_index]}'"
            echo "       last path: '${RESULTS_PATHS[$last_index]}'"
            
            # Verify the key was stored correctly
            if [[ "${RESULTS_KEYS[$last_index]}" != "$key" ]]; then
                echo "    ERROR: Key mismatch! Expected '$key', got '${RESULTS_KEYS[$last_index]}'"
            fi
            
            if [[ "${RESULTS_PATHS[$last_index]}" != "$results_dir" ]]; then
                echo "    ERROR: Path mismatch! Expected '$results_dir', got '${RESULTS_PATHS[$last_index]}'"
            fi
            
            # Show all current array contents
            echo "    -> Current complete arrays:"
            for (( j=0; j<${#RESULTS_KEYS[@]}; j++ )); do
                echo "       [$j] '${RESULTS_KEYS[j]}' -> '${RESULTS_PATHS[j]}'"
            done
            
            # Verify the arrays are synchronized
            if [ ${#RESULTS_KEYS[@]} -ne ${#RESULTS_PATHS[@]} ]; then
                echo "    ERROR: Array mismatch! Keys: ${#RESULTS_KEYS[@]}, Paths: ${#RESULTS_PATHS[@]}"
                exit 1
            fi
            
            # Run individual visualization
            if ! run_individual_visualization "$results_dir" "$sim_type" "$event_type"; then
                echo "    WARNING: Individual visualization failed for $sim_type $event_type"
            fi
        done
    done
    
    echo ""
    echo "=== SIMULATION RESULTS SUMMARY ==="
    echo "Raw array lengths: keys=${#RESULTS_KEYS[@]}, paths=${#RESULTS_PATHS[@]}"
    
    # Clean up any empty elements before proceeding
    echo "Cleaning arrays..."
    local clean_keys=()
    local clean_paths=()
    
    for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
        local key="${RESULTS_KEYS[i]}"
        local path="${RESULTS_PATHS[i]}"
        
        echo "  Raw [$i]: key='$key', path='$path'"
        
        # Only keep non-empty entries
        if [[ -n "$key" && -n "$path" ]]; then
            clean_keys+=("$key")
            clean_paths+=("$path")
            echo "    ✓ Kept"
        else
            echo "    ✗ Removed (empty)"
        fi
    done
    
    # Replace with cleaned arrays
    RESULTS_KEYS=("${clean_keys[@]}")
    RESULTS_PATHS=("${clean_paths[@]}")
    
    echo ""
    echo "After cleaning - Successfully completed simulations:"
    for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
        echo "  [$i] ${RESULTS_KEYS[i]}: ${RESULTS_PATHS[i]}"
    done
    echo ""
    
    # Create pairwise comparisons for each event type
    echo "=== STEP 2: Creating Pairwise Comparisons ==="
    
    if [ "$ENABLE_PAIRWISE_VIZ" = false ]; then
        echo "Skipping pairwise comparisons (disabled)"
    else
        for event_type in "${EVENT_TYPES[@]}"; do
            echo ""
            echo "Creating pairwise comparison for: $event_type"
            
            local memilio_key="memilio_${event_type}"
            local panvadere_key="panvadere_${event_type}"
            
            local memilio_results=$(find_result_path "$memilio_key")
            local panvadere_results=$(find_result_path "$panvadere_key")
            
            echo "  -> Looking for: $memilio_key and $panvadere_key"
            echo "  -> Found memilio: '$memilio_results'"
            echo "  -> Found panvadere: '$panvadere_results'"
            
            if [[ -n "$memilio_results" && -n "$panvadere_results" ]]; then
                echo "  -> SUCCESS: Found both results, creating comparison..."
                if ! run_pairwise_comparison "$memilio_results" "$panvadere_results" "$event_type" "$viz_output_dir"; then
                    echo "  WARNING: Some pairwise visualizations may have failed for $event_type"
                fi
            else
                echo "  WARNING: Missing results for $event_type"
                echo "    memilio: '$memilio_results'"
                echo "    panvadere: '$panvadere_results'"
            fi
        done
    fi
    
    # Create comprehensive comparison
    echo ""
    echo "=== STEP 3: Creating Comprehensive Comparison ==="
    
    if [ "$ENABLE_COMPREHENSIVE_VIZ" = false ]; then
        echo "Skipping comprehensive comparison (disabled)"
    else
        # Prepare arrays for the Python script
        local all_paths=()
        local all_labels=()
        
        for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
            local key="${RESULTS_KEYS[i]}"
            local path="${RESULTS_PATHS[i]}"
            if [[ -n "$path" && -d "$path" ]]; then
                all_paths+=("$path")
                all_labels+=("$key")
            fi
        done
        
        if [ ${#all_paths[@]} -gt 0 ]; then
            echo "Running comprehensive comparison with ${#all_paths[@]} scenarios..."
            echo "Labels: ${all_labels[*]}"
            
            $PYTHON3_DIR "$COMPARISON_PYTHON_SCRIPT" \
                --results-paths "${all_paths[@]}" \
                --labels "${all_labels[@]}" \
                --output-dir "$viz_output_dir" \
                --timestamp "$timestamp" \
                $VIZ_OPTIONS
        else
            echo "WARNING: No valid results found for comprehensive comparison"
        fi
    fi
    
    # Copy results to last_result folders
    echo ""
    echo "=== STEP 4: Copying Results to Last Result Folders ==="
    
    for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
        local key="${RESULTS_KEYS[i]}"
        local results_path="${RESULTS_PATHS[i]}"
        
        # Extract sim_type and event_type from key
        local sim_type="${key%%_*}"  # Everything before first underscore
        local event_type="${key#*_}" # Everything after first underscore
        
        if [[ -n "$results_path" && -d "$results_path/infection_state_per_age_group" ]]; then
            copy_results_to_last_result "$results_path" "$sim_type" "$event_type"
        else
            echo "  -> No infection state results found for $key"
        fi
    done
    
    echo ""
    echo "=== COMPREHENSIVE SIMULATION ANALYSIS COMPLETED SUCCESSFULLY ==="
    echo "Results summary:"
    for (( i=0; i<${#RESULTS_KEYS[@]}; i++ )); do
        echo "  ${RESULTS_KEYS[i]}: ${RESULTS_PATHS[i]}"
    done
    echo ""
    echo "Comprehensive visualizations: $viz_output_dir"
    echo "Individual pairwise comparisons available in subfolders of: $viz_output_dir"
}


# Run main function
main "$@"