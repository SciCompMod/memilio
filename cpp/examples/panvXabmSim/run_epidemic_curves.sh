#!/bin/zsh

# =============================================================================
# EPIDEMIC CURVE DYNAMICS ANALYSIS SCRIPT
# =============================================================================
# This script runs all 4 scenarios with both initialization types 
# (transmission-informed vs uniform) and generates epidemic curve comparisons.
# 
# Based on the paper's methodology:
# - "panvadere" = transmission-informed initialization (detailed Vadere data)
# - "memilio" = uniform initialization (case number informed)
# =============================================================================

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
COMPARISON_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/infection_cumulated_comp.py"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
VIZ_OUTPUT_DIR="$MAIN_PATH/examples/panvXabmSim/results/epidemic_curves_analysis"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"

# Event types mapping to paper scenarios
# R1 = restaurant_table_equals_household (strong household clustering)
# R2 = restaurant_table_equals_half_household (weaker household clustering) 
# W1 = work_meeting_baseline (few meetings - sequential clustering)
# W2 = work_meeting_many (many meetings - high mixing)
EVENT_TYPES=(
    "restaurant_table_equals_household"
    "restaurant_table_equals_half_household" 
    "work_meeting_baseline"
    "work_meeting_many"
)

# Paper scenario labels for output
SCENARIO_LABELS=(
    "R1_restaurant_strong_clustering"
    "R2_restaurant_weaker_clustering"
    "W1_workplace_few_meetings"
    "W2_workplace_many_meetings"
)

# Initialization types
# panvadere = transmission-informed (preserves social clustering from Vadere)
# memilio = uniform initialization (distributes infections uniformly)
INIT_TYPES=("panvadere" "memilio")
INIT_LABELS=("transmission_informed" "uniform_initialized")

# Simulation parameters
NUM_DAYS=10
NUM_PERSONS=1000
RUNS=100  # Increased for statistical significance
INFECTION_K=22.6
VIZ_OPTIONS="--s90percentile"

# Results storage
RESULTS_STORAGE_FILE="/tmp/epidemic_curves_results_$$"
rm -f "$RESULTS_STORAGE_FILE"
touch "$RESULTS_STORAGE_FILE"

# Helper functions
store_result() {
    local key=$1
    local path=$2
    echo "$key|$path" >> "$RESULTS_STORAGE_FILE"
}

find_result_path() {
    local search_key=$1
    while IFS='|' read -r key path; do
        if [[ "$key" == "$search_key" ]]; then
            echo "$path"
            return 0
        fi
    done < "$RESULTS_STORAGE_FILE"
    echo ""
    return 1
}

create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local scenario_label=$1
    local init_type=$2
    local results_dir="$RESULTS_BASE_DIR/epidemic_curves_${timestamp}_${scenario_label}_${init_type}"
    echo "$results_dir"
}

create_viz_output_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local viz_dir="$VIZ_OUTPUT_DIR/epidemic_analysis_$timestamp"
    mkdir -p "$viz_dir"
    echo "$viz_dir"
}

run_simulation() {
    local results_dir=$1
    local event_type=$2
    local init_type=$3
    local scenario_label=$4
    
    echo "    -> Running $init_type initialization for $scenario_label..."
    echo "       Event: $event_type, Output: $results_dir"
    
    $MAIN_EXECUTABLE \
        --event "$event_type" \
        --sim "$init_type" \
        --output_dir "$results_dir" \
        --days "$NUM_DAYS" \
        --n_persons "$NUM_PERSONS" \
        --runs "$RUNS" \
        --infection_k "$INFECTION_K"
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "       ✓ Completed successfully"
        return 0
    else
        echo "       ✗ Failed with exit code $exit_code"
        return $exit_code
    fi
}

create_epidemic_curve_comparison() {
    local transmission_informed_dir=$1
    local uniform_dir=$2
    local scenario_label=$3
    local viz_output_dir=$4
    
    echo "    -> Creating epidemic curve comparison for $scenario_label..."
    
    # Create scenario-specific output directory
    local scenario_output_dir="$viz_output_dir/$scenario_label"
    mkdir -p "$scenario_output_dir"
    
    # Generate comparison plot
    $PYTHON3_DIR "$COMPARISON_VIZ_SCRIPT" \
        --path-to-memilio-sim "$uniform_dir/amount_of_infections" \
        --path-to-panvXabmSim "$transmission_informed_dir/amount_of_infections" \
        --output-path "$scenario_output_dir" \
        $VIZ_OPTIONS
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "       ✓ Comparison plot created in: $scenario_output_dir"
        return 0
    else
        echo "       ✗ Failed to create comparison plot"
        return $exit_code
    fi
}

validate_prerequisites() {
    if [ ! -f "$MAIN_EXECUTABLE" ]; then
        echo "ERROR: Main executable not found: $MAIN_EXECUTABLE"
        echo "Please build the project first."
        exit 1
    fi

    if [ ! -f "$PYTHON3_DIR" ]; then
        echo "ERROR: Python3 not found: $PYTHON3_DIR"
        echo "Please check your Python virtual environment setup."
        exit 1
    fi

    if [ ! -f "$COMPARISON_VIZ_SCRIPT" ]; then
        echo "ERROR: Visualization script not found: $COMPARISON_VIZ_SCRIPT"
        exit 1
    fi
}

# Main execution
main() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    echo "=== Epidemic Curve Dynamics Analysis ==="
    echo "Analyzing transmission-informed vs. uniform initialization"
    echo "Timestamp: $timestamp"
    echo "Scenarios: ${SCENARIO_LABELS[*]}"
    echo "Runs per scenario: $RUNS"
    echo "Days: $NUM_DAYS, Persons: $NUM_PERSONS"
    echo ""

    # Validate prerequisites
    validate_prerequisites
    
    # Create main visualization output directory
    local viz_output_dir=$(create_viz_output_dir)
    echo "Visualization output directory: $viz_output_dir"
    echo ""
    
    # Track total simulations
    local total_scenarios=${#EVENT_TYPES[@]}
    local total_simulations=$((total_scenarios * 2))  # 2 init types per scenario
    local current_sim=0
    
    # Run simulations for all scenarios and initialization types
    echo "=== STEP 1: Running All Simulations ==="
    
    for ((i=1; i<=${#EVENT_TYPES[@]}; i++)); do
        local event_type="${EVENT_TYPES[$i]}"
        local scenario_label="${SCENARIO_LABELS[$i]}"
        
        echo ""
        echo "[$i/$total_scenarios] Processing Scenario: $scenario_label"
        echo "Event Type: $event_type"
        
        # Storage for this scenario's results
        local transmission_informed_dir=""
        local uniform_dir=""
        
        # Run both initialization types for this scenario
        for ((j=1; j<=${#INIT_TYPES[@]}; j++)); do
            local init_type="${INIT_TYPES[$j]}"
            local init_label="${INIT_LABELS[$j]}"
            
            current_sim=$((current_sim + 1))
            echo "  [$current_sim/$total_simulations] Running $init_label initialization"
            
            # Create results directory
            local results_dir=$(create_results_dir "$scenario_label" "$init_label")
            
            # Run simulation
            if ! run_simulation "$results_dir" "$event_type" "$init_type" "$scenario_label"; then
                echo "    ERROR: Simulation failed for $scenario_label with $init_label initialization"
                continue
            fi
            
            # Verify results exist
            if [ ! -d "$results_dir/amount_of_infections" ]; then
                echo "    ERROR: Results directory incomplete: $results_dir"
                continue
            fi
            
            # Store results path
            store_result "${scenario_label}_${init_label}" "$results_dir"
            
            # Keep track for comparison
            if [ "$init_type" = "panvadere" ]; then
                transmission_informed_dir="$results_dir"
            else
                uniform_dir="$results_dir"
            fi
            
            echo "    ✓ Results stored: $results_dir"
        done
        
        # Create comparison plot for this scenario
        if [ -n "$transmission_informed_dir" ] && [ -n "$uniform_dir" ]; then
            echo ""
            echo "  Creating epidemic curve comparison for $scenario_label..."
            if create_epidemic_curve_comparison "$transmission_informed_dir" "$uniform_dir" "$scenario_label" "$viz_output_dir"; then
                echo "  ✓ Comparison completed for $scenario_label"
            else
                echo "  ✗ Comparison failed for $scenario_label"
            fi
        else
            echo "  ✗ Missing results for $scenario_label - skipping comparison"
        fi
    done
    
    echo ""
    echo "=== STEP 2: Creating Summary Overview ==="
    
    # Create a summary report
    local summary_file="$viz_output_dir/analysis_summary.txt"
    cat > "$summary_file" << EOF
Epidemic Curve Dynamics Analysis Summary
========================================
Timestamp: $timestamp
Analysis Type: Transmission-Informed vs. Uniform Initialization
Total Scenarios: $total_scenarios
Runs per Scenario: $RUNS
Simulation Days: $NUM_DAYS
Population: $NUM_PERSONS

Scenarios Analyzed:
EOF
    
    for ((i=1; i<=${#SCENARIO_LABELS[@]}; i++)); do
        echo "  ${i}. ${SCENARIO_LABELS[$i]} (${EVENT_TYPES[$i]})" >> "$summary_file"
    done
    
    cat >> "$summary_file" << EOF

Results Location: $viz_output_dir

Individual scenario comparisons are available in:
EOF
    
    for scenario_label in "${SCENARIO_LABELS[@]}"; do
        echo "  - $viz_output_dir/$scenario_label/" >> "$summary_file"
    done
    
    echo ""
    echo "=== EPIDEMIC CURVE ANALYSIS COMPLETED ==="
    echo "Total simulations run: $total_simulations"
    echo "Results directory: $viz_output_dir"
    echo "Summary report: $summary_file"
    echo ""
    echo "Next steps:"
    echo "1. Review individual scenario comparisons in subdirectories"
    echo "2. Analyze differences in epidemic curves between initialization types"
    echo "3. Quantify impact on final epidemic size and peak timing"
    echo "=================================="
}

# Run main function
main "$@"
