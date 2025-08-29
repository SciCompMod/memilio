#!/bin/zsh

# =============================================================================
# MULTI-SEED SIMULATION SCRIPT
# =============================================================================
# This script runs simulations with 20 different seeds for both Memilio and 
# Panvadere simulations, then creates a comprehensive visualization showing
# all median runs for each seed.
# =============================================================================

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
VIZ_OUTPUT_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"
MULTI_SEED_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/multi_seed_comparison.py"

# Get event label function
get_event_label() {
    case "$1" in
        "restaurant_table_equals_household")
            echo "R1"
            ;;
        "restaurant_table_equals_half_household")
            echo "R2"
            ;;
        "work_meeting_many")
            echo "W1"
            ;;
        "work_meeting_baseline")
            echo "W2"
            ;;
        *)
            echo "UNKNOWN"
            ;;
    esac
}

# Event types to run
EVENT_TYPES=(
    "restaurant_table_equals_household"
    "restaurant_table_equals_half_household" 
    "work_meeting_many"
    "work_meeting_baseline"
)

# Simulation parameters
NUM_DAYS=10
NUM_PERSONS=1000
RUNS=20  # Number of runs per seed (to get median)
INFECTION_K=22.6



# Array of 20 different seeds
SEEDS=(
    1402121
    35897932
    27182818
    18284590
    45235360
    28747135
    99887766
    55443322
    11223344
    66778899
    44556677
    88990011
    22334455
    77889900
    33445566
    00112233
    55667788
    99001122
    33447799
    66880022
)

# Both simulation types
SIM_TYPES=("memilio" "panvadere")

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local sim_type=$1
    local seed=$2
    local event_type=$3
    local results_dir="$RESULTS_BASE_DIR/multiseed_${timestamp}_${sim_type}_${event_type}_seed${seed}"
    echo "$results_dir"
}

create_viz_output_dir() {
    local event_label=$1
    local viz_dir="$VIZ_OUTPUT_BASE_DIR/seeds_${event_label}"
    
    # Try to create the directory, with fallback options
    if mkdir -p "$viz_dir" 2>/dev/null; then
        echo "$viz_dir"
    else
        echo "Warning: Cannot create $viz_dir (permission denied)" >&2
        # Fallback to current user's home directory
        local fallback_dir="$HOME/seeds_${event_label}"
        mkdir -p "$fallback_dir"
        echo "Using fallback directory: $fallback_dir" >&2
        echo "$fallback_dir"
    fi
}

run_simulation_with_seed() {
    local results_dir=$1
    local sim_type=$2
    local seed=$3
    local event_type=$4
    
    echo "  -> Running $sim_type simulation with seed $seed for $event_type..."
    $MAIN_EXECUTABLE --event "$event_type" --sim "$sim_type" --output_dir "$results_dir" \
        --days "$NUM_DAYS" --n_persons "$NUM_PERSONS" --runs "$RUNS" \
        --infection_k "$INFECTION_K" --seed "$seed"
    return $?
}

# Main execution
main() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    echo "=== Multi-Seed Simulation Analysis for All Scenarios ==="
    echo "Event Types: ${EVENT_TYPES[*]}"
    echo "Number of seeds: ${#SEEDS[@]}"
    echo "Simulation types: ${SIM_TYPES[*]}"
    echo "Days: $NUM_DAYS, Persons: $NUM_PERSONS, Runs per seed: $RUNS"
    echo ""
    
    # Process each event type separately
    for event_type in "${EVENT_TYPES[@]}"; do
        local event_label=$(get_event_label "$event_type")
        echo ""
        echo "=== PROCESSING EVENT TYPE: $event_type (Label: $event_label) ==="
        
        # Create visualization output directory for this event type
        local viz_output_dir=$(create_viz_output_dir "$event_label")
        echo "Visualization output directory: $viz_output_dir"
        echo ""
        
        # Results storage arrays for this event type
        local -a results_paths=()
        local -a results_labels=()
        
        # Run simulations for all seeds and simulation types for this event type
        echo "=== STEP 1: Running All Seed Simulations for $event_type ==="
        
        local total_simulations=$((${#SEEDS[@]} * ${#SIM_TYPES[@]}))
        local current_sim=0
        
        for seed in "${SEEDS[@]}"; do
            echo ""
            echo "Processing seed: $seed for $event_type"
            
            for sim_type in "${SIM_TYPES[@]}"; do
                current_sim=$((current_sim + 1))
                echo "  [$current_sim/$total_simulations] Running $sim_type with seed $seed"
                
                # Create results directory for this combination
                local results_dir=$(create_results_dir "$sim_type" "$seed" "$event_type")
                
                # Run simulation
                if ! run_simulation_with_seed "$results_dir" "$sim_type" "$seed" "$event_type"; then
                    echo "    ERROR: $sim_type simulation failed for seed $seed"
                    echo "    -> Skipping this result"
                    continue
                fi
                
                # Verify results exist
                if [ ! -d "$results_dir" ] || [ ! -d "$results_dir/amount_of_infections" ]; then
                    echo "    ERROR: Results directory incomplete: $results_dir"
                    continue
                fi
                
                echo "    SUCCESS: Results saved in $results_dir"
                
                # Store results for visualization
                results_paths+=("$results_dir")
                results_labels+=("${sim_type}_seed${seed}")
                
                echo "    -> Added to visualization queue: ${sim_type}_seed${seed}"
            done
        done
        
        echo ""
        echo "=== STEP 2: Creating Multi-Seed Visualization for $event_type ==="
        echo "Total successful simulations: ${#results_paths[@]}"
        
        if [ ${#results_paths[@]} -eq 0 ]; then
            echo "ERROR: No successful simulations found for $event_type. Cannot create visualization."
            continue
        fi
        
        echo "Creating comprehensive multi-seed comparison for $event_label..."
        echo "Results paths: ${#results_paths[@]} entries"
        echo "Labels: ${#results_labels[@]} entries"
        
        # Call Python script for multi-seed visualization
        echo "Running visualization script..."
        $PYTHON3_DIR "$MULTI_SEED_VIZ_SCRIPT" \
            --results-paths "${results_paths[@]}" \
            --labels "${results_labels[@]}" \
            --output-dir "$viz_output_dir" \
            --event-type "$event_label" \
            --num-seeds "${#SEEDS[@]}"
        
        if [ $? -eq 0 ]; then
            echo "✓ Multi-seed visualization completed successfully for $event_label!"
        else
            echo "WARNING: Multi-seed visualization may have encountered issues for $event_label."
        fi
        
        echo ""
        echo "=== COMPLETED ANALYSIS FOR $event_type (Label: $event_label) ==="
        echo "Total simulations run: ${#results_paths[@]}"
        echo "Visualization output: $viz_output_dir"
        echo "=========================="
    done
    
    echo ""
    echo "=== ALL MULTI-SEED ANALYSES COMPLETED ==="
    echo "All event types processed: ${EVENT_TYPES[*]}"
    echo "=========================="
}

# Validate prerequisites
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

# Run main function
main "$@"
