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
VIZ_OUTPUT_DIR="$MAIN_PATH/examples/panvXabmSim/results/multi_seed_viz"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"
MULTI_SEED_VIZ_SCRIPT="$MAIN_PATH/examples/panvXabmSim/viz/multi_seed_comparison.py"

# Simulation parameters
EVENT_TYPE="work_meeting_baseline"  # Can be changed as needed
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

# Results storage arrays
declare -a ALL_RESULTS_PATHS=()
declare -a ALL_RESULTS_LABELS=()

# Functions
create_results_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local sim_type=$1
    local seed=$2
    local results_dir="$RESULTS_BASE_DIR/multiseed_${timestamp}_${sim_type}_seed${seed}"
    echo "$results_dir"
}

create_viz_output_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local viz_dir="$VIZ_OUTPUT_DIR/multi_seed_analysis_$timestamp"
    
    # Try to create the directory, with fallback options
    if mkdir -p "$viz_dir" 2>/dev/null; then
        echo "$viz_dir"
    else
        echo "Warning: Cannot create $viz_dir (permission denied)" >&2
        # Fallback to current user's home directory
        local fallback_dir="$HOME/multi_seed_analysis_$timestamp"
        mkdir -p "$fallback_dir"
        echo "Using fallback directory: $fallback_dir" >&2
        echo "$fallback_dir"
    fi
}

run_simulation_with_seed() {
    local results_dir=$1
    local sim_type=$2
    local seed=$3
    
    echo "  -> Running $sim_type simulation with seed $seed..."
    $MAIN_EXECUTABLE --event "$EVENT_TYPE" --sim "$sim_type" --output_dir "$results_dir" \
        --days "$NUM_DAYS" --n_persons "$NUM_PERSONS" --runs "$RUNS" \
        --infection_k "$INFECTION_K" --seed "$seed"
    return $?
}

# Main execution
main() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    echo "=== Multi-Seed Simulation Analysis ==="
    echo "Event Type: $EVENT_TYPE"
    echo "Number of seeds: ${#SEEDS[@]}"
    echo "Simulation types: ${SIM_TYPES[*]}"
    echo "Days: $NUM_DAYS, Persons: $NUM_PERSONS, Runs per seed: $RUNS"
    echo ""
    
    # Create main visualization output directory
    local viz_output_dir=$(create_viz_output_dir)
    echo "Visualization output directory: $viz_output_dir"
    echo ""
    
    # Run simulations for all seeds and simulation types
    echo "=== STEP 1: Running All Seed Simulations ==="
    
    local total_simulations=$((${#SEEDS[@]} * ${#SIM_TYPES[@]}))
    local current_sim=0
    
    for seed in "${SEEDS[@]}"; do
        echo ""
        echo "Processing seed: $seed"
        
        for sim_type in "${SIM_TYPES[@]}"; do
            current_sim=$((current_sim + 1))
            echo "  [$current_sim/$total_simulations] Running $sim_type with seed $seed"
            
            # Create results directory for this combination
            local results_dir=$(create_results_dir "$sim_type" "$seed")
            
            # Run simulation
            if ! run_simulation_with_seed "$results_dir" "$sim_type" "$seed"; then
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
            ALL_RESULTS_PATHS+=("$results_dir")
            ALL_RESULTS_LABELS+=("${sim_type}_seed${seed}")
            
            echo "    -> Added to visualization queue: ${sim_type}_seed${seed}"
        done
    done
    
    echo ""
    echo "=== STEP 2: Creating Multi-Seed Visualization ==="
    echo "Total successful simulations: ${#ALL_RESULTS_PATHS[@]}"
    
    if [ ${#ALL_RESULTS_PATHS[@]} -eq 0 ]; then
        echo "ERROR: No successful simulations found. Cannot create visualization."
        exit 1
    fi
    
    echo "Creating comprehensive multi-seed comparison..."
    echo "Results paths: ${#ALL_RESULTS_PATHS[@]} entries"
    echo "Labels: ${#ALL_RESULTS_LABELS[@]} entries"
    
    # Call Python script for multi-seed visualization
    echo "Running visualization script..."
    $PYTHON3_DIR "$MULTI_SEED_VIZ_SCRIPT" \
        --results-paths "${ALL_RESULTS_PATHS[@]}" \
        --labels "${ALL_RESULTS_LABELS[@]}" \
        --output-dir "$viz_output_dir" \
        --event-type "$EVENT_TYPE" \
        --num-seeds "${#SEEDS[@]}"
    
    if [ $? -eq 0 ]; then
        echo "✓ Multi-seed visualization completed successfully!"
    else
        echo "WARNING: Multi-seed visualization may have encountered issues."
    fi
    
    echo ""
    echo "=== MULTI-SEED ANALYSIS COMPLETED ==="
    echo "Total simulations run: ${#ALL_RESULTS_PATHS[@]}"
    echo "Visualization output: $viz_output_dir"
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
