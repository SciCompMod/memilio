#!/bin/zsh

# Simple Memilio vs Panvadere Contact Analysis Runner
# ===================================================

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
MAIN_EXECUTABLE="$MAIN_PATH/build/bin/panvXabm"
PYTHON3_PATH="$MAIN_PATH/v_m/bin/python3"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
ANALYSIS_SCRIPT="$MAIN_PATH/examples/panvXabmSim/analyze_contacts.py"

# Simulation parameters
EVENT_TYPE="work_meeting_baseline"
NUM_DAYS=10
NUM_PERSONS=1000
RUNS=20
INFECTION_K=22.6

# Generate timestamp for this run
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "=== Contact Analysis Runner ==="
echo "Event: $EVENT_TYPE | Days: $NUM_DAYS | Persons: $NUM_PERSONS | Runs: $RUNS"
echo "Timestamp: $TIMESTAMP"
echo ""

# Create result directories
MEMILIO_DIR="$RESULTS_BASE_DIR/contact_analysis_${TIMESTAMP}_memilio"
PANVADERE_DIR="$RESULTS_BASE_DIR/contact_analysis_${TIMESTAMP}_panvadere"

echo "Results will be saved to:"
echo "  Memilio:   $MEMILIO_DIR"
echo "  Panvadere: $PANVADERE_DIR"
echo ""

# Function to run simulation
run_simulation() {
    local sim_type=$1
    local results_dir=$2
    
    echo "=== Running $sim_type Simulation ==="
    echo "Output directory: $results_dir"
    
    $MAIN_EXECUTABLE \
        --event "$EVENT_TYPE" \
        --sim "$sim_type" \
        --output_dir "$results_dir" \
        --days "$NUM_DAYS" \
        --n_persons "$NUM_PERSONS" \
        --runs "$RUNS" \
        --infection_k "$INFECTION_K"
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "✓ $sim_type simulation completed successfully"
    else
        echo "✗ $sim_type simulation failed with exit code $exit_code"
        return $exit_code
    fi
    echo ""
    return 0
}

# Run Memilio simulation
if ! run_simulation "memilio" "$MEMILIO_DIR"; then
    echo "Failed to run Memilio simulation. Exiting."
    exit 1
fi

# Run Panvadere simulation
if ! run_simulation "panvadere" "$PANVADERE_DIR"; then
    echo "Failed to run Panvadere simulation. Exiting."
    exit 1
fi

# Check if required CSV files exist
echo "=== Checking Results ==="
required_files=("best_run_detailed_infection.csv" "best_run_contact_data.csv")

for sim_dir in "$MEMILIO_DIR" "$PANVADERE_DIR"; do
    sim_name=$(basename "$sim_dir" | cut -d'_' -f4)
    echo "Checking $sim_name results..."
    
    for file in "${required_files[@]}"; do
        if [ -f "$sim_dir/$file" ]; then
            echo "  ✓ $file found"
        else
            echo "  ✗ $file missing"
        fi
    done
done
echo ""

# Run contact analysis
echo "=== Running Contact Analysis ==="
echo "Analyzing potential contacts for infectious persons..."

if [ -f "$ANALYSIS_SCRIPT" ]; then
    # Activate virtual environment and run analysis with rolling averages
    source venv/bin/activate
    python "$ANALYSIS_SCRIPT" \
        --memilio-dir "$MEMILIO_DIR" \
        --panvadere-dir "$PANVADERE_DIR" \
        --output-file "contact_analysis_${TIMESTAMP}.png" \
        --title "Contact Analysis with 24h Rolling Avg - $TIMESTAMP" \
        --rolling-window 24
    
    analysis_exit_code=$?
    if [ $analysis_exit_code -eq 0 ]; then
        echo "✓ Contact analysis completed successfully"
        echo "✓ Visualization saved as: contact_analysis_${TIMESTAMP}.png"
    else
        echo "✗ Contact analysis failed with exit code $analysis_exit_code"
        exit 1
    fi
else
    echo "✗ Analysis script not found: $ANALYSIS_SCRIPT"
    echo "Please run this script from the panvXabmSim directory"
    exit 1
fi

echo ""
echo "=== Analysis Complete ==="
echo "Memilio results:   $MEMILIO_DIR"
echo "Panvadere results: $PANVADERE_DIR"
echo "Visualization:     contact_analysis_${TIMESTAMP}.png"
echo ""
