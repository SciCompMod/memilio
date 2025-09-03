#!/bin/zsh

# =============================================================================
# COMPREHENSIVE EPIDEMIC VISUALIZATION SCRIPT
# =============================================================================
# This script uses the output from run_epidemic_curves.sh and creates:
# 1. Multi-panel epidemic curve comparisons (2x2 grid)
# 2. Temporal heatmap visualizations for each scenario
# 3. Comparative temporal heatmaps showing spatial-temporal patterns
# 4. Organized output folder structure
# =============================================================================

# Configuration
MAIN_PATH="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp"
PYTHON3_DIR="$MAIN_PATH/v_m/bin/python3"
RESULTS_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/results"
VIZ_BASE_DIR="$MAIN_PATH/examples/panvXabmSim/viz"

# Visualization Control Flags (set to true/false to enable/disable)
ENABLE_EPIDEMIC_CURVES=true
ENABLE_COMPARATIVE_HEATMAPS=true   
ENABLE_INFECTION_TREES=true
ENABLE_LOCATION_PIE_CHARTS=true

# Visualization scripts
SIMPLE_VIZ_SCRIPT="$VIZ_BASE_DIR/simple_multi_panel.py"
COMPARATIVE_HEATMAP_SCRIPT="$VIZ_BASE_DIR/comparative_temporal_heatmap.py"
TEMPORAL_HEATMAP_SCRIPT="$VIZ_BASE_DIR/temporal_infection_heatmap.py"
COMPARATIVE_TREES_SCRIPT="$VIZ_BASE_DIR/comparative_infection_trees.py"
LOCATION_PIE_SCRIPT="$VIZ_BASE_DIR/location_infection_pie_chart.py"

# Create comprehensive output directory
create_output_dir() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local output_dir="$RESULTS_BASE_DIR/comprehensive_visualization_$timestamp"
    mkdir -p "$output_dir"
    echo "$output_dir"
}

# Find scenario results based on pattern matching
find_scenario_results() {
    local scenario_pattern=$1
    local init_type=$2
    
    # Look for directories matching the pattern and sort by date (newest first)
    local results_dir=$(find "$RESULTS_BASE_DIR" -type d -name "*${scenario_pattern}*${init_type}*" | sort -r | head -1)
    
    if [ -n "$results_dir" ] && [ -d "$results_dir" ]; then
        echo "$results_dir"
        return 0
    else
        echo ""
        return 1
    fi
}

# Create epidemic curves visualization
create_epidemic_curves() {
    local output_dir=$1
    
    if [ "$ENABLE_EPIDEMIC_CURVES" != "true" ]; then
        echo "⏭️  Skipping epidemic curves (disabled)"
        return 0
    fi
    
    echo "=== Creating Multi-Panel Epidemic Curves ==="
    
    # Create subdirectory for epidemic curves
    local curves_dir="$output_dir/epidemic_curves"
    mkdir -p "$curves_dir"
    
    # Run the simple multi-panel visualization
    echo "Generating 2x2 epidemic curve comparison..."
    cd "$curves_dir"
    $PYTHON3_DIR "$SIMPLE_VIZ_SCRIPT" --results-dir "$RESULTS_BASE_DIR" --s90percentile
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "✅ Epidemic curves visualization completed"
        echo "   Output: $curves_dir"
    else
        echo "❌ Failed to create epidemic curves visualization"
        return $exit_code
    fi
    
    return 0
}

# Create heatmap visualizations for a scenario pair
create_scenario_heatmaps() {
    local scenario_label=$1
    local scenario_pattern=$2
    local output_dir=$3
    
    echo ""
    echo "=== Creating Heatmaps for $scenario_label ==="
    
    # Find transmission-informed and uniform results
    local transmission_dir=$(find_scenario_results "$scenario_pattern" "transmission_informed")
    local uniform_dir=$(find_scenario_results "$scenario_pattern" "uniform_initialized")
    
    if [ -z "$transmission_dir" ] || [ -z "$uniform_dir" ]; then
        echo "❌ Could not find complete results for $scenario_label"
        echo "   Transmission-informed: $transmission_dir"
        echo "   Uniform initialized: $uniform_dir"
        return 1
    fi
    
    echo "Found results:"
    echo "  Transmission-informed: $(basename "$transmission_dir")"
    echo "  Uniform initialized: $(basename "$uniform_dir")"
    
    # Create subdirectory for this scenario
    local scenario_dir="$output_dir/heatmaps/$scenario_label"
    mkdir -p "$scenario_dir"
    
    # Extract scenario code (R1, R2, W1, W2) from scenario label
    local scenario_code="${scenario_label:0:2}"
    
    # Determine if this is a workplace scenario (W1, W2) or restaurant scenario (R1, R2)
    local use_workplaces=""
    local view_type="households"
    local viz_style="rectangles"
    if [[ "$scenario_code" == "W1" ]] || [[ "$scenario_code" == "W2" ]]; then
        use_workplaces="--use-workplaces"
        view_type="workplaces"
        viz_style="hexagons"
        echo "Using workplace-based visualization for $scenario_code with hexagon style"
    else
        echo "Using household-based visualization for $scenario_code with rectangle style"
    fi
    
    # Create comparative temporal heatmap (average across runs only)
    if [ "$ENABLE_COMPARATIVE_HEATMAPS" = "true" ]; then
        echo "Creating comparative temporal heatmap (averaged across runs) - $view_type view..."
        cd "$scenario_dir"
        $PYTHON3_DIR "$COMPARATIVE_HEATMAP_SCRIPT" \
            --data-dir-method1 "$transmission_dir" \
            --data-dir-method2 "$uniform_dir" \
            --method1-name "Transmission-Informed" \
            --method2-name "Uniform" \
            --time-points 0 24 72 240 \
            --output-path "$scenario_dir/comparative_temporal_heatmap_averaged_${view_type}_${scenario_code}.png" \
            --scenario-name "$scenario_code" \
            --viz-style "$viz_style" \
            --use-median \
            --use-size-based-colors \
            $use_workplaces
        
        if [ $? -eq 0 ]; then
            echo "✅ Averaged heatmap ($view_type view) created"
        else
            echo "❌ Failed to create averaged heatmap ($view_type view)"
        fi
    else
        echo "⏭️  Skipping comparative heatmaps (disabled)"
    fi
    
    # Create comparative infection trees if script exists
    if [ "$ENABLE_INFECTION_TREES" = "true" ] && [ -f "$COMPARATIVE_TREES_SCRIPT" ]; then
        echo "Creating comparative infection trees..."
        cd "$scenario_dir"
        $PYTHON3_DIR "$COMPARATIVE_TREES_SCRIPT" \
            --data-dir-method1 "$transmission_dir" \
            --data-dir-method2 "$uniform_dir" \
            --method1-name "Transmission-Informed" \
            --method2-name "Uniform" \
            --output-path "$scenario_dir/comparative_infection_trees_${scenario_code}.png" \
            --scenario-name "$scenario_code" \
            --connection-style "current" \
            --layout-method "depth-first"

        if [ $? -eq 0 ]; then
            echo "✅ Infection trees created for $scenario_code"
        else
            echo "❌ Failed to create infection trees for $scenario_code"
        fi
    elif [ "$ENABLE_INFECTION_TREES" != "true" ]; then
        echo "⏭️  Skipping infection trees (disabled)"
    fi
    
    # Create location infection pie chart if script exists
    if [ "$ENABLE_LOCATION_PIE_CHARTS" = "true" ] && [ -f "$LOCATION_PIE_SCRIPT" ]; then
        echo "Creating location infection pie chart..."
        cd "$scenario_dir"
        
        # Create comparative pie chart
        $PYTHON3_DIR "$LOCATION_PIE_SCRIPT" \
            --comparative \
            --data-dir "$transmission_dir" \
            --data-dir-2 "$uniform_dir" \
            --scenario-name "Transmission-Informed" \
            --scenario-name-2 "Uniform" \
            --output "$scenario_dir/comparative_location_pie_${scenario_code}.png" \
            --title-suffix "$scenario_code"

        if [ $? -eq 0 ]; then
            echo "✅ Location infection pie charts created for $scenario_code"
        else
            echo "❌ Failed to create location pie charts for $scenario_code"
        fi
    elif [ "$ENABLE_LOCATION_PIE_CHARTS" != "true" ]; then
        echo "⏭️  Skipping location pie charts (disabled)"
    fi
    
    echo "✅ All visualizations completed for $scenario_label"
    echo "   Output: $scenario_dir"
    
    return 0
}

# Validate prerequisites
validate_prerequisites() {
    echo "=== Validating Prerequisites ==="
    
    # Check if epidemic curve results exist
    if [ ! -d "$RESULTS_BASE_DIR" ]; then
        echo "❌ Results directory not found: $RESULTS_BASE_DIR"
        echo "Please run ./run_epidemic_curves.sh first to generate the data."
        exit 1
    fi
    
    # Count available result directories
    local epidemic_dirs=$(ls "$RESULTS_BASE_DIR" | grep "epidemic_curves_.*_transmission_informed" | wc -l)
    if [ "$epidemic_dirs" -eq 0 ]; then
        echo "❌ No epidemic curve results found in $RESULTS_BASE_DIR"
        echo "Please run ./run_epidemic_curves.sh first to generate the data."
        exit 1
    fi
    
    echo "✅ Found $epidemic_dirs scenario results"
    
    # Check visualization scripts
    if [ ! -f "$SIMPLE_VIZ_SCRIPT" ]; then
        echo "❌ Simple visualization script not found: $SIMPLE_VIZ_SCRIPT"
        exit 1
    fi
    
    if [ ! -f "$COMPARATIVE_HEATMAP_SCRIPT" ]; then
        echo "❌ Comparative heatmap script not found: $COMPARATIVE_HEATMAP_SCRIPT"
        exit 1
    fi
    
    # Check Python environment
    if [ ! -f "$PYTHON3_DIR" ]; then
        echo "❌ Python3 not found: $PYTHON3_DIR"
        echo "Please check your Python virtual environment setup."
        exit 1
    fi
    
    echo "✅ All prerequisites validated"
    echo ""
}

# Main execution function
main() {
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    
    echo "=== Comprehensive Epidemic Visualization Generator ==="
    echo "Timestamp: $timestamp"
    echo "This script generates epidemic curves and heatmap visualizations"
    echo "from the results of run_epidemic_curves.sh"
    echo ""
    
    # Show visualization settings
    echo "=== Visualization Settings ==="
    echo "📊 Epidemic Curves:       $([ "$ENABLE_EPIDEMIC_CURVES" = "true" ] && echo "✅ ENABLED" || echo "❌ DISABLED")"
    echo "🗺️  Comparative Heatmaps:  $([ "$ENABLE_COMPARATIVE_HEATMAPS" = "true" ] && echo "✅ ENABLED" || echo "❌ DISABLED")"
    echo "🌳 Infection Trees:       $([ "$ENABLE_INFECTION_TREES" = "true" ] && echo "✅ ENABLED" || echo "❌ DISABLED")"
    echo "🥧 Location Pie Charts:   $([ "$ENABLE_LOCATION_PIE_CHARTS" = "true" ] && echo "✅ ENABLED" || echo "❌ DISABLED")"
    echo ""
    
    # Validate prerequisites
    validate_prerequisites
    
    # Create main output directory
    local output_dir=$(create_output_dir)
    echo "=== Creating Comprehensive Visualizations ==="
    echo "Output directory: $output_dir"
    echo ""
    
    # Create epidemic curves
    if ! create_epidemic_curves "$output_dir"; then
        echo "❌ Failed to create epidemic curves"
        exit 1
    fi
    
    # Define scenarios to process (using explicit declaration to avoid zsh array indexing issues)
    declare -a scenario_labels
    declare -a scenario_patterns
    scenario_labels=("R1_restaurant_strong" "R2_restaurant_weak" "W1_workplace_few" "W2_workplace_many")
    scenario_patterns=("R1_restaurant_strong_clustering" "R2_restaurant_weaker_clustering" "W1_workplace_few_meetings" "W2_workplace_many_meetings")
    
    # Create heatmaps directory
    mkdir -p "$output_dir/heatmaps"
    
    # Process each scenario
    local successful_scenarios=0
    local total_scenarios=${#scenario_labels[@]}
    
    for ((i=1; i<=${#scenario_labels[@]}; i++)); do
        local scenario_label="${scenario_labels[$i]}"
        local scenario_pattern="${scenario_patterns[$i]}"
        
        if create_scenario_heatmaps "$scenario_label" "$scenario_pattern" "$output_dir"; then
            successful_scenarios=$((successful_scenarios + 1))
        fi
    done
}

# Run main function
main "$@"
