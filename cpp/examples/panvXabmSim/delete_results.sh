#!/bin/zsh

# Configuration
RESULTS_BASE_DIR="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results"

# Function to delete all run_* folders
delete_run_folders() {
    echo "Deleting all run_* and results_* folders in: $RESULTS_BASE_DIR"

    if [ ! -d "$RESULTS_BASE_DIR" ]; then
        echo "Results directory does not exist: $RESULTS_BASE_DIR"
        return 1
    fi

    # Find and delete all run_* directories or results_* directories
    local run_folders=("$RESULTS_BASE_DIR"/run_*)
    run_folders+=("$RESULTS_BASE_DIR"/results_*)
    run_folders+=("$RESULTS_BASE_DIR"/results_viz/comparison_*)
    run_folders+=("$RESULTS_BASE_DIR"/debug_*)
    run_folders+=("$RESULTS_BASE_DIR"/multiseed_*)

    if [ ${#run_folders[@]} -eq 1 ] && [ ! -d "${run_folders[0]}" ]; then
        echo "No run_* or results_* folders found to delete."
        return 0
    fi
    
    for folder in "${run_folders[@]}"; do
        if [ -d "$folder" ]; then
            echo "Deleting: $folder"
            rm -rf "$folder"
        fi
    done
    
    echo "Cleanup completed."
}

# Main execution
main() {
    echo "Starting cleanup of simulation results..."
    delete_run_folders
}

# Run main function
main "$@"