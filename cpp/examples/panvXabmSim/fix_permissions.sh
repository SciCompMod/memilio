#!/bin/zsh

# =============================================================================
# PERMISSION FIX SCRIPT
# =============================================================================
# This script fixes permission issues that may occur when directories are
# created with sudo and then accessed by regular user scripts.
# =============================================================================

RESULTS_DIR="/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results"

echo "=== Fixing Permissions for Multi-Seed Analysis ==="

if [ -d "$RESULTS_DIR" ]; then
    echo "Checking permissions on results directory..."
    
    # Check if there are any root-owned directories
    ROOT_DIRS=$(find "$RESULTS_DIR" -user root 2>/dev/null)
    
    if [ -n "$ROOT_DIRS" ]; then
        echo "Found directories owned by root:"
        echo "$ROOT_DIRS"
        echo ""
        echo "Fixing permissions..."
        
        # Fix ownership
        sudo chown -R saschakorf:staff "$RESULTS_DIR"
        
        # Ensure proper permissions
        chmod -R u+rwx,g+rx,o+rx "$RESULTS_DIR"
        
        echo "✓ Permissions fixed successfully!"
    else
        echo "✓ No permission issues found."
    fi
else
    echo "Results directory does not exist: $RESULTS_DIR"
    echo "Creating it with correct permissions..."
    mkdir -p "$RESULTS_DIR"
    echo "✓ Results directory created."
fi

echo ""
echo "=== Permission Fix Complete ==="
echo "You can now run multi-seed scripts without sudo."
