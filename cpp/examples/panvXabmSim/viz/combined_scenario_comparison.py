import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
from datetime import datetime
from pathlib import Path


def load_h5_results(base_path, percentile="p50"):
    """Load HDF5 results for a given group and percentile."""
    file_path = os.path.join(
        base_path, "amount_of_infections", percentile, "Results.h5")
    try:
        with h5py.File(file_path, 'r') as f:
            data = {k: v[()] for k, v in f['0'].items()}
        return data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None


def calculate_cumulative_infections(data):
    """Calculate cumulative infections from infection state data."""
    if data is None:
        return None
    current_infections = data.sum(axis=1)
    return current_infections


def parse_label(label):
    """Parse label to extract simulation type and seed."""
    parts = label.split('_')
    if len(parts) >= 2:
        sim_type = parts[0]
        seed_part = parts[1]
        if seed_part.startswith('seed'):
            seed = seed_part[4:]
            return sim_type, seed
    return label, "unknown"


def group_results_by_simulation_type(results_paths, labels):
    """Group results by simulation type and extract seed information."""
    grouped_results = {'memilio': [], 'panvadere': []}

    for path, label in zip(results_paths, labels):
        sim_type, seed = parse_label(label)
        data = load_h5_results(path, "p50")
        if data is not None:
            cumulative = calculate_cumulative_infections(data['Total'])
            if cumulative is not None:
                grouped_results[sim_type].append({
                    'seed': seed,
                    'time': data['Time'],
                    'cumulative': cumulative,
                    'path': path,
                    'label': label
                })

    return grouped_results


def get_scenario_full_name(event_label):
    """Convert event label to full scenario name."""
    scenario_names = {
        'R1': 'R1: Limited Inter-Household Mixing',
        'R2': 'R2: More Inter-Household Mixing',
        'W1': 'W1: Few Meetings, Limited Mixing',
        'W2': 'W2: More Meetings, More Mixing'
    }
    return scenario_names.get(event_label, event_label)


def create_combined_median_comparison(scenario_data_list, output_dir):
    """
    Create a 2x2 grid showing median comparisons for all four scenarios.
    
    Args:
        scenario_data_list: List of tuples (event_label, grouped_results)
        output_dir: Directory to save the output
    """
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 12
    })

    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    fig.suptitle('Multi-Seed Median Comparison Across All Scenarios',
                 fontsize=28, fontweight='bold', y=0.98)

    # Flatten axes for easier iteration
    axes_flat = axes.flatten()

    # Find global y-axis max for consistent scaling
    global_y_max = 0
    for event_label, grouped_results in scenario_data_list:
        for sim_type in ['memilio', 'panvadere']:
            if grouped_results[sim_type]:
                for result in grouped_results[sim_type]:
                    local_max = np.max(result['cumulative'])
                    global_y_max = max(global_y_max, local_max)

    # Plot each scenario
    for idx, (event_label, grouped_results) in enumerate(scenario_data_list):
        ax = axes_flat[idx]
        
        # Set subplot title (scenario name)
        scenario_name = get_scenario_full_name(event_label)
        ax.set_title(scenario_name, fontsize=18, fontweight='bold', pad=10)

        # Calculate statistics for Memilio
        if grouped_results['memilio']:
            # Find minimum length to handle different array sizes
            min_len = min(len(result['cumulative']) for result in grouped_results['memilio'])
            memilio_curves = np.array([result['cumulative'][:min_len]
                                      for result in grouped_results['memilio']])
            memilio_median = np.median(memilio_curves, axis=0)
            memilio_q25 = np.percentile(memilio_curves, 25, axis=0)
            memilio_q75 = np.percentile(memilio_curves, 75, axis=0)
            time_axis = grouped_results['memilio'][0]['time'][:min_len]

            # Plot Memilio
            ax.plot(time_axis, memilio_median, color='blue',
                   linewidth=3, label='Uniform (median)')
            ax.fill_between(time_axis, memilio_q25, memilio_q75,
                           color='blue', alpha=0.3, label='Uniform (25-75%)')

        # Calculate statistics for Panvadere
        if grouped_results['panvadere']:
            # Find minimum length to handle different array sizes
            min_len = min(len(result['cumulative']) for result in grouped_results['panvadere'])
            panvadere_curves = np.array([result['cumulative'][:min_len]
                                        for result in grouped_results['panvadere']])
            panvadere_median = np.median(panvadere_curves, axis=0)
            panvadere_q25 = np.percentile(panvadere_curves, 25, axis=0)
            panvadere_q75 = np.percentile(panvadere_curves, 75, axis=0)
            time_axis = grouped_results['panvadere'][0]['time'][:min_len]

            # Plot Panvadere
            ax.plot(time_axis, panvadere_median, color='red',
                   linewidth=3, label='Transmission-Informed (median)')
            ax.fill_between(time_axis, panvadere_q25, panvadere_q75,
                           color='red', alpha=0.3, label='Transmission-Informed (25-75%)')

        # Set consistent y-axis limits
        ax.set_ylim(0, global_y_max * 1.05)
        
        # Labels and grid
        ax.set_xlabel('Days', fontsize=16)
        ax.set_ylabel('Cumulative Infections', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12, loc='upper left')
        ax.tick_params(labelsize=14)

    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(output_dir, "combined_median_comparison_all_scenarios.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Combined median comparison plot saved to {output_file}")
    plt.close()

    return output_file


def main():
    """Main function for CLI usage."""
    parser = argparse.ArgumentParser(
        description="Create combined median comparison across all scenarios")
    parser.add_argument("--scenario-dirs", nargs='+', required=True,
                        help="Directories containing results for each scenario")
    parser.add_argument("--scenario-labels", nargs='+', required=True,
                        help="Labels for each scenario (e.g., R1, R2, W1, W2)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for the combined visualization")

    args = parser.parse_args()

    # Validate inputs
    if len(args.scenario_dirs) != len(args.scenario_labels):
        print("Error: Number of scenario directories must match number of labels")
        return 1

    # Create output directory
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"✓ Output directory ready: {args.output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        return 1

    print(f"Processing {len(args.scenario_dirs)} scenarios...")

    # Load data for each scenario
    scenario_data_list = []
    
    for scenario_dir, event_label in zip(args.scenario_dirs, args.scenario_labels):
        print(f"\nLoading scenario {event_label} from {scenario_dir}")
        
        # The scenario_dir is the viz output dir (seeds_XX), but we need to look in the parent
        # directory for the actual multiseed results
        results_base = os.path.dirname(scenario_dir)
        
        if not os.path.exists(results_base):
            print(f"Warning: Results directory not found: {results_base}")
            continue
            
        # Look for subdirectories with results matching this event label
        result_paths = []
        result_labels = []
        
        # Map event labels to their event type strings
        event_type_map = {
            'R1': 'restaurant_table_equals_household',
            'R2': 'restaurant_table_equals_half_household',
            'W1': 'work_meeting_baseline',
            'W2': 'work_meeting_many'
        }
        event_type_str = event_type_map.get(event_label, '')
        
        # Collect all matching directories and extract their timestamps
        # Format: multiseed_YYYYMMDD_HHMMSS_simtype_eventtype_seed
        timestamp_dirs = {}
        for item in os.listdir(results_base):
            item_path = os.path.join(results_base, item)
            if os.path.isdir(item_path) and 'multiseed' in item.lower() and event_type_str in item.lower():
                # Extract timestamp (format: multiseed_20260102_002425_...)
                parts = item.split('_')
                if len(parts) >= 3:
                    timestamp = f"{parts[1]}_{parts[2]}"  # YYYYMMDD_HHMMSS
                    if timestamp not in timestamp_dirs:
                        timestamp_dirs[timestamp] = []
                    timestamp_dirs[timestamp].append(item)
        
        if not timestamp_dirs:
            print(f"  No matching results found")
            continue
        
        # Get the latest timestamp
        latest_timestamp = max(timestamp_dirs.keys())
        print(f"  Using results from latest run: {latest_timestamp}")
        
        # Use only results from the latest run
        for item in timestamp_dirs[latest_timestamp]:
            item_path = os.path.join(results_base, item)
            # Extract simulation type and seed from directory name
            if 'memilio' in item.lower():
                sim_type = 'memilio'
            elif 'panvadere' in item.lower():
                sim_type = 'panvadere'
            else:
                continue
            
            # Extract seed
            seed_part = [part for part in item.split('_') if 'seed' in part.lower()]
            if seed_part:
                seed = seed_part[0].replace('seed', '')
                result_paths.append(item_path)
                result_labels.append(f"{sim_type}_seed{seed}")
        
        if not result_paths:
            print(f"Warning: No valid results found in {scenario_dir}")
            continue
        
        print(f"  Found {len(result_paths)} result sets")
        
        # Group results
        grouped_results = group_results_by_simulation_type(result_paths, result_labels)
        print(f"  Memilio results: {len(grouped_results['memilio'])}")
        print(f"  Panvadere results: {len(grouped_results['panvadere'])}")
        
        if grouped_results['memilio'] or grouped_results['panvadere']:
            scenario_data_list.append((event_label, grouped_results))

    if not scenario_data_list:
        print("Error: No valid scenario data found")
        return 1

    # Create combined visualization
    print(f"\nCreating combined visualization for {len(scenario_data_list)} scenarios...")
    create_combined_median_comparison(scenario_data_list, args.output_dir)
    
    print("\n✓ Combined median comparison completed successfully!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
