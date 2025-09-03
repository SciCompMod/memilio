#!/usr/bin/env python3
"""
Multi-Seed Simulation Comparison Visualization

This script creates comprehensive visualizations comparing simulation results
across multiple seeds for both Memilio and Panvadere simulations.
"""

import sys
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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

    # Sum all infection states (excluding susceptible state 0)
    current_infections = data.sum(axis=1)
    return current_infections


def parse_label(label):
    """Parse label to extract simulation type and seed."""
    # Expected format: "memilio_seed1402121" or "panvadere_seed35897932"
    parts = label.split('_')
    if len(parts) >= 2:
        sim_type = parts[0]
        seed_part = parts[1]
        if seed_part.startswith('seed'):
            seed = seed_part[4:]  # Remove 'seed' prefix
            return sim_type, seed
    return label, "unknown"


def group_results_by_simulation_type(results_paths, labels):
    """Group results by simulation type and extract seed information."""
    grouped_results = {'memilio': [], 'panvadere': []}

    for path, label in zip(results_paths, labels):
        sim_type, seed = parse_label(label)

        # Load median results
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


def create_multi_seed_comparison_plot(grouped_results, output_dir, event_type, num_seeds):
    """Create comprehensive multi-seed comparison visualization."""

    plt.style.use('default')
    # Set larger font sizes globally
    plt.rcParams.update({'font.size': 18, 'axes.titlesize': 24, 'axes.labelsize': 20,
                         'xtick.labelsize': 16, 'ytick.labelsize': 16, 'legend.fontsize': 16})

    fig, axes = plt.subplots(2, 2, figsize=(20, 16))
    fig.suptitle(f'Scenario {event_type} - Multi-Seed Simulation',
                 fontsize=32, fontweight='bold', y=0.98)

    # Color schemes
    memilio_colors = plt.cm.Blues(np.linspace(
        0.3, 0.8, len(grouped_results['memilio'])))
    panvadere_colors = plt.cm.Reds(np.linspace(
        0.3, 0.8, len(grouped_results['panvadere'])))

    # Plot 1: All Memilio runs
    ax1 = axes[0, 0]
    ax1.set_title('Uniform - All Seeds', fontsize=24, fontweight='bold')

    memilio_finals = []
    for i, result in enumerate(grouped_results['memilio']):
        ax1.plot(result['time'], result['cumulative'],
                 color=memilio_colors[i], alpha=0.6, linewidth=2)
        memilio_finals.append(result['cumulative'][-1])

    ax1.set_xlabel('Days', fontsize=20)
    ax1.set_ylabel('Cumulative Infections', fontsize=20)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(labelsize=16)

    # Plot 2: All Panvadere runs
    ax2 = axes[0, 1]
    ax2.set_title('Transmission-Informed - All Seeds',
                  fontsize=24, fontweight='bold')

    panvadere_finals = []
    for i, result in enumerate(grouped_results['panvadere']):
        ax2.plot(result['time'], result['cumulative'],
                 color=panvadere_colors[i], alpha=0.6, linewidth=2)
        panvadere_finals.append(result['cumulative'][-1])

    ax2.set_xlabel('Days', fontsize=20)
    ax2.set_ylabel('Cumulative Infections', fontsize=20)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(labelsize=16)

    # Make y-axis limits the same for both plots
    if memilio_finals and panvadere_finals:
        # Get all cumulative values from both plots
        all_memilio_values = np.concatenate(
            [result['cumulative'] for result in grouped_results['memilio']])
        all_panvadere_values = np.concatenate(
            [result['cumulative'] for result in grouped_results['panvadere']])
        y_max = max(np.max(all_memilio_values), np.max(all_panvadere_values))
        y_min = 0  # Start from 0 for infections

        ax1.set_ylim(y_min, y_max * 1.05)  # Add 5% padding
        ax2.set_ylim(y_min, y_max * 1.05)

    # Plot 3: Median comparison with confidence intervals
    ax3 = axes[1, 0]
    ax3.set_title('Median Comparison with Confidence Intervals',
                  fontsize=24, fontweight='bold')

    if memilio_finals and panvadere_finals:
        # Calculate statistics for Memilio
        memilio_curves = np.array([result['cumulative']
                                  for result in grouped_results['memilio']])
        memilio_median = np.median(memilio_curves, axis=0)
        memilio_q25 = np.percentile(memilio_curves, 25, axis=0)
        memilio_q75 = np.percentile(memilio_curves, 75, axis=0)

        # Calculate statistics for Panvadere
        panvadere_curves = np.array([result['cumulative']
                                    for result in grouped_results['panvadere']])
        panvadere_median = np.median(panvadere_curves, axis=0)
        panvadere_q25 = np.percentile(panvadere_curves, 25, axis=0)
        panvadere_q75 = np.percentile(panvadere_curves, 75, axis=0)

        # Use time from first result (should be the same for all)
        time_axis = grouped_results['memilio'][0]['time'] if grouped_results[
            'memilio'] else grouped_results['panvadere'][0]['time']

        # Plot Memilio
        ax3.plot(time_axis, memilio_median, color='blue',
                 linewidth=4, label='Uniform (median)')
        ax3.fill_between(time_axis, memilio_q25, memilio_q75,
                         color='blue', alpha=0.3, label='Uniform (25-75%)')

        # Plot Panvadere
        ax3.plot(time_axis, panvadere_median, color='red',
                 linewidth=4, label='Transmission-Informed (median)')
        ax3.fill_between(time_axis, panvadere_q25, panvadere_q75,
                         color='red', alpha=0.3, label='Transmission-Informed (25-75%)')

    ax3.set_xlabel('Days', fontsize=22)
    ax3.set_ylabel('Cumulative Infections', fontsize=22)
    ax3.grid(True, alpha=0.3)
    ax3.legend(fontsize=16)
    ax3.tick_params(labelsize=16)

    ax4 = axes[1, 1]
    ax4.set_title('Final Infection Counts Distribution',
                  fontsize=24, fontweight='bold')

    if memilio_finals and panvadere_finals:
        # Create box plot
        box_data = [memilio_finals, panvadere_finals]
        box_labels = ['Memilio', 'Panvadere']

        bp = ax4.boxplot(box_data, tick_labels=box_labels, patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')

        # Make median lines thick and black
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(4)

        for i, (data, label) in enumerate(zip(box_data, box_labels)):
            y = data
            x = np.full(len(y), i+1)
            ax4.scatter(x, y, alpha=0.7, s=50,
                        color='darkblue' if i == 0 else 'darkred', zorder=10)  # High zorder to be on top

        # Add statistics text with larger font
        stats_text = f"Uniform: μ={np.mean(memilio_finals):.0f}, σ={np.std(memilio_finals):.0f}\n"
        stats_text += f"Transmission-Informed: μ={np.mean(panvadere_finals):.0f}, σ={np.std(panvadere_finals):.0f}"
        ax4.text(0.02, 0.1, stats_text, transform=ax4.transAxes, fontsize=18,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    ax4.set_ylabel('Final Cumulative Infections', fontsize=22)
    ax4.grid(True, alpha=0.3)
    ax4.tick_params(labelsize=16)

    plt.tight_layout()

    # Save the plot, name is multi_seed_comparison_{event_type}
    output_file = os.path.join(
        output_dir, f"multi_seed_comparison_{event_type}.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Multi-seed comparison plot saved to {output_file}")

    return output_file


def create_seed_variation_analysis(grouped_results, output_dir):
    """Create analysis of variation across seeds."""

    # Set larger font sizes
    plt.rcParams.update({'font.size': 18, 'axes.titlesize': 24, 'axes.labelsize': 20,
                         'xtick.labelsize': 16, 'ytick.labelsize': 16, 'legend.fontsize': 16})

    fig, axes = plt.subplots(2, 1, figsize=(14, 12))
    fig.suptitle('Seed Variation Analysis', fontsize=28, fontweight='bold')

    # Extract final infection counts
    memilio_finals = [result['cumulative'][-1]
                      for result in grouped_results['memilio']]
    memilio_seeds = [result['seed'] for result in grouped_results['memilio']]

    panvadere_finals = [result['cumulative'][-1]
                        for result in grouped_results['panvadere']]
    panvadere_seeds = [result['seed']
                       for result in grouped_results['panvadere']]

    # Plot 1: Final infection counts by seed
    ax1 = axes[0]
    ax1.set_title('Final Infection Counts by Seed',
                  fontsize=24, fontweight='bold')

    if memilio_finals:
        ax1.scatter(range(len(memilio_seeds)), memilio_finals,
                    color='blue', s=80, alpha=0.7, label='Memilio')
    if panvadere_finals:
        ax1.scatter(range(len(panvadere_seeds)), panvadere_finals,
                    color='red', s=80, alpha=0.7, label='Panvadere')

    ax1.set_xlabel('Seed Index', fontsize=20)
    ax1.set_ylabel('Final Cumulative Infections', fontsize=20)
    ax1.legend(fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(labelsize=16)

    # Plot 2: Coefficient of variation over time
    ax2 = axes[1]
    ax2.set_title('Coefficient of Variation Over Time',
                  fontsize=24, fontweight='bold')

    if grouped_results['memilio']:
        # Calculate CV for Memilio
        memilio_curves = np.array([result['cumulative']
                                  for result in grouped_results['memilio']])
        memilio_mean = np.mean(memilio_curves, axis=0)
        memilio_std = np.std(memilio_curves, axis=0)
        # Avoid division by zero
        memilio_cv = memilio_std / (memilio_mean + 1e-10)

        time_axis = grouped_results['memilio'][0]['time']
        ax2.plot(time_axis, memilio_cv, color='blue',
                 linewidth=3, label='Memilio CV')

    if grouped_results['panvadere']:
        # Calculate CV for Panvadere
        panvadere_curves = np.array([result['cumulative']
                                    for result in grouped_results['panvadere']])
        panvadere_mean = np.mean(panvadere_curves, axis=0)
        panvadere_std = np.std(panvadere_curves, axis=0)
        panvadere_cv = panvadere_std / \
            (panvadere_mean + 1e-10)  # Avoid division by zero

        time_axis = grouped_results['panvadere'][0]['time']
        ax2.plot(time_axis, panvadere_cv, color='red',
                 linewidth=3, label='Panvadere CV')

    ax2.set_xlabel('Days', fontsize=20)
    ax2.set_ylabel('Coefficient of Variation', fontsize=20)
    ax2.legend(fontsize=16)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(labelsize=16)

    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(output_dir, "seed_variation_analysis.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Seed variation analysis plot saved to {output_file}")

    return output_file


def create_summary_statistics(grouped_results, output_dir, event_type, num_seeds):
    """Create summary statistics report."""

    # Prepare statistics
    stats = []

    for sim_type in ['memilio', 'panvadere']:
        if grouped_results[sim_type]:
            finals = [result['cumulative'][-1]
                      for result in grouped_results[sim_type]]
            seeds = [result['seed'] for result in grouped_results[sim_type]]

            stats.append({
                'Simulation Type': sim_type.capitalize(),
                'Number of Seeds': len(finals),
                'Mean Final Infections': f"{np.mean(finals):.1f}",
                'Std Final Infections': f"{np.std(finals):.1f}",
                'Min Final Infections': f"{np.min(finals):.1f}",
                'Max Final Infections': f"{np.max(finals):.1f}",
                'CV Final Infections': f"{np.std(finals)/np.mean(finals):.3f}",
                'Seeds': ', '.join(seeds[:5]) + ('...' if len(seeds) > 5 else '')
            })

    # Create DataFrame and save to CSV
    df = pd.DataFrame(stats)
    csv_file = os.path.join(output_dir, "multi_seed_statistics.csv")
    df.to_csv(csv_file, index=False)
    print(f"Summary statistics saved to {csv_file}")

    # Create text report
    txt_file = os.path.join(output_dir, "multi_seed_report.txt")
    with open(txt_file, 'w') as f:
        f.write(f"Multi-Seed Simulation Analysis Report\n")
        f.write(f"=====================================\n\n")
        f.write(f"Event Type: {event_type}\n")
        f.write(f"Total Seeds Analyzed: {num_seeds}\n")
        f.write(
            f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("Summary Statistics:\n")
        f.write("-" * 50 + "\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")

        # Add interpretation
        if len(stats) == 2:
            memilio_mean = float(stats[0]['Mean Final Infections'])
            panvadere_mean = float(stats[1]['Mean Final Infections'])
            diff_pct = abs(memilio_mean - panvadere_mean) / \
                max(memilio_mean, panvadere_mean) * 100

            f.write("Key Findings:\n")
            f.write("-" * 20 + "\n")
            f.write(
                f"• Mean difference between simulations: {diff_pct:.1f}%\n")
            f.write(
                f"• {'Memilio' if memilio_mean > panvadere_mean else 'Panvadere'} shows higher infection rates on average\n")

            memilio_cv = float(stats[0]['CV Final Infections'])
            panvadere_cv = float(stats[1]['CV Final Infections'])
            f.write(
                f"• {'Memilio' if memilio_cv < panvadere_cv else 'Panvadere'} shows more consistent results across seeds\n")

    print(f"Detailed report saved to {txt_file}")

    return csv_file, txt_file


def main():
    """Main function for CLI usage."""
    parser = argparse.ArgumentParser(
        description="Multi-seed simulation comparison visualization")
    parser.add_argument("--results-paths", nargs='+', required=True,
                        help="Paths to simulation results directories")
    parser.add_argument("--labels", nargs='+', required=True,
                        help="Labels for each result (format: simtype_seedN)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for visualizations")
    parser.add_argument("--event-type", default="unknown",
                        help="Event type being analyzed")
    parser.add_argument("--num-seeds", type=int, default=20,
                        help="Number of seeds used")

    args = parser.parse_args()

    # Validate inputs
    if len(args.results_paths) != len(args.labels):
        print("Error: Number of results paths must match number of labels")
        return 1

    # Create output directory with fallback options
    original_output_dir = args.output_dir
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"✓ Output directory created: {args.output_dir}")
    except PermissionError:
        print(f"Warning: Permission denied creating {args.output_dir}")
        # Fallback 1: Try creating in user's home directory
        fallback_dir = os.path.expanduser(
            "~/multi_seed_analysis_" + datetime.now().strftime('%Y%m%d_%H%M%S'))
        try:
            os.makedirs(fallback_dir, exist_ok=True)
            args.output_dir = fallback_dir
            print(f"✓ Using fallback directory: {args.output_dir}")
        except Exception as e:
            # Fallback 2: Use current directory
            args.output_dir = os.path.join(
                os.getcwd(), "multi_seed_analysis_" + datetime.now().strftime('%Y%m%d_%H%M%S'))
            os.makedirs(args.output_dir, exist_ok=True)
            print(f"✓ Using current directory fallback: {args.output_dir}")
    except Exception as e:
        print(f"Error creating output directory: {e}")
        return 1

    print(f"Processing {len(args.results_paths)} simulation results...")
    print(f"Output directory: {args.output_dir}")

    # Group results by simulation type
    grouped_results = group_results_by_simulation_type(
        args.results_paths, args.labels)

    print(f"Found {len(grouped_results['memilio'])} Memilio results")
    print(f"Found {len(grouped_results['panvadere'])} Panvadere results")

    if not grouped_results['memilio'] and not grouped_results['panvadere']:
        print("Error: No valid results found to visualize")
        return 1

    # Create visualizations
    try:
        # Main comparison plot
        create_multi_seed_comparison_plot(grouped_results, args.output_dir,
                                          args.event_type, args.num_seeds)

        # Variation analysis
        create_seed_variation_analysis(grouped_results, args.output_dir)

        # Summary statistics
        create_summary_statistics(grouped_results, args.output_dir,
                                  args.event_type, args.num_seeds)

        print("✓ All visualizations completed successfully!")
        return 0

    except Exception as e:
        print(f"Error creating visualizations: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
