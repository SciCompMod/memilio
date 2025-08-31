#!/usr/bin/env python3
"""
Location Infection Distribution Pie Chart Generator

This script analyzes infection data across all runs to create pie charts showing
the average distribution of infections by location type.

Usage:
    python location_infection_pie_chart.py --data-dir <path> [options]

Author: Generated for epidemic simulation visualization
Date: August 2025
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import glob
from pathlib import Path
from collections import defaultdict
import seaborn as sns


def load_run_infections(run_file_path):
    """
    Load infection data from a single run file.

    Args:
        run_file_path (str): Path to the run's detailed infection CSV file

    Returns:
        pd.DataFrame: DataFrame with infection data
    """
    try:
        df = pd.read_csv(run_file_path)
        return df
    except Exception as e:
        print(f"Warning: Could not load {run_file_path}: {e}")
        return pd.DataFrame()


def analyze_location_infections(data_dir, exclude_negative_time=True, exclude_event_panvadere=True):
    """
    Analyze infection distributions by location type across all runs.

    Args:
        data_dir (str): Directory containing all_runs_detailed_infections folder
        exclude_negative_time (bool): Whether to exclude infections with negative timesteps
        exclude_event_panvadere (bool): Whether to exclude EventPanvadere locations

    Returns:
        dict: Location type distribution statistics
    """
    infections_dir = os.path.join(data_dir, "all_runs_detailed_infections")

    if not os.path.exists(infections_dir):
        raise FileNotFoundError(f"Directory not found: {infections_dir}")

    # Find all run files
    run_files = glob.glob(os.path.join(
        infections_dir, "run_*_detailed_infection.csv"))

    if not run_files:
        raise FileNotFoundError(f"No run files found in {infections_dir}")

    print(f"Found {len(run_files)} run files")

    # Dictionary to store location type counts for each run
    run_location_counts = {}
    all_location_types = set()

    for i, run_file in enumerate(sorted(run_files)):
        if i % 10 == 0:
            print(f"Processing run {i}...")

        df = load_run_infections(run_file)

        if df.empty:
            continue

        # Filter data if requested
        if exclude_negative_time:
            df = df[df['Timestep'] >= 0]

        if exclude_event_panvadere:
            df = df[df['Location_Type'] != 'EventPanvadere']

        # Count infections by location type for this run
        location_counts = df['Location_Type'].value_counts()

        # Extract run number from filename
        run_num = int(os.path.basename(run_file).split('_')[1])
        run_location_counts[run_num] = location_counts.to_dict()

        # Keep track of all location types we've seen
        all_location_types.update(location_counts.index)

    # Calculate statistics across all runs
    location_stats = {}

    for location_type in all_location_types:
        counts = []
        for run_counts in run_location_counts.values():
            counts.append(run_counts.get(location_type, 0))

        location_stats[location_type] = {
            'mean': np.mean(counts),
            'std': np.std(counts),
            'median': np.median(counts),
            'min': np.min(counts),
            'max': np.max(counts),
            'total_runs': len(counts),
            'runs_with_infections': sum(1 for c in counts if c > 0)
        }

    return location_stats, run_location_counts


def create_pie_chart(location_stats, output_path, scenario_name="", title_suffix=""):
    """
    Create a pie chart showing average infection distribution by location type.

    Args:
        location_stats (dict): Location type statistics
        output_path (str): Path for output file
        scenario_name (str): Name of the scenario for labeling
        title_suffix (str): Additional text for title
    """
    # Extract mean values and sort by size
    location_means = {loc: stats['mean']
                      for loc, stats in location_stats.items()}
    sorted_locations = sorted(location_means.items(),
                              key=lambda x: x[1], reverse=True)

    # Prepare data for pie chart
    labels = [loc for loc, _ in sorted_locations]
    sizes = [count for _, count in sorted_locations]
    total_infections = sum(sizes)

    # Calculate percentages
    percentages = [size/total_infections * 100 for size in sizes]

    # Define colors for different location types
    color_map = {
        'Home': '#FF9999',           # Light red
        'Work': '#66B2FF',           # Light blue
        'School': '#99FF99',         # Light green
        'BasicsShop': '#FFCC99',     # Light orange
        'SocialEvent': '#FF99CC',    # Light pink
        'Restaurant': '#FFD700',     # Gold
        'Hospital': '#FF6666',       # Red
        'ICU': '#CC0000',           # Dark red
        'EventPanvadere': '#CCCCCC',  # Gray
        'Unknown': '#999999'         # Dark gray
    }

    colors = [color_map.get(label, plt.cm.Set3(i/len(labels)))
              for i, label in enumerate(labels)]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    # Create pie chart
    wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                                      startangle=90, textprops={'fontsize': 24})

    # Customize the appearance
    plt.setp(autotexts, size=20, weight="bold")

    # Create title
    title = f"Average Infection Distribution by Location Type"
    if scenario_name:
        title += f"\n{scenario_name}"
    if title_suffix:
        title += f" - {title_suffix}"

    ax.set_title(title, fontsize=26, fontweight='bold', pad=20)

    # Add summary statistics as text
    stats_text = f"Total Infections (avg): {total_infections:.1f}\n"
    stats_text += f"Number of Location Types: {len(labels)}\n"
    stats_text += f"Runs Analyzed: {next(iter(location_stats.values()))['total_runs']}"

    # Add statistics box
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(1.15, 0.95, stats_text, transform=ax.transAxes, fontsize=16,
            verticalalignment='top', bbox=props)

    # Add detailed breakdown
    breakdown_text = "Detailed Breakdown:\n"
    for label, size, pct in zip(labels, sizes, percentages):
        std = location_stats[label]['std']
        breakdown_text += f"{label}: {size:.1f} ± {std:.1f} ({pct:.1f}%)\n"

    ax.text(1.15, 0.75, breakdown_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    # Adjust layout to prevent text cutoff
    plt.tight_layout()
    plt.subplots_adjust(right=0.75)

    # Save with high DPI
    plt.savefig(output_path, dpi=500, bbox_inches='tight', facecolor='white')
    print(f"Pie chart saved to: {output_path}")

    return fig


def create_comparative_pie_charts(data_dirs, scenario_names, output_path):
    """
    Create side-by-side pie charts for comparing scenarios.

    Args:
        data_dirs (list): List of data directories
        scenario_names (list): List of scenario names
        output_path (str): Path for output file
    """
    if len(data_dirs) != 2 or len(scenario_names) != 2:
        raise ValueError("Comparative charts require exactly 2 scenarios")

    # Analyze both scenarios
    stats_list = []
    for data_dir in data_dirs:
        stats, _ = analyze_location_infections(data_dir)
        stats_list.append(stats)

    # Create figure with subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Common color scheme for consistency
    all_locations = set()
    for stats in stats_list:
        all_locations.update(stats.keys())

    color_map = {
        'Home': '#FF9999',           # Light red
        'Work': '#66B2FF',           # Light blue
        'School': '#99FF99',         # Light green
        'BasicsShop': '#FFCC99',     # Light orange
        'SocialEvent': '#FF99CC',    # Light pink
        'Restaurant': '#FFD700',     # Gold
        'Hospital': '#FF6666',       # Red
        'ICU': '#CC0000',           # Dark red
        'EventPanvadere': '#CCCCCC',  # Gray
        'Unknown': '#999999'         # Dark gray
    }

    axes = [ax1, ax2]

    for i, (stats, scenario_name, ax) in enumerate(zip(stats_list, scenario_names, axes)):
        # Prepare data
        location_means = {loc: stat['mean'] for loc, stat in stats.items()}
        sorted_locations = sorted(
            location_means.items(), key=lambda x: x[1], reverse=True)

        labels = [loc for loc, _ in sorted_locations]
        sizes = [count for _, count in sorted_locations]
        colors = [color_map.get(label, plt.cm.Set3(j/len(all_locations)))
                  for j, label in enumerate(labels)]

        # Create pie chart
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, colors=colors,
                                          autopct='%1.1f%%', startangle=90,
                                          textprops={'fontsize': 20})

        plt.setp(autotexts, size=18, weight="bold")
        ax.set_title(f"{scenario_name}\nAvg Total: {sum(sizes):.1f} infections",
                     fontsize=22, fontweight='bold')

    plt.suptitle("Infection Distribution by Location Type - Comparison",
                 fontsize=28, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=500, bbox_inches='tight', facecolor='white')
    print(f"Comparative pie charts saved to: {output_path}")

    return fig


def main():
    parser = argparse.ArgumentParser(
        description="Generate pie charts for infection distribution by location type")
    parser.add_argument("--data-dir", required=True,
                        help="Directory containing all_runs_detailed_infections folder")
    parser.add_argument("--output",
                        help="Output file path (default: location_infection_distribution.png)")
    parser.add_argument("--scenario-name", default="",
                        help="Scenario name for labeling")
    parser.add_argument("--title-suffix", default="",
                        help="Additional text for title")
    parser.add_argument("--include-negative-time", action="store_true",
                        help="Include infections with negative timesteps")
    parser.add_argument("--include-event-panvadere", action="store_true",
                        help="Include EventPanvadere locations")
    parser.add_argument("--comparative", action="store_true",
                        help="Create comparative charts (requires --data-dir-2 and --scenario-name-2)")
    parser.add_argument("--data-dir-2",
                        help="Second data directory for comparison")
    parser.add_argument("--scenario-name-2", default="",
                        help="Second scenario name for comparison")

    args = parser.parse_args()

    # Set default output path
    if not args.output:
        if args.comparative:
            args.output = "comparative_location_infection_distribution.png"
        else:
            scenario_part = f"_{args.scenario_name}" if args.scenario_name else ""
            args.output = f"location_infection_distribution{scenario_part}.png"

    try:
        if args.comparative:
            if not args.data_dir_2:
                raise ValueError(
                    "--data-dir-2 is required for comparative charts")

            create_comparative_pie_charts(
                [args.data_dir, args.data_dir_2],
                [args.scenario_name or "Scenario 1",
                    args.scenario_name_2 or "Scenario 2"],
                args.output
            )
        else:
            # Analyze single scenario
            print(f"Analyzing infection data from: {args.data_dir}")
            location_stats, run_counts = analyze_location_infections(
                args.data_dir,
                exclude_negative_time=not args.include_negative_time,
                exclude_event_panvadere=not args.include_event_panvadere
            )

            print("\nLocation Type Statistics:")
            print("-" * 50)
            for location, stats in sorted(location_stats.items(),
                                          key=lambda x: x[1]['mean'], reverse=True):
                print(f"{location:15} | Mean: {stats['mean']:6.1f} | "
                      f"Std: {stats['std']:5.1f} | "
                      f"Runs: {stats['runs_with_infections']}/{stats['total_runs']}")

            # Create pie chart
            create_pie_chart(location_stats, args.output,
                             args.scenario_name, args.title_suffix)

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
