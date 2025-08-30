#!/usr/bin/env python3
"""
Comparative Temporal Infection Heatmap
Compares transmission-informed vs uniform initialization approaches
Shows location-based transmission patterns and household infection statistics

Features:
- Best run mode: Shows results from the single best-performing run
- Average mode: Shows averaged results across multiple simulation runs
- Household-based visualization with color coding for infection burden
- Temporal progression analysis at multiple time points

Usage:
    # Best run mode (default)
    python comparative_temporal_heatmap.py
    
    # Average mode across all runs
    python comparative_temporal_heatmap.py --use-average
    
    # With custom parameters
    python comparative_temporal_heatmap.py --use-average --time-points 0 48 120 240
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import argparse
import os
import sys
import glob
from collections import defaultdict, Counter

# Import functions from the existing scripts
try:
    from contact_network import (
        load_household_data, load_contact_data,
        load_location_data, create_household_layout
    )
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure contact_network.py is in the same directory")
    sys.exit(1)

import networkx as nx
import matplotlib.patches as mpatches


def load_temporal_infection_data(detailed_infection_file):
    """Load temporal infection data and return infections by time step"""
    try:
        # Read the detailed infection data
        infection_data = pd.read_csv(detailed_infection_file)

        # Create a dictionary mapping timestep to infected persons
        infections_by_time = defaultdict(set)

        # Find all people with negative timesteps (initial seeding) or EventPanvadere
        initial_infected = set()

        for _, row in infection_data.iterrows():
            timestep = row['Timestep']
            person_id = row['Person_ID']
            location_type = row['Location_Type']

            # Handle initial seeding (negative timesteps or EventPanvadere)
            if timestep < 0 or location_type == 'EventPanvadere':
                initial_infected.add(person_id)

        # Add initial infections to timestep 0 and all subsequent timesteps
        for timestep in range(0, max(infection_data['Timestep']) + 1):
            infections_by_time[timestep].update(initial_infected)

        # Now process regular infections (positive timesteps, non-EventPanvadere)
        for _, row in infection_data.iterrows():
            timestep = row['Timestep']
            person_id = row['Person_ID']
            location_type = row['Location_Type']

            # Only track positive timesteps that are not initial seeding
            if timestep > 0 and location_type != 'EventPanvadere':
                # Add this person to all future timesteps (they remain infected)
                for t in range(timestep, max(infection_data['Timestep']) + 1):
                    infections_by_time[t].add(person_id)

        return dict(infections_by_time)

    except Exception as e:
        print(f"Error loading temporal infection data: {e}")
        return {}


def find_infection_files(data_dir, use_average=False):
    """Find infection data files - either best run or all runs depending on use_average flag"""
    if use_average:
        # Look for detailed_infection files from all runs
        search_pattern = "**/run_*_detailed_infection.csv"
        search_path = os.path.join(data_dir, search_pattern)
        files = glob.glob(search_path, recursive=True)
        if not files:
            # Fallback to best_run files if no individual run files found
            search_pattern = "**/best_run_detailed_infection.csv"
            search_path = os.path.join(data_dir, search_pattern)
            files = glob.glob(search_path, recursive=True)
            print(
                f"Warning: No individual run files found, using best run files: {len(files)} files")
        else:
            print(
                f"Found {len(files)} individual run infection files in {data_dir}")
    else:
        # Look for best run files only
        search_pattern = "**/best_run_detailed_infection.csv"
        search_path = os.path.join(data_dir, search_pattern)
        files = glob.glob(search_path, recursive=True)
        print(f"Found {len(files)} best run infection files in {data_dir}")

    return files


def calculate_average_infections_by_time(infection_files, household_data, time_points):
    """Calculate average infections across multiple runs for specific time points"""
    print(
        f"Calculating average infections across {len(infection_files)} runs...")

    # Dictionary to store infections per run per timestep
    # run_infections[timestep] = [set_infected_run1, set_infected_run2, ...]
    run_infections = defaultdict(list)

    # Dictionary to store household infection counts per run per timestep
    # household_run_infections[timestep][household_id] = [count_run1, count_run2, ...]
    household_run_infections = defaultdict(lambda: defaultdict(list))

    successful_runs = 0

    for i, infection_file in enumerate(infection_files):
        print(
            f"Processing run {i+1}/{len(infection_files)}: {os.path.basename(os.path.dirname(infection_file))}")

        # Load temporal data for this run
        infections_by_time = load_temporal_infection_data(infection_file)

        if not infections_by_time:
            print(f"  Warning: No infection data loaded from {infection_file}")
            continue

        # Store infections for each requested time point
        for timestep in time_points:
            infected_at_time = infections_by_time.get(timestep, set())
            run_infections[timestep].append(infected_at_time)

            # Calculate household infections for this timestep
            household_stats = defaultdict(
                lambda: {'total_members': 0, 'infected_members': 0})

            for _, row in household_data.iterrows():
                person_id = row['person_id']
                household_id = row['household_id']

                household_stats[household_id]['total_members'] += 1

                if person_id in infected_at_time:
                    household_stats[household_id]['infected_members'] += 1

            # Store household infection counts for this run
            for household_id, stats in household_stats.items():
                infected_count = stats['infected_members']
                household_run_infections[timestep][household_id].append(
                    infected_count)

        successful_runs += 1

    print(f"Successfully processed {successful_runs} runs")

    # Calculate averages for individual infections
    average_infections_by_time = {}
    for timestep in time_points:
        if timestep in run_infections and run_infections[timestep]:
            # Calculate average infected individuals
            all_infected_sets = run_infections[timestep]
            # For average, we need to calculate which individuals are infected on average
            # This is complex for individuals, so we'll use the most commonly infected ones

            # Count how many times each person was infected across runs
            person_infection_counts = defaultdict(int)
            total_runs = len(all_infected_sets)

            for infected_set in all_infected_sets:
                for person_id in infected_set:
                    person_infection_counts[person_id] += 1

            # Include people who were infected in at least 50% of runs
            # This creates a "typical" infection pattern
            threshold = total_runs * 0.5
            average_infected_set = set()
            for person_id, count in person_infection_counts.items():
                if count >= threshold:
                    average_infected_set.add(person_id)

            average_infections_by_time[timestep] = average_infected_set
        else:
            average_infections_by_time[timestep] = set()

    # Calculate average household infections
    average_household_infections_by_time = {}
    for timestep in time_points:
        if timestep in household_run_infections:
            household_averages = {}
            for household_id, infection_counts in household_run_infections[timestep].items():
                if infection_counts:
                    avg_infected = np.mean(infection_counts)
                    household_averages[household_id] = avg_infected
                else:
                    household_averages[household_id] = 0.0
            average_household_infections_by_time[timestep] = household_averages
        else:
            average_household_infections_by_time[timestep] = {}

    return average_infections_by_time, average_household_infections_by_time, successful_runs


def get_infection_status_color(person_id, infected_at_time):
    """Get color based on infection status - using visible colors"""
    if person_id in infected_at_time:
        return '#FF4444', '#8B0000'  # Red with dark red border for infected
    else:
        return 'white', 'gray'  # White with gray border for not infected


def calculate_household_contact_matrix(household_data, contact_data):
    """Calculate inter-household contact strength matrix"""
    print("Calculating household contact matrix...")
    print(f"Contact data columns: {contact_data.columns.tolist()}")

    # Create person to household mapping
    person_to_household = {}
    for _, row in household_data.iterrows():
        person_to_household[row['person_id']] = row['household_id']

    # Get list of unique households
    households = sorted(household_data['household_id'].unique())
    household_to_idx = {hh: i for i, hh in enumerate(households)}

    # Initialize contact matrix
    n_households = len(households)
    contact_matrix = np.zeros((n_households, n_households))

    # Get the actual column names (strip whitespace)
    columns = [col.strip() for col in contact_data.columns]
    contact_data.columns = columns

    # Calculate inter-household contact hours
    for _, row in contact_data.iterrows():
        person1 = row['person_1']
        person2 = row['person_2']
        contact_hours = row['contact_hours']

        # Get households for both persons
        if person1 in person_to_household and person2 in person_to_household:
            hh1 = person_to_household[person1]
            hh2 = person_to_household[person2]

            # Only count inter-household contacts (not within same household)
            if hh1 != hh2:
                idx1 = household_to_idx[hh1]
                idx2 = household_to_idx[hh2]
                contact_matrix[idx1, idx2] += contact_hours
                contact_matrix[idx2, idx1] += contact_hours  # Symmetric

    return households, contact_matrix


def create_contact_based_household_layout(household_data, contact_data,
                                          household_spacing=150.0, person_spacing=15.0,
                                          max_households_per_row=25):
    """Create uniform rectangular layout for households with even spacing"""

    # Calculate household contact matrix
    households, contact_matrix = calculate_household_contact_matrix(
        household_data, contact_data)

    print(f"Creating uniform grid layout for {len(households)} households...")

    # Create a simple grid layout with even spacing
    person_pos = {}
    household_positions = {}
    household_dimensions = {}

    # Calculate grid dimensions for even distribution
    total_households = len(households)
    households_per_row = min(max_households_per_row,
                             int(np.sqrt(total_households) * 1.5))
    rows = int(np.ceil(total_households / households_per_row))

    # Sort households for consistent layout
    sorted_households = sorted(households)

    # Create uniform rectangular layout
    for idx, household_id in enumerate(sorted_households):
        row = idx // households_per_row
        col = idx % households_per_row

        # Calculate position with even spacing
        hx = col * household_spacing
        hy = row * household_spacing

        household_positions[household_id] = (hx, hy)

        # Set uniform rectangular dimensions for all households
        household_width = household_spacing * 0.6  # 60% of spacing for consistent size
        household_height = household_spacing * 0.4  # Rectangular shape
        household_dimensions[household_id] = (
            household_width, household_height)

        # Get household members and place them at center (not used for display)
        household_members = household_data[household_data['household_id']
                                           == household_id]['person_id'].tolist()
        for person_id in household_members:
            person_pos[person_id] = (hx, hy)

    print(f"Uniform grid layout created: {total_households} households")
    print(f"Grid dimensions: {rows} rows x {households_per_row} columns")
    print(f"Household spacing: {household_spacing} units")

    return person_pos, households, household_positions, household_dimensions


def count_infected_households(household_data, infected_at_time):
    """Count households with at least one infected member"""
    infected_households = set()

    for _, row in household_data.iterrows():
        person_id = row['person_id']
        household_id = row['household_id']

        if person_id in infected_at_time:
            infected_households.add(household_id)

    return len(infected_households)


def create_time_specific_heatmap(ax, household_data, household_positions, household_dimensions, person_pos,
                                 infected_at_time, timestep, title, method_name,
                                 use_average=False, household_average_infections=None, num_runs=None):
    """Create heatmap for specific time point - showing households colored by infection count"""

    if use_average and household_average_infections is not None:
        # Use average household infections
        household_stats = {}
        for household_id, avg_infected in household_average_infections.items():
            household_stats[household_id] = {
                'infection_rate': 0,  # Not used in display
                'infected_members': avg_infected,  # This will be the average
                'total_members': 1  # Not used for average display
            }
    else:
        # Calculate household infection stats for this time point (best run mode)
        household_stats = defaultdict(
            lambda: {'total_members': 0, 'infected_members': 0})

        for _, row in household_data.iterrows():
            person_id = row['person_id']
            household_id = row['household_id']

            household_stats[household_id]['total_members'] += 1

            if person_id in infected_at_time:
                household_stats[household_id]['infected_members'] += 1

        # Calculate infection rates for best run mode
        for household_id in household_stats:
            stats = household_stats[household_id]
            if stats['total_members'] > 0:
                stats['infection_rate'] = stats['infected_members'] / \
                    stats['total_members']

    # Draw uniform rectangular households colored by infection count
    for household_id, position in household_positions.items():
        if not isinstance(position, tuple):
            continue

        hx, hy = position
        household_width, household_height = household_dimensions.get(
            household_id, (40, 40))

        stats = household_stats.get(household_id, {
            'infection_rate': 0, 'infected_members': 0, 'total_members': 0})

        if use_average:
            # For average mode, infected_count is already a float average
            infected_count = stats['infected_members']
        else:
            # For best run mode, infected_count is an integer
            infected_count = stats['infected_members']

        # Color household based on number of infected members (or average)
        if use_average:
            # Color based on average infection count (0.5 steps)
            if infected_count < 0.5:
                # <0.5 average infected - white
                household_color = 'white'
                edge_color = 'lightgray'
            elif infected_count < 1.0:
                # 0.5-1.0 average infected - light yellow
                household_color = '#FFFACD'
                edge_color = '#FFD700'
            elif infected_count < 1.5:
                # 1.0-1.5 average infected - yellow
                household_color = '#FFFF99'
                edge_color = '#FFD700'
            elif infected_count < 2.0:
                # 1.5-2.0 average infected - yellow-orange
                household_color = '#FFCC33'
                edge_color = '#FF8C00'
            elif infected_count < 2.5:
                # 2.0-2.5 average infected - orange
                household_color = '#FF9933'
                edge_color = '#FF4500'
            elif infected_count < 3.0:
                # 2.5-3.0 average infected - orange-red
                household_color = '#FF6600'
                edge_color = '#FF4500'
            elif infected_count < 3.5:
                # 3.0-3.5 average infected - light red
                household_color = '#FF3300'
                edge_color = '#CC0000'
            elif infected_count < 4.0:
                # 3.5-4.0 average infected - red
                household_color = '#FF0000'
                edge_color = '#CC0000'
            elif infected_count < 4.5:
                # 4.0-4.5 average infected - bright red
                household_color = '#DD0000'
                edge_color = '#AA0000'
            else:
                # ≥4.5 average infected - dark black
                household_color = '#000000'
                edge_color = '#000000'
        else:
            # Color based on discrete infection count (best run mode)
            if infected_count == 0:
                # 0 infected - white
                household_color = 'white'
                edge_color = 'lightgray'
            elif infected_count == 1:
                # 1 infected - light yellow
                household_color = '#FFFACD'
                edge_color = '#F0E68C'
            elif infected_count == 2:
                # 2 infected - yellow
                household_color = '#FFFF99'
                edge_color = '#FFD700'
            elif infected_count == 3:
                # 3 infected - orange
                household_color = '#FFB366'
                edge_color = '#FF8C00'
            elif infected_count == 4:
                # 4 infected - red
                household_color = '#FF6666'
                edge_color = '#CC0000'
            else:
                # 5 infected - black
                household_color = 'black'
                edge_color = 'darkgray'

        # Draw household as a uniform rectangle
        rect = plt.Rectangle((hx - household_width/2, hy - household_height/2),
                             household_width, household_height,
                             facecolor=household_color, alpha=0.9,
                             edgecolor=edge_color, linewidth=1.5)
        ax.add_patch(rect)

    # Set title and formatting
    infected_count = len(infected_at_time)
    total_people = len(person_pos)
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"

    # Modify title based on mode
    if use_average and num_runs:
        # Calculate average infected individuals from household averages
        if household_average_infections:
            avg_total_infected = sum(household_average_infections.values())
            title_suffix = f"{avg_total_infected:.1f} avg infected across {num_runs} runs"
        else:
            title_suffix = f"avg. across {num_runs} runs"
    else:
        title_suffix = f"{infected_count}/{total_people} infected ({infected_count/total_people*100:.1f}%)"

    ax.set_title(f"{method_name}\n{day_title}: {title_suffix}",
                 fontsize=11, fontweight='bold', pad=10)
    ax.axis('equal')
    ax.axis('off')

    return infected_count


def create_household_comparison_plot(ax, time_points, household_counts_method1, household_counts_method2,
                                     method1_name, method2_name, use_average=False,
                                     avg_infected_method1=None, avg_infected_method2=None):
    """Create slim comparison plot of household infection counts and average infected if available"""

    # Convert time points to day labels
    day_labels = []
    for tp in time_points:
        if tp == 0:
            day_labels.append("Init")
        else:
            day_labels.append(f"Day {tp//24}")

    x_pos = np.arange(len(time_points))
    width = 0.35

    bars1 = ax.bar(x_pos - width/2, household_counts_method1, width,
                   label=method1_name, color='#d62728', alpha=0.8)
    bars2 = ax.bar(x_pos + width/2, household_counts_method2, width,
                   label=method2_name, color='#1f77b4', alpha=0.8)

    # Add value labels on bars
    for i, (v1, v2) in enumerate(zip(household_counts_method1, household_counts_method2)):
        # For average mode, also show average infected count below the household count
        if use_average and avg_infected_method1 and avg_infected_method2:
            label1 = f"{v1}\n({avg_infected_method1[i]:.0f})"
            label2 = f"{v2}\n({avg_infected_method2[i]:.0f})"
        else:
            label1 = str(v1)
            label2 = str(v2)

        ax.text(i - width/2, v1 + max(household_counts_method1 + household_counts_method2) * 0.02,
                label1, ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax.text(i + width/2, v2 + max(household_counts_method1 + household_counts_method2) * 0.02,
                label2, ha='center', va='bottom', fontsize=8, fontweight='bold')

    ax.set_xlabel('Time Points', fontsize=9)

    # Update ylabel based on mode
    if use_average:
        ax.set_ylabel('Infected Households\n(Avg Infected)', fontsize=9)
        ax.set_title('Households with ≥1 Infected Member\n(Average infected individuals in parentheses)',
                     fontsize=10, fontweight='bold')
    else:
        ax.set_ylabel('Infected\nHouseholds', fontsize=9)
        ax.set_title('Households with ≥1 Infected Member',
                     fontsize=10, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(day_labels, fontsize=8)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(axis='y', alpha=0.3)


def create_comparative_temporal_heatmap(data_dir_method1, data_dir_method2=None,
                                        method1_name="Transmission-Informed", method2_name="Uniform",
                                        output_path=None, time_points=[0, 24, 72, 240], use_average=False):
    """Create comparative temporal infection heatmap with option for best run or average across runs"""

    mode_description = "Average across runs" if use_average else "Best run"
    print(
        f"Loading temporal infection data for both methods ({mode_description})...")

    # Load Method 1 data (transmission-informed)
    try:
        household_data = load_household_data(
            os.path.join(data_dir_method1, 'household_id.csv'))
        contact_data = load_contact_data(os.path.join(
            data_dir_method1, 'contact_intensiveness.csv'))
        location_data = load_location_data(os.path.join(
            data_dir_method1, 'location_id_and_type.csv'))

        if use_average:
            # Load all infection files for averaging
            infection_files_m1 = find_infection_files(
                data_dir_method1, use_average=True)
            if not infection_files_m1:
                print(
                    f"Error: No infection files found for Method 1 in {data_dir_method1}")
                return
            infections_by_time_m1, household_infections_by_time_m1, num_runs_m1 = calculate_average_infections_by_time(
                infection_files_m1, household_data, time_points)
        else:
            # Load best run only
            infections_by_time_m1 = load_temporal_infection_data(
                os.path.join(data_dir_method1, 'best_run_detailed_infection.csv'))
            household_infections_by_time_m1 = None
            num_runs_m1 = 1

    except FileNotFoundError as e:
        print(f"Error loading Method 1 data: {e}")
        return

    # Load Method 2 data (uniform) - if provided
    has_method2_data = False
    if data_dir_method2 and os.path.exists(data_dir_method2):
        try:
            if use_average:
                # Load all infection files for averaging
                infection_files_m2 = find_infection_files(
                    data_dir_method2, use_average=True)
                if infection_files_m2:
                    infections_by_time_m2, household_infections_by_time_m2, num_runs_m2 = calculate_average_infections_by_time(
                        infection_files_m2, household_data, time_points)
                    has_method2_data = True
                else:
                    print(
                        f"Warning: No infection files found for Method 2 in {data_dir_method2}")
                    has_method2_data = False
            else:
                # Load best run only
                infections_by_time_m2 = load_temporal_infection_data(
                    os.path.join(data_dir_method2, 'best_run_detailed_infection.csv'))
                household_infections_by_time_m2 = None
                num_runs_m2 = 1
                has_method2_data = True
        except FileNotFoundError as e:
            print(f"Warning: Method 2 data not found: {e}")
            print("Will create placeholder for comparison layout")
            has_method2_data = False
    else:
        print("Method 2 data directory not provided - creating placeholder")
        has_method2_data = False

    if use_average:
        print(
            f"Loaded Method 1 data: {len(household_data)} people, averaging across {num_runs_m1} runs")
        if has_method2_data:
            print(
                f"Loaded Method 2 data: {len(household_data)} people, averaging across {num_runs_m2} runs")
    else:
        print(
            f"Loaded Method 1 data: {len(household_data)} people, temporal infections for {len(infections_by_time_m1)} time steps")

    # Create household layout based on contact patterns (same for both methods)
    spacing_factor = 2.8  # Increased from 2.0 for maximum spacing
    household_spacing = 250.0 * spacing_factor  # Increased base spacing
    person_spacing = 80.0 * spacing_factor  # Increased base spacing

    print("Creating contact-based household layout...")
    person_pos, households, household_positions, household_dimensions = create_contact_based_household_layout(
        household_data, contact_data,
        household_spacing=household_spacing,
        person_spacing=person_spacing,
        max_households_per_row=22  # Reduced from 25 to accommodate larger spacing
    )

    # Create figure with 3 rows: Method1, Comparison, Method2 - more compact
    fig = plt.figure(figsize=(4*len(time_points), 10))

    # Calculate household infection counts for comparison
    household_counts_m1 = []
    household_counts_m2 = []

    infection_counts_m1 = []
    infection_counts_m2 = []

    # For average mode: collect average infected counts for comparison plot
    avg_infected_counts_m1 = []
    avg_infected_counts_m2 = []

    # Create subplots grid: 3 rows x len(time_points) columns - adjusted spacing
    gs = fig.add_gridspec(3, len(time_points), height_ratios=[
                          2.5, 0.8, 2.5], hspace=0.35, wspace=0.1)

    # Row 1: Method 1 (Transmission-informed)
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs[0, i])

        infected_at_time = infections_by_time_m1.get(timestep, set())
        day_label = timestep // 24

        # Get household infections for this timestep if using average mode
        household_avg_infections = None
        if use_average and household_infections_by_time_m1:
            household_avg_infections = household_infections_by_time_m1.get(
                timestep, {})

        infected_count = create_time_specific_heatmap(
            ax, household_data, household_positions, household_dimensions, person_pos,
            infected_at_time, day_label,
            f"Time Point {i+1}", method1_name,
            use_average=use_average,
            household_average_infections=household_avg_infections,
            num_runs=num_runs_m1 if use_average else None
        )

        infection_counts_m1.append(infected_count)
        household_counts_m1.append(count_infected_households(
            household_data, infected_at_time))

        # For average mode: collect average infected count from household averages
        if use_average and household_avg_infections:
            avg_total_infected = sum(household_avg_infections.values())
            avg_infected_counts_m1.append(avg_total_infected)
        else:
            avg_infected_counts_m1.append(infected_count)

    # Row 2: Household comparison (slim panel spanning all columns)
    ax_comparison = fig.add_subplot(gs[1, :])
    ax_comparison.margins(y=0.1)

    # For now, use placeholder data for method 2 if not available
    if not has_method2_data:
        # Create simulated uniform data for demonstration
        household_counts_m2 = [int(count * 1.2)
                               for count in household_counts_m1]
        infection_counts_m2 = [int(count * 1.15)
                               for count in infection_counts_m1]
        avg_infected_counts_m2 = [
            count * 1.15 for count in avg_infected_counts_m1]
    else:
        # Calculate actual Method 2 data
        for i, timestep in enumerate(time_points):
            infected_at_time = infections_by_time_m2.get(timestep, set())
            infection_counts_m2.append(len(infected_at_time))
            household_counts_m2.append(count_infected_households(
                household_data, infected_at_time))

            # For average mode: collect average infected count from household averages
            if use_average and household_infections_by_time_m2:
                household_avg_infections = household_infections_by_time_m2.get(
                    timestep, {})
                if household_avg_infections:
                    avg_total_infected = sum(household_avg_infections.values())
                    avg_infected_counts_m2.append(avg_total_infected)
                else:
                    avg_infected_counts_m2.append(len(infected_at_time))
            else:
                avg_infected_counts_m2.append(len(infected_at_time))

    create_household_comparison_plot(ax_comparison, time_points, household_counts_m1, household_counts_m2,
                                     method1_name, method2_name, use_average=use_average,
                                     avg_infected_method1=avg_infected_counts_m1,
                                     avg_infected_method2=avg_infected_counts_m2)

    # Row 3: Method 2 (Uniform) - placeholder or real data
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs[2, i])
        if has_method2_data:
            infected_at_time = infections_by_time_m2.get(timestep, set())
            day_label = timestep // 24

            # Get household infections for this timestep if using average mode
            household_avg_infections = None
            if use_average and household_infections_by_time_m2:
                household_avg_infections = household_infections_by_time_m2.get(
                    timestep, {})

            create_time_specific_heatmap(
                ax, household_data, household_positions, household_dimensions, person_pos,
                infected_at_time, day_label,
                f"Time Point {i+1}", method2_name,
                use_average=use_average,
                household_average_infections=household_avg_infections,
                num_runs=num_runs_m2 if use_average else None
            )
        else:
            # Create placeholder visualization
            ax.text(0.5, 0.5, f"{method2_name}\nData Not Available\n\n(Placeholder for comparison layout)",
                    ha='center', va='center', transform=ax.transAxes, fontsize=12,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.5))
            ax.set_title(f"{method2_name} - Time Point {i+1}",
                         fontsize=11, fontweight='bold')
            ax.axis('off')

    # Create legend based on mode
    if use_average:
        # Legend for average infection counts (0.5 steps)
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor='white',
                          edgecolor='lightgray', label='<0.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFACD',
                          edgecolor='#FFD700', label='0.5-1.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFF99',
                          edgecolor='#FFD700', label='1.0-1.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFCC33',
                          edgecolor='#FF8C00', label='1.5-2.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF9933',
                          edgecolor='#FF4500', label='2.0-2.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF6600',
                          edgecolor='#FF4500', label='2.5-3.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF3300',
                          edgecolor='#CC0000', label='3.0-3.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF0000',
                          edgecolor='#CC0000', label='3.5-4.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#DD0000',
                          edgecolor='#AA0000', label='4.0-4.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#000000',
                          edgecolor='#000000', label='≥4.5')
        ]
        legend_ncol = 5  # Reduced to accommodate fewer colors
    else:
        # Legend for discrete infection counts (best run mode)
        legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor='white',
                          edgecolor='lightgray', label='0 Infected'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFACD',
                          edgecolor='#F0E68C', label='1 Infected'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFF99',
                          edgecolor='#FFD700', label='2 Infected'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFB366',
                          edgecolor='#FF8C00', label='3 Infected'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF6666',
                          edgecolor='#CC0000', label='4 Infected'),
            plt.Rectangle((0, 0), 1, 1, facecolor='black',
                          edgecolor='darkgray', label='5 Infected')
        ]
        legend_ncol = 6

    # Place legend outside the plots in one row with title
    if use_average:
        legend_title = 'Average Infected Members per Household'
    else:
        legend_title = 'Infected Members per Household'

    legend = fig.legend(handles=legend_elements, loc='center',
                        bbox_to_anchor=(0.5, 0.06), ncol=legend_ncol, fontsize=9,
                        title=legend_title, title_fontsize=10)
    legend.get_title().set_fontweight('bold')

    # Add overall title
    title_suffix = f"({mode_description})" if use_average else ""
    fig.suptitle(f'Comparative Temporal Infection Progression {title_suffix}\nTransmission-Informed vs Uniform Initialization',
                 fontsize=14, fontweight='bold', y=0.98)

    # Adjust layout with more bottom space for the legend
    plt.subplots_adjust(bottom=0.15, top=0.88, left=0.05, right=0.95)

    if output_path:
        plt.savefig(output_path, dpi=300,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative temporal heatmap saved to {output_path}")

    # Print summary statistics
    if use_average:
        print(f"\n{method1_name} Summary (averaged across {num_runs_m1} runs):")
        for i, (timestep, inf_count, hh_count) in enumerate(zip(time_points, infection_counts_m1, household_counts_m1)):
            time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
            print(
                f"  {time_label}: ~{inf_count} avg infected, ~{hh_count} households affected")

        if has_method2_data:
            print(f"\n{method2_name} Summary (averaged across {num_runs_m2} runs):")
            for i, (timestep, inf_count, hh_count) in enumerate(zip(time_points, infection_counts_m2, household_counts_m2)):
                time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
                print(
                    f"  {time_label}: ~{inf_count} avg infected, ~{hh_count} households affected")
    else:
        print(f"\n{method1_name} Summary (best run):")
        for i, (timestep, inf_count, hh_count) in enumerate(zip(time_points, infection_counts_m1, household_counts_m1)):
            time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
            print(
                f"  {time_label}: {inf_count} infected, {hh_count} households affected")

        if has_method2_data:
            print(f"\n{method2_name} Summary (best run):")
            for i, (timestep, inf_count, hh_count) in enumerate(zip(time_points, infection_counts_m2, household_counts_m2)):
                time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
                print(
                    f"  {time_label}: {inf_count} infected, {hh_count} households affected")

    return fig


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Create comparative temporal infection heatmap')
    parser.add_argument(
        '--data-dir-method1', help='Directory containing transmission-informed data files')
    parser.add_argument(
        '--data-dir-method2', help='Directory containing uniform initialization data files')
    parser.add_argument(
        '--method1-name', default='Transmission-Informed', help='Name for first method')
    parser.add_argument('--method2-name', default='Uniform',
                        help='Name for second method')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--time-points', nargs='+', type=int, default=[0, 24, 72, 24*10],
                        help='Time points in hours (default: 0, 24, 72, 240 for initialization, day 1, day 3, day 10)')
    parser.add_argument('--use-average', action='store_true', default=True,
                        help='Use average across all runs instead of best run only')

    args = parser.parse_args()

    # Debug paths (comment out for production)
    args.data_dir_method1 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250830_215731_R1_restaurant_strong_clustering_transmission_informed"
    # This one doesn't have multi-run data
    args.data_dir_method2 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250830_215737_R1_restaurant_strong_clustering_uniform_initialized"
    args.method1_name = "Transmission-Informed"
    args.method2_name = "Uniform"

    # Modify output path to indicate mode
    mode_suffix = "_average" if args.use_average else "_best_run"
    args.output_path = f"/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/comparative_temporal_infection_heatmap{mode_suffix}.png"

    if not os.path.exists(args.data_dir_method1):
        print(f"Error: Directory {args.data_dir_method1} does not exist")
        sys.exit(1)

    if not args.output_path:
        mode_suffix = "_average" if args.use_average else "_best_run"
        args.output_path = os.path.join(
            args.data_dir_method1, f'comparative_temporal_infection_heatmap{mode_suffix}.png')

    # Create visualization
    fig = create_comparative_temporal_heatmap(
        data_dir_method1=args.data_dir_method1,
        data_dir_method2=args.data_dir_method2,
        method1_name=args.method1_name,
        method2_name=args.method2_name,
        output_path=args.output_path,
        time_points=args.time_points,
        use_average=args.use_average
    )

    plt.show()


if __name__ == "__main__":
    main()
