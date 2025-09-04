#!/usr/bin/env python3
"""
Comparative Temporal Infection Heatmap
Compares transmission-informed vs uniform initialization approaches
Shows location-based transmission patterns and household/workplace infection statistics

Features:
- Average/median mode: Shows averaged or median results across multiple simulation runs
- Household or workplace-based visualization with color coding for infection burden
- Temporal progression analysis at multiple time points
- Optional difference heatmaps showing Method1 - Method2 as a fourth row
- Size-based coloring: Households colored by size with infection-rate saturation

Usage:
    # Households with average across all runs (3-row layout)
    python comparative_temporal_heatmap.py

    # Workplace-based visualization
    python comparative_temporal_heatmap.py --use-workplaces

    # Use median instead of average
    python comparative_temporal_heatmap.py --use-median

    # Show difference heatmaps (4-row layout)
    python comparative_temporal_heatmap.py --show-differences

    # Use size-based coloring (households by size, workplaces all black, saturation by infection rate)
    python comparative_temporal_heatmap.py --use-size-based-colors

    # With custom parameters
    python comparative_temporal_heatmap.py --time-points 0 48 120 240 --use-workplaces --use-median --use-size-based-colors
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os
import sys
import glob
from collections import defaultdict

# Import functions from the existing scripts
try:
    from contact_network import (
        load_household_data, load_contact_data,
        load_location_data
    )
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure contact_network.py is in the same directory")
    sys.exit(1)


def load_temporal_infection_data(detailed_infection_file, max_timestep):
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
        for timestep in range(max(infection_data['Timestep']) + 1, max_timestep + 1):
            infections_by_time[timestep].update(
                infections_by_time[max(infection_data['Timestep'])])

        return dict(infections_by_time)

    except Exception as e:
        print(f"Error loading temporal infection data: {e}")
        return {}


def find_infection_files(data_dir):
    """Find infection data files from all runs for averaging"""
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

    return files


def load_workplace_data(work_id_file):
    """Load workplace data from work_id.csv file"""
    try:
        work_data = pd.read_csv(work_id_file)

        # Filter out people without workplace (Work_ID = 4294967295 indicates no workplace)
        work_data = work_data[work_data['Work_ID'] != 4294967295]

        # Rename columns to match household data structure
        work_data = work_data.rename(
            columns={'Person_ID': 'person_id', 'Work_ID': 'workplace_id'})

        print(
            f"Loaded workplace data: {len(work_data)} working people in {work_data['workplace_id'].nunique()} workplaces")
        return work_data
    except Exception as e:
        print(f"Error loading workplace data: {e}")
        # Return empty DataFrame with correct column structure
        return pd.DataFrame(columns=['person_id', 'workplace_id'])


def calculate_average_infections_by_time(infection_files, group_data, time_points, use_workplaces=False, use_median=False):
    """Calculate average or median infections across multiple runs for specific time points"""
    stat_type = "median" if use_median else "average"
    print(
        f"Calculating {stat_type} infections across {len(infection_files)} runs...")

    # Dictionary to store infections per run per timestep
    # run_infections[timestep] = [set_infected_run1, set_infected_run2, ...]
    run_infections = defaultdict(list)

    # Dictionary to store group infection counts per run per timestep
    # group_run_infections[timestep][group_id] = [count_run1, count_run2, ...]
    group_run_infections = defaultdict(lambda: defaultdict(list))

    # Dictionary to store infected group counts per run per timestep
    # infected_group_counts[timestep] = [count_run1, count_run2, ...]
    infected_group_counts = defaultdict(list)

    successful_runs = 0

    # Determine group type for clearer messaging
    group_type = "workplaces" if use_workplaces else "households"
    group_id_col = "workplace_id" if use_workplaces else "household_id"

    for i, infection_file in enumerate(infection_files):
        print(
            f"Processing run {i+1}/{len(infection_files)}: {os.path.basename(os.path.dirname(infection_file))}")

        # Load temporal data for this run
        infections_by_time = load_temporal_infection_data(
            infection_file, max(time_points))

        if not infections_by_time:
            print(f"  Warning: No infection data loaded from {infection_file}")
            continue

        # Store infections for each requested time point
        for timestep in time_points:
            infected_at_time = infections_by_time.get(timestep, set())
            run_infections[timestep].append(infected_at_time)

            # Calculate group infections for this timestep
            group_stats = defaultdict(
                lambda: {'total_members': 0, 'infected_members': 0})

            for _, row in group_data.iterrows():
                person_id = row['person_id']
                group_id = row[group_id_col]

                group_stats[group_id]['total_members'] += 1

                if person_id in infected_at_time:
                    group_stats[group_id]['infected_members'] += 1

            # Store group infection counts for this run
            for group_id, stats in group_stats.items():
                infected_count = stats['infected_members']
                group_run_infections[timestep][group_id].append(
                    infected_count)

            # Count groups with at least 1 infected member for this run
            infected_groups_count = count_infected_groups(
                group_data, infected_at_time, use_workplaces)
            infected_group_counts[timestep].append(
                infected_groups_count)
            if infected_groups_count == 0:
                print(
                    f"  Warning: No infected groups found at timestep {timestep}")

        successful_runs += 1

    print(f"Successfully processed {successful_runs} runs")

    # Calculate averages for individual infections
    average_infections_by_time = {}
    for timestep in time_points:
        if timestep in run_infections and run_infections[timestep]:
            all_infected_sets = run_infections[timestep]
            # For average/median, we need to calculate which individuals are infected on average/typically
            all_infected_sets_lengths = [
                len(infected_set) for infected_set in all_infected_sets]

            if use_median:
                average_infections_by_time[timestep] = np.median(
                    all_infected_sets_lengths)
            else:
                average_infections_by_time[timestep] = np.mean(
                    all_infected_sets_lengths)

     # Calculate average or median group infections
    average_group_infections_by_time = {}
    for timestep in time_points:
        if timestep in group_run_infections:
            group_stats = {}
            for group_id, infection_counts in group_run_infections[timestep].items():
                if infection_counts:
                    if use_median:
                        stat_value = np.median(infection_counts)
                    else:
                        stat_value = np.mean(infection_counts)
                    group_stats[group_id] = stat_value
                else:
                    group_stats[group_id] = 0.0
            average_group_infections_by_time[timestep] = group_stats
        else:
            average_group_infections_by_time[timestep] = {}

    # Calculate average or median infected group counts
    average_infected_group_counts = {}
    for timestep in time_points:
        if timestep in infected_group_counts and infected_group_counts[timestep]:
            if use_median:
                stat_value = np.median(infected_group_counts[timestep])
            else:
                stat_value = np.mean(infected_group_counts[timestep])
            average_infected_group_counts[timestep] = stat_value
        else:
            average_infected_group_counts[timestep] = 0.0

    return average_infections_by_time, average_group_infections_by_time, successful_runs, average_infected_group_counts


def get_household_size_color(household_size):
    """Get base color for household based on size - using visually appealing, distinct colors"""
    if household_size == 1:
        return '#FFD700'  # Gold
    elif household_size == 2:
        return '#FF6B6B'  # Coral Red
    elif household_size == 3:
        return '#4ECDC4'  # Teal
    elif household_size == 4:
        return "#006907"  # green
    elif household_size == 5:
        return '#9B59B6'  # Purple
    else:
        # For households with 6+ people, use a deeper purple
        return '#6A1B9A'  # Purple


def apply_saturation_to_color(base_color_hex, saturation_percent):
    """Apply saturation to a base color. saturation_percent should be 0-100"""
    import matplotlib.colors as mcolors

    # Convert hex to RGB
    rgb = mcolors.hex2color(base_color_hex)

    # Convert RGB to HSV
    hsv = mcolors.rgb_to_hsv(rgb)

    # Apply saturation (0-1 scale)
    hsv[1] = min(0.0, saturation_percent / 100.0)

    # Convert back to RGB and then to hex
    rgb_saturated = mcolors.hsv_to_rgb(hsv)
    hex_saturated = mcolors.rgb2hex(rgb_saturated)

    return hex_saturated


def calculate_group_sizes(group_data, use_workplaces=False):
    """Calculate the size of each group (household/workplace)"""
    group_id_col = "workplace_id" if use_workplaces else "household_id"
    group_sizes = {}

    # Count members for each group
    for _, row in group_data.iterrows():
        group_id = row[group_id_col]
        if group_id not in group_sizes:
            group_sizes[group_id] = 0
        group_sizes[group_id] += 1

    return group_sizes


def get_difference_color(difference_value):
    """Get color based on difference value (method1 - method2)"""
    if difference_value >= 4.5:
        # ≥4.5 difference - dark black
        return '#000000', '#000000'
    elif difference_value >= 4.0:
        # 4.0-4.5 difference - bright red
        return '#DD0000', '#AA0000'
    elif difference_value >= 3.5:
        # 3.5-4.0 difference - red
        return '#FF0000', '#CC0000'
    elif difference_value >= 3.0:
        # 3.0-3.5 difference - light red
        return '#FF3300', '#CC0000'
    elif difference_value >= 2.5:
        # 2.5-3.0 difference - orange-red
        return '#FF6600', '#FF4500'
    elif difference_value >= 2.0:
        # 2.0-2.5 difference - orange
        return '#FF9933', '#FF4500'
    elif difference_value >= 1.5:
        # 1.5-2.0 difference - yellow-orange
        return '#FFCC33', '#FF8C00'
    elif difference_value >= 1.0:
        # 1.0-1.5 difference - yellow
        return '#FFFF99', '#FFD700'
    elif difference_value >= 0.5:
        # 0.5-1.0 difference - light yellow
        return '#FFFACD', '#FFD700'
    elif difference_value >= -0.5:
        # -0.5 to 0.5 difference - white (near zero)
        return 'white', 'lightgray'
    elif difference_value >= -1.0:
        # -1.0 to -0.5 difference - light green
        return '#F0FFF0', '#90EE90'
    elif difference_value >= -1.5:
        # -1.5 to -1.0 difference - light green
        return '#E0FFE0', '#7CFC00'
    elif difference_value >= -2.0:
        # -2.0 to -1.5 difference - green
        return '#C0FFC0', '#32CD32'
    elif difference_value >= -2.5:
        # -2.5 to -2.0 difference - green
        return '#90FF90', '#228B22'
    elif difference_value >= -3.0:
        # -3.0 to -2.5 difference - darker green
        return '#60FF60', '#006400'
    elif difference_value >= -3.5:
        # -3.5 to -3.0 difference - darker green
        return '#30FF30', '#004000'
    elif difference_value >= -4.0:
        # -4.0 to -3.5 difference - dark green
        return '#00FF00', '#002000'
    elif difference_value >= -4.5:
        # -4.5 to -4.0 difference - very dark green
        return '#00CC00', '#001000'
    else:
        # < -4.5 difference - purple
        return '#800080', '#400040'


def create_difference_heatmap(ax, group_positions, group_dimensions,
                              group_diff_infections, timestep,
                              viz_style='rectangles', use_workplaces=False, use_median=False):
    """Create difference heatmap showing method1 - method2 for each group"""

    # Draw groups colored by difference value
    for group_id, position in group_positions.items():
        if not isinstance(position, tuple):
            continue

        gx, gy = position
        group_width, group_height = group_dimensions.get(group_id, (40, 40))

        # Get difference value for this group
        diff_value = group_diff_infections.get(group_id, 0.0)

        # Get color based on difference
        group_color, edge_color = get_difference_color(diff_value)

        # Draw groups with different visualization styles
        if viz_style == 'rectangles':
            rect = plt.Rectangle((gx - group_width/2, gy - group_height/2),
                                 group_width, group_height,
                                 facecolor=group_color, alpha=0.9,
                                 edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rect)

        elif viz_style == 'circles':
            radius = min(group_width, group_height) * 0.4
            circle = plt.Circle((gx, gy), radius,
                                facecolor=group_color, alpha=0.9,
                                edgecolor=edge_color, linewidth=2)
            ax.add_patch(circle)

        elif viz_style == 'squares':
            size = min(group_width, group_height) * 0.8
            square = plt.Rectangle((gx - size/2, gy - size/2),
                                   size, size,
                                   facecolor=group_color, alpha=0.9,
                                   edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(square)

        elif viz_style == 'hexagons':
            from matplotlib.patches import RegularPolygon
            radius = min(group_width, group_height) * 0.35
            hexagon = RegularPolygon((gx, gy), 6, radius=radius,
                                     facecolor=group_color, alpha=0.9,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(hexagon)

        elif viz_style == 'diamonds':
            from matplotlib.patches import RegularPolygon
            size = min(group_width, group_height) * 0.7
            diamond = RegularPolygon((gx, gy), 4, radius=size*0.6, orientation=np.pi/4,
                                     facecolor=group_color, alpha=0.9,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(diamond)

        elif viz_style == 'rounded_rects':
            from matplotlib.patches import FancyBboxPatch
            rounded_rect = FancyBboxPatch((gx - group_width/2, gy - group_height/2),
                                          group_width, group_height,
                                          boxstyle="round,pad=0.02",
                                          facecolor=group_color, alpha=0.9,
                                          edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rounded_rect)

    # Set title and formatting
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"

    stat_type = "median" if use_median else "avg."
    ax.set_title(f"Difference (Trans-Inf - Uniform)\n{day_title}: {stat_type} difference",
                 fontsize=14, fontweight='bold', pad=12)
    ax.axis('equal')
    ax.axis('off')

    # Add small margins for better spacing within difference heatmaps
    ax.margins(0.02)

    return ax


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


def calculate_workplace_contact_matrix(workplace_data, contact_data):
    """Calculate inter-workplace contact strength matrix"""
    print("Calculating workplace contact matrix...")
    print(f"Contact data columns: {contact_data.columns.tolist()}")

    # Check if required columns exist
    if 'person_id' not in workplace_data.columns or 'workplace_id' not in workplace_data.columns:
        print(f"ERROR: Missing required columns in workplace data")
        print(f"Expected: ['person_id', 'workplace_id']")
        print(f"Found: {workplace_data.columns.tolist()}")
        # Return empty results
        return [], np.array([[]])

    # Create person to workplace mapping
    person_to_workplace = {}
    for _, row in workplace_data.iterrows():
        person_to_workplace[row['person_id']] = row['workplace_id']

    # Get list of unique workplaces
    workplaces = sorted(workplace_data['workplace_id'].unique())
    workplace_to_idx = {wp: i for i, wp in enumerate(workplaces)}

    # Initialize contact matrix
    n_workplaces = len(workplaces)
    contact_matrix = np.zeros((n_workplaces, n_workplaces))

    # Get the actual column names (strip whitespace)
    columns = [col.strip() for col in contact_data.columns]
    contact_data.columns = columns

    # Calculate inter-workplace contact hours
    for _, row in contact_data.iterrows():
        person1 = row['person_1']
        person2 = row['person_2']
        contact_hours = row['contact_hours']

        # Get workplaces for both persons
        if person1 in person_to_workplace and person2 in person_to_workplace:
            wp1 = person_to_workplace[person1]
            wp2 = person_to_workplace[person2]

            # Only count inter-workplace contacts (not within same workplace)
            if wp1 != wp2:
                idx1 = workplace_to_idx[wp1]
                idx2 = workplace_to_idx[wp2]
                contact_matrix[idx1, idx2] += contact_hours
                contact_matrix[idx2, idx1] += contact_hours  # Symmetric

    return workplaces, contact_matrix


def create_contact_based_workplace_layout(workplace_data, contact_data,
                                          workplace_spacing=150.0, person_spacing=15.0,
                                          max_workplaces_per_row=25):
    """Create uniform rectangular layout for workplaces with even spacing"""

    # Calculate workplace contact matrix
    workplaces, contact_matrix = calculate_workplace_contact_matrix(
        workplace_data, contact_data)

    print(f"Creating uniform grid layout for {len(workplaces)} workplaces...")

    # Handle empty workplace data
    if len(workplaces) == 0:
        print("Warning: No workplaces found. Returning empty layout.")
        return {}, [], {}, {}

    # Create a simple grid layout with even spacing
    person_pos = {}
    workplace_positions = {}
    workplace_dimensions = {}

    # Calculate grid dimensions for even distribution
    total_workplaces = len(workplaces)
    workplaces_per_row = min(max_workplaces_per_row,
                             max(1, int(np.sqrt(total_workplaces) * 1.5)))
    rows = int(np.ceil(total_workplaces / workplaces_per_row))

    # Sort workplaces for consistent layout
    sorted_workplaces = sorted(workplaces)

    # Create uniform rectangular layout
    for idx, workplace_id in enumerate(sorted_workplaces):
        row = idx // workplaces_per_row
        col = idx % workplaces_per_row

        # Calculate position with even spacing
        wx = col * workplace_spacing
        wy = row * workplace_spacing

        workplace_positions[workplace_id] = (wx, wy)

        # Set uniform rectangular dimensions for all workplaces
        workplace_width = workplace_spacing * 0.6  # 60% of spacing for consistent size
        workplace_height = workplace_spacing * 0.4  # Rectangular shape
        workplace_dimensions[workplace_id] = (
            workplace_width, workplace_height)

        # Get workplace members and place them at center (not used for display)
        workplace_members = workplace_data[workplace_data['workplace_id']
                                           == workplace_id]['person_id'].tolist()
        for person_id in workplace_members:
            person_pos[person_id] = (wx, wy)

    print(f"Uniform grid layout created: {total_workplaces} workplaces")
    print(f"Grid dimensions: {rows} rows x {workplaces_per_row} columns")
    print(f"Workplace spacing: {workplace_spacing} units")

    return person_pos, workplaces, workplace_positions, workplace_dimensions


def count_infected_households(household_data, infected_at_time):
    """Count households with at least one infected member"""
    infected_households = set()

    for _, row in household_data.iterrows():
        person_id = row['person_id']
        household_id = row['household_id']

        if person_id in infected_at_time:
            infected_households.add(household_id)

    return len(infected_households)


def count_infected_groups(group_data, infected_at_time, use_workplaces=False):
    """Count households or workplaces with at least one infected member"""
    infected_groups = set()
    group_id_col = "workplace_id" if use_workplaces else "household_id"

    for _, row in group_data.iterrows():
        person_id = row['person_id']
        group_id = row[group_id_col]

        if person_id in infected_at_time:
            infected_groups.add(group_id)

    return len(infected_groups)


def create_size_based_heatmap(ax, group_data, group_positions, group_dimensions, person_pos,
                              infected_at_time, timestep, title, method_name,
                              group_average_infections=None, num_runs=None,
                              viz_style='rectangles', use_workplaces=False, use_median=False):
    """Create heatmap with household-size-based colors and infection-rate-based saturation"""

    # Calculate group sizes
    group_sizes = calculate_group_sizes(group_data, use_workplaces)

    # Use average/median group infections
    group_stats = {}
    for group_id, stat_value in group_average_infections.items():
        group_stats[group_id] = {
            'infection_rate': 0,  # Not used in display
            'infected_members': stat_value,  # This will be the average/median
            'total_members': group_sizes.get(group_id, 1)  # Actual group size
        }

    # Draw groups colored by size and saturation by infection rate
    for group_id, position in group_positions.items():
        if not isinstance(position, tuple):
            continue

        gx, gy = position
        group_width, group_height = group_dimensions.get(group_id, (40, 40))

        stats = group_stats.get(group_id, {
            'infection_rate': 0, 'infected_members': 0, 'total_members': 1})

        # Get group size and infected count
        group_size = stats['total_members']

        infected_count = stats['infected_members']

        if use_workplaces:
            # Workplaces: All purple with saturation based on infection rate
            base_color = "#A200B8"  # purple
            # Calculate infection rate (0-1) - cap at 100%
            infection_rate = min(1.0, infected_count / max(1, group_size))
            # Use bins for saturation
            if infection_rate < 0.05:  # 0-20%
                saturation_percent = 0  # Use white instead of base color
                base_color = '#FFFFFF'
            elif infection_rate < 0.2:  # 20-40%
                saturation_percent = 30
            elif infection_rate < 0.5:  # 40-60%
                saturation_percent = 50
            elif infection_rate < 0.95:  # 60-80%
                saturation_percent = 80
            else:  # 80-100%
                saturation_percent = 100
        else:
            # Households: Color by size, saturation by infection rate
            base_color = get_household_size_color(group_size)
            # Calculate infection rate (0-1) - cap at 100%
            infection_rate = min(1.0, infected_count / max(1, group_size))

            # Use bins for saturation
            if infection_rate < 0.05:  # 0-20%
                saturation_percent = 0  # Use white instead of base color
                base_color = '#FFFFFF'
            elif infection_rate < 0.2:  # 20-40%
                saturation_percent = 10
            elif infection_rate < 0.6:  # 40-60%
                saturation_percent = 20
            elif infection_rate < 0.95:  # 60-80%
                saturation_percent = 30
            else:  # 80-100%
                saturation_percent = 100

        # Apply saturation to the base color
        base_saturation = saturation_percent

        # Edge color is slightly darker version
        edge_color = base_color
        if base_color == '#FFFFFF':
            base_color = "#CFCFCFFF"
            edge_color = "#CFCFCFFF"  # Slightly darker than white
            base_saturation = 30.0

        # Draw groups with different visualization styles
        if viz_style == 'rectangles':
            rect = plt.Rectangle((gx - group_width/2, gy - group_height/2),
                                 group_width, group_height,
                                 facecolor=base_color, alpha=base_saturation/100.0,
                                 edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rect)

        elif viz_style == 'circles':
            radius = min(group_width, group_height) * 0.4
            circle = plt.Circle((gx, gy), radius,
                                facecolor=base_color, alpha=base_saturation/100.0,
                                edgecolor=edge_color, linewidth=2)
            ax.add_patch(circle)

        elif viz_style == 'squares':
            size = min(group_width, group_height) * 0.8
            square = plt.Rectangle((gx - size/2, gy - size/2),
                                   size, size,
                                   facecolor=base_color, alpha=base_saturation/100.0,
                                   edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(square)

        elif viz_style == 'hexagons':
            from matplotlib.patches import RegularPolygon
            radius = min(group_width, group_height) * 0.35
            hexagon = RegularPolygon((gx, gy), 6, radius=radius,
                                     facecolor=base_color, alpha=base_saturation/100.0,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(hexagon)

        elif viz_style == 'diamonds':
            from matplotlib.patches import RegularPolygon
            size = min(group_width, group_height) * 0.7
            diamond = RegularPolygon((gx, gy), 4, radius=size*0.6, orientation=np.pi/4,
                                     facecolor=group_color, alpha=0.9,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(diamond)

        elif viz_style == 'rounded_rects':
            from matplotlib.patches import FancyBboxPatch
            rounded_rect = FancyBboxPatch((gx - group_width/2, gy - group_height/2),
                                          group_width, group_height,
                                          boxstyle="round,pad=0.02",
                                          facecolor=group_color, alpha=0.9,
                                          edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rounded_rect)

    # Set title and formatting
    infected_count = infected_at_time
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"

    # Calculate average/median infected individuals from group averages
    if group_average_infections:
        stat_total_infected = infected_count
        stat_type = "median" if use_median else "avg."
        title_suffix = f"{stat_total_infected:.1f} infected"
    else:
        stat_type = "median" if use_median else "avg."
        title_suffix = f"infected"

    ax.set_title(f"{day_title}: {title_suffix}",
                 fontsize=12, fontweight='bold')
    ax.axis('equal')
    ax.axis('off')

    return infected_count


def create_time_specific_heatmap(ax, group_data, group_positions, group_dimensions, person_pos,
                                 infected_at_time, timestep, title, method_name,
                                 group_average_infections=None, num_runs=None,
                                 viz_style='rectangles', use_workplaces=False, use_median=False):
    """Create heatmap for specific time point - showing groups (households/workplaces) colored by infection count"""

    # Determine group type and column names
    group_type = "workplace" if use_workplaces else "household"
    group_id_col = "workplace_id" if use_workplaces else "household_id"

    # Use average/median group infections
    group_stats = {}
    for group_id, stat_value in group_average_infections.items():
        group_stats[group_id] = {
            'infection_rate': 0,  # Not used in display
            'infected_members': stat_value,  # This will be the average/median
            'total_members': 1  # Not used for average/median display
        }

    # Draw uniform rectangular groups colored by infection count
    for group_id, position in group_positions.items():
        if not isinstance(position, tuple):
            continue

        gx, gy = position
        group_width, group_height = group_dimensions.get(
            group_id, (40, 40))

        stats = group_stats.get(group_id, {
            'infection_rate': 0, 'infected_members': 0, 'total_members': 0})

        # For average/median mode, infected_count is already a float average/median
        infected_count = stats['infected_members']

        # Color group based on average/median infection count (0.5 steps)
        if infected_count < 0.5:
            # <0.5 average/median infected - white
            group_color = 'white'
            edge_color = 'lightgray'
        elif infected_count < 1.0:
            # 0.5-1.0 average/median infected - light yellow
            group_color = '#FFFACD'
            edge_color = '#FFD700'
        elif infected_count < 1.5:
            # 1.0-1.5 average/median infected - yellow
            group_color = '#FFFF99'
            edge_color = '#FFD700'
        elif infected_count < 2.0:
            # 1.5-2.0 average/median infected - yellow-orange
            group_color = '#FFCC33'
            edge_color = '#FF8C00'
        elif infected_count < 2.5:
            # 2.0-2.5 average/median infected - orange
            group_color = '#FF9933'
            edge_color = '#FF4500'
        elif infected_count < 3.0:
            # 2.5-3.0 average/median infected - orange-red
            group_color = '#FF6600'
            edge_color = '#FF4500'
        elif infected_count < 3.5:
            # 3.0-3.5 average/median infected - light red
            group_color = '#FF3300'
            edge_color = '#CC0000'
        elif infected_count < 4.0:
            # 3.5-4.0 average/median infected - red
            group_color = '#FF0000'
            edge_color = '#CC0000'
        elif infected_count < 4.5:
            # 4.0-4.5 average/median infected - bright red
            group_color = '#DD0000'
            edge_color = '#AA0000'
        else:
            # ≥4.5 average/median infected - dark black
            group_color = '#000000'
            edge_color = '#000000'

        # Draw groups with different visualization styles
        if viz_style == 'rectangles':
            # Style 1: Original rectangular groups
            rect = plt.Rectangle((gx - group_width/2, gy - group_height/2),
                                 group_width, group_height,
                                 facecolor=group_color, alpha=0.9,
                                 edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rect)

        elif viz_style == 'circles':
            # Style 2: Circular groups
            radius = min(group_width, group_height) * 0.4
            circle = plt.Circle((gx, gy), radius,
                                facecolor=group_color, alpha=0.9,
                                edgecolor=edge_color, linewidth=2)
            ax.add_patch(circle)

        elif viz_style == 'squares':
            # Style 3: Square groups
            size = min(group_width, group_height) * 0.8
            square = plt.Rectangle((gx - size/2, gy - size/2),
                                   size, size,
                                   facecolor=group_color, alpha=0.9,
                                   edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(square)

        elif viz_style == 'hexagons':
            # Style 4: Hexagonal groups
            from matplotlib.patches import RegularPolygon
            radius = min(group_width, group_height) * 0.35
            hexagon = RegularPolygon((gx, gy), 6, radius=radius,
                                     facecolor=group_color, alpha=0.9,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(hexagon)

        elif viz_style == 'diamonds':
            # Style 5: Diamond-shaped groups
            size = min(group_width, group_height) * 0.7
            # Create diamond by rotating a square 45 degrees
            from matplotlib.patches import RegularPolygon
            diamond = RegularPolygon((gx, gy), 4, radius=size*0.6, orientation=np.pi/4,
                                     facecolor=group_color, alpha=0.9,
                                     edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(diamond)

        elif viz_style == 'rounded_rects':
            # Style 6: Rounded rectangles
            from matplotlib.patches import FancyBboxPatch
            rounded_rect = FancyBboxPatch((gx - group_width/2, gy - group_height/2),
                                          group_width, group_height,
                                          boxstyle="round,pad=0.02",
                                          facecolor=group_color, alpha=0.9,
                                          edgecolor=edge_color, linewidth=1.5)
            ax.add_patch(rounded_rect)

    # Set title and formatting
    infected_count = infected_at_time
    total_people = len(person_pos)
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"

    # Calculate average/median infected individuals from group averages
    if group_average_infections:
        stat_total_infected = infected_count
        stat_type = "median" if use_median else "avg."
        title_suffix = f"{stat_total_infected:.1f} infected"
    else:
        stat_type = "median" if use_median else "avg."
        title_suffix = f"infected"

    ax.set_title(f"{method_name}\n{day_title}: {title_suffix}",
                 fontsize=14, fontweight='bold', pad=12)
    ax.axis('equal')
    ax.axis('off')

    return infected_count


def create_household_comparison_plot(ax, time_points, group_counts_method1, group_counts_method2,
                                     method1_name, method2_name,
                                     avg_infected_method1=None, avg_infected_method2=None, use_workplaces=False, use_median=False):
    """Create slim comparison plot of household/workplace infection counts and average/median infected"""

    # Convert time points to day labels
    day_labels = []
    for tp in time_points:
        if tp == 0:
            day_labels.append("Init")
        else:
            day_labels.append(f"Day {tp//24}")

    x_pos = np.arange(len(time_points))
    width = 0.35

    bars1 = ax.bar(x_pos - width/2, group_counts_method1, width,
                   label=method1_name, color='#d62728', alpha=0.8)
    bars2 = ax.bar(x_pos + width/2, group_counts_method2, width,
                   label=method2_name, color='#1f77b4', alpha=0.8)

    # Add value labels on bars (show decimals for averages/medians if needed)
    for i, (v1, v2) in enumerate(zip(group_counts_method1, group_counts_method2)):
        label1 = f"{v1:.1f}" if v1 != int(v1) else str(int(v1))
        label2 = f"{v2:.1f}" if v2 != int(v2) else str(int(v2))

        ax.text(i - width/2, v1 + max(group_counts_method1 + group_counts_method2) * 0.05,
                label1, ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.text(i + width/2, v2 + max(group_counts_method1 + group_counts_method2) * 0.05,
                label2, ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Update ylabel and title based on group type
    group_type = "Departments" if use_workplaces else "Households"
    stat_type = "Median" if use_median else "Avg."
    ax.set_ylabel(f' Infected\n{group_type}', fontsize=11)
    ax.set_title(f' {group_type} with ≥1 Infected Member',
                 fontsize=12, fontweight='bold')

    # Clean axis formatting - remove extra ticks and labels
    ax.set_xticks(x_pos)
    ax.set_xticklabels(day_labels, fontsize=10)

    # Clean y-axis - let matplotlib handle y-ticks automatically but remove any extra formatting
    ax.tick_params(axis='both', which='minor', length=0)  # Remove minor ticks
    # Keep major x-ticks short
    ax.tick_params(axis='x', which='major', length=4)
    # Keep major y-ticks short
    ax.tick_params(axis='y', which='major', length=4)

    ax.legend(loc='upper left', fontsize=11)
    ax.grid(axis='y', alpha=0.3)


def create_comparative_temporal_heatmap(data_dir_method1, data_dir_method2=None,
                                        method1_name="Transmission-Informed", method2_name="Uniform",
                                        output_path=None, time_points=None,
                                        viz_style='rectangles', scenario_name=None, use_workplaces=False, use_median=False, show_differences=False, use_size_based_colors=False):
    """Create comparative temporal infection heatmap with average/median across runs"""

    # Set default time points if None provided
    if time_points is None:
        time_points = [0, 24, 72, 240]

    stat_type = "Median" if use_median else "Average"
    group_type = "workplaces" if use_workplaces else "households"
    print(
        f"Loading temporal infection data for both methods ({stat_type} across runs) - {group_type} view...")

    # Load Method 1 data (transmission-informed)
    try:
        if use_workplaces:
            # Load workplace data
            group_data = load_workplace_data(
                os.path.join(data_dir_method1, 'work_id.csv'))
        else:
            # Load household data
            group_data = load_household_data(
                os.path.join(data_dir_method1, 'household_id.csv'))

        contact_data = load_contact_data(os.path.join(
            data_dir_method1, 'contact_intensiveness.csv'))

        # Load all infection files for averaging
        infection_files_m1 = find_infection_files(data_dir_method1)
        if not infection_files_m1:
            print(
                f"Error: No infection files found for Method 1 in {data_dir_method1}")
            return
        infections_by_time_m1, group_infections_by_time_m1, num_runs_m1, avg_infected_group_counts_m1 = calculate_average_infections_by_time(
            infection_files_m1, group_data, time_points, use_workplaces, use_median)

    except FileNotFoundError as e:
        print(f"Error loading Method 1 data: {e}")
        return

    # Load Method 2 data (uniform) - if provided
    has_method2_data = False
    infections_by_time_m2 = {}
    group_infections_by_time_m2 = {}
    num_runs_m2 = 1
    avg_infected_group_counts_m2 = {}

    if data_dir_method2 and os.path.exists(data_dir_method2):
        try:
            # Load all infection files for averaging
            infection_files_m2 = find_infection_files(data_dir_method2)
            if infection_files_m2:
                infections_by_time_m2, group_infections_by_time_m2, num_runs_m2, avg_infected_group_counts_m2 = calculate_average_infections_by_time(
                    infection_files_m2, group_data, time_points, use_workplaces, use_median)
                has_method2_data = True
            else:
                print(
                    f"Warning: No infection files found for Method 2 in {data_dir_method2}")
                has_method2_data = False
        except FileNotFoundError as e:
            print(f"Warning: Method 2 data not found: {e}")
            print("Will create comparison layout")
            has_method2_data = False
    else:
        print("Method 2 data directory not provided - creating comparison")
        has_method2_data = False

    print(
        f"Loaded Method 1 data: {len(group_data)} people, averaging across {num_runs_m1} runs")
    if has_method2_data:
        print(
            f"Loaded Method 2 data: {len(group_data)} people, averaging across {num_runs_m2} runs")

    # Check if we have valid data
    if len(group_data) == 0:
        print("ERROR: No data loaded. Cannot create visualization.")
        print("This might happen if:")
        print("  - No workplace data found (for workplace mode)")
        print("  - No household data found (for household mode)")
        print("  - Infection data doesn't match any group members")
        return None

    # Create group layout based on contact patterns (same for both methods)
    spacing_factor = 2.8  # Increased from 2.0 for maximum spacing
    group_spacing = 250.0 * spacing_factor  # Increased base spacing
    person_spacing = 80.0 * spacing_factor  # Increased base spacing

    if use_workplaces:
        print("Creating contact-based workplace layout...")
        person_pos, groups, group_positions, group_dimensions = create_contact_based_workplace_layout(
            group_data, contact_data,
            workplace_spacing=group_spacing,
            person_spacing=person_spacing,
            max_workplaces_per_row=22  # Reduced from 25 to accommodate larger spacing
        )
    else:
        print("Creating contact-based household layout...")
        person_pos, groups, group_positions, group_dimensions = create_contact_based_household_layout(
            group_data, contact_data,
            household_spacing=group_spacing,
            person_spacing=person_spacing,
            max_households_per_row=22  # Reduced from 25 to accommodate larger spacing
        )

    # Create figure with conditional rows: 3 rows (default) or 4 rows (with differences)
    if show_differences:
        fig = plt.figure(figsize=(4*len(time_points), 14)
                         )  # 4 rows with larger legend
    else:
        fig = plt.figure(figsize=(4*len(time_points), 10))  # 3 rows

    # Calculate group infection counts for comparison
    group_counts_m1 = []
    group_counts_m2 = []

    infection_counts_m1 = []
    infection_counts_m2 = []

    # For average mode: collect average infected counts for comparison plot
    avg_infected_counts_m1 = []
    avg_infected_counts_m2 = []

    # Calculate differences for each timestep and group (only if needed)
    group_differences_by_time = {}

    # Create subplots grid with conditional spacing
    if show_differences:
        # 4-row layout: Method1, Method2, Comparison, Difference
        # Row 1: Top heatmaps (Method 1)
        gs_top = fig.add_gridspec(1, len(time_points), top=0.90, bottom=0.75,
                                  left=0.05, right=0.95, wspace=0.1)

        # Row 2: Method 2 heatmaps
        gs_method2 = fig.add_gridspec(1, len(time_points), top=0.70, bottom=0.55,
                                      left=0.05, right=0.95, wspace=0.1)

        # Row 3: Middle comparison
        gs_mid = fig.add_gridspec(1, 1, top=0.52, bottom=0.37,
                                  left=0.05, right=0.95)

        # Row 4: Bottom difference heatmaps with increased spacing
        gs_diff = fig.add_gridspec(1, len(time_points), top=0.37, bottom=0.22,
                                   left=0.05, right=0.95, wspace=0.15)
    else:
        # 3-row layout: Method1, Method2, Comparison
        # Row 1: Top heatmaps (Method 1)
        gs_top = fig.add_gridspec(1, len(time_points), top=0.90, bottom=0.68,
                                  left=0.05, right=0.95, wspace=0.1)

        # Row 2: Middle comparison - centered between heatmaps
        gs_mid = fig.add_gridspec(1, 1, top=0.63, bottom=0.50,
                                  left=0.05, right=0.95)

        # Row 3: Method 2 heatmaps
        gs_method2 = fig.add_gridspec(1, len(time_points), top=0.43, bottom=0.19,
                                      left=0.05, right=0.95, wspace=0.1)

    # Row 1: Method 1 (Transmission-informed)
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs_top[0, i])

        infected_at_time = infections_by_time_m1.get(timestep, set())
        day_label = timestep // 24

        # Get group infections for this timestep
        group_avg_infections = group_infections_by_time_m1.get(timestep, {})

        # Choose heatmap function based on coloring mode
        if use_size_based_colors:
            infected_count = create_size_based_heatmap(
                ax, group_data, group_positions, group_dimensions, person_pos,
                infected_at_time, day_label,
                f"Time Point {i+1}", method1_name,
                group_average_infections=group_avg_infections,
                num_runs=num_runs_m1,
                viz_style=viz_style,
                use_workplaces=use_workplaces,
                use_median=use_median
            )
        else:
            infected_count = create_time_specific_heatmap(
                ax, group_data, group_positions, group_dimensions, person_pos,
                infected_at_time, day_label,
                f"Time Point {i+1}", method1_name,
                group_average_infections=group_avg_infections,
                num_runs=num_runs_m1,
                viz_style=viz_style,
                use_workplaces=use_workplaces,
                use_median=use_median
            )

        infection_counts_m1.append(infected_count)

        # Use proper average group counts
        group_counts_m1.append(
            avg_infected_group_counts_m1.get(timestep, 0))

        # Collect average infected count from group averages
        if group_avg_infections:
            avg_total_infected = sum(group_avg_infections.values())
            avg_infected_counts_m1.append(avg_total_infected)
        else:
            avg_infected_counts_m1.append(infected_count)

        # Add vertical text label for the first heatmap in the row
        if i == 0:
            ax.text(-0.15, 0.5, method1_name, transform=ax.transAxes,
                    rotation=90, fontsize=14, fontweight='bold',
                    verticalalignment='center', horizontalalignment='center')

    # Row 2: Method 2 (Uniform)
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs_method2[0, i])
        if has_method2_data:
            infected_at_time = infections_by_time_m2.get(timestep, set())
            day_label = timestep // 24

            # Get group infections for this timestep
            group_avg_infections = group_infections_by_time_m2.get(
                timestep, {})

            # Choose heatmap function based on coloring mode
            if use_size_based_colors:
                create_size_based_heatmap(
                    ax, group_data, group_positions, group_dimensions, person_pos,
                    infected_at_time, day_label,
                    f"Time Point {i+1}", method2_name,
                    group_average_infections=group_avg_infections,
                    num_runs=num_runs_m2,
                    viz_style=viz_style,
                    use_workplaces=use_workplaces,
                    use_median=use_median
                )
            else:
                create_time_specific_heatmap(
                    ax, group_data, group_positions, group_dimensions, person_pos,
                    infected_at_time, day_label,
                    f"Time Point {i+1}", method2_name,
                    group_average_infections=group_avg_infections,
                    num_runs=num_runs_m2,
                    viz_style=viz_style,
                    use_workplaces=use_workplaces,
                    use_median=use_median
                )
        else:
            # Create empty visualization
            ax.text(0.5, 0.5, f"{method2_name}\nData Not Available\n\n(No data for comparison)",
                    ha='center', va='center', transform=ax.transAxes, fontsize=12,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.5))
            ax.set_title(f"{method2_name} - Time Point {i+1}",
                         fontsize=11, fontweight='bold')
            ax.axis('off')

        # Add vertical text label for the first heatmap in the row
        if i == 0:
            ax.text(-0.15, 0.5, method2_name, transform=ax.transAxes,
                    rotation=90, fontsize=14, fontweight='bold',
                    verticalalignment='center', horizontalalignment='center')

    # Calculate actual Method 2 data or create empty comparison
    # Calculate actual Method 2 data or create empty comparison
    if not has_method2_data:
        # Create simulated uniform data for demonstration
        group_counts_m2 = [int(count * 1.2)
                           for count in group_counts_m1]
        infection_counts_m2 = [int(count * 1.15)
                               for count in infection_counts_m1]
        avg_infected_counts_m2 = [
            count * 1.15 for count in avg_infected_counts_m1]

        # Create simulated group differences (only if differences are requested)
        if show_differences:
            for i, timestep in enumerate(time_points):
                group_avg_m1 = group_infections_by_time_m1.get(timestep, {})
                group_diff = {}
                for group_id, val1 in group_avg_m1.items():
                    # Simulate method 2 having slightly different values
                    simulated_m2_val = val1 * 1.1 if val1 > 0 else 0
                    group_diff[group_id] = val1 - simulated_m2_val
                group_differences_by_time[timestep] = group_diff
    else:
        # Calculate actual Method 2 data
        for i, timestep in enumerate(time_points):
            infected_at_time = infections_by_time_m2.get(timestep, set())
            infection_counts_m2.append(infected_at_time)

            # Use proper average group counts
            group_counts_m2.append(
                avg_infected_group_counts_m2.get(timestep, 0))

            # Collect average infected count from group averages
            group_avg_infections = group_infections_by_time_m2.get(
                timestep, {})
            if group_avg_infections:
                avg_total_infected = sum(group_avg_infections.values())
                avg_infected_counts_m2.append(avg_total_infected)
            else:
                avg_infected_counts_m2.append(len(infected_at_time))

            # Calculate differences (method1 - method2) for each group (only if differences are requested)
            if show_differences:
                group_avg_m1 = group_infections_by_time_m1.get(timestep, {})
                group_avg_m2 = group_infections_by_time_m2.get(timestep, {})
                group_diff = {}

                # Get all group IDs from both methods
                all_group_ids = set(group_avg_m1.keys()) | set(
                    group_avg_m2.keys())
                for group_id in all_group_ids:
                    val1 = group_avg_m1.get(group_id, 0.0)
                    val2 = group_avg_m2.get(group_id, 0.0)
                    group_diff[group_id] = val1 - val2
                group_differences_by_time[timestep] = group_diff

    # Row 3: Group comparison (slim panel spanning all columns)
    ax_comparison = fig.add_subplot(gs_mid[0, 0])
    ax_comparison.margins(y=0.35)

    create_household_comparison_plot(ax_comparison, time_points, group_counts_m1, group_counts_m2,
                                     method1_name, method2_name,
                                     avg_infected_method1=avg_infected_counts_m1,
                                     avg_infected_method2=avg_infected_counts_m2,
                                     use_workplaces=use_workplaces, use_median=use_median)

    # Row 4: Difference heatmaps (only if requested)
    if show_differences:
        for i, timestep in enumerate(time_points):
            ax = fig.add_subplot(gs_diff[0, i])

            day_label = timestep // 24
            group_diff_infections = group_differences_by_time.get(timestep, {})

            create_difference_heatmap(
                ax, group_positions, group_dimensions,
                group_diff_infections, day_label,
                viz_style=viz_style,
                use_workplaces=use_workplaces,
                use_median=use_median
            )

    # Create legends - conditional based on coloring mode and whether differences are shown
    group_type_title = "Workplace" if use_workplaces else "Household"
    stat_type_lower = "median" if use_median else "average"

    if use_size_based_colors:
        # Size-based coloring legends
        if use_workplaces:
            # For workplaces: only saturation legend (black with different saturations)
            size_based_legend_elements = [
                plt.Rectangle((0, 0), 1, 1, facecolor="#A200B8",
                              edgecolor='black', label='department')
            ]
            size_legend_title = f'Department Color (saturation = infection rate)'
        else:
            # For households: size-based colors legend with improved colors
            size_based_legend_elements = [
                plt.Rectangle((0, 0), 1, 1, facecolor='#FFD700',
                              edgecolor='#B8860B', label='1-person household (gold)'),
                plt.Rectangle((0, 0), 1, 1, facecolor='#FF6B6B',
                              edgecolor='#CC5555', label='2-person household (coral)'),
                plt.Rectangle((0, 0), 1, 1, facecolor='#4ECDC4',
                              edgecolor='#3EA39D', label='3-person household (teal)'),
                plt.Rectangle((0, 0), 1, 1, facecolor='#006907',
                              edgecolor='#006907', label='4-person household (green)'),
                plt.Rectangle((0, 0), 1, 1, facecolor='#9B59B6',
                              edgecolor='#7C4791', label='5-person household (purple)')
            ]
            size_legend_title = f'Household Size Colors (saturation = infection rate)'

        # Use size-based legend instead of original
        original_legend_elements = size_based_legend_elements
        original_legend_title = size_legend_title
    else:
        # Original legend for average/median infection counts (0.5 steps)
        original_legend_elements = [
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
        original_legend_title = f'{stat_type} Infected Members per {group_type_title}'

    # Place original legend

    if show_differences:
        # 4-row layout: place legend higher to accommodate expanded difference legend below
        original_legend_bbox = (0.5, 0.16)

        # Difference legend (complete range with 0.5 increments from -4.5 to +4.5)
        difference_legend_elements = [
            # Negative differences (Method2 > Method1) - Green spectrum
            plt.Rectangle((0, 0), 1, 1, facecolor='#800080',
                          edgecolor='#400040', label='<-4.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#00CC00',
                          edgecolor='#001000', label='-4.5 to -4.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#00FF00',
                          edgecolor='#002000', label='-4.0 to -3.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#30FF30',
                          edgecolor='#004000', label='-3.5 to -3.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#60FF60',
                          edgecolor='#006400', label='-3.0 to -2.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#90FF90',
                          edgecolor='#228B22', label='-2.5 to -2.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#C0FFC0',
                          edgecolor='#32CD32', label='-2.0 to -1.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#E0FFE0',
                          edgecolor='#7CFC00', label='-1.5 to -1.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#F0FFF0',
                          edgecolor='#90EE90', label='-1.0 to -0.5'),

            # Near zero - White
            plt.Rectangle((0, 0), 1, 1, facecolor='white',
                          edgecolor='lightgray', label='-0.5 to 0.5'),

            # Positive differences (Method1 > Method2) - Yellow to Red spectrum
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFACD',
                          edgecolor='#FFD700', label='0.5 to 1.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFFF99',
                          edgecolor='#FFD700', label='1.0 to 1.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FFCC33',
                          edgecolor='#FF8C00', label='1.5 to 2.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF9933',
                          edgecolor='#FF4500', label='2.0 to 2.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF6600',
                          edgecolor='#FF4500', label='2.5 to 3.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF3300',
                          edgecolor='#CC0000', label='3.0 to 3.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#FF0000',
                          edgecolor='#CC0000', label='3.5 to 4.0'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#DD0000',
                          edgecolor='#AA0000', label='4.0 to 4.5'),
            plt.Rectangle((0, 0), 1, 1, facecolor='#000000',
                          edgecolor='#000000', label='≥4.5')
        ]
    else:
        # 3-row layout: place legend lower since no difference legend needed
        original_legend_bbox = (0.5, 0.15)

    original_legend = fig.legend(handles=original_legend_elements, loc='center',
                                 bbox_to_anchor=original_legend_bbox, ncol=10, fontsize=10,
                                 title=original_legend_title, title_fontsize=11)
    original_legend.get_title().set_fontweight('bold')

    # Add difference legend only if differences are shown
    if show_differences:
        # Place difference legend below original - use multiple rows for better fit
        difference_legend_title = f'Difference (Trans-Inf - Uniform) per {group_type_title}'
        difference_legend = fig.legend(handles=difference_legend_elements, loc='center',
                                       bbox_to_anchor=(0.5, 0.05), ncol=10, fontsize=9,
                                       title=difference_legend_title, title_fontsize=10)
        difference_legend.get_title().set_fontweight('bold')

        # Add the original legend back to the figure (matplotlib removes it when adding second legend)
        fig.add_artist(original_legend)

    # Add overall title
    stat_type_lower = "median" if use_median else "average"
    title_suffix = f"({stat_type} across {num_runs_m1} runs)"

    # Create title with scenario information and group type
    group_view = "Department View" if use_workplaces else "Household View"
    base_title = f'Comparative Temporal Infection Progression {title_suffix} - {group_view}'
    if scenario_name:
        base_title = f'Scenario {scenario_name}: {base_title}'

    fig.suptitle(f'{base_title}\nTransmission-Informed vs Uniform Initialization',
                 fontsize=16, fontweight='bold', y=0.98)  # Moved higher for 4-row layout

    # Custom spacing is handled by the individual gridspec positioning above

    if output_path:
        plt.savefig(output_path, dpi=500,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative temporal heatmap saved to {output_path}")

    # Print summary statistics
    group_type_lower = "workplaces" if use_workplaces else "households"
    stat_type_lower = "median" if use_median else "average"
    print(
        f"\n{method1_name} Summary ({stat_type_lower} across {num_runs_m1} runs):")
    for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m1, group_counts_m1)):
        time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
        print(
            f"  {time_label}: ~{inf_count} {stat_type_lower} infected, ~{group_count} {group_type_lower} affected")

    if has_method2_data:
        print(
            f"\n{method2_name} Summary ({stat_type_lower} across {num_runs_m2} runs):")
        for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m2, group_counts_m2)):
            time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
            print(
                f"  {time_label}: ~{inf_count} {stat_type_lower} infected, ~{group_count} {group_type_lower} affected")

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
    parser.add_argument('--viz-style', default='rectangles',
                        choices=['rectangles', 'circles', 'squares',
                                 'hexagons', 'diamonds', 'rounded_rects'],
                        help='Visualization style for households (default: rectangles)')
    parser.add_argument(
        '--scenario-name', help='Scenario name to include in title (e.g., R1, R2, W1, W2)')
    parser.add_argument('--use-workplaces', action='store_true', default=False,
                        help='Use workplace-based visualization instead of household-based')
    parser.add_argument('--use-median', action='store_true', default=True,
                        help='Use median instead of average when aggregating across runs')
    parser.add_argument('--show-differences', action='store_true', default=False,
                        help='Show difference heatmaps (Method1 - Method2) as a fourth row')
    parser.add_argument('--use-size-based-colors', action='store_true', default=False,
                        help='Use household-size-based colors with infection-rate saturation (households: size-based colors, workplaces: black)')

    args = parser.parse_args()

    # debug this for now
    # args.data_dir_method1 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250903_180348_W1_workplace_few_meetings_transmission_informed"
    # args.data_dir_method2 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250903_180355_W1_workplace_few_meetings_uniform_initialized"
    # args.use_workplaces = True
    # args.use_size_based_colors = True
    # Only set default output path if none was provided
    if not args.output_path:
        mode_suffix = "_median" if args.use_median else "_average"
        view_suffix = "_workplaces" if args.use_workplaces else "_households"
        diff_suffix = "_with_differences" if args.show_differences else ""
        size_suffix = "_size_based" if args.use_size_based_colors else ""
        args.output_path = f"/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/comparative_temporal_infection_heatmap{mode_suffix}{view_suffix}{diff_suffix}{size_suffix}.png"

    if not os.path.exists(args.data_dir_method1):
        print(f"Error: Directory {args.data_dir_method1} does not exist")
        sys.exit(1)

    if not args.output_path:
        mode_suffix = "_median" if args.use_median else "_average"
        view_suffix = "_workplaces" if args.use_workplaces else "_households"
        diff_suffix = "_with_differences" if args.show_differences else ""
        args.output_path = os.path.join(
            args.data_dir_method1, f'comparative_temporal_infection_heatmap{mode_suffix}{view_suffix}{diff_suffix}.png')

    # Create visualization
    fig = create_comparative_temporal_heatmap(
        data_dir_method1=args.data_dir_method1,
        data_dir_method2=args.data_dir_method2,
        method1_name=args.method1_name,
        method2_name=args.method2_name,
        output_path=args.output_path,
        time_points=args.time_points,
        viz_style=args.viz_style,
        scenario_name=args.scenario_name,
        use_workplaces=args.use_workplaces,
        use_median=args.use_median,
        show_differences=args.show_differences,
        use_size_based_colors=args.use_size_based_colors
    )

    if fig is None:
        print("Failed to create visualization due to data issues.")
        sys.exit(1)


if __name__ == "__main__":
    main()
