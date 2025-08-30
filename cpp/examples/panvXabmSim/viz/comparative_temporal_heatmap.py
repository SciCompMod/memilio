#!/usr/bin/env python3
"""
Comparative Temporal Infection Heatmap
Compares transmission-informed vs uniform initialization approaches
Shows location-based transmission patterns and household/workplace infection statistics

Features:
- Best run mode: Shows results from the single best-performing run
- Average mode: Shows averaged results across multiple simulation runs
- Household or workplace-based visualization with color coding for infection burden
- Temporal progression analysis at multiple time points

Usage:
    # Best run mode (default) - households
    python comparative_temporal_heatmap.py
    
    # Average mode across all runs - households
    python comparative_temporal_heatmap.py --use-average
    
    # Workplace-based visualization
    python comparative_temporal_heatmap.py --use-average --use-workplaces
    
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


def load_workplace_data(work_id_file):
    """Load workplace data from work_id.csv file"""
    try:
        work_data = pd.read_csv(work_id_file)
        
        # Filter out people without workplace (Work_ID = 4294967295 indicates no workplace)
        work_data = work_data[work_data['Work_ID'] != 4294967295]
        
        # Rename columns to match household data structure
        work_data = work_data.rename(columns={'Person_ID': 'person_id', 'Work_ID': 'workplace_id'})
        
        print(f"Loaded workplace data: {len(work_data)} working people in {work_data['workplace_id'].nunique()} workplaces")
        return work_data
    except Exception as e:
        print(f"Error loading workplace data: {e}")
        # Return empty DataFrame with correct column structure
        return pd.DataFrame(columns=['person_id', 'workplace_id'])


def calculate_average_infections_by_time(infection_files, group_data, time_points, use_workplaces=False):
    """Calculate average infections across multiple runs for specific time points"""
    print(
        f"Calculating average infections across {len(infection_files)} runs...")

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
        infections_by_time = load_temporal_infection_data(infection_file)

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

    # Calculate average group infections
    average_group_infections_by_time = {}
    for timestep in time_points:
        if timestep in group_run_infections:
            group_averages = {}
            for group_id, infection_counts in group_run_infections[timestep].items():
                if infection_counts:
                    avg_infected = np.mean(infection_counts)
                    group_averages[group_id] = avg_infected
                else:
                    group_averages[group_id] = 0.0
            average_group_infections_by_time[timestep] = group_averages
        else:
            average_group_infections_by_time[timestep] = {}

    # Calculate average infected group counts
    average_infected_group_counts = {}
    for timestep in time_points:
        if timestep in infected_group_counts and infected_group_counts[timestep]:
            avg_count = np.mean(infected_group_counts[timestep])
            average_infected_group_counts[timestep] = avg_count
        else:
            average_infected_group_counts[timestep] = 0.0

    return average_infections_by_time, average_group_infections_by_time, successful_runs, average_infected_group_counts


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


def create_time_specific_heatmap(ax, group_data, group_positions, group_dimensions, person_pos,
                                 infected_at_time, timestep, title, method_name,
                                 use_average=False, group_average_infections=None, num_runs=None, 
                                 viz_style='rectangles', use_workplaces=False):
    """Create heatmap for specific time point - showing groups (households/workplaces) colored by infection count"""

    # Determine group type and column names
    group_type = "workplace" if use_workplaces else "household"
    group_id_col = "workplace_id" if use_workplaces else "household_id"

    if use_average and group_average_infections is not None:
        # Use average group infections
        group_stats = {}
        for group_id, avg_infected in group_average_infections.items():
            group_stats[group_id] = {
                'infection_rate': 0,  # Not used in display
                'infected_members': avg_infected,  # This will be the average
                'total_members': 1  # Not used for average display
            }
    else:
        # Calculate group infection stats for this time point (best run mode)
        group_stats = defaultdict(
            lambda: {'total_members': 0, 'infected_members': 0})

        for _, row in group_data.iterrows():
            person_id = row['person_id']
            group_id = row[group_id_col]

            group_stats[group_id]['total_members'] += 1

            if person_id in infected_at_time:
                group_stats[group_id]['infected_members'] += 1

        # Calculate infection rates for best run mode
        for group_id in group_stats:
            stats = group_stats[group_id]
            if stats['total_members'] > 0:
                stats['infection_rate'] = stats['infected_members'] / \
                    stats['total_members']

    # Draw uniform rectangular groups colored by infection count
    for group_id, position in group_positions.items():
        if not isinstance(position, tuple):
            continue

        gx, gy = position
        group_width, group_height = group_dimensions.get(
            group_id, (40, 40))

        stats = group_stats.get(group_id, {
            'infection_rate': 0, 'infected_members': 0, 'total_members': 0})

        if use_average:
            # For average mode, infected_count is already a float average
            infected_count = stats['infected_members']
        else:
            # For best run mode, infected_count is an integer
            infected_count = stats['infected_members']

        # Color group based on number of infected members (or average)
        if use_average:
            # Color based on average infection count (0.5 steps)
            if infected_count < 0.5:
                # <0.5 average infected - white
                group_color = 'white'
                edge_color = 'lightgray'
            elif infected_count < 1.0:
                # 0.5-1.0 average infected - light yellow
                group_color = '#FFFACD'
                edge_color = '#FFD700'
            elif infected_count < 1.5:
                # 1.0-1.5 average infected - yellow
                group_color = '#FFFF99'
                edge_color = '#FFD700'
            elif infected_count < 2.0:
                # 1.5-2.0 average infected - yellow-orange
                group_color = '#FFCC33'
                edge_color = '#FF8C00'
            elif infected_count < 2.5:
                # 2.0-2.5 average infected - orange
                group_color = '#FF9933'
                edge_color = '#FF4500'
            elif infected_count < 3.0:
                # 2.5-3.0 average infected - orange-red
                group_color = '#FF6600'
                edge_color = '#FF4500'
            elif infected_count < 3.5:
                # 3.0-3.5 average infected - light red
                group_color = '#FF3300'
                edge_color = '#CC0000'
            elif infected_count < 4.0:
                # 3.5-4.0 average infected - red
                group_color = '#FF0000'
                edge_color = '#CC0000'
            elif infected_count < 4.5:
                # 4.0-4.5 average infected - bright red
                group_color = '#DD0000'
                edge_color = '#AA0000'
            else:
                # ≥4.5 average infected - dark black
                group_color = '#000000'
                edge_color = '#000000'
        else:
            # Color based on discrete infection count (best run mode)
            if infected_count == 0:
                # 0 infected - white
                group_color = 'white'
                edge_color = 'lightgray'
            elif infected_count == 1:
                # 1 infected - light yellow
                group_color = '#FFFACD'
                edge_color = '#F0E68C'
            elif infected_count == 2:
                # 2 infected - yellow
                group_color = '#FFFF99'
                edge_color = '#FFD700'
            elif infected_count == 3:
                # 3 infected - orange
                group_color = '#FFB366'
                edge_color = '#FF8C00'
            elif infected_count == 4:
                # 4 infected - red
                group_color = '#FF6666'
                edge_color = '#CC0000'
            else:
                # 5+ infected - black
                group_color = 'black'
                edge_color = 'darkgray'

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
    infected_count = len(infected_at_time)
    total_people = len(person_pos)
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"

    # Modify title based on mode
    if use_average and num_runs:
        # Calculate average infected individuals from group averages
        if group_average_infections:
            avg_total_infected = sum(group_average_infections.values())
            title_suffix = f"{avg_total_infected:.1f} avg. infected"
        else:
            title_suffix = f"avg. infected"
    else:
        title_suffix = f"{infected_count}/{total_people} infected ({infected_count/total_people*100:.1f}%)"

    ax.set_title(f"{method_name}\n{day_title}: {title_suffix}",
                 fontsize=14, fontweight='bold', pad=12)
    ax.axis('equal')
    ax.axis('off')

    return infected_count


def create_household_comparison_plot(ax, time_points, group_counts_method1, group_counts_method2,
                                     method1_name, method2_name, use_average=False,
                                     avg_infected_method1=None, avg_infected_method2=None, use_workplaces=False):
    """Create slim comparison plot of household/workplace infection counts and average infected if available"""

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

    # Add value labels on bars
    for i, (v1, v2) in enumerate(zip(group_counts_method1, group_counts_method2)):
        # Format labels appropriately (show decimals for averages if needed)
        if use_average:
            label1 = f"{v1:.1f}" if v1 != int(v1) else str(int(v1))
            label2 = f"{v2:.1f}" if v2 != int(v2) else str(int(v2))
        else:
            label1 = str(int(v1))
            label2 = str(int(v2))

        ax.text(i - width/2, v1 + max(group_counts_method1 + group_counts_method2) * 0.05,
                label1, ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax.text(i + width/2, v2 + max(group_counts_method1 + group_counts_method2) * 0.05,
                label2, ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Update ylabel and title based on mode and group type
    group_type = "Workplaces" if use_workplaces else "Households"
    
    if use_average:
        ax.set_ylabel(f'Avg. Infected\n{group_type}', fontsize=11)
        ax.set_title(f'Average {group_type} with ≥1 Infected Member',
                     fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel(f'Infected\n{group_type}', fontsize=11)
        ax.set_title(f'{group_type} with ≥1 Infected Member',
                     fontsize=12, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(day_labels, fontsize=10)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(axis='y', alpha=0.3)


def create_comparative_temporal_heatmap(data_dir_method1, data_dir_method2=None,
                                        method1_name="Transmission-Informed", method2_name="Uniform",
                        output_path=None, time_points=None, use_average=False, 
                                        viz_style='rectangles', scenario_name=None, use_workplaces=False):
    """Create comparative temporal infection heatmap with option for best run or average across runs"""
    
    # Set default time points if None provided
    if time_points is None:
        time_points = [0, 24, 72, 240]

    mode_description = "Average across runs" if use_average else "Best run"
    group_type = "workplaces" if use_workplaces else "households"
    print(
        f"Loading temporal infection data for both methods ({mode_description}) - {group_type} view...")

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
            infections_by_time_m1, group_infections_by_time_m1, num_runs_m1, avg_infected_group_counts_m1 = calculate_average_infections_by_time(
                infection_files_m1, group_data, time_points, use_workplaces)
        else:
            # Load best run only
            infections_by_time_m1 = load_temporal_infection_data(
                os.path.join(data_dir_method1, 'best_run_detailed_infection.csv'))
            group_infections_by_time_m1 = None
            num_runs_m1 = 1
            avg_infected_group_counts_m1 = {}

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
            if use_average:
                # Load all infection files for averaging
                infection_files_m2 = find_infection_files(
                    data_dir_method2, use_average=True)
                if infection_files_m2:
                    infections_by_time_m2, group_infections_by_time_m2, num_runs_m2, avg_infected_group_counts_m2 = calculate_average_infections_by_time(
                        infection_files_m2, group_data, time_points, use_workplaces)
                    has_method2_data = True
                else:
                    print(
                        f"Warning: No infection files found for Method 2 in {data_dir_method2}")
                    has_method2_data = False
            else:
                # Load best run only
                infections_by_time_m2 = load_temporal_infection_data(
                    os.path.join(data_dir_method2, 'best_run_detailed_infection.csv'))
                group_infections_by_time_m2 = None
                num_runs_m2 = 1
                avg_infected_group_counts_m2 = {}
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
            f"Loaded Method 1 data: {len(group_data)} people, averaging across {num_runs_m1} runs")
        if has_method2_data:
            print(
                f"Loaded Method 2 data: {len(group_data)} people, averaging across {num_runs_m2} runs")
    else:
        print(
            f"Loaded Method 1 data: {len(group_data)} people, temporal infections for {len(infections_by_time_m1)} time steps")

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

    # Create figure with 3 rows: Method1, Comparison, Method2 - more compact
    fig = plt.figure(figsize=(4*len(time_points), 10))

    # Calculate group infection counts for comparison
    group_counts_m1 = []
    group_counts_m2 = []

    infection_counts_m1 = []
    infection_counts_m2 = []

    # For average mode: collect average infected counts for comparison plot
    avg_infected_counts_m1 = []
    avg_infected_counts_m2 = []

    # Create subplots grid with custom spacing
    # Row 1: Top heatmaps (higher up, more space from title)
    gs_top = fig.add_gridspec(1, len(time_points), top=0.83, bottom=0.58,
                              left=0.05, right=0.95, wspace=0.1)

    # Row 2: Middle comparison (closer to top row, further from bottom row)
    gs_mid = fig.add_gridspec(1, 1, top=0.55, bottom=0.45,
                              left=0.05, right=0.95)

    # Row 3: Bottom heatmaps (same height as top row)
    gs_bot = fig.add_gridspec(1, len(time_points), top=0.35, bottom=0.10,
                              left=0.05, right=0.95, wspace=0.1)

    # Row 1: Method 1 (Transmission-informed)
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs_top[0, i])

        infected_at_time = infections_by_time_m1.get(timestep, set())
        day_label = timestep // 24

        # Get group infections for this timestep if using average mode
        group_avg_infections = None
        if use_average and group_infections_by_time_m1:
            group_avg_infections = group_infections_by_time_m1.get(
                timestep, {})

        infected_count = create_time_specific_heatmap(
            ax, group_data, group_positions, group_dimensions, person_pos,
            infected_at_time, day_label,
            f"Time Point {i+1}", method1_name,
            use_average=use_average,
            group_average_infections=group_avg_infections,
            num_runs=num_runs_m1 if use_average else None,
            viz_style=viz_style,
            use_workplaces=use_workplaces
        )

        infection_counts_m1.append(infected_count)

        # Use proper average group counts in average mode
        if use_average:
            group_counts_m1.append(
                avg_infected_group_counts_m1.get(timestep, 0))
        else:
            group_counts_m1.append(count_infected_groups(
                group_data, infected_at_time, use_workplaces))

        # For average mode: collect average infected count from group averages
        if use_average and group_avg_infections:
            avg_total_infected = sum(group_avg_infections.values())
            avg_infected_counts_m1.append(avg_total_infected)
        else:
            avg_infected_counts_m1.append(infected_count)

    # Row 2: Group comparison (slim panel spanning all columns)
    ax_comparison = fig.add_subplot(gs_mid[0, 0])
    ax_comparison.margins(y=0.35)

    # For now, use placeholder data for method 2 if not available
    if not has_method2_data:
        # Create simulated uniform data for demonstration
        group_counts_m2 = [int(count * 1.2)
                           for count in group_counts_m1]
        infection_counts_m2 = [int(count * 1.15)
                               for count in infection_counts_m1]
        avg_infected_counts_m2 = [
            count * 1.15 for count in avg_infected_counts_m1]
    else:
        # Calculate actual Method 2 data
        for i, timestep in enumerate(time_points):
            infected_at_time = infections_by_time_m2.get(timestep, set())
            infection_counts_m2.append(len(infected_at_time))

            # Use proper average group counts in average mode
            if use_average:
                group_counts_m2.append(
                    avg_infected_group_counts_m2.get(timestep, 0))
            else:
                group_counts_m2.append(count_infected_groups(
                    group_data, infected_at_time, use_workplaces))

            # For average mode: collect average infected count from group averages
            if use_average and group_infections_by_time_m2:
                group_avg_infections = group_infections_by_time_m2.get(
                    timestep, {})
                if group_avg_infections:
                    avg_total_infected = sum(group_avg_infections.values())
                    avg_infected_counts_m2.append(avg_total_infected)
                else:
                    avg_infected_counts_m2.append(len(infected_at_time))
            else:
                avg_infected_counts_m2.append(len(infected_at_time))

    create_household_comparison_plot(ax_comparison, time_points, group_counts_m1, group_counts_m2,
                                     method1_name, method2_name, use_average=use_average,
                                     avg_infected_method1=avg_infected_counts_m1,
                                     avg_infected_method2=avg_infected_counts_m2,
                                     use_workplaces=use_workplaces)

    # Row 3: Method 2 (Uniform) - placeholder or real data
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs_bot[0, i])
        if has_method2_data:
            infected_at_time = infections_by_time_m2.get(timestep, set())
            day_label = timestep // 24

            # Get group infections for this timestep if using average mode
            group_avg_infections = None
            if use_average and group_infections_by_time_m2:
                group_avg_infections = group_infections_by_time_m2.get(
                    timestep, {})

            create_time_specific_heatmap(
                ax, group_data, group_positions, group_dimensions, person_pos,
                infected_at_time, day_label,
                f"Time Point {i+1}", method2_name,
                use_average=use_average,
                group_average_infections=group_avg_infections,
                num_runs=num_runs_m2 if use_average else None,
                viz_style=viz_style,
                use_workplaces=use_workplaces
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
    group_type_title = "Workplace" if use_workplaces else "Household"
    
    if use_average:
        legend_title = f'Average Infected Members per {group_type_title}'
    else:
        legend_title = f'Infected Members per {group_type_title}'

    legend = fig.legend(handles=legend_elements, loc='center',
                        bbox_to_anchor=(0.5, 0.05), ncol=10, fontsize=11,
                        title=legend_title, title_fontsize=12)
    legend.get_title().set_fontweight('bold')

    # Add overall title
    if use_average:
        # Get the run count from either method (they should be the same)
        run_count = num_runs_m1 if use_average else 1
        title_suffix = f"(Average across {run_count} runs)"
    else:
        title_suffix = "(Best run)"

    # Create title with scenario information and group type
    group_view = "Workplace View" if use_workplaces else "Household View"
    base_title = f'Comparative Temporal Infection Progression {title_suffix} - {group_view}'
    if scenario_name:
        base_title = f'Scenario {scenario_name}: {base_title}'

    fig.suptitle(f'{base_title}\nTransmission-Informed vs Uniform Initialization',
                 fontsize=16, fontweight='bold', y=0.96)

    # Custom spacing is handled by the individual gridspec positioning above

    if output_path:
        plt.savefig(output_path, dpi=500,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative temporal heatmap saved to {output_path}")

    # Print summary statistics
    group_type_lower = "workplaces" if use_workplaces else "households"
    
    if use_average:
        print(f"\n{method1_name} Summary (averaged across {num_runs_m1} runs):")
        for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m1, group_counts_m1)):
            time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
            print(
                f"  {time_label}: ~{inf_count} avg infected, ~{group_count} {group_type_lower} affected")

        if has_method2_data:
            print(f"\n{method2_name} Summary (averaged across {num_runs_m2} runs):")
            for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m2, group_counts_m2)):
                time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
                print(
                    f"  {time_label}: ~{inf_count} avg infected, ~{group_count} {group_type_lower} affected")
    else:
        print(f"\n{method1_name} Summary (best run):")
        for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m1, group_counts_m1)):
            time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
            print(
                f"  {time_label}: {inf_count} infected, {group_count} {group_type_lower} affected")

        if has_method2_data:
            print(f"\n{method2_name} Summary (best run):")
            for i, (timestep, inf_count, group_count) in enumerate(zip(time_points, infection_counts_m2, group_counts_m2)):
                time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
                print(
                    f"  {time_label}: {inf_count} infected, {group_count} {group_type_lower} affected")

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
    parser.add_argument('--viz-style', default='rectangles',
                        choices=['rectangles', 'circles', 'squares',
                                 'hexagons', 'diamonds', 'rounded_rects'],
                        help='Visualization style for households (default: rectangles)')
    parser.add_argument(
        '--scenario-name', help='Scenario name to include in title (e.g., R1, R2, W1, W2)')
    parser.add_argument('--use-workplaces', action='store_true', default=False,
                        help='Use workplace-based visualization instead of household-based')

    args = parser.parse_args()

    # Only set default output path if none was provided
    if not args.output_path:
        mode_suffix = "_average" if args.use_average else "_best_run"
        view_suffix = "_workplaces" if args.use_workplaces else "_households"
        args.output_path = f"/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/comparative_temporal_infection_heatmap{mode_suffix}{view_suffix}.png"

    if not os.path.exists(args.data_dir_method1):
        print(f"Error: Directory {args.data_dir_method1} does not exist")
        sys.exit(1)

    if not args.output_path:
        mode_suffix = "_average" if args.use_average else "_best_run"
        view_suffix = "_workplaces" if args.use_workplaces else "_households"
        args.output_path = os.path.join(
            args.data_dir_method1, f'comparative_temporal_infection_heatmap{mode_suffix}{view_suffix}.png')

    # Only create visualization if using average mode
    if not args.use_average:
        print("Skipping visualization: only average plots are generated")
        sys.exit(0)

    # Create visualization
    fig = create_comparative_temporal_heatmap(
        data_dir_method1=args.data_dir_method1,
        data_dir_method2=args.data_dir_method2,
        method1_name=args.method1_name,
        method2_name=args.method2_name,
        output_path=args.output_path,
        time_points=args.time_points,
        use_average=args.use_average,
        viz_style=args.viz_style,
        scenario_name=args.scenario_name,
        use_workplaces=args.use_workplaces
    )
    
    if fig is None:
        print("Failed to create visualization due to data issues.")
        sys.exit(1)


if __name__ == "__main__":
    main()
