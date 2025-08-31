#!/usr/bin/env python3
"""
Comparative Infection Tree Visualization
Creates side-by-side transmission tree comparisons for different initialization strategies
Shows transmission chain differences and generation patterns with opportunity analysis

Usage:
    python comparative_infection_trees.py --data-dir-method1 <transmission_informed_dir> --data-dir-method2 <uniform_dir>
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
import numpy as np
from infection_timeline import (
    load_simulation_data, transform_infection_data, create_contact_data_for_analysis,
    build_transmission_tree, calculate_tree_positions, has_contact_at_location
)


def analyze_transmission_opportunities(infection_events, contact_data, transmission_tree):
    """Analyze transmission opportunities and calculate missed/realized transmissions"""

    opportunities = {
        'realized_transmissions': 0,
        'missed_opportunities': 0,
        'total_contacts': 0,
        'transmission_efficiency': 0.0
    }

    # Count realized transmissions
    for infector, infected_list in transmission_tree.items():
        opportunities['realized_transmissions'] += len(infected_list)

    # Estimate missed opportunities by looking at contacts that didn't lead to transmission
    infectious_people = set()
    for event in infection_events:
        infectious_people.add(event['person_id'])

    # Count total potential transmission contacts by looking at locations where infectious people were present
    # Group by timestep and location to find people at same place/time
    potential_transmissions = 0
    for (timestep, location), group in contact_data.groupby(['Timestep', 'Location_ID']):
        people = group['Person_ID'].tolist()

        # Check if any infectious people were at this location at this time
        infectious_present = any(
            person_id in infectious_people for person_id in people)

        if infectious_present and len(people) > 1:
            # Count potential transmission opportunities (n choose 2 for n people at same location/time)
            n = len(people)
            potential_transmissions += n * (n - 1) // 2  # All possible pairs

    opportunities['total_contacts'] = potential_transmissions
    opportunities['missed_opportunities'] = max(
        0, potential_transmissions - opportunities['realized_transmissions'])

    if potential_transmissions > 0:
        opportunities['transmission_efficiency'] = opportunities['realized_transmissions'] / \
            potential_transmissions

    return opportunities


def create_comparative_infection_trees(data_dir_method1, data_dir_method2,
                                       method1_name="Transmission-Informed",
                                       method2_name="Uniform",
                                       max_infections=50, output_path=None,
                                       scenario_name=None):
    """Create side-by-side infection tree comparison"""

    print(
        f"Creating comparative infection trees for {scenario_name or 'scenario'}...")

    # Load data for both methods
    contact_file_1 = os.path.join(
        data_dir_method1, 'best_run_contact_data.csv')
    infection_file_1 = os.path.join(
        data_dir_method1, 'best_run_detailed_infection.csv')
    infected_status_file_1 = os.path.join(
        data_dir_method1, 'best_run_infected_status.csv')

    contact_file_2 = os.path.join(
        data_dir_method2, 'best_run_contact_data.csv')
    infection_file_2 = os.path.join(
        data_dir_method2, 'best_run_detailed_infection.csv')
    infected_status_file_2 = os.path.join(
        data_dir_method2, 'best_run_infected_status.csv')

    # Check if files exist
    for file_path, name in [(contact_file_1, f"{method1_name} contact"),
                            (infection_file_1, f"{method1_name} infection"),
                            (infected_status_file_1,
                             f"{method1_name} infected status"),
                            (contact_file_2, f"{method2_name} contact"),
                            (infection_file_2, f"{method2_name} infection"),
                            (infected_status_file_2, f"{method2_name} infected status")]:
        if not os.path.exists(file_path):
            print(f"Error: {name} file not found: {file_path}")
            return

    # Load and process data for method 1
    print(f"Loading {method1_name} data...")
    contact_df_1, infection_df_1 = load_simulation_data(
        contact_file_1, infection_file_1)
    infection_events_1 = transform_infection_data(infection_df_1)
    contact_analysis_df_1 = create_contact_data_for_analysis(
        contact_df_1, infection_df_1)

    # Load infected status data
    infected_status_df_1 = pd.read_csv(infected_status_file_1)
    print(f"Loaded infected status data: {len(infected_status_df_1)} records")

    # Load and process data for method 2
    print(f"Loading {method2_name} data...")
    contact_df_2, infection_df_2 = load_simulation_data(
        contact_file_2, infection_file_2)
    infection_events_2 = transform_infection_data(infection_df_2)
    contact_analysis_df_2 = create_contact_data_for_analysis(
        contact_df_2, infection_df_2)

    # Load infected status data
    infected_status_df_2 = pd.read_csv(infected_status_file_2)
    print(f"Loaded infected status data: {len(infected_status_df_2)} records")
    contact_analysis_df_2 = create_contact_data_for_analysis(
        contact_df_2, infection_df_2)

    # Limit infections for visualization clarity
    if len(infection_events_1) > max_infections:
        infection_events_1 = infection_events_1[:max_infections]
    if len(infection_events_2) > max_infections:
        infection_events_2 = infection_events_2[:max_infections]

    # Update infection events to rename EventPanvadere based on scenario
    def update_event_location_names(events, scenario_name):
        updated_events = []
        for event in events:
            updated_event = event.copy()
            if event['location_type'] == 'EventPanvadere':
                if scenario_name and ('W' in scenario_name or 'work' in scenario_name.lower()):
                    updated_event['location_type'] = 'Workplace Event'
                else:
                    updated_event['location_type'] = 'Restaurant Event'
            updated_events.append(updated_event)
        return updated_events

    infection_events_1 = update_event_location_names(
        infection_events_1, scenario_name)
    infection_events_2 = update_event_location_names(
        infection_events_2, scenario_name)

    # Build transmission trees using contact data AND infected status data
    print("Building transmission trees...")
    tree_1, generations_1, potential_1 = build_transmission_tree(
        infection_events_1, contact_analysis_df_1, infected_status_df_1)
    tree_2, generations_2, potential_2 = build_transmission_tree(
        infection_events_2, contact_analysis_df_2, infected_status_df_2)

    # Calculate positions
    positions_1 = calculate_tree_positions(infection_events_1, generations_1)
    positions_2 = calculate_tree_positions(infection_events_2, generations_2)

    # Create visualization with better proportions
    fig = plt.figure(figsize=(20, 18))

    # Create a 3-row layout: trees on top, opportunity analysis below (more space for bottom)
    gs = fig.add_gridspec(3, 2, height_ratios=[
                          2.5, 2.5, 1.5], hspace=0.3, wspace=0.2)

    # Location colors - will be dynamically updated based on scenario
    base_location_colors = {
        'Work': '#4169E1',        # Royal Blue
        'Home': '#2E8B57',        # Sea Green
        'School': '#FF6347',      # Tomato Red
        'SocialEvent': '#9370DB',  # Medium Purple
        'BasicsShop': '#FF8C00',  # Dark Orange
        'Patient Zero': "#FF0000",  # Red
        # Default Deep Pink - will be overridden based on scenario
        'EventPanvadere': '#FF1493'
    }

    # Update EventPanvadere color based on scenario type
    location_colors = base_location_colors.copy()
    if scenario_name and ('W' in scenario_name or 'work' in scenario_name.lower()):
        # Blue Violet for Workplace
        location_colors['EventPanvadere'] = '#8A2BE2'
        location_colors['Workplace Event Cluster'] = '#8A2BE2'
    else:
        # Deep Pink for Restaurant
        location_colors['EventPanvadere'] = '#FF1493'
        location_colors['Restaurant Event Cluster'] = '#FF1493'

    # Plot Method 1 tree
    ax1 = fig.add_subplot(gs[0:2, 0])
    plot_infection_tree(ax1, infection_events_1, tree_1, positions_1, location_colors,
                        f"{method1_name}\n{len(infection_events_1)} infections", generations_1, infected_status_df_1, contact_analysis_df_1)

    # Plot Method 2 tree
    ax2 = fig.add_subplot(gs[0:2, 1])
    plot_infection_tree(ax2, infection_events_2, tree_2, positions_2, location_colors,
                        f"{method2_name}\n{len(infection_events_2)} infections", generations_2, infected_status_df_2, contact_analysis_df_2)

    # Create transmission opportunities over time analysis (bottom panel spanning both columns)
    ax3 = fig.add_subplot(gs[2, :])
    plot_transmission_opportunities_over_time(
        ax3, infection_events_1, infection_events_2, method1_name, method2_name,
        contact_analysis_df_1, contact_analysis_df_2, infected_status_df_1, infected_status_df_2)

    # Add overall title
    title = "Comparative Infection Transmission Trees"
    if scenario_name:
        title = f"Scenario {scenario_name}: {title}"
    fig.suptitle(f"{title}\n{method1_name} vs {method2_name} Initialization",
                 fontsize=16, fontweight='bold', y=0.95)

    # Create legend (use the first subplot for positioning)
    legend_elements = []
    for loc_type, color in location_colors.items():
        # Skip EventPanvadere, Patient Zero, and specific cluster entries from legend
        if loc_type not in ['EventPanvadere', 'Workplace Event Cluster', 'Restaurant Event Cluster', 'Patient Zero']:
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=color, markersize=10,
                                              label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='black', lw=2, linestyle='-',
                   label='Certain transmission', alpha=0.9),
        plt.Line2D([0], [0], color='black', lw=2, linestyle='--',
                   label='Uncertain transmission ', alpha=0.9)
    ])

    ax1.legend(handles=legend_elements, loc='upper left', title='Legend',
               title_fontsize=10, bbox_to_anchor=(0, 1), fontsize=14)

    plt.tight_layout()

    if output_path:
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=500,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative infection trees saved to {output_path}")

    # Print comparison summary with simplified metrics
    print(f"\n=== COMPARATIVE TRANSMISSION ANALYSIS ===")
    print(f"{method1_name}:")
    print(f"  Total infections: {len(infection_events_1)}")

    # Calculate clustering metrics
    clustering_1 = calculate_clustering_metrics(infection_events_1, tree_1)
    print(
        f"  Same-location transmissions: {clustering_1['same_location_transmissions']}")
    print(
        f"  Clustering coefficient: {clustering_1['clustering_coefficient']:.3f}")

    print(f"\n{method2_name}:")
    print(f"  Total infections: {len(infection_events_2)}")

    clustering_2 = calculate_clustering_metrics(infection_events_2, tree_2)
    print(
        f"  Same-location transmissions: {clustering_2['same_location_transmissions']}")
    print(
        f"  Clustering coefficient: {clustering_2['clustering_coefficient']:.3f}")

    # Add difference analysis
    infection_diff = len(infection_events_2) - len(infection_events_1)
    clustering_diff = clustering_2['clustering_coefficient'] - \
        clustering_1['clustering_coefficient']
    print(f"\nDIFFERENCES ({method2_name} - {method1_name}):")
    print(f"  Total infections: {infection_diff:+d}")
    print(f"  Clustering coefficient: {clustering_diff:+.3f}")

    return fig


def calculate_clustering_metrics(infection_events, transmission_tree):
    """Calculate clustering-related metrics for transmission patterns"""

    metrics = {
        'same_location_transmissions': 0,
        'total_transmissions': 0,
        'clustering_coefficient': 0.0,
        'location_clusters': {}
    }

    # Count same-location transmissions
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_location = event_dict.get(
            infector_id, {}).get('location_type', 'Unknown')

        for infected_id in infected_list:
            metrics['total_transmissions'] += 1
            infected_location = event_dict.get(
                infected_id, {}).get('location_type', 'Unknown')

            if infector_location == infected_location and infector_location != 'Unknown':
                metrics['same_location_transmissions'] += 1

            # Count location clusters
            if infected_location not in metrics['location_clusters']:
                metrics['location_clusters'][infected_location] = 0
            metrics['location_clusters'][infected_location] += 1

    # Calculate clustering coefficient
    if metrics['total_transmissions'] > 0:
        metrics['clustering_coefficient'] = metrics['same_location_transmissions'] / \
            metrics['total_transmissions']

    return metrics


def plot_infection_tree(ax, infection_events, transmission_tree, positions, location_colors, title, generations=None, infected_status_df=None, contact_data=None):
    """Plot a single infection tree on the given axes with enhanced clustering visualization"""

    times = [event['time'] for event in infection_events]

    # Set up y-axis values early for use throughout the function
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values) if y_values else 1

    # Group infections by location type to show clustering
    location_groups = {}
    for event in infection_events:
        loc_type = event['location_type']
        if loc_type not in location_groups:
            location_groups[loc_type] = []
        location_groups[loc_type].append(event)

    # # Draw background rectangles to highlight location clusters
    # for loc_type, events in location_groups.items():
    #     if len(events) > 5 and loc_type != 'Patient Zero':  # Only show clusters with >5 infections
    #         y_positions = [positions[e['person_id']] for e in events]
    #         times_loc = [e['time'] for e in events]

    #         # Create a subtle background rectangle for each cluster
    #         # Ensure cluster uses the same color as EventPanvadere when it's the location type
    #         if loc_type == 'EventPanvadere':
    #             rect_color = location_colors.get('EventPanvadere', '#FF1493')
    #         else:
    #             rect_color = location_colors.get(
    #                 loc_type, '#FF1493')  # Use same default as points
    #         ax.axhspan(min(y_positions) - 0.8, max(y_positions) + 0.8,
    #                    xmin=(min(times_loc) - min(times)) /
    #                    (max(times) - min(times)) - 0.02,
    #                    xmax=(max(times_loc) - min(times)) /
    #                    (max(times) - min(times)) + 0.02,
    #                    alpha=0.15, color=rect_color, zorder=0)

    #         # Add cluster label for significant clusters - positioned to the right
    #         ax.text(max(times_loc) + (max(times) - min(times)) * 0.05, np.mean(y_positions),
    #                 f"{loc_type} Cluster\n({len(events)} cases)",
    #                 ha='left', va='center', fontsize=12,
    #                 bbox=dict(boxstyle='round,pad=0.3',
    #                           facecolor=rect_color, alpha=0.3),
    #                 fontweight='bold')

    # Draw transmission tree connections with enhanced styling
    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        # Find infector event and position
        infector_event = next(
            e for e in infection_events if e['person_id'] == infector_id)
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for infected_id in infected_list:
            infected_event = next(
                e for e in infection_events if e['person_id'] == infected_id)
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Use thicker lines for same-location transmissions
            same_location = (
                infector_event['location_type'] == infected_event['location_type'])

            # Check if transmission could be from multiple sources (uncertain)
            # Count how many people had contact AND were infectious at the time of infection
            potential_infectors = 0
            if infected_status_df is not None and contact_data is not None:
                # Find patient zero time for special case handling
                patient_zero_time = min(e['time'] for e in infection_events)

                for e in infection_events:
                    if e['time'] < infected_event['time']:  # Infected before current person

                        # Check if they had contact at the location
                        if has_contact_at_location(e['person_id'], infected_event['person_id'],
                                                   infected_event['location'], infected_event['time'], contact_data):

                            # Special case: assume patient zero is always infectious, even before timestep 0
                            is_patient_zero = (e['time'] == patient_zero_time)

                            if is_patient_zero:
                                # Patient zero is assumed infectious from before simulation start
                                potential_infectors += 1
                            else:
                                # For timesteps < 0, only patient zero can be infectious
                                if infected_event['time'] < 0:
                                    continue  # Only patient zero infectious before timestep 0

                                # Check if this person was actually infectious at time of transmission
                                infectious_check = infected_status_df[
                                    (infected_status_df['Person_ID'] == e['person_id']) &
                                    (infected_status_df['Timestep'] == int(infected_event['time'])) &
                                    (infected_status_df['Infected'] == 1)
                                ]

                                if len(infectious_check) > 0:
                                    potential_infectors += 1
            else:
                # Fallback to location-type based logic if no data available
                potential_infectors = sum(1 for e in infection_events
                                          if e['time'] < infected_event['time']
                                          and e['location_type'] == infected_event['location_type'])

            line_width = 2  # Use consistent line width for all transmissions
            line_alpha = 0.9 if same_location else 0.6
            line_color = 'black'  # Always use black for transmission lines

            # Use dashed line for uncertain transmissions (multiple potential infectors)
            line_style = '--' if potential_infectors > 1 else '-'

            # Draw vertical first, then horizontal connection
            ax.plot([infector_time, infector_time], [infector_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)
            ax.plot([infector_time, infected_time], [infected_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)

            # Add enhanced arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time - (infected_time -
                                infector_time) * 0.1, infected_y),
                        arrowprops=dict(arrowstyle='->', color=line_color, lw=line_width, alpha=line_alpha))

    # Plot infection points with simplified labeling
    for event in infection_events:
        y_pos = positions[event['person_id']]
        # Better color lookup with specific handling for EventPanvadere
        location_type = event['location_type']
        if location_type == 'EventPanvadere':
            # Ensure EventPanvadere gets the right color
            color = location_colors.get(
                'EventPanvadere', '#FF1493')  # Default to Deep Pink
        else:
            # Use Deep Pink instead of gray as default
            color = location_colors.get(location_type, '#FF1493')
        size = 300 if event['time'] == min(times) else 200

        # Use consistent border width
        border_width = 1

        ax.scatter(event['time'], y_pos, s=size, c=color, alpha=1.0,
                   zorder=5, edgecolors='black', linewidth=border_width)

    # Calculate y limits first
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values) if y_values else 1

    # Add timeline markers for key events - all labels at same height
    outbreak_time = min(times)  # Patient zero time
    label_y_pos = -0.5  # Fixed position slightly below y=0

    # Add simulation start line (timestep 0)
    ax.axvline(x=0, color='blue', linestyle='-',
               alpha=0.8, linewidth=3, zorder=1)
    ax.text(0, label_y_pos, 'Simulation Start',
            ha='center', va='top', fontsize=14, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='blue', alpha=0.7, edgecolor='darkblue'))

    if len(times) > 1:
        # Find a reasonable "outbreak event" time (e.g., when multiple infections start)
        sorted_times = sorted(times)
        outbreak_event_time = sorted_times[min(
            1, len(sorted_times)-1)]  # 4th infection or last

        # Add vertical line for patient zero
        ax.axvline(x=outbreak_time, color='red', linestyle='--',
                   alpha=0.7, linewidth=3, zorder=1)
        ax.text(outbreak_time, label_y_pos, 'Patient Zero',
                ha='center', va='top', fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.7, edgecolor='darkred'))

        # Always plot outbreak event line and label, even if close to simulation start
        outbreak_event = next(
            (e for e in infection_events if e['time'] == outbreak_event_time), None)
        if outbreak_event:
            loc_type = outbreak_event.get('location_type', 'EventPanvadere')
            outbreak_color = location_colors.get(loc_type, '#FF1493')
        else:
            outbreak_color = '#FF1493'

        ax.axvline(x=outbreak_event_time, color=outbreak_color,
                   linestyle='--', alpha=0.7, linewidth=3, zorder=1)
        ax.text(outbreak_event_time, label_y_pos, 'Outbreak Event',
                ha='center', va='top', fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=outbreak_color, alpha=0.7, edgecolor=outbreak_color))

    # Formatting with bigger fonts
    ax.set_xlabel('Timestep', fontsize=18)
    ax.set_ylabel('Transmission Order', fontsize=18)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Increase tick label sizes drastically
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=14)

    # Set limits with padding
    ax.set_ylim(min(y_values) - y_range * 0.05, max(y_values) + y_range * 0.1)

    # Statistics box removed per user request


def plot_transmission_opportunities_over_time(ax, infection_events_1, infection_events_2, method1_name, method2_name,
                                              contact_data_1, contact_data_2, infected_status_df_1, infected_status_df_2):
    """
    Plot transmission opportunities over time based on actual contact data.

    Algorithm:
    1. At each timepoint, find all infectious agents
    2. Look up where these infectious agents are located
    3. Count how many other (non-infectious) agents are at those same locations
    4. Sum gives the total transmission opportunities at that timepoint
    """

    # Get the full time range from both methods
    times_1 = [event['time'] for event in infection_events_1]
    times_2 = [event['time'] for event in infection_events_2]
    all_times = sorted(set(times_1 + times_2))

    if not all_times:
        return

    min_time = min(all_times)
    max_time = max(all_times)

    # Create time range for full simulation duration (240 timesteps)
    time_range = np.arange(0, 240, 1)

    def calculate_opportunities(time_point, infection_events, contact_data, infected_status_df):
        """Calculate transmission opportunities at a specific timepoint"""

        # Step 1: Find all infectious agents at this timepoint
        infectious_agents = set()

        # Get agents who were infected before or at this timepoint
        for event in infection_events:
            if event['time'] <= time_point:
                person_id = event['person_id']

                # Check if this person is actually infectious at this timepoint
                infectious_check = infected_status_df[
                    (infected_status_df['Person_ID'] == person_id) &
                    (infected_status_df['Timestep'] == int(time_point)) &
                    (infected_status_df['Infected'] == 1)
                ]

                if len(infectious_check) > 0:
                    infectious_agents.add(person_id)

        if not infectious_agents:
            return 0

        # Step 2: Find locations where infectious agents are present at this timepoint
        infectious_locations = set()
        infectious_agents_at_locations = {}

        current_contacts = contact_data[contact_data['Timestep'] == int(
            time_point)]

        for _, contact in current_contacts.iterrows():
            person_id = contact['Person_ID']
            location_id = contact['Location_ID']

            if person_id in infectious_agents:
                infectious_locations.add(location_id)
                if location_id not in infectious_agents_at_locations:
                    infectious_agents_at_locations[location_id] = set()
                infectious_agents_at_locations[location_id].add(person_id)

        # Step 3: Count non-infectious agents at these locations
        transmission_opportunities = 0

        for location_id in infectious_locations:
            # Get all people at this location at this timepoint
            people_at_location = current_contacts[
                current_contacts['Location_ID'] == location_id
            ]['Person_ID'].tolist()

            # Count non-infectious people (susceptible to infection)
            susceptible_at_location = [
                person_id for person_id in people_at_location
                if person_id not in infectious_agents
            ]
            if (len(susceptible_at_location) > 50):
                print(susceptible_at_location)

            transmission_opportunities += len(susceptible_at_location)

        return transmission_opportunities

    # Calculate opportunities for both methods
    opportunities_1 = []
    opportunities_2 = []

    for t in time_range:
        opp_1 = calculate_opportunities(
            t, infection_events_1, contact_data_1, infected_status_df_1)
        opp_2 = calculate_opportunities(
            t, infection_events_2, contact_data_2, infected_status_df_2)

        opportunities_1.append(opp_1)
        opportunities_2.append(opp_2)

    # Plot the curves
    ax.plot(time_range, opportunities_1, 'o-', color='#d62728', linewidth=3,
            markersize=4, alpha=0.8, label=f'{method1_name}')
    ax.plot(time_range, opportunities_2, 's-', color='#1f77b4', linewidth=3,
            markersize=4, alpha=0.8, label=f'{method2_name}')

    # Fill area between curves to show difference
    ax.fill_between(time_range, opportunities_1, opportunities_2,
                    alpha=0.3, color='lightgray', label='Difference')

    # Calculate cumulative opportunities for day differences
    timesteps_per_day = 24  # 1 timestep = 1 hour
    day_1_idx = min(timesteps_per_day, len(
        opportunities_1) - 1)  # 1 day = 24 timesteps
    # 3 days = 72 timesteps
    day_3_idx = min(3 * timesteps_per_day, len(opportunities_1) - 1)
    # 5 days = 120 timesteps
    day_5_idx = min(5 * timesteps_per_day, len(opportunities_1) - 1)
    # 7 days = 168 timesteps
    day_7_idx = min(7 * timesteps_per_day, len(opportunities_1) - 1)

    cum_opp_1 = np.cumsum(opportunities_1)
    cum_opp_2 = np.cumsum(opportunities_2)
    cum_opp_rel = (cum_opp_2 - cum_opp_1)
    cum_opp_perc = (cum_opp_rel / cum_opp_1 * 100)

    cum_diff_1_day = cum_opp_2[day_1_idx] - \
        cum_opp_1[day_1_idx] if day_1_idx < len(cum_opp_1) else 0
    cum_diff_3_days = cum_opp_2[day_3_idx] - \
        cum_opp_1[day_3_idx] if day_3_idx < len(cum_opp_1) else 0
    cum_diff_5_days = cum_opp_2[day_5_idx] - \
        cum_opp_1[day_5_idx] if day_5_idx < len(cum_opp_1) else 0
    cum_diff_7_days = cum_opp_2[day_7_idx] - \
        cum_opp_1[day_7_idx] if day_7_idx < len(cum_opp_1) else 0

    # Add day boundary markers with differences
    day_boundaries = [timesteps_per_day, 3 * timesteps_per_day,
                      5 * timesteps_per_day, 7 * timesteps_per_day]
    day_numbers = [1, 3, 5, 7]
    day_diffs = [cum_diff_1_day, cum_diff_3_days,
                 cum_diff_5_days, cum_diff_7_days]
    colors = ['green', 'orange', 'purple', 'blue']

    # All day labels at same height
    max_opp = max(max(opportunities_1), max(opportunities_2)
                  ) if opportunities_1 and opportunities_2 else 100
    y_positions = [max_opp * 0.9] * 4

    for i, (boundary, day_num, diff, color) in enumerate(zip(day_boundaries, day_numbers, day_diffs, colors)):
        if boundary <= max(time_range):
            # Add vertical line with thicker width
            ax.axvline(x=boundary, color=color,
                       linestyle=':', alpha=0.7, linewidth=3)

            # Add compact annotation with bigger font
            if opportunities_1 and opportunities_2:
                diff_text = f"Day {day_num}: {diff:+.0f}"
                ax.text(boundary - 10, y_positions[i], diff_text, rotation=0, fontsize=12,
                        bbox=dict(boxstyle='round,pad=0.3',
                                  facecolor=color, alpha=0.4),
                        ha='left', va='center', fontweight='bold')

    # Formatting with bigger fonts
    ax.set_xlabel(
        'Timestep', fontsize=18)
    ax.set_ylabel(
        'Susceptible People at Risk', fontsize=18)
    ax.set_title(
        (
            f'Transmission Opportunities Over Time\n'
            f'Day X: Cumulative Difference in At-Risk Population\n'
            f'Final Cumulative Difference: {cum_opp_rel[-1]:.0f} ({cum_opp_perc[-1]:.1f}%)'
        ),
        fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=18)

    # Increase tick label sizes drastically
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=14)

    # Set x-axis limits to start from 0
    ax.set_xlim(0, max(time_range) + 1)


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description='Create comparative infection tree visualization')
    parser.add_argument('--data-dir-method1',
                        help='Directory containing transmission-informed data files')
    parser.add_argument('--data-dir-method2',
                        help='Directory containing uniform initialization data files')
    parser.add_argument('--method1-name', default='Transmission-Informed',
                        help='Name for first method')
    parser.add_argument('--method2-name', default='Uniform',
                        help='Name for second method')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--scenario-name',
                        help='Scenario name for titles (e.g., R1, W2)')
    parser.add_argument('--max-infections', type=int, default=20,
                        help='Maximum infections to display for clarity')

    args = parser.parse_args()
    # for debug set filepaths directly
    # args.data_dir_method1 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250831_222421_R1_restaurant_strong_clustering_transmission_informed"
    # args.data_dir_method2 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/epidemic_curves_20250831_222428_R1_restaurant_strong_clustering_uniform_initialized"
    # args.method1_name = "Method 1"
    # args.method2_name = "Method 2"
    # args.output_path = "/path/to/output.png"
    # args.scenario_name = "Scenario 1"
    # args.max_infections = 100

    # Check if directories exist
    for dir_path, name in [(args.data_dir_method1, "Method 1"),
                           (args.data_dir_method2, "Method 2")]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
            sys.exit(1)

    # Set default output path if none provided
    if not args.output_path:
        base_name = f"comparative_infection_trees"
        if args.scenario_name:
            base_name += f"_{args.scenario_name}"
        args.output_path = f"{base_name}.png"

    # Create visualization
    fig = create_comparative_infection_trees(
        args.data_dir_method1, args.data_dir_method2,
        args.method1_name, args.method2_name,
        args.max_infections, args.output_path, args.scenario_name
    )

    if fig is None:
        print("Failed to create visualization due to data issues.")
        sys.exit(1)

    print("Comparative infection tree analysis complete!")


if __name__ == "__main__":
    main()
