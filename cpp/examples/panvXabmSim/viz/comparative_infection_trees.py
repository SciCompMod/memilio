#!/usr/bin/env python3
"""
Enhanced Comparative Infection Tree Visualization
Creates side-by-side transmission tree comparisons for different initialization strategies
with multiple layout and connection style options

Usage:
    python enhanced_comparative_infection_trees.py --data-dir-method1 <transmission_informed_dir> --data-dir-method2 <uniform_dir>
                                                   --layout-method <current|hierarchical|location|temporal>
                                                   --connection-style <current|L-shaped|curved|straight>
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys
import numpy as np
from matplotlib.patches import FancyArrowPatch
from infection_timeline import (
    load_simulation_data, transform_infection_data, create_contact_data_for_analysis,
    build_transmission_tree, has_contact_at_location
)


# Layout calculation functions
def calculate_current_positions(infection_events, transmission_tree, infection_generations):
    """Original tree-based layout from the existing code"""
    generations = {}
    for person_id, generation in infection_generations.items():
        if generation not in generations:
            generations[generation] = []
        generations[generation].append(person_id)

    person_positions = {}
    current_y = 0

    for gen in sorted(generations.keys()):
        people_in_gen = generations[gen]
        people_with_times = []
        for person_id in people_in_gen:
            event = next(
                e for e in infection_events if e['person_id'] == person_id)
            people_with_times.append((person_id, event['time']))

        people_with_times.sort(key=lambda x: x[1])

        for i, (person_id, _) in enumerate(people_with_times):
            person_positions[person_id] = current_y + i

        current_y += len(people_in_gen) + 2

    return person_positions


def calculate_transmission_order_positions(infection_events, transmission_tree, infection_generations):
    """Calculate positions based on transmission order to naturally avoid crossings"""
    # Sort all events by time to get chronological order
    sorted_events = sorted(infection_events, key=lambda x: x['time'])

    positions = {}
    current_y = 0

    # Process events in chronological order
    for event in sorted_events:
        person_id = event['person_id']
        positions[person_id] = current_y
        current_y += 1

    return positions


def calculate_depth_first_positions(infection_events, transmission_tree, infection_generations):
    """Calculate positions using depth-first traversal to minimize crossings"""
    positions = {}
    current_y = 0
    visited = set()

    # Find patient zero
    patient_zero = min(infection_events, key=lambda x: x['time'])
    patient_zero_id = patient_zero['person_id']

    def dfs_assign_positions(person_id):
        nonlocal current_y
        if person_id in visited:
            return

        visited.add(person_id)
        positions[person_id] = current_y
        current_y += 1

        # Get children and sort by infection time
        if person_id in transmission_tree:
            children = transmission_tree[person_id]
            if children:
                # Sort children by their infection time
                event_dict = {e['person_id']: e for e in infection_events}
                children_with_times = [
                    (child, event_dict[child]['time']) for child in children]
                children_with_times.sort(key=lambda x: x[1])

                # Process children in order
                for child_id, _ in children_with_times:
                    dfs_assign_positions(child_id)

    # Start DFS from patient zero
    dfs_assign_positions(patient_zero_id)

    # Handle any unvisited nodes (shouldn't happen in a proper tree)
    for event in infection_events:
        if event['person_id'] not in visited:
            dfs_assign_positions(event['person_id'])

    return positions


def calculate_current_positions_no_crossings(infection_events, transmission_tree, infection_generations):
    """Original tree-based layout but with y-positions reordered to minimize crossings"""
    generations = {}
    for person_id, generation in infection_generations.items():
        if generation not in generations:
            generations[generation] = []
        generations[generation].append(person_id)

    person_positions = {}
    current_y = 0

    for gen in sorted(generations.keys()):
        people_in_gen = generations[gen]
        people_with_times = []
        for person_id in people_in_gen:
            event = next(
                e for e in infection_events if e['person_id'] == person_id)
            people_with_times.append((person_id, event['time']))

        # For each generation, we need to consider transmission relationships
        # to minimize crossings when drawing vertical-then-horizontal lines
        if gen == 0:
            # Patient zero - just place at current position
            people_with_times.sort(key=lambda x: x[1])  # Sort by time
            for i, (person_id, _) in enumerate(people_with_times):
                person_positions[person_id] = current_y + i
        else:
            # For other generations, group by parent and sort to minimize crossings
            parent_groups = {}
            for person_id, time in people_with_times:
                # Find parent in transmission tree
                parent = None
                for potential_parent, children in transmission_tree.items():
                    if person_id in children:
                        parent = potential_parent
                        break

                if parent not in parent_groups:
                    parent_groups[parent] = []
                parent_groups[parent].append((person_id, time))

            # Sort parents by their y-position to maintain order
            sorted_parents = []
            for parent in parent_groups.keys():
                if parent in person_positions:
                    sorted_parents.append((parent, person_positions[parent]))
            # Sort by parent y-position
            sorted_parents.sort(key=lambda x: x[1])

            # Assign positions to children based on parent order and infection time
            for parent, parent_y in sorted_parents:
                if parent in parent_groups:
                    children = parent_groups[parent]
                    # Sort children by infection time
                    children.sort(key=lambda x: x[1])

                    for person_id, time in children:
                        person_positions[person_id] = current_y
                        current_y += 1

        # Add spacing between generations only after processing the generation
        if gen > 0 or len(people_with_times) > 1:
            current_y += 2

    return person_positions


def calculate_hierarchical_positions(infection_events, transmission_tree, infection_generations):
    """Hierarchical layout that minimizes edge crossings"""
    positions = {}
    event_dict = {e['person_id']: e for e in infection_events}

    # Start with patient zero at center
    patient_zero = min(infection_events, key=lambda x: x['time'])
    positions[patient_zero['person_id']] = 0

    # Process each generation
    for gen in sorted(set(infection_generations.values()))[1:]:
        gen_members = [p for p, g in infection_generations.items() if g == gen]

        # Group children by their parents
        parent_children = {}
        for person_id in gen_members:
            parent = None
            for potential_parent, children in transmission_tree.items():
                if person_id in children:
                    parent = potential_parent
                    break

            if parent not in parent_children:
                parent_children[parent] = []
            parent_children[parent].append(person_id)

        # Position children around their parents
        for parent, children in parent_children.items():
            if parent in positions:
                parent_pos = positions[parent]
                # Sort children by infection time
                children_with_times = [
                    (child, event_dict[child]['time']) for child in children]
                children_with_times.sort(key=lambda x: x[1])

                # Spread children symmetrically around parent
                num_children = len(children)
                for i, (child, _) in enumerate(children_with_times):
                    offset = (i - (num_children - 1) / 2) * 2
                    positions[child] = parent_pos + offset

    return positions


def calculate_location_positions(infection_events, transmission_tree, infection_generations):
    """Group by location type, then sort by time within each location"""
    location_groups = {}
    for event in infection_events:
        loc_type = event['location_type']
        if loc_type not in location_groups:
            location_groups[loc_type] = []
        location_groups[loc_type].append(event)

    positions = {}
    current_y = 0

    # Sort location groups by predefined order
    location_order = ['Patient Zero', 'EventPanvadere', 'Workplace Event', 'Restaurant Event',
                      'Work', 'Home', 'School', 'SocialEvent', 'BasicsShop']
    ordered_locations = []

    for loc in location_order:
        if loc in location_groups:
            ordered_locations.append((loc, location_groups[loc]))

    # Add any remaining locations not in the predefined order
    for loc, events in location_groups.items():
        if loc not in location_order:
            ordered_locations.append((loc, events))

    for loc_type, events in ordered_locations:
        events.sort(key=lambda x: x['time'])

        for i, event in enumerate(events):
            positions[event['person_id']] = current_y + i

        current_y += len(events) + 3  # Add spacing between location clusters

    return positions


def calculate_temporal_positions(infection_events, transmission_tree, infection_generations):
    """Temporal-spatial hybrid: group by time bins and location"""
    time_location_groups = {}

    for event in infection_events:
        time_bin = int(event['time'] // 24)  # Group by day
        loc_type = event['location_type']
        key = (time_bin, loc_type)

        if key not in time_location_groups:
            time_location_groups[key] = []
        time_location_groups[key].append(event)

    positions = {}
    current_y = 0

    for (time_bin, loc_type), events in sorted(time_location_groups.items()):
        events.sort(key=lambda x: x['time'])

        for i, event in enumerate(events):
            positions[event['person_id']] = current_y + i

        current_y += len(events) + 1

    return positions


# Connection drawing functions
def draw_l_shaped_connections(ax, transmission_tree, positions, infection_events, line_color, line_width, line_alpha, line_style):
    """Draw L-shaped connections (horizontal then vertical)"""
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_event = event_dict[infector_id]
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for infected_id in infected_list:
            infected_event = event_dict[infected_id]
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Draw horizontal first, then vertical
            ax.plot([infector_time, infected_time], [infector_y, infector_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)
            ax.plot([infected_time, infected_time], [infector_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)

            # Add arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time, infected_y -
                                (infected_y - infector_y) * 0.1),
                        arrowprops=dict(arrowstyle='->', color=line_color, lw=line_width, alpha=line_alpha))


def draw_current_connections(ax, transmission_tree, positions, infection_events, line_color, line_width, line_alpha, line_style):
    """Draw current connections (vertical then horizontal) - original method"""
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_event = event_dict[infector_id]
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for infected_id in infected_list:
            infected_event = event_dict[infected_id]
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Draw vertical first, then horizontal (original approach)
            ax.plot([infector_time, infector_time], [infector_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)
            ax.plot([infector_time, infected_time], [infected_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)

            # Add arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time - (infected_time -
                                infector_time) * 0.1, infected_y),
                        arrowprops=dict(arrowstyle='->', color=line_color, lw=line_width, alpha=line_alpha))


def draw_no_crossing_vertical_connections(ax, transmission_tree, positions, infection_events, line_color, line_width, line_alpha, line_style):
    """Draw simple vertical-then-horizontal connections with optimized y-positions to avoid crossings"""
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_event = event_dict[infector_id]
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for infected_id in infected_list:
            infected_event = event_dict[infected_id]
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Simple vertical then horizontal - no crossings due to optimized positioning
            ax.plot([infector_time, infector_time], [infector_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)
            ax.plot([infector_time, infected_time], [infected_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)

            # Add arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time - (infected_time -
                                infector_time) * 0.1, infected_y),
                        arrowprops=dict(arrowstyle='->', color=line_color, lw=line_width, alpha=line_alpha))


def draw_curved_connections(ax, transmission_tree, positions, infection_events, line_color, line_width, line_alpha, line_style):
    """Draw curved connections using Bezier curves"""
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_event = event_dict[infector_id]
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for i, infected_id in enumerate(infected_list):
            infected_event = event_dict[infected_id]
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Vary curve radius to avoid overlaps
            curve_rad = 0.3 + (i * 0.1)

            arrow = FancyArrowPatch(
                (infector_time, infector_y),
                (infected_time, infected_y),
                connectionstyle=f"arc3,rad={curve_rad}",
                arrowstyle='->',
                mutation_scale=20,
                color=line_color,
                alpha=line_alpha,
                linewidth=line_width,
                zorder=2
            )
            ax.add_patch(arrow)


def draw_straight_connections(ax, transmission_tree, positions, infection_events, line_color, line_width, line_alpha, line_style):
    """Draw straight diagonal connections"""
    event_dict = {e['person_id']: e for e in infection_events}

    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        infector_event = event_dict[infector_id]
        infector_y = positions[infector_id]
        infector_time = infector_event['time']

        for infected_id in infected_list:
            infected_event = event_dict[infected_id]
            infected_y = positions[infected_id]
            infected_time = infected_event['time']

            # Draw straight line
            ax.plot([infector_time, infected_time], [infector_y, infected_y],
                    color=line_color, linewidth=line_width, alpha=line_alpha,
                    linestyle=line_style, zorder=2)

            # Add arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infector_time + (infected_time - infector_time) * 0.9,
                                infector_y + (infected_y - infector_y) * 0.9),
                        arrowprops=dict(arrowstyle='->', color=line_color, lw=line_width, alpha=line_alpha))


def add_location_pie_inset(ax, infection_events, location_colors):
    """Add a small pie chart in the lower right corner of the tree plot"""

    # Count infections by location type (excluding Patient Zero)
    counts = {}
    for event in infection_events:
        loc_type = event['location_type']
        if loc_type != 'Patient Zero':  # Skip Patient Zero
            counts[loc_type] = counts.get(loc_type, 0) + 1

    if not counts:
        return  # No data to show

    total = sum(counts.values())

    # Prepare data for pie chart
    sizes = []
    colors = []
    labels = []

    # Define consistent order: Home, Work, School, SocialEvent, BasicsShop
    location_order = ['Home', 'Work', 'School', 'SocialEvent', 'BasicsShop']

    # Add locations in consistent order (only if they have counts)
    for loc_type in location_order:
        if loc_type in counts:
            count = counts[loc_type]
            percentage = (count / total) * 100
            sizes.append(percentage)
            colors.append(location_colors.get(loc_type, '#999999'))
            labels.append(f'{percentage:.1f}%')

    # Add any remaining locations not in the predefined order
    remaining_locations = set(counts.keys()) - set(location_order)
    for loc_type in sorted(remaining_locations):
        count = counts[loc_type]
        percentage = (count / total) * 100
        sizes.append(percentage)
        colors.append(location_colors.get(loc_type, '#999999'))
        labels.append(f'{percentage:.1f}%')

    # Create inset axes for pie chart in lower right corner (moved lower)
    # Position: [left, bottom, width, height] in axes coordinates
    pie_ax = ax.inset_axes([0.75, -0.05, 0.23, 0.25])

    # Create pie chart with no labels
    wedges = pie_ax.pie(sizes, colors=colors, startangle=90,
                        labels=None, labeldistance=0)

    # Remove axes and make it clean
    pie_ax.set_aspect('equal')


def plot_infection_tree_enhanced(ax, infection_events, transmission_tree, positions, location_colors,
                                 title, generations, infected_status_df, contact_analysis_df,
                                 layout_method='current', connection_style='current'):
    """Enhanced plot_infection_tree function with layout and connection options"""

    times = [event['time'] for event in infection_events]
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values) if y_values else 1

    # Connection drawing parameters
    connection_methods = {
        'current': draw_current_connections,
        'no-crossings': draw_no_crossing_vertical_connections,
        'L-shaped': draw_l_shaped_connections,
        'curved': draw_curved_connections,
        'straight': draw_straight_connections
    }

    # Draw transmission tree connections
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
            potential_infectors = 0
            if infected_status_df is not None and contact_analysis_df is not None:
                patient_zero_time = min(e['time'] for e in infection_events)

                for e in infection_events:
                    if e['time'] < infected_event['time']:
                        if has_contact_at_location(e['person_id'], infected_event['person_id'],
                                                   infected_event['location'], infected_event['time'], contact_analysis_df):
                            is_patient_zero = (e['time'] == patient_zero_time)

                            if is_patient_zero:
                                potential_infectors += 1
                            else:
                                if infected_event['time'] < 0:
                                    continue

                                infectious_check = infected_status_df[
                                    (infected_status_df['Person_ID'] == e['person_id']) &
                                    (infected_status_df['Timestep'] == int(infected_event['time'])) &
                                    (infected_status_df['Infected'] == 1)
                                ]

                                if len(infectious_check) > 0:
                                    potential_infectors += 1
            else:
                potential_infectors = sum(1 for e in infection_events
                                          if e['time'] < infected_event['time']
                                          and e['location_type'] == infected_event['location_type'])

            line_width = 2
            line_alpha = 0.9 if same_location else 0.6
            line_color = 'black'
            line_style = '--' if potential_infectors > 1 else '-'

            # Draw connections using selected method
            mini_tree = {infector_id: [infected_id]}
            connection_methods[connection_style](ax, mini_tree, positions, infection_events,
                                                 line_color, line_width, line_alpha, line_style)

    # Plot infection points
    for event in infection_events:
        y_pos = positions[event['person_id']]
        location_type = event['location_type']
        if location_type == 'EventPanvadere':
            color = location_colors.get('EventPanvadere', '#FF1493')
        else:
            color = location_colors.get(location_type, '#FF1493')
        size = 300 if event['time'] == min(times) else 200
        border_width = 1

        ax.scatter(event['time'], y_pos, s=size, c=color, alpha=1.0,
                   zorder=5, edgecolors='black', linewidth=border_width)

    # Add timeline markers
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values) if y_values else 1
    outbreak_time = min(times)

    # Add simulation start line
    ax.axvline(x=0, color='blue', linestyle='-',
               alpha=0.8, linewidth=3, zorder=1)
    ax.text(0, -y_range * 0.12, 'Simulation Start',
            ha='center', va='top', fontsize=14, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='blue', alpha=0.7, edgecolor='darkblue'))

    if len(times) > 1:
        sorted_times = sorted(times)
        outbreak_event_time = sorted_times[min(1, len(sorted_times)-1)]

        # Add vertical line for patient zero
        ax.axvline(x=outbreak_time, color='red', linestyle='--',
                   alpha=0.7, linewidth=3, zorder=1)
        ax.text(outbreak_time, -y_range * 0.02, 'Patient Zero',
                ha='center', va='top', fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.7, edgecolor='darkred'))

        # Outbreak event line
        outbreak_event = next(
            (e for e in infection_events if e['time'] == outbreak_event_time), None)
        if outbreak_event:
            loc_type = outbreak_event.get('location_type', 'EventPanvadere')
            outbreak_color = location_colors.get(loc_type, '#FF1493')
        else:
            outbreak_color = '#FF1493'

        ax.axvline(x=outbreak_event_time, color=outbreak_color,
                   linestyle='--', alpha=0.7, linewidth=3, zorder=1)
        ax.text(outbreak_event_time, -y_range * 0.07, 'Outbreak Event',
                ha='center', va='top', fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=outbreak_color, alpha=0.7, edgecolor=outbreak_color))

    # Formatting
    ax.set_xlabel('Timestep', fontsize=18)
    ax.set_ylabel('Transmission Order', fontsize=18)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=14)

    # Set limits with padding
    ax.set_ylim(min(y_values) - y_range * 0.15, max(y_values) + y_range * 0.1)

    # Hide negative x-axis tick labels but keep the axis range
    x_ticks = ax.get_xticks()
    tick_labels = [str(int(tick)) if tick >= 0 else '' for tick in x_ticks]
    ax.set_xticklabels(tick_labels)

    # Add location distribution pie chart in lower right corner
    add_location_pie_inset(ax, infection_events, location_colors)


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

    # Estimate missed opportunities
    infectious_people = set()
    for event in infection_events:
        infectious_people.add(event['person_id'])

    potential_transmissions = 0
    for (timestep, location), group in contact_data.groupby(['Timestep', 'Location_ID']):
        people = group['Person_ID'].tolist()
        infectious_present = any(
            person_id in infectious_people for person_id in people)

        if infectious_present and len(people) > 1:
            n = len(people)
            potential_transmissions += n * (n - 1) // 2

    opportunities['total_contacts'] = potential_transmissions
    opportunities['missed_opportunities'] = max(
        0, potential_transmissions - opportunities['realized_transmissions'])

    if potential_transmissions > 0:
        opportunities['transmission_efficiency'] = opportunities['realized_transmissions'] / \
            potential_transmissions

    return opportunities


def calculate_clustering_metrics(infection_events, transmission_tree):
    """Calculate clustering-related metrics for transmission patterns"""
    metrics = {
        'same_location_transmissions': 0,
        'total_transmissions': 0,
        'clustering_coefficient': 0.0,
        'location_clusters': {}
    }

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

            if infected_location not in metrics['location_clusters']:
                metrics['location_clusters'][infected_location] = 0
            metrics['location_clusters'][infected_location] += 1

    if metrics['total_transmissions'] > 0:
        metrics['clustering_coefficient'] = metrics['same_location_transmissions'] / \
            metrics['total_transmissions']

    return metrics


def create_comparative_infection_trees_enhanced(data_dir_method1, data_dir_method2,
                                                method1_name="Transmission-Informed",
                                                method2_name="Uniform",
                                                max_infections=50, output_path=None,
                                                scenario_name=None, layout_method='current',
                                                connection_style='current'):
    """Create side-by-side infection tree comparison with enhanced layout options"""

    print(
        f"Creating enhanced comparative infection trees for {scenario_name or 'scenario'}...")
    print(f"Layout method: {layout_method}")
    print(f"Connection style: {connection_style}")

    # Load data for both methods (same as original)
    # Load data for both methods (same as original)
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
    infected_status_df_1 = pd.read_csv(infected_status_file_1)

    # Load and process data for method 2
    print(f"Loading {method2_name} data...")
    contact_df_2, infection_df_2 = load_simulation_data(
        contact_file_2, infection_file_2)
    infection_events_2 = transform_infection_data(infection_df_2)
    contact_analysis_df_2 = create_contact_data_for_analysis(
        contact_df_2, infection_df_2)
    infected_status_df_2 = pd.read_csv(infected_status_file_2)

    # Limit infections to first 5 days for visualization clarity
    max_timesteps = 120
    inf_events_for_opp_1 = infection_events_1
    inf_events_for_opp_2 = infection_events_2
    infection_events_1 = [
        event for event in infection_events_1 if event['time'] <= max_timesteps]
    infection_events_2 = [
        event for event in infection_events_2 if event['time'] <= max_timesteps]

    # Update location names based on scenario
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

    # Build transmission trees
    print("Building transmission trees...")
    tree_1, generations_1, potential_1 = build_transmission_tree(
        infection_events_1, contact_analysis_df_1, infected_status_df_1)
    tree_2, generations_2, potential_2 = build_transmission_tree(
        infection_events_2, contact_analysis_df_2, infected_status_df_2)

    # Calculate positions using selected layout method
    layout_methods = {
        'current': calculate_current_positions,
        'no-crossings': calculate_current_positions_no_crossings,
        'transmission-order': calculate_transmission_order_positions,
        'depth-first': calculate_depth_first_positions,
        'hierarchical': calculate_hierarchical_positions,
        'location': calculate_location_positions,
        'temporal': calculate_temporal_positions
    }

    if layout_method in ['current', 'no-crossings', 'transmission-order', 'depth-first', 'hierarchical']:
        positions_1 = layout_methods[layout_method](
            infection_events_1, tree_1, generations_1)
        positions_2 = layout_methods[layout_method](
            infection_events_2, tree_2, generations_2)
    else:
        positions_1 = layout_methods[layout_method](
            infection_events_1, tree_1, generations_1)
        positions_2 = layout_methods[layout_method](
            infection_events_2, tree_2, generations_2)

    # Create visualization
    fig = plt.figure(figsize=(20, 18))
    gs = fig.add_gridspec(3, 2, height_ratios=[
                          2.5, 2.5, 1.5], hspace=0.3, wspace=0.2)

    # Location colors
    base_location_colors = {
        'Work': '#4169E1',
        'Home': '#2E8B57',
        'School': '#FF6347',
        'SocialEvent': '#9370DB',
        'BasicsShop': '#FF8C00',
        'Patient Zero': "#FF0000",
        'EventPanvadere': '#FF1493'
    }

    location_colors = base_location_colors.copy()
    if scenario_name and ('W' in scenario_name or 'work' in scenario_name.lower()):
        location_colors['EventPanvadere'] = '#8A2BE2'
        location_colors['Workplace Event'] = '#8A2BE2'
    else:
        location_colors['EventPanvadere'] = '#FF1493'
        location_colors['Restaurant Event'] = '#FF1493'

    # Plot trees with enhanced layout
    ax1 = fig.add_subplot(gs[0:2, 0])
    plot_infection_tree_enhanced(ax1, infection_events_1, tree_1, positions_1, location_colors,
                                 f"{method1_name}\n{len(infection_events_1)} infections",
                                 generations_1, infected_status_df_1, contact_analysis_df_1,
                                 layout_method, connection_style)

    ax2 = fig.add_subplot(gs[0:2, 1])
    plot_infection_tree_enhanced(ax2, infection_events_2, tree_2, positions_2, location_colors,
                                 f"{method2_name}\n{len(infection_events_2)} infections",
                                 generations_2, infected_status_df_2, contact_analysis_df_2,
                                 layout_method, connection_style)

    # Bottom panel - transmission opportunities analysis
    ax3 = fig.add_subplot(gs[2, :])
    plot_transmission_opportunities_over_time(ax3, inf_events_for_opp_1, inf_events_for_opp_2,
                                              method1_name, method2_name,
                                              contact_analysis_df_1, contact_analysis_df_2,
                                              infected_status_df_1, infected_status_df_2)

    # Add overall title

    title = f"Comparative Infection Transmission Trees"
    scenarios = {
        'R1_restaurant_strong_clustering': 'R1',
        'R2_restaurant_weaker_clustering': 'R2',
        'W1_workplace_few_meetings': 'W1',
        'W2_workplace_many_meetings': 'W2'
    }
    if scenario_name:
        title = f"Scenario {scenario_name}: {title}"

    fig.suptitle(f"{title}\n{method1_name} vs {method2_name} Initialization\n",
                 fontsize=16, fontweight='bold', y=0.95)

    # Create legend
    legend_elements = []
    for loc_type, color in location_colors.items():
        if loc_type not in ['EventPanvadere', 'Workplace Event', 'Restaurant Event', 'Patient Zero']:
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=color, markersize=10,
                                              label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='black', lw=2, linestyle='-',
                   label='Certain transmission', alpha=0.9),
        plt.Line2D([0], [0], color='black', lw=2, linestyle='--',
                   label='Uncertain transmission', alpha=0.9)
    ])

    ax1.legend(handles=legend_elements, loc='upper left',
               title_fontsize=10, bbox_to_anchor=(0, 1), fontsize=14)

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=500,
                    bbox_inches='tight', facecolor='white')
        print(f"Enhanced comparative infection trees saved to {output_path}")

    # Print comparison summary
    print(
        f"\n=== COMPARATIVE TRANSMISSION ANALYSIS ({layout_method.upper()} LAYOUT) ===")
    print(f"{method1_name}:")
    print(f"  Total infections: {len(infection_events_1)}")

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

    infection_diff = len(infection_events_2) - len(infection_events_1)
    clustering_diff = clustering_2['clustering_coefficient'] - \
        clustering_1['clustering_coefficient']
    print(f"\nDIFFERENCES ({method2_name} - {method1_name}):")
    print(f"  Total infections: {infection_diff:+d}")
    print(f"  Clustering coefficient: {clustering_diff:+.3f}")

    return fig


def plot_transmission_opportunities_over_time(ax, infection_events_1, infection_events_2, method1_name, method2_name,
                                              contact_data_1, contact_data_2, infected_status_df_1, infected_status_df_2):
    """
    Plot transmission opportunities over time based on actual contact data.
    """
    times_1 = [event['time'] for event in infection_events_1]
    times_2 = [event['time'] for event in infection_events_2]
    all_times = sorted(set(times_1 + times_2))

    if not all_times:
        return

    time_range = np.arange(0, 240, 1)

    def calculate_opportunities(time_point, infection_events, contact_data, infected_status_df):
        """Calculate transmission opportunities at a specific timepoint"""
        infectious_agents = set()

        for event in infection_events:
            if event['time'] <= time_point:
                person_id = event['person_id']
                infectious_check = infected_status_df[
                    (infected_status_df['Person_ID'] == person_id) &
                    (infected_status_df['Timestep'] == int(time_point)) &
                    (infected_status_df['Infected'] == 1)
                ]

                if len(infectious_check) > 0:
                    infectious_agents.add(person_id)

        if not infectious_agents:
            return 0

        infectious_locations = set()
        current_contacts = contact_data[contact_data['Timestep'] == int(
            time_point)]

        for _, contact in current_contacts.iterrows():
            person_id = contact['Person_ID']
            location_id = contact['Location_ID']

            if person_id in infectious_agents:
                infectious_locations.add(location_id)

        transmission_opportunities = 0
        for location_id in infectious_locations:
            people_at_location = current_contacts[
                current_contacts['Location_ID'] == location_id
            ]['Person_ID'].tolist()

            susceptible_at_location = [
                person_id for person_id in people_at_location
                if person_id not in infectious_agents
            ]

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

    # Fill area between curves
    ax.fill_between(time_range, opportunities_1, opportunities_2,
                    alpha=0.3, color='lightgray', label='Difference')

    # Calculate cumulative differences
    timesteps_per_day = 24
    day_indices = [min(i * timesteps_per_day, len(opportunities_1) - 1)
                   for i in [1, 3, 5, 7]]

    cum_opp_1 = np.cumsum(opportunities_1)
    cum_opp_2 = np.cumsum(opportunities_2)
    cum_opp_rel = (cum_opp_2 - cum_opp_1)
    cum_opp_perc = (cum_opp_rel / cum_opp_1 * 100)

    day_diffs = [cum_opp_2[idx] - cum_opp_1[idx] if idx < len(cum_opp_1) else 0
                 for idx in day_indices]

    # Add day boundary markers
    day_boundaries = [24, 72, 120, 168]
    day_numbers = [1, 3, 5, 7]
    colors = ['green', 'orange', 'purple', 'blue']

    max_opp = max(max(opportunities_1), max(opportunities_2)
                  ) if opportunities_1 and opportunities_2 else 100
    y_positions = [max_opp * 0.9] * 4

    for i, (boundary, day_num, diff, color) in enumerate(zip(day_boundaries, day_numbers, day_diffs, colors)):
        if boundary <= max(time_range):
            ax.axvline(x=boundary, color=color,
                       linestyle=':', alpha=0.7, linewidth=3)

            if opportunities_1 and opportunities_2:
                diff_text = f"Day {day_num}: {diff:+.0f}"
                ax.text(boundary - 16, y_positions[i], diff_text, rotation=0, fontsize=16,
                        bbox=dict(boxstyle='round,pad=0.3',
                                  facecolor=color, alpha=0.4),
                        ha='left', va='center', fontweight='bold')

    # Formatting
    ax.set_xlabel('Timestep', fontsize=18)
    ax.set_ylabel('Susceptible People at Risk', fontsize=18)
    ax.set_title(
        f'Transmission Opportunities Over Time\n'
        f'Day X: Cumulative Difference in At-Risk Population\n'
        f'Final Cumulative Difference: {cum_opp_rel[-1]:.0f} ({cum_opp_perc[-1]:.1f}%) ({cum_opp_1[-1]:.0f} vs {cum_opp_2[-1]:.0f})',
        fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.7)
    ax.legend(loc='upper right', fontsize=18)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.set_xlim(0, max(time_range) + 1)


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description='Create enhanced comparative infection tree visualization')
    parser.add_argument('--data-dir-method1', required=True,
                        help='Directory containing transmission-informed data files')
    parser.add_argument('--data-dir-method2', required=True,
                        help='Directory containing uniform initialization data files')
    parser.add_argument('--method1-name', default='Transmission-Informed',
                        help='Name for first method')
    parser.add_argument('--method2-name', default='Uniform',
                        help='Name for second method')
    parser.add_argument('--layout-method', choices=['current', 'no-crossings', 'transmission-order', 'depth-first', 'hierarchical', 'location', 'temporal'],
                        default='current', help='Y-axis layout method')
    parser.add_argument('--connection-style', choices=['current', 'no-crossings', 'L-shaped', 'curved', 'straight'],
                        default='current', help='Connection drawing style')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--scenario-name',
                        help='Scenario name for titles (e.g., R1, W2)')
    parser.add_argument('--max-infections', type=int, default=150,
                        help='Maximum infections to display for clarity')

    args = parser.parse_args()

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
        base_name += f"_{args.layout_method}_{args.connection_style.replace('-', '_')}"
        args.output_path = f"{base_name}.png"

    # Create visualization
    fig = create_comparative_infection_trees_enhanced(
        args.data_dir_method1, args.data_dir_method2,
        args.method1_name, args.method2_name,
        args.max_infections, args.output_path, args.scenario_name,
        args.layout_method, args.connection_style
    )

    if fig is None:
        print("Failed to create visualization due to data issues.")
        sys.exit(1)

    print("Enhanced comparative infection tree analysis complete!")


if __name__ == "__main__":
    main()
