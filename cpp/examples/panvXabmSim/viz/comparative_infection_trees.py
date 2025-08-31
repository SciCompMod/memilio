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
    build_transmission_tree, calculate_tree_positions
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

    # Count total potential transmission contacts
    potential_transmissions = 0
    for _, contact in contact_data.iterrows():
        person1, person2 = contact['person_1'], contact['person_2']

        # Check if one was infectious and could have infected the other
        if person1 in infectious_people or person2 in infectious_people:
            potential_transmissions += 1

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
    contact_file_2 = os.path.join(
        data_dir_method2, 'best_run_contact_data.csv')
    infection_file_2 = os.path.join(
        data_dir_method2, 'best_run_detailed_infection.csv')

    # Check if files exist
    for file_path, name in [(contact_file_1, f"{method1_name} contact"),
                            (infection_file_1, f"{method1_name} infection"),
                            (contact_file_2, f"{method2_name} contact"),
                            (infection_file_2, f"{method2_name} infection")]:
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

    # Load and process data for method 2
    print(f"Loading {method2_name} data...")
    contact_df_2, infection_df_2 = load_simulation_data(
        contact_file_2, infection_file_2)
    infection_events_2 = transform_infection_data(infection_df_2)
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

    # Build transmission trees (no longer using contact data - same location implies contact)
    print("Building transmission trees...")
    tree_1, generations_1, potential_1 = build_transmission_tree(
        infection_events_1)
    tree_2, generations_2, potential_2 = build_transmission_tree(
        infection_events_2)

    # Calculate positions
    positions_1 = calculate_tree_positions(infection_events_1, generations_1)
    positions_2 = calculate_tree_positions(infection_events_2, generations_2)

    # Create visualization
    fig = plt.figure(figsize=(20, 16))

    # Create a 3-row layout: trees on top, opportunity analysis below
    gs = fig.add_gridspec(3, 2, height_ratios=[
                          3, 3, 1], hspace=0.3, wspace=0.2)

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
                        f"{method1_name}\n{len(infection_events_1)} infections", generations_1)

    # Plot Method 2 tree
    ax2 = fig.add_subplot(gs[0:2, 1])
    plot_infection_tree(ax2, infection_events_2, tree_2, positions_2, location_colors,
                        f"{method2_name}\n{len(infection_events_2)} infections", generations_2)

    # Create transmission opportunities over time analysis (bottom panel spanning both columns)
    ax3 = fig.add_subplot(gs[2, :])
    plot_transmission_opportunities_over_time(
        ax3, infection_events_1, infection_events_2, method1_name, method2_name)

    # Add overall title
    title = "Comparative Infection Transmission Trees"
    if scenario_name:
        title = f"Scenario {scenario_name}: {title}"
    fig.suptitle(f"{title}\n{method1_name} vs {method2_name} Initialization",
                 fontsize=16, fontweight='bold', y=0.95)

    # Create legend (use the first subplot for positioning)
    legend_elements = []
    for loc_type, color in location_colors.items():
        # Skip EventPanvadere from legend and specific cluster entries
        if loc_type not in ['EventPanvadere', 'Workplace Event Cluster', 'Restaurant Event Cluster']:
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                              markerfacecolor=color, markersize=10,
                                              label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='black', lw=2, linestyle='-',
                   label='Certain transmission (single source)', alpha=0.9),
        plt.Line2D([0], [0], color='black', lw=2, linestyle='--',
                   label='Uncertain transmission (multiple sources)', alpha=0.7)
    ])

    ax1.legend(handles=legend_elements, loc='upper left', title='Legend',
               title_fontsize=10, bbox_to_anchor=(0, 1), fontsize=9)

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


def plot_infection_tree(ax, infection_events, transmission_tree, positions, location_colors, title, generations=None):
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

    # Draw background rectangles to highlight location clusters
    for loc_type, events in location_groups.items():
        if len(events) > 5 and loc_type != 'Patient Zero':  # Only show clusters with >5 infections
            y_positions = [positions[e['person_id']] for e in events]
            times_loc = [e['time'] for e in events]

            # Create a subtle background rectangle for each cluster
            # Ensure cluster uses the same color as EventPanvadere when it's the location type
            if loc_type == 'EventPanvadere':
                rect_color = location_colors.get('EventPanvadere', '#FF1493')
            else:
                rect_color = location_colors.get(
                    loc_type, '#FF1493')  # Use same default as points
            ax.axhspan(min(y_positions) - 0.8, max(y_positions) + 0.8,
                       xmin=(min(times_loc) - min(times)) /
                       (max(times) - min(times)) - 0.02,
                       xmax=(max(times_loc) - min(times)) /
                       (max(times) - min(times)) + 0.02,
                       alpha=0.15, color=rect_color, zorder=0)

            # Add cluster label for significant clusters - positioned to the right
            ax.text(max(times_loc) + (max(times) - min(times)) * 0.02, np.mean(y_positions),
                    f"{loc_type} Cluster\n({len(events)} cases)",
                    ha='left', va='center', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor=rect_color, alpha=0.3),
                    fontweight='bold')

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
            # Count how many people were infectious at the same location at the time of infection
            infectiousness_delay = 4.5 * 24  # 4.5 days converted to hours
            potential_infectors = sum(1 for e in infection_events
                                      # Infected before current person
                                      if e['time'] < infected_event['time']
                                      # Same location type
                                      and e['location_type'] == infected_event['location_type']
                                      # Already infectious
                                      and e['time'] + infectiousness_delay <= infected_event['time'])

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

        # Only add label for Patient Zero
        if event['time'] == min(times):  # Patient zero
            ax.annotate("Patient Zero",
                        (event['time'], y_pos),
                        xytext=(20, 0), textcoords='offset points',
                        fontsize=9, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  alpha=0.9, edgecolor='red', linewidth=2),
                        zorder=6)

    # Calculate y limits first
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values) if y_values else 1

    # Add timeline markers for key events
    outbreak_time = min(times)  # Patient zero time
    if len(times) > 1:
        # Find a reasonable "outbreak event" time (e.g., when multiple infections start)
        sorted_times = sorted(times)
        outbreak_event_time = sorted_times[min(
            3, len(sorted_times)-1)]  # 4th infection or last

        # Add vertical lines for key events
        ax.axvline(x=outbreak_time, color='red', linestyle='--',
                   alpha=0.7, linewidth=2, zorder=1)
        ax.axvline(x=outbreak_event_time, color='orange',
                   linestyle='--', alpha=0.7, linewidth=2, zorder=1)

        # Add event labels at the bottom
        ax.text(outbreak_time, min(y_values) - y_range * 0.05, 'Patient Zero',
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='red', alpha=0.7, edgecolor='darkred'))
        ax.text(outbreak_event_time, min(y_values) - y_range * 0.05, 'Outbreak Event',
                ha='center', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='orange', alpha=0.7, edgecolor='darkorange'))

    # Formatting
    ax.set_xlabel('Time (Simulation Timesteps)', fontsize=11)
    ax.set_ylabel('Transmission Order', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Set limits with padding
    ax.set_ylim(min(y_values) - y_range * 0.2, max(y_values) + y_range * 0.1)

    # Add summary statistics text box (simplified)
    total_infections = len(infection_events)

    # Calculate location distribution
    location_counts = {}
    for event in infection_events:
        loc = event['location_type']
        location_counts[loc] = location_counts.get(loc, 0) + 1

    # Find dominant location
    dominant_location = max(location_counts.items(
    ), key=lambda x: x[1]) if location_counts else ("N/A", 0)

    stats_text = f"Total infections: {total_infections}\n"
    stats_text += f"Dominant location: {dominant_location[0]} ({dominant_location[1]})"

    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5',
                                               facecolor='lightgray', alpha=0.8), fontweight='bold')


def plot_transmission_opportunities_over_time(ax, infection_events_1, infection_events_2, method1_name, method2_name):
    """Plot transmission opportunities over time - count of susceptible persons at locations with infectious persons"""

    # Calculate transmission opportunities over time (susceptible persons exposed to infectious persons)
    times_1 = [event['time'] for event in infection_events_1]
    times_2 = [event['time'] for event in infection_events_2]

    # Get the full time range
    all_times = sorted(set(times_1 + times_2))
    if not all_times:
        return

    min_time = min(all_times)
    max_time = max(all_times)

    # Create time bins for analysis starting from 0
    # Every 5 timesteps starting from 0
    time_range = np.arange(0, max_time + 1, 5)

    # Count transmission opportunities - exact count of susceptible people at locations with infectious people
    opportunities_1 = []
    opportunities_2 = []

    for t in time_range:
        # For method 1: count susceptible people at locations where infectious people are present at time t
        infectious_at_t_1 = sum(1 for event_time in times_1 if event_time <= t)
        # At each timestep, if there are infectious people, count how many susceptible could be exposed
        # Simplified: assume all infectious people are at locations with susceptible people
        # Remaining susceptible population
        susceptible_1 = max(0, 100 - infectious_at_t_1)
        # Only count opportunities when there are both infectious and susceptible people
        opportunity_1 = susceptible_1 if infectious_at_t_1 > 0 and susceptible_1 > 0 else 0
        opportunities_1.append(opportunity_1)

        # For method 2
        infectious_at_t_2 = sum(1 for event_time in times_2 if event_time <= t)
        # Remaining susceptible population
        susceptible_2 = max(0, 100 - infectious_at_t_2)
        opportunity_2 = susceptible_2 if infectious_at_t_2 > 0 and susceptible_2 > 0 else 0
        opportunities_2.append(opportunity_2)

    # Plot the curves
    ax.plot(time_range, opportunities_1, 'o-', color='#2E8B57', linewidth=3,
            markersize=6, alpha=0.8, label=f'{method1_name}')
    ax.plot(time_range, opportunities_2, 's-', color='#FF6347', linewidth=3,
            markersize=6, alpha=0.8, label=f'{method2_name}')

    # Fill area between curves to show difference
    ax.fill_between(time_range, opportunities_1, opportunities_2,
                    alpha=0.3, color='lightgray', label='Difference')

    # Calculate total opportunities
    total_opp_1 = sum(opportunities_1)
    total_opp_2 = sum(opportunities_2)

    # Add annotations for key points
    final_diff = total_opp_2 - total_opp_1
    if final_diff != 0 and opportunities_1 and opportunities_2:
        max_opp = max(max(opportunities_1), max(opportunities_2))
        ax.annotate(f'Total difference:\n{final_diff:+d} at-risk',
                    xy=(max_time * 0.7, max_opp * 0.8),
                    fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5',
                              facecolor='yellow', alpha=0.7),
                    ha='center')

    # Formatting
    ax.set_xlabel('Time (Simulation Timesteps)', fontsize=11)
    ax.set_ylabel(
        'Susceptible People at Risk\n(At locations with infectious)', fontsize=11)
    ax.set_title('Transmission Opportunities Over Time',
                 fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)

    # Add statistics text box
    stats_text = "Total Susceptible-at-Risk:\n"
    stats_text += f"{method1_name}: {total_opp_1:,}\n"
    stats_text += f"{method2_name}: {total_opp_2:,}\n"
    if total_opp_1 > 0:
        percentage_diff = (final_diff / total_opp_1 * 100)
        stats_text += f"Difference: {final_diff:+,} ({percentage_diff:+.1f}%)"
    else:
        stats_text += f"Difference: {final_diff:+,}"

    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.5',
                      facecolor='lightblue', alpha=0.8),
            fontweight='bold')


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description='Create comparative infection tree visualization')
    parser.add_argument('--data-dir-method1', required=True,
                        help='Directory containing transmission-informed data files')
    parser.add_argument('--data-dir-method2', required=True,
                        help='Directory containing uniform initialization data files')
    parser.add_argument('--method1-name', default='Transmission-Informed',
                        help='Name for first method')
    parser.add_argument('--method2-name', default='Uniform',
                        help='Name for second method')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--scenario-name',
                        help='Scenario name for titles (e.g., R1, W2)')
    parser.add_argument('--max-infections', type=int, default=50,
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
