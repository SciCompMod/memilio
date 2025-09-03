

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import argparse
import os
import sys
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


def get_infection_status_color(person_id, infected_at_time):
    """Get color based on infection status at specific time

    Args:
        person_id: ID of the person
        infected_at_time: Set of person IDs infected at this time point
    """
    if person_id in infected_at_time:
        return '#FF4444', '#8B0000'  # Red with dark red border for infected
    else:
        return '#FFFFFF', '#CCCCCC'  # White with light gray border for not infected


def create_time_specific_heatmap(ax, household_data, household_positions, person_pos,
                                 infected_at_time, timestep, title):
    """Create heatmap for specific time point"""

    # Calculate household infection stats for this time point
    household_stats = defaultdict(
        lambda: {'total_members': 0, 'infected_members': 0})

    for _, row in household_data.iterrows():
        person_id = row['person_id']
        household_id = row['household_id']

        household_stats[household_id]['total_members'] += 1

        if person_id in infected_at_time:
            household_stats[household_id]['infected_members'] += 1

    # Calculate infection rates
    for household_id in household_stats:
        stats = household_stats[household_id]
        if stats['total_members'] > 0:
            stats['infection_rate'] = stats['infected_members'] / \
                stats['total_members']

    # Create the graph
    G = nx.Graph()
    for person_id in person_pos.keys():
        G.add_node(person_id, node_type='person')

    # Create node colors based on infection status at this time
    node_colors = []
    edge_colors = []
    node_sizes = []

    for person_id in person_pos.keys():
        node_color, edge_color = get_infection_status_color(
            person_id, infected_at_time)

        node_colors.append(node_color)
        edge_colors.append(edge_color)

        # Slightly larger nodes for infected individuals
        if person_id in infected_at_time:
            node_sizes.append(100)
        else:
            node_sizes.append(80)

    # Draw household boundaries
    for household_id, position in household_positions.items():
        if not isinstance(position, tuple):
            continue

        hx, hy = position
        household_members = household_data[household_data['household_id']
                                           == household_id]['person_id'].tolist()

        if len(household_members) > 0:
            member_positions = [person_pos[pid]
                                for pid in household_members if pid in person_pos]
            if member_positions:
                xs, ys = zip(*member_positions)
                min_x, max_x = min(xs) - 25, max(xs) + 25
                min_y, max_y = min(ys) - 25, max(ys) + 25

                stats = household_stats.get(household_id, {
                                            'infection_rate': 0, 'infected_members': 0, 'total_members': len(household_members)})

                # Color household boundary based on infection rate at this time point only
                if stats['infection_rate'] == 0:
                    boundary_color = '#FFFFFF'  # White for not infected
                    alpha = 0.2
                elif stats['infection_rate'] < 1.0:
                    boundary_color = '#FFA500'  # Orange for partially infected
                    alpha = 0.4
                else:
                    boundary_color = '#FF4444'  # Red for fully infected
                    alpha = 0.5

                # Draw household boundary
                width = max_x - min_x
                height = max_y - min_y
                rect = plt.Rectangle((min_x, min_y), width, height,
                                     facecolor=boundary_color, alpha=alpha, zorder=0,
                                     edgecolor='gray', linewidth=1)
                ax.add_patch(rect)

    # Draw person nodes
    nx.draw_networkx_nodes(G, person_pos,
                           node_color=node_colors,
                           node_size=node_sizes,
                           alpha=0.8,
                           edgecolors=edge_colors,
                           linewidths=1.5,
                           ax=ax)

    # Set title and formatting
    infected_count = len(infected_at_time)
    total_people = len(person_pos)
    if timestep == 0:
        day_title = "Initialization"
    else:
        day_title = f"Day {timestep}"
    ax.set_title(f"{title}\n{day_title}: {infected_count}/{total_people} infected ({infected_count/total_people*100:.1f}%)",
                 fontsize=12, fontweight='bold', pad=10)
    ax.axis('equal')
    ax.axis('off')

    return infected_count


# 0, 5, 10 days in hours
def create_temporal_infection_heatmap(data_dir, output_path=None, time_points=[0, 120, 240]):
    """Create time-progressive household infection heatmap"""
    print("Loading temporal infection data...")

    # Load data
    try:
        household_data = load_household_data(
            os.path.join(data_dir, 'household_id.csv'))
        contact_data = load_contact_data(
            os.path.join(data_dir, 'contact_intensiveness.csv'))
        location_data = load_location_data(
            os.path.join(data_dir, 'location_id_and_type.csv'))
        infections_by_time = load_temporal_infection_data(
            os.path.join(data_dir, 'best_run_detailed_infection.csv'))
    except FileNotFoundError as e:
        print(f"Error: Required file not found: {e}")
        return

    print(
        f"Loaded data: {len(household_data)} people, temporal infections for {len(infections_by_time)} time steps")

    # Create household layout (same for all time points)
    spacing_factor = 0.2
    household_spacing = 300.0 * spacing_factor
    person_spacing = 50.0 * spacing_factor

    person_pos, households, household_positions = create_household_layout(
        household_data,
        household_spacing=household_spacing,
        person_spacing=person_spacing,
        max_households_per_row=25
    )

    # Create figure with subplots for each time point
    fig, axes = plt.subplots(
        1, len(time_points), figsize=(6*len(time_points), 8))
    if len(time_points) == 1:
        axes = [axes]

    infection_counts = []

    for i, timestep in enumerate(time_points):
        ax = axes[i]

        # Get infected people at this time point
        infected_at_time = infections_by_time.get(timestep, set())

        # Create heatmap for this time point
        day_label = timestep // 24  # Convert hours to days
        infected_count = create_time_specific_heatmap(
            ax, household_data, household_positions, person_pos,
            infected_at_time, day_label,
            f"Time Point {i+1}"
        )
        infection_counts.append(infected_count)

    # Create shared legend
    legend_elements = [
        plt.Circle((0, 0), 1, facecolor='#FF4444',
                   edgecolor='#8B0000', linewidth=2, label='Infected'),
        plt.Circle((0, 0), 1, facecolor='#FFFFFF',
                   edgecolor='#CCCCCC', linewidth=2, label='Not Infected'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FF4444', alpha=0.5,
                      edgecolor='gray', label='Fully Infected Household'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FFA500', alpha=0.4,
                      edgecolor='gray', label='Partially Infected Household'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FFFFFF', alpha=0.2,
                      edgecolor='gray', label='Not Infected Household')
    ]

    # Place legend outside the plots
    fig.legend(handles=legend_elements, loc='center',
               bbox_to_anchor=(0.5, 0.02), ncol=3, fontsize=10)

    # Add overall title
    fig.suptitle('Household Infection Progression Over Time\n(Best Run from Multi-Run Simulation)',
                 fontsize=16, fontweight='bold', y=0.95)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, top=0.85)

    if output_path:
        plt.savefig(output_path, dpi=300,
                    bbox_inches='tight', facecolor='white')
        print(f"Temporal infection heatmap saved to {output_path}")

    # Print summary statistics
    print(f"\nInfection Progression Summary:")
    for i, (timestep, count) in enumerate(zip(time_points, infection_counts)):
        if timestep == 0:
            time_label = "Initialization"
        else:
            day = timestep // 24
            time_label = f"Day {day}"
        print(
            f"  {time_label}: {count}/{len(person_pos)} infected ({count/len(person_pos)*100:.1f}%)")

    return fig


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Create temporal household infection heatmap')
    parser.add_argument('--data-dir', help='Directory containing data files')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--time-points', nargs='+', type=int, default=[0, 120, 240],
                        help='Time points in hours (default: 0, 120, 240 for initialization, day 5, day 10)')

    args = parser.parse_args()

    # Debug paths (comment out for production)
    # args.data_dir = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250825_231431_work_meeting_many_panvadere"
    # args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/temporal_infection_heatmap.png"

    if not os.path.exists(args.data_dir):
        print(f"Error: Directory {args.data_dir} does not exist")
        sys.exit(1)

    if not args.output_path:
        args.output_path = os.path.join(
            args.data_dir, 'temporal_infection_heatmap.png')

    # Create visualization
    fig = create_temporal_infection_heatmap(
        data_dir=args.data_dir,
        output_path=args.output_path,
        time_points=args.time_points
    )

    plt.show()


if __name__ == "__main__":
    main()
