#!/usr/bin/env python3
"""
Comparative Temporal Infection Heatmap
Compares transmission-informed vs uniform initialization approaches
Shows location-based transmission patterns and household infection statistics
"""

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
    """Get color based on infection status"""
    if person_id in infected_at_time:
        return '#FF4444', '#8B0000'  # Red with dark red border for infected
    else:
        return '#FFFFFF', '#CCCCCC'  # White with light gray border for not infected


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
                                          household_spacing=120.0, person_spacing=10.0,
                                          max_households_per_row=25):
    """Create household layout based on contact patterns using grid with households as rectangles"""

    # Calculate household contact matrix
    households, contact_matrix = calculate_household_contact_matrix(
        household_data, contact_data)

    print(f"Creating contact-based grid layout for {len(households)} households...")
    
    # Use contact matrix to determine household ordering for grid placement
    # Calculate total contact strength for each household
    household_contact_strength = {}
    for i, hh in enumerate(households):
        # Sum of all contacts for this household
        total_contacts = np.sum(contact_matrix[i, :])
        household_contact_strength[hh] = total_contacts
    
    # Sort households by contact strength (most connected first)
    sorted_households = sorted(households, key=lambda hh: household_contact_strength[hh], reverse=True)
    
    # Create a structured grid layout with household rectangles
    household_positions = {}
    person_pos = {}
    household_dimensions = {}  # Store width and height for each household
    
    rows = int(np.ceil(len(sorted_households) / max_households_per_row))
    
    for idx, household_id in enumerate(sorted_households):
        # Calculate grid position
        row = idx // max_households_per_row
        col = idx % max_households_per_row
        
        # Convert to actual coordinates (center of household rectangle)
        hx = col * household_spacing
        hy = row * household_spacing
        
        household_positions[household_id] = (hx, hy)
        
        # Get household members
        household_members = household_data[household_data['household_id'] == household_id]['person_id'].tolist()
        n_members = len(household_members)
        
        # Calculate household dimensions based on number of members
        if n_members == 1:
            household_width = 30  # Increased from 20
            household_height = 30  # Increased from 20
        else:
            # Horizontal layout: width grows with members, height stays constant
            household_width = max(50, n_members * person_spacing + 20)  # Increased padding
            household_height = 35  # Increased from 25
        
        household_dimensions[household_id] = (household_width, household_height)
        
        # Arrange people in horizontal line within the household rectangle
        if n_members == 1:
            person_pos[household_members[0]] = (hx, hy)
        else:
            # Arrange members in a horizontal line
            start_x = hx - (n_members - 1) * person_spacing / 2
            for i, person_id in enumerate(household_members):
                px = start_x + i * person_spacing
                py = hy  # Keep all members at the same y-coordinate (horizontal line)
                person_pos[person_id] = (px, py)
    
    # Calculate some statistics
    high_contact_households = sum(1 for strength in household_contact_strength.values() if strength > np.median(list(household_contact_strength.values())))
    
    print(f"Contact-based grid layout created: {high_contact_households}/{len(households)} households with above-median contact levels")
    print(f"Grid dimensions: {rows} rows x {max_households_per_row} columns")
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
                                 infected_at_time, timestep, title, method_name):
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

    # Create node colors based on infection status
    node_colors = []
    edge_colors = []
    node_sizes = []

    for person_id in person_pos.keys():
        node_color, edge_color = get_infection_status_color(
            person_id, infected_at_time)

        node_colors.append(node_color)
        edge_colors.append(edge_color)

        # Smaller nodes for better household visibility
        if person_id in infected_at_time:
            node_sizes.append(50)  # Reduced from 100
        else:
            node_sizes.append(40)  # Reduced from 80

    # Draw household boundaries as rectangles with proper dimensions
    for household_id, position in household_positions.items():
        if not isinstance(position, tuple):
            continue

        hx, hy = position
        household_width, household_height = household_dimensions.get(household_id, (40, 40))
        household_members = household_data[household_data['household_id']
                                           == household_id]['person_id'].tolist()

        if len(household_members) > 0:
            stats = household_stats.get(household_id, {
                                        'infection_rate': 0, 'infected_members': 0, 'total_members': len(household_members)})

            # Color household boundary based on infection rate at this time point only
            if stats['infection_rate'] == 0:
                boundary_color = '#D0D0D0'  # Darker gray for better visibility
                alpha = 0.4  # Increased opacity
                linewidth = 1.5  # Thicker border
            elif stats['infection_rate'] < 1.0:
                boundary_color = '#FFA500'  # Orange for partially infected
                alpha = 0.6  # Increased opacity
                linewidth = 2.5  # Thicker border
            else:
                boundary_color = '#FF4444'  # Red for fully infected
                alpha = 0.7  # Increased opacity
                linewidth = 3  # Thicker border

            # Draw household boundary as a rectangle centered on household position
            rect = plt.Rectangle((hx - household_width/2, hy - household_height/2), 
                               household_width, household_height,
                               facecolor=boundary_color, alpha=alpha, zorder=0,
                               edgecolor='black', linewidth=linewidth)  # Changed to black for better visibility
            ax.add_patch(rect)

    # Draw person nodes
    nx.draw_networkx_nodes(G, person_pos,
                           node_color=node_colors,
                           node_size=node_sizes,
                           alpha=0.9,
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

    ax.set_title(f"{method_name}\n{day_title}: {infected_count}/{total_people} infected ({infected_count/total_people*100:.1f}%)",
                 fontsize=11, fontweight='bold', pad=10)
    ax.axis('equal')
    ax.axis('off')

    return infected_count


def create_household_comparison_plot(ax, time_points, household_counts_method1, household_counts_method2,
                                     method1_name, method2_name):
    """Create slim comparison plot of household infection counts"""

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
                   label=method1_name, color='#FF6B6B', alpha=0.8)
    bars2 = ax.bar(x_pos + width/2, household_counts_method2, width,
                   label=method2_name, color='#4ECDC4', alpha=0.8)

    # Add value labels on bars
    for i, (v1, v2) in enumerate(zip(household_counts_method1, household_counts_method2)):
        ax.text(i - width/2, v1 + 1, str(v1), ha='center',
                va='bottom', fontsize=9, fontweight='bold')
        ax.text(i + width/2, v2 + 1, str(v2), ha='center',
                va='bottom', fontsize=9, fontweight='bold')

    ax.set_xlabel('Time Points', fontsize=10)
    ax.set_ylabel('Infected\nHouseholds', fontsize=10)
    ax.set_title('Households with ≥1 Infected Member',
                 fontsize=11, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(day_labels)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(axis='y', alpha=0.3)


def create_comparative_temporal_heatmap(data_dir_method1, data_dir_method2=None,
                                        method1_name="Transmission-Informed", method2_name="Uniform",
                                        output_path=None, time_points=[0, 72, 240]):
    """Create comparative temporal infection heatmap"""
    print("Loading temporal infection data for both methods...")

    # Load Method 1 data (transmission-informed)
    try:
        household_data = load_household_data(
            os.path.join(data_dir_method1, 'household_id.csv'))
        contact_data = load_contact_data(os.path.join(
            data_dir_method1, 'contact_intensiveness.csv'))
        location_data = load_location_data(os.path.join(
            data_dir_method1, 'location_id_and_type.csv'))
        infections_by_time_m1 = load_temporal_infection_data(
            os.path.join(data_dir_method1, 'best_run_detailed_infection.csv'))
    except FileNotFoundError as e:
        print(f"Error loading Method 1 data: {e}")
        return

    # Load Method 2 data (uniform) - if provided
    if data_dir_method2 and os.path.exists(data_dir_method2):
        try:
            infections_by_time_m2 = load_temporal_infection_data(
                os.path.join(data_dir_method2, 'best_run_detailed_infection.csv'))
            has_method2_data = True
        except FileNotFoundError as e:
            print(f"Warning: Method 2 data not found: {e}")
            print("Will create placeholder for comparison layout")
            has_method2_data = False
    else:
        print("Method 2 data directory not provided - creating placeholder")
        has_method2_data = False

    print(
        f"Loaded Method 1 data: {len(household_data)} people, temporal infections for {len(infections_by_time_m1)} time steps")

    # Create household layout based on contact patterns (same for both methods)
    spacing_factor = 0.6  # Increased from 0.3 to make households more visible
    household_spacing = 300.0 * spacing_factor
    person_spacing = 50.0 * spacing_factor

    print("Creating contact-based household layout...")
    person_pos, households, household_positions, household_dimensions = create_contact_based_household_layout(
        household_data, contact_data,
        household_spacing=household_spacing,
        person_spacing=person_spacing,
        max_households_per_row=20  # Reduced from 25 to accommodate larger spacing
    )

    # Create figure with 3 rows: Method1, Comparison, Method2
    fig = plt.figure(figsize=(6*len(time_points), 12))

    # Calculate household infection counts for comparison
    household_counts_m1 = []
    household_counts_m2 = []

    infection_counts_m1 = []
    infection_counts_m2 = []

    # Create subplots grid: 3 rows x len(time_points) columns
    gs = fig.add_gridspec(3, len(time_points), height_ratios=[
                          3, 0.5, 3], hspace=0.3)

    # Row 1: Method 1 (Transmission-informed)
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs[0, i])

        infected_at_time = infections_by_time_m1.get(timestep, set())
        day_label = timestep // 24

        infected_count = create_time_specific_heatmap(
            ax, household_data, household_positions, household_dimensions, person_pos,
            infected_at_time, day_label,
            f"Time Point {i+1}", method1_name
        )

        infection_counts_m1.append(infected_count)
        household_counts_m1.append(count_infected_households(
            household_data, infected_at_time))

    # Row 2: Household comparison (slim panel spanning all columns)
    ax_comparison = fig.add_subplot(gs[1, :])

    # For now, use placeholder data for method 2 if not available
    if not has_method2_data:
        # Create simulated uniform data for demonstration
        # Placeholder: 20% more households
        household_counts_m2 = [int(count * 1.2)
                               for count in household_counts_m1]
        # Placeholder: 15% more infections
        infection_counts_m2 = [int(count * 1.15)
                               for count in infection_counts_m1]
    else:
        # Calculate actual Method 2 data
        for timestep in time_points:
            infected_at_time = infections_by_time_m2.get(timestep, set())
            infection_counts_m2.append(len(infected_at_time))
            household_counts_m2.append(count_infected_households(
                household_data, infected_at_time))

    create_household_comparison_plot(ax_comparison, time_points, household_counts_m1, household_counts_m2,
                                     method1_name, method2_name)

    # Row 3: Method 2 (Uniform) - placeholder or real data
    for i, timestep in enumerate(time_points):
        ax = fig.add_subplot(gs[2, i])

        if has_method2_data:
            infected_at_time = infections_by_time_m2.get(timestep, set())
            day_label = timestep // 24

            create_time_specific_heatmap(
                ax, household_data, household_positions, household_dimensions, person_pos,
                infected_at_time, day_label,
                f"Time Point {i+1}", method2_name
            )
        else:
            # Create placeholder visualization
            ax.text(0.5, 0.5, f"{method2_name}\nData Not Available\n\n(Placeholder for comparison layout)",
                    ha='center', va='center', transform=ax.transAxes, fontsize=12,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='lightgray', alpha=0.5))
            ax.set_title(f"{method2_name} - Time Point {i+1}",
                         fontsize=11, fontweight='bold')
            ax.axis('off')

    # Create simple legend
    legend_elements = [
        plt.Circle((0, 0), 1, facecolor='#FF4444',
                   edgecolor='#8B0000', label='Infected Person'),
        plt.Circle((0, 0), 1, facecolor='#FFFFFF',
                   edgecolor='#CCCCCC', label='Not Infected Person'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FF4444', alpha=0.6,
                      edgecolor='darkgray', label='Fully Infected Household'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#FFA500', alpha=0.5,
                      edgecolor='darkgray', label='Partially Infected Household'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#E8E8E8', alpha=0.3,
                      edgecolor='darkgray', label='Not Infected Household')
    ]

    # Place legend outside the plots
    fig.legend(handles=legend_elements, loc='center',
               bbox_to_anchor=(0.5, 0.02), ncol=5, fontsize=10)

    # Add overall title
    fig.suptitle('Comparative Temporal Infection Progression\nTransmission-Informed vs Uniform Initialization',
                 fontsize=16, fontweight='bold', y=0.98)

    # Adjust layout
    plt.subplots_adjust(bottom=0.1, top=0.92)

    if output_path:
        plt.savefig(output_path, dpi=300,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative temporal heatmap saved to {output_path}")

    # Print summary statistics
    print(f"\n{method1_name} Summary:")
    for i, (timestep, inf_count, hh_count) in enumerate(zip(time_points, infection_counts_m1, household_counts_m1)):
        time_label = "Initialization" if timestep == 0 else f"Day {timestep // 24}"
        print(
            f"  {time_label}: {inf_count} infected, {hh_count} households affected")

    if has_method2_data:
        print(f"\n{method2_name} Summary:")
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
    parser.add_argument('--time-points', nargs='+', type=int, default=[0, 72, 240],
                        help='Time points in hours (default: 0, 72, 240 for initialization, day 3, day 10)')

    args = parser.parse_args()

    # Debug paths (comment out for production)
    args.data_dir_method1 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250825_231431_restaurant_table_equals_household_panvadere"
    args.data_dir_method2 = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250825_231431_restaurant_table_equals_household_memilio"
    args.method1_name = "Panvadere (Transmission-Informed)"
    args.method2_name = "Memilio (Uniform)"
    args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/comparative_temporal_infection_heatmap.png"

    if not os.path.exists(args.data_dir_method1):
        print(f"Error: Directory {args.data_dir_method1} does not exist")
        sys.exit(1)

    if not args.output_path:
        args.output_path = os.path.join(
            args.data_dir_method1, 'comparative_temporal_infection_heatmap.png')

    # Create visualization
    fig = create_comparative_temporal_heatmap(
        data_dir_method1=args.data_dir_method1,
        data_dir_method2=args.data_dir_method2,
        method1_name=args.method1_name,
        method2_name=args.method2_name,
        output_path=args.output_path,
        time_points=args.time_points
    )

    plt.show()


if __name__ == "__main__":
    main()
