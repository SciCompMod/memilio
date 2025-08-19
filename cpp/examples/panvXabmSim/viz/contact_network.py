import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
import numpy as np
from collections import defaultdict
import seaborn as sns
from matplotlib.patches import Circle
import matplotlib.patches as mpatches
import argparse
import os
import sys


def load_household_data(household_file):
    """Load household assignment data from CSV file"""
    df = pd.read_csv(household_file)
    # Rename columns to match expected format
    df.columns = ['person_id', 'household_id']
    return df


def load_infection_count_data(infection_file):
    """Load infection count data from CSV file"""
    df = pd.read_csv(infection_file)
    # Rename columns to match expected format
    df.columns = ['person_id', 'infection_count']
    return df


def load_contact_data(contact_file):
    """Load contact intensiveness data from CSV file"""
    df = pd.read_csv(contact_file)
    # Rename columns to match expected format
    df.columns = ['person_1', 'person_2', 'location_id', 'contact_hours']
    return df


def load_location_data(location_file):
    """Load location type data from CSV file"""
    df = pd.read_csv(location_file)
    # Rename columns to match expected format
    df.columns = ['location_id', 'location_type']
    return df


def calculate_infection_rates(infection_count_data, total_runs=5):
    """
    Calculate infection rates for each person based on infection counts

    Parameters:
    infection_count_data: DataFrame with columns ['person_id', 'infection_count']
    total_runs: total number of simulation runs
    """
    infection_rates = {}

    for _, row in infection_count_data.iterrows():
        person_id = row['person_id']
        infection_count = row['infection_count']
        infection_rate = (infection_count / total_runs) * 100
        infection_rates[person_id] = infection_rate

    return infection_rates


def prepare_contact_data_with_locations(contact_data, location_data):
    """
    Merge contact data with location types

    Parameters:
    contact_data: DataFrame with contact intensiveness
    location_data: DataFrame with location types
    """
    # Merge contact data with location types
    merged_data = contact_data.merge(
        location_data, on='location_id', how='left')

    # Handle missing location types
    merged_data['location_type'] = merged_data['location_type'].fillna(
        'Unknown')

    return merged_data


def calculate_contact_strength_from_hours(contact_data_with_types):
    """
    Calculate contact strength between person pairs based on contact hours
    Returns dict with (person1, person2, location_type) -> total_contact_hours
    """
    contact_strength = defaultdict(float)
    contact_locations = defaultdict(set)

    for _, row in contact_data_with_types.iterrows():
        person1, person2 = sorted([row['person_1'], row['person_2']])
        location_type = row['location_type']
        location_id = row['location_id']
        contact_hours = row['contact_hours']

        # Sum contact hours by location type
        key = (person1, person2, location_type)
        contact_strength[key] += contact_hours
        contact_locations[key].add(location_id)

    return contact_strength, contact_locations


def create_household_layout(household_data, household_spacing=15.0, person_spacing=2.0, max_households_per_row=15):
    """
    Create compact household layout with multiple rows
    """
    # Group persons by household
    households = household_data.groupby(
        'household_id')['person_id'].apply(list).to_dict()

    pos = {}
    household_positions = {}

    num_households = len(households)

    # Calculate number of rows needed
    households_per_row = min(max_households_per_row, num_households)
    num_rows = int(np.ceil(num_households / households_per_row))
    row_spacing = 12.0

    household_items = list(households.items())

    for i, (household_id, members) in enumerate(household_items):
        row = i // households_per_row
        col = i % households_per_row

        # Calculate position for this row
        actual_households_in_row = min(
            households_per_row, num_households - row * households_per_row)
        total_width_this_row = (
            actual_households_in_row - 1) * household_spacing
        start_x = -total_width_this_row / 2

        household_center_x = start_x + col * household_spacing
        household_center_y = row * row_spacing - \
            ((num_rows - 1) * row_spacing / 2)

        household_positions[household_id] = (
            household_center_x, household_center_y)

        # Arrange household members
        n_members = len(members)
        if n_members == 1:
            pos[members[0]] = (household_center_x, household_center_y)
        elif n_members <= 6:
            member_width = (n_members - 1) * person_spacing
            member_start_x = household_center_x - member_width / 2

            for j, person_id in enumerate(members):
                x_pos = member_start_x + j * person_spacing
                pos[person_id] = (x_pos, household_center_y)
        # else:
        #     # For large households, show representative members only
        #     visible_members = members[:3] + members[-2:]
        #     for j, person_id in enumerate(visible_members):
        #         if j < 3:
        #             x_offset = (j - 1) * person_spacing
        #         else:
        #             x_offset = (j - 1.5) * person_spacing
        #         pos[person_id] = (household_center_x +
        #                           x_offset, household_center_y)

            household_positions[f"{household_id}_hidden"] = len(members) - 5

    return pos, households, household_positions


def create_location_layout(contact_data_with_types, household_positions, area_distance=80.0):
    """
    Create layout for location nodes
    """
    location_types = contact_data_with_types['location_type'].unique()

    # Get the range of household positions
    if household_positions:
        household_x_positions = [
            pos[0] for pos in household_positions.values() if isinstance(pos, tuple)]
        household_y_positions = [
            pos[1] for pos in household_positions.values() if isinstance(pos, tuple)]

        min_x = min(household_x_positions) if household_x_positions else 0
        max_x = max(household_x_positions) if household_x_positions else 0
        min_y = min(household_y_positions) if household_y_positions else 0
        max_y = max(household_y_positions) if household_y_positions else 0

        center_x = (min_x + max_x) / 2
        width = max_x - min_x
        household_area_height = max_y - min_y + 10
    else:
        center_x = 0
        width = 100
        household_area_height = 5

    # Define location areas around the household area
    location_area_mapping = {
        'Work': (center_x - width/2, -(area_distance + household_area_height/4)),
        'School': (center_x + width/2, -(area_distance + household_area_height/4)),
        'SocialEvent': (center_x - width/2, area_distance + household_area_height / 4),
        'BasicsShop': (center_x + width/2, area_distance + household_area_height/4),
        'Home': (center_x, area_distance + household_area_height/2 + 20),
        'Unknown': (center_x, area_distance + household_area_height/2 + 20)
    }

    location_positions = {}
    location_nodes = set()
    location_areas = {}

    for location_type in location_types:
        # Get area position for this type
        if location_type in location_area_mapping:
            area_pos = location_area_mapping[location_type]
        else:
            # Place unknown types in a default area
            print(f"Unknown location type: {location_type}")
            area_pos = location_area_mapping['Unknown']

        location_areas[location_type] = area_pos

        # Get unique location IDs for this type
        type_locations = contact_data_with_types[
            contact_data_with_types['location_type'] == location_type
        ]['location_id'].unique()

        area_x, area_y = area_pos

        # Arrange locations with spacing in multiple rows
        n_locations = len(type_locations)
        if n_locations == 1:
            location_node = f"{location_type}_{type_locations[0]}"
            location_positions[location_node] = (area_x, area_y)
            location_nodes.add(location_node)
        else:
            max_locations_per_row = 16  # Maximum locations per row
            spacing = 6.0
            row_spacing = 8.0

            # Calculate number of rows needed
            locations_per_row = min(max_locations_per_row, n_locations)
            num_rows = int(np.ceil(n_locations / locations_per_row))

            for i, location_id in enumerate(type_locations):
                row = i // locations_per_row
                col = i % locations_per_row

                # Calculate position for this row
                actual_locations_in_row = min(
                    locations_per_row, n_locations - row * locations_per_row)
                total_width_this_row = (actual_locations_in_row - 1) * spacing
                start_x = area_x - total_width_this_row / 2

                x_pos = start_x + col * spacing
                y_pos = area_y + row * row_spacing - \
                    ((num_rows - 1) * row_spacing / 2)

                location_node = f"{location_type}_{location_id}"
                location_positions[location_node] = (x_pos, y_pos)
                location_nodes.add(location_node)

    return location_positions, location_nodes, location_areas


def get_location_edge_colors():
    """Define colors for different location types"""
    return {
        'Home': '#2E8B57',
        'Work': '#4169E1',
        'School': '#FF6347',
        'SocialEvent': '#9370DB',
        'BasicsShop': '#FF8C00',
        'Unknown': '#808080'
    }


def get_infection_color(infection_rate):
    """Color mapping for infection rates: light yellow → deep red"""
    colors = ['#FFF8DC', '#FFE4B5', '#FFA500', '#FF6347', '#B22222']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('infection', colors, N=n_bins)
    return cmap(infection_rate / 100)


def get_infection_outline(infection_rate):
    """ Outer ring of person is either white for an infection rate below 50% or red for above """
    if infection_rate < 10:
        return '#FFFFFF'  # White
    else:
        return '#FF0000'  # Red


def get_edge_width_from_hours(contact_hours, max_hours=None):
    """Convert contact hours to edge width"""
    if max_hours is None:
        max_hours = 100  # Default max hours

    # Normalize to 0.1-1.0 range for edge width
    min_width = 0.1
    max_width = 1.0

    normalized = min(contact_hours / max_hours, 1.0)
    width = min_width + (max_width - min_width) * normalized

    return width


def visualize_contact_network(household_data, infection_count_data, contact_data, location_data,
                              total_runs=5, figsize=(24, 18), save_path=None,
                              min_contact_hours=5, max_persons=200, min_infection_rate=0,
                              viz_home_contacts=False):
    """
    Create the contact network visualization based on your simulation data

    Parameters:
    household_data: DataFrame with person_id and household_id
    infection_count_data: DataFrame with person_id and infection_count
    contact_data: DataFrame with person_1, person_2, location_id, contact_hours
    location_data: DataFrame with location_id and location_type
    total_runs: number of simulation runs
    min_contact_hours: minimum contact hours to show edge
    max_persons: maximum persons to display (for performance)
    min_infection_rate: minimum infection rate (%) to display person (0-100)
    """

    print("Processing data...")

    # Calculate infection rates
    infection_rates = calculate_infection_rates(
        infection_count_data, total_runs)

    # Filter persons by minimum infection rate
    if min_infection_rate > 0:
        print(
            f"Filtering persons with infection rate >= {min_infection_rate}%...")
        selected_by_infection = [person_id for person_id, rate in infection_rates.items()
                                 if rate >= min_infection_rate]

        household_data = household_data[household_data['person_id'].isin(
            selected_by_infection)]
        infection_count_data = infection_count_data[infection_count_data['person_id'].isin(
            selected_by_infection)]

        # Recalculate infection rates for filtered data
        infection_rates = calculate_infection_rates(
            infection_count_data, total_runs)

        print(f"Persons after infection rate filter: {len(household_data)}")

    # Merge contact data with location types
    contact_data_with_types = prepare_contact_data_with_locations(
        contact_data, location_data)

    # Filter for performance if needed (after infection rate filtering)
    if len(household_data) > max_persons:
        print(
            f"Reducing dataset from {len(household_data)} to {max_persons} persons for visualization...")
        selected_persons = household_data['person_id'].sample(
            n=max_persons, random_state=42).tolist()

        household_data = household_data[household_data['person_id'].isin(
            selected_persons)]
        infection_count_data = infection_count_data[infection_count_data['person_id'].isin(
            selected_persons)]

        # Recalculate infection rates for filtered data
        infection_rates = calculate_infection_rates(
            infection_count_data, total_runs)

    # Filter contact data to only include remaining persons
    remaining_persons = household_data['person_id'].tolist()
    contact_data_with_types = contact_data_with_types[
        (contact_data_with_types['person_1'].isin(remaining_persons)) &
        (contact_data_with_types['person_2'].isin(remaining_persons))
    ]

    # Filter contacts by minimum hours
    contact_data_with_types = contact_data_with_types[
        contact_data_with_types['contact_hours'] >= min_contact_hours
    ]

    # Filter out Hospital and ICU for location types
    contact_data_with_types = contact_data_with_types[
        contact_data_with_types['location_type'].isin(
            ['Hospital', 'ICU']) == False
    ]

    # Filter out homes, if not visualizing home contacts
    if not viz_home_contacts:
        contact_data_with_types = contact_data_with_types[
            contact_data_with_types['location_type'] != 'Home'
        ]

    print(
        f"Visualizing {len(household_data)} persons with {len(contact_data_with_types)} significant contacts...")
    if min_infection_rate > 0:
        print(f"(Filtered for infection rate >= {min_infection_rate}%)")

    # Calculate contact strengths
    contact_strength, contact_locations = calculate_contact_strength_from_hours(
        contact_data_with_types)

    # Create layouts
    person_pos, households, household_positions = create_household_layout(
        household_data, household_spacing=12.0, person_spacing=1.8, max_households_per_row=15)

    location_pos, location_nodes, location_areas = create_location_layout(
        contact_data_with_types, household_positions, area_distance=60.0)

    # Combine all positions
    all_pos = {**person_pos, **location_pos}

    # Create the graph
    plt.figure(figsize=figsize)
    G = nx.Graph()

    # Add person nodes
    for person_id in person_pos.keys():
        G.add_node(person_id, node_type='person')

    # Add location nodes
    for location_node in location_nodes:
        G.add_node(location_node, node_type='location')

    # Draw household areas
    for household_id, household_pos in household_positions.items():
        if '_hidden' in str(household_id):
            continue

        household_center_x, household_center_y = household_pos
        members = households.get(household_id, [])

        if len(members) > 1:
            fancy_box = FancyBboxPatch(
                (household_center_x - 3, household_center_y - 1), 6, 2,
                boxstyle="round,pad=0.1",
                facecolor='lightgray',
                alpha=0.25,
                edgecolor='darkgray',
                linewidth=0.8,
                zorder=0
            )
            plt.gca().add_patch(fancy_box)

            if len(members) > 6:
                plt.text(household_center_x + 5, household_center_y, f"({len(members)})",
                         fontsize=8, alpha=0.7, fontweight='bold')

    # Draw location areas
    for location_type, (area_x, area_y) in location_areas.items():
        if any(location_type in node for node in location_nodes):
            # Count locations of this type to determine rectangle size
            type_locations = contact_data_with_types[
                contact_data_with_types['location_type'] == location_type
            ]['location_id'].unique()
            n_locations = len(type_locations)

            # Calculate dynamic rectangle size based on location layout
            max_locations_per_row = 16
            spacing = 6.0
            row_spacing = 8.0

            if n_locations <= max_locations_per_row:
                # Single row
                rect_width = max(16, n_locations * spacing + 4)
                rect_height = 16
            else:
                # Multiple rows
                num_rows = int(np.ceil(n_locations / max_locations_per_row))
                rect_width = max_locations_per_row * spacing + 4
                rect_height = max(16, num_rows * row_spacing + 4)

            rect = plt.Rectangle((area_x - rect_width/2, area_y - rect_height/2),
                                 rect_width, rect_height,
                                 facecolor='lightblue', alpha=0.1, zorder=0,
                                 edgecolor='blue', linewidth=1)
            plt.gca().add_patch(rect)

            label_y = area_y + rect_height/2 + 5 if area_y > 0 else area_y - rect_height/2 - 5
            plt.text(area_x, label_y, location_type.replace('_', ' '),
                     ha='center', va='center', fontsize=10, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Draw person nodes with infection rate colors
    person_colors = [get_infection_color(
        infection_rates.get(p, 0)) for p in person_pos.keys()]
    person_outlines = [get_infection_outline(
        infection_rates.get(p, 0)) for p in person_pos.keys()]

    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=person_colors,
                           node_size=120,
                           alpha=0.9,
                           edgecolors=person_outlines,
                           linewidths=1.5)

    # Draw location nodes
    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(location_nodes),
                           node_color='lightsteelblue',
                           node_size=80,
                           alpha=0.7,
                           node_shape='s',
                           edgecolors='navy',
                           linewidths=0.5)

    # Get edge colors and max contact hours for width scaling
    edge_colors = get_location_edge_colors()
    max_contact_hours = max(contact_data_with_types['contact_hours']) if len(
        contact_data_with_types) > 0 else 100

    # Draw edges based on contact intensity
    edges_by_type = defaultdict(list)
    edge_widths_by_type = defaultdict(list)

    # # Add person-to-location edges
    # for (person1, person2, location_type), total_hours in contact_strength.items():
    #     if total_hours >= min_contact_hours:
    #         locations = contact_locations[(person1, person2, location_type)]

    #         for location_id in locations:
    #             location_node = f"{location_type}_{location_id}"
    #             if location_node in location_nodes:
    #                 # Add edges from both persons to location
    #                 edge_width = get_edge_width_from_hours(
    #                     total_hours, max_contact_hours)

    #                 edges_by_type[location_type].append(
    #                     (person1, location_node))
    #                 edges_by_type[location_type].append(
    #                     (person2, location_node))
    #                 edge_widths_by_type[location_type].extend(
    #                     [edge_width, edge_width])

    #                 G.add_edge(person1, location_node,
    #                            weight=total_hours, opacity=edge_width)
    #                 G.add_edge(person2, location_node,
    #                            weight=total_hours, opacity=edge_width)

    # # Draw edges by type with appropriate colors and widths
    # for location_type, edges in edges_by_type.items():
    #     if edges:
    #         color = edge_colors.get(location_type, 'gray')
    #         widths = edge_widths_by_type[location_type]

    #         # Draw edges with varying widths
    #         for i, (person, location) in enumerate(edges):
    #             width = widths[i] if i < len(widths) else 0.5
    #             nx.draw_networkx_edges(G, all_pos,
    #                                    edgelist=[(person, location)],
    #                                    width=width,
    #                                    edge_color=color,
    #                                    alpha=0.6,
    #                                    connectionstyle="arc3,rad=0.1")

    # Draw location labels
    location_labels = {loc: loc.split('_')[-1] for loc in location_nodes}
    nx.draw_networkx_labels(G, all_pos,
                            labels=location_labels,
                            font_size=6)

    # Create infection rate colorbar
    sm = plt.cm.ScalarMappable(
        cmap='YlOrRd', norm=plt.Normalize(vmin=0, vmax=100))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca(), shrink=0.6, aspect=30)
    cbar.set_label('Infection Rate (%)', rotation=270, labelpad=15)

    # Create legend
    legend_elements = []

    # Infection rate samples
    for rate, label in [(0, 'Never infected'), (50, 'Sometimes infected'), (100, 'Always infected')]:
        color = get_infection_color(rate)
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          markeredgecolor='black', markeredgewidth=0.5,
                                          label=label, linestyle='None'))

    legend_elements.append(plt.Line2D([0], [0], color='none', label=''))

    # Location type edges
    for loc_type, color in edge_colors.items():
        if loc_type in [edge_type for edge_type in edges_by_type.keys()]:
            legend_elements.append(plt.Line2D([0], [0], color=color, lw=2,
                                              label=f'{loc_type} Contacts'))

    # Edge width explanation
    legend_elements.append(plt.Line2D([0], [0], color='none', label=''))
    legend_elements.append(plt.Line2D([0], [0], color='gray', lw=1,
                                      label=f'Thin: {min_contact_hours}h contact'))
    legend_elements.append(plt.Line2D([0], [0], color='gray', lw=3,
                                      label=f'Thick: {max_contact_hours:.0f}h contact'))

    plt.legend(handles=legend_elements, loc='center left',
               bbox_to_anchor=(0.02, 0.98),
               frameon=True, fancybox=True, shadow=True,
               fontsize=9, title='Legend', title_fontsize=10)

    # Add statistics
    num_persons = len(person_pos)
    num_households = len(
        [h for h in household_positions.keys() if '_hidden' not in str(h)])
    avg_infection_rate = np.mean(list(infection_rates.values()))
    total_contact_hours = contact_data_with_types['contact_hours'].sum()

    title_text = f'Contact Network Analysis\n' + \
        f'{num_persons} persons, {num_households} households, {total_runs} simulation runs\n' + \
        f'Avg infection rate: {avg_infection_rate:.1f}%, Total contact hours: {total_contact_hours:,}'

    if min_infection_rate > 0:
        title_text += f'\n(Showing only persons with ≥{min_infection_rate}% infection rate)'

    # Amount of persons with over 50% infection rate:
    num_high_infection = sum(
        1 for rate in infection_rates.values() if rate > 99)
    title_text += f'\n{num_high_infection} persons with >50% infection rate'

    plt.title(title_text, fontsize=14, pad=20)

    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Contact network visualization saved to {save_path}")

    # Print summary statistics
    print("\n=== Network Summary ===")
    print(f"Total persons: {num_persons}")
    print(f"Total households: {num_households}")
    print(f"Average infection rate: {avg_infection_rate:.2f}%")
    print(f"Total contact hours: {total_contact_hours:,}")
    print(
        f"Contact types: {', '.join(contact_data_with_types['location_type'].unique())}")
    if min_infection_rate > 0:
        print(f"Minimum infection rate filter: {min_infection_rate}%")


def main():
    """Main function to run the analysis"""
    parser = argparse.ArgumentParser(
        description='Generate contact network visualization from simulation data')
    parser.add_argument('--data-dir',
                        help='Directory containing CSV files')
    parser.add_argument(
        '--output-path', help='Output path for visualization (optional)')
    parser.add_argument('--scenario-name', default='Simulation',
                        help='Scenario name for titles')
    parser.add_argument('--total-runs', type=int, default=5,
                        help='Total number of simulation runs')
    parser.add_argument('--min-contact-hours', type=int,
                        default=10, help='Minimum contact hours to show')
    parser.add_argument('--max-persons', type=int,
                        default=200, help='Maximum persons to display')
    parser.add_argument('--min-infection-rate', type=float,
                        default=0, help='Minimum infection rate (%) to display')
    parser.add_argument('--viz_home_contacts', action='store_true',
                        help='Visualize home contacts')

    args = parser.parse_args()

    # # For debug, set a path:
    args.data_dir = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/last_result"
    args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/last_result/contact_network.png"
    args.scenario_name = "Test Scenario"
    args.total_runs = 10

    # Check if directory exists
    if not os.path.exists(args.data_dir):
        print(f"Error: Directory {args.data_dir} does not exist")
        sys.exit(1)

    # Load data files
    print("Loading data...")
    try:
        household_data = load_household_data(
            os.path.join(args.data_dir, 'household_id.csv'))
        infection_count_data = load_infection_count_data(
            os.path.join(args.data_dir, 'infection_count.csv'))
        contact_data = load_contact_data(os.path.join(
            args.data_dir, 'contact_intensiveness.csv'))
        location_data = load_location_data(os.path.join(
            args.data_dir, 'location_id_and_type.csv'))
    except FileNotFoundError as e:
        print(f"Error: Required CSV file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

    print(f"Loaded {len(household_data)} household assignments")
    print(f"Loaded {len(infection_count_data)} infection records")
    print(f"Loaded {len(contact_data)} contact records")
    print(f"Loaded {len(location_data)} location types")

    # Set output path
    if args.output_path:
        save_path = args.output_path
    else:
        save_path = os.path.join(
            args.data_dir, f'contact_network_{args.scenario_name.lower()}.png')

    # Create visualization
    visualize_contact_network(
        household_data=household_data,
        infection_count_data=infection_count_data,
        contact_data=contact_data,
        location_data=location_data,
        total_runs=args.total_runs,
        figsize=(24, 18),
        save_path=save_path,
        min_contact_hours=args.min_contact_hours,
        max_persons=args.max_persons,
        min_infection_rate=args.min_infection_rate,
        viz_home_contacts=args.viz_home_contacts
    )


if __name__ == "__main__":
    main()
