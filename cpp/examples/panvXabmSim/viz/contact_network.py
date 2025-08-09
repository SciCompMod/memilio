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


def analyze_infection_clustering(household_data, infection_rates):
    """Analyze and visualize infection clustering patterns"""

    # Calculate household infection rates
    household_infections = {}
    for _, row in household_data.iterrows():
        hh_id = row['household_id']
        person_rate = infection_rates.get(row['person_id'], 0)

        if hh_id not in household_infections:
            household_infections[hh_id] = []
        household_infections[hh_id].append(person_rate)

    # Calculate clustering metrics
    household_avg_rates = {hh: np.mean(rates)
                           for hh, rates in household_infections.items()}
    household_variance = {hh: np.var(rates)
                          for hh, rates in household_infections.items()}

    # High-infection households (>70% average)
    high_infection_households = [
        hh for hh, avg in household_avg_rates.items() if avg > 70]

    return household_avg_rates, household_variance, high_infection_households


def create_comparison_visualization(contact_data, household_data,
                                    infection_rates_sim1, infection_rates_sim2,
                                    sim1_name="Widespread", sim2_name="Clustered",
                                    save_path=None):
    """Create side-by-side comparison of two simulations with clustering highlights"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(32, 16))

    # Simulation 1 (left)
    plt.sca(ax1)

    # Create layouts (reuse from main function)
    contact_strength, contact_locations = calculate_contact_strength(
        contact_data)
    person_pos, households, household_positions = create_household_layout(
        household_data, household_spacing=12.0, person_spacing=1.8, max_households_per_row=15)
    location_pos, location_nodes, location_areas = create_location_layout(
        contact_data, household_positions, area_distance=60.0)
    all_pos = {**person_pos, **location_pos}

    # Analyze clustering for sim1
    _, _, high_infection_hh_sim1 = analyze_infection_clustering(
        household_data, infection_rates_sim1)

    # Highlight infection clusters for sim1
    for hh_id in high_infection_hh_sim1:
        hh_members = household_data[household_data['household_id']
                                    == hh_id]['person_id'].tolist()
        hh_positions = [all_pos[person]
                        for person in hh_members if person in all_pos]

        if len(hh_positions) >= 2:
            center_x = np.mean([pos[0] for pos in hh_positions])
            center_y = np.mean([pos[1] for pos in hh_positions])

            circle = Circle((center_x, center_y), radius=10,
                            facecolor='red', alpha=0.15, zorder=-1)
            ax1.add_patch(circle)

            ax1.text(center_x, center_y + 12, "CLUSTER",
                     ha='center', va='center', fontsize=8, fontweight='bold',
                     color='red', alpha=0.8)

    # Draw network for sim1
    G1 = nx.Graph()
    for person_id in person_pos.keys():
        G1.add_node(person_id)

    # Draw nodes
    nx.draw_networkx_nodes(G1, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=[get_infection_color_v2(
                               infection_rates_sim1.get(p, 0)) for p in person_pos.keys()],
                           node_size=120,
                           alpha=0.9,
                           edgecolors='white',
                           linewidths=1.5)

    # Add basic edges and areas (simplified)
    draw_basic_network_elements(
        ax1, all_pos, location_nodes, location_areas, household_positions, households)

    ax1.set_title(f'{sim1_name} Infection Pattern\n({len(high_infection_hh_sim1)} clustered households)',
                  fontsize=14, fontweight='bold')
    ax1.axis('equal')
    ax1.axis('off')

    # Simulation 2 (right)
    plt.sca(ax2)

    # Analyze clustering for sim2
    _, _, high_infection_hh_sim2 = analyze_infection_clustering(
        household_data, infection_rates_sim2)

    # Highlight infection clusters for sim2
    for hh_id in high_infection_hh_sim2:
        hh_members = household_data[household_data['household_id']
                                    == hh_id]['person_id'].tolist()
        hh_positions = [all_pos[person]
                        for person in hh_members if person in all_pos]

        if len(hh_positions) >= 2:
            center_x = np.mean([pos[0] for pos in hh_positions])
            center_y = np.mean([pos[1] for pos in hh_positions])

            circle = Circle((center_x, center_y), radius=10,
                            facecolor='red', alpha=0.15, zorder=-1)
            ax2.add_patch(circle)

            ax2.text(center_x, center_y + 12, "CLUSTER",
                     ha='center', va='center', fontsize=8, fontweight='bold',
                     color='red', alpha=0.8)

    # Draw network for sim2
    G2 = nx.Graph()
    for person_id in person_pos.keys():
        G2.add_node(person_id)

    # Draw nodes
    nx.draw_networkx_nodes(G2, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=[get_infection_color_v2(
                               infection_rates_sim2.get(p, 0)) for p in person_pos.keys()],
                           node_size=120,
                           alpha=0.9,
                           edgecolors='white',
                           linewidths=1.5)

    # Add basic edges and areas
    draw_basic_network_elements(
        ax2, all_pos, location_nodes, location_areas, household_positions, households)

    ax2.set_title(f'{sim2_name} Infection Pattern\n({len(high_infection_hh_sim2)} clustered households)',
                  fontsize=14, fontweight='bold')
    ax2.axis('equal')
    ax2.axis('off')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def draw_basic_network_elements(ax, all_pos, location_nodes, location_areas, household_positions, households):
    """Helper function to draw basic network elements"""

    # Draw location areas
    for location_type, (area_x, area_y) in location_areas.items():
        if any(location_type in node for node in location_nodes):
            rect = plt.Rectangle((area_x - 8, area_y - 8), 16, 16,
                                 facecolor='lightblue', alpha=0.1, zorder=0,
                                 edgecolor='blue', linewidth=1)
            ax.add_patch(rect)

            if area_y > 0:
                label_y = area_y + 10
            else:
                label_y = area_y - 10

            ax.text(area_x, label_y, location_type.replace('_', ' '),
                    ha='center', va='center', fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Draw household areas
    for household_id, household_pos in household_positions.items():
        if '_hidden' in household_id:
            continue

        household_center_x, household_center_y = household_pos
        members = households[household_id]

        if len(members) > 1:
            fancy_box = FancyBboxPatch(
                (household_center_x - 6, household_center_y - 1), 12, 2,
                boxstyle="round,pad=0.1",
                facecolor='lightgray',
                alpha=0.25,
                edgecolor='darkgray',
                linewidth=0.8,
                zorder=0
            )
            ax.add_patch(fancy_box)


def add_statistics_panel(infection_rates, household_data, position='right', panel_title="SIMULATION"):
    """Add statistical comparison panel to current plot"""

    # Calculate key metrics
    total_infected = sum(1 for rate in infection_rates.values() if rate > 0)
    total_persons = len(infection_rates)
    avg_rate = np.mean(list(infection_rates.values()))

    # Clustering metrics
    household_avg_rates, household_variance, high_infection_hh = analyze_infection_clustering(
        household_data, infection_rates)
    clustering_index = len(
        high_infection_hh) / len(household_avg_rates) * 100 if household_avg_rates else 0

    # Calculate infection spread (standard deviation across households)
    household_spread = np.std(
        list(household_avg_rates.values())) if household_avg_rates else 0

    # Determine pattern type
    if clustering_index > 30 and household_spread > 25:
        pattern_type = "HIGHLY CLUSTERED"
        pattern_color = 'red'
    elif clustering_index > 15:
        pattern_type = "MODERATELY CLUSTERED"
        pattern_color = 'orange'
    else:
        pattern_type = "WIDESPREAD"
        pattern_color = 'green'

    # Create text box
    stats_text = f"""{panel_title} ANALYSIS

Total Infected: {total_infected}/{total_persons} ({total_infected/total_persons*100:.1f}%)
Average Rate: {avg_rate:.1f}%

CLUSTERING METRICS:
- High-infection households: {len(high_infection_hh)}
- Clustering index: {clustering_index:.1f}%
- Household spread: {household_spread:.1f}%
- Max household rate: {max(household_avg_rates.values()):.1f}%

PATTERN: {pattern_type}"""

    # Position the text box
    if position == 'right':
        x, y = 0.98, 0.98
        ha, va = 'right', 'top'
    elif position == 'left':
        x, y = 0.02, 0.98
        ha, va = 'left', 'top'
    elif position == 'bottom-right':
        x, y = 0.98, 0.02
        ha, va = 'right', 'bottom'
    else:  # bottom-left
        x, y = 0.02, 0.02
        ha, va = 'left', 'bottom'

    plt.text(x, y, stats_text, transform=plt.gca().transAxes,
             fontsize=9, ha=ha, va=va, fontfamily='monospace',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='white', alpha=0.95,
                       edgecolor=pattern_color, linewidth=2))


def create_comparison_with_statistics(contact_data, household_data,
                                      infection_rates_sim1, infection_rates_sim2,
                                      sim1_name="Widespread", sim2_name="Clustered",
                                      save_path=None):
    """Create comparison with statistical overlays"""

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(32, 16))

    # Plot simulation 1
    plt.sca(ax1)
    visualize_household_network(contact_data, household_data, infection_rates_sim1,
                                figsize=(16, 16), enable_reduction=False)
    add_statistics_panel(infection_rates_sim1, household_data,
                         position='left', panel_title=sim1_name.upper())

    # Plot simulation 2
    plt.sca(ax2)
    visualize_household_network(contact_data, household_data, infection_rates_sim2,
                                figsize=(16, 16), enable_reduction=False)
    add_statistics_panel(infection_rates_sim2, household_data,
                         position='right', panel_title=sim2_name.upper())

    # Add main title
    fig.suptitle('Infection Pattern Comparison',
                 fontsize=18, fontweight='bold', y=0.95)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def create_difference_visualization(contact_data, household_data,
                                    infection_rates_sim1, infection_rates_sim2,
                                    sim1_name="Simulation 1", sim2_name="Simulation 2",
                                    save_path=None):
    """Show differences between two simulations"""

    # Calculate differences
    differences = {}
    all_persons = set(infection_rates_sim1.keys()) | set(
        infection_rates_sim2.keys())

    for person in all_persons:
        rate1 = infection_rates_sim1.get(person, 0)
        rate2 = infection_rates_sim2.get(person, 0)
        differences[person] = rate2 - rate1  # Positive = more infected in sim2

    def get_difference_color(diff):
        """Blue for negative (less in sim2), Red for positive (more in sim2), White for no change"""
        if abs(diff) < 5:  # Minimal difference
            return '#F0F0F0'  # Light gray
        elif diff > 0:  # Higher in sim2
            # Scale red intensity based on difference magnitude
            intensity = min(abs(diff) / 100, 1.0)
            return plt.cm.Reds(0.3 + 0.7 * intensity)
        else:  # Higher in sim1
            # Scale blue intensity based on difference magnitude
            intensity = min(abs(diff) / 100, 1.0)
            return plt.cm.Blues(0.3 + 0.7 * intensity)

    plt.figure(figsize=(20, 14))

    # Create layouts
    contact_strength, contact_locations = calculate_contact_strength(
        contact_data)
    person_pos, households, household_positions = create_household_layout(
        household_data, household_spacing=12.0, person_spacing=1.8, max_households_per_row=15)
    location_pos, location_nodes, location_areas = create_location_layout(
        contact_data, household_positions, area_distance=60.0)
    all_pos = {**person_pos, **location_pos}

    # Create graph
    G = nx.Graph()
    for person_id in person_pos.keys():
        G.add_node(person_id)

    # Draw basic network elements
    draw_basic_network_elements(
        plt.gca(), all_pos, location_nodes, location_areas, household_positions, households)

    # Draw nodes with difference colors
    node_colors = [get_difference_color(
        differences.get(p, 0)) for p in person_pos.keys()]

    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=node_colors,
                           node_size=120,
                           alpha=0.9,
                           edgecolors='black',
                           linewidths=1.5)

    # Add difference annotations for significant changes
    for person_id, diff in differences.items():
        if abs(diff) > 30 and person_id in person_pos:  # Only annotate large differences
            x, y = person_pos[person_id]
            sign = "+" if diff > 0 else ""
            plt.text(x, y + 3, f"{sign}{diff:.0f}%",
                     ha='center', va='center', fontsize=6, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # Create custom legend for differences
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=plt.cm.Blues(0.7),
                   markersize=12, label=f'Higher in {sim1_name}', linestyle='None'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#F0F0F0',
                   markersize=12, label='No significant difference', linestyle='None'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=plt.cm.Reds(0.7),
                   markersize=12, label=f'Higher in {sim2_name}', linestyle='None')
    ]

    plt.legend(handles=legend_elements, loc='upper left',
               title='Infection Rate Differences', title_fontsize=12,
               bbox_to_anchor=(0.02, 0.98))

    # Add statistics
    total_diff = sum(abs(d) for d in differences.values())
    major_increases = sum(1 for d in differences.values() if d > 20)
    major_decreases = sum(1 for d in differences.values() if d < -20)

    diff_stats = f"""DIFFERENCE ANALYSIS
    
Total absolute difference: {total_diff:.1f}%
Major increases ({sim2_name}): {major_increases} persons
Major decreases ({sim2_name}): {major_decreases} persons

Average difference: {np.mean(list(differences.values())):.1f}%
Max increase: {max(differences.values()):.1f}%
Max decrease: {min(differences.values()):.1f}%"""

    plt.text(0.98, 0.02, diff_stats, transform=plt.gca().transAxes,
             fontsize=10, ha='right', va='bottom', fontfamily='monospace',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='white', alpha=0.95,
                       edgecolor='gray', linewidth=1))

    plt.title(f'Infection Rate Differences: {sim2_name} vs {sim1_name}\n' +
              f'(Red = Higher in {sim2_name}, Blue = Higher in {sim1_name})',
              fontsize=14, fontweight='bold', pad=20)
    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def load_contact_data(contact_file):
    """Load contact data from CSV file"""
    return pd.read_csv(contact_file)


def load_infection_data(infection_file):
    """Load infection data from CSV file"""
    return pd.read_csv(infection_file)


def load_household_data(household_file):
    """Load household assignment data from CSV file"""
    return pd.read_csv(household_file)


def calculate_infection_rates(infection_data, total_runs=25):
    """
    Calculate infection rates for each person across multiple runs

    Parameters:
    infection_data: DataFrame with columns ['person_id', 'run_id', 'infected', 'infection_time']
    total_runs: total number of simulation runs
    """
    # Count infections per person
    infection_counts = infection_data.groupby('person_id')['infected'].sum()

    # Calculate infection rates (percentage)
    infection_rates = (infection_counts / total_runs) * 100

    # Ensure all persons are included (even those never infected)
    all_persons = set(infection_data['person_id'].unique())
    for person in all_persons:
        if person not in infection_rates:
            infection_rates[person] = 0.0

    return infection_rates


def calculate_contact_strength(contact_data):
    """
    Calculate contact strength between person pairs
    Returns dict with (person1, person2, location_type) -> contact_count
    """
    contact_strength = defaultdict(int)
    contact_locations = defaultdict(set)

    for _, row in contact_data.iterrows():
        person1, person2 = sorted([row['person_1'], row['person_2']])
        location_type = row['location_type']
        location_id = row['location_id']

        # Count contacts by location type
        key = (person1, person2, location_type)
        contact_strength[key] += 1
        contact_locations[key].add(location_id)

    return contact_strength, contact_locations


def create_household_layout(household_data, household_spacing=15.0, person_spacing=2.0, max_households_per_row=15):
    """
    Create compact household layout with multiple rows and data reduction

    Parameters:
    household_data: DataFrame with columns ['person_id', 'household_id']
    household_spacing: horizontal spacing between households
    person_spacing: spacing between persons within household
    max_households_per_row: maximum households per row before wrapping
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
    row_spacing = 12.0  # Vertical spacing between household rows

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

        # Store as tuple
        household_positions[household_id] = (
            household_center_x, household_center_y)

        # Arrange household members (reduce crowding by limiting visible members)
        n_members = len(members)
        if n_members == 1:
            pos[members[0]] = (household_center_x, household_center_y)
        elif n_members <= 6:  # Show all members for small households
            member_width = (n_members - 1) * person_spacing
            member_start_x = household_center_x - member_width / 2

            for j, person_id in enumerate(members):
                x_pos = member_start_x + j * person_spacing
                pos[person_id] = (x_pos, household_center_y)
        else:  # For large households, show representative members only
            # Show first 3, last 2, and indicate there are more
            visible_members = members[:3] + members[-2:]
            for j, person_id in enumerate(visible_members):
                if j < 3:
                    x_offset = (j - 1) * person_spacing
                else:  # last 2 members
                    x_offset = (j - 1.5) * person_spacing
                pos[person_id] = (household_center_x +
                                  x_offset, household_center_y)

            # Store info about hidden members for later use
            household_positions[f"{household_id}_hidden"] = len(members) - 5

    return pos, households, household_positions


def create_location_layout(contact_data, household_positions, area_distance=80.0, area_size=15.0):
    """
    Create layout for location nodes with adjusted distance for multi-row households

    Parameters:
    contact_data: DataFrame with contact information
    household_positions: dict mapping household_id to (x, y) position tuples
    area_distance: distance from household area to location areas
    area_size: size of each location area
    """
    location_types = ['Work', 'School', 'Social_Event', 'Shopping']

    # Get the range of household positions - now always tuples
    if household_positions:
        # Extract x and y coordinates from tuples with safety checks
        household_x_positions = [pos[0]
                                 for pos in household_positions.values()]
        household_y_positions = [pos[1]
                                 for pos in household_positions.values()]

        min_x = min(household_x_positions) if household_x_positions else 0
        max_x = max(household_x_positions) if household_x_positions else 0
        min_y = min(household_y_positions) if household_y_positions else 0
        max_y = max(household_y_positions) if household_y_positions else 0

        center_x = (min_x + max_x) / 2
        width = max_x - min_x
        household_area_height = max_y - min_y + 10  # Add padding
    else:
        center_x = 0
        width = 100
        household_area_height = 5

    # Define location areas above and below the household area
    location_areas = {
        # Below left
        'Work': (center_x - width/4, -(area_distance + household_area_height/2)),
        # Below right
        'School': (center_x + width/4, -(area_distance + household_area_height/2)),
        # Above left
        'Social_Event': (center_x - width/4, area_distance + household_area_height/2),
        # Above right
        'Shopping': (center_x + width/4, area_distance + household_area_height/2)
    }

    # Get unique locations for each type
    location_positions = {}
    location_nodes = set()

    for location_type in location_types:
        if location_type in contact_data['location_type'].values:
            # Get unique location IDs for this type
            type_locations = contact_data[
                contact_data['location_type'] == location_type
            ]['location_id'].unique()

            area_x, area_y = location_areas[location_type]

            # Arrange locations with spacing
            n_locations = len(type_locations)
            if n_locations == 1:
                # Single location at area center
                location_node = f"{location_type}_{type_locations[0]}"
                location_positions[location_node] = (area_x, area_y)
                location_nodes.add(location_node)
            else:
                # Multiple locations in a row
                spacing = 6.0  # Spacing between location nodes

                for i, location_id in enumerate(type_locations):
                    x_offset = (i - (n_locations - 1) / 2) * spacing
                    y_offset = 0  # Keep all at same height for each type

                    location_node = f"{location_type}_{location_id}"
                    location_positions[location_node] = (
                        area_x + x_offset, area_y + y_offset)
                    location_nodes.add(location_node)

    return location_positions, location_nodes, location_areas


def apply_data_reduction(contact_data, infection_data, household_data,
                         reduction_strategy='smart', max_persons=200,
                         min_infection_rate=10, min_contacts=2):
    """
    Apply data reduction techniques to make visualization manageable

    Parameters:
    contact_data: DataFrame with contact information
    infection_data: DataFrame with infection data
    household_data: DataFrame with household assignments
    reduction_strategy: 'smart', 'random', 'high_risk', or 'representative'
    max_persons: maximum number of persons to include
    min_infection_rate: minimum infection rate (%) to include person
    min_contacts: minimum number of contacts to include person

    Returns:
    Filtered versions of all three dataframes
    """

    # Calculate person-level statistics
    person_stats = {}

    # Infection rates
    infection_rates = infection_data.groupby(
        'person_id')['infected'].mean() * 100

    # Contact counts
    contact_counts = {}
    for _, row in contact_data.iterrows():
        person1, person2 = row['person_1'], row['person_2']
        contact_counts[person1] = contact_counts.get(person1, 0) + 1
        contact_counts[person2] = contact_counts.get(person2, 0) + 1

    # Combine statistics
    all_persons = set(household_data['person_id'].unique())
    for person_id in all_persons:
        person_stats[person_id] = {
            'infection_rate': infection_rates.get(person_id, 0),
            'contact_count': contact_counts.get(person_id, 0),
            'household_id': household_data[household_data['person_id'] == person_id]['household_id'].iloc[0]
        }

    # Apply reduction strategy
    if reduction_strategy == 'smart':
        selected_persons = _smart_selection(
            person_stats, household_data, max_persons)
    elif reduction_strategy == 'random':
        selected_persons = _random_selection(all_persons, max_persons)
    elif reduction_strategy == 'high_risk':
        selected_persons = _high_risk_selection(
            person_stats, max_persons, min_infection_rate)
    elif reduction_strategy == 'representative':
        selected_persons = _representative_selection(
            person_stats, household_data, max_persons)
    else:
        selected_persons = list(all_persons)[:max_persons]

    # Filter data
    filtered_household = household_data[household_data['person_id'].isin(
        selected_persons)]

    # Keep only households that still have members
    remaining_households = filtered_household['household_id'].unique()
    filtered_household = filtered_household[filtered_household['household_id'].isin(
        remaining_households)]

    # Filter contacts (only between selected persons)
    filtered_contact = contact_data[
        (contact_data['person_1'].isin(selected_persons)) &
        (contact_data['person_2'].isin(selected_persons))
    ]

    # Filter infections
    filtered_infection = infection_data[infection_data['person_id'].isin(
        selected_persons)]

    print(
        f"Data reduction: {len(all_persons)} → {len(selected_persons)} persons")
    print(
        f"Households: {household_data['household_id'].nunique()} → {filtered_household['household_id'].nunique()}")
    print(f"Contacts: {len(contact_data)} → {len(filtered_contact)}")

    return filtered_contact, filtered_infection, filtered_household


def _smart_selection(person_stats, household_data, max_persons):
    """Smart selection prioritizing interesting persons and household diversity"""
    # Score each person
    scored_persons = []

    for person_id, stats in person_stats.items():
        # Score based on multiple factors
        infection_score = stats['infection_rate'] / 100  # 0-1
        contact_score = min(stats['contact_count'] /
                            50, 1.0)  # Cap at 50 contacts

        # Bonus for extreme cases (very high or very low infection rates)
        extremity_bonus = abs(stats['infection_rate'] - 50) / 50

        total_score = infection_score * 0.4 + contact_score * 0.4 + extremity_bonus * 0.2
        scored_persons.append((person_id, total_score, stats['household_id']))

    # Sort by score
    scored_persons.sort(key=lambda x: x[1], reverse=True)

    # Select persons ensuring household diversity
    selected_persons = []
    selected_households = set()
    household_counts = {}

    for person_id, score, household_id in scored_persons:
        if len(selected_persons) >= max_persons:
            break

        # Limit persons per household to maintain diversity
        household_count = household_counts.get(household_id, 0)
        max_per_household = 4  # Limit to 4 persons per household in visualization

        if household_count < max_per_household:
            selected_persons.append(person_id)
            selected_households.add(household_id)
            household_counts[household_id] = household_count + 1

    return selected_persons


def _random_selection(all_persons, max_persons):
    """Random selection for baseline comparison"""
    import random
    return random.sample(list(all_persons), min(max_persons, len(all_persons)))


def _high_risk_selection(person_stats, max_persons, min_infection_rate):
    """Select high-risk persons (high infection rates)"""
    high_risk = [pid for pid, stats in person_stats.items()
                 if stats['infection_rate'] >= min_infection_rate]
    return high_risk[:max_persons]


def _representative_selection(person_stats, household_data, max_persons):
    """Select representative sample maintaining household structure"""
    # Group by household size
    household_sizes = household_data.groupby('household_id').size()

    # Sample households proportionally by size
    selected_households = []

    # Select some small, medium, and large households
    small_hh = household_sizes[household_sizes <= 2].index[:5]
    medium_hh = household_sizes[(household_sizes > 2) & (
        household_sizes <= 4)].index[:8]
    large_hh = household_sizes[household_sizes > 4].index[:4]

    selected_households = list(small_hh) + list(medium_hh) + list(large_hh)

    # Get all persons from selected households
    selected_persons = household_data[
        household_data['household_id'].isin(selected_households)
    ]['person_id'].tolist()

    return selected_persons[:max_persons]


def aggregate_location_contacts(contact_data, min_contacts_per_location=5):
    """
    Aggregate low-frequency location contacts to reduce visual clutter

    Parameters:
    contact_data: DataFrame with contact data
    min_contacts_per_location: minimum contacts needed to show specific location

    Returns:
    Modified contact_data with aggregated locations
    """
    # Count contacts per location
    location_counts = contact_data.groupby(
        ['location_type', 'location_id']).size()

    # Identify low-frequency locations to aggregate
    low_freq_locations = location_counts[location_counts <
                                         min_contacts_per_location].index

    # Aggregate low-frequency locations by type
    modified_contact_data = contact_data.copy()

    for location_type, location_id in low_freq_locations:
        mask = ((modified_contact_data['location_type'] == location_type) &
                (modified_contact_data['location_id'] == location_id))

        # Replace with aggregated location
        modified_contact_data.loc[mask,
                                  'location_id'] = f"{location_type}_Other"

    print(
        f"Location aggregation: {contact_data['location_id'].nunique()} → {modified_contact_data['location_id'].nunique()} unique locations")

    return modified_contact_data


def get_location_edge_colors():
    """Define colors for different location types for edges"""
    return {
        'Home': '#2E8B57',        # Sea Green for household
        'Work': '#4169E1',        # Royal Blue
        'School': '#FF6347',      # Tomato Red
        'Social_Event': '#9370DB',  # Medium Purple
        'Shopping': '#FF8C00'     # Dark Orange
    }


def get_infection_color_v2(infection_rate):
    """Better color palette: light yellow → warm orange → deep red"""
    # Create custom colormap with more pleasing colors
    colors = ['#FFF8DC', '#FFE4B5', '#FFA500',
              '#FF6347', '#B22222']  # Cornsilk to FireBrick
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('infection', colors, N=n_bins)
    return cmap(infection_rate / 100)


def draw_curved_edges(G, all_pos, edges, width, color, alpha):
    """Draw edges with slight curves for more organic look"""
    nx.draw_networkx_edges(G, all_pos,
                           edgelist=edges,
                           width=width,
                           edge_color=color,
                           alpha=alpha,
                           connectionstyle="arc3,rad=0.1")  # Slight curve


def create_enhanced_legend(infection_rates, edge_colors):
    """Create a more informative legend with sample nodes and better styling"""
    legend_elements = []

    # Infection rate samples with actual colored nodes
    for rate, label in [(0, 'Never infected'), (50, 'Sometimes infected'), (100, 'Always infected')]:
        color = get_infection_color_v2(rate)
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          markeredgecolor='black', markeredgewidth=0.5,
                                          label=label, linestyle='None'))

    # Add separator
    legend_elements.append(plt.Line2D([0], [0], color='none', label=''))

    # Location type lines with better styling
    for loc_type, color in edge_colors.items():
        if loc_type != 'Home':  # Skip home since it's rarely visible
            legend_elements.append(plt.Line2D([0], [0], color=color, lw=3,
                                              label=f'{loc_type.replace("_", " ")} Contacts'))

    # Create legend with enhanced styling
    legend = plt.legend(handles=legend_elements, loc='upper left',
                        bbox_to_anchor=(0.02, 0.98),
                        frameon=True, fancybox=True, shadow=True,
                        fontsize=9, title='Legend', title_fontsize=10)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_edgecolor('gray')


def visualize_household_network(contact_data, household_data, infection_rates,
                                figsize=(28, 20), total_runs=25, save_path=None,
                                enable_reduction=True, max_persons=150,
                                reduction_strategy='smart'):
    """
    Create the household-centered network visualization with enhanced styling
    """

    # Apply data reduction if requested
    if enable_reduction and len(household_data) > max_persons:
        print("Applying data reduction...")
        # Reconstruct infection_data from infection_rates for the reduction function
        infection_data_for_reduction = []
        for person_id, rate in infection_rates.items():
            # Create mock infection data based on rates
            for run_id in range(1, total_runs + 1):
                infected = 1 if (rate / 100) > 0.5 else 0  # Simple threshold
                infection_data_for_reduction.append({
                    'person_id': person_id,
                    'run_id': run_id,
                    'infected': infected,
                    'infection_time': ''
                })

        infection_data_df = pd.DataFrame(infection_data_for_reduction)

        # Apply the reduction
        contact_data, infection_data_filtered, household_data = apply_data_reduction(
            contact_data,
            infection_data_df,
            household_data,
            reduction_strategy=reduction_strategy,
            max_persons=max_persons
        )

        # Aggregate low-frequency locations
        contact_data = aggregate_location_contacts(
            contact_data, min_contacts_per_location=3)

    plt.figure(figsize=figsize)

    # Calculate contact strengths
    contact_strength, contact_locations = calculate_contact_strength(
        contact_data)

    # Create layouts with improved parameters
    person_pos, households, household_positions = create_household_layout(
        household_data, household_spacing=12.0, person_spacing=1.8, max_households_per_row=15)

    location_pos, location_nodes, location_areas = create_location_layout(
        contact_data, household_positions, area_distance=60.0)

    # Combine all positions
    all_pos = {**person_pos, **location_pos}

    # Create the graph
    G = nx.Graph()

    # Add person nodes
    for person_id in person_pos.keys():
        G.add_node(person_id, node_type='person')

    # Add location nodes
    for location_node in location_nodes:
        G.add_node(location_node, node_type='location')

    # Add household edges (direct connections) - but limit to avoid clutter
    household_edges = []
    for household_id, members in households.items():
        if len(members) <= 2:  # Only show direct connections for small households
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    person1, person2 = members[i], members[j]
                    strength = contact_strength.get(
                        (min(person1, person2), max(person1, person2), 'Home'), 0)
                    if strength > 0:
                        G.add_edge(person1, person2,
                                   edge_type='household', weight=strength)
                        household_edges.append((person1, person2))
        elif len(members) <= 6:  # For medium households, show some connections
            # Connect each person to household "center" person (first person)
            center_person = members[0]
            for member in members[1:]:
                strength = contact_strength.get(
                    (min(center_person, member), max(center_person, member), 'Home'), 0)
                if strength > 0:
                    G.add_edge(center_person, member,
                               edge_type='household', weight=strength)
                    household_edges.append((center_person, member))

    # Add location-mediated edges (but filter for important ones only)
    location_edges = []
    person_location_edges = []

    for (person1, person2, location_type), strength in contact_strength.items():
        if location_type != 'Home' and strength >= 2:  # Only show frequent contacts
            locations = contact_locations[(person1, person2, location_type)]

            for location_id in locations:
                location_node = f"{location_type}_{location_id}"
                if location_node in location_nodes:
                    G.add_edge(person1, location_node,
                               edge_type='person_to_location', weight=strength)
                    G.add_edge(person2, location_node,
                               edge_type='person_to_location', weight=strength)

                    person_location_edges.append((person1, location_node))
                    person_location_edges.append((person2, location_node))

    # Draw household areas with rounded corners
    for household_id, household_pos in household_positions.items():
        if '_hidden' in household_id:
            continue

        household_center_x, household_center_y = household_pos
        members = households[household_id]

        if len(members) > 1:
            # Rounded rectangles with subtle shadow
            fancy_box = FancyBboxPatch(
                (household_center_x - 6, household_center_y - 1), 12, 2,
                boxstyle="round,pad=0.1",
                facecolor='lightgray',
                alpha=0.25,
                edgecolor='darkgray',
                linewidth=0.8,
                zorder=0
            )
            plt.gca().add_patch(fancy_box)

            # Add household size indicator for large households
            if len(members) > 6:
                plt.text(household_center_x + 5, household_center_y, f"({len(members)})",
                         fontsize=8, alpha=0.7, fontweight='bold')

    # Draw location areas (smaller)
    for location_type, (area_x, area_y) in location_areas.items():
        if any(location_type in node for node in location_nodes):
            rect = plt.Rectangle((area_x - 8, area_y - 8), 16, 16,
                                 facecolor='lightblue', alpha=0.1, zorder=0,
                                 edgecolor='blue', linewidth=1)
            plt.gca().add_patch(rect)

            if area_y > 0:  # Above household line
                label_y = area_y + 10
            else:  # Below household line
                label_y = area_y - 10

            plt.text(area_x, label_y, location_type.replace('_', ' '),
                     ha='center', va='center', fontsize=11, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    # Get edge colors
    edge_colors = get_location_edge_colors()

    # Draw person nodes with enhanced styling
    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=[get_infection_color_v2(
                               infection_rates.get(p, 0)) for p in person_pos.keys()],
                           node_size=120,  # Slightly larger
                           alpha=0.9,
                           edgecolors='white',  # White borders for better contrast
                           linewidths=1.5)

    # Draw location nodes (smaller)
    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(location_nodes),
                           node_color='lightsteelblue',
                           node_size=80,
                           alpha=0.7,
                           node_shape='s',
                           edgecolors='navy',
                           linewidths=0.5)

    # Draw household edges with curves
    if household_edges:
        draw_curved_edges(G, all_pos, household_edges,
                          0.8, edge_colors['Home'], 0.6)

    # Draw person-to-location edges with curves, colored by type
    if person_location_edges:
        edges_by_type = defaultdict(list)

        for person, location_node in person_location_edges:
            location_type = location_node.split('_')[0]
            edges_by_type[location_type].append((person, location_node))

        for location_type, edges in edges_by_type.items():
            color = edge_colors.get(location_type, 'gray')
            draw_curved_edges(G, all_pos, edges, 0.2, color, 0.3)

    # Only draw labels for locations (no person labels for cleaner look)
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

    # Enhanced legend
    edge_colors = get_location_edge_colors()
    create_enhanced_legend(infection_rates, edge_colors)

    # Add title with data info
    num_persons = len(person_pos)
    num_households = len(
        [h for h in household_positions.keys() if '_hidden' not in h])

    plt.title(f'Household-Centered Contact Network\n' +
              f'{num_persons} persons, {num_households} households ({total_runs} simulation runs)',
              fontsize=14, pad=20)
    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()


def print_network_stats(contact_data, infection_data, household_data, total_runs=25):
    """Print comprehensive network statistics"""

    print("=== Network Statistics ===")
    print(f"Total persons: {len(household_data)}")
    print(f"Total households: {household_data['household_id'].nunique()}")
    print(f"Total contacts: {len(contact_data)}")
    print(f"Simulation runs: {total_runs}")

    # Household statistics
    household_sizes = household_data.groupby('household_id').size()
    print(f"\n=== Household Statistics ===")
    print(f"Average household size: {household_sizes.mean():.2f}")
    print(
        f"Household size range: {household_sizes.min()} - {household_sizes.max()}")

    # Contact statistics by location type
    print(f"\n=== Contact Statistics by Location ===")
    location_stats = contact_data.groupby('location_type').size()
    for location_type, count in location_stats.items():
        print(f"{location_type}: {count} contacts")

    # Infection statistics
    infection_rates = calculate_infection_rates(infection_data, total_runs)
    if hasattr(infection_rates, 'values'):
        rates = list(infection_rates.values)
    else:
        rates = list(infection_rates.values())

    print(f"\n=== Infection Statistics ===")
    print(f"Average infection rate: {np.mean(rates):.2f}%")
    print(f"Max infection rate: {max(rates):.2f}%")
    print(f"Persons never infected: {sum(1 for r in rates if r == 0)}")
    print(f"Persons always infected: {sum(1 for r in rates if r == 100)}")


def main():
    """Main function to run the analysis"""

    # Load data
    base_dir = '/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/viz'
    contact_data = load_contact_data(f'{base_dir}/medium_contact_data.csv')
    infection_data = load_infection_data(
        f'{base_dir}/medium_infection_data.csv')
    household_data = load_household_data(
        f'{base_dir}/medium_household_data.csv')

    contact_data = load_contact_data(f'{base_dir}/medium_contact_data.csv')
    infection_data2 = load_infection_data(
        f'{base_dir}/medium_infection_data2.csv')
    household_data2 = load_household_data(
        f'{base_dir}/medium_household_data2.csv')

    total_runs = 25

    print(f"Loaded {len(contact_data)} contact records")
    print(f"Loaded {len(infection_data)} infection records")
    print(f"Loaded {len(household_data)} household assignments")

    # Calculate infection rates
    print("Calculating infection rates...")
    infection_rates = calculate_infection_rates(infection_data, total_runs)

    # Print statistics
    print_network_stats(contact_data, infection_data,
                        household_data, total_runs)

    # Create visualization with enhanced styling
    print("\nCreating enhanced household-centered visualization...")
    visualize_household_network(contact_data, household_data, infection_rates,
                                figsize=(20, 14),
                                total_runs=total_runs,
                                enable_reduction=True,
                                max_persons=150,
                                reduction_strategy='smart',
                                save_path='household_network.png')

    # # Create time-filtered visualization (example: first 24 hours)
    # print("\nCreating time-filtered visualization (first 24 hours)...")
    # contact_data_24h = contact_data[contact_data['time'] <= 24]
    # visualize_household_network(contact_data_24h, household_data, infection_rates,
    #                             figsize=(18, 12),
    #                             total_runs=total_runs,
    #                             enable_reduction=True,
    #                             max_persons=100,
    #                             save_path='household_network_24h.png')

    # print("Analysis complete!")


if __name__ == "__main__":
    main()
