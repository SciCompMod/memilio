
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import argparse
import os
import sys
from collections import defaultdict, Counter

# Import functions from the existing scripts
# Make sure contact_network.py and time_loc_type.py are in the same directory
try:
    from contact_network import (
        load_household_data, load_infection_count_data, load_contact_data,
        load_location_data, calculate_infection_rates, prepare_contact_data_with_locations,
        create_household_layout, create_location_layout, get_location_edge_colors,
        get_infection_color, get_infection_outline, visualize_contact_network,
        calculate_contact_strength_from_hours, get_edge_width_from_hours
    )
    from time_loc_type import (
        read_and_analyze_data, calculate_group_averages, create_pie_charts,
        print_detailed_summary
    )
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure contact_network.py and time_loc_type.py are in the same directory")
    sys.exit(1)

import networkx as nx
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch


def parse_city_config(config_file):
    """Parse the city configuration file to extract demographic data"""
    config_data = {}
    try:
        with open(config_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = line.strip()
            if ',' in line and not line.startswith('POPULATION') and not line.startswith('---'):
                key, value = line.split(',', 1)
                try:
                    if '.' in value:
                        config_data[key] = float(value)
                    else:
                        config_data[key] = int(value)
                except ValueError:
                    config_data[key] = value
    except FileNotFoundError:
        print(
            f"Warning: City config file {config_file} not found. Using simulation data only.")
        return None
    except Exception as e:
        print(f"Error reading city config: {e}")
        return None
    return config_data


def create_age_distribution_pie():
    """Create age distribution pie chart from city config"""
    return {
        'categories': ['0-4y\n(44)', '5-14y\n(94)', '15-34y\n(222)',
                       '35-59y\n(334)', '60-79y\n(233)', '80+y\n(73)'],
        'counts': [44, 94, 222, 334, 233, 73],
        'colors': ['#FFA07A', '#20B2AA', '#DA70D6', '#32CD32', "#989606", "#FE8515"]
    }


def create_household_distribution_pie():
    """Create household size distribution pie chart from city config"""
    return {
        'categories': ['1p\n(205)', '2p\n(172)', '3p\n(59)', '4p\n(46)', '5p\n(18)'],
        'counts': [205, 172, 59, 46, 18],
        'colors': ['#FFA07A', '#20B2AA', '#DA70D6', '#32CD32', "#989606"]
    }


def create_occupation_distribution_pie():
    """Create occupation distribution pie chart from city config"""
    return {
        'categories': ['Workers\n(431)', 'Pupils\n(105)', 'Others\n(464)'],
        'counts': [431, 105, 464],
        'colors': ['#4169E1', '#FF5805', '#BA55D3']
    }


def create_contact_network_subplot(ax, household_data, infection_count_data, contact_data,
                                   location_data, total_runs=5, min_contact_hours=5,
                                   max_persons=100, viz_home_contacts=False):
    """Create contact network visualization in a subplot by calling the core logic directly"""

    # Performance optimization: check dataset size and apply intelligent filtering
    num_persons = len(household_data)
    num_contacts = len(contact_data)

    print(
        f"Processing contact network: {num_persons} persons, {num_contacts} contacts...")

    if num_persons > 1000:
        print(
            f"Large dataset detected ({num_persons} persons). Applying heavy filtering for clarity...")
        # For any dataset over 300, sample much fewer households
        unique_households = household_data['household_id'].unique()
        if len(unique_households) > 50:
            sampled_households = np.random.choice(
                unique_households, 50, replace=False)
            household_data = household_data[household_data['household_id'].isin(
                sampled_households)]
            print(
                f"Sampled {len(sampled_households)} households, {len(household_data)} persons")

    # Merge contact data with location types
    contact_data_with_types = prepare_contact_data_with_locations(
        contact_data, location_data)

    # Filter contact data to only include remaining persons
    remaining_persons = household_data['person_id'].tolist()
    contact_data_with_types = contact_data_with_types[
        (contact_data_with_types['person_1'].isin(remaining_persons)) &
        (contact_data_with_types['person_2'].isin(remaining_persons))
    ]

    # Performance optimization: aggressive contact filtering for large datasets
    if len(contact_data_with_types) > 5000:
        print(
            f"Large contact dataset ({len(contact_data_with_types)} contacts). Increasing filter threshold...")
        min_contact_hours = max(min_contact_hours, 15)

    # Filter contacts by minimum hours
    contact_data_with_types = contact_data_with_types[
        ((contact_data_with_types['contact_hours'] >= 50)
         & (contact_data_with_types['location_type'].isin(['Work', 'School', 'SocialEvent'])))
        | ((contact_data_with_types['contact_hours'] >= 5)
           & (contact_data_with_types['location_type'].isin(['BasicsShop'])))
    ]

    # Filter out Hospital, ICU, and Unknown location types
    contact_data_with_types = contact_data_with_types[
        ~contact_data_with_types['location_type'].isin(
            ['Hospital', 'ICU', 'Unknown'])
    ]

    # Filter out homes, if not visualizing home contacts
    if not viz_home_contacts:
        contact_data_with_types = contact_data_with_types[
            contact_data_with_types['location_type'] != 'Home'
        ]

    print(
        f"Filtered to {len(household_data)} persons with {len(contact_data_with_types)} significant contacts")

    # Performance optimization: skip edge drawing for very large datasets
    draw_edges = len(contact_data_with_types) < 20000
    if not draw_edges:
        print("Skipping edge drawing for performance (too many contacts)")
    else:
        print(f"Will draw edges for {len(contact_data_with_types)} contacts")

    # Calculate contact strengths only if we're drawing edges
    if draw_edges:
        contact_strength, contact_locations = calculate_contact_strength_from_hours(
            contact_data_with_types)
    else:
        contact_strength, contact_locations = {}, {}

    spacing_factor = 0.22

    household_spacing = 500.0 * spacing_factor
    person_spacing = 80.0 * spacing_factor

    # Use the original create_household_layout but with much cleaner parameters
    person_pos, households, household_positions = create_household_layout(
        household_data,
        household_spacing=household_spacing,
        person_spacing=person_spacing,
        # Much fewer households per row for clarity
        max_households_per_row=16
    )

    location_pos, location_nodes, location_areas = create_location_layout(
        # Adjusted for better spacing
        contact_data_with_types, household_positions, area_distance=3500.0 * spacing_factor)

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

    # Draw location areas with proper bounding boxes - more prominent
    for location_type, (area_x, area_y) in location_areas.items():
        if any(location_type in node for node in location_nodes):
            # Count locations for sizing
            type_locations = contact_data_with_types[
                contact_data_with_types['location_type'] == location_type
            ]['location_id'].unique()
            n_locations = len(type_locations)

            # More prominent rectangles
            rect_width = max(12, min(n_locations * 4 + 4, 25))
            rect_height = max(8, min(10, 15))

            # rect = plt.Rectangle((area_x - rect_width/2, area_y - rect_height/2),
            #                      rect_width, rect_height,
            #                      facecolor='lightblue', alpha=0.2, zorder=0,  # More visible
            #                      edgecolor='blue', linewidth=1.5)  # Thicker border
            # ax.add_patch(rect)

            # More prominent labels
            label_y = area_y + rect_height/2 + \
                200 if area_y > 0 else area_y - rect_height/2 - 150
            ax.text(area_x, label_y, location_type.replace('_', ' '),
                    ha='center', va='center', fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9,
                              edgecolor='blue', linewidth=1))

    # Draw person nodes with bigger size
    person_colors = ['lightblue' for _ in person_pos.keys()]
    person_outlines = ['darkblue' for _ in person_pos.keys()]

    # Bigger nodes that scale better
    node_size = max(60, min(120, 8000 / len(person_pos)))  # Larger base size

    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(person_pos.keys()),
                           node_color=person_colors,
                           node_size=node_size,
                           alpha=0.8,
                           edgecolors=person_outlines,
                           linewidths=1.5,  # Thicker outline
                           ax=ax)

    # Draw location nodes - also bigger
    location_node_size = max(
        40, min(80, 3000 / len(location_nodes)) if location_nodes else 80)
    nx.draw_networkx_nodes(G, all_pos,
                           nodelist=list(location_nodes),
                           node_color='lightsteelblue',
                           node_size=location_node_size*1.2,
                           alpha=0.8,  # More opaque
                           node_shape='s',
                           edgecolors='navy',
                           linewidths=1.0,
                           ax=ax)

    # Draw edges only for manageable datasets
    if draw_edges and contact_strength:
        # Get edge colors and max contact hours for width scaling
        edge_colors = get_location_edge_colors()
        max_contact_hours = max(contact_data_with_types['contact_hours']) if len(
            contact_data_with_types) > 0 else 100

        # Track edges per location type
        edges_per_location_type = defaultdict(int)
        max_edges_per_type = 1000

        # Sort contacts by total hours (strongest first)
        sorted_contacts = sorted(
            contact_strength.items(), key=lambda x: x[1], reverse=True)

        for (person1, person2, location_type), total_hours in sorted_contacts:
            # Check if we've reached the limit for this location type
            if edges_per_location_type[location_type] >= max_edges_per_type:
                continue

            locations = contact_locations[(person1, person2, location_type)]

            for location_id in locations:
                if edges_per_location_type[location_type] >= max_edges_per_type:
                    break

                location_node = f"{location_type}_{location_id}"
                if location_node in location_nodes:
                    # Slightly more visible edges while still subtle
                    edge_width = max(
                        0.3, min(1.2, (total_hours / max_contact_hours) * 3.0))
                    color = edge_colors.get(location_type, 'gray')

                    nx.draw_networkx_edges(G, all_pos,
                                           edgelist=[
                                               (person1, location_node), (person2, location_node)],
                                           width=edge_width,
                                           edge_color=color,
                                           alpha=0.25,
                                           ax=ax)
                    edges_per_location_type[location_type] += 2

    if len(location_nodes) < 100:
        location_labels = {loc: loc.split('_')[-1] for loc in location_nodes}
        nx.draw_networkx_labels(G, all_pos,
                                labels=location_labels,
                                font_size=max(4, min(6, 150 / len(location_nodes))), ax=ax)

    # Set title and formatting
    num_persons = len(person_pos)
    num_households = len(
        [h for h in household_positions.keys() if '_hidden' not in str(h)])

    ax.set_title(f'Contact Network\n{num_persons} agents, {num_households} households',
                 fontsize=11, fontweight='bold', pad=8)
    ax.axis('equal')
    ax.axis('off')

    print(
        f"Contact network visualization complete: {num_persons} persons, {len(location_nodes)} locations")


def create_single_pie_chart(ax, person_totals, group_people, location_types,
                            num_days, title, group_name):
    """Create a single pie chart with improved formatting"""
    # Define consistent color order matching location order: Home, Work, School, SocialEvent, BasicsShop
    color_map = {
        'Home': '#2E8B57',        # Green
        'Work': '#4169E1',        # Blue
        'School': "#FF5805",      # Red
        'SocialEvent': '#9370DB',  # Purple
        'BasicsShop': "#FFCB2E"   # Orange
    }

    # Define consistent order: Home, Work, School, SocialEvent, BasicsShop
    desired_order = ['Home', 'Work', 'School', 'SocialEvent', 'BasicsShop']

    if not group_people:
        ax.text(0.5, 0.5, f'No {group_name}', ha='center', va='center',
                transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontsize=10, fontweight='bold', pad=8)
        ax.axis('off')
        return

    # Calculate averages
    averages = calculate_group_averages(person_totals, group_people, num_days)

    # Create location to average mapping
    location_to_avg = dict(zip(location_types, averages))

    # Filter out categories with 0 hours and organize in desired order
    non_zero_data = []
    non_zero_labels = []
    non_zero_colors = []

    # Add data in desired order (only if > 0)
    for loc_type in desired_order:
        if loc_type in location_to_avg and location_to_avg[loc_type] > 0:
            avg = location_to_avg[loc_type]
            non_zero_data.append(avg)

            # Shorter labels for better fit
            if loc_type == 'SocialEvent':
                short_label = 'Social'
            elif loc_type == 'BasicsShop':
                short_label = 'Shop'
            else:
                short_label = loc_type
            non_zero_labels.append(f"{short_label}\n{avg:.1f}h")
            non_zero_colors.append(color_map.get(loc_type, '#999999'))

    # Add any remaining location types not in desired order
    remaining_locations = set(location_types) - set(desired_order)
    for i, loc_type in enumerate(location_types):
        if loc_type in remaining_locations and location_to_avg[loc_type] > 0:
            avg = location_to_avg[loc_type]
            non_zero_data.append(avg)
            short_label = loc_type
            non_zero_labels.append(f"{short_label}\n{avg:.1f}h")
            # Gray for unknown location types
            non_zero_colors.append('#999999')

    if non_zero_data:
        # Calculate percentages for positioning logic
        total = sum(non_zero_data)
        percentages = [(val/total)*100 for val in non_zero_data]

        # Create pie chart with compact spacing
        wedges, texts, autotexts = ax.pie(non_zero_data,
                                          labels=non_zero_labels,
                                          colors=non_zero_colors,
                                          autopct='%1.1f%%',
                                          startangle=90,
                                          textprops={'fontsize': 11},
                                          labeldistance=1.2,  # Closer to pie
                                          pctdistance=0.7,
                                          radius=0.9)          # Larger pie radius

        # Adjust positioning for small slices (< 5%)
        small_slice_threshold = 5.0  # Slices smaller than 10%
        very_small_threshold = 0.0   # Slices smaller than 5%

        # Count small slices for progressive positioning
        small_slice_count = 0
        very_small_slice_count = 0

        for i, (text, autotext, percentage) in enumerate(zip(texts, autotexts, percentages)):
            if percentage < small_slice_threshold:
                # Adjust percentage text position progressively for small slices
                if percentage < very_small_threshold:
                    # Very small slices: move percentage text further from center
                    base_distance = 0.75
                    distance_increment = 0.15 * very_small_slice_count
                    pct_distance = base_distance + distance_increment
                    very_small_slice_count += 1
                else:
                    # Small slices: move percentage text moderately from center
                    base_distance = 0.75
                    distance_increment = 0.15 * small_slice_count
                    pct_distance = base_distance + distance_increment
                    # small_slice_count += 1

                # Get the wedge center angle for positioning
                wedge = wedges[i]
                theta = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)

                # Calculate new position for percentage text
                new_x = pct_distance * np.cos(theta)
                new_y = pct_distance * np.sin(theta)
                autotext.set_position((new_x, new_y))

                # Also handle outer label positioning for very small slices
                if percentage < 8:
                    next_percentage = percentages[(
                        i + 1) % len(percentages)] if len(percentages) > 1 else 100
                    if next_percentage < 8:
                        pos = text.get_position()
                        new_x = pos[0] * 1.1
                        new_y = pos[1] * 1.1 - 0.25
                        text.set_position((new_x, new_y))

        # Improved text formatting
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(10)

        # Adjust label positioning to avoid overlap
        for text in texts:
            text.set_fontsize(11)
            text.set_fontweight('normal')

    ax.set_title(title, fontsize=11, fontweight='bold', pad=5)


def create_config_pie_chart(ax, pie_data, title):
    """Create pie chart from config data with improved formatting"""
    if pie_data['counts']:
        # Calculate percentages for positioning logic
        total = sum(pie_data['counts'])
        percentages = [(val/total)*100 for val in pie_data['counts']]

        wedges, texts, autotexts = ax.pie(pie_data['counts'],
                                          labels=pie_data['categories'],
                                          colors=pie_data['colors'],
                                          autopct='%1.1f%%',
                                          startangle=90,
                                          textprops={'fontsize': 11},
                                          labeldistance=1.2,  # Closer to pie
                                          pctdistance=0.55,
                                          radius=0.9)          # Larger pie radius

        # Adjust positioning for small slices (< 5%) - compact version
        small_slice_threshold = 5.0  # Slices smaller than 10%
        very_small_threshold = 0.0   # Slices smaller than 5%

        # Count small slices for progressive positioning
        small_slice_count = 0
        very_small_slice_count = 0

        for i, (text, autotext, percentage) in enumerate(zip(texts, autotexts, percentages)):
            if percentage < small_slice_threshold:
                # Adjust percentage text position progressively for small slices
                if percentage < very_small_threshold:
                    # Very small slices: move percentage text further from center
                    base_distance = 0.7
                    distance_increment = 0.15 * very_small_slice_count
                    pct_distance = base_distance + distance_increment
                    # very_small_slice_count += 1
                else:
                    # Small slices: move percentage text moderately from center
                    base_distance = 0.65
                    distance_increment = 0.05 * small_slice_count
                    pct_distance = base_distance + distance_increment
                    # small_slice_count += 1

                # Get the wedge center angle for positioning
                wedge = wedges[i]
                theta = np.deg2rad((wedge.theta1 + wedge.theta2) / 2)

                # Calculate new position for percentage text
                new_x = pct_distance * np.cos(theta)
                new_y = pct_distance * np.sin(theta)
                autotext.set_position((new_x, new_y))

            # Special handling for the specific percentage range you mentioned
            if (percentage > 4.3 and percentage < 4.5):
                autotext_pos = autotext.get_position()
                autotext.set_position(
                    (autotext_pos[0] * 1.0, autotext_pos[1] * 1.05))
            elif percentage < 8:
                # Check if next slice is also small
                next_percentage = percentages[(
                    i + 1) % len(percentages)] if len(percentages) > 1 else 100
                if next_percentage < 5:
                    # Move outer text further away
                    pos = text.get_position()
                    new_x = pos[0] * 1.0
                    new_y = pos[1] * 1.0
                    text.set_position((new_x, new_y))

        # Improved text formatting
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(10)

        for text in texts:
            text.set_fontsize(11)
            text.set_fontweight('normal')

    ax.set_title(title, fontsize=11, fontweight='bold', pad=5)


def calculate_degree_distribution(contact_data):
    """Calculate the degree distribution from contact data"""
    print("  Calculating degree distribution...")

    # Count contacts per person
    degree_counts = Counter()

    # Process each contact record
    for _, row in contact_data.iterrows():
        person1 = row['Person 1'] if 'Person 1' in contact_data.columns else row.iloc[0]
        person2 = row['Person 2'] if 'Person 2' in contact_data.columns else row.iloc[1]

        degree_counts[person1] += 1
        degree_counts[person2] += 1

    # Convert to list of degrees
    degrees = list(degree_counts.values())

    print(f"  Calculated degrees for {len(degree_counts)} people")
    print(f"  Mean degree: {np.mean(degrees):.1f}")
    print(f"  Max degree: {max(degrees) if degrees else 0}")

    return degrees, degree_counts


def create_degree_distribution_plot(ax, contact_data, title="Potential Contact Distribution"):
    """Create degree distribution histogram"""
    degrees, degree_counts = calculate_degree_distribution(contact_data)

    if not degrees:
        ax.text(0.5, 0.5, 'No contact data', ha='center', va='center',
                transform=ax.transAxes, fontsize=12)
        ax.set_title(title, fontsize=11, fontweight='bold', pad=8)
        ax.axis('off')
        return

    # Create histogram
    bins = min(30, max(10, len(set(degrees))))  # Adaptive binning
    counts, bin_edges, patches = ax.hist(degrees, bins=bins, alpha=0.7, color='skyblue',
                                         edgecolor='navy', linewidth=0.5)

    # Color bars by frequency (darker = more frequent)
    max_count = max(counts)
    for i, (count, patch) in enumerate(zip(counts, patches)):
        # Create gradient from light blue to dark blue
        intensity = count / max_count
        patch.set_facecolor(plt.cm.Blues(0.3 + 0.7 * intensity))

    # Add statistics text box
    mean_deg = np.mean(degrees)
    std_deg = np.std(degrees)
    max_deg = max(degrees)
    median_deg = np.median(degrees)

    # Identify potential super-connectors (top 5% of degrees)
    threshold_95 = np.percentile(degrees, 95)
    super_connectors = sum(1 for d in degrees if d >= threshold_95)

    stats_text = f'Mean: {mean_deg:.1f}\nStd: {std_deg:.1f}\nMax: {max_deg}\nMedian: {median_deg:.1f}\nTop 5%: {super_connectors} people'
    ax.text(0.98, 0.9, stats_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      alpha=0.9, edgecolor='gray'),
            fontsize=12)

    # Formatting
    ax.set_xlabel('Number of Potential Contacts', fontsize=12)
    ax.set_ylabel('Number of People', fontsize=12)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=5)
    ax.grid(True, alpha=0.3)

    # Highlight super-connectors region if it exists
    if threshold_95 < max_deg:
        ax.axvline(threshold_95, color='red',
                   linestyle='--', alpha=0.7, linewidth=2)
        ax.text(threshold_95, ax.get_ylim()[1] * 0.8, f'95th percentile\n({threshold_95:.0f} contacts)',
                ha='left' if threshold_95 < max_deg * 0.7 else 'right',
                va='center', fontsize=10, color='black',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=1.0))

    print(
        f"  Degree distribution plot created: {len(degrees)} people, mean degree {mean_deg:.1f}")
    return degrees, degree_counts


def create_combined_visualization(data_dir, location_file, city_config_file=None, output_path=None,
                                  total_runs=5, min_contact_hours=5, max_persons=100):
    """Create the improved combined 6-panel visualization"""
    print("Loading data...")

    # Load contact network data
    try:
        household_data = load_household_data(
            os.path.join(data_dir, 'household_id.csv'))
        infection_count_data = load_infection_count_data(
            os.path.join(data_dir, 'infection_count.csv'))
        contact_data = load_contact_data(
            os.path.join(data_dir, 'contact_intensiveness.csv'))
        location_data = load_location_data(
            os.path.join(data_dir, 'location_id_and_type.csv'))
    except FileNotFoundError as e:
        print(f"Error: Required CSV file not found: {e}")
        return

    # Load time allocation data
    person_totals, location_types, num_days, workers, school_attendees = read_and_analyze_data(
        location_file)

    if person_totals is None:
        print("Failed to load time allocation data")
        return

    # Load city config data if provided
    config_data = None
    if city_config_file:
        config_data = parse_city_config(city_config_file)

    print(
        f"Loaded data for {len(household_data)} people across {num_days} days")
    print(f"Workers: {len(workers)}, Students: {len(school_attendees)}")

    # Print detailed summary
    print_detailed_summary(person_totals, location_types,
                           num_days, workers, school_attendees)

    # Create figure with compact layout - expanded for degree distribution
    # Slightly wider to accommodate new panel
    fig = plt.figure(figsize=(18, 12))

    # Create grid layout with minimal spacing - now 4x3 for degree distribution
    gs = fig.add_gridspec(4, 3, width_ratios=[1.8, 1, 1],
                          hspace=0.18, wspace=0.15,  # Slightly more vertical space
                          left=0.02, right=0.98, top=0.90, bottom=0.05)

    # Left panel: Contact network (spans first column, first 3 rows)
    ax_network = fig.add_subplot(gs[:4, 0])
    create_contact_network_subplot(ax_network, household_data, infection_count_data,
                                   contact_data, location_data, total_runs,
                                   min_contact_hours, max_persons)

    # Right panels: Pie charts in a 3x2 grid on the right side
    ax_all = fig.add_subplot(gs[0, 1])
    ax_workers = fig.add_subplot(gs[1, 1])
    ax_students = fig.add_subplot(gs[2, 1])
    ax_age = fig.add_subplot(gs[0, 2])
    ax_household = fig.add_subplot(gs[1, 2])
    ax_occupation = fig.add_subplot(gs[2, 2])

    ax_degree = fig.add_subplot(gs[3, 1:])  # Spans columns 1 and 2

    # Create time allocation pie charts
    all_people = list(person_totals.keys())

    # More compact titles
    create_single_pie_chart(ax_all, person_totals, all_people, location_types,
                            num_days, f'Average Time Spent (All) \n({len(all_people)} Agents)', 'Agents')

    create_single_pie_chart(ax_workers, person_totals, workers, location_types,
                            num_days, f'Average Time Spent (Workers) \n({len(workers)} Workers)', 'Workers')

    create_single_pie_chart(ax_students, person_totals, school_attendees, location_types,
                            num_days, f'Average Time Spent (Students) \n({len(school_attendees)} Students)', 'Students')

    # Create demographic pie charts from city config
    if config_data:
        age_pie = create_age_distribution_pie()
        household_pie = create_household_distribution_pie()
        occupation_pie = create_occupation_distribution_pie()

        create_config_pie_chart(
            ax_age, age_pie, 'Age Distribution\n(1000 agents)')
        create_config_pie_chart(
            ax_household, household_pie, 'Household Sizes\n(500 households)')
        create_config_pie_chart(
            ax_occupation, occupation_pie, 'Occupation Distribution\n(1000 agents)')
    else:
        # Fallback: create occupation chart from simulation data
        all_people = set(person_totals.keys())
        worker_set = set(workers)
        student_set = set(school_attendees)

        workers_only = worker_set - student_set
        students_only = student_set - worker_set
        both = worker_set & student_set
        others = all_people - worker_set - student_set

        # Create fallback occupation pie with matching colors
        sim_categories = []
        sim_counts = []
        sim_colors = ['#4169E1', '#FF5805', '#9370DB', '#2E8B57']

        if len(workers_only) > 0:
            sim_categories.append(f'Workers\n({len(workers_only)})')
            sim_counts.append(len(workers_only))
        if len(students_only) > 0:
            sim_categories.append(f'Students\n({len(students_only)})')
            sim_counts.append(len(students_only))
        if len(both) > 0:
            sim_categories.append(f'Worker-Students\n({len(both)})')
            sim_counts.append(len(both))
        if len(others) > 0:
            sim_categories.append(f'Others\n({len(others)})')
            sim_counts.append(len(others))

        # Create simulation-based occupation chart
        if sim_counts:
            occupation_sim_pie = {
                'categories': sim_categories,
                'counts': sim_counts,
                'colors': sim_colors[:len(sim_categories)]
            }
            create_config_pie_chart(ax_occupation, occupation_sim_pie,
                                    'Occupation Distribution\n(From simulation)')

        # Show empty charts for missing data
        for ax, title in [(ax_household, 'Household Sizes\n(No config data)')]:
            ax.text(0.5, 0.5, 'Config data\nnot available', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12)
            ax.set_title(title, fontsize=12, fontweight='bold', pad=20)
            ax.axis('off')

    # Create degree distribution plot (Panel H)
    print("Creating degree distribution visualization...")
    create_degree_distribution_plot(ax_degree, contact_data)

    # Add overall title with compact spacing
    fig.suptitle('Population Structure and Contact Network Analysis\n 7-day Baseline Simulation (No Infections)',
                 fontsize=16, fontweight='bold', y=1.0)

    # Add panel labels with compact positioning
    label_props = {'fontsize': 26,
                   'fontweight': 'bold', 'va': 'top', 'ha': 'left'}

    ax_network.text(0.01, 0.99, 'A',
                    transform=ax_network.transAxes, **label_props)
    ax_all.text(-0.35, 1.2, 'B', transform=ax_all.transAxes,
                **label_props)
    ax_workers.text(-0.35, 1.2, 'C',
                    transform=ax_workers.transAxes, **label_props)
    ax_students.text(
        -0.35, 1.2, 'D', transform=ax_students.transAxes, **label_props)
    ax_age.text(-0.35, 1.2, 'E', transform=ax_age.transAxes, **label_props)
    ax_household.text(
        -0.35, 1.2, 'F', transform=ax_household.transAxes, **label_props)
    ax_occupation.text(
        -0.35, 1.2, 'G', transform=ax_occupation.transAxes, **label_props)
    ax_degree.text(
        -0.0, 1.2, 'H', transform=ax_degree.transAxes, **label_props)

    if output_path:
        plt.savefig(output_path, dpi=400,
                    bbox_inches='tight', facecolor='white')
        print(f"Combined visualization saved to {output_path}")

    return fig


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Create improved combined population structure visualization')
    parser.add_argument('--data-dir', help='Directory containing CSV files')
    parser.add_argument('--output-path', help='Output path for visualization')
    parser.add_argument('--total-runs', type=int, default=5,
                        help='Total simulation runs')
    parser.add_argument('--min-contact-hours', type=int,
                        default=40, help='Minimum contact hours')
    parser.add_argument('--max-persons', type=int,
                        default=1000, help='Maximum persons to display')

    args = parser.parse_args()

    # # Debug paths (comment out for production)
    args.data_dir = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_default_2025-09-09121414"
    args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/population_structure_with_degree_distribution.png"

    location_file = args.data_dir + "/location_type_and_id.txt"
    city_config_file = args.data_dir + "/city_config.txt"

    if not os.path.exists(args.data_dir):
        print(f"Error: Directory {args.data_dir} does not exist")
        sys.exit(1)

    if not args.output_path:
        args.output_path = os.path.join(
            args.data_dir, 'population_structure_improved.png')

    # Create visualization
    fig = create_combined_visualization(
        data_dir=args.data_dir,
        location_file=location_file,
        output_path=args.output_path,
        total_runs=args.total_runs,
        min_contact_hours=args.min_contact_hours,
        max_persons=args.max_persons,
        city_config_file=city_config_file
    )

    plt.show()


if __name__ == "__main__":
    main()
