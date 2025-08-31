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

    # Build transmission trees
    print("Building transmission trees...")
    tree_1, generations_1, potential_1 = build_transmission_tree(
        infection_events_1, contact_analysis_df_1)
    tree_2, generations_2, potential_2 = build_transmission_tree(
        infection_events_2, contact_analysis_df_2)

    # Calculate positions
    positions_1 = calculate_tree_positions(infection_events_1, generations_1)
    positions_2 = calculate_tree_positions(infection_events_2, generations_2)

    # Analyze transmission opportunities
    opportunities_1 = analyze_transmission_opportunities(
        infection_events_1, contact_analysis_df_1, tree_1)
    opportunities_2 = analyze_transmission_opportunities(
        infection_events_2, contact_analysis_df_2, tree_2)

    # Create visualization
    fig = plt.figure(figsize=(20, 16))

    # Create a 3-row layout: trees on top, opportunity analysis below
    gs = fig.add_gridspec(3, 2, height_ratios=[
                          3, 3, 1], hspace=0.3, wspace=0.2)

    # Location colors
    location_colors = {
        'Work': '#4169E1',        # Royal Blue
        'Home': '#2E8B57',        # Sea Green
        'School': '#FF6347',      # Tomato Red
        'SocialEvent': '#9370DB',  # Medium Purple
        'BasicsShop': '#FF8C00',  # Dark Orange
        'EventPanvadere': '#8A2BE2',  # Blue Violet
        'Patient Zero': "#FF0000"  # Red
    }

    # Plot Method 1 tree
    ax1 = fig.add_subplot(gs[0:2, 0])
    plot_infection_tree(ax1, infection_events_1, tree_1, positions_1, location_colors,
                        f"{method1_name}\n{len(infection_events_1)} infections")

    # Plot Method 2 tree
    ax2 = fig.add_subplot(gs[0:2, 1])
    plot_infection_tree(ax2, infection_events_2, tree_2, positions_2, location_colors,
                        f"{method2_name}\n{len(infection_events_2)} infections")

    # Create opportunity analysis comparison (bottom panel spanning both columns)
    ax3 = fig.add_subplot(gs[2, :])
    plot_opportunity_comparison(
        ax3, opportunities_1, opportunities_2, method1_name, method2_name)

    # Add overall title
    title = "Comparative Infection Transmission Trees"
    if scenario_name:
        title = f"Scenario {scenario_name}: {title}"
    fig.suptitle(f"{title}\n{method1_name} vs {method2_name} Initialization",
                 fontsize=16, fontweight='bold', y=0.95)

    # Create legend (use the first subplot for positioning)
    legend_elements = []
    for loc_type, color in location_colors.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='black', lw=2, label='Transmission Chain'),
        plt.Line2D([0], [0], color='black', lw=0, marker='>', markersize=8,
                   label='Transmission Direction')
    ])

    ax1.legend(handles=legend_elements, loc='upper left', title='Infection Locations',
               title_fontsize=10, bbox_to_anchor=(0, 1))

    plt.tight_layout()

    if output_path:
        # Ensure directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=500,
                    bbox_inches='tight', facecolor='white')
        print(f"Comparative infection trees saved to {output_path}")

    # Print comparison summary
    print(f"\n=== COMPARATIVE TRANSMISSION ANALYSIS ===")
    print(f"{method1_name}:")
    print(
        f"  Realized transmissions: {opportunities_1['realized_transmissions']}")
    print(f"  Total potential contacts: {opportunities_1['total_contacts']}")
    print(
        f"  Transmission efficiency: {opportunities_1['transmission_efficiency']:.2%}")
    print(f"  Generations: {max(generations_1.values()) + 1}")

    print(f"\n{method2_name}:")
    print(
        f"  Realized transmissions: {opportunities_2['realized_transmissions']}")
    print(f"  Total potential contacts: {opportunities_2['total_contacts']}")
    print(
        f"  Transmission efficiency: {opportunities_2['transmission_efficiency']:.2%}")
    print(f"  Generations: {max(generations_2.values()) + 1}")

    return fig


def plot_infection_tree(ax, infection_events, transmission_tree, positions, location_colors, title):
    """Plot a single infection tree on the given axes"""

    times = [event['time'] for event in infection_events]

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

            # Draw L-shaped connection
            ax.plot([infector_time, infected_time], [infector_y, infector_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)
            ax.plot([infected_time, infected_time], [infector_y, infected_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)

            # Add arrow
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time, infected_y -
                                (infected_y - infector_y) * 0.1),
                        arrowprops=dict(arrowstyle='->', color='black', lw=2))

    # Plot infection points
    for event in infection_events:
        y_pos = positions[event['person_id']]
        color = location_colors.get(event['location_type'], 'gray')
        size = 300 if event['time'] == min(times) else 200

        ax.scatter(event['time'], y_pos, s=size, c=color, alpha=1.0,
                   zorder=5, edgecolors='black', linewidth=2)

        # Add person labels
        if event['time'] == min(times):  # Patient zero
            ax.annotate(f"P{event['person_id']} (P0)",
                        (event['time'], y_pos),
                        xytext=(20, 0), textcoords='offset points',
                        fontsize=9, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  alpha=0.9, edgecolor='red', linewidth=2),
                        zorder=6)
        else:
            ax.annotate(f"P{event['person_id']}",
                        (event['time'], y_pos),
                        xytext=(15, 0), textcoords='offset points',
                        fontsize=8, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                  alpha=0.9, edgecolor='gray'),
                        zorder=6)

    # Formatting
    ax.set_xlabel('Time (Simulation Timesteps)', fontsize=11)
    ax.set_ylabel('Transmission Order', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Set limits with padding
    y_values = list(positions.values())
    y_range = max(y_values) - min(y_values)
    ax.set_ylim(min(y_values) - y_range * 0.1, max(y_values) + y_range * 0.1)


def plot_opportunity_comparison(ax, opportunities_1, opportunities_2, method1_name, method2_name):
    """Plot transmission opportunity analysis comparison"""

    methods = [method1_name, method2_name]
    realized = [opportunities_1['realized_transmissions'],
                opportunities_2['realized_transmissions']]
    efficiency = [opportunities_1['transmission_efficiency'],
                  opportunities_2['transmission_efficiency']]

    # Create bar plots
    x_pos = np.arange(len(methods))
    width = 0.35

    bars1 = ax.bar(x_pos - width/2, realized, width, label='Realized Transmissions',
                   color='#2E8B57', alpha=0.8)

    # Create secondary y-axis for efficiency
    ax2 = ax.twinx()
    bars2 = ax2.bar(x_pos + width/2, [e*100 for e in efficiency], width,
                    label='Transmission Efficiency (%)', color='#FF6347', alpha=0.8)

    # Add value labels on bars
    for i, (v1, v2) in enumerate(zip(realized, efficiency)):
        ax.text(i - width/2, v1 + max(realized) * 0.01, str(v1),
                ha='center', va='bottom', fontweight='bold', fontsize=10)
        ax2.text(i + width/2, v2*100 + max(efficiency)*100 * 0.01, f'{v2:.1%}',
                 ha='center', va='bottom', fontweight='bold', fontsize=10)

    # Formatting
    ax.set_xlabel('Initialization Method', fontsize=11)
    ax.set_ylabel('Realized Transmissions', fontsize=11)
    ax2.set_ylabel('Transmission Efficiency (%)', fontsize=11)
    ax.set_title('Transmission Opportunity Analysis',
                 fontsize=12, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(methods)
    ax.grid(True, alpha=0.3, axis='y')

    # Combine legends
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')


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
