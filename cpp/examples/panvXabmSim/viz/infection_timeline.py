import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import sys


def load_simulation_data(contact_file, infection_file):
    """Load and preprocess simulation data"""
    print("Loading simulation data...")

    # Load infection data
    infection_df = pd.read_csv(infection_file)
    print(f"Loaded {len(infection_df)} infection events")

    # Load contact data (this is large, so we'll be selective)
    print("Loading contact data (this may take a moment)...")
    contact_df = pd.read_csv(contact_file)
    print(f"Loaded {len(contact_df)} contact records")

    return contact_df, infection_df


def transform_infection_data(infection_df):
    """Transform infection data to the format expected by the visualization"""
    infection_events = []

    for _, row in infection_df.iterrows():
        infection_events.append({
            'person_id': int(row['Person_ID']),
            'time': float(row['Timestep']),
            'location': str(row['Location_ID']),
            'location_type': row['Location_Type']
        })

    # First person is patient zero
    sorted_infection_events = sorted(infection_events, key=lambda x: x['time'])
    sorted_infection_events[0]['location_type'] = 'Patient Zero'
    sorted_infection_events[0]['location'] = 'Unknown'
    return sorted(infection_events, key=lambda x: x['time'])


def create_contact_data_for_analysis(contact_df, infection_df, time_window=5):
    """
    Create contact data focusing on relevant timeframes around infections
    Since the contact data is huge, we'll focus on times around infections
    """
    print("Processing contact data for infection analysis...")

    # Get all infection timesteps and create a window around them
    infection_times = set(infection_df['Timestep'])
    relevant_times = set()

    for inf_time in infection_times:
        for t in range(inf_time - time_window, inf_time + time_window + 1):
            relevant_times.add(t)

    print(f"Filtering to {len(relevant_times)} relevant timesteps...")

    # Filter contact data to relevant timesteps
    relevant_contact_df = contact_df[contact_df['Timestep'].isin(
        relevant_times)].copy()
    print(
        f"Reduced contact data from {len(contact_df)} to {len(relevant_contact_df)} records")

    # Transform to expected format
    contact_records = []

    # Group by timestep and location to find people at same place/time
    for (timestep, location), group in relevant_contact_df.groupby(['Timestep', 'Location_ID']):
        people = group['Person_ID'].tolist()

        # Create pairwise contacts for people at same location/time
        for i, person1 in enumerate(people):
            for person2 in people[i+1:]:
                contact_records.append({
                    'person_1': int(person1),
                    'person_2': int(person2),
                    'time': float(timestep),
                    'location_type': 'Unknown',  # We'll infer this
                    'location_id': str(location)
                })

    print(f"Generated {len(contact_records)} pairwise contact records")
    return pd.DataFrame(contact_records)


def has_contact_at_location(person1, person2, location, time, contact_data, time_tolerance=2):
    """Check if two people had contact at specific location around specific time"""
    relevant_contacts = contact_data[
        (((contact_data['person_1'] == person1) & (contact_data['person_2'] == person2)) |
         ((contact_data['person_1'] == person2) & (contact_data['person_2'] == person1))) &
        (contact_data['location_id'] == location) &
        (abs(contact_data['time'] - time) <= time_tolerance)
    ]
    return len(relevant_contacts) > 0


def get_contact_strength_between(person1, person2, location, contact_data):
    """Get contact strength between two specific people at location"""
    contacts = contact_data[
        (((contact_data['person_1'] == person1) & (contact_data['person_2'] == person2)) |
         ((contact_data['person_1'] == person2) & (contact_data['person_2'] == person1))) &
        (contact_data['location_id'] == location)
    ]
    return len(contacts)


def find_potential_infectors(infected_person, infection_time, location,
                             infection_events, contact_data, time_window=8):
    """Find who could have infected this person and calculate probabilities"""

    # Find people who were infectious at the same location within time window
    potential_infectors = []

    for event in infection_events:
        if (event['person_id'] != infected_person and
            event['time'] < infection_time and
                event['time'] + time_window < infection_time):

            # Check if they were at same location
            if has_contact_at_location(event['person_id'], infected_person,
                                       location, infection_time, contact_data):

                # Calculate probability based on time proximity and contact strength
                time_diff = infection_time - event['time']
                contact_strength = get_contact_strength_between(
                    event['person_id'], infected_person, location, contact_data)
                probability = 1

                if probability > 0:
                    potential_infectors.append(
                        (event['person_id'], probability))

    # Normalize probabilities
    total_prob = sum(prob for _, prob in potential_infectors)
    if total_prob > 0:
        potential_infectors = [(pid, prob/total_prob)
                               for pid, prob in potential_infectors]

    return potential_infectors


def create_timeline_with_probability_bars(infection_events, contact_data,
                                          max_infections=50, figsize=(20, 12), save_path=None):
    """
    Create timeline showing infection events with clear potential infector visualization
    Limited to first max_infections for readability
    """
    # Limit to first N infections for visualization clarity
    if len(infection_events) > max_infections:
        infection_events = infection_events[:max_infections]
        print(
            f"Limiting visualization to first {max_infections} infections for clarity")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1])

    # Location colors - updated for simulation data
    location_colors = {
        'Work': '#4169E1',        # Royal Blue
        'Home': '#2E8B57',        # Sea Green
        'School': '#FF6347',      # Tomato Red
        'SocialEvent': '#9370DB',  # Medium Purple
        'BasicsShop': '#FF8C00',  # Dark Orange
        'EventPanvadere': '#8A2BE2',  # Blue Violet
        'Patient Zero': "#FF0000"  # Gray
    }

    # Main timeline (top)
    times = [event['time'] for event in infection_events]
    persons = [event['person_id'] for event in infection_events]

    # Plot infection points with location-based colors
    for event in infection_events:
        color = location_colors.get(event['location_type'], 'gray')
        size = 200 if event['time'] == min(
            times) else 120  # Larger for first infection
        alpha = 1.0 if event['time'] == min(times) else 0.8

        ax1.scatter(event['time'], event['person_id'], s=size, c=color,
                    alpha=alpha, zorder=3, edgecolors='black', linewidth=2)

    # Add person labels
    for i, event in enumerate(infection_events):
        if event['time'] < 0.0:
            ax1.annotate(f"{event['location_type']}",
                         (event['time'], event['person_id']),
                         xytext=(12, 1), textcoords='offset points',
                         fontsize=8, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.3',
                                   facecolor='white', alpha=0.9, edgecolor='gray'))

    # Draw transmission probability analysis
    for i, event in enumerate(infection_events[1:], 1):  # Skip first infection
        infected_person = event['person_id']
        potential_infectors = find_potential_infectors(
            infected_person, event['time'], event['location'],
            # Only consider earlier infections
            infection_events[:i], contact_data)

        if len(potential_infectors) > 1:
            # Multiple infectors - draw direct arrows from each infector to infected person
            for j, (infector, probability) in enumerate(potential_infectors):
                # Find infector position
                infector_time = None
                for e in infection_events:
                    if e['person_id'] == infector:
                        infector_time = e['time']
                        break

                if infector_time:
                    # Calculate arrow width based on probability (thinner for uncertain)
                    arrow_width = 0.8 + 1.2 * probability  # Range: 0.8 to 2.0
                    alpha_intensity = 0.4 + 0.4 * probability  # Range: 0.4 to 0.8

                    # Use orange color for uncertain transmissions
                    arrow_color = 'orange'

                    # Add slight curve to avoid overlapping arrows
                    connection_style = f"arc3,rad={0.1 + j * 0.05}"

                    # Draw direct arrow from infector to infected person
                    ax1.annotate('', xy=(event['time'], infected_person),
                                 xytext=(infector_time, infector),
                                 arrowprops=dict(arrowstyle='->',
                                                 color=arrow_color,
                                                 alpha=alpha_intensity,
                                                 lw=arrow_width,
                                                 connectionstyle=connection_style))

                    # Add probability label near the arrow
                    mid_x = (infector_time + event['time']) / 2
                    mid_y = (infector + infected_person) / \
                        2 + j * max(persons) * 0.015
                    ax1.text(mid_x, mid_y,
                             f'{probability:.0%}',
                             ha='center', va='center', fontsize=6, fontweight='bold',
                             bbox=dict(boxstyle='round,pad=0.2',
                                       facecolor='orange', alpha=0.6, edgecolor='darkorange'))

        elif len(potential_infectors) == 1:
            # Single infector - draw simple arrow
            infector, probability = potential_infectors[0]

            # Find infector position
            infector_time = None
            for e in infection_events:
                if e['person_id'] == infector:
                    infector_time = e['time']
                    break

            if infector_time:
                # Draw direct arrow with thicker line for certain transmission
                ax1.annotate('', xy=(event['time'], infected_person),
                             xytext=(infector_time, infector),
                             arrowprops=dict(arrowstyle='->',
                                             color='darkgreen',
                                             alpha=0.8,
                                             lw=3,
                                             connectionstyle="arc3,rad=0.1"))

                # Add certainty label
                mid_x = (infector_time + event['time']) / 2
                mid_y = (infector + infected_person) / 2
                ax1.text(mid_x, mid_y + max(persons) * 0.02,
                         f'P{infector}→P{infected_person}',
                         ha='center', va='center', fontsize=6, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.2',
                                   facecolor='lightgreen', alpha=0.8))

    ax1.set_xlabel('Time (Simulation Timesteps)', fontsize=12)
    ax1.set_ylabel('Person ID', fontsize=12)
    ax1.set_title(f'Infection Timeline - Simulation Data\n' +
                  f'First {len(infection_events)} infections (Timesteps {min(times)}-{max(times)})',
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Set y-axis limits with some padding
    person_range = max(persons) - min(persons)
    ax1.set_ylim(min(persons) - person_range * 0.1,
                 max(persons) + person_range * 0.1)

    # Create legend
    legend_elements = []
    for loc_type, color in location_colors.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='darkgreen', lw=3,
                   label='Certain Transmission'),
        plt.Line2D([0], [0], color='orange', lw=2,
                   label='Uncertain Transmission')
    ])

    ax1.legend(handles=legend_elements, loc='upper left',
               title='Legend', title_fontsize=10)

    # Uncertainty summary (bottom)
    uncertainty_scores = []
    times_uncertain = []

    for i, event in enumerate(infection_events[1:], 1):
        potential_infectors = find_potential_infectors(
            event['person_id'], event['time'], event['location'],
            infection_events[:i], contact_data)

        if potential_infectors:
            probs = [prob for _, prob in potential_infectors]
            if len(probs) > 1:
                entropy = -sum(p * np.log2(p) if p > 0 else 0 for p in probs)
            else:
                entropy = 0
            uncertainty_scores.append(entropy)
            times_uncertain.append(event['time'])

    if uncertainty_scores:
        ax2.plot(times_uncertain, uncertainty_scores, 'o-', color='red',
                 linewidth=2, markersize=6, alpha=0.7)
        ax2.fill_between(times_uncertain, uncertainty_scores,
                         alpha=0.3, color='red')

        if max(uncertainty_scores) > 0:
            high_uncertainty_threshold = max(uncertainty_scores) * 0.7
            ax2.axhline(y=high_uncertainty_threshold, color='orange',
                        linestyle='--', alpha=0.7,
                        label=f'High Uncertainty (>{high_uncertainty_threshold:.1f})')
            ax2.legend()

    ax2.set_xlabel('Time (Simulation Timesteps)', fontsize=12)
    ax2.set_ylabel('Uncertainty\n(Entropy)', fontsize=12)
    ax2.set_title('Transmission Uncertainty Over Time',
                  fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    if uncertainty_scores:
        ax2.set_ylim(0, max(uncertainty_scores) * 1.1)

    # Add summary statistics
    if uncertainty_scores:
        avg_uncertainty = np.mean(uncertainty_scores)
        max_uncertainty = max(uncertainty_scores)
        high_uncertainty_count = sum(
            1 for u in uncertainty_scores if u > max_uncertainty * 0.7)

        stats_text = f"""UNCERTAINTY SUMMARY:
                        Average: {avg_uncertainty:.2f}
                        Maximum: {max_uncertainty:.2f}
                        High uncertainty events: {high_uncertainty_count}"""

        ax2.text(0.98, 0.98, stats_text, transform=ax2.transAxes,
                 fontsize=10, ha='right', va='top',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Infection timeline saved to {save_path}")


def analyze_simulation_infections(contact_file, infection_file, max_display=50, save_path=None, scenario_name="Simulation"):
    """Main function to analyze and visualize simulation infections"""

    # Load data
    contact_df, infection_df = load_simulation_data(
        contact_file, infection_file)

    # Transform infection data
    infection_events = transform_infection_data(infection_df)

    # Create relevant contact data
    contact_analysis_df = create_contact_data_for_analysis(
        contact_df, infection_df)

    # Print summary
    print(f"\n=== {scenario_name.upper()} ANALYSIS SUMMARY ===")
    print(f"Total infections: {len(infection_events)}")
    print(
        f"Timeline: Timesteps {infection_events[0]['time']:.0f} to {infection_events[-1]['time']:.0f}")
    print(
        f"Duration: {infection_events[-1]['time'] - infection_events[0]['time']:.0f} timesteps")

    # Location type breakdown
    location_counts = {}
    for event in infection_events:
        loc_type = event['location_type']
        location_counts[loc_type] = location_counts.get(loc_type, 0) + 1

    print(f"\nInfections by location type:")
    for loc_type, count in sorted(location_counts.items(), key=lambda x: -x[1]):
        percentage = count / len(infection_events) * 100
        print(f"  {loc_type}: {count} ({percentage:.1f}%)")

    # Create visualization
    print(f"\nCreating visualization...")
    create_timeline_with_probability_bars(
        infection_events, contact_analysis_df,
        max_infections=max_display, save_path=save_path
    )

    if save_path:
        print(f"Visualization saved as '{save_path}'")

    return infection_events, contact_analysis_df


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description='Analyze infection timeline from simulation data')
    parser.add_argument('--contact-file',
                        help='Path to contact data CSV file')
    parser.add_argument('--infection-file',
                        help='Path to infection data CSV file')
    parser.add_argument(
        '--output-path', help='Output path for visualization (optional)')
    parser.add_argument('--scenario-name', default='Simulation',
                        help='Scenario name for titles')
    parser.add_argument('--max-display', type=int, default=50,
                        help='Maximum infections to display')

    args = parser.parse_args()

    # for debugging hardcode filepath
    args.contact_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/best_run_contact_data.csv"
    args.infection_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/best_run_detailed_infection.csv"
    args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/infection_timeline_simulation.png"

    # Check if files exist
    if not os.path.exists(args.contact_file):
        print(f"Error: Contact file {args.contact_file} does not exist")
        sys.exit(1)

    if not os.path.exists(args.infection_file):
        print(f"Error: Infection file {args.infection_file} does not exist")
        sys.exit(1)

    # Set output path
    if args.output_path:
        save_path = args.output_path
    else:
        base_dir = os.path.dirname(args.contact_file)
        save_path = os.path.join(
            base_dir, f'infection_timeline_{args.scenario_name.lower()}.png')

    # Analyze the simulation data
    infection_events, contact_data = analyze_simulation_infections(
        args.contact_file,
        args.infection_file,
        max_display=args.max_display,
        save_path=save_path,
        scenario_name=args.scenario_name
    )


if __name__ == "__main__":
    main()
