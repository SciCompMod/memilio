import pandas as pd
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
                             infection_events, contact_data=None, time_window=24, infectiousness_delay=4.5):
    """Find who could have infected this person based on location and timing
    
    Args:
        infected_person: ID of person who got infected
        infection_time: Time when infection occurred
        location: Location where infection occurred
        infection_events: List of all infection events
        contact_data: Not used anymore - presence at same location implies contact
        time_window: Maximum time window for potential transmission (hours)
        infectiousness_delay: Days after infection when person becomes infectious
    """

    # Find people who were infectious at the same location within time window
    potential_infectors = []

    for event in infection_events:
        if (event['person_id'] != infected_person and
            event['time'] < infection_time):

            # Check if they became infectious before the infection occurred
            # Person becomes infectious after infectiousness_delay days (converted to hours)
            infectious_time = event['time'] + (infectiousness_delay * 24)  # Convert days to hours
            
            # They must be infectious and the infection must occur within the time window
            if (infectious_time <= infection_time and 
                infection_time <= infectious_time + time_window):

                # Check if they were at same location type (same location implies contact)
                if event['location_type'] == location.get('type', location) if isinstance(location, dict) else event['location_type'] == location:

                    # Calculate probability based on how long they were infectious
                    # Longer infectious period = higher probability
                    days_infectious = (infection_time - infectious_time) / 24
                    probability = min(1.0, days_infectious / 7.0)  # Max probability after 7 days infectious

                    if probability > 0:
                        potential_infectors.append(
                            (event['person_id'], probability))

    # Normalize probabilities
    total_prob = sum(prob for _, prob in potential_infectors)
    if total_prob > 0:
        potential_infectors = [(pid, prob/total_prob)
                               for pid, prob in potential_infectors]

    return potential_infectors


def build_transmission_tree(infection_events, contact_data=None, show_all_potential_infectors=False):
    """Build a transmission tree to organize the layout hierarchically
    
    Args:
        infection_events: List of infection events with time, location_type, person_id
        contact_data: Not used anymore - same location implies contact
        show_all_potential_infectors: Whether to show all potential transmission paths
    """
    print("Building transmission tree...")

    # Build transmission relationships
    transmission_tree = {}
    infection_generations = {}
    all_potential_transmissions = {}  # Store all potential infectors for visualization

    # Patient zero is generation 0
    patient_zero = min(infection_events, key=lambda x: x['time'])
    transmission_tree[patient_zero['person_id']] = []
    infection_generations[patient_zero['person_id']] = 0

    # For each subsequent infection, find most likely infector
    for i, event in enumerate(infection_events[1:], 1):
        infected_person = event['person_id']
        # Pass location_type instead of full location object
        location = event['location_type']
        potential_infectors = find_potential_infectors(
            infected_person, event['time'], location,
            infection_events[:i], contact_data)

        if potential_infectors:
            # Store all potential infectors for optional visualization
            all_potential_transmissions[infected_person] = potential_infectors

            # Choose most likely infector for primary tree structure
            most_likely_infector = max(
                potential_infectors, key=lambda x: x[1])[0]

            # Add to transmission tree
            if most_likely_infector not in transmission_tree:
                transmission_tree[most_likely_infector] = []
            transmission_tree[most_likely_infector].append(infected_person)

            # Set generation (one more than infector)
            if most_likely_infector in infection_generations:
                infection_generations[infected_person] = infection_generations[most_likely_infector] + 1
            else:
                infection_generations[infected_person] = 1
        else:
            # Orphan infection - assign to generation 1
            infection_generations[infected_person] = 1

    return transmission_tree, infection_generations, all_potential_transmissions


def calculate_tree_positions(infection_events, infection_generations):
    """Calculate Y positions based on transmission tree structure with proper tree layout"""
    print("Calculating tree-based positions...")

    # Group people by generation
    generations = {}
    for person_id, generation in infection_generations.items():
        if generation not in generations:
            generations[generation] = []
        generations[generation].append(person_id)

    # Create tree positions with proper spacing
    person_positions = {}
    current_y = 0

    for gen in sorted(generations.keys()):
        people_in_gen = generations[gen]

        # Get infection times for sorting within generation
        people_with_times = []
        for person_id in people_in_gen:
            event = next(
                e for e in infection_events if e['person_id'] == person_id)
            people_with_times.append((person_id, event['time']))

        # Sort by time within generation
        people_with_times.sort(key=lambda x: x[1])

        # Assign positions with consistent spacing
        for i, (person_id, _) in enumerate(people_with_times):
            person_positions[person_id] = current_y + i

        # Move to next generation with spacing
        current_y += len(people_in_gen) + 2  # Add spacing between generations

    return person_positions


def create_timeline_with_probability_bars(infection_events, contact_data,
                                          max_infections=50, figsize=(20, 12), save_path=None, show_all_potential_infectors=False):
    """
    Create timeline showing infection events with tree-based layout and clear transmission lines
    Limited to first max_infections for readability
    """
    # Limit to first N infections for visualization clarity
    if len(infection_events) > max_infections:
        infection_events = infection_events[:max_infections]
        print(
            f"Limiting visualization to first {max_infections} infections for clarity")

    # Build transmission tree for better layout (no longer using contact data)
    transmission_tree, infection_generations, all_potential_transmissions = build_transmission_tree(
        infection_events, None, show_all_potential_infectors)
    person_positions = calculate_tree_positions(
        infection_events, infection_generations)

    # Create single plot instead of subplots
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Location colors - updated for simulation data
    location_colors = {
        'Work': '#4169E1',        # Royal Blue
        'Home': '#2E8B57',        # Sea Green
        'School': '#FF6347',      # Tomato Red
        'SocialEvent': '#9370DB',  # Medium Purple
        'BasicsShop': '#FF8C00',  # Dark Orange
        'EventPanvadere': '#8A2BE2',  # Blue Violet
        'Patient Zero': "#FF0000"  # Red
    }

    # Main timeline (top) - using tree-based positioning
    times = [event['time'] for event in infection_events]

    # Draw tree structure with straight vertical and horizontal lines
    # First, remove the old spine drawing

    # Draw the tree connections using straight lines only
    for infector_id, infected_list in transmission_tree.items():
        if not infected_list:
            continue

        # Find infector event and position
        infector_event = next(
            e for e in infection_events if e['person_id'] == infector_id)
        infector_y = person_positions[infector_id]
        infector_time = infector_event['time']

        if len(infected_list) == 1:
            # Single infection - use L-shaped path (horizontal then vertical)
            infected_id = infected_list[0]
            infected_event = next(
                e for e in infection_events if e['person_id'] == infected_id)
            infected_y = person_positions[infected_id]
            infected_time = infected_event['time']

            # Draw L-shaped connection: horizontal then vertical
            # Step 1: Horizontal line from infector to infected person's time
            ax.plot([infector_time, infected_time], [infector_y, infector_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)

            # Step 2: Vertical line from infector's Y to infected person's Y
            ax.plot([infected_time, infected_time], [infector_y, infected_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)

            # Add small arrow at the end
            ax.annotate('', xy=(infected_time, infected_y),
                        xytext=(infected_time, infected_y -
                                (infected_y - infector_y) * 0.1),
                        arrowprops=dict(arrowstyle='->', color='black', lw=2))

        else:
            # Multiple infections - use tree structure with vertical and horizontal lines

            # Find the horizontal position for the vertical line (between infector and earliest infected)
            earliest_infected_time = min(
                next(e for e in infection_events if e['person_id'] == iid)[
                    'time']
                for iid in infected_list
            )
            vertical_line_x = infector_time + \
                (earliest_infected_time - infector_time) * 0.3

            # Draw horizontal line from infector to vertical line position
            ax.plot([infector_time, vertical_line_x], [infector_y, infector_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)

            # Find Y range for vertical line
            infected_y_positions = [person_positions[iid]
                                    for iid in infected_list]
            min_y = min(infected_y_positions + [infector_y])
            max_y = max(infected_y_positions + [infector_y])

            # Draw vertical line connecting all infected people
            ax.plot([vertical_line_x, vertical_line_x], [min_y, max_y],
                    color='black', linewidth=2, alpha=0.8, zorder=2)

            # Draw horizontal lines from vertical line to each infected person
            for infected_id in infected_list:
                infected_event = next(
                    e for e in infection_events if e['person_id'] == infected_id)
                infected_y = person_positions[infected_id]
                infected_time = infected_event['time']

                # Horizontal line from vertical line to infected person
                ax.plot([vertical_line_x, infected_time], [infected_y, infected_y],
                        color='black', linewidth=2, alpha=0.8, zorder=2)

                # Small arrow at the end
                ax.annotate('', xy=(infected_time, infected_y),
                            xytext=(infected_time - (infected_time -
                                                     vertical_line_x) * 0.1, infected_y),
                            arrowprops=dict(arrowstyle='->', color='black', lw=2))

    # Draw all potential infectors if requested (with less saturated colors)
    if show_all_potential_infectors:
        for infected_person, potential_infectors in all_potential_transmissions.items():
            if len(potential_infectors) > 1:  # Only show if there are multiple options
                # Find infected person's position and time
                infected_event = next(
                    e for e in infection_events if e['person_id'] == infected_person)
                infected_y = person_positions[infected_person]
                infected_time = infected_event['time']

                # Draw lines to all potential infectors (except the most likely one)
                most_likely_infector = max(
                    potential_infectors, key=lambda x: x[1])[0]

                for infector_id, probability in potential_infectors:
                    # Skip the most likely one (already drawn)
                    if infector_id != most_likely_infector:
                        # Find infector position
                        infector_event = next(
                            e for e in infection_events if e['person_id'] == infector_id)
                        infector_y = person_positions[infector_id]
                        infector_time = infector_event['time']

                        # Draw with more prominent styling for better readability
                        # Increased base alpha for better visibility
                        alpha_val = 0.6 + (probability * 0.3)
                        line_color = 'darkorange'  # More visible color than lightgray
                        line_width = 2.5  # Thicker lines for better readability

                        # Use L-shaped path like the main tree structure
                        # Step 1: Horizontal line from infector towards infected person's time
                        intermediate_x = infector_time + \
                            (infected_time - infector_time) * 0.7
                        ax.plot([infector_time, intermediate_x], [infector_y, infector_y],
                                color=line_color, linewidth=line_width, alpha=alpha_val,
                                linestyle='--', zorder=2)  # Dashed instead of dotted, higher zorder

                        # Step 2: Vertical line from infector's Y to infected person's Y
                        ax.plot([intermediate_x, intermediate_x], [infector_y, infected_y],
                                color=line_color, linewidth=line_width, alpha=alpha_val,
                                linestyle='--', zorder=2)

                        # Step 3: Short horizontal line to infected person
                        ax.plot([intermediate_x, infected_time], [infected_y, infected_y],
                                color=line_color, linewidth=line_width, alpha=alpha_val,
                                linestyle='--', zorder=2)

                        # End with a more prominent dot to show uncertainty
                        ax.scatter(infected_time, infected_y, s=50, c='red',
                                   alpha=alpha_val, zorder=5, marker='o', edgecolors='darkred', linewidth=2)

                        # Add probability label for uncertain transmissions with better styling
                        mid_x = (infector_time + infected_time) / 2
                        mid_y = (infector_y + infected_y) / 2
                        ax.text(mid_x, mid_y, f'{probability:.0%}',
                                ha='center', va='center', fontsize=7, fontweight='bold', alpha=0.9,
                                bbox=dict(boxstyle='round,pad=0.2',
                                          facecolor='yellow', alpha=0.8, edgecolor='darkorange'),
                                zorder=6)

    # Plot infection points AFTER drawing the tree lines so they appear on top
    for event in infection_events:
        y_pos = person_positions[event['person_id']]
        color = location_colors.get(event['location_type'], 'gray')
        size = 300 if event['time'] == min(
            times) else 200  # Larger for patient zero
        alpha = 1.0

        # Plot the infection point
        ax.scatter(event['time'], y_pos, s=size, c=color,
                   alpha=alpha, zorder=5, edgecolors='black', linewidth=2)

    # Add person labels for key events
    for event in infection_events:
        y_pos = person_positions[event['person_id']]
        if event['time'] == min(times):  # Patient zero
            ax.annotate(f"P{event['person_id']} (Patient Zero)",
                        (event['time'], y_pos),
                        xytext=(20, 0), textcoords='offset points',
                        fontsize=10, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.4',
                                  facecolor='white', alpha=0.9, edgecolor='red', linewidth=2),
                        zorder=6)
        else:
            # Add person ID labels
            ax.annotate(f"P{event['person_id']}",
                        (event['time'], y_pos),
                        xytext=(15, 0), textcoords='offset points',
                        fontsize=9, fontweight='bold', alpha=1.0,
                        bbox=dict(boxstyle='round,pad=0.2',
                                  facecolor='white', alpha=0.9, edgecolor='gray'),
                        zorder=6)

    ax.set_xlabel('Time (Simulation Timesteps)', fontsize=12)
    ax.set_ylabel('Transmission Order', fontsize=12)
    ax.set_title('Infection Transmission Tree\n' +
                 f'First {len(infection_events)} infections (Timesteps {min(times):.0f}-{max(times):.0f})',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Set y-axis limits with some padding based on tree positions
    y_positions_values = list(person_positions.values())
    y_range = max(y_positions_values) - min(y_positions_values)
    ax.set_ylim(min(y_positions_values) - y_range * 0.1,
                max(y_positions_values) + y_range * 0.1)

    # Create legend
    legend_elements = []
    for loc_type, color in location_colors.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          label=loc_type))

    legend_elements.extend([
        plt.Line2D([0], [0], color='black', lw=2,
                   label='Certain Transmission'),
        plt.Line2D([0], [0], color='black', lw=0, marker='>', markersize=8,
                   label='Transmission Direction')
    ])

    if show_all_potential_infectors:
        legend_elements.extend([
            plt.Line2D([0], [0], color='darkorange', lw=2.5, linestyle='--',
                       label='Uncertain Transmission'),
            plt.Line2D([0], [0], color='red', lw=0, marker='o', markersize=8,
                       label='Uncertainty Indicator')
        ])

    ax.legend(handles=legend_elements, loc='upper left',
              title='Legend', title_fontsize=10)

    plt.tight_layout()

    if save_path:
        # Ensure the directory exists before saving
        save_dir = os.path.dirname(save_path)
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir)

        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Infection timeline saved to {save_path}")


def analyze_simulation_infections(contact_file, infection_file, max_display=50, save_path=None, scenario_name="Simulation", show_all_potential_infectors=False):
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

    print("\nInfections by location type:")
    for loc_type, count in sorted(location_counts.items(), key=lambda x: -x[1]):
        percentage = count / len(infection_events) * 100
        print(f"  {loc_type}: {count} ({percentage:.1f}%)")

    # Create visualization
    print("\nCreating visualization...")
    create_timeline_with_probability_bars(
        infection_events, contact_analysis_df,
        max_infections=max_display, save_path=save_path,
        show_all_potential_infectors=show_all_potential_infectors
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
    parser.add_argument('--show-all-potential-infectors', action="store_true",
                        help='Show all potential infectors with less saturated colors')
    parser.add_argument('--display-one-infection', action="store_true",
                        help='Whether to display only one infection')

    args = parser.parse_args()

    # for debugging hardcode filepath
    # args.contact_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/best_run_contact_data.csv"
    # args.infection_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/best_run_detailed_infection.csv"
    # args.output_path = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/run_20250813_144038_memilio/infection_timeline_simulation.png"

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
    analyze_simulation_infections(
        args.contact_file,
        args.infection_file,
        max_display=args.max_display,
        save_path=save_path,
        scenario_name=args.scenario_name,
        show_all_potential_infectors=args.show_all_potential_infectors
    )


if __name__ == "__main__":
    main()
