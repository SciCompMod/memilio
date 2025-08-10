import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
                             infection_events, contact_data, time_window=24):
    """Find who could have infected this person and calculate probabilities"""

    # Find people who were infectious at the same location within time window
    potential_infectors = []

    for event in infection_events:
        if (event['person_id'] != infected_person and
            event['time'] < infection_time and
                event['time'] > infection_time - time_window):

            # Check if they were at same location
            if has_contact_at_location(event['person_id'], infected_person,
                                       location, infection_time, contact_data):

                # Calculate probability based on time proximity and contact strength
                time_diff = infection_time - event['time']
                contact_strength = get_contact_strength_between(
                    event['person_id'], infected_person, location, contact_data)

                # Simple probability model
                time_factor = max(0, 1 - time_diff / time_window)
                # Normalize contact strength
                contact_factor = min(contact_strength / 3, 1.0)
                probability = time_factor * contact_factor

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
                                          figsize=(20, 12), save_path=None):
    """
    Create timeline showing infection events with clear potential infector visualization
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, height_ratios=[3, 1])

    # Location colors
    location_colors = {
        'Work': '#4169E1',        # Royal Blue
        'Home': '#2E8B57',        # Sea Green
        'School': '#FF6347',      # Tomato Red
        'Social_Event': '#9370DB',  # Medium Purple
        'Shopping': '#FF8C00'     # Dark Orange
    }

    # Sort by time
    infection_events = sorted(infection_events, key=lambda x: x['time'])

    # Main timeline (top)
    times = [event['time'] for event in infection_events]
    persons = [event['person_id'] for event in infection_events]

    # Plot infection points with location-based colors
    for event in infection_events:
        color = location_colors.get(event['location_type'], 'gray')
        size = 200 if event['time'] == 0 else 120  # Larger for patient zero
        alpha = 1.0 if event['time'] == 0 else 0.8

        ax1.scatter(event['time'], event['person_id'], s=size, c=color,
                    alpha=alpha, zorder=3, edgecolors='black', linewidth=2)

    # Add person labels
    for event in infection_events:
        label_text = f"P{event['person_id']}"
        if event['time'] == 0:
            label_text = f"P{event['person_id']}\n(Patient Zero)"

        ax1.annotate(f"{label_text}\n{event['location']}",
                     (event['time'], event['person_id']),
                     xytext=(10, 0), textcoords='offset points',
                     fontsize=9, fontweight='bold',
                     bbox=dict(boxstyle='round,pad=0.3',
                               facecolor='white', alpha=0.9, edgecolor='gray'))

    # Draw enhanced probability visualization
    for event in infection_events[1:]:  # Skip patient zero
        infected_person = event['person_id']
        potential_infectors = find_potential_infectors(
            infected_person, event['time'], event['location'],
            infection_events, contact_data)

        if len(potential_infectors) > 1:
            # Create a grouped bar chart for potential infectors
            bar_width = 3.0
            bar_spacing = 0.5
            total_width = len(potential_infectors) * (bar_width + bar_spacing)
            start_x = event['time'] - total_width - 2

            # Draw connection lines from each potential infector to infected person
            for j, (infector, probability) in enumerate(potential_infectors):
                bar_x = start_x + j * (bar_width + bar_spacing)
                bar_y = infected_person - 0.4
                bar_height = probability * 0.8

                # Color bar based on probability (gradient)
                alpha_intensity = 0.5 + 0.5 * probability  # 50-100% alpha based on probability
                bar_color = location_colors.get(event['location_type'], 'gray')

                # Draw probability bar
                bar = ax1.barh(bar_y, bar_width, height=bar_height,
                               left=bar_x, alpha=alpha_intensity,
                               color=bar_color, edgecolor='black', linewidth=1)

                # Add infector label on the bar
                ax1.text(bar_x + bar_width/2, bar_y + bar_height/2,
                         f'P{infector}\n{probability:.0%}',
                         ha='center', va='center', fontsize=8, fontweight='bold',
                         color='white')

                # Draw arrow from bar to infected person
                arrow_start = (bar_x + bar_width, bar_y + bar_height/2)
                arrow_end = (event['time'] - 1, infected_person)

                # Arrow thickness based on probability
                arrow_width = 1 + 3 * probability

                ax1.annotate('', xy=arrow_end, xytext=arrow_start,
                             arrowprops=dict(arrowstyle='->',
                                             color=bar_color,
                                             alpha=alpha_intensity,
                                             lw=arrow_width))

            # Add uncertainty summary box
            infector_list = ", ".join(
                [f"P{pid}({prob:.0%})" for pid, prob in potential_infectors])
            uncertainty_text = f"Potential infectors:\n{infector_list}"

            ax1.text(start_x + total_width/2, infected_person + 0.6, uncertainty_text,
                     ha='center', va='bottom', fontsize=7,
                     bbox=dict(boxstyle='round,pad=0.3',
                               facecolor='yellow', alpha=0.7, edgecolor='orange'))

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
                # Draw direct arrow
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
                ax1.text(mid_x, mid_y + 0.3, f'P{infector}→P{infected_person}\n(Certain)',
                         ha='center', va='center', fontsize=7, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.2',
                                   facecolor='lightgreen', alpha=0.8))

    ax1.set_xlabel('Time (hours)', fontsize=12)
    ax1.set_ylabel('Person ID', fontsize=12)
    ax1.set_title('Infection Timeline with Clear Transmission Attribution\n' +
                  '(Bars show potential infectors with probabilities, arrows show transmission paths)',
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(persons) + 2)  # Extra space for labels

    # Create legend
    legend_elements = []
    for loc_type, color in location_colors.items():
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w',
                                          markerfacecolor=color, markersize=10,
                                          label=loc_type.replace('_', ' ')))

    # Add transmission type legend
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

    for event in infection_events[1:]:
        potential_infectors = find_potential_infectors(
            event['person_id'], event['time'], event['location'],
            infection_events, contact_data)

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
                 linewidth=3, markersize=8, alpha=0.7)
        ax2.fill_between(times_uncertain, uncertainty_scores,
                         alpha=0.3, color='red')

        # Add uncertainty threshold line
        if max(uncertainty_scores) > 0:
            high_uncertainty_threshold = max(uncertainty_scores) * 0.7
            ax2.axhline(y=high_uncertainty_threshold, color='orange',
                        linestyle='--', alpha=0.7,
                        label=f'High Uncertainty (>{high_uncertainty_threshold:.1f})')
            ax2.legend()

    ax2.set_xlabel('Time (hours)', fontsize=12)
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


# Sample data for testing
sample_infection_events = [
    {'person_id': 1, 'time': 0.0, 'location': 'W1', 'location_type': 'Work'},
    {'person_id': 4, 'time': 8.5, 'location': 'W1', 'location_type': 'Work'},
    {'person_id': 7, 'time': 9.2, 'location': 'W1', 'location_type': 'Work'},
    {'person_id': 2, 'time': 14.3, 'location': 'H1', 'location_type': 'Home'},
    {'person_id': 3, 'time': 15.1, 'location': 'H1', 'location_type': 'Home'},
    {'person_id': 5, 'time': 18.7, 'location': 'H2', 'location_type': 'Home'},
    {'person_id': 8, 'time': 24.2, 'location': 'E1',
        'location_type': 'Social_Event'},
    {'person_id': 9, 'time': 24.8, 'location': 'E1',
        'location_type': 'Social_Event'},
    {'person_id': 6, 'time': 25.5, 'location': 'H2', 'location_type': 'Home'},
    {'person_id': 10, 'time': 32.1, 'location': 'S1', 'location_type': 'School'},
    {'person_id': 11, 'time': 33.4, 'location': 'S1', 'location_type': 'School'},
    {'person_id': 12, 'time': 36.7, 'location': 'SH1', 'location_type': 'Shopping'},
    {'person_id': 13, 'time': 37.2, 'location': 'SH1', 'location_type': 'Shopping'},
]

sample_contact_data = [
    {'person_1': 1, 'person_2': 4, 'time': 8.0,
        'location_type': 'Work', 'location_id': 'W1'},
    {'person_1': 1, 'person_2': 7, 'time': 8.5,
        'location_type': 'Work', 'location_id': 'W1'},
    {'person_1': 4, 'person_2': 7, 'time': 9.0,
        'location_type': 'Work', 'location_id': 'W1'},
    {'person_1': 1, 'person_2': 2, 'time': 14.0,
        'location_type': 'Home', 'location_id': 'H1'},
    {'person_1': 1, 'person_2': 3, 'time': 14.5,
        'location_type': 'Home', 'location_id': 'H1'},
    {'person_1': 2, 'person_2': 3, 'time': 15.0,
        'location_type': 'Home', 'location_id': 'H1'},
    {'person_1': 4, 'person_2': 5, 'time': 18.0,
        'location_type': 'Home', 'location_id': 'H2'},
    {'person_1': 4, 'person_2': 6, 'time': 24.0,
        'location_type': 'Home', 'location_id': 'H2'},
    {'person_1': 5, 'person_2': 6, 'time': 25.0,
        'location_type': 'Home', 'location_id': 'H2'},
    {'person_1': 1, 'person_2': 8, 'time': 23.5,
        'location_type': 'Social_Event', 'location_id': 'E1'},
    {'person_1': 4, 'person_2': 8, 'time': 23.8,
        'location_type': 'Social_Event', 'location_id': 'E1'},
    {'person_1': 7, 'person_2': 8, 'time': 24.0,
        'location_type': 'Social_Event', 'location_id': 'E1'},
    {'person_1': 8, 'person_2': 9, 'time': 24.5,
        'location_type': 'Social_Event', 'location_id': 'E1'},
    {'person_1': 3, 'person_2': 10, 'time': 31.5,
        'location_type': 'School', 'location_id': 'S1'},
    {'person_1': 10, 'person_2': 11, 'time': 33.0,
        'location_type': 'School', 'location_id': 'S1'},
    {'person_1': 2, 'person_2': 12, 'time': 36.0,
        'location_type': 'Shopping', 'location_id': 'SH1'},
    {'person_1': 5, 'person_2': 12, 'time': 36.3,
        'location_type': 'Shopping', 'location_id': 'SH1'},
    {'person_1': 8, 'person_2': 12, 'time': 36.5,
        'location_type': 'Shopping', 'location_id': 'SH1'},
    {'person_1': 12, 'person_2': 13, 'time': 37.0,
        'location_type': 'Shopping', 'location_id': 'SH1'},
]


def run_example():
    """Run the timeline visualization with sample data"""

    print("Creating enhanced timeline with clear transmission attribution...")

    # Convert to DataFrame
    contact_df = pd.DataFrame(sample_contact_data)

    # Create the visualization
    create_timeline_with_probability_bars(
        sample_infection_events,
        contact_df,
        figsize=(22, 12),
        save_path='infection_timeline_enhanced.png'
    )

    print("Visualization saved as 'infection_timeline_enhanced.png'")

    # Print summary
    print(f"\nSample scenario summary:")
    print(f"- {len(sample_infection_events)} infection events")
    print(f"- {len(sample_contact_data)} contact records")
    print(
        f"- Timeline spans {max(e['time'] for e in sample_infection_events):.1f} hours")

    # Analyze uncertainty cases
    contact_df = pd.DataFrame(sample_contact_data)
    uncertain_cases = 0
    certain_cases = 0

    for event in sample_infection_events[1:]:
        potential_infectors = find_potential_infectors(
            event['person_id'], event['time'], event['location'],
            sample_infection_events, contact_df)

        if len(potential_infectors) > 1:
            uncertain_cases += 1
            infector_info = ", ".join(
                [f"P{pid}({prob:.0%})" for pid, prob in potential_infectors])
            print(
                f"  UNCERTAIN - Person {event['person_id']}: {infector_info}")
        elif len(potential_infectors) == 1:
            certain_cases += 1
            infector, prob = potential_infectors[0]
            print(
                f"  CERTAIN   - Person {event['person_id']}: P{infector} (100%)")

    print(f"\nTransmission certainty:")
    print(f"- Certain transmissions: {certain_cases}")
    print(f"- Uncertain transmissions: {uncertain_cases}")
    print(
        f"- Uncertainty rate: {uncertain_cases/(certain_cases+uncertain_cases)*100:.1f}%")


if __name__ == "__main__":
    run_example()
