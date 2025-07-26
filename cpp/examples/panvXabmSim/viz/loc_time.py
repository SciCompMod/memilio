import matplotlib.pyplot as plt
import numpy as np


def read_and_analyze_data(filename):
    """
    Read location data and analyze by person groups
    """
    # Define location types in order
    location_types = ['Home', 'Work', 'BasicsShop', 'School', 'SocialEvent']
    location_to_index = {location: i for i,
                         location in enumerate(location_types)}

    # maximum number of hours:
    max_hours = 24 * 6  # 7 days of data, 24 hours each

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        print(f"Successfully read {len(lines)} lines from {filename}")

        # Initialize tracking variables
        current_hour = 0
        # {person_id: [Home, Work, BasicsShop, School, SocialEvent]}
        person_totals = {}

        # Process each line
        for line_num, line in enumerate(lines):
            line = line.strip()

            # Skip header line
            if 'Person ID' in line and 'Location Type' in line:
                continue

            # Blank line indicates new hour
            if not line:
                current_hour += 1
                if current_hour >= max_hours:
                    print("Reached maximum hours limit, stopping processing.")
                    break
                continue

            # Parse person ID and location
            if ',' in line:
                try:
                    person_id_str, location_type = line.split(',', 1)
                    person_id = int(person_id_str.strip())
                    location_type = location_type.strip()

                    # Initialize person if not exists
                    if person_id not in person_totals:
                        # [Home, Work, BasicsShop, School, SocialEvent]
                        person_totals[person_id] = [0, 0, 0, 0, 0]

                    # Find location index and increment counter
                    if location_type in location_to_index:
                        location_index = location_to_index[location_type]
                        person_totals[person_id][location_index] += 1
                    else:
                        print(
                            f"Warning: Unknown location type '{location_type}' on line {line_num + 1}")

                except ValueError as e:
                    print(
                        f"Warning: Could not parse line {line_num + 1}: '{line}' - {e}")

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        return None, None, None, None, None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None, None, None, None, None

    if not person_totals:
        print("No valid data found!")
        return None, None, None, None, None

    # Calculate number of days
    total_hours = current_hour
    num_days = total_hours // 24
    if total_hours % 24 > 0:
        print(
            f"Warning: Incomplete day detected with {total_hours % 24} hours.")

    print(f"Processed {total_hours} hours of data across {num_days} days")
    print(f"Number of people: {len(person_totals)}")

    # Identify different groups
    all_people = list(person_totals.keys())
    workers = [pid for pid in all_people if person_totals[pid]
               [1] > 0]  # Index 1 = Work
    # Index 3 = School
    school_attendees = [pid for pid in all_people if person_totals[pid][3] > 0]
    # CCheck if no worker is also a school attendee and vice versa
    no_worker_school_attendee = all(
        pid not in school_attendees for pid in workers)
    no_school_attendee_worker = all(
        pid not in workers for pid in school_attendees)
    if not no_worker_school_attendee:
        print("Warning: Some workers are also school attendees!")
    if not no_school_attendee_worker:
        print("Warning: Some school attendees are also workers!")

    print(
        f"Workers (people with Work hours): {len(workers)} people - IDs: {workers}")
    print(
        f"School attendees (people with School hours): {len(school_attendees)} people - IDs: {school_attendees}")

    return person_totals, location_types, num_days, workers, school_attendees


def calculate_group_averages(person_totals, group_people, num_days):
    """
    Calculate average hours per day for a group of people
    """
    if not group_people:
        return [0, 0, 0, 0, 0]

    # Sum totals across all people in the group
    group_totals = [0, 0, 0, 0, 0]
    for person_id in group_people:
        person_data = person_totals[person_id]
        for i in range(5):
            group_totals[i] += person_data[i]

    # Calculate average per day per person
    num_people = len(group_people)
    averages = [total / (num_days * num_people) for total in group_totals]

    return averages


def create_pie_charts(person_totals, location_types, num_days, workers, school_attendees):
    """
    Create 3 pie charts: All people, Workers, School attendees
    """
    colors = ['#FF9999', '#66B2FF', '#99FF99', '#FFCC99', '#FF99CC']

    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

    # 1. All people average
    all_people = list(person_totals.keys())
    all_averages = calculate_group_averages(
        person_totals, all_people, num_days)

    # Filter out categories with 0 hours for all people
    non_zero_data1 = []
    non_zero_labels1 = []
    non_zero_colors1 = []

    for i, (avg, label) in enumerate(zip(all_averages, location_types)):
        if avg > 0:
            non_zero_data1.append(avg)
            non_zero_labels1.append(f"{label}\n({avg:.1f}h)")
            non_zero_colors1.append(colors[i])

    if non_zero_data1:
        wedges1, texts1, autotexts1 = ax1.pie(non_zero_data1,
                                              labels=non_zero_labels1,
                                              colors=non_zero_colors1,
                                              autopct='%1.1f%%',
                                              startangle=90)

        for autotext in autotexts1:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(10)

    ax1.set_title(f'All People\nAverage per day per person\n({len(all_people)} people)',
                  fontsize=12, fontweight='bold')

    # 2. Workers only
    if workers:
        worker_averages = calculate_group_averages(
            person_totals, workers, num_days)

        non_zero_data2 = []
        non_zero_labels2 = []
        non_zero_colors2 = []

        for i, (avg, label) in enumerate(zip(worker_averages, location_types)):
            if avg > 0:
                non_zero_data2.append(avg)
                non_zero_labels2.append(f"{label}\n({avg:.1f}h)")
                non_zero_colors2.append(colors[i])

        if non_zero_data2:
            wedges2, texts2, autotexts2 = ax2.pie(non_zero_data2,
                                                  labels=non_zero_labels2,
                                                  colors=non_zero_colors2,
                                                  autopct='%1.1f%%',
                                                  startangle=90)

            for autotext in autotexts2:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
                autotext.set_fontsize(10)
    else:
        ax2.text(0.5, 0.5, 'No Workers', ha='center', va='center',
                 transform=ax2.transAxes, fontsize=14)
        ax2.axis('off')

    ax2.set_title(f'Workers Only\nAverage per day per person\n({len(workers)} people)',
                  fontsize=12, fontweight='bold')

    # 3. School attendees only
    if school_attendees:
        school_averages = calculate_group_averages(
            person_totals, school_attendees, num_days)

        non_zero_data3 = []
        non_zero_labels3 = []
        non_zero_colors3 = []

        for i, (avg, label) in enumerate(zip(school_averages, location_types)):
            if avg > 0:
                non_zero_data3.append(avg)
                non_zero_labels3.append(f"{label}\n({avg:.1f}h)")
                non_zero_colors3.append(colors[i])

        if non_zero_data3:
            wedges3, texts3, autotexts3 = ax3.pie(non_zero_data3,
                                                  labels=non_zero_labels3,
                                                  colors=non_zero_colors3,
                                                  autopct='%1.1f%%',
                                                  startangle=90)

            for autotext in autotexts3:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
                autotext.set_fontsize(10)
    else:
        ax3.text(0.5, 0.5, 'No School Attendees', ha='center', va='center',
                 transform=ax3.transAxes, fontsize=14)
        ax3.axis('off')

    ax3.set_title(f'School Attendees Only\nAverage per day per person\n({len(school_attendees)} people)',
                  fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.suptitle('Average Daily Hours by Group',
                 fontsize=16, fontweight='bold', y=1.05)
    plt.show()


def print_detailed_summary(person_totals, location_types, num_days, workers, school_attendees):
    """
    Print detailed summary statistics
    """
    print("\n" + "="*70)
    print("DETAILED SUMMARY")
    print("="*70)

    all_people = list(person_totals.keys())

    # Calculate averages for each group
    all_averages = calculate_group_averages(
        person_totals, all_people, num_days)
    worker_averages = calculate_group_averages(
        person_totals, workers, num_days) if workers else [0]*5
    school_averages = calculate_group_averages(
        person_totals, school_attendees, num_days) if school_attendees else [0]*5

    print(f"Analysis period: {num_days} days")
    print(f"Total people: {len(all_people)}")
    print(f"Workers: {len(workers)} people {workers}")
    print(
        f"School attendees: {len(school_attendees)} people {school_attendees}")

    print(f"\nAverage hours per day per person:")
    print("Group\t\t\t" +
          "\t".join(f"{loc[:8]}" for loc in location_types) + "\tTotal")
    print("-" * 70)

    # All people
    all_total = sum(all_averages)
    print(f"All People\t\t" +
          "\t".join(f"{avg:.1f}" for avg in all_averages) + f"\t{all_total:.1f}")

    # Workers
    if workers:
        worker_total = sum(worker_averages)
        print(f"Workers\t\t\t" +
              "\t".join(f"{avg:.1f}" for avg in worker_averages) + f"\t{worker_total:.1f}")

    # School attendees
    if school_attendees:
        school_total = sum(school_averages)
        print(f"School Attendees\t" +
              "\t".join(f"{avg:.1f}" for avg in school_averages) + f"\t{school_total:.1f}")


def main():
    """
    Main function to create group-based pie charts
    """
    # Use the uploaded file directly
    filename = "/Users/saschakorf/Nosynch/Arbeit/memilio/cpp/examples/panvXabmSim/results/results_2025-07-23203930/location_type_and_id.txt"

    print(f"Reading location data from {filename}...")
    print("Creating pie charts for different groups...")

    # Read and analyze data
    person_totals, location_types, num_days, workers, school_attendees = read_and_analyze_data(
        filename)

    if person_totals is not None:
        # Create pie charts
        create_pie_charts(person_totals, location_types,
                          num_days, workers, school_attendees)

        # Print detailed summary
        print_detailed_summary(person_totals, location_types,
                               num_days, workers, school_attendees)
    else:
        print("Failed to process data. Please check your file format.")


if __name__ == "__main__":
    main()
