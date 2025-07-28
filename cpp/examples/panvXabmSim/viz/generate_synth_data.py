import pandas as pd
import numpy as np
import random
from collections import defaultdict
import itertools


class ContactNetworkDataGenerator:
    def __init__(self, num_persons=100, num_runs=25, seed=42):
        """
        Initialize the data generator

        Parameters:
        num_persons: Total number of people in the simulation
        num_runs: Number of simulation runs
        seed: Random seed for reproducibility
        """
        self.num_persons = num_persons
        self.num_runs = num_runs
        self.seed = seed

        # Set random seeds
        np.random.seed(seed)
        random.seed(seed)

        # Define reasonable infrastructure numbers based on population
        # Average 4 people per household
        self.num_households = max(10, num_persons // 4)
        # Average 20 people per workplace
        self.num_workplaces = max(5, num_persons // 20)
        # Average 50 people per school
        self.num_schools = max(2, num_persons // 50)
        # Various social event locations
        self.num_social_venues = max(3, num_persons // 15)
        self.num_shopping_centers = max(
            2, num_persons // 25)  # Shopping centers

        # Age groups and their probabilities
        self.age_groups = {
            'child': 0.20,      # 0-17 (go to school)
            'adult': 0.65,      # 18-64 (go to work)
            'elderly': 0.15     # 65+ (mostly home/shopping)
        }

        # Contact probabilities by location type
        self.contact_probabilities = {
            'Home': 0.95,           # Very high within household
            'Work': 0.3,            # Moderate at work
            'School': 0.4,          # Higher at school
            'Social_Event': 0.6,    # High at social events
            'Shopping': 0.1         # Low at shopping
        }

        # Initialize data structures
        self.persons = []
        self.households = {}
        self.workplaces = {}
        self.schools = {}

    def generate_households(self):
        """Generate household assignments"""
        print("Generating households...")

        household_data = []
        person_id = 1

        for household_id in range(1, self.num_households + 1):
            # Determine household size (weighted towards smaller households)
            household_size_probs = [0.15, 0.25,
                                    0.25, 0.20, 0.10, 0.05]  # 1-6 people
            household_size = np.random.choice(
                range(1, 7), p=household_size_probs)

            household_members = []

            for _ in range(household_size):
                if person_id <= self.num_persons:
                    # Assign age group based on household composition
                    if len(household_members) == 0:  # First person (head of household)
                        age_group = np.random.choice(
                            ['adult', 'elderly'], p=[0.8, 0.2])
                    elif len(household_members) == 1 and household_size > 2:  # Spouse
                        # Similar age
                        age_group = household_members[0]['age_group']
                    else:  # Children or other family members
                        age_group = np.random.choice(
                            ['child', 'adult'], p=[0.7, 0.3])

                    person = {
                        'person_id': person_id,
                        'household_id': f'H{household_id}',
                        'age_group': age_group
                    }

                    household_members.append(person)
                    household_data.append({
                        'person_id': person_id,
                        'household_id': f'H{household_id}'
                    })

                    person_id += 1

            self.households[f'H{household_id}'] = household_members

        return pd.DataFrame(household_data)

    def assign_locations(self):
        """Assign people to workplaces and schools"""
        print("Assigning locations...")

        # Initialize location assignments
        workplace_assignments = defaultdict(list)
        school_assignments = defaultdict(list)

        for household_id, members in self.households.items():
            for person in members:
                person_id = person['person_id']
                age_group = person['age_group']

                # Assign workplace to adults
                if age_group == 'adult':
                    # 90% of adults work
                    if np.random.random() < 0.9:
                        workplace_id = f'W{np.random.randint(1, self.num_workplaces + 1)}'
                        workplace_assignments[workplace_id].append(person_id)
                        person['workplace'] = workplace_id

                # Assign school to children
                elif age_group == 'child':
                    school_id = f'S{np.random.randint(1, self.num_schools + 1)}'
                    school_assignments[school_id].append(person_id)
                    person['school'] = school_id

        self.workplaces = dict(workplace_assignments)
        self.schools = dict(school_assignments)

    def generate_contacts(self):
        """Generate contact data across all location types"""
        print("Generating contact data...")

        contact_data = []
        time_counter = 0.0

        # Generate contacts for multiple days
        for day in range(7):  # One week of contacts
            daily_time_offset = day * 24

            # Morning household contacts (6-9 AM)
            contact_data.extend(self._generate_household_contacts(
                daily_time_offset + 6, daily_time_offset + 9, 'morning'))

            # Work/School contacts (9 AM - 5 PM)
            contact_data.extend(self._generate_workplace_contacts(
                daily_time_offset + 9, daily_time_offset + 17))
            contact_data.extend(self._generate_school_contacts(
                daily_time_offset + 8, daily_time_offset + 15))

            # Shopping contacts (throughout day)
            contact_data.extend(self._generate_shopping_contacts(
                daily_time_offset + 10, daily_time_offset + 20))

            # Evening household contacts (5-10 PM)
            contact_data.extend(self._generate_household_contacts(
                daily_time_offset + 17, daily_time_offset + 22, 'evening'))

            # Social events (weekends and evenings)
            if day >= 5 or np.random.random() < 0.3:  # Weekend or weekday evening event
                contact_data.extend(self._generate_social_contacts(
                    daily_time_offset + 18, daily_time_offset + 23))

        return pd.DataFrame(contact_data)

    def _generate_household_contacts(self, start_time, end_time, period):
        """Generate contacts within households"""
        contacts = []

        for household_id, members in self.households.items():
            if len(members) < 2:
                continue

            # Generate contacts between household members
            for person1, person2 in itertools.combinations(members, 2):
                # Multiple contacts during this period
                # Average 3 contacts per period
                num_contacts = np.random.poisson(3)

                for _ in range(num_contacts):
                    if np.random.random() < self.contact_probabilities['Home']:
                        contact_time = np.random.uniform(start_time, end_time)

                        contacts.append({
                            'person_1': person1['person_id'],
                            'person_2': person2['person_id'],
                            'time': round(contact_time, 2),
                            'location_type': 'Home',
                            'location_id': household_id
                        })

        return contacts

    def _generate_workplace_contacts(self, start_time, end_time):
        """Generate contacts in workplaces"""
        contacts = []

        for workplace_id, workers in self.workplaces.items():
            if len(workers) < 2:
                continue

            # Generate contacts between coworkers
            for person1, person2 in itertools.combinations(workers, 2):
                # Fewer contacts at work than at home
                num_contacts = np.random.poisson(1.5)

                for _ in range(num_contacts):
                    if np.random.random() < self.contact_probabilities['Work']:
                        contact_time = np.random.uniform(start_time, end_time)

                        contacts.append({
                            'person_1': person1,
                            'person_2': person2,
                            'time': round(contact_time, 2),
                            'location_type': 'Work',
                            'location_id': workplace_id
                        })

        return contacts

    def _generate_school_contacts(self, start_time, end_time):
        """Generate contacts in schools"""
        contacts = []

        for school_id, students in self.schools.items():
            if len(students) < 2:
                continue

            # Generate contacts between students
            for person1, person2 in itertools.combinations(students, 2):
                # Higher contact rate at school
                num_contacts = np.random.poisson(2)

                for _ in range(num_contacts):
                    if np.random.random() < self.contact_probabilities['School']:
                        contact_time = np.random.uniform(start_time, end_time)

                        contacts.append({
                            'person_1': person1,
                            'person_2': person2,
                            'time': round(contact_time, 2),
                            'location_type': 'School',
                            'location_id': school_id
                        })

        return contacts

    def _generate_shopping_contacts(self, start_time, end_time):
        """Generate contacts at shopping centers"""
        contacts = []

        # Random subset of people go shopping each day
        all_persons = [p['person_id']
                       for household in self.households.values() for p in household]
        shoppers_per_center = {}

        for shopping_center in range(1, self.num_shopping_centers + 1):
            # 20-30% of people might shop on any given day
            num_shoppers = int(len(all_persons) *
                               np.random.uniform(0.05, 0.15))
            shoppers = np.random.choice(
                all_persons, size=num_shoppers, replace=False)
            shoppers_per_center[f'SH{shopping_center}'] = shoppers

        for location_id, shoppers in shoppers_per_center.items():
            if len(shoppers) < 2:
                continue

            # Generate random encounters between shoppers
            for person1, person2 in itertools.combinations(shoppers, 2):
                if np.random.random() < self.contact_probabilities['Shopping']:
                    contact_time = np.random.uniform(start_time, end_time)

                    contacts.append({
                        'person_1': person1,
                        'person_2': person2,
                        'time': round(contact_time, 2),
                        'location_type': 'Shopping',
                        'location_id': location_id
                    })

        return contacts

    def _generate_social_contacts(self, start_time, end_time):
        """Generate contacts at social events"""
        contacts = []

        # Generate 1-3 social events
        num_events = np.random.randint(1, 4)

        for event_num in range(1, num_events + 1):
            event_id = f'E{event_num}'

            # Random subset of people attend each social event
            all_persons = [p['person_id']
                           for household in self.households.values() for p in household]
            num_attendees = int(len(all_persons) * np.random.uniform(0.1, 0.3))
            attendees = np.random.choice(
                all_persons, size=num_attendees, replace=False)

            # Generate contacts between attendees
            for person1, person2 in itertools.combinations(attendees, 2):
                if np.random.random() < self.contact_probabilities['Social_Event']:
                    contact_time = np.random.uniform(start_time, end_time)

                    contacts.append({
                        'person_1': person1,
                        'person_2': person2,
                        'time': round(contact_time, 2),
                        'location_type': 'Social_Event',
                        'location_id': event_id
                    })

        return contacts

    def generate_infections(self, contact_data):
        """Generate infection data based on contacts and transmission probabilities"""
        print("Generating infection data...")

        infection_data = []

        # Transmission probabilities by location type
        transmission_probs = {
            'Home': 0.15,
            'Work': 0.05,
            'School': 0.08,
            'Social_Event': 0.12,
            'Shopping': 0.02
        }

        for run_id in range(1, self.num_runs + 1):
            # Initialize infection status for this run
            infected_persons = set()
            infection_times = {}
            infection_locations = {}

            # Choose 1-3 initial infected persons (patient zero)
            initial_infected = np.random.choice(
                range(1, self.num_persons + 1),
                size=np.random.randint(1, 4),
                replace=False
            )

            for person_id in initial_infected:
                infected_persons.add(person_id)
                infection_times[person_id] = 0.0
                infection_locations[person_id] = 'External'

            # Sort contacts by time to simulate disease spread
            run_contacts = contact_data.copy()
            run_contacts = run_contacts.sort_values('time')

            # Simulate transmission through contacts
            for _, contact in run_contacts.iterrows():
                person1, person2 = contact['person_1'], contact['person_2']
                location_type = contact['location_type']
                contact_time = contact['time']

                # Check if one person is infected and the other is not
                if person1 in infected_persons and person2 not in infected_persons:
                    if np.random.random() < transmission_probs.get(location_type, 0.05):
                        infected_persons.add(person2)
                        infection_times[person2] = contact_time
                        infection_locations[person2] = location_type

                elif person2 in infected_persons and person1 not in infected_persons:
                    if np.random.random() < transmission_probs.get(location_type, 0.05):
                        infected_persons.add(person1)
                        infection_times[person1] = contact_time
                        infection_locations[person1] = location_type

            # Record infection status for all persons in this run
            for person_id in range(1, self.num_persons + 1):
                if person_id in infected_persons:
                    infection_data.append({
                        'person_id': person_id,
                        'run_id': run_id,
                        'infected': 1,
                        'infection_time': infection_times[person_id],
                        'infection_location': infection_locations[person_id]
                    })
                else:
                    infection_data.append({
                        'person_id': person_id,
                        'run_id': run_id,
                        'infected': 0,
                        'infection_time': '',
                        'infection_location': ''
                    })

        return pd.DataFrame(infection_data)

    def generate_all_data(self, output_prefix='generated'):
        """Generate all data files"""
        print(
            f"Generating data for {self.num_persons} persons and {self.num_runs} runs...")
        print(
            f"Infrastructure: {self.num_households} households, {self.num_workplaces} workplaces, {self.num_schools} schools")

        # Generate household data
        household_data = self.generate_households()

        # Assign locations
        self.assign_locations()

        # Generate contact data
        contact_data = self.generate_contacts()

        # Generate infection data
        infection_data = self.generate_infections(contact_data)

        # Save to CSV files
        household_file = f'{output_prefix}_household_data.csv'
        contact_file = f'{output_prefix}_contact_data.csv'
        infection_file = f'{output_prefix}_infection_data.csv'

        household_data.to_csv(household_file, index=False)
        contact_data.to_csv(contact_file, index=False)
        infection_data.to_csv(infection_file, index=False)

        print(f"\nGenerated files:")
        print(f"- {household_file}: {len(household_data)} household assignments")
        print(f"- {contact_file}: {len(contact_data)} contact records")
        print(f"- {infection_file}: {len(infection_data)} infection records")

        # Print summary statistics
        print(f"\nSummary statistics:")
        print(
            f"- Average household size: {household_data.groupby('household_id').size().mean():.2f}")
        print(f"- Total contacts: {len(contact_data)}")
        print(f"- Contacts by location type:")
        for loc_type, count in contact_data['location_type'].value_counts().items():
            print(f"  - {loc_type}: {count}")

        infection_rates = infection_data.groupby(
            'person_id')['infected'].mean() * 100
        print(f"- Average infection rate: {infection_rates.mean():.1f}%")
        print(
            f"- Infection rate range: {infection_rates.min():.1f}% - {infection_rates.max():.1f}%")

        return household_data, contact_data, infection_data


def main():
    """Example usage of the data generator"""

    # Small dataset for testing
    print("=== Generating Small Dataset ===")
    small_generator = ContactNetworkDataGenerator(
        num_persons=5, num_runs=1, seed=42)
    small_generator.generate_all_data('small')

    print("\n" + "="*50 + "\n")

    # Medium dataset
    print("=== Generating Medium Dataset ===")
    medium_generator = ContactNetworkDataGenerator(
        num_persons=200, num_runs=25, seed=123)
    medium_generator.generate_all_data('medium')

    print("\n" + "="*50 + "\n")

    # Large dataset
    # print("=== Generating Large Dataset ===")
    # large_generator = ContactNetworkDataGenerator(
    #     num_persons=500, num_runs=50, seed=456)
    # large_generator.generate_all_data('large')


if __name__ == "__main__":
    main()
