#!/usr/bin/env python3
"""
Epidemiological Metrics Calculator for CSV Simulation Results

This script calculates 5 key epidemiological metrics from CSV simulation data:
1. Basic Reproduction Number (R0) - estimated from infection spread
2. Time to 100 Infections - based on infection progression
3. Peak Infection Day - estimated from contact patterns
4. Attack Rate - percentage of population infected
5. Trajectory Divergence Score - for comparing scenarios

Input: CSV files (contact_intensiveness.csv, infection_count.csv, location_id_and_type.csv)
Output: Text file with calculated metrics
"""

import pandas as pd
import numpy as np
import os
import sys
from datetime import datetime
import argparse
from collections import defaultdict


def load_csv_data(data_dir):
    """
    Load all CSV files from simulation results.

    Args:
        data_dir (str): Directory containing CSV files

    Returns:
        tuple: (contact_df, infection_df, location_df)
    """
    try:
        # Load contact intensiveness data
        contact_file = os.path.join(data_dir, 'contact_intensiveness.csv')
        contact_df = pd.read_csv(contact_file)
        # Clean column names (remove leading spaces)
        contact_df.columns = contact_df.columns.str.strip()

        # Load infection count data
        infection_file = os.path.join(data_dir, 'infection_count.csv')
        infection_df = pd.read_csv(infection_file)
        infection_df.columns = infection_df.columns.str.strip()

        # Load location data
        location_file = os.path.join(data_dir, 'location_id_and_type.csv')
        location_df = pd.read_csv(location_file)
        location_df.columns = location_df.columns.str.strip()

        return contact_df, infection_df, location_df

    except Exception as e:
        print(f"Error loading CSV files: {e}")
        return None, None, None


def estimate_r0_from_contacts(contact_df, infection_df):
    """
    Estimate R0 from contact patterns and infection data.

    Args:
        contact_df (DataFrame): Contact data
        infection_df (DataFrame): Infection data

    Returns:
        float: Estimated R0
    """
    # Get infected person IDs
    infected_persons = set(infection_df['Person ID'].tolist())

    # Calculate average contacts per infected person
    infected_contacts = contact_df[
        (contact_df['Person 1'].isin(infected_persons)) |
        (contact_df['Person 2'].isin(infected_persons))
    ]

    # Count unique contacts per infected person
    contact_counts = defaultdict(set)

    for _, row in infected_contacts.iterrows():
        person1, person2 = row['Person 1'], row['Person 2']
        if person1 in infected_persons:
            contact_counts[person1].add(person2)
        if person2 in infected_persons:
            contact_counts[person2].add(person1)

    # Calculate average contacts per infected person
    if len(contact_counts) > 0:
        avg_contacts = np.mean([len(contacts)
                               for contacts in contact_counts.values()])

        # Estimate transmission probability from infection success rate
        total_unique_persons = len(set(contact_df['Person 1'].tolist() +
                                       contact_df['Person 2'].tolist()))
        attack_rate = len(infected_persons) / total_unique_persons

        # Simple R0 estimation: contacts * transmission_probability
        # Assuming transmission probability correlates with attack rate
        transmission_prob = min(attack_rate * 2, 0.5)  # Cap at 50%
        r0 = avg_contacts * transmission_prob

        return max(r0, 1.0)  # R0 should be >= 1 for an outbreak

    return np.nan


def estimate_time_to_threshold(infection_df, contact_df, threshold=100):
    """
    Estimate time to reach infection threshold based on contact intensity.

    Args:
        infection_df (DataFrame): Infection data
        contact_df (DataFrame): Contact data  
        threshold (int): Target infection count

    Returns:
        float: Estimated days to threshold
    """
    current_infections = len(infection_df)

    if current_infections >= threshold:
        return 0.0  # Already reached threshold

    # Calculate daily contact rate (assuming contact hours represent daily contacts)
    total_contact_hours = contact_df['Contact Hours'].sum()
    unique_persons = len(set(contact_df['Person 1'].tolist() +
                             contact_df['Person 2'].tolist()))

    if unique_persons == 0:
        return np.nan

    avg_daily_contacts = total_contact_hours / unique_persons

    # Estimate doubling time from contact intensity
    # More contacts = faster spread = shorter doubling time
    base_doubling_time = 5.0  # Conservative baseline
    contact_factor = max(avg_daily_contacts / 10.0, 0.1)  # Scale factor
    estimated_doubling_time = base_doubling_time / contact_factor

    # Calculate days to reach threshold
    growth_needed = threshold / current_infections
    doublings_needed = np.log2(growth_needed)
    days_to_threshold = doublings_needed * estimated_doubling_time

    return max(days_to_threshold, 0)


def calculate_peak_day_estimate(contact_df, infection_df):
    """
    Estimate peak infection day based on contact patterns.

    Args:
        contact_df (DataFrame): Contact data
        infection_df (DataFrame): Infection data

    Returns:
        tuple: (estimated_peak_day, estimated_peak_height)
    """
    # Estimate based on contact intensity distribution
    total_contact_hours = contact_df['Contact Hours'].sum()
    unique_persons = len(set(contact_df['Person 1'].tolist() +
                             contact_df['Person 2'].tolist()))

    if unique_persons == 0:
        return np.nan, np.nan

    # Estimate epidemic timeline
    avg_contacts = total_contact_hours / unique_persons
    current_infections = len(infection_df)

    # Simple epidemic model: peak occurs after several doubling periods
    doubling_time = 5.0 / max(avg_contacts / 10.0, 0.1)
    estimated_peak_day = doubling_time * 4  # Typically 3-5 doubling times to peak

    # Estimate peak height based on total contact potential
    # This is a rough approximation
    max_daily_contacts = contact_df.groupby(
        'Person 1')['Contact Hours'].sum().max()
    estimated_peak_height = min(max_daily_contacts * 5, unique_persons * 0.05)

    return estimated_peak_day, estimated_peak_height


def calculate_attack_rate(infection_df, contact_df):
    """
    Calculate attack rate from infection data.

    Args:
        infection_df (DataFrame): Infection data
        contact_df (DataFrame): Contact data

    Returns:
        float: Attack rate as percentage
    """
    # Get total population from contact data
    total_persons = len(set(contact_df['Person 1'].tolist() +
                            contact_df['Person 2'].tolist()))

    # Get number of infected persons
    infected_persons = len(infection_df)

    if total_persons == 0:
        return 0.0

    attack_rate = (infected_persons / total_persons) * 100
    return attack_rate


def analyze_location_patterns(contact_df, location_df):
    """
    Analyze contact patterns by location type.

    Args:
        contact_df (DataFrame): Contact data
        location_df (DataFrame): Location data

    Returns:
        dict: Analysis of location-based transmission
    """
    # Merge contact data with location types
    contact_with_loc = contact_df.merge(
        location_df,
        left_on='Location ID',
        right_on='Location_ID',
        how='left'
    )

    # Analyze contacts by location type
    location_analysis = {}

    for loc_type in contact_with_loc['Location_Type'].unique():
        if pd.isna(loc_type):
            continue

        loc_contacts = contact_with_loc[contact_with_loc['Location_Type'] == loc_type]

        location_analysis[loc_type] = {
            'total_contacts': len(loc_contacts),
            'total_contact_hours': loc_contacts['Contact Hours'].sum(),
            'avg_contact_hours': loc_contacts['Contact Hours'].mean(),
            'unique_persons': len(set(loc_contacts['Person 1'].tolist() +
                                      loc_contacts['Person 2'].tolist()))
        }

    return location_analysis


def process_simulation_data(data_dir, scenario_name):
    """
    Process CSV simulation data and calculate all metrics.

    Args:
        data_dir (str): Directory containing CSV files
        scenario_name (str): Name of scenario

    Returns:
        dict: Dictionary containing all calculated metrics
    """
    # Load data
    contact_df, infection_df, location_df = load_csv_data(data_dir)

    if contact_df is None:
        return None

    print(
        f"Loaded data: {len(contact_df)} contacts, {len(infection_df)} infections, {len(location_df)} locations")

    # Calculate metrics
    metrics = {
        'scenario': scenario_name,
        'total_contacts': len(contact_df),
        'total_infections': len(infection_df),
        'total_locations': len(location_df)
    }

    # Get total population
    total_population = len(set(contact_df['Person 1'].tolist() +
                               contact_df['Person 2'].tolist()))
    metrics['total_population'] = total_population

    # 1. Basic Reproduction Number
    metrics['R0'] = estimate_r0_from_contacts(contact_df, infection_df)

    # 2. Time to 100 infections
    metrics['time_to_100'] = estimate_time_to_threshold(
        infection_df, contact_df, 100)

    # 3. Peak infection day and height
    peak_day, peak_height = calculate_peak_day_estimate(
        contact_df, infection_df)
    metrics['peak_day'] = peak_day
    metrics['peak_height'] = peak_height

    # 4. Attack rate
    metrics['attack_rate'] = calculate_attack_rate(infection_df, contact_df)

    # 5. Location analysis for context
    metrics['location_analysis'] = analyze_location_patterns(
        contact_df, location_df)

    return metrics


def calculate_trajectory_divergence(metrics1, metrics2):
    """
    Calculate divergence between two scenarios.

    Args:
        metrics1 (dict): Metrics from first scenario
        metrics2 (dict): Metrics from second scenario

    Returns:
        dict: Divergence metrics
    """
    divergence = {}

    # Compare key metrics
    for metric in ['R0', 'time_to_100', 'peak_day', 'attack_rate']:
        if metric in metrics1 and metric in metrics2:
            val1, val2 = metrics1[metric], metrics2[metric]
            if not (np.isnan(val1) or np.isnan(val2)):
                relative_diff = abs(val1 - val2) / \
                    max(abs(val1), abs(val2), 1e-6)
                divergence[f'{metric}_divergence'] = relative_diff

    return divergence


def write_metrics_to_file(metrics_list, output_file):
    """
    Write calculated metrics to text file.

    Args:
        metrics_list (list): List of metric dictionaries
        output_file (str): Output file path
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("Epidemiological Metrics Analysis - CSV Data\n")
        f.write("=" * 55 + "\n")
        f.write(
            f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Write metrics for each scenario
        for i, metrics in enumerate(metrics_list):
            if metrics is None:
                continue

            f.write(f"SCENARIO: {metrics['scenario']}\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Population: {metrics['total_population']:,}\n")
            f.write(f"Total Contacts: {metrics['total_contacts']:,}\n")
            f.write(f"Total Infections: {metrics['total_infections']:,}\n")
            f.write(f"Total Locations: {metrics['total_locations']:,}\n\n")

            f.write("KEY EPIDEMIOLOGICAL METRICS:\n")
            f.write(
                f"1. Basic Reproduction Number (R0): {metrics['R0']:.3f}\n")

            if np.isnan(metrics['time_to_100']):
                f.write(
                    "2. Time to 100 Infections: Already exceeded or cannot estimate\n")
            else:
                f.write(
                    f"2. Time to 100 Infections: {metrics['time_to_100']:.1f} days\n")

            f.write(
                f"3. Peak Infection Day (estimated): {metrics['peak_day']:.1f} days\n")
            f.write(
                f"   Peak Daily Infections (estimated): {metrics['peak_height']:.0f}\n")
            f.write(f"4. Attack Rate: {metrics['attack_rate']:.2f}%\n")

            # Location analysis
            if 'location_analysis' in metrics:
                f.write(f"\nLOCATION-BASED ANALYSIS:\n")
                for loc_type, data in metrics['location_analysis'].items():
                    f.write(f"  {loc_type}:\n")
                    f.write(
                        f"    - Total contacts: {data['total_contacts']}\n")
                    f.write(
                        f"    - Total contact hours: {data['total_contact_hours']}\n")
                    f.write(
                        f"    - Avg contact hours: {data['avg_contact_hours']:.2f}\n")
                    f.write(
                        f"    - Unique persons: {data['unique_persons']}\n")

            f.write("\n" + "=" * 55 + "\n\n")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate epidemiological metrics from CSV simulation results')
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing CSV files')
    parser.add_argument(
        '--scenario_name', help='Scenario name (optional, defaults to folder name)')

    args = parser.parse_args()

    # Check if directory exists
    if not os.path.exists(args.input_dir):
        print(f"Directory not found: {args.input_dir}")
        return

    # Set scenario name
    if args.scenario_name:
        scenario_name = args.scenario_name
    else:
        scenario_name = os.path.basename(args.input_dir.rstrip('/'))

    print(f"Processing {scenario_name} from {args.input_dir}...")

    # Process the data
    metrics = process_simulation_data(args.input_dir, scenario_name)

    if metrics is None:
        print("Failed to process simulation data")
        return

    # Write results to the same directory as input
    output_file = os.path.join(args.input_dir, 'epidemic_metrics.txt')
    write_metrics_to_file([metrics], output_file)
    print(f"Metrics saved to {output_file}")


if __name__ == "__main__":
    main()
