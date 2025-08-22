#!/usr/bin/env python3
"""
Comprehensive Epidemiological Analysis for Agent-Based Model Simulation Results

This script performs advanced epidemiological analysis on ABM simulation data including:
1. Temporal transmission dynamics analysis
2. Contact network structure analysis  
3. Household and location-based transmission patterns
4. Transmission tree analysis (generation times, serial intervals)
5. Superspreading event detection
6. Secondary attack rates by context
7. Spatial clustering analysis
8. Epidemic curve characteristics
9. Network centrality and clustering metrics
10. Comparative scenario analysis

Input: Multiple CSV files from ABM simulation
Output: Comprehensive epidemiological report
"""

import pandas as pd
import numpy as np
import os
import sys
from datetime import datetime
import argparse
from collections import defaultdict, Counter
import networkx as nx
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
import warnings
warnings.filterwarnings('ignore')


def load_comprehensive_data(data_dir):
    """
    Load all available CSV files from simulation results.

    Args:
        data_dir (str): Directory containing CSV files

    Returns:
        dict: Dictionary containing all loaded DataFrames
    """
    data = {}

    try:
        # Core files
        files_to_load = {
            'contact_intensiveness': 'contact_intensiveness.csv',
            'infection_count': 'infection_count.csv',
            'location_data': 'location_id_and_type.csv',
            'household_data': 'household_id.csv',
            'detailed_infections': 'best_run_detailed_infection.csv',
            'contact_data': 'best_run_contact_data.csv'
        }

        for key, filename in files_to_load.items():
            filepath = os.path.join(data_dir, filename)
            if os.path.exists(filepath):
                data[key] = pd.read_csv(filepath)
                # Clean column names
                data[key].columns = data[key].columns.str.strip()
                print(f"Loaded {key}: {len(data[key])} records")
            else:
                print(f"Warning: {filename} not found")
                data[key] = None

        return data

    except Exception as e:
        print(f"Error loading data: {e}")
        return None


def build_transmission_network(detailed_infections_df, contact_data_df, time_window=24):
    """
    Build transmission network from detailed infection and contact data.

    Args:
        detailed_infections_df: DataFrame with infection events
        contact_data_df: DataFrame with contact data
        time_window: Time window to look for contacts before infection

    Returns:
        tuple: (transmission_graph, infection_timeline, transmission_pairs)
    """
    if detailed_infections_df is None or contact_data_df is None:
        return None, None, None

    # Sort infections by time
    infections = detailed_infections_df.sort_values('Timestep').copy()

    # Build infection timeline
    infection_timeline = {}
    for _, row in infections.iterrows():
        person_id = row['Person_ID']
        timestep = row['Timestep']
        location_id = row['Location_ID']
        location_type = row['Location_Type']

        infection_timeline[person_id] = {
            'time': timestep,
            'location_id': location_id,
            'location_type': location_type
        }

    # Build transmission network
    transmission_graph = nx.DiGraph()
    transmission_pairs = []

    # Patient zero (earliest infection)
    patient_zero = infections.iloc[0]['Person_ID']
    transmission_graph.add_node(
        patient_zero, generation=0, time=infections.iloc[0]['Timestep'])

    for _, infected_row in infections.iloc[1:].iterrows():
        infected_person = infected_row['Person_ID']
        infection_time = infected_row['Timestep']
        infection_location = infected_row['Location_ID']

        # Find potential infectors (people at same location within time window)
        potential_infectors = []

        # Look for contacts at the infection location within time window
        relevant_contacts = contact_data_df[
            (contact_data_df['Location_ID'] == infection_location) &
            (contact_data_df['Timestep'] >= infection_time - time_window) &
            (contact_data_df['Timestep'] <= infection_time) &
            (contact_data_df['Person_ID'] != infected_person)
        ]

        # Filter to only previously infected people
        for _, contact_row in relevant_contacts.iterrows():
            contact_person = contact_row['Person_ID']
            if contact_person in infection_timeline:
                contact_infection_time = infection_timeline[contact_person]['time']
                if contact_infection_time < infection_time:
                    time_diff = infection_time - contact_infection_time
                    potential_infectors.append((contact_person, time_diff))

        # Choose most likely infector (most recent infection)
        if potential_infectors:
            infector = min(potential_infectors, key=lambda x: x[1])[0]
            transmission_graph.add_edge(infector, infected_person)
            transmission_pairs.append(
                (infector, infected_person, infection_time, infection_location))

            # Set generation
            if infector in transmission_graph.nodes:
                infector_gen = transmission_graph.nodes[infector].get(
                    'generation', 0)
                transmission_graph.add_node(
                    infected_person, generation=infector_gen + 1, time=infection_time)
        else:
            # No clear infector found - might be external
            transmission_graph.add_node(
                infected_person, generation=1, time=infection_time)

    return transmission_graph, infection_timeline, transmission_pairs


def calculate_reproduction_numbers(transmission_graph):
    """
    Calculate various reproduction number metrics.

    Args:
        transmission_graph: NetworkX directed graph of transmissions

    Returns:
        dict: Dictionary with reproduction number metrics
    """
    if transmission_graph is None:
        return {}

    # Calculate individual reproduction numbers
    reproduction_numbers = []
    for node in transmission_graph.nodes():
        offspring_count = len(list(transmission_graph.successors(node)))
        reproduction_numbers.append(offspring_count)

    reproduction_numbers = np.array(reproduction_numbers)

    metrics = {
        'R_mean': np.mean(reproduction_numbers),
        'R_median': np.median(reproduction_numbers),
        'R_std': np.std(reproduction_numbers),
        'R_max': np.max(reproduction_numbers),
        'R_min': np.min(reproduction_numbers),
        'R0_estimate': np.mean([reproduction_numbers[i] for i, n in enumerate(transmission_graph.nodes())
                               if transmission_graph.nodes[n].get('generation', 0) == 0]) if len(transmission_graph.nodes()) > 0 else 0,
        'superspreader_threshold': np.percentile(reproduction_numbers, 90) if len(reproduction_numbers) > 0 else 0,
        'superspreader_count': np.sum(reproduction_numbers >= np.percentile(reproduction_numbers, 90)) if len(reproduction_numbers) > 0 else 0
    }

    return metrics


def analyze_generation_times(transmission_graph, transmission_pairs):
    """
    Analyze generation times and serial intervals.

    Args:
        transmission_graph: NetworkX transmission graph
        transmission_pairs: List of transmission pairs with timing

    Returns:
        dict: Generation time analysis
    """
    if not transmission_pairs:
        return {}

    generation_times = []
    serial_intervals = []

    for infector, infected, infection_time, location in transmission_pairs:
        if infector in transmission_graph.nodes and infected in transmission_graph.nodes:
            infector_time = transmission_graph.nodes[infector].get('time', 0)
            serial_interval = infection_time - infector_time
            serial_intervals.append(serial_interval)

            # Generation time approximation (serial interval)
            generation_times.append(serial_interval)

    if len(generation_times) == 0:
        return {}

    return {
        'generation_time_mean': np.mean(generation_times),
        'generation_time_median': np.median(generation_times),
        'generation_time_std': np.std(generation_times),
        'serial_interval_mean': np.mean(serial_intervals),
        'serial_interval_median': np.median(serial_intervals),
        'serial_interval_std': np.std(serial_intervals),
        'min_serial_interval': np.min(serial_intervals),
        'max_serial_interval': np.max(serial_intervals)
    }


def analyze_household_transmission(detailed_infections_df, household_data_df):
    """
    Analyze household transmission patterns.

    Args:
        detailed_infections_df: Infection events
        household_data_df: Household assignments

    Returns:
        dict: Household transmission analysis
    """
    if detailed_infections_df is None or household_data_df is None:
        return {}

    # Merge infection data with household data
    infection_households = detailed_infections_df.merge(
        household_data_df, on='Person_ID', how='left'
    )

    # Analyze household clusters
    household_infections = infection_households.groupby('Household_ID').agg({
        'Person_ID': 'count',
        'Timestep': ['min', 'max', 'mean']
    }).round(2)

    household_infections.columns = [
        'infected_count', 'first_infection', 'last_infection', 'mean_infection_time']

    # Calculate household secondary attack rates
    total_households = len(household_data_df['Household_ID'].unique())
    households_with_infections = len(
        household_infections[household_infections['infected_count'] > 0])
    households_with_multiple = len(
        household_infections[household_infections['infected_count'] > 1])

    # Household size distribution
    household_sizes = household_data_df.groupby('Household_ID').size()

    metrics = {
        'total_households': total_households,
        'households_with_infections': households_with_infections,
        'households_with_multiple_infections': households_with_multiple,
        'household_attack_rate': households_with_infections / total_households * 100,
        'household_secondary_transmission_rate': households_with_multiple / households_with_infections * 100 if households_with_infections > 0 else 0,
        'mean_household_size': household_sizes.mean(),
        'median_household_size': household_sizes.median(),
        'max_household_cluster_size': household_infections['infected_count'].max() if len(household_infections) > 0 else 0,
        'mean_household_cluster_size': household_infections['infected_count'].mean() if len(household_infections) > 0 else 0
    }

    return metrics


def analyze_location_transmission_risk(detailed_infections_df, contact_data_df, location_data_df):
    """
    Analyze transmission risk by location type.

    Args:
        detailed_infections_df: Infection events
        contact_data_df: Contact data
        location_data_df: Location type data

    Returns:
        dict: Location-based transmission risk analysis
    """
    if detailed_infections_df is None:
        return {}

    # Check if Location_Type is already in detailed_infections_df
    if 'Location_Type' in detailed_infections_df.columns:
        infections_with_location = detailed_infections_df.copy()
    else:
        # Merge with location types if needed
        if location_data_df is None:
            return {}
        infections_with_location = detailed_infections_df.merge(
            location_data_df,
            left_on='Location_ID',
            right_on='Location_ID',
            how='left'
        )

    # Count infections by location type
    location_infections = infections_with_location.groupby('Location_Type').agg({
        'Person_ID': 'count',
        'Timestep': ['min', 'max', 'mean']
    })

    # If we have contact data, calculate exposure rates
    location_metrics = {}

    for location_type in infections_with_location['Location_Type'].unique():
        if pd.isna(location_type):
            continue

        type_infections = len(
            infections_with_location[infections_with_location['Location_Type'] == location_type])

        # Calculate person-hours of exposure if contact data available
        if contact_data_df is not None and location_data_df is not None:
            type_locations = location_data_df[location_data_df['Location_Type']
                                              == location_type]['Location_ID'].tolist()
            type_contacts = contact_data_df[contact_data_df['Location_ID'].isin(
                type_locations)]
            total_person_visits = len(type_contacts)
            unique_persons_exposed = len(
                type_contacts['Person_ID'].unique()) if len(type_contacts) > 0 else 0

            location_metrics[location_type] = {
                'total_infections': type_infections,
                'total_person_visits': total_person_visits,
                'unique_persons_exposed': unique_persons_exposed,
                'transmission_rate_per_visit': type_infections / total_person_visits * 100 if total_person_visits > 0 else 0,
                'attack_rate_per_person': type_infections / unique_persons_exposed * 100 if unique_persons_exposed > 0 else 0
            }
        else:
            location_metrics[location_type] = {
                'total_infections': type_infections
            }

    return location_metrics


def analyze_contact_network_structure(contact_data_df, detailed_infections_df):
    """
    Analyze the structure of the contact network.

    Args:
        contact_data_df: Contact data
        detailed_infections_df: Infection data

    Returns:
        dict: Network structure metrics
    """
    if contact_data_df is None:
        return {}

    # Build contact network (undirected)
    contact_network = nx.Graph()

    # Add edges for contacts
    contact_pairs = contact_data_df.groupby(
        ['Person_ID', 'Location_ID']).size().reset_index(name='contact_count')

    # Create person-person contacts at same location-time
    for location_id in contact_data_df['Location_ID'].unique():
        for timestep in contact_data_df['Timestep'].unique():
            people_at_location = contact_data_df[
                (contact_data_df['Location_ID'] == location_id) &
                (contact_data_df['Timestep'] == timestep)
            ]['Person_ID'].tolist()

            # Add edges between all pairs at same location-time
            for i, person1 in enumerate(people_at_location):
                for person2 in people_at_location[i+1:]:
                    if contact_network.has_edge(person1, person2):
                        contact_network[person1][person2]['weight'] += 1
                    else:
                        contact_network.add_edge(person1, person2, weight=1)

    if len(contact_network.nodes()) == 0:
        return {}

    # Calculate network metrics
    metrics = {
        'total_nodes': len(contact_network.nodes()),
        'total_edges': len(contact_network.edges()),
        'density': nx.density(contact_network),
        'average_degree': np.mean([d for n, d in contact_network.degree()]),
        'max_degree': max([d for n, d in contact_network.degree()]) if len(contact_network.nodes()) > 0 else 0
    }

    # Connected components
    if len(contact_network.nodes()) > 0:
        components = list(nx.connected_components(contact_network))
        metrics['num_connected_components'] = len(components)
        metrics['largest_component_size'] = len(
            max(components, key=len)) if components else 0
        metrics['largest_component_fraction'] = metrics['largest_component_size'] / \
            len(contact_network.nodes())

    # Clustering coefficient
    if len(contact_network.nodes()) > 2:
        metrics['average_clustering'] = nx.average_clustering(contact_network)

    # Centrality measures for infected individuals
    if detailed_infections_df is not None:
        infected_persons = set(detailed_infections_df['Person_ID'].tolist())
        infected_in_network = [
            p for p in infected_persons if p in contact_network.nodes()]

        if infected_in_network:
            degree_centrality = nx.degree_centrality(contact_network)
            infected_centralities = [degree_centrality[p]
                                     for p in infected_in_network]

            metrics['infected_mean_centrality'] = np.mean(
                infected_centralities)
            metrics['infected_max_centrality'] = np.max(infected_centralities)

    return metrics


def analyze_epidemic_curve(detailed_infections_df):
    """
    Analyze the epidemic curve characteristics.

    Args:
        detailed_infections_df: Infection events with timestamps

    Returns:
        dict: Epidemic curve metrics
    """
    if detailed_infections_df is None:
        return {}

    # Daily infection counts
    daily_infections = detailed_infections_df.groupby('Timestep').size()

    if len(daily_infections) == 0:
        return {}

    # Find peak
    peak_day = daily_infections.idxmax()
    peak_infections = daily_infections.max()

    # Growth phase analysis
    cumulative_infections = daily_infections.cumsum()

    # Exponential growth phase (first half to peak)
    growth_phase = daily_infections[daily_infections.index <= peak_day]

    metrics = {
        'total_days': len(daily_infections),
        'peak_day': peak_day,
        'peak_infections': peak_infections,
        'days_to_peak': len(growth_phase),
        'mean_daily_infections': daily_infections.mean(),
        'std_daily_infections': daily_infections.std(),
        'total_infections': len(detailed_infections_df),
        'attack_rate_time_series': list(cumulative_infections.values)
    }

    # Doubling time calculation (exponential phase)
    if len(growth_phase) > 3:
        try:
            # Fit exponential to growth phase
            growth_values = growth_phase.values
            growth_times = np.arange(len(growth_values))

            # Log-linear regression
            log_values = np.log(np.maximum(growth_values, 1))
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                growth_times, log_values)

            doubling_time = np.log(2) / slope if slope > 0 else np.inf
            metrics['doubling_time'] = doubling_time
            metrics['exponential_fit_r_squared'] = r_value**2

        except:
            metrics['doubling_time'] = np.nan
            metrics['exponential_fit_r_squared'] = np.nan

    return metrics


def calculate_superspreading_metrics(transmission_graph, percentile_threshold=80):
    """
    Identify and analyze superspreading events.

    Args:
        transmission_graph: NetworkX transmission graph
        percentile_threshold: Percentile threshold for superspreading

    Returns:
        dict: Superspreading analysis
    """
    if transmission_graph is None:
        return {}

    # Get reproduction numbers
    reproduction_numbers = []
    superspreaders = []

    for node in transmission_graph.nodes():
        offspring_count = len(list(transmission_graph.successors(node)))
        reproduction_numbers.append(offspring_count)
        if offspring_count >= np.percentile(reproduction_numbers, percentile_threshold):
            superspreaders.append(node)

    if len(reproduction_numbers) == 0:
        return {}

    threshold = np.percentile(reproduction_numbers, percentile_threshold)
    superspreader_events = np.sum(np.array(reproduction_numbers) >= threshold)

    metrics = {
        'superspreader_threshold': threshold,
        'num_superspreaders': superspreader_events,
        'superspreader_fraction': superspreader_events / len(reproduction_numbers) * 100,
        'infections_from_superspreaders': np.sum([r for r in reproduction_numbers if r >= threshold]),
        'superspreader_contribution': np.sum([r for r in reproduction_numbers if r >= threshold]) / np.sum(reproduction_numbers) * 100 if np.sum(reproduction_numbers) > 0 else 0
    }

    return metrics


def comprehensive_epidemiological_analysis(data_dir, scenario_name):
    """
    Perform comprehensive epidemiological analysis on simulation data.

    Args:
        data_dir (str): Directory containing simulation CSV files
        scenario_name (str): Name of the scenario

    Returns:
        dict: Comprehensive epidemiological metrics
    """
    print(f"Loading data from {data_dir}...")
    data = load_comprehensive_data(data_dir)

    if data is None:
        return None

    # Initialize results
    results = {
        'scenario': scenario_name,
        'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    }

    # Basic data summary
    results['data_summary'] = {}
    for key, df in data.items():
        if df is not None:
            results['data_summary'][key] = {
                'records': len(df),
                'columns': list(df.columns)
            }

    # 1. Build transmission network and analyze transmission dynamics
    print("Building transmission network...")
    transmission_graph, infection_timeline, transmission_pairs = build_transmission_network(
        data['detailed_infections'], data['contact_data']
    )

    # 2. Calculate reproduction numbers
    print("Calculating reproduction numbers...")
    results['reproduction_analysis'] = calculate_reproduction_numbers(
        transmission_graph)

    # 3. Analyze generation times and serial intervals
    print("Analyzing generation times...")
    results['generation_analysis'] = analyze_generation_times(
        transmission_graph, transmission_pairs)

    # 4. Household transmission analysis
    print("Analyzing household transmission...")
    results['household_analysis'] = analyze_household_transmission(
        data['detailed_infections'], data['household_data']
    )

    # 5. Location-based transmission risk
    print("Analyzing location transmission risk...")
    results['location_risk_analysis'] = analyze_location_transmission_risk(
        data['detailed_infections'], data['contact_data'], data['location_data']
    )

    # 6. Contact network structure analysis
    print("Analyzing contact network structure...")
    results['network_analysis'] = analyze_contact_network_structure(
        data['contact_data'], data['detailed_infections']
    )

    # 7. Epidemic curve analysis
    print("Analyzing epidemic curve...")
    results['epidemic_curve_analysis'] = analyze_epidemic_curve(
        data['detailed_infections'])

    # 8. Superspreading analysis
    print("Analyzing superspreading events...")
    results['superspreading_analysis'] = calculate_superspreading_metrics(
        transmission_graph)

    # 9. Traditional epidemiological metrics (legacy support)
    print("Calculating traditional metrics...")
    if data['contact_intensiveness'] is not None and data['infection_count'] is not None:
        total_population = len(set(
            data['contact_intensiveness']['Person 1'].tolist() +
            data['contact_intensiveness']['Person 2'].tolist()
        )) if data['contact_intensiveness'] is not None else 0

        results['traditional_metrics'] = {
            'total_population': total_population,
            'total_infections': len(data['infection_count']) if data['infection_count'] is not None else 0,
            'attack_rate': len(data['infection_count']) / total_population * 100 if total_population > 0 else 0,
            'total_contacts': len(data['contact_intensiveness']) if data['contact_intensiveness'] is not None else 0
        }

    # 10. Transmission chain analysis
    if transmission_graph is not None:
        print("Analyzing transmission chains...")
        results['transmission_chain_analysis'] = {
            'total_transmission_events': len(transmission_pairs),
            'max_generation': max([transmission_graph.nodes[n].get('generation', 0)
                                   for n in transmission_graph.nodes()]) if len(transmission_graph.nodes()) > 0 else 0,
            'transmission_tree_depth': nx.dag_longest_path_length(transmission_graph) if len(transmission_graph.nodes()) > 0 else 0,
            'number_of_chains': len(list(nx.weakly_connected_components(transmission_graph))) if len(transmission_graph.nodes()) > 0 else 0
        }

    return results


def write_comprehensive_report(results, output_file):
    """
    Write comprehensive epidemiological analysis report.

    Args:
        results (dict): Analysis results
        output_file (str): Output file path
    """
    with open(output_file, 'w') as f:
        # Header
        f.write("╔" + "═" * 78 + "╗\n")
        f.write("║" + " " * 15 +
                "COMPREHENSIVE EPIDEMIOLOGICAL ANALYSIS REPORT" + " " * 15 + "║\n")
        f.write("╚" + "═" * 78 + "╝\n\n")

        f.write(f"Scenario: {results['scenario']}\n")
        f.write(f"Generated: {results['timestamp']}\n")
        f.write("─" * 80 + "\n\n")

        # Data Summary
        f.write("📊 DATA SUMMARY\n")
        f.write("─" * 20 + "\n")
        if 'data_summary' in results:
            for data_type, info in results['data_summary'].items():
                f.write(
                    f"  {data_type.replace('_', ' ').title()}: {info['records']:,} records\n")
        f.write("\n")

        # Reproduction Number Analysis
        if 'reproduction_analysis' in results and results['reproduction_analysis']:
            f.write("🦠 REPRODUCTION NUMBER ANALYSIS\n")
            f.write("─" * 35 + "\n")
            ra = results['reproduction_analysis']
            f.write(
                f"  Mean Reproduction Number (R): {ra.get('R_mean', 'N/A'):.3f}\n")
            f.write(
                f"  Median Reproduction Number: {ra.get('R_median', 'N/A'):.3f}\n")
            f.write(f"  R0 Estimate: {ra.get('R0_estimate', 'N/A'):.3f}\n")
            f.write(f"  Standard Deviation: {ra.get('R_std', 'N/A'):.3f}\n")
            f.write(
                f"  Range: {ra.get('R_min', 'N/A'):.1f} - {ra.get('R_max', 'N/A'):.1f}\n")
            f.write(
                f"  Superspreader Threshold (90th percentile): {ra.get('superspreader_threshold', 'N/A'):.1f}\n")
            f.write(
                f"  Number of Superspreaders: {ra.get('superspreader_count', 'N/A')}\n")
            f.write("\n")

        # Generation Time Analysis
        if 'generation_analysis' in results and results['generation_analysis']:
            f.write("⏱️  GENERATION TIME & SERIAL INTERVAL ANALYSIS\n")
            f.write("─" * 50 + "\n")
            ga = results['generation_analysis']
            f.write(
                f"  Mean Generation Time: {ga.get('generation_time_mean', 'N/A'):.1f} timesteps\n")
            f.write(
                f"  Median Generation Time: {ga.get('generation_time_median', 'N/A'):.1f} timesteps\n")
            f.write(
                f"  Mean Serial Interval: {ga.get('serial_interval_mean', 'N/A'):.1f} timesteps\n")
            f.write(
                f"  Serial Interval Range: {ga.get('min_serial_interval', 'N/A'):.1f} - {ga.get('max_serial_interval', 'N/A'):.1f}\n")
            f.write(
                f"  Serial Interval Std Dev: {ga.get('serial_interval_std', 'N/A'):.1f}\n")
            f.write("\n")

        # Household Analysis
        if 'household_analysis' in results and results['household_analysis']:
            f.write("🏠 HOUSEHOLD TRANSMISSION ANALYSIS\n")
            f.write("─" * 38 + "\n")
            ha = results['household_analysis']
            f.write(
                f"  Total Households: {ha.get('total_households', 'N/A'):,}\n")
            f.write(
                f"  Households with Infections: {ha.get('households_with_infections', 'N/A'):,}\n")
            f.write(
                f"  Households with Multiple Infections: {ha.get('households_with_multiple_infections', 'N/A'):,}\n")
            f.write(
                f"  Household Attack Rate: {ha.get('household_attack_rate', 'N/A'):.2f}%\n")
            f.write(
                f"  Household Secondary Transmission Rate: {ha.get('household_secondary_transmission_rate', 'N/A'):.2f}%\n")
            f.write(
                f"  Mean Household Size: {ha.get('mean_household_size', 'N/A'):.1f}\n")
            f.write(
                f"  Max Household Cluster Size: {ha.get('max_household_cluster_size', 'N/A')}\n")
            f.write("\n")

        # Location Risk Analysis
        if 'location_risk_analysis' in results and results['location_risk_analysis']:
            f.write("📍 LOCATION-BASED TRANSMISSION RISK\n")
            f.write("─" * 40 + "\n")
            for location_type, metrics in results['location_risk_analysis'].items():
                f.write(f"  {location_type}:\n")
                f.write(
                    f"    Total Infections: {metrics.get('total_infections', 'N/A')}\n")
                if 'transmission_rate_per_visit' in metrics:
                    f.write(
                        f"    Transmission Rate per Visit: {metrics['transmission_rate_per_visit']:.4f}%\n")
                    f.write(
                        f"    Attack Rate per Person: {metrics['attack_rate_per_person']:.2f}%\n")
                    f.write(
                        f"    Total Person-Visits: {metrics['total_person_visits']:,}\n")
                f.write("\n")

        # Network Analysis
        if 'network_analysis' in results and results['network_analysis']:
            f.write("🌐 CONTACT NETWORK STRUCTURE\n")
            f.write("─" * 32 + "\n")
            na = results['network_analysis']
            f.write(
                f"  Total Nodes (People): {na.get('total_nodes', 'N/A'):,}\n")
            f.write(
                f"  Total Edges (Contacts): {na.get('total_edges', 'N/A'):,}\n")
            f.write(f"  Network Density: {na.get('density', 'N/A'):.6f}\n")
            f.write(
                f"  Average Degree: {na.get('average_degree', 'N/A'):.2f}\n")
            f.write(f"  Max Degree: {na.get('max_degree', 'N/A')}\n")
            f.write(
                f"  Connected Components: {na.get('num_connected_components', 'N/A')}\n")
            f.write(
                f"  Largest Component Size: {na.get('largest_component_size', 'N/A'):,}\n")
            if 'average_clustering' in na:
                f.write(
                    f"  Average Clustering Coefficient: {na['average_clustering']:.4f}\n")
            if 'infected_mean_centrality' in na:
                f.write(
                    f"  Mean Centrality of Infected: {na['infected_mean_centrality']:.4f}\n")
            f.write("\n")

        # Epidemic Curve Analysis
        if 'epidemic_curve_analysis' in results and results['epidemic_curve_analysis']:
            f.write("📈 EPIDEMIC CURVE CHARACTERISTICS\n")
            f.write("─" * 37 + "\n")
            eca = results['epidemic_curve_analysis']
            f.write(f"  Total Duration: {eca.get('total_days', 'N/A')} days\n")
            f.write(f"  Peak Day: {eca.get('peak_day', 'N/A')}\n")
            f.write(
                f"  Peak Infections: {eca.get('peak_infections', 'N/A')}\n")
            f.write(f"  Days to Peak: {eca.get('days_to_peak', 'N/A')}\n")
            f.write(
                f"  Mean Daily Infections: {eca.get('mean_daily_infections', 'N/A'):.2f}\n")
            if 'doubling_time' in eca and not np.isnan(eca['doubling_time']):
                f.write(f"  Doubling Time: {eca['doubling_time']:.2f} days\n")
                f.write(
                    f"  Exponential Fit R²: {eca.get('exponential_fit_r_squared', 'N/A'):.3f}\n")
            f.write("\n")

        # Superspreading Analysis
        if 'superspreading_analysis' in results and results['superspreading_analysis']:
            f.write("💥 SUPERSPREADING EVENT ANALYSIS\n")
            f.write("─" * 37 + "\n")
            sa = results['superspreading_analysis']
            f.write(
                f"  Superspreader Threshold: {sa.get('superspreader_threshold', 'N/A'):.1f} secondary infections\n")
            f.write(
                f"  Number of Superspreaders: {sa.get('num_superspreaders', 'N/A')}\n")
            f.write(
                f"  Superspreader Fraction: {sa.get('superspreader_fraction', 'N/A'):.2f}%\n")
            f.write(
                f"  Infections from Superspreaders: {sa.get('infections_from_superspreaders', 'N/A')}\n")
            f.write(
                f"  Superspreader Contribution: {sa.get('superspreader_contribution', 'N/A'):.2f}%\n")
            f.write("\n")

        # Transmission Chain Analysis
        if 'transmission_chain_analysis' in results and results['transmission_chain_analysis']:
            f.write("🔗 TRANSMISSION CHAIN ANALYSIS\n")
            f.write("─" * 34 + "\n")
            tca = results['transmission_chain_analysis']
            f.write(
                f"  Total Transmission Events: {tca.get('total_transmission_events', 'N/A')}\n")
            f.write(
                f"  Maximum Generation: {tca.get('max_generation', 'N/A')}\n")
            f.write(
                f"  Transmission Tree Depth: {tca.get('transmission_tree_depth', 'N/A')}\n")
            f.write(
                f"  Number of Transmission Chains: {tca.get('number_of_chains', 'N/A')}\n")
            f.write("\n")

        # Traditional Metrics
        if 'traditional_metrics' in results and results['traditional_metrics']:
            f.write("📊 TRADITIONAL EPIDEMIOLOGICAL METRICS\n")
            f.write("─" * 43 + "\n")
            tm = results['traditional_metrics']
            f.write(
                f"  Total Population: {tm.get('total_population', 'N/A'):,}\n")
            f.write(
                f"  Total Infections: {tm.get('total_infections', 'N/A'):,}\n")
            f.write(
                f"  Overall Attack Rate: {tm.get('attack_rate', 'N/A'):.2f}%\n")
            f.write(
                f"  Total Contact Events: {tm.get('total_contacts', 'N/A'):,}\n")
            f.write("\n")

        f.write("─" * 80 + "\n")
        f.write("Report completed successfully.\n")


def calculate_trajectory_divergence(results1, results2):
    """
    Calculate divergence between two scenarios.

    Args:
        results1 (dict): Results from first scenario
        results2 (dict): Results from second scenario

    Returns:
        dict: Divergence metrics
    """
    divergence = {}

    # Compare reproduction numbers
    if ('reproduction_analysis' in results1 and 'reproduction_analysis' in results2 and
            results1['reproduction_analysis'] and results2['reproduction_analysis']):

        r1 = results1['reproduction_analysis'].get('R_mean', 0)
        r2 = results2['reproduction_analysis'].get('R_mean', 0)
        if r1 > 0 and r2 > 0:
            divergence['R_mean_divergence'] = abs(r1 - r2) / max(r1, r2)

    # Compare attack rates
    if ('traditional_metrics' in results1 and 'traditional_metrics' in results2 and
            results1['traditional_metrics'] and results2['traditional_metrics']):

        a1 = results1['traditional_metrics'].get('attack_rate', 0)
        a2 = results2['traditional_metrics'].get('attack_rate', 0)
        if a1 > 0 and a2 > 0:
            divergence['attack_rate_divergence'] = abs(a1 - a2) / max(a1, a2)

    return divergence


def write_metrics_to_file(metrics_list, output_file):
    """
    Legacy function - writes simplified metrics for backward compatibility.

    Args:
        metrics_list (list): List of metric dictionaries
        output_file (str): Output file path
    """
    # Use new comprehensive report if we have new-style results
    if len(metrics_list) > 0 and 'reproduction_analysis' in metrics_list[0]:
        write_comprehensive_report(metrics_list[0], output_file)
        return

    # Otherwise use legacy format
    with open(output_file, 'w') as f:
        f.write("Epidemiological Metrics Analysis - Legacy Format\n")
        f.write("=" * 55 + "\n")
        f.write(
            f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for metrics in metrics_list:
            if metrics is None:
                continue

            f.write(f"SCENARIO: {metrics['scenario']}\n")
            f.write("-" * 40 + "\n")
            f.write("Legacy metrics format\n\n")


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive epidemiological analysis of ABM simulation results')
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing CSV files')
    parser.add_argument('--output_dir',
                        help='Directory to save analysis reports (optional, defaults to current directory)')
    parser.add_argument('--scenario_name',
                        help='Scenario name (optional, defaults to folder name)')
    parser.add_argument('--legacy_mode', action='store_true',
                        help='Use legacy analysis mode (simplified metrics)')

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

    # Set output directory
    if args.output_dir:
        output_dir = args.output_dir
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
    else:
        output_dir = "."  # Current directory

    print(f"Processing {scenario_name} from {args.input_dir}...")
    print(f"Output will be saved to {output_dir}")

    # Perform comprehensive analysis
    results = comprehensive_epidemiological_analysis(
        args.input_dir, scenario_name)

    if results is None:
        print("Failed to process simulation data")
        return

    # Write results to the specified output directory
    output_file = os.path.join(
        output_dir, f'comprehensive_epidemic_analysis_{scenario_name}.txt')
    write_comprehensive_report(results, output_file)
    print(f"Comprehensive analysis saved to {output_file}")

    # Also create a simple summary file
    summary_file = os.path.join(
        output_dir, f'epidemic_metrics_summary_{scenario_name}.txt')
    with open(summary_file, 'w') as f:
        f.write(f"EPIDEMIC ANALYSIS SUMMARY - {scenario_name}\n")
        f.write("=" * 50 + "\n\n")

        if 'reproduction_analysis' in results and results['reproduction_analysis']:
            ra = results['reproduction_analysis']
            f.write(f"Mean R: {ra.get('R_mean', 'N/A'):.3f}\n")
            f.write(f"R0 Estimate: {ra.get('R0_estimate', 'N/A'):.3f}\n")

        if 'traditional_metrics' in results and results['traditional_metrics']:
            tm = results['traditional_metrics']
            f.write(f"Attack Rate: {tm.get('attack_rate', 'N/A'):.2f}%\n")
            f.write(
                f"Total Infections: {tm.get('total_infections', 'N/A'):,}\n")

        if 'epidemic_curve_analysis' in results and results['epidemic_curve_analysis']:
            eca = results['epidemic_curve_analysis']
            f.write(f"Peak Day: {eca.get('peak_day', 'N/A')}\n")
            f.write(f"Peak Infections: {eca.get('peak_infections', 'N/A')}\n")

    print(f"Summary saved to {summary_file}")


if __name__ == "__main__":
    main()
