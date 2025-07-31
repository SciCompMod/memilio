#!/usr/bin/env python3
"""
Epidemiological Metrics Calculator for ABM Simulation Results

This script calculates 5 key epidemiological metrics from agent-based model simulation data:
1. Basic Reproduction Number (R0)
2. Time to 100 Infections
3. Peak Infection Day
4. Attack Rate
5. Trajectory Divergence Score

Input: HDF5 files containing simulation results
Output: Text file with calculated metrics
"""

import h5py
import numpy as np
import pandas as pd
import os
import sys
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.stats import linregress
import argparse


def load_h5_results(file_path, dataset_key='0'):
    """
    Load simulation results from HDF5 file.

    Args:
        file_path (str): Path to HDF5 results file
        dataset_key (str): Key for dataset in HDF5 file (default '0')

    Returns:
        dict: Dictionary containing time series data
    """
    try:
        with h5py.File(file_path, 'r') as f:
            data = {k: v[()] for k, v in f[dataset_key].items()}
        return data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None


def calculate_daily_infections(infection_state_data, time_data):
    """
    Calculate daily new infections from infection state time series.

    Args:
        infection_state_data (np.array): Infection state data over time
        time_data (np.array): Time points

    Returns:
        tuple: (time_array, daily_infections_array)
    """
    # Sum all infected states (Exposed, I_Asymp, I_Symp, I_Severe, I_Critical)
    # Assuming states 1-5 are infected states based on your code structure
    # Exposed, I_Asymp, I_Symp, I_Severe, I_Critical
    infected_states = [1, 2, 3, 4, 5]

    total_infected = np.sum(infection_state_data[:, infected_states], axis=1)

    # Calculate daily new infections (difference between consecutive time points)
    daily_new = np.diff(total_infected)
    daily_new = np.maximum(daily_new, 0)  # Ensure non-negative

    return time_data[1:], daily_new


def calculate_r0(time_data, daily_infections, generation_time=5.0):
    """
    Calculate basic reproduction number R0 from early exponential growth.

    Args:
        time_data (np.array): Time points
        daily_infections (np.array): Daily new infections
        generation_time (float): Average generation time in days

    Returns:
        float: Estimated R0 value
    """
    # Find early exponential growth phase (first 14 days with >5 infections)
    early_mask = (time_data <= 14) & (daily_infections >= 5)

    if np.sum(early_mask) < 5:  # Need at least 5 points
        return np.nan

    early_time = time_data[early_mask]
    early_infections = daily_infections[early_mask]

    # Fit exponential growth: log(infections) = r*t + c
    try:
        log_infections = np.log(early_infections + 1)  # Add 1 to avoid log(0)
        slope, intercept, r_value, p_value, std_err = linregress(
            early_time, log_infections)

        # R0 = 1 + r * generation_time
        r0 = 1 + slope * generation_time
        return max(r0, 0)  # Ensure non-negative
    except:
        return np.nan


def calculate_time_to_threshold(time_data, cumulative_infections, threshold=100):
    """
    Calculate time to reach infection threshold.

    Args:
        time_data (np.array): Time points
        cumulative_infections (np.array): Cumulative infection count
        threshold (int): Infection threshold

    Returns:
        float: Days to reach threshold (NaN if never reached)
    """
    threshold_idx = np.where(cumulative_infections >= threshold)[0]

    if len(threshold_idx) == 0:
        return np.nan

    return time_data[threshold_idx[0]]


def calculate_peak_day(time_data, daily_infections):
    """
    Calculate peak infection day and height.

    Args:
        time_data (np.array): Time points
        daily_infections (np.array): Daily new infections

    Returns:
        tuple: (peak_day, peak_height)
    """
    peak_idx = np.argmax(daily_infections)
    peak_day = time_data[peak_idx]
    peak_height = daily_infections[peak_idx]

    return peak_day, peak_height


def calculate_attack_rate(infection_state_data, total_population):
    """
    Calculate final attack rate.

    Args:
        infection_state_data (np.array): Infection state data over time
        total_population (int): Total population size

    Returns:
        float: Attack rate as percentage
    """
    # Get final time point data
    final_data = infection_state_data[-1, :]

    # Sum recovered and dead (assuming these are in the final states)
    # States: 0=Susceptible, 1=Exposed, 2=I_Asymp, 3=I_Symp, 4=I_Severe, 5=I_Critical, 6=Recovered, 7=Dead
    total_ever_infected = np.sum(final_data[1:])  # All non-susceptible

    attack_rate = (total_ever_infected / total_population) * 100
    return attack_rate


def calculate_trajectory_divergence(model1_infections, model2_infections):
    """
    Calculate Root Mean Square Error between two infection trajectories.

    Args:
        model1_infections (np.array): Infection trajectory from model 1
        model2_infections (np.array): Infection trajectory from model 2

    Returns:
        float: RMSE divergence score
    """
    # Ensure same length
    min_len = min(len(model1_infections), len(model2_infections))
    traj1 = model1_infections[:min_len]
    traj2 = model2_infections[:min_len]

    # Calculate RMSE
    rmse = np.sqrt(np.mean((traj1 - traj2) ** 2))

    # Normalize by peak infection count for comparability
    peak_max = max(np.max(traj1), np.max(traj2))
    if peak_max > 0:
        normalized_rmse = rmse / peak_max
    else:
        normalized_rmse = 0

    return normalized_rmse


def process_simulation_data(data_path, scenario_name, total_population=10000):
    """
    Process simulation data and calculate all metrics.

    Args:
        data_path (str): Path to HDF5 file
        scenario_name (str): Name of scenario for output
        total_population (int): Total population size

    Returns:
        dict: Dictionary containing all calculated metrics
    """
    # Load data
    data = load_h5_results(data_path)
    if data is None:
        return None

    time_data = data['Time']

    # Assuming 'Total' contains infection state data
    # Shape should be (time_points, infection_states)
    if 'Total' in data:
        infection_state_data = data['Total']
    else:
        print(f"Warning: 'Total' key not found in {data_path}")
        return None

    # Calculate daily infections
    time_daily, daily_infections = calculate_daily_infections(
        infection_state_data, time_data)
    cumulative_infections = np.cumsum(daily_infections)

    # Calculate all metrics
    metrics = {
        'scenario': scenario_name,
        'total_population': total_population,
        'simulation_days': len(time_data),
    }

    # 1. Basic Reproduction Number
    metrics['R0'] = calculate_r0(time_daily, daily_infections)

    # 2. Time to 100 infections
    metrics['time_to_100'] = calculate_time_to_threshold(
        time_daily, cumulative_infections, 100)

    # 3. Peak infection day and height
    peak_day, peak_height = calculate_peak_day(time_daily, daily_infections)
    metrics['peak_day'] = peak_day
    metrics['peak_height'] = peak_height

    # 4. Attack rate
    metrics['attack_rate'] = calculate_attack_rate(
        infection_state_data, total_population)

    return metrics


def write_metrics_to_file(metrics_list, output_file):
    """
    Write calculated metrics to text file.

    Args:
        metrics_list (list): List of metric dictionaries
        output_file (str): Output file path
    """
    with open(output_file, 'w') as f:
        # Write header
        f.write("Epidemiological Metrics Analysis\n")
        f.write("=" * 50 + "\n")
        f.write(
            f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Write metrics for each scenario
        for i, metrics in enumerate(metrics_list):
            if metrics is None:
                continue

            f.write(f"SCENARIO {i+1}: {metrics['scenario']}\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total Population: {metrics['total_population']:,}\n")
            f.write(f"Simulation Days: {metrics['simulation_days']}\n\n")

            f.write("KEY METRICS:\n")
            f.write(
                f"1. Basic Reproduction Number (R0): {metrics['R0']:.3f}\n")

            if np.isnan(metrics['time_to_100']):
                f.write("2. Time to 100 Infections: Never reached\n")
            else:
                f.write(
                    f"2. Time to 100 Infections: {metrics['time_to_100']:.1f} days\n")

            f.write(f"3. Peak Infection Day: {metrics['peak_day']:.1f} days\n")
            f.write(
                f"   Peak Daily Infections: {metrics['peak_height']:.0f}\n")
            f.write(f"4. Attack Rate: {metrics['attack_rate']:.2f}%\n")

            f.write("\n" + "=" * 50 + "\n\n")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate epidemiological metrics from ABM simulation results')
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing HDF5 result files')
    parser.add_argument(
        '--output_file', default='epidemic_metrics.txt', help='Output text file')
    parser.add_argument('--population', type=int,
                        default=10000, help='Total population size')
    parser.add_argument('--scenarios', nargs='+',
                        help='Scenario names (optional)')

    args = parser.parse_args()

    # Find HDF5 files in input directory
    h5_files = []
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.endswith('.h5') or file.endswith('.hdf5'):
                h5_files.append(os.path.join(root, file))

    if not h5_files:
        print(f"No HDF5 files found in {args.input_dir}")
        return

    print(f"Found {len(h5_files)} HDF5 files")

    # Process each file
    all_metrics = []
    for i, file_path in enumerate(h5_files):
        if args.scenarios and i < len(args.scenarios):
            scenario_name = args.scenarios[i]
        else:
            scenario_name = os.path.basename(file_path).replace(
                '.h5', '').replace('.hdf5', '')

        print(f"Processing {scenario_name}...")
        metrics = process_simulation_data(
            file_path, scenario_name, args.population)
        all_metrics.append(metrics)

    # Write results
    write_metrics_to_file(all_metrics, args.output_file)
    print(f"Metrics saved to {args.output_file}")


if __name__ == "__main__":
    main()
