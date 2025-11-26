#############################################################################
# Copyright (C) 2020-2025 MEmilio
#
# Authors: Agatha Schmidt, Henrik Zunker
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################
"""
Data generation module for GNN-based surrogate models.

This module provides functionality to generate training data for Graph Neural Network
surrogate models by running a ODE-SECIR model based simulations across multiple regions and
time periods with varying damping interventions.
"""

import copy
import os
import pickle
import random
import time
from enum import Enum
from typing import Dict, List, Tuple

import numpy as np
from progress.bar import Bar

import memilio.simulation as mio
import memilio.simulation.osecir as osecir
from memilio.simulation import AgeGroup, set_log_level, Damping
from memilio.simulation.osecir import (
    Index_InfectionState, InfectionState, GraphParameterStudy,
    interpolate_simulation_result)

from memilio.surrogatemodel.GNN.GNN_utils import (
    scale_data,
    remove_confirmed_compartments
)
import memilio.surrogatemodel.utils.dampings as dampings


class Location(Enum):
    """Contact location types for the model."""
    Home = 0
    School = 1
    Work = 2
    Other = 3


# Default number of age groups
DEFAULT_NUM_AGE_GROUPS = 6

# Contact location names corresponding to provided contact files
CONTACT_LOCATIONS = ["home", "school_pf_eig", "work", "other"]


def set_covid_parameters(
        model, start_date, num_groups=DEFAULT_NUM_AGE_GROUPS):
    """Sets COVID-19 specific parameters for all age groups.

    Parameters are based on Kühn et al. (2021): https://doi.org/10.1016/j.mbs.2021.108648
    The function sets age-stratified parameters (when available) the following age groups:
    0-4, 5-14, 15-34, 35-59, 60-79, 80+

    :param model: MEmilio ODE SECIR model to configure.
    :param start_date: Start date of the simulation (used to set StartDay parameter).
    :param num_groups: Number of age groups (default: 6).

    """
    # Age-specific transmission and progression parameters
    transmission_probability = [0.03, 0.06, 0.06, 0.06, 0.09, 0.175]
    recovered_per_infected_no_symptoms = [0.25, 0.25, 0.2, 0.2, 0.2, 0.2]
    severe_per_infected_symptoms = [
        0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225]
    critical_per_severe = [0.075, 0.075, 0.075, 0.15, 0.3, 0.4]
    deaths_per_critical = [0.05, 0.05, 0.14, 0.14, 0.4, 0.6]

    # Age-specific compartment transition times (in days)
    time_infected_no_symptoms = [2.74, 2.74, 2.565, 2.565, 2.565, 2.565]
    time_infected_symptoms = [7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775]
    time_infected_severe = [5, 5, 5.925, 7.55, 8.5, 11]
    time_infected_critical = [6.95, 6.95, 6.86, 17.36, 17.1, 11.6]

    # Apply parameters for each age group
    for i in range(num_groups):
        age_group = AgeGroup(i)

        model.parameters.TimeExposed[age_group] = 3.335
        model.parameters.TimeInfectedNoSymptoms[age_group] = time_infected_no_symptoms[i]
        model.parameters.TimeInfectedSymptoms[age_group] = time_infected_symptoms[i]
        model.parameters.TimeInfectedSevere[age_group] = time_infected_severe[i]
        model.parameters.TimeInfectedCritical[age_group] = time_infected_critical[i]

        model.parameters.RelativeTransmissionNoSymptoms[age_group] = 1.0
        model.parameters.TransmissionProbabilityOnContact[age_group] = transmission_probability[i]
        model.parameters.RiskOfInfectionFromSymptomatic[age_group] = 0.25
        model.parameters.MaxRiskOfInfectionFromSymptomatic[age_group] = 0.5
        model.parameters.RecoveredPerInfectedNoSymptoms[age_group] = recovered_per_infected_no_symptoms[i]
        model.parameters.SeverePerInfectedSymptoms[age_group] = severe_per_infected_symptoms[i]
        model.parameters.CriticalPerSevere[age_group] = critical_per_severe[i]
        model.parameters.DeathsPerCritical[age_group] = deaths_per_critical[i]

    # Set simulation start day
    model.parameters.StartDay = start_date.day_in_year


def set_contact_matrices(model, data_dir, num_groups=DEFAULT_NUM_AGE_GROUPS):
    """Loads and configures contact matrices for different location types.

    Contact matrices are loaded for the locations defined in CONTACT_LOCATIONS.

    :param model: MEmilio ODE SECIR model to configure.
    :param data_dir: Root directory containing contact matrix data (should contain Germany/contacts/ subdirectory).
    :param num_groups: Number of age groups (default: 6).

    """
    contact_matrices = mio.ContactMatrixGroup(
        len(CONTACT_LOCATIONS), num_groups)

    # Load contact matrices for each location
    for location_idx, location_name in enumerate(CONTACT_LOCATIONS):
        contact_file = os.path.join(
            data_dir, "Germany", "contacts", f"baseline_{location_name}.txt"
        )

        if not os.path.exists(contact_file):
            raise FileNotFoundError(
                f"Contact matrix file not found: {contact_file}"
            )

        contact_matrices[location_idx] = mio.ContactMatrix(
            mio.read_mobility_plain(contact_file)
        )

    model.parameters.ContactPatterns.cont_freq_mat = contact_matrices


def get_graph(num_groups, data_dir, mobility_directory, start_date, end_date):
    """Creates a graph with mobility connections.

    Creates a graph where each node represents a geographic region (here county) with its own ODE model, 
    and edges represent mobility/commuter connections between these regions. Each node is initialized with 
    population data and COVID parameters. Edges are weighted by commuter mobility patterns. 

    :param num_groups: Number of age groups to model.
    :param data_dir: Root directory containing population and contact data.
    :param mobility_directory: Path to mobility/commuter data file.
    :param start_date: Simulation start date.
    :param end_date: Simulation end date (used for data loading).
    :returns: Configured ModelGraph with nodes for each region and mobility edges.

    """
    model = osecir.Model(num_groups)
    set_covid_parameters(model, start_date, num_groups)
    set_contact_matrices(model, data_dir, num_groups)

    # Initialize empty graph
    graph = osecir.ModelGraph()

    # To account for underreporting
    scaling_factor_infected = [2.5] * num_groups
    scaling_factor_icu = 1.0
    # Test & trace capacity as in Kühn et al. (2021): https://doi.org/10.1016/j.mbs.2021.108648
    tnt_capacity_factor = 7.5 / 100000.0

    # Paths to population data
    pydata_dir = os.path.join(data_dir, "Germany", "pydata")
    path_population_data = os.path.join(
        pydata_dir, "county_current_population.json")

    # Verify data files exist
    if not os.path.exists(path_population_data):
        raise FileNotFoundError(
            f"Population data not found: {path_population_data}"
        )

    # Create one node per county
    is_node_for_county = True

    # Populate graph nodes with county level data
    osecir.set_nodes(
        model.parameters,
        start_date,
        end_date,
        pydata_dir,
        path_population_data,
        is_node_for_county,
        graph,
        scaling_factor_infected,
        scaling_factor_icu,
        tnt_capacity_factor,
        0,  # Zero days extrapolating the data
        False  # No extrapolation of data
    )

    # Add mobility edges between regions
    osecir.set_edges(mobility_directory, graph, len(CONTACT_LOCATIONS))

    return graph


def get_compartment_factors():
    """Draw base factors for compartment initialization.

    Factors follow the sampling strategy from Schmidt et al.
    (2024): symptomatic individuals between 0.01% and 5% of the population,
    exposed and asymptomatic compartments proportional to that symptomatic
    proportion, and hospital/ICU/deaths sampled hierarchically.
    Recovered is taken from the remaining feasible proportion and susceptibles
    fill the residual.
    """
    p_infected = random.uniform(0.0001, 0.05)
    p_exposed = p_infected * random.uniform(0.1, 5.0)
    p_ins = p_infected * random.uniform(0.1, 5.0)
    p_hosp = p_infected * random.uniform(0.001, 1.0)
    p_icu = p_hosp * random.uniform(0.001, 1.0)
    p_dead = p_icu * random.uniform(0.001, 1.0)

    sum_randoms = (
        p_infected + p_exposed + p_ins + p_hosp + p_icu + p_dead
    )

    if sum_randoms >= 1.0:
        raise RuntimeError(
            "Sampled compartment factors exceed total population. Adjust bounds."
        )

    p_recovered = random.uniform(0.0, 1.0 - sum_randoms)
    p_susceptible = max(1.0 - (sum_randoms + p_recovered), 0.0)

    return {
        "infected": p_infected,
        "exposed": p_exposed,
        "infected_no_symptoms": p_ins,
        "hospitalized": p_hosp,
        "critical": p_icu,
        "dead": p_dead,
        "recovered": p_recovered,
        "susceptible": p_susceptible
    }


def _initialize_compartments_for_node(
        model, factors, num_groups, within_group_variation):
    """Initializes epidemic compartments using shared base factors.

    :param model: Model instance for the specific node.
    :param factors: Compartment factors obtained from get_compartment_factors.
    :param num_groups: Number of age groups.
    :param within_group_variation: Whether to apply additional random scaling per age group.
    """

    def _variation():
        return random.uniform(0.1, 1.0) if within_group_variation else 1.0

    for age_idx in range(num_groups):
        age_group = AgeGroup(age_idx)
        total_population = model.populations.get_group_total_AgeGroup(
            age_group)

        infected_symptoms = total_population * \
            factors["infected"] * _variation()

        exposed = infected_symptoms * factors["exposed"] * _variation()
        infected_no_symptoms = infected_symptoms * \
            factors["infected_no_symptoms"] * _variation()
        infected_severe = infected_symptoms * \
            factors["hospitalized"] * _variation()
        infected_critical = infected_severe * \
            factors["critical"] * _variation()
        dead = infected_critical * factors["dead"] * _variation()
        recovered = total_population * factors["recovered"] * _variation()

        model.populations[age_group, Index_InfectionState(
            InfectionState.InfectedSymptoms)] = infected_symptoms
        model.populations[age_group, Index_InfectionState(
            InfectionState.Exposed)] = exposed
        model.populations[age_group, Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = infected_no_symptoms
        model.populations[age_group, Index_InfectionState(
            InfectionState.InfectedSevere)] = infected_severe
        model.populations[age_group, Index_InfectionState(
            InfectionState.InfectedCritical)] = infected_critical
        model.populations[age_group, Index_InfectionState(
            InfectionState.Dead)] = dead
        model.populations[age_group, Index_InfectionState(
            InfectionState.Recovered)] = recovered

        # Susceptibles are the remainder
        model.populations.set_difference_from_group_total_AgeGroup(
            (age_group, InfectionState.Susceptible), total_population
        )


def _apply_dampings_to_model(model, damping_days, damping_factors, num_groups):
    """Applies contact dampings (NPIs) to model at specified days.

    Currently only supports global dampings (same for all age groups and spatial units).

    :param model: Model to apply dampings to.
    :param damping_days: Days at which to apply dampings.
    :param damping_factors: Multiplicative factors for contact reduction.
    :param num_groups: Number of age groups.
    :returns: Tuple of (damped_contact_matrices, damping_coefficients).

    """
    damped_matrices = []
    damping_coefficients = []

    for day, factor in zip(damping_days, damping_factors):
        # Create uniform damping matrix for all age groups
        damping_matrix = np.ones((num_groups, num_groups)) * factor

        # Add damping to model
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(
            Damping(coeffs=damping_matrix, t=day, level=0, type=0)
        )

        # Store resulting contact matrix and coefficients
        damped_matrices.append(
            model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
                day + 1))
        damping_coefficients.append(damping_matrix)

    return damped_matrices, damping_coefficients


def run_secir_groups_simulation(
        days, damping_days, damping_factors, graph, within_group_variation,
        num_groups=DEFAULT_NUM_AGE_GROUPS):
    """Runs a multi-region ODE SECIR simulation with age groups and NPIs.

    Performs the following steps:
    1. Initialize each node with compartment specific (random) factors
    2. Apply contact dampings (NPIs) at specified days
    3. Run the simulation
    4. Post-process and return results

    :param days: Total number of days to simulate.
    :param damping_days: List of days when NPIs are applied.
    :param damping_factors: List of contact reduction factors for each damping.
    :param graph: Pre-configured ModelGraph with nodes and edges.
    :param within_group_variation: Whether to apply per-age random variation during initialization.
    :param num_groups: Number of age groups (default: 6).
    :returns: Tuple containing dataset_entry (simulation results for each node [num_nodes, time_steps, compartments]), 
              damped_matrices (contact matrices at each damping day), damping_coefficients (damping coefficient matrices), 
              and runtime (execution time in seconds).
    :raises ValueError: If lengths of damping_days and damping_factors don't match.

    """
    if len(damping_days) != len(damping_factors):
        raise ValueError(
            f"Length mismatch: damping_days ({len(damping_days)}) != "
            f"damping_factors ({len(damping_factors)})"
        )

    # Initialize each node in the graph
    damped_matrices = []
    damping_coefficients = []

    factors = get_compartment_factors()

    for node_idx in range(graph.num_nodes):
        model = graph.get_node(node_idx).property

        # Initialize compartment populations
        _initialize_compartments_for_node(
            model, factors, num_groups, within_group_variation)

        # Apply dampings/NPIs
        if damping_days:
            node_damped_mats, node_damping_coeffs = _apply_dampings_to_model(
                model, damping_days, damping_factors, num_groups
            )
            # Store only from first node, since dampings are global at the moment
            if node_idx == 0:
                damped_matrices = node_damped_mats
                damping_coefficients = node_damping_coeffs

        model.apply_constraints()

        # Update graph with initialized populations
        graph.get_node(node_idx).property.populations = model.populations

    # Run simulation and measure runtime
    study = GraphParameterStudy(graph, 0, days, 0.5, 1)

    start_time = time.perf_counter()
    study_results = study.run()
    runtime = time.perf_counter() - start_time

    # Interpolate results to daily values
    graph_run = study_results[0]
    results = interpolate_simulation_result(graph_run)

    # Remove confirmed compartments (not used in GNN)
    for result_idx in range(len(results)):
        results[result_idx] = remove_confirmed_compartments(
            np.asarray(results[result_idx]), num_groups
        )

    dataset_entry = copy.deepcopy(results)

    return dataset_entry, damped_matrices, damping_coefficients, runtime


def generate_data(num_runs, data_dir, output_path, input_width, label_width,
                  start_date, end_date, save_data=True, transform=True,
                  damping_method="classic", max_number_damping=3,
                  mobility_file="commuter_mobility.txt",
                  num_groups=DEFAULT_NUM_AGE_GROUPS,
                  within_group_variation=True):
    """Generates training dataset for GNN surrogate model.

    Runs num_runs-ODE SECIR simulations with random initial conditions and damping patterns to create a training dataset.
    If save_data=True, saves a pickle file named:
    'GNN_data_{label_width}days_{max_number_damping}dampings_{damping_method}{num_runs}.pickle'

    :param num_runs: Number of simulation runs to generate.
    :param data_dir: Root directory containing all required input data (population, contacts, mobility).
    :param output_path: Directory where generated dataset will be saved.
    :param input_width: Number of time steps for model input.
    :param label_width: Number of time steps for model output/labels.
    :param start_date: Simulation start date.
    :param end_date: Simulation end date (used for set_nodes function).
    :param save_data: Whether to save dataset (default: True).
    :param transform: Whether to apply scaling transformation to data (default: True).
    :param damping_method: Method for generating damping patterns: "classic", "active", "random".
    :param max_number_damping: Maximum number of damping events per simulation.
    :param mobility_file: Filename of mobility file (in data_dir/Germany/mobility/).
    :param num_groups: Number of age groups (default: 6).
    :param within_group_variation: Whether to apply random variation per spatial unit/age group when initializing compartments.
    :returns: Dictionary with keys: "inputs" ([num_runs, input_width, num_nodes, features]), 
              "labels" ([num_runs, label_width, num_nodes, features]), 
              "contact_matrix" (List of damped contact matrices), "damping_days" (List of damping day arrays), 
              "damping_factors" (List of damping factor arrays).

    """
    set_log_level(mio.LogLevel.Error)

    # Calculate total simulation days
    total_days = label_width + input_width - 1

    # Initialize output dictionary
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": [],
        "damping_days": [],
        "damping_factors": []
    }

    # Build mobility file path
    mobility_path = os.path.join(
        data_dir, "Germany", "mobility", mobility_file)

    # Verify mobility file exists
    if not os.path.exists(mobility_path):
        raise FileNotFoundError(f"Mobility file not found: {mobility_path}")

    # Create graph (reused for all runs with different initial conditions)
    graph = get_graph(num_groups, data_dir,
                      mobility_path, start_date, end_date)

    print(f"\nGenerating {num_runs} simulation runs...")
    bar = Bar(
        'Progress', max=num_runs,
        suffix='%(percent)d%% [%(elapsed_td)s / %(eta_td)s]')

    runtimes = []

    for _ in range(num_runs):
        # Generate random damping pattern
        if max_number_damping > 0:
            damping_days, damping_factors = dampings.generate_dampings(
                total_days,
                max_number_damping,
                method=damping_method,
                min_distance=2,
                min_damping_day=2
            )
        else:
            damping_days = []
            damping_factors = []

        # Run simulation
        simulation_result, damped_mats, damping_coeffs, runtime = \
            run_secir_groups_simulation(
                total_days, damping_days, damping_factors, graph,
                within_group_variation, num_groups
            )

        runtimes.append(runtime)

        # Split into inputs and labels
        # Shape: [num_nodes, time_steps, features] -> transpose to [time_steps, num_nodes, features]
        result_transposed = np.asarray(simulation_result).transpose(1, 0, 2)
        inputs = result_transposed[:input_width]
        labels = result_transposed[input_width:]

        # Store results
        data["inputs"].append(inputs)
        data["labels"].append(labels)
        data["contact_matrix"].append(np.array(damped_mats))
        data["damping_factors"].append(damping_coeffs)
        data["damping_days"].append(damping_days)

        bar.next()

    bar.finish()

    # Print performance statistics
    print(f"\nSimulation Statistics:")
    print(f"  Total days simulated: {total_days}")
    print(f"  Average runtime: {np.mean(runtimes):.3f}s")
    print(f"  Median runtime: {np.median(runtimes):.3f}s")
    print(f"  Total time: {np.sum(runtimes):.1f}s")

    # Save dataset if requested
    if save_data:
        # Apply scaling transformation
        inputs_scaled, labels_scaled = scale_data(data, transform)

        all_data = {
            "inputs": inputs_scaled,
            "labels": labels_scaled,
            "damping_day": data["damping_days"],
            "contact_matrix": data["contact_matrix"],
            "damping_coeff": data["damping_factors"]
        }

        # Create output directory if needed
        os.makedirs(output_path, exist_ok=True)

        # Generate filename
        if num_runs < 1000:
            filename = f'GNN_data_{label_width}days_{max_number_damping}dampings_{damping_method}{num_runs}.pickle'
        else:
            filename = f'GNN_data_{label_width}days_{max_number_damping}dampings_{damping_method}{num_runs//1000}k.pickle'

        # Save to pickle file
        output_file = os.path.join(output_path, filename)
        with open(output_file, 'wb') as f:
            pickle.dump(all_data, f)

        print(f"\nDataset saved to: {output_file}")

    return data


def main():
    """Main function for dataset generation.

    Example configuration for generating GNN training data.

    """
    # Set random seed for reproducibility
    random.seed(10)

    # Configuration
    data_dir = os.path.join(os.getcwd(), 'data')
    output_path = os.path.join(os.getcwd(), 'generated_datasets')

    # Simulation parameters
    input_width = 5  # Days of history used as input
    num_runs = 1  # Number of simulation runs
    max_dampings = 0  # Number of NPI dampings per simulation

    # Prediction horizons to generate data for
    prediction_horizons = [30]  # Days to predict into the future

    # Simulation time period
    start_date = mio.Date(2020, 10, 1)
    end_date = mio.Date(2021, 10, 31)

    # Generate datasets
    print("=" * 70)
    print("GNN Surrogate Model - Dataset Generation")
    print("=" * 70)
    print(f"Data directory: {data_dir}")
    print(f"Output directory: {output_path}")
    print(f"Simulation period: {start_date} to {end_date}")
    print(f"Number of runs per configuration: {num_runs}")
    print(f"Input width: {input_width} days")
    print(f"Max dampings: {max_dampings}")
    print("=" * 70)

    for label_width in prediction_horizons:
        print(f"\n{'='*70}")
        print(f"Generating data for {label_width}-day predictions")
        print(f"{'='*70}")

        generate_data(
            num_runs=num_runs,
            data_dir=data_dir,
            output_path=output_path,
            input_width=input_width,
            label_width=label_width,
            start_date=start_date,
            end_date=end_date,
            save_data=True,
            damping_method="active",
            max_number_damping=max_dampings
        )

    print(f"\n{'='*70}")
    print("Dataset generation complete!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
