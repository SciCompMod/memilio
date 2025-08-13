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

import copy
import os
import pickle
import random
import json
from datetime import date

import numpy as np
import tensorflow as tf
from progress.bar import Bar
from sklearn.preprocessing import FunctionTransformer

from memilio.simulation import (AgeGroup, Damping, LogLevel, set_log_level)
from memilio.simulation.osecir import (Index_InfectionState,
                                       InfectionState, Model,
                                       interpolate_simulation_result, simulate)
import memilio.surrogatemodel.utils.dampings as dampings
from memilio.surrogatemodel.utils.helper_functions import (
    interpolate_age_groups, remove_confirmed_compartments, normalize_simulation_data)
import memilio.simulation as mio
import memilio.simulation.osecir as osecir


def run_secir_groups_simulation(days, damping_days, damping_factors, populations):
    """ Uses an ODE SECIR model allowing for asymptomatic infection with 6 different age groups. The model is not stratified by region.
    Virus-specific parameters are fixed and initial number of persons in the particular infection states are chosen randomly from defined ranges.

    :param days: Describes how many days we simulate within a single run.
    :param damping_days: The days when damping is applied.
    :param damping_factors: damping factors associated to the damping days.
    :param populations: List containing the population in each age group.
    :returns: Tuple of lists (list_of_simulation_results, list_of_damped_matrices), the first containing the simulation results, the second list containing the 
        damped contact matrices. 

    """
    # Collect indices of confirmed compartments
    del_indices = []

    if len(damping_days) != len(damping_factors):
        raise ValueError("Length of damping_days and damping_factors differ!")

    set_log_level(LogLevel.Off)

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1

    age_groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    # Get number of age groups
    num_groups = len(age_groups)

    # Initialize Parameters
    model = Model(num_groups)

    # Set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.TimeExposed[AgeGroup(i)] = 3.2
        model.parameters.TimeInfectedNoSymptoms[AgeGroup(i)] = 2.
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

        # Initial number of people in each compartment with random numbers
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = random.uniform(
            0.00025, 0.005) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = random.uniform(
            0.0001, 0.0035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptomsConfirmed)] = 0
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = random.uniform(
            0.00007, 0.001) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptomsConfirmed)] = 0
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = random.uniform(
            0.00003, 0.0006) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = random.uniform(
            0.00001, 0.0002) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(
            0.002, 0.08) * populations[i]
        model.populations[AgeGroup(i),
                          Index_InfectionState(InfectionState.Dead)] = 0
        model.populations.set_difference_from_group_total_AgeGroup(
            (AgeGroup(i), Index_InfectionState(InfectionState.Susceptible)), populations[i])

        # Compartment transition propabilities
        model.parameters.RelativeTransmissionNoSymptoms[AgeGroup(i)] = 0.5
        model.parameters.TransmissionProbabilityOnContact[AgeGroup(i)] = 0.1
        model.parameters.RecoveredPerInfectedNoSymptoms[AgeGroup(i)] = 0.09
        model.parameters.RiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.25
        model.parameters.SeverePerInfectedSymptoms[AgeGroup(i)] = 0.2
        model.parameters.CriticalPerSevere[AgeGroup(i)] = 0.25
        model.parameters.DeathsPerCritical[AgeGroup(i)] = 0.3
        # twice the value of RiskOfInfectionFromSymptomatic
        model.parameters.MaxRiskOfInfectionFromSymptomatic[AgeGroup(i)] = 0.5

        # Collecting deletable indices
        index_no_sym_conf = model.populations.get_flat_index(
            osecir.MultiIndex_PopulationsArray(mio.AgeGroup(i), osecir.InfectionState.InfectedNoSymptomsConfirmed))
        index_sym_conf = model.populations.get_flat_index(
            osecir.MultiIndex_PopulationsArray(mio.AgeGroup(i), osecir.InfectionState.InfectedSymptomsConfirmed))
        del_indices.append(index_no_sym_conf)
        del_indices.append(index_sym_conf)

    # StartDay is the n-th day of the year
    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # Load baseline and minimum contact matrix and assign them to the model
    baseline = getBaselineMatrix()
    minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum

    # Generate a damping matrix and assign it to the model
    damped_matrices = []

    for i in np.arange(len(damping_days)):
        damping = np.ones((num_groups, num_groups)
                          ) * damping_factors[i]
        day = damping_days[i]
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=(damping), t=day, level=0, type=0))
        damped_matrices.append(model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
            day+1))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)

    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    # Omit first column, as the time points are not of interest here.
    result_array = remove_confirmed_compartments(
        np.transpose(result.as_ndarray()[1:, :]), del_indices)

    dataset_entries = copy.deepcopy(result_array)
    dataset_entries = np.nan_to_num(dataset_entries)

    return dataset_entries.tolist(), damped_matrices


def generate_data(
        num_runs, path_out, path_population, input_width, label_width,
        normalize=True, save_data=True, damping_method="random", number_dampings=5):
    """ Generate data sets of num_runs many equation-based model simulations and possibly transforms the computed results by a log(1+x) transformation.
    Divides the results in input and label data sets and returns them as a dictionary of two TensorFlow Stacks.
    In general, we have 8 different compartments and 6 age groups.  If we choose
    input_width = 5 and label_width = 20, the dataset has
    - input with dimension 5 x 8 x 6
    - labels with dimension 20 x 8 x 6

    :param num_runs: Number of times, the function run_secir_groups_simulation is called.
    :param path_out: Path, where the dataset is saved to.
    :param path_population: Path, where we try to read the population data.
    :param input_width: Int value that defines the number of time series used for the input.
    :param label_width: Int value that defines the size of the labels.
    :param normalize: Default: true Option to transform dataset by logarithmic normalization.
    :param save_data: Default: true Option to save the dataset.
    :param damping_method: String specifying the damping method, that should be used. Possible values "classic", "active", "random".
    :param number_dampings: Maximal number of possible dampings. 
    :returns: Data dictionary of input and label data sets.

    """
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrices": [],
        "damping_factors": [],
        "damping_days": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1

    # Load population data
    population = get_population(path_population)

    # show progess in terminal for longer runs
    # Due to the random structure, there's currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):

        # Generate random damping days
        damping_days, damping_factors = dampings.generate_dampings(
            days, number_dampings, method=damping_method, min_distance=2,
            min_damping_day=2)

        data_run, damped_matrices = run_secir_groups_simulation(
            days, damping_days, damping_factors, population[random.randint(0, len(population) - 1)])
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        data['contact_matrices'].append(damped_matrices)
        data['damping_factors'].append(damping_factors)
        data['damping_days'].append(damping_days)
        bar.next()
    bar.finish()

    if normalize:
        # logarithmic normalization
        transformer = FunctionTransformer(np.log1p, validate=True)

        # transform inputs and labels
        data['inputs'] = normalize_simulation_data(
            data['inputs'], transformer, num_runs)
        data['labels'] = normalize_simulation_data(
            data['labels'], transformer, num_runs)
    else:
        data['inputs'] = tf.convert_to_tensor(data['inputs'])
        data['labels'] = tf.convert_to_tensor(data['labels'])

    if save_data:
        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path_out):
            os.mkdir(path_out)

        # save dict to pickle file
        if num_runs < 1000:
            filename = 'data_secir_groups_%ddays_%d_' % (
                label_width, num_runs) + damping_method+'.pickle'
        else:
            filename = 'data_secir_groups_%ddays_%dk_' % (
                label_width, num_runs//1000) + damping_method+'.pickle'

        with open(os.path.join(path_out, filename), 'wb') as f:
            pickle.dump(data, f)
    return data


def getBaselineMatrix():
    """ loads the baselinematrix"""

    baseline_contact_matrix0 = os.path.join(
        "./data/Germany/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/Germany/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/Germany/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/Germany/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


def getMinimumMatrix():
    """ loads the minimum matrix"""

    minimum_contact_matrix0 = os.path.join(
        "./data/Germany/contacts/minimum_home.txt")
    minimum_contact_matrix1 = os.path.join(
        "./data/Germany/contacts/minimum_school_pf_eig.txt")
    minimum_contact_matrix2 = os.path.join(
        "./data/Germany/contacts/minimum_work.txt")
    minimum_contact_matrix3 = os.path.join(
        "./data/Germany/contacts/minimum_other.txt")

    minimum = np.loadtxt(minimum_contact_matrix0) \
        + np.loadtxt(minimum_contact_matrix1) + \
        np.loadtxt(minimum_contact_matrix2) + \
        np.loadtxt(minimum_contact_matrix3)

    return minimum


def get_population(path):
    """ read population data in list from dataset

    :param path: Path to the dataset containing the population 
    :returns: List of interpolated age grouped population data

    """

    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population.append(interpolate_age_groups(data_entry))
    return population


if __name__ == "__main__":
    # Store data relative to current file two levels higher.
    path = os.path.dirname(os.path.realpath(__file__))
    path_output = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    path_population = os.path.abspath(
        r"data//Germany//pydata//county_current_population.json")

    input_width = 5
    label_width = 90
    num_runs = 100
    data = generate_data(num_runs, path_output, path_population, input_width,
                         label_width, damping_method="active")
