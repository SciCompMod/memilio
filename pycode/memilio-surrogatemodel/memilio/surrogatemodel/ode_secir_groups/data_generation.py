#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

from memilio.simulation import (Damping, LogLevel, set_log_level)
from memilio.simulation.secir import (AgeGroup, Index_InfectionState,
                                      InfectionState, Model,
                                      interpolate_simulation_result, simulate)


def run_secir_simulation(days, damping_day, populations):
    """! Uses an ODE SECIR model allowing for asymptomatic infection with 6 different age groups. The model is not stratified by region. 
    Virus-specific parameters are fixed and initial number of persons in the particular infection states are chosen randomly from defined ranges.

    @param Days Describes how many days we simulate within a single run.
    @return List containing the populations in each compartment for each day of the simulation.
   """
    set_log_level(LogLevel.Off)

    start_day = 1
    start_month = 1
    start_year = 2019
    dt = 0.1

    # Define age Groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']
    num_groups = len(groups)

    # Initialize Parameters
    model = Model(num_groups)

    # Set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.IncubationTime[AgeGroup(i)] = 5.2
        model.parameters.TimeInfectedSymptoms[AgeGroup(i)] = 6.
        model.parameters.SerialInterval[AgeGroup(i)] = 4.2
        model.parameters.TimeInfectedSevere[AgeGroup(i)] = 12.
        model.parameters.TimeInfectedCritical[AgeGroup(i)] = 8.

        # Initial number of people in each compartment with random numbers
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Exposed)] = random.uniform(
            0.00025, 0.0005) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedNoSymptoms)] = random.uniform(
            0.0001, 0.00035) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSymptoms)] = random.uniform(
            0.00007, 0.0001) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedSevere)] = random.uniform(
            0.00003, 0.00006) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.InfectedCritical)] = random.uniform(
            0.00001, 0.00002) * populations[i]
        model.populations[AgeGroup(i), Index_InfectionState(
            InfectionState.Recovered)] = random.uniform(
            0.002, 0.008) * populations[i]
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

    model.parameters.StartDay = (
        date(start_year, start_month, start_day) - date(start_year, 1, 1)).days

    # Load baseline and minimum contact matrix and assign them to the model
    baseline = getBaselineMatrix()
    minimum = getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum

    # Generate a damping matrix and assign it to the model
    damping = np.ones((num_groups, num_groups)
                      ) * np.float16(random.uniform(0, 0.5))

    model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=(damping), t=damping_day, level=0, type=0))

    damped_contact_matrix = model.parameters.ContactPatterns.cont_freq_mat.get_matrix_at(
        damping_day+1)

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    # Run Simulation
    result = simulate(0, days, dt, model)
    # Interpolate simulation result on days time scale
    result = interpolate_simulation_result(result)

    # Using an array instead of a list to avoid problems with references
    result_array = result.as_ndarray()

    # Omit first column, as the time points are not of interest here.
    dataset_entries = copy.deepcopy(result_array[1:, :].transpose())

    return dataset_entries.tolist(), damped_contact_matrix


def generate_data(
        num_runs, path, input_width, label_width, normalize=True,
        save_data=True):
    """! Generate data sets of num_runs many equation-based model simulations and transforms the computed results by a log(1+x) transformation.
    Divides the results in input and label data sets and returns them as a dictionary of two TensorFlow Stacks.

    In general, we have 8 different compartments and 6 age groups.  If we choose, 
    input_width = 5 and label_width = 20, the dataset has 
    - input with dimension 5 x 8 x 6
    - labels with dimension 20 x 8 x 6

   @param num_runs Number of times, the function run_secir_simulation is called.
   @param path Path, where the dataset is saved to.
   @param input_width Int value that defines the number of time series used for the input.
   @param label_width Int value that defines the size of the labels.
   @param normalize [Default: true] Option to transform dataset by logarithmic normalization.
   @param save_data [Default: true] Option to save the dataset.
   @return Data dictionary of input and label data sets.
   """
    data = {
        "inputs": [],
        "labels": [],
        "contact_matrix": [],
        "damping_day": []
    }

    # The number of days is the same as the sum of input and label width.
    # Since the first day of the input is day 0, we still need to subtract 1.
    days = input_width + label_width - 1

    # Load population data
    population = get_population()

    # show progess in terminal for longer runs
    # Due to the random structure, theres currently no need to shuffle the data
    bar = Bar('Number of Runs done', max=num_runs)
    for _ in range(0, num_runs):

        # Generate a random damping day
        # The damping has to take place at least 5 days before the end of the simulation and should not be part of the input
        damping_day = random.randrange(
            input_width, (input_width+label_width)-5)

        data_run, damped_contact_matrix = run_secir_simulation(
            days, damping_day, population[random.randint(0, len(population) - 1)])
        data['inputs'].append(data_run[:input_width])
        data['labels'].append(data_run[input_width:])
        data['contact_matrix'].append(np.array(damped_contact_matrix))
        data['damping_day'].append(damping_day)
        bar.next()
    bar.finish()

    if normalize:
        # logarithmic normalization
        transformer = FunctionTransformer(np.log1p, validate=True)
        inputs = np.asarray(data['inputs']).transpose(2, 0, 1).reshape(8, -1)
        scaled_inputs = transformer.transform(inputs)
        scaled_inputs = scaled_inputs.transpose().reshape(num_runs, input_width, 48)
        scaled_inputs_list = scaled_inputs.tolist()

        labels = np.asarray(data['labels']).transpose(2, 0, 1).reshape(8, -1)
        scaled_labels = transformer.transform(labels)
        scaled_labels = scaled_labels.transpose().reshape(num_runs, label_width, 48)
        scaled_labels_list = scaled_labels.tolist()

        # cast dfs to tensors
        data['inputs'] = tf.stack(scaled_inputs_list)
        data['labels'] = tf.stack(scaled_labels_list)

    if save_data:
        # check if data directory exists. If necessary, create it.
        if not os.path.isdir(path):
            os.mkdir(path)

        # save dict to json file
        with open(os.path.join(path, 'data_secir_simple.pickle'), 'wb') as f:
            pickle.dump(data, f)
    return data


def getBaselineMatrix():
    """! loads the baselinematrix
    """

    baseline_contact_matrix0 = os.path.join(
        "./data/contacts/baseline_home.txt")
    baseline_contact_matrix1 = os.path.join(
        "./data/contacts/baseline_school_pf_eig.txt")
    baseline_contact_matrix2 = os.path.join(
        "./data/contacts/baseline_work.txt")
    baseline_contact_matrix3 = os.path.join(
        "./data/contacts/baseline_other.txt")

    baseline = np.loadtxt(baseline_contact_matrix0) \
        + np.loadtxt(baseline_contact_matrix1) + \
        np.loadtxt(baseline_contact_matrix2) + \
        np.loadtxt(baseline_contact_matrix3)

    return baseline


def getMinimumMatrix():
    """! loads the minimum matrix
    """

    minimum_contact_matrix0 = os.path.join(
        "./data/contacts/minimum_home.txt")
    minimum_contact_matrix1 = os.path.join(
        "./data/contacts/minimum_school_pf_eig.txt")
    minimum_contact_matrix2 = os.path.join(
        "./data/contacts/minimum_work.txt")
    minimum_contact_matrix3 = os.path.join(
        "./data/contacts/minimum_other.txt")

    minimum = np.loadtxt(minimum_contact_matrix0) \
        + np.loadtxt(minimum_contact_matrix1) + \
        np.loadtxt(minimum_contact_matrix2) + \
        np.loadtxt(minimum_contact_matrix3)

    return minimum


def get_population(path="/home/schm_a45/Documents/Code/memilio/memilio/data/pydata/Germany/county_population_dim401.json"):

    with open(path) as f:
        data = json.load(f)
    population = []
    for data_entry in data:
        population_county = []
        population_county.append(
            data_entry['<3 years'] + data_entry['3-5 years'] / 2)
        population_county.append(data_entry['6-14 years'])
        population_county.append(
            data_entry['15-17 years'] + data_entry['18-24 years'] +
            data_entry['25-29 years'] + data_entry['30-39 years'] / 2)
        population_county.append(
            data_entry['30-39 years'] / 2 + data_entry['40-49 years'] +
            data_entry['50-64 years'] * 2 / 3)
        population_county.append(
            data_entry['65-74 years'] + data_entry['>74 years'] * 0.2 +
            data_entry['50-64 years'] * 1 / 3)
        population_county.append(
            data_entry['>74 years'] * 0.8)

        population.append(population_county)
    return population


if __name__ == "__main__":
    # Store data relative to current file two levels higher.
    path = os.path.dirname(os.path.realpath(__file__))
    path_data = os.path.join(os.path.dirname(os.path.realpath(
        os.path.dirname(os.path.realpath(path)))), 'data')

    input_width = 5
    label_width = 30
    num_runs = 10000
    data = generate_data(num_runs, path_data, input_width,
                         label_width)
